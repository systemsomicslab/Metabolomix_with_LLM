"""注入ごとの RT / m/z 証拠（spec §9.1 `assay-feature-evidence.v1`）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

**なぜ「注入ごと」でなければならないか。** アライメント成果物は、特徴あたり1つの
代表 RT / m/z（SMF の `retention_time_in_seconds`、`.arf` スポットの
`MassCenter`）を持つ。これを各注入へ複製すると、RT のばらつきは定義上ゼロになり、
「全注入で標準品の RT が一致した」というQC合格が**計算せずに**出てしまう。実測値は
`AlignedPeakProperties` の行（= 注入1件）にしか無いので、読むのはそこだけにする。

**2つの軸をどう結ぶか。** どちらも名前では結ばない:

- feature 軸（`.arf` スポット ↔ mzTab SMF_ID）は、件数一致に加えて m/z の数値一致で
  確かめる（`mztab/evidence.py` の `compare_feature_mz` と同じ判定を共有する。
  1つずらした誤接合は Da オーダーに爆発する）。
- assay 軸（`.arf` の注入行 ↔ mzTab の assay）は、job の source→assay 対応
  (`ds.assay_sources`: abundance 列 → 実行時に予定した raw のパス) だけを使う。
  対応表が無いときに assay 表示名の一致へ落とさない——落とすと、v2 が要求する
  厳密な ID/source 照合が、既存の緩い推定で黙って代用される。

**照合できないときは表を作らない。** 部分的に埋まった証拠表は「取れた分だけ」で
QC を計算させてしまう。`availability=False` と `reasons` を返し、`rows` は空にする
（spec §9.1「対応成果物がない、または該当値を取得できない場合はmetricを
not_evaluableにする」）。値を取得できることと同定が確定していることは別で、
ここが答えるのは前者だけ。
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

from metabolomix.arf import reader as arf_reader
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.version import server_version
from metabolomix.mztab import evidence as mztab_evidence

__all__ = [
    "EVIDENCE_SCHEMA",
    "build_assay_evidence",
    "normalize_cell",
]

#: spec §9.1 の長形式表の schema 名。
EVIDENCE_SCHEMA = "assay-feature-evidence.v1"

#: この reader 自身の版。行ごとに保存する（同じ成果物から別版の reader が
#: 別の値を出したとき、どちらの表かを後から言えるようにする）。
_READER_NAME = "arf-assay-evidence"

#: 対応済みの RT 元単位 → 分への係数。ここに無い単位は推測しない。
_RT_FACTORS = {"minute": 1.0, "second": 1.0 / 60.0}

#: この module が読める `adapter["evidence_reader"]`。adapter が別の reader を
#: 宣言している場合は、その形式用の実装が要る（推測で `.arf` を読まない）。
_SUPPORTED_READER = "arf"


def normalize_cell(cell: dict, source_rt_unit: str) -> dict:
    """注入1件分のセルを、分に揃えた観測値へ落とす（純関数）。

    `source_rt_unit` が未対応なら `ValueError`。この層は例外を投げ、公開境界
    (`build_assay_evidence`) が `EVIDENCE_UNIT_UNSUPPORTED` の DomainError へ
    変換する——単位を既定値で補うと、秒の値が分として QC 窓へ入る。

    欠損には理由を付ける。`None` をそのまま返すだけだと、下流は「取得できなかった」
    と「測定されなかった」を区別できない。
    """
    factor = _RT_FACTORS.get(source_rt_unit)
    if factor is None:
        raise ValueError(
            f"unsupported RT unit: {source_rt_unit!r} "
            f"(supported: {sorted(_RT_FACTORS)})")

    rt = cell.get("rt")
    mz = cell.get("m_z")
    reasons = []
    if rt is None:
        reasons.append("rt_absent_in_artifact")
    if mz is None:
        reasons.append("mz_absent_in_artifact")

    return {
        "observed_rt_min": None if rt is None else float(rt) * factor,
        "observed_mz": None if mz is None else float(mz),
        "missing_reason": "+".join(reasons) if reasons else None,
    }


# ---------- 軸の照合 ----------

def _normalize_stem(value: str) -> str:
    """raw のパス→ファイル名（拡張子なし）の比較用正規化。

    `metabolomix/analysis/sample_manifest.py` と同じく `os.path.normcase` を使う
    （Windows の大文字小文字差でだけ対応が切れるのを避ける）。
    """
    return os.path.normcase(Path(value).stem)


def _assay_targets(ds) -> tuple[dict, list]:
    """assay_id → 予定 raw の stem。job の source→assay 対応だけから作る。"""
    assay_ids = list(getattr(ds, "sample_assay_ids", []) or [])
    sources = dict(getattr(ds, "assay_sources", {}) or {})
    targets: dict[str, str] = {}
    unmapped: list[str] = []
    for assay_id in assay_ids:
        raw = sources.get(f"abundance_{assay_id}")
        if not raw:
            unmapped.append(assay_id)
            continue
        targets[assay_id] = _normalize_stem(raw)
    return targets, unmapped


def _resolve_assay_columns(cells: list[dict], targets: dict) -> tuple[dict, list]:
    """1スポットのセル列を assay_id へ解決する。

    `.arf` の行順（FileID 順）は mzTab の assay 順と一致しないので、位置では
    結ばない。同じ raw に2行が対応したら `duplicate_assay_key` で止める
    ——feature×assay が一意でなければ、どちらの値が採られたか誰にも言えない。
    """
    by_stem: dict[str, list[int]] = {}
    for index, cell in enumerate(cells):
        name = cell.get("name")
        if not name:
            continue
        by_stem.setdefault(_normalize_stem(str(name)), []).append(index)

    columns: dict[str, int] = {}
    problems: list[dict] = []
    for assay_id, stem in targets.items():
        found = by_stem.get(stem, [])
        if len(found) > 1:
            problems.append({"code": "duplicate_assay_key",
                             "assay_id": assay_id, "stem": stem,
                             "rows": found})
        elif not found:
            problems.append({"code": "assay_source_unresolved",
                             "assay_id": assay_id, "stem": stem})
        else:
            columns[assay_id] = found[0]
    return columns, problems


def _reason(code: str, message: str, **details) -> dict:
    return {"code": code, "message": message, "details": details}


def _sha256(path: Path) -> str | None:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def _envelope(ds, artifact: Path, *, availability: bool, reasons: list[dict],
              rows: list[dict], artifact_hash: str | None, rt_unit: str,
              detail: dict | None = None) -> dict:
    n_detected = sum(1 for r in rows if r["detection_status"] == "detected")
    n_gap = sum(1 for r in rows if r["detection_status"] == "gap_filled")
    return {
        "schema": EVIDENCE_SCHEMA,
        "dataset_id": getattr(ds, "dataset_id", None),
        "availability": availability,
        "reasons": reasons,
        "source_artifact_hash": artifact_hash,
        "source_locator": str(artifact),
        "reader_version": f"{_READER_NAME}/{server_version()}",
        "rt_unit": {"source": rt_unit, "normalized": "minute"},
        "rows": rows,
        "coverage": {
            "features": len(getattr(ds, "feature_ids", []) or []),
            "assays": len(getattr(ds, "sample_assay_ids", []) or []),
            "rows": len(rows),
            "detected": n_detected,
            "gap_filled": n_gap,
            "values_missing": sum(1 for r in rows if r["observed_rt_min"] is None
                                  or r["observed_mz"] is None),
        },
        "detail": detail or {},
    }


def build_assay_evidence(ds, artifact: Path, adapter: dict) -> dict:
    """`.arf` 1本から `assay-feature-evidence.v1` を組む。

    `adapter` は `metabolomix.console.profile_adapter.adapter_capabilities()` の戻り値。
    どの形式を・どの単位で読んでよいかは adapter が宣言した範囲だけに限る
    （spec §5.1 と同じ規則。未確認の版の Key 配置を推測で読まない）。

    照合できない場合は例外ではなく封筒で返す——別アライメントの `.arf` が同じ
    フォルダに居るのは異常事態ではない。単位・reader の未対応だけは DomainError
    にする。これは入力の問題ではなく「この組み合わせ用の実装が無い」という
    計画不足で、黙って既定値へ落とすと数値が静かに狂うため。
    """
    artifact = Path(artifact)
    reader = adapter.get("evidence_reader")
    if reader != _SUPPORTED_READER:
        raise DomainError(
            "EVIDENCE_READER_UNSUPPORTED",
            f"このadapterが宣言するevidence_readerに対応する実装がありません: "
            f"{reader!r}（対応: {_SUPPORTED_READER!r}）。",
            {"evidence_reader": reader, "supported": [_SUPPORTED_READER],
             "adapter_id": adapter.get("adapter_id")},
        )
    rt_unit = adapter.get("evidence_rt_unit")
    if rt_unit not in _RT_FACTORS:
        raise DomainError(
            "EVIDENCE_UNIT_UNSUPPORTED",
            f"adapterが宣言するRT単位に対応していません: {rt_unit!r}"
            f"（対応: {sorted(_RT_FACTORS)}）。",
            {"evidence_rt_unit": rt_unit, "supported": sorted(_RT_FACTORS),
             "adapter_id": adapter.get("adapter_id")},
        )

    artifact_hash = _sha256(artifact)

    def fail(reasons, detail=None):
        return _envelope(ds, artifact, availability=False, reasons=reasons,
                         rows=[], artifact_hash=artifact_hash, rt_unit=rt_unit,
                         detail=detail)

    try:
        with open(artifact, "rb") as handle:
            spots = arf_reader.deserialize(handle)
    except Exception as exc:                      # noqa: BLE001 - 破損/別形式を包む
        return fail([_reason("artifact_unreadable",
                             "成果物を `.arf` として読めませんでした。",
                             path=str(artifact), error=str(exc))])

    normalized = mztab_evidence.normalize_arf_spots(spots)
    feature_ids = list(getattr(ds, "feature_ids", []) or [])

    if not normalized:
        # `deserialize` は壊れたファイルでも例外を出さず空を返す（別形式の
        # トップレベルオブジェクトは読み飛ばす仕様）。ブロック数で
        # 「読めなかった」と「読めたが注入行が無い」を分ける——前者は入力の
        # 取り違え、後者は成果物の種類違いで、読み手の次の一手が変わる。
        try:
            with open(artifact, "rb") as handle:
                blocks = arf_reader.readable_block_count(handle)
        except Exception:                         # noqa: BLE001
            blocks = 0
        if blocks == 0:
            return fail([_reason(
                "artifact_unreadable",
                "成果物を `.arf` として読めませんでした"
                "（LZ4 MessagePack ブロックが1件もありません）。",
                path=str(artifact))])
        return fail([_reason(
            "per_injection_rows_missing",
            "注入ごとのピーク行（AlignedPeakProperties）が1件もありません。"
            "代表値だけの成果物から各注入の値を複製することはしません。",
            n_spots=len(spots), n_blocks=blocks)])

    if len(normalized) != len(feature_ids):
        return fail([_reason(
            "feature_count_mismatch",
            "`.arf` のスポット数と mzTab の特徴数が一致しません。",
            n_arf_spots=len(normalized), n_mztab_features=len(feature_ids))])

    feature_meta = getattr(ds, "feature_metadata", {}) or {}
    mz_ok, mz_detail = mztab_evidence.compare_feature_mz(
        [spot.get("mz") for spot in normalized],
        [(feature_meta.get(fid) or {}).get("mz") for fid in feature_ids])
    if not mz_ok:
        return fail([_reason(
            "feature_axis_mz_mismatch",
            "スポット順と SMF_ID の対応が m/z で確認できません（誤接合の疑い）。",
            **mz_detail)], detail=mz_detail)

    targets, unmapped = _assay_targets(ds)
    if unmapped or not targets:
        return fail([_reason(
            "assay_source_unavailable",
            "job の source→assay 対応（assay_sources）が無いため、注入を"
            "assay へ結べません。assay 表示名の一致では代用しません。",
            unmapped_assays=unmapped)])

    # 列の解決は**全スポットで同じ**であることを要求する。スポットごとに
    # 別の行が同じ assay に当たるなら、その `.arf` の注入軸は一貫していない。
    base_columns, problems = _resolve_assay_columns(normalized[0]["samples"], targets)
    if problems:
        codes = sorted({p["code"] for p in problems})
        return fail([_reason(
            code,
            "注入行を assay へ一意に解決できません。",
            problems=[p for p in problems if p["code"] == code])
            for code in codes])

    rows: list[dict] = []
    for fid, spot in zip(feature_ids, normalized):
        cells = spot["samples"]
        columns, problems = _resolve_assay_columns(cells, targets)
        if problems:
            codes = sorted({p["code"] for p in problems})
            return fail([_reason(
                code,
                f"特徴 {fid!r} で注入行を assay へ一意に解決できません。",
                feature_id=fid,
                problems=[p for p in problems if p["code"] == code])
                for code in codes])
        for assay_id, index in columns.items():
            cell = cells[index]
            values = normalize_cell(cell, rt_unit)
            rows.append({
                "dataset_id": getattr(ds, "dataset_id", None),
                "feature_id": fid,
                "assay_id": assay_id,
                "observed_rt_min": values["observed_rt_min"],
                "observed_mz": values["observed_mz"],
                # 検出状態は gap-fill フラグそのもの。値が非ゼロであることを
                # 検出の代用にしない（mzTab の非ゼロの 70% が gap-fill だった）。
                "detection_status": ("gap_filled" if cell.get("is_gap_filled")
                                     else "detected"),
                "missing_reason": values["missing_reason"],
                "source_artifact_hash": artifact_hash,
                "source_locator": str(artifact),
                "reader_version": f"{_READER_NAME}/{server_version()}",
            })

    return _envelope(
        ds, artifact, availability=True, reasons=[], rows=rows,
        artifact_hash=artifact_hash, rt_unit=rt_unit,
        detail={**mz_detail,
                "feature_axis": "arf_spot_index == mztab_smf_id (m/z verified)",
                "assay_axis": "arf_file_name == job source raw stem",
                "assay_columns": base_columns})
