"""実験情報シート（sample-manifest.v1 / v2）の厳密読込・対応・出所・原子的適用（spec §7）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。ここで縛るのは
3 段階の責務分離と、その上に載る**統計母集団の選択**:

``parse_manifest``
    シート自体の厳密検証。列の過不足・重複ID・未知source_file・行欠落・
    親パスへの脱出・不正role・不正注入順・batch内の注入順重複を**1件でも**
    見つけたら、部分的な行リストを返さずに例外にする。明示されたセルだけへ
    `provenance={value, source: "user_manifest", confidence: "confirmed"}` を
    付ける（空欄セルは `{value: None, source: None, confidence: None}` のまま
    `resolve_metadata` へ回す）。

    ``source_file`` は解決後の絶対パス文字列として保持する（シート上は
    `source_root` 基準の相対パスだが、下流の `resolve_metadata` は
    `ds.assay_sources`（絶対パス）とだけ照合するため、ここで解決しておかないと
    `source_root` を持ち回る必要が生じる）。

``resolve_metadata``
    明示シート（あれば）と dataset 由来の推定を統合し、``ds.sample_names`` と
    同じ順に並べて返す。シートの行順・mzTab の assay 列順は一致していなくて
    よい——``ds.assay_sources``（`assay[N]-ms_run_ref` → `ms_run[N]-location` を
    辿った「どの列がどの生データか」の対応表、Task 7）で結ぶ。

    role・batch・injection_order・qc_pool のうち明示値が無いセルだけ推定で
    埋める。ただし既存 `detect_sample_roles`（ファイル名トークン判定）は
    **シートが全く無い場合にだけ**使う——シートがある以上、blank セルは
    トークン推定へ回さず「sample」扱いに留める（部分入力シートでの
    per-セルトークン推定は本 task の範囲外。7.3 のフル運用は次版の課題）。
    明示 `role=unknown` は正当な列挙値であり、sample へ再推定しない。
    group・qc_pool は常に明示値のみを使い、欠落を推測で埋めない
    （「全QCは同じpool」に化けさせない）。batch・injection_order は
    mzTab MTD の `assay[N]-custom[...]`（run_order/batch）を候補として使えるが、
    出所を `source="mztab", confidence="unverified"` に留める——Console の
    読込順由来か実注入順由来かをここでは確認できないため、自動ドリフト補正の
    根拠にしない。明示値と異なる候補は上書きせず ``conflicts`` に残す。

``apply_metadata``
    ``rows`` が ``ds.sample_names`` と同じ順・同じ件数で並んでいる前提で、
    対応（件数・sample_id の一意性）と型を全件検証してから、Task 6
    (`result_state.metadata_fingerprints`) の指紋比較で「どの列が変わったか」を
    求め、`invalidate_results` を一度だけ呼んでから ``ds`` へ一括代入する。
    検証で1件でも落ちれば ``ds`` には一切触れない（半端な適用を作らない）。

``select_statistical_samples`` / ``validate_independent_samples`` / ``validate_standard_assays``
    v2 が足す2点——``role=standard``（標準品注入は生物試料ではない）と
    ``biological_sample_id``（反復注入を束ねる生物単位）——を、検定の**直前**に
    効かせる。台帳としてのシートは反復注入も標準品も書けてよく、止めるのは
    「その集合をnと数えてよいか」を問う場所だけ。選択規則を1か所に集めるのは、
    統計・QC・標準注入解決が各々書き直すと片方だけ ``include=false`` を
    見落とすため。
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path

from metabolomix.analysis import preprocessing
from metabolomix.analysis.result_state import invalidate_results, metadata_fingerprints
from metabolomix.core.atomic_io import DomainError

__all__ = [
    "FIELDS",
    "FIELDS_V2",
    "apply_metadata",
    "parse_injection_order",
    "parse_manifest",
    "resolve_metadata",
    "select_statistical_samples",
    "validate_independent_samples",
    "validate_standard_assays",
]

#: sample-manifest.v1 が持つ8列（この順で並ぶことを要求する）。
FIELDS = ("sample_id", "source_file", "role", "group", "batch",
          "injection_order", "qc_pool", "include")

#: sample-manifest.v2。v1の8列順をそのまま保ち、末尾に生物単位を足すだけ
#: （spec §7）。列を挿入せず末尾に足すのは、v1の読み書きをする既存経路と
#: 見出し行の突き合わせが位置で壊れないようにするため。
FIELDS_V2 = FIELDS + ("biological_sample_id",)

_SCHEMA_HEADERS = {
    "# schema = sample-manifest.v1": FIELDS,
    "# schema = sample-manifest.v2": FIELDS_V2,
}
_ROLE_VALUES = frozenset({"sample", "qc", "blank", "unknown"})
#: v2で足す `standard`（標準品注入）は測定の道具であって生物試料ではない。
#: v1のシートには書けない——v1の読み手はこの値の意味を知らないため、
#: 「未知のrole」として素通りさせず読込時に止める。
_ROLE_VALUES_V2 = _ROLE_VALUES | {"standard"}
_INJECTION_ORDER_RE = re.compile(r"[1-9][0-9]*")
_ASSAY_ID_RE = re.compile(r"^assay\[(\d+)\]$")

# csv.DictReader の余剰列を捨てずに検知するための番兵キー。
_EXTRA_KEY = "__extra_columns__"


def parse_injection_order(text: str) -> int | None:
    """注入順セルを厳密パースする。

    空文字はNone（欠落値）。それ以外は`[1-9][0-9]*`のみを受け付け、"NaN"
    "Infinity"「1.5」「0」「-1」のような数値『らしい』文字列を浮動小数点や
    0/負数として黙って通さない——batch内一意性チェックが意味を失うため。
    """
    if text == "":
        return None
    if not _INJECTION_ORDER_RE.fullmatch(text):
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            f"注入順は正の整数です: {text!r}",
            {"value": text},
        )
    return int(text)


def _normalize_path(value: str) -> str:
    """パス比較専用の正規化。`metabolomix.console.validation._normalize_path`と同じ規則
    （層を跨いだprivate参照を避けるため複製する）。
    """
    return os.path.normcase(os.path.normpath(value))


def parse_manifest(path: Path, *, source_root: Path,
                    expected_sources: list[str]) -> list[dict]:
    """sample-manifest を厳密に読み、列+provenance+conflictsの行リストを返す。

    版は**先頭のschema行**で決まる（v1=8列、v2=9列）。列見出しから版を推測しない
    ——推測すると、v1と書いてあるシートへ9列目を足しただけで v2 の意味論
    （`standard` role・生物単位）が黙って有効になる。

    `expected_sources`は`source_root`基準の相対パス一覧（実行時に予定した
    raw測定単位）。全件に1行ずつ対応することを要求し、未知のsource_file・
    行欠落・重複を許さない。1件でも不正なら例外にし、部分的な行リストを
    返さない。
    """
    path = Path(path)
    source_root = Path(source_root).resolve()
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    schema_line = lines[0].rstrip("\n") if lines else ""
    fields = _SCHEMA_HEADERS.get(schema_line)
    if fields is None:
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            f"先頭行が既知のschema宣言ではありません: {schema_line!r}"
            f"（対応: {', '.join(sorted(_SCHEMA_HEADERS))}）。",
            {"path": str(path), "schema_line": schema_line,
             "supported": sorted(_SCHEMA_HEADERS)},
        )
    role_values = _ROLE_VALUES_V2 if fields is FIELDS_V2 else _ROLE_VALUES

    reader = csv.DictReader(lines[1:], delimiter="\t",
                            restkey=_EXTRA_KEY, restval=None)
    fieldnames = list(reader.fieldnames or [])
    if fieldnames != list(fields):
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            "列見出しが仕様と一致しません（過不足または順序違い）。",
            {"expected": list(fields), "found": fieldnames},
        )

    # 予定raw（相対パス）を正規化して索引化。全件へ1行ずつ来ることを後で照合する。
    expected_by_norm: dict[str, str] = {}
    for rel in expected_sources:
        resolved = (source_root / rel).resolve()
        expected_by_norm[_normalize_path(str(resolved))] = rel

    rows: list[dict] = []
    seen_ids: set[str] = set()
    seen_raw_norm: dict[str, str] = {}

    for line_no, raw_row in enumerate(reader, start=3):  # 1行目schema, 2行目見出し
        if raw_row.get(_EXTRA_KEY) is not None:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"{line_no}行目に列数超過があります。",
                {"line": line_no},
            )
        if any(raw_row.get(field) is None for field in fields):
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"{line_no}行目の列数が不足しています。",
                {"line": line_no},
            )

        sample_id = raw_row["sample_id"]
        if not sample_id:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"{line_no}行目: sample_id が空です。",
                {"line": line_no},
            )
        if sample_id in seen_ids:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"sample_id が重複しています: {sample_id!r}",
                {"sample_id": sample_id},
            )
        seen_ids.add(sample_id)

        source_file_text = raw_row["source_file"]
        if not source_file_text:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"{line_no}行目: source_file が空です。",
                {"line": line_no, "sample_id": sample_id},
            )
        resolved = (source_root / source_file_text).resolve()
        try:
            resolved.relative_to(source_root)
        except ValueError:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"source_file が元フォルダの外を指しています: {source_file_text!r}",
                {"sample_id": sample_id, "source_file": source_file_text},
            ) from None
        norm = _normalize_path(str(resolved))
        if norm not in expected_by_norm:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"未知の source_file です（予定raw一覧にありません）: "
                f"{source_file_text!r}",
                {"sample_id": sample_id, "source_file": source_file_text},
            )
        if norm in seen_raw_norm:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"同一rawへ複数行があります: {source_file_text!r} と "
                f"{seen_raw_norm[norm]!r}",
                {"source_file": source_file_text},
            )
        seen_raw_norm[norm] = source_file_text

        role = raw_row["role"] or None
        if role is not None and role not in role_values:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"role が不正です: {role!r}",
                {"sample_id": sample_id, "role": role},
            )

        group = raw_row["group"] or None
        batch = raw_row["batch"] or None
        qc_pool = raw_row["qc_pool"] or None
        injection_order = parse_injection_order(raw_row["injection_order"])

        include_text = raw_row["include"]
        if include_text == "":
            include, include_confirmed = True, False
        elif include_text == "true":
            include, include_confirmed = True, True
        elif include_text == "false":
            include, include_confirmed = False, True
        else:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"include は true/false のみです: {include_text!r}",
                {"sample_id": sample_id, "include": include_text},
            )

        values = {
            "sample_id": sample_id, "source_file": str(resolved), "role": role,
            "group": group, "batch": batch, "injection_order": injection_order,
            "qc_pool": qc_pool, "include": include,
        }
        if fields is FIELDS_V2:
            # 空欄は「生物学的独立性が未確認」（spec §7）。ここでは欠落として
            # 通し、検定の直前 `validate_independent_samples` で止める——
            # シートは標準品・blankのように生物単位を持たない行も載せるため。
            values["biological_sample_id"] = raw_row["biological_sample_id"] or None
        provenance: dict[str, dict] = {}
        for field in fields:
            if field == "include":
                provenance[field] = (
                    {"value": include, "source": "user_manifest", "confidence": "confirmed"}
                    if include_confirmed else
                    {"value": True, "source": "default", "confidence": "inferred"}
                )
                continue
            value = values[field]
            provenance[field] = (
                {"value": None, "source": None, "confidence": None} if value is None else
                {"value": value, "source": "user_manifest", "confidence": "confirmed"}
            )

        row = dict(values)
        row["provenance"] = provenance
        row["conflicts"] = []
        rows.append(row)

    missing = [rel for norm, rel in expected_by_norm.items() if norm not in seen_raw_norm]
    if missing:
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            f"予定raw {len(missing)} 件に行がありません: {', '.join(missing)}",
            {"missing_sources": missing},
        )

    # batch内（blank=未確認batchも1グループとして扱う）の注入順重複を拒否する。
    by_batch: dict[str | None, list[int]] = {}
    for row in rows:
        if row["injection_order"] is None:
            continue
        by_batch.setdefault(row["batch"], []).append(row["injection_order"])
    for batch_key, orders in by_batch.items():
        if len(orders) != len(set(orders)):
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"batch内で注入順が重複しています（batch={batch_key!r}）。",
                {"batch": batch_key},
            )

    return rows


# ---------- resolve_metadata ----------

def _abundance_key(assay_id: str | None) -> str | None:
    """`assay[N]` → `abundance_assay[N]`（`ds.assay_sources`のキー形式）。"""
    m = _ASSAY_ID_RE.fullmatch(assay_id or "")
    return f"abundance_assay[{m.group(1)}]" if m else None


def _raw_path_for_sample(ds, index: int) -> str | None:
    """i番目サンプルへ実行時に予定されたrawの絶対パス（分かれば）。"""
    assay_ids = list(getattr(ds, "sample_assay_ids", []) or [])
    if index >= len(assay_ids):
        return None
    key = _abundance_key(assay_ids[index])
    if key is None:
        return None
    return (getattr(ds, "assay_sources", {}) or {}).get(key)


def _mztab_batch_and_order(ds, index: int) -> tuple[str | None, int | None]:
    """mzTab MTDのassay[N]-custom[...]から拾ったbatch/注入順候補（未確認）。"""
    assay_ids = list(getattr(ds, "sample_assay_ids", []) or [])
    if index >= len(assay_ids):
        return None, None
    meta = (getattr(ds, "assay_metadata", {}) or {}).get(assay_ids[index]) or {}
    return meta.get("batch"), meta.get("run_order")


def _unverified_mztab_prov(value) -> dict:
    if value is None:
        return {"value": None, "source": None, "confidence": None}
    return {"value": value, "source": "mztab", "confidence": "unverified"}


def _default_row(ds, index: int, name: str) -> dict:
    """明示シートが全く無いときの1行。既存のトークン/mzTab推定経路を使う。"""
    source_file = _raw_path_for_sample(ds, index)
    role = preprocessing.detect_sample_roles([name]).get(name, "sample")
    batch, injection_order = _mztab_batch_and_order(ds, index)

    values = {
        "sample_id": name, "source_file": source_file, "role": role,
        "group": None, "batch": batch, "injection_order": injection_order,
        "qc_pool": None, "include": True,
    }
    provenance = {
        "sample_id": {"value": name, "source": "default", "confidence": "inferred"},
        "source_file": (_unverified_mztab_prov(source_file) if source_file
                        else {"value": None, "source": None, "confidence": None}),
        "role": {"value": role, "source": "filename_token", "confidence": "inferred"},
        "group": {"value": None, "source": None, "confidence": None},
        "batch": _unverified_mztab_prov(batch),
        "injection_order": _unverified_mztab_prov(injection_order),
        "qc_pool": {"value": None, "source": None, "confidence": None},
        "include": {"value": True, "source": "default", "confidence": "inferred"},
    }
    row = dict(values)
    row["provenance"] = provenance
    row["conflicts"] = []
    return row


def _merge_row(ds, index: int, name: str, row: dict) -> dict:
    """明示シートの1行へ、埋まっていない列だけdataset由来の候補で補う。"""
    values = dict(row)
    provenance = {field: dict(prov) for field, prov in row["provenance"].items()}
    conflicts = list(row.get("conflicts", []))

    if values["role"] is None:
        # シートがある以上、blankをトークン推定へは回さず「sample」扱いに留める
        # （明示role=unknownは既に非Noneなのでここへは来ない＝再推定されない）。
        values["role"] = "sample"
        provenance["role"] = {"value": "sample", "source": "default", "confidence": "inferred"}

    mz_batch, mz_order = _mztab_batch_and_order(ds, index)
    if values["batch"] is None:
        if mz_batch is not None:
            values["batch"] = mz_batch
            provenance["batch"] = _unverified_mztab_prov(mz_batch)
    elif mz_batch is not None and mz_batch != values["batch"]:
        conflicts.append({"field": "batch", **_unverified_mztab_prov(mz_batch)})

    if values["injection_order"] is None:
        if mz_order is not None:
            values["injection_order"] = mz_order
            provenance["injection_order"] = _unverified_mztab_prov(mz_order)
    elif mz_order is not None and mz_order != values["injection_order"]:
        conflicts.append({"field": "injection_order", **_unverified_mztab_prov(mz_order)})

    merged = dict(values)
    merged["provenance"] = provenance
    merged["conflicts"] = conflicts
    return merged


def resolve_metadata(ds, rows: list[dict] | None) -> list[dict]:
    """明示シート（あれば）とdataset由来の既定値を統合し、ds.sample_namesの順で返す。

    `rows`が`None`なら明示シートが無いケース（既存のdetect_sample_roles/mzTab
    推定だけを使う）。`rows`があれば`ds.assay_sources`の対応表で各行をサンプルへ
    結び付ける——シートの行順やdataset内の並びは一致していなくてよい。
    """
    sample_names = list(getattr(ds, "sample_names", []) or [])

    if rows is None:
        return [_default_row(ds, i, name) for i, name in enumerate(sample_names)]

    by_norm_path = {_normalize_path(row["source_file"]): row for row in rows}

    resolved: list[dict] = []
    for i, name in enumerate(sample_names):
        raw_path = _raw_path_for_sample(ds, i)
        row = by_norm_path.get(_normalize_path(raw_path)) if raw_path else None
        if row is None:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"サンプル {name!r} に対応するシート行が見つかりません"
                "（raw対応表で解決できませんでした）。",
                {"sample_name": name, "raw_path": raw_path},
            )
        resolved.append(_merge_row(ds, i, name, row))
    return resolved


# ---------- apply_metadata ----------

_TYPE_CHECKS: dict[str, tuple[type, ...]] = {
    "sample_id": (str,), "source_file": (str, type(None)),
    "role": (str, type(None)), "group": (str, type(None)),
    "batch": (str, type(None)), "injection_order": (int, type(None)),
    "qc_pool": (str, type(None)), "include": (bool,),
    "biological_sample_id": (str, type(None)),
}


def _validate_rows(rows: list[dict], sample_names: list[str]) -> None:
    """対応（件数・ID一意性）と型を全件検証する。1件でも落ちればここで止める。"""
    if len(rows) != len(sample_names):
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            f"行数がサンプル数と一致しません（rows={len(rows)}, "
            f"samples={len(sample_names)}）。",
            {"rows": len(rows), "samples": len(sample_names)},
        )
    seen_ids: set[str] = set()
    for index, row in enumerate(rows):
        missing_keys = [f for f in FIELDS if f not in row]
        if missing_keys or "provenance" not in row or "conflicts" not in row:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"{index}行目の形が不正です（欠けているキー: {missing_keys}）。",
                {"index": index, "missing": missing_keys},
            )
        for field, allowed_types in _TYPE_CHECKS.items():
            if field not in row:  # biological_sample_id はv1行に無い。
                continue
            value = row[field]
            if not isinstance(value, allowed_types):
                raise DomainError(
                    "SAMPLE_MANIFEST_INVALID",
                    f"{index}行目の {field} の型が不正です: {value!r}",
                    {"index": index, "field": field},
                )
        allowed_roles = (_ROLE_VALUES_V2 if "biological_sample_id" in row
                         else _ROLE_VALUES)
        if row["role"] is not None and row["role"] not in allowed_roles:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"{index}行目のroleが不正です: {row['role']!r}",
                {"index": index, "role": row["role"]},
            )
        sample_id = row["sample_id"]
        if not sample_id or sample_id in seen_ids:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"sample_id が空または重複しています: {sample_id!r}",
                {"index": index, "sample_id": sample_id},
            )
        seen_ids.add(sample_id)


def apply_metadata(ds, rows: list[dict]) -> dict:
    """厳密検証してから一括反映する（部分適用なし）。

    `rows`は`ds.sample_names`と同じ順・同じ件数で並んでいる前提
    （`resolve_metadata`の出力、またはそれと同じ契約の合成データ）。対応と型を
    全件検証してから、Task 6の指紋比較（`metadata_fingerprints`）で
    「何が変わったか」を求め、`invalidate_results`を一度だけ呼ぶ。検証で
    1件でも落ちれば`ds`には一切触れない。
    """
    sample_names = list(getattr(ds, "sample_names", []) or [])
    _validate_rows(rows, sample_names)

    old_rows = getattr(ds, "sample_metadata_rows", None) or []
    if not old_rows:
        # 初回適用は「全列が確定した」ことそのものが前処理入力の変化にあたる。
        changed = {field for field in FIELDS_V2
                   if field in FIELDS or any(field in row for row in rows)}
    else:
        before = metadata_fingerprints(old_rows)
        after = metadata_fingerprints(rows)
        # 比較対象はv1の8列ではなくv2の全列。v1の列だけで畳むと、v2で足した
        # `biological_sample_id` の付け替えが「何も変わっていない」と読まれ、
        # 旧IDのnで出した検定結果が生き残る（列がv1行に無い場合は新旧とも
        # 欠落＝同値なので、v1経路の changed_fields は従来どおり）。
        changed = {field for field in FIELDS_V2
                   if before.get(field) != after.get(field)}

    # ここまでの検証が全部通ってから、まとめて反映する（半端な適用を作らない）。
    invalidate_results(ds, changed)
    ds.sample_metadata_rows = rows
    ds.sample_ids = [row["sample_id"] for row in rows]
    ds.metadata_revision = getattr(ds, "metadata_revision", 0) + 1

    return {
        "metadata_revision": ds.metadata_revision,
        "changed_fields": sorted(changed),
        "sample_ids": list(ds.sample_ids),
    }


# ---------- 統計母集団の選択（spec §7・§9） ----------

def select_statistical_samples(rows: list[dict], groups: list[str]) -> list[dict]:
    """検定の対象になる行だけを、`rows`の順のまま返す。

    条件は3つのANDで、ここが**唯一の定義**（spec §9「全検定はincluded
    role=sampleのみを対象とし」）:

    - ``include`` が真——除外した注入はQCからも統計からも消える（raw対応一覧
      からは消さない。それは`parse_manifest`が持つ台帳の役目）。
    - ``role`` が ``"sample"``——``standard``（標準品注入）・``qc``・``blank``・
      ``unknown`` は生物試料ではない。standardを群へ自動投入しない、が §7 の要求。
    - ``group`` が引数 ``groups`` に在る——空欄の行は群未確定であって「その他群」
      ではないので、どの比較にも入れない。

    同じ3条件を統計・QC・標準注入解決が各々書き直すと、片方だけ
    ``include=false`` を見落として母集団が食い違う。呼び出し側はこの関数を通す。
    """
    wanted = set(groups)
    return [row for row in rows
            if row.get("include", True)
            and row.get("role") == "sample"
            and row.get("group") in wanted]


def validate_independent_samples(rows: list[dict], groups: list[str]) -> None:
    """検定対象の各注入が、独立した生物試料であることを確認する。

    注入数をそのままnと読むと、同一個体の2注入が「n=2」に化けて群内分散を
    過小評価する。初期版は自動平均も混合モデルも持たないので、対応できない
    ことを ``REPEATED_MEASURES_UNSUPPORTED`` で明示して止める（spec §7）。

    ``biological_sample_id`` の空欄・列自体が無いv1行は「独立性が未確認」で
    あって「独立」ではない。ここで推測して通すと、同じ個体の反復注入が
    区別されないまま n に化ける——止めて、シートで宣言させる。

    重複を見るのは**検定対象だけ**。除外行・standard・対象外群のIDが重なって
    いても検定のnには入らないので止めない。
    """
    selected = select_statistical_samples(rows, groups)

    missing = [row["sample_id"] for row in selected
               if not row.get("biological_sample_id")]
    if missing:
        raise DomainError(
            "BIOLOGICAL_SAMPLE_ID_REQUIRED",
            f"検定対象 {len(missing)} 件に biological_sample_id がありません"
            f"（sample-manifest.v2 で明示してください）: {', '.join(missing)}",
            {"sample_ids": missing, "groups": list(groups)},
        )

    by_id: dict[str, list[str]] = {}
    for row in selected:
        by_id.setdefault(row["biological_sample_id"], []).append(row["sample_id"])
    repeated = {bio_id: ids for bio_id, ids in by_id.items() if len(ids) > 1}
    if repeated:
        raise DomainError(
            "REPEATED_MEASURES_UNSUPPORTED",
            "同一 biological_sample_id の注入が検定対象に複数あります"
            "（反復測定モデル・自動平均は初期版の対象外です）: "
            + ", ".join(f"{bio_id}={'/'.join(ids)}"
                        for bio_id, ids in sorted(repeated.items())),
            {"repeated": {k: sorted(v) for k, v in sorted(repeated.items())},
             "groups": list(groups)},
        )


def validate_standard_assays(rows: list[dict],
                             standard_assays: dict[str, list[str]]) -> None:
    """``standard_assays`` の各sample_idが、実在する採用済み標準注入か確認する。

    要求側（`pipeline/request_v2.py`）が見られるのは「target_idが在るか」「形が
    正しいか」までで、sample_idが本当に標準品注入かはシートと突き合わせないと
    分からない（spec §6.2「role=standardかつinclude=trueであることをmanifestと
    照合する」）。ここを緩めると、生物試料を内部標準の分母に使った比が
    ``internal_standard_ratio`` の名前で出て行く。
    """
    by_id = {row["sample_id"]: row for row in rows}
    bad: dict[str, str] = {}
    for target_id, sample_ids in (standard_assays or {}).items():
        for sample_id in sample_ids:
            row = by_id.get(sample_id)
            if row is None:
                bad[f"{target_id}:{sample_id}"] = "not_in_manifest"
            elif row.get("role") != "standard":
                bad[f"{target_id}:{sample_id}"] = f"role={row.get('role')!r}"
            elif not row.get("include", True):
                bad[f"{target_id}:{sample_id}"] = "include=false"
    if bad:
        raise DomainError(
            "STANDARD_ASSAY_INVALID",
            "standard_assays が標準品注入（role=standard・include=true）を"
            "指していません: "
            + ", ".join(f"{key}（{reason}）" for key, reason in sorted(bad.items())),
            {"invalid": dict(sorted(bad.items()))},
        )
