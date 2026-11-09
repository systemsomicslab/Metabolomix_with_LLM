"""完了した run の成果物を、対話セッションが使える DatasetState へ戻す層。

`load_dataset(pipeline_path=..., polarity=...)`（`dataset.opening`）が `resolve_unit` で
極性の単位を 1 つ選び、その単位の解析フォルダ・Console のジョブ・試料の対応表を引く。
v3 の run は最後の工程が書いた受領（`receipt.read_receipt`）から読み、受領が無い・照合できない
ときだけ run の記録の単位から読む。旧 v1 の run は読むだけの記録（`store.is_legacy`）で、
1 単位として今までどおり開く。`source="mztab"` で読むとき（`load_run_dataset`）、差次的結果を
戻すのは旧 v1 の run だけ。

グローバル session は import しない（この層は DatasetState を作って返すだけで、
どこへ置くかは呼び出し側＝MCP 公開層が決める）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from metabolomix.analysis import result_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline import recovery, store

__all__ = ["UnitRun", "load_run_dataset", "resolve_unit", "run_manifest_rows"]

#: `service._persist_differential` が差次的結果に付ける output_name の接頭辞。
_DIFFERENTIAL_PREFIX = "differential:"

#: `service._handle_preprocess` が前処理の成果物に付ける output_name。
_PREPROCESS_OUTPUT = "preprocess"

_INPUT_SUBDIR = "input"


@dataclass(frozen=True)
class UnitRun:
    """run の極性の単位 1 つを開くのに要るもの（パスはすべて絶対）。"""

    pipeline_root: Path
    record: dict
    polarity: str
    legacy: bool
    unit_root: Path
    #: MS-DIAL が `.arf` 等を書いた解析フォルダ（受領の `analysis_dir`。無ければ `<単位>/input`）。
    analysis_dir: Path
    job_path: Path
    #: 受領の時点で解いた試料の対応表の行（無ければ None）。
    manifest_rows: list[dict] | None
    notes: list[str] = field(default_factory=list)


def _polarity_options(record: dict, units: dict) -> list[dict]:
    return [{"value": p, "label": f"{p}（{store.unit_status(record, p)}）"}
            for p in units if store.unit_status(record, p) == "validated"]


def resolve_unit(pipeline_path, polarity: str | None = None, *,
                 allow_incomplete: bool = False) -> UnitRun:
    """run の極性の単位を 1 つ選び、その解析フォルダ・ジョブ・試料の対応表を返す。

    `pipeline_path` は `pipeline-run.json` でもその親ディレクトリでもよい。`polarity` を省いて
    単位が 1 つならそれ、2 つ以上なら `PIPELINE_POLARITY_REQUIRED`（`details.options` が選べる
    極性）。v3 の単位は検証を通ったもの（`validated`）だけ開き、それ以外は
    `PIPELINE_UNIT_NOT_VALIDATED`（`allow_incomplete` なら探索用に読む）。
    """
    root = recovery.normalize_pipeline_root(Path(pipeline_path).expanduser())
    record = store.load_run(root)
    legacy = store.is_legacy(record)
    units = store.run_units(record)
    if not units:
        raise DomainError(
            "PIPELINE_UPSTREAM_INCOMPLETE",
            "この run はまだ上流（Console）の成果物を持っていません。"
            "pipeline_status で状態を確認してください。", {"pipeline_root": str(root)})

    if polarity is None:
        if len(units) > 1:
            options = _polarity_options(record, units)
            if not options and not allow_incomplete:
                raise _not_validated(root, record, units, None)
            raise DomainError(
                "PIPELINE_POLARITY_REQUIRED",
                "この run には極性の単位が複数あります。どれを開くか polarity で選んでください。",
                {"options": options or [{"value": p, "label": p} for p in units]})
        polarity = next(iter(units))
    elif polarity not in units:
        raise DomainError(
            "PIPELINE_POLARITY_NOT_FOUND",
            f"この run に極性 {polarity!r} の単位はありません（{', '.join(units)}）。",
            {"polarity": polarity, "available": list(units)})

    if not legacy and store.unit_status(record, polarity) != "validated" and not allow_incomplete:
        raise _not_validated(root, record, units, polarity)

    unit = units[polarity]
    receipt = None if legacy else receipt_mod.read_receipt(root, record)
    unit_receipt = ((receipt or {}).get("units") or {}).get(polarity)
    notes: list[str] = []
    if not legacy and unit_receipt is None:
        notes.append("run の受領（results/receipt.json）が無い・照合できないので、run の記録の"
                     "単位から解析フォルダを引いた（試料の対応表は載らない）。")

    unit_root = root / (unit.get("unit_root") or ".")
    job_rel = ((unit_receipt or {}).get("console_job")
               or (unit.get("upstream") or {}).get("console_job_path"))
    if not job_rel:
        raise DomainError(
            "PIPELINE_UPSTREAM_INCOMPLETE",
            "この run はまだ上流（Console）の成果物を持っていません。"
            "pipeline_status で状態を確認してください。",
            {"pipeline_root": str(root), "polarity": polarity})
    if unit_receipt is not None and unit_receipt.get("analysis_dir"):
        analysis_dir = root / unit_receipt["analysis_dir"]
    else:
        analysis_dir = unit_root / _INPUT_SUBDIR

    if legacy:
        rows = run_manifest_rows(root, record)
    else:
        rows = _receipt_rows(root, unit_receipt)
    return UnitRun(pipeline_root=root, record=record, polarity=polarity, legacy=legacy,
                   unit_root=unit_root, analysis_dir=analysis_dir, job_path=root / job_rel,
                   manifest_rows=rows, notes=notes)


def _not_validated(root: Path, record: dict, units: dict, polarity: str | None) -> DomainError:
    statuses = {p: store.unit_status(record, p) for p in units}
    which = f"極性 {polarity!r} の単位は" if polarity else "どの単位も"
    return DomainError(
        "PIPELINE_UNIT_NOT_VALIDATED",
        f"{which}検証を通っていません（{statuses}）。pipeline_status で確かめてください。"
        "中断・失敗した実行の出力を探索で読むなら allow_incomplete=True。",
        {"pipeline_root": str(root), "polarity": polarity, "units": statuses,
         "validated": [p for p, s in statuses.items() if s == "validated"]})


def _receipt_rows(root: Path, unit_receipt: dict | None) -> list[dict] | None:
    """受領の単位の試料の対応表の行。`source_file` は pipeline_root 相対なので絶対にして返す。"""
    rows = ((unit_receipt or {}).get("samples") or {}).get("rows")
    if not isinstance(rows, list) or not rows:
        return None
    out = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        row = dict(row)
        source = row.get("source_file")
        if isinstance(source, str) and source and not Path(source).is_absolute():
            row["source_file"] = str(root / source)
        out.append(row)
    return out or None


def load_run_dataset(pipeline_path, *, polarity: str | None = None,
                     allow_incomplete: bool = False):
    """run の成果物から DatasetState を組み立てて返す（保存も session 更新もしない）。

    `pipeline_path` は `pipeline-run.json` でもその親ディレクトリでもよい
    （`pipeline_status` と同じ規約）。v3 の run は選んだ単位のジョブと受領の試料の対応表だけを
    載せる（差次的結果は無い）。旧 v1 の run だけ、差次的結果を指紋（データ・群）が一致する
    ものに限って `ds.results` へ戻し、戻したもの・戻さなかったものを `ds.handoff_report` に残す
    （`restore_differentials`）。
    """
    unit = resolve_unit(pipeline_path, polarity, allow_incomplete=allow_incomplete)
    ds = mztab_loading.load_dataset_state(job_path=unit.job_path,
                                          allow_incomplete=allow_incomplete)

    # 群が無ければ、載せた差次的結果の群の指紋を照合できない。
    if unit.manifest_rows is not None:
        ds.sample_metadata_rows = unit.manifest_rows

    if unit.legacy:
        # 群の指紋は試料対応表を戻した後でなければ run と同じ規則で計算できない。
        ds.handoff_report = restore_differentials(ds, unit.pipeline_root, unit.record)
    return ds


def run_manifest_rows(pipeline_root: Path, record: dict) -> list[dict] | None:
    """run が解決した試料対応表の行を返す。どこにも無ければ None。

    `resolve_metadata` が `record["inputs"]["manifest"]` にだけ書く行
    （`resolve_sample_metadata` の行）を返す。
    これを戻さないと群の指紋が「群なし」になり、差次的結果が全部
    `group_mismatch` で戻らない。`source_file` の絶対/相対の違いは群の指紋
    （`group` と `biological_sample_id` だけを見る）に効かない。
    """
    rows = (record.get("inputs") or {}).get("manifest")
    if isinstance(rows, list) and rows:
        return [dict(row) for row in rows if isinstance(row, dict)]
    return None


def _handoff_rejection(ds, provenance: dict) -> str | None:
    """引き継げない理由を返す（引き継げるなら None）。

    判定の規則そのものは `result_state.is_current` の handoff 分岐と同じ
    （kind はここでは見ない——呼び出し側が `differential:` の ref だけを渡す）。
    ここでは「なぜ戻さなかったか」を利用者に言えるよう、理由を分けて返す。
    """
    recorded_data = provenance.get("dataset_hash")
    if recorded_data is None:
        return "legacy_result_without_dataset_hash"
    if recorded_data != result_state.dataset_fingerprint(ds):
        if not _computed_before_evidence(ds, provenance):
            return "dataset_mismatch"
    recorded_groups = provenance.get("group_fingerprint")
    if recorded_groups is not None and recorded_groups != result_state.group_fingerprint(ds):
        return "group_mismatch"
    return None


def _computed_before_evidence(ds, provenance: dict) -> bool:
    """結果が、検出状態（`.arf` の gap-fill）の取り込み前に計算されたものか。

    取り込みが入ると検出マスクが指紋に入るので、取り込み前に計算された完了済みの
    run の結果は指紋が必ずずれる。検出マスクを除いて計算した今の指紋と一致し、かつ
    今の ds に検出マスクがあるときだけ、取り込みが後から付いただけとみなす
    （行列・特徴・検体など他のどこかが違えば一致しない）。
    """
    if getattr(ds, "detected_mask", None) is None:
        return False
    return provenance.get("dataset_hash") == result_state.dataset_fingerprint(
        ds, include_detected=False)


def _run_preprocess_parameters(pipeline_root: Path, record: dict,
                               parent_ids: list) -> dict | None:
    """差次的結果の親（run の前処理）が実際に使ったレシピを返す。読めなければ None。

    書き出しのメタ行 `# preprocess = ...` は本来 `ds.results` の親結果から引くが、
    run の前処理結果は session へ戻さない（session の前処理ではないので、戻すと
    「今の前処理」に見えてしまう。`is_current` の handoff 分岐も差次的結果以外は
    通さない）。そこでレシピだけを
    引き継ぎ結果の来歴へ写す。写さないと、書き出しは session の
    `preprocessing_recipe`（空、または session で後からやり直した別のレシピ）を
    書き、正しい数字に違う条件のラベルが付く。

    ここは `current_result_refs` ではなく**履歴を result_id で引く**。前処理は
    resume のたびに作り直されて新しい result_id で永続化される一方、成果物で
    検証が通った差次的工程は skip されるので、差次的結果の親は最後の前処理とは
    限らない。前処理の成果物は result_id ごとに別ファイル（`results/<id>.json`）
    なので、古い ref も hash 照合できる。照合できなければ None（不明）を返し、
    読み込み自体は止めない——数字は差次的結果側で検証済みで、欠けるのはラベルだけ。
    """
    wanted = set(parent_ids or [])
    for ref in reversed(record.get("results") or []):
        if not isinstance(ref, dict) or ref.get("output_name") != _PREPROCESS_OUTPUT:
            continue
        if ref.get("result_id") not in wanted:
            continue
        try:
            data = store.read_result_data(pipeline_root, [ref], _PREPROCESS_OUTPUT)[0]
        except (DomainError, OSError, ValueError):
            return None
        # `preprocess_auto` が `_commit_preprocess` に渡したレシピ＝来歴の
        # `effective_parameters` と同じもの（`_handle_preprocess` は provenance を
        # 落として plan を保存するので、ここから読む）。
        recipe = (data or {}).get("resolved_recipe")
        return dict(recipe) if isinstance(recipe, dict) else None
    return None


def restore_differentials(ds, pipeline_root: Path, record: dict) -> dict:
    """run の差次的結果のうち、指紋が今の ds と一致するものを `ds.results` へ戻す。

    読むのは `store.current_result_refs` が選ぶ**今の ref だけ**——同じ比較が
    export で再計算されると `differential:<cid>` が履歴に何本も並び、上書きされた
    古い ref の hash はもうファイルと一致しない（それを「壊れた」と読まない）。
    今の ref の成果物が改変されていれば、行列と同じく `RESULT_INTEGRITY_MISMATCH`
    で読み込みごと止める。

    result_id は付け替えない（pipeline の TSV のメタ行 `# result_id` と同じ ID で
    引けるように）。`ds.last_differential` は触らない——「直近の結果」は session で
    計算したものに限り、戻した結果は result_id で名指しして使う。
    """
    restored: list[dict] = []
    rejected: list[dict] = []
    for ref in store.current_result_refs(record):
        output_name = ref.get("output_name") or ""
        if not output_name.startswith(_DIFFERENTIAL_PREFIX):
            continue
        comparison_id = output_name.split(":", 1)[1]
        result = store.read_result_data(pipeline_root, [ref], output_name)[0] or {}
        provenance = dict(result.get("provenance") or {})
        result_id = provenance.get("result_id") or ref.get("result_id")
        reason = _handoff_rejection(ds, provenance)
        if reason is not None:
            rejected.append({"comparison_id": comparison_id, "result_id": result_id,
                             "reason": reason})
            continue
        provenance["handoff"] = {
            "pipeline_root": str(pipeline_root),
            "comparison_id": comparison_id,
            # None は「run の前処理結果を読めなかった」。書き出しは不明と書く。
            "preprocess": _run_preprocess_parameters(
                pipeline_root, record, provenance.get("parent_ids")),
        }
        # 戻した時点の試料情報・除外の版。以後に変わったら `is_current` が古い結果として扱う。
        provenance["metadata_revision"] = getattr(ds, "metadata_revision", 0)
        if result_state.dataset_fingerprint(ds) != provenance.get("dataset_hash"):
            # 検出状態の取り込み前に計算された結果。`is_current` が同じ規則で通せるよう
            # 来歴に印を付ける（結果の数字は変えない。検出率による絞り込みは、その計算では
            # 掛かっていなかった——利用者へは `restored` の印で知らせる）。
            provenance["evidence_attached_after"] = True
            result["evidence_attached_after"] = True
        result["provenance"] = provenance
        ds.results[result_id] = result
        entry = {"comparison_id": comparison_id, "result_id": result_id,
                 "group_a": result.get("a"), "group_b": result.get("b")}
        if result.get("evidence_attached_after"):
            entry["evidence_attached_after"] = True
        restored.append(entry)
    return {"restored": restored, "rejected": rejected}
