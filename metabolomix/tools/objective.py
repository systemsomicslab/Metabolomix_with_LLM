"""objective レコード（解析フォルダ配下 reports/<analysis_id>.objective.md）のツール。

`record_objective` が作成と更新の両方を担う（目的が無ければ作り、あれば渡された欄だけを直す）。
deps: mcp_core / session_state / notes（下位レイヤ）。tools_* / server は import しない。
"""
from pathlib import Path

from metabolomix.tools import notes
from metabolomix.core import mcp_core
from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import config_invalid_as_envelope
from metabolomix.tools._holder import resolve_holder

__all__ = [
    "record_objective",
]

#: 目的を新しく作るときに要る引数（既存の目的を直すときは要らない）。
_REQUIRED_FOR_CREATE = ("dataset", "polarity", "groups", "comparison", "sub_questions")


def _resolve_objective_file(analysis_id: str, holder) -> Path | None:
    """レポート候補ディレクトリから objective を探す（ファイル名優先、無ければ frontmatter）。

    objective はレポートと同じ解析フォルダ配下 reports/ に住むので、別の解析フォルダを
    load_dataset した後は見えない（記録はデータの側に付いて回る契約）。`holder` は
    読み書きするハンドル（または対話のセッション）。
    """
    candidates = mcp_core._report_dir_candidates(holder)
    for directory in candidates:
        try:
            direct = notes.note_path(directory, notes.objective_slug(analysis_id))
        except ValueError:
            return None
        if direct.is_file():
            return direct
    for directory in candidates:
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.md")):
            meta, _ = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
            if meta.get("type") == "objective" and str(meta.get("analysis_id", "")) == analysis_id:
                return path
    return None


# 目的を作る呼び出しは同じパスへ上書きするので冪等だが、既存の目的に add_subquestions を渡すと
# 既存 Q の最大+1 で採番して追記するため、ツール全体としては冪等でない。
@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=False), structured_output=False)
@config_invalid_as_envelope
def record_objective(
    analysis_id: str,
    dataset: str | None = None,
    polarity: str | None = None,
    groups: list[str] | None = None,
    comparison: str | None = None,
    sub_questions: list[str] | None = None,
    biological_context: str | None = None,
    inferred_objective: str | None = None,
    confirmed_objective: str | None = None,
    expected_biology: list[str] | None = None,
    assay_kind: str | None = None,
    status: str | None = None,
    add_subquestions: list[str] | None = None,
    dataset_handle: str | None = None,
) -> str:
    """実験目的を解析フォルダ配下 reports/<analysis_id>.objective.md に記録する（作成と更新）。

    その analysis_id の目的がまだ無ければ作る（dataset / polarity / groups / comparison /
    sub_questions が要る。sub_questions は Q1..Qn の本文）。あれば**渡した欄だけ**を直す
    （confirmed_objective・biological_context・status・assay_kind など。創発した小問は
    add_subquestions に渡すと既存 Q の最大+1 で採番して追記する。既存の小問は書き換えない。
    既存の目的に sub_questions を渡しても無視される（追記は add_subquestions を使う））。
    データから推測した目的はユーザーに確認してから confirmed_objective に渡す。
    biological_context は対象系（生物種/細胞/処理）。assay_kind は lipid / metabolite /
    unknown で、確定した種別を渡すと以後のツール出力の解釈規則（脂質名文法か候補集合か）が
    サーバ側で切り替わる（作成時に省けば unknown）。
    dataset_handle: 置き場所（そのデータセットの解析フォルダ配下 reports/）を決めるハンドル
    （省略時は直近。設定の `[paths] reports_dir` があればそちらが先）。どちらも無ければ
    missing_state（load_dataset）。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    if assay_kind is not None:
        try:
            assay_kind = session_state.normalize_assay_kind(assay_kind)
        except ValueError as exc:
            return str(exc)
    given = {"dataset": dataset, "polarity": polarity, "groups": groups,
             "comparison": comparison, "sub_questions": sub_questions}
    # 置き場所を先に決める。設定ファイルが読めない（CONFIG_INVALID）・書けないときに、
    # assay_kind だけが切り替わった半端な状態を残さない。
    try:
        existing = _resolve_objective_file(analysis_id, holder)
        report_dir = mcp_core._resolve_report_dir(holder) if existing is None else None
    except OSError as exc:
        return str(exc)

    if existing is None:
        missing = [name for name in _REQUIRED_FOR_CREATE if given[name] is None]
        if missing:
            return (f"objective が見つかりません: {analysis_id}。新しく作るには "
                    f"{', '.join(missing)} も渡してください。")
        meta_fields = {
            "assay_kind": assay_kind or "unknown",
            "dataset": dataset,
            "polarity": polarity,
            "groups": groups,
            "comparison": comparison,
            "biological_context": biological_context or "",
            "inferred_objective": inferred_objective or "",
            "confirmed_objective": confirmed_objective or "",
            "expected_biology": expected_biology or [],
        }
        if status is not None:
            meta_fields["status"] = status
        try:
            path = notes.write_objective(report_dir, analysis_id, meta_fields, sub_questions)
            if add_subquestions:
                notes.add_subquestions(path, add_subquestions)
        except (ValueError, OSError) as exc:
            return str(exc)
        session_state.session.set_assay_kind(assay_kind or "unknown")
        total = len(sub_questions) + len(add_subquestions or [])
        return (f"objective を記録: {path.name}（confirmed={bool(confirmed_objective)}, "
                f"小問{total}件）。")

    updates = {name: value for name, value in (
        ("dataset", dataset), ("polarity", polarity), ("groups", groups),
        ("comparison", comparison), ("biological_context", biological_context),
        ("inferred_objective", inferred_objective),
        ("confirmed_objective", confirmed_objective),
        ("expected_biology", expected_biology), ("status", status)) if value is not None}
    if assay_kind is not None:
        updates["assay_kind"] = assay_kind
    try:
        if updates:
            notes.update_objective_meta(existing, updates)
        if add_subquestions:
            notes.add_subquestions(existing, add_subquestions)
    except (ValueError, OSError) as exc:
        return str(exc)
    if assay_kind is not None:
        session_state.session.set_assay_kind(assay_kind)
    note = ("。sub_questions は既存の目的には使わない（追記は add_subquestions）"
            if sub_questions is not None else "")
    return (f"objective を更新: {existing.name}（更新フィールド={list(updates) or 'なし'}, "
            f"追加小問={len(add_subquestions or [])}件）{note}")
