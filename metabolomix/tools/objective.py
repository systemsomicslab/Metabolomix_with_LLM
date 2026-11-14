"""objective レコード（解析フォルダ配下 reports/<analysis_id>.objective.md）のツール。

`record_objective` が作成と更新の両方を担う（目的が無ければ渡した欄だけで作り、あれば渡された欄だけを直す）。
アッセイ種別（`assay_kind`）はハンドル（データセット）ごとに持ち、目的ファイルにも残す。`load_dataset` は
`recall_assay_kind` で解析フォルダの目的ファイルから読み戻す。
deps: mcp_core / session_state / notes（下位レイヤ）。tools_* / server は import しない。
"""
import os
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

#: 目的を新しく作るときの欄の既定（渡さなかった欄は空のまま作り、後から `record_objective` で足せる）。
_CREATE_DEFAULTS = {"dataset": "", "polarity": "", "groups": [], "comparison": ""}


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


def _with_caveat(text: str, assay_kind: str | None) -> str:
    """種別を渡した呼び出しの戻り値の先頭へ、確定した種別の意味論ダイジェストを（1 回だけ）前置する。

    入口の `load_dataset` はたいてい種別が `unknown` のうちに走るので、確定後の規則
    （脂質名文法か候補集合か）はここで届く。同じ種別のままなら何も付かない
    （`session.take_caveat` の規則: リソース既読・`caveat_mode=off` でも付かない）。
    """
    if assay_kind is None:
        return text
    caveat = session_state.session.take_caveat(assay_kind)
    return f"{caveat}\n\n{text}" if caveat else text


def _analysis_dir(holder) -> str | None:
    return getattr(getattr(holder, "dataset", None), "analysis_dir", None) or None


def _same_dir(a, b) -> bool:
    try:
        return (os.path.normcase(str(Path(a).resolve()))
                == os.path.normcase(str(Path(b).resolve())))
    except (OSError, ValueError, TypeError):
        return False


def recall_assay_kind(ds, report_dirs) -> tuple[str | None, str | None]:
    """解析フォルダの目的ファイル（`reports/*.objective.md`）に記録した種別を読み戻す。(種別, 注記)。

    読むのは `report_dirs`（記録の置き場所の候補。`mcp_core._report_dir_candidates`）の目的ファイルの
    うち、そのデータセットの解析フォルダ配下 `reports/` にあるもの全部と、ほかの置き場所（設定の
    `[paths] reports_dir`）にあって frontmatter の `analysis_dir` がこの解析フォルダのもの。
    `unknown` でない種別が 1 つに決まればそれ、食い違えば None と理由（どれも採らない）、
    無ければ (None, None)。読めないファイルは飛ばす。
    """
    analysis_dir = getattr(ds, "analysis_dir", None)
    if not analysis_dir:
        return None, None
    own = Path(analysis_dir) / "reports"
    found: dict[str, list[str]] = {}
    seen: set[str] = set()
    for directory in report_dirs:
        directory = Path(directory)
        if not directory.is_dir():
            continue
        mine = _same_dir(directory, own)
        for path in sorted(directory.glob(f"*{notes.OBJECTIVE_SUFFIX}.md")):
            key = os.path.normcase(str(path.resolve()))
            if key in seen:
                continue
            seen.add(key)
            try:
                meta, _ = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError):
                continue
            if meta.get("type") != "objective":
                continue
            if not mine and not (meta.get("analysis_dir")
                                 and _same_dir(meta["analysis_dir"], analysis_dir)):
                continue
            kind = str(meta.get("assay_kind") or "").strip().lower()
            if kind in ("lipid", "metabolite"):
                found.setdefault(kind, []).append(path.name)
    if not found:
        return None, None
    if len(found) > 1:
        listed = "、".join(f"{kind}（{', '.join(names)}）" for kind, names in sorted(found.items()))
        return None, (f"目的ファイルの assay_kind が食い違う（{listed}）ので読み戻していない"
                      "（unknown のまま。record_objective(assay_kind=...) で確定する）。")
    kind, names = next(iter(found.items()))
    return kind, f"assay_kind={kind} を目的ファイル（{', '.join(names)}）から読み戻した。"


def _kind_note(holder, assay_kind: str | None) -> str:
    """種別を書いたときの戻り値の 1 句（どのハンドルの種別が切り替わったか）。"""
    if assay_kind is None:
        return ""
    handle_id = getattr(holder, "handle_id", None) or getattr(holder, "current_handle_id", None)
    where = f"ハンドル {handle_id}" if handle_id else "セッション"
    return f", assay_kind={assay_kind}（{where} のデータセットに記録）"


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

    その analysis_id の目的がまだ無ければ、**渡した欄だけで**作る（欠けた欄は空のまま。後から
    この関数で足せる）。sub_questions は Q1..Qn の本文。あれば**渡した欄だけ**を直す
    （confirmed_objective・biological_context・status・assay_kind など。創発した小問は
    add_subquestions に渡すと既存 Q の最大+1 で採番して追記する。既存の小問は書き換えない。
    小問の無い目的に sub_questions を渡すと Q1..Qn として足し、小問のある目的に渡しても
    無視される（追記は add_subquestions を使う））。
    biological_context は対象系（生物種/細胞/処理）。assay_kind は lipid / metabolite /
    unknown で、**ハンドル（データセット）ごとに持つ**: 渡すと書いたハンドルの種別が切り替わり、
    以後そのハンドルの解釈規則（ダイジェストの脂質名文法か候補集合か・概観と spec v2 の組成の
    ラベル）が変わる。目的ファイルにも残り、同じ解析フォルダを `load_dataset` で開き直すと
    読み戻す（作成時に省けば unknown で、ハンドルの種別は動かさない）。
    dataset_handle: 置き場所（そのデータセットの解析フォルダ配下 reports/）と種別を書くハンドル
    （省略時は直近。設定の `[paths] reports_dir` があれば置き場所はそちらが先）。どちらも無ければ
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
    # 置き場所を先に決める。設定ファイルが読めない（CONFIG_INVALID）・書けないときに、
    # assay_kind だけが切り替わった半端な状態を残さない。
    try:
        existing = _resolve_objective_file(analysis_id, holder)
        report_dir = mcp_core._resolve_report_dir(holder) if existing is None else None
    except OSError as exc:
        return str(exc)

    if existing is None:
        meta_fields = {
            "assay_kind": assay_kind or "unknown",
            "dataset": dataset if dataset is not None else _CREATE_DEFAULTS["dataset"],
            "polarity": polarity if polarity is not None else _CREATE_DEFAULTS["polarity"],
            "groups": groups if groups is not None else list(_CREATE_DEFAULTS["groups"]),
            "comparison": comparison if comparison is not None else _CREATE_DEFAULTS["comparison"],
            "biological_context": biological_context or "",
            "inferred_objective": inferred_objective or "",
            "confirmed_objective": confirmed_objective or "",
            "expected_biology": expected_biology or [],
        }
        if _analysis_dir(holder):
            # 読み戻し（`recall_assay_kind`）が、共有の置き場所（`[paths] reports_dir`）でも
            # どの解析フォルダの目的かを見分ける。
            meta_fields["analysis_dir"] = _analysis_dir(holder)
        if status is not None:
            meta_fields["status"] = status
        try:
            path = notes.write_objective(report_dir, analysis_id, meta_fields, sub_questions or [])
            if add_subquestions:
                notes.add_subquestions(path, add_subquestions)
        except (ValueError, OSError) as exc:
            return str(exc)
        if assay_kind is not None:   # 渡さない作成は、ハンドルの種別を動かさない
            holder.assay_kind = assay_kind
        total = len(sub_questions or []) + len(add_subquestions or [])
        empty = [name for name, value in (("dataset", dataset), ("polarity", polarity),
                                          ("groups", groups), ("comparison", comparison),
                                          ("sub_questions", sub_questions)) if value is None]
        rest = (f"。空のまま作った欄: {', '.join(empty)}（後から record_objective で足せる）"
                if empty else "")
        return _with_caveat(f"objective を記録: {path.name}（confirmed={bool(confirmed_objective)}, "
                            f"小問{total}件{_kind_note(holder, assay_kind)}）{rest}。", assay_kind)

    updates = {name: value for name, value in (
        ("dataset", dataset), ("polarity", polarity), ("groups", groups),
        ("comparison", comparison), ("biological_context", biological_context),
        ("inferred_objective", inferred_objective),
        ("confirmed_objective", confirmed_objective),
        ("expected_biology", expected_biology), ("status", status)) if value is not None}
    if assay_kind is not None:
        updates["assay_kind"] = assay_kind
    first_subquestions = False
    try:
        meta, subqs, _body = notes.parse_objective(existing)
        if updates and not meta.get("analysis_dir") and _analysis_dir(holder):
            updates["analysis_dir"] = _analysis_dir(holder)
        if updates:
            notes.update_objective_meta(existing, updates)
        # 小問の無い目的（欄を省いて作ったもの）にだけ、sub_questions を Q1..Qn として足す。
        first_subquestions = bool(sub_questions) and not subqs
        if first_subquestions:
            notes.add_subquestions(existing, sub_questions, initial=True)
        if add_subquestions:
            notes.add_subquestions(existing, add_subquestions)
    except (ValueError, OSError) as exc:
        return str(exc)
    if assay_kind is not None:
        holder.assay_kind = assay_kind
    note = ("。sub_questions は小問のある目的には使わない（追記は add_subquestions）"
            if sub_questions is not None and not first_subquestions else "")
    added = len(add_subquestions or []) + (len(sub_questions) if first_subquestions else 0)
    fields = [*[k for k in updates if k != "analysis_dir"],
              *(["sub_questions"] if first_subquestions else [])]
    return _with_caveat(f"objective を更新: {existing.name}（更新フィールド={fields or 'なし'}, "
                        f"追加小問={added}件{_kind_note(holder, assay_kind)}）{note}", assay_kind)
