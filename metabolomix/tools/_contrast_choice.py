"""対比の選び方の誤りを、ツールの封筒にする（ツールではない共有の部品）。

`stats_results.as_two_group` は対比を省いて読める対比が 2 つ以上あると `CONTRAST_REQUIRED`
（details の `ambiguous: true` と選択肢 `choices`）を送出する。これは誤りではなく利用者に選ばせる問いなので、
対比を 1 つ選ぶツール（`dataset_plot(kind="volcano")`・`save_figure`・`dataset_export_differential`・
`dataset_export_evidence_bundle`）は `needs_input`（`field` / `argument` = `contrast_id`、選択肢は対比の一覧）で
返す。spec で対比を選ぶツール（`dataset_novel_candidates` の `statistics.result_id`）は `spec_field` を渡し、
`field` をその spec の路に、選択肢の値を `result_id#contrast_id` の形にする（`argument` は載せない）。
未知の ID・交互作用の名指しなどはそのまま `mztab_error` の封筒にする。
"""
from __future__ import annotations

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error, needs_input

__all__ = ["contrast_envelope"]


def contrast_envelope(exc: DomainError, tool: str, *, spec_field: str | None = None) -> str:
    """`DomainError` を封筒にする。対比を省いた曖昧さは `needs_input`、それ以外は `mztab_error`。

    `spec_field`（例 ``"statistics.result_id"``）を渡すと、答えはその spec の欄に `result_id#contrast_id` の形で
    入れて呼び直す問いになる（`argument` を持たないツール向け）。
    """
    details = dict(exc.details or {})
    choices = details.get("choices") or []
    if (exc.code == "CONTRAST_REQUIRED" and details.get("ambiguous") and choices
            and (spec_field is None or details.get("result_id"))):
        if spec_field is None:
            options = [{"value": c["id"], "label": c["label"]} for c in choices]
            return needs_input("contrast_id", exc.message, options, argument="contrast_id",
                               retry_tools=[tool])
        rid = details.get("result_id")
        options = [{"value": f"{rid}#{c['id']}", "label": c["label"]} for c in choices]
        reason = (f"この比較（result_id={rid}）には対比が {len(choices)} つあります。spec の {spec_field} を"
                  " `result_id#contrast_id` の形にして、どの対比を使うか選んでください。")
        return needs_input(spec_field, reason, options, retry_tools=[tool])
    return mztab_error(exc.code, exc.message, details or None)
