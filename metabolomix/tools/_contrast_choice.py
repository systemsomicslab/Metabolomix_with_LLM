"""対比の選び方の誤りを、ツールの封筒にする（ツールではない共有の部品）。

`stats_results.as_two_group` は対比を省いて読める対比が 2 つ以上あると `CONTRAST_REQUIRED`
（details の `ambiguous: true` と選択肢 `choices`）を送出する。これは誤りではなく利用者に選ばせる問いなので、
対比を 1 つ選ぶツール（`dataset_plot(kind="volcano")`・`save_figure`・`dataset_export_differential`・
`dataset_export_evidence_bundle`）は `needs_input`（`field` / `argument` = `contrast_id`、選択肢は対比の一覧）で
返す。未知の ID・交互作用の名指しなどはそのまま `mztab_error` の封筒にする。
"""
from __future__ import annotations

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error, needs_input

__all__ = ["contrast_envelope"]


def contrast_envelope(exc: DomainError, tool: str) -> str:
    """`DomainError` を封筒にする。対比を省いた曖昧さは `needs_input`、それ以外は `mztab_error`。"""
    details = dict(exc.details or {})
    choices = details.get("choices") or []
    if exc.code == "CONTRAST_REQUIRED" and details.get("ambiguous") and choices:
        options = [{"value": c["id"], "label": c["label"]} for c in choices]
        return needs_input("contrast_id", exc.message, options, argument="contrast_id",
                           retry_tools=[tool])
    return mztab_error(exc.code, exc.message, details or None)
