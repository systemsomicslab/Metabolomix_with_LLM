"""ビューア HTML の組み立て。テンプレートは同じフォルダの viewer.html。

データを埋め込むときは `</` を `<\\/` にして script 要素から抜け出せないようにする。
None を渡すと MCP Apps 用の空テンプレート（データは curation_view_data で取る）。
"""
from __future__ import annotations

import json
from pathlib import Path

from lipidmix.curation.flags import SUBMISSION_PREFIX

_TEMPLATE = Path(__file__).with_name("viewer.html")
_PLACEHOLDER = "/*__CURATION_DATA__*/null"


def render_html(review: dict | None) -> str:
    template = _TEMPLATE.read_text(encoding="utf-8")
    if review is None:
        data = "null"
    else:
        data = json.dumps(review, ensure_ascii=False, separators=(",", ":"), allow_nan=False,
                          default=str).replace("</", "<\\/")
    html = template.replace(_PLACEHOLDER, data, 1)
    return html.replace("__SUBMISSION_PREFIX__", SUBMISSION_PREFIX.strip())
