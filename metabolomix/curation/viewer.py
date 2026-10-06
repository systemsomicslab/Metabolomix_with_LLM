"""ビューア HTML の組み立て。テンプレートは同じフォルダの viewer.html。

置換の順序: 先に `__SUBMISSION_PREFIX__` をテンプレートへ入れ、その**後で**データを埋め込む
（逆にすると、スポット名などに同じ文字列があったとき埋め込んだデータまで書き換わる）。
埋め込む JSON の `<` はすべて `\\u003c` にする——`</script>` だけでなく `<!--` などでも
script 要素の解析状態が変わりうるため。`<` は JSON では文字列の中にしか現れないので、
`\\u003c` は同じ文字として読み戻される。
None を渡すと MCP Apps 用の空テンプレート（データは curation_view_data で取る）。
"""
from __future__ import annotations

import json
from pathlib import Path

from metabolomix.curation.flags import SUBMISSION_PREFIX

_TEMPLATE = Path(__file__).with_name("viewer.html")
_SUGGEST_TEMPLATE = Path(__file__).with_name("suggest_viewer.html")
_COMMON_JS = Path(__file__).with_name("viewer_common.js")
_PLACEHOLDER = "/*__CURATION_DATA__*/null"
_COMMON_PLACEHOLDER = "/*__VIEWER_COMMON__*/"


def _embed(data) -> str:
    if data is None:
        return "null"
    return json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False,
                      default=str).replace("<", "\\u003c")


def _render(template_path: Path, data) -> str:
    template = template_path.read_text(encoding="utf-8")
    template = template.replace(_COMMON_PLACEHOLDER, _COMMON_JS.read_text(encoding="utf-8"), 1)
    template = template.replace("__SUBMISSION_PREFIX__", SUBMISSION_PREFIX.strip())
    return template.replace(_PLACEHOLDER, _embed(data), 1)


def render_html(review: dict | None) -> str:
    return _render(_TEMPLATE, review)


def render_suggest_html(suggestion: dict | None) -> str:
    return _render(_SUGGEST_TEMPLATE, suggestion)
