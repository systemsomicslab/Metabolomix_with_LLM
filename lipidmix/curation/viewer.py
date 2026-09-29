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

from lipidmix.curation.flags import SUBMISSION_PREFIX

_TEMPLATE = Path(__file__).with_name("viewer.html")
_PLACEHOLDER = "/*__CURATION_DATA__*/null"


def render_html(review: dict | None) -> str:
    template = _TEMPLATE.read_text(encoding="utf-8")
    template = template.replace("__SUBMISSION_PREFIX__", SUBMISSION_PREFIX.strip())
    if review is None:
        data = "null"
    else:
        data = json.dumps(review, ensure_ascii=False, separators=(",", ":"), allow_nan=False,
                          default=str).replace("<", "\\u003c")
    return template.replace(_PLACEHOLDER, data, 1)


def render_suggest_html(suggestion: dict | None) -> str:
    return "<!doctype html><title>suggest</title>"   # Task 7 で本実装に置き換える
