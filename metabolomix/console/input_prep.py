"""測定ファイルの随伴ファイルの規則（pipeline の入力配置が使う）。

`.wiff` は `.wiff.scan` が無いと読めない。MS-DIAL の `SupportMsRawDataExtension` は
`wiff` と `wiff2` を別フォーマットとして数えるので、随伴は明示規則の完全一致だけを返す。
"""
from __future__ import annotations

from pathlib import Path


#: 随伴ファイルの明示規則（形式ごと）。要素は (基準, 追加suffix)。
#: 基準 "name" は primary のファイル名全体に、"stem" は主拡張子を除いた部分に
#: suffix を足した名前を候補にし、その名前がフォルダに実在するときだけ随伴として扱う。
#:
#: 以前は「primary の先頭トークンで始まり、計測拡張子でない」という前方一致で
#: 判定していたが、これは同じ stem を持つ PAI2/DCL/ARF/タグファイル
#: （例: `sample_a.pai2`）まで随伴として誤って拾ってしまう欠陥だった。
#: 明示登録した完全一致の名前だけを随伴とすることでこれを防ぐ。
COMPANION_RULES: dict[str, tuple[tuple[str, str], ...]] = {
    "wiff": (("name", ".scan"), ("stem", ".timeseries.data")),
}


def _companions_of(entries: list[Path], primary: Path, ext: str) -> list[Path]:
    """primary と同じ測定に属する随伴ファイルを、`COMPANION_RULES` の完全一致でだけ返す。

    `.wiff` は `.wiff.scan` が無いと読めず、`.timeseries.data` も同じ測定の一部。
    一方 `.wiff2` は**別フォーマット**なので連れて行ってはいけない —— それを
    置いていくことがこの関数の目的そのものだから（`COMPANION_RULES` に無い名前は
    そもそも候補に上がらない）。
    """
    by_name = {p.name: p for p in entries}
    out = []
    for base, suffix in COMPANION_RULES.get(ext, ()):
        name = (primary.name if base == "name" else primary.stem) + suffix
        found = by_name.get(name)
        if found is not None and found != primary:
            out.append(found)
    return out
