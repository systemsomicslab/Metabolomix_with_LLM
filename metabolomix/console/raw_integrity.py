"""開いた形式の計測ファイル（mzML・imzML・mzXML）が途中で切れていないかの安い検査。

ダウンロードや変換が途中で止まった mzML は、拡張子も名前も正しいまま末尾が欠ける。MS-DIAL
Console はそれを読むところで落ち、どのファイルが原因かをログに残さない（ROSMAP の 4 KB の
mzML。2026-11-13 の事実確認 B7）。そこで Console を起こす前に、ファイルの末尾が根の要素の
閉じタグで終わっているか（`</indexedmzML>` / `</mzML>` / `</mzXML>`）と、サイズが下限
（`MIN_OPEN_FORMAT_BYTES`）以上かを確かめる。読むのは末尾の数 KB だけで、全体の検証
（XML として読める・スペクトルが壊れていない）はしない。

imzML は閉じタグだけを見て、サイズの下限は当てない。imzML が持つのはメタデータだけで、
スペクトルは隣の `.ibd`（バイナリ）にあるので、画素の少ない正しい imzML は小さくてよい。
`.ibd` の中身（途中で切れていないか）は imzML の中の位置と突き合わせないと分からないので見ない。

ベンダーの形式（`.wiff` `.raw` `.d` ほか）は中身の決まりが公開されていないので見ない。
依存は stdlib だけ（pipeline の worker と受付が読む）。
"""
from __future__ import annotations

from pathlib import Path

__all__ = [
    "MIN_OPEN_FORMAT_BYTES",
    "OPEN_FORMAT_EXTENSIONS",
    "broken_open_format_files",
    "open_format_problem",
]

#: 末尾の閉じタグ（拡張子 → 受ける閉じタグ）。imzML は mzML の上に乗る形式で、根は `mzML`。
_TRAILERS: dict[str, tuple[bytes, ...]] = {
    "mzml": (b"</indexedmzML>", b"</mzML>"),
    "imzml": (b"</indexedmzML>", b"</mzML>"),
    "mzxml": (b"</mzXML>",),
}

OPEN_FORMAT_EXTENSIONS = frozenset(_TRAILERS)

#: これより小さい開いた形式のファイルは壊れているとみなす（バイト）。LC-MS の 1 注入の
#: mzML は数 MB 以上になる。ヘッダだけ・スペクトル 0 件のファイル（閉じタグはある）を捕まえる。
#: スペクトルを中に持たない imzML（`_NO_SIZE_FLOOR`）には当てない。
MIN_OPEN_FORMAT_BYTES = 32 * 1024

#: サイズの下限を当てない拡張子。imzML のスペクトルは隣の `.ibd` にある。
_NO_SIZE_FLOOR = frozenset({"imzml"})

#: 末尾を読む量。閉じタグの後ろの空白・改行と、UTF-16 の 2 倍の幅を見込む。
_TAIL_BYTES = 4096


def _extension(path: Path) -> str:
    return path.suffix.lower().lstrip(".")


def open_format_problem(path: Path) -> dict | None:
    """開いた形式のファイル 1 つを確かめ、問題があれば `{"name", "size", "reason"}` を返す。

    `reason` は `too_small`（`MIN_OPEN_FORMAT_BYTES` 未満。imzML には当てない）・`truncated`
    （末尾が閉じタグで終わらない）・`unreadable`（読めない）。開いた形式でない・問題が無ければ None。
    """
    path = Path(path)
    extension = _extension(path)
    trailers = _TRAILERS.get(extension)
    if trailers is None or not path.is_file():
        return None
    try:
        size = path.stat().st_size
        with path.open("rb") as fh:
            fh.seek(max(0, size - _TAIL_BYTES))
            tail = fh.read()
    except OSError as exc:
        return {"name": path.name, "size": None, "reason": "unreadable", "detail": str(exc)}
    # UTF-16 で書かれていても ASCII の閉じタグは NUL を抜けば同じバイト列になる。
    tail = tail.replace(b"\x00", b"").rstrip()
    if not tail.endswith(trailers):
        return {"name": path.name, "size": size, "reason": "truncated"}
    if size < MIN_OPEN_FORMAT_BYTES and extension not in _NO_SIZE_FLOOR:
        return {"name": path.name, "size": size, "reason": "too_small"}
    return None


def broken_open_format_files(paths) -> list[dict]:
    """`paths` のうち開いた形式で問題のあるファイル（`open_format_problem` の結果）の一覧。"""
    found = []
    for path in paths:
        problem = open_format_problem(Path(path))
        if problem is not None:
            found.append(problem)
    return found
