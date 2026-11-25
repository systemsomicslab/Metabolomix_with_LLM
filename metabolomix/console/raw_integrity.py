"""開いた形式の計測ファイル（mzML・imzML・mzXML）が途中で切れていないかの安い検査。

ダウンロードや変換が途中で止まった mzML は、拡張子も名前も正しいまま末尾が欠ける。MS-DIAL
Console はそれを読むところで落ち、どのファイルが原因かをログに残さない（ROSMAP の 4 KB の
mzML。2026-11-13 の事実確認 B7）。そこで Console を起こす前に、ファイルの末尾が根の要素の
閉じタグで終わっているか（`</indexedmzML>` / `</mzML>` / `</mzXML>`）と、サイズが下限
（`MIN_OPEN_FORMAT_BYTES`）以上かを確かめる。読むのは末尾の数 KB（mzML は見出しも）だけで、全体の検証
（XML として読める・スペクトルが壊れていない）はしない。

mzML は見出し（`<run` の前まで）も読み、中身の無い自己終了の `<referenceableParamGroup .../>` など
（`_EMPTY_HEADER_ELEMENTS`）を断る（2026-11-20）。
登録の前に装置の情報を抜いて書き直された mzML（Metabolomics Workbench の ST004581 の 480 本）で、MS-DIAL の読み込みは
これでスペクトルを 1 本も読まず、Console は exit 0 のままピーク 0（`.pai2` が 1 バイト）で終わる。末尾の閉じタグは
正しいので、切れた mzML の検査は通ってしまう。開きタグと閉じタグの形に直せば読める（索引の offset と `fileChecksum`
も作り直す）。

imzML は閉じタグだけを見て、サイズの下限は当てない。imzML が持つのはメタデータだけで、
スペクトルは隣の `.ibd`（バイナリ）にあるので、画素の少ない正しい imzML は小さくてよい。
`.ibd` の中身（途中で切れていないか）は imzML の中の位置と突き合わせないと分からないので見ない。

ベンダーの形式（`.wiff` `.raw` `.d` ほか）は中身の決まりが公開されていないので見ない。
依存は stdlib だけ（pipeline の worker と受付が読む）。
"""
from __future__ import annotations

import re
from pathlib import Path
from xml.parsers import expat

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

#: 見出しを確かめる拡張子と、MS-DIAL が中身の無い自己終了の形で読めなくなる見出しの要素。
#: RawDataHandler 1.3.9725.140（MSDIALCUI の Debug ビルドが `1.3.*-*` で引く版）の mzML の読み込みは、これらの要素で
#: 閉じタグが来るまで読み進むループを開き、`IsEmptyElement` を見ない。自己終了だと閉じタグが来ないのでファイルの
#: 終わりまで読み、`run` もスペクトルも飛ばしてスペクトル 0 で終わる（2026-11-20、IL で確かめた）。見出しのほかの
#: 要素（softwareList・instrumentConfiguration・dataProcessing・sourceFileList など）は自己終了でも害が無い。
_HEADER_CHECKED = frozenset({"mzml"})
_EMPTY_HEADER_ELEMENTS: tuple[str, ...] = (
    "referenceableParamGroupList", "referenceableParamGroup", "fileDescription", "fileContent", "sampleList")
# Expat が認識した開始タグだけを調べる。属性値の中の > や /> はタグの終端ではない。
_SELF_CLOSING = re.compile(r"""<(?:[^"'<>]|"[^"]*"|'[^']*')*/>""")
#: 見出しを読む上限と 1 回に読む量（ふつうの見出しは数十 KB 以下）。
_HEADER_MAX_BYTES = 1024 * 1024
_HEADER_CHUNK = 64 * 1024


def _extension(path: Path) -> str:
    return path.suffix.lower().lstrip(".")


def _empty_header_element(fh) -> str | None:
    """mzML の見出し（`<run` の前まで）にある、中身の無い自己終了の要素の名前（無ければ None）。"""
    fh.seek(0)
    first = fh.read(_HEADER_CHUNK)
    encoding = None
    if first.startswith((b"\xff\xfe", b"<\x00")):
        encoding = "UTF-16LE"
    elif first.startswith((b"\xfe\xff", b"\x00<")):
        encoding = "UTF-16BE"
    parser = expat.ParserCreate(encoding)
    found = None

    class HeaderDone(Exception):
        pass

    def start(name, attrs):
        nonlocal found
        local = name.rsplit(":", 1)[-1]
        if local in ("run", "spectrumList", "chromatogramList"):
            raise HeaderDone
        if local in _EMPTY_HEADER_ELEMENTS:
            # UTF-16 の漢字にも引用符と同じバイトが現れる。開始タグは完全だが、
            # 続く読み取り範囲の末尾は文字の途中で切れることがある。
            tag = parser.GetInputContext().decode(encoding or "latin-1", errors="ignore")
            if _SELF_CLOSING.match(tag):
                found = local
                raise HeaderDone

    parser.StartElementHandler = start
    # DTD の実体を展開しない。目的はファイルに書かれたタグの確認だけ。
    parser.DefaultHandler = lambda _text: None
    consumed = 0
    chunk = first
    try:
        while chunk:
            parser.Parse(chunk, False)
            consumed += len(chunk)
            chunk = fh.read(min(_HEADER_CHUNK, _HEADER_MAX_BYTES - consumed))
        # 長いコメントなどで Expat が保留したトークンも、読み取り済みの範囲で処理する。
        parser.Parse(b"", True)
    except (HeaderDone, expat.ExpatError, LookupError, ValueError):
        pass
    return found


def open_format_problem(path: Path) -> dict | None:
    """開いた形式のファイル 1 つを確かめ、問題があれば `{"name", "size", "reason"}` を返す。

    `reason` は `too_small`（`MIN_OPEN_FORMAT_BYTES` 未満。imzML には当てない）・`truncated`
    （末尾が閉じタグで終わらない）・`empty_header_element`（mzML の見出しに中身の無い自己終了の
    `referenceableParamGroup` など（`_EMPTY_HEADER_ELEMENTS`）がある。`element` に名前）・`unreadable`（読めない）。開いた形式でない・
    問題が無ければ None。
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
            element = _empty_header_element(fh) if extension in _HEADER_CHECKED else None
    except OSError as exc:
        return {"name": path.name, "size": None, "reason": "unreadable", "detail": str(exc)}
    # UTF-16 で書かれていても ASCII の閉じタグは NUL を抜けば同じバイト列になる。
    tail = tail.replace(b"\x00", b"").rstrip()
    if not tail.endswith(trailers):
        return {"name": path.name, "size": size, "reason": "truncated"}
    if size < MIN_OPEN_FORMAT_BYTES and extension not in _NO_SIZE_FLOOR:
        return {"name": path.name, "size": size, "reason": "too_small"}
    if element is not None:
        return {"name": path.name, "size": size, "reason": "empty_header_element", "element": element}
    return None


def broken_open_format_files(paths) -> list[dict]:
    """`paths` のうち開いた形式で問題のあるファイル（`open_format_problem` の結果）の一覧。"""
    found = []
    for path in paths:
        problem = open_format_problem(Path(path))
        if problem is not None:
            found.append(problem)
    return found
