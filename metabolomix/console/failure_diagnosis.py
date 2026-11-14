"""MS-DIAL Console が完了しなかったときの原因の見立て（ログと入力から）。

Console は LC-MS の処理でファイルごとの名前をログに出さないことが多く、落ちても「どの計測
ファイルで」を自分では言わない（2026-11-13 の事実確認 B7）。ここでは手掛かりを 3 つ集め、
`MSDIAL_EXECUTION_FAILED` の details に `diagnosis` として返す:

1. **資源不足** — 終了コード 1450（Windows の `ERROR_NO_SYSTEM_RESOURCES`。commit の上限）か、
   ログに同じ意味の文面・`OutOfMemoryException` があれば `cause: "system_resources"`。
2. **壊れた入力** — 配置済みの開いた形式（mzML ほか）を `raw_integrity` でもう一度確かめる
   （受付の後に切れたものも拾う）。見つかれば `cause: "broken_input"`。
3. **ログに名前の出た入力** — エラーの行に計測ファイルの名前があればそれ、無ければ最初の
   エラーより前に最後に名前の出たファイルを「疑わしい」として挙げる（見立てで、確定ではない）。

ログの行を返すときは、絶対パスをファイル名にする（研究室の参照ライブラリの置き場所を戻り値に
出さない。CLAUDE.md）。依存は stdlib と同じパッケージの `raw_integrity` だけ。
"""
from __future__ import annotations

import re
from pathlib import Path, PureWindowsPath

from metabolomix.console.raw_integrity import broken_open_format_files

__all__ = ["SYSTEM_RESOURCES_EXIT_CODES", "diagnose_failure"]

#: 資源不足と読む終了コード。1450 = ERROR_NO_SYSTEM_RESOURCES（commit の上限に当たった）。
SYSTEM_RESOURCES_EXIT_CODES = frozenset({1450})

_RESOURCE_PATTERNS = re.compile(
    r"ERROR_NO_SYSTEM_RESOURCES|Insufficient system resources|OutOfMemoryException"
    r"|システム リソースが不足|システムリソースが不足", re.IGNORECASE)
_ERROR_PATTERN = re.compile(r"exception|\berror\b|\bfail(?:ed|ure)?\b|エラー|例外|失敗",
                            re.IGNORECASE)
#: パスらしい語（バックスラッシュを含む語と、`C:/` で始まる語）。ファイル名だけに縮める。
_PATH_TOKEN = re.compile(r"[^\s\"'<>|]*\\[^\s\"'<>|]*|[A-Za-z]:/[^\s\"'<>|]*")

#: ログを読む量（末尾から。落ちた理由は終わりの近くにある）。
_LOG_TAIL_BYTES = 512 * 1024
_MAX_LINES = 3
_MAX_LINE_CHARS = 240
_MAX_SUSPECTS = 10


def _read_log_tail(log_path: Path) -> list[str]:
    try:
        with log_path.open("rb") as fh:
            fh.seek(0, 2)
            size = fh.tell()
            start = max(0, size - _LOG_TAIL_BYTES)
            fh.seek(start)
            data = fh.read()
    except OSError:
        return []
    if start > 0:
        # 途中から読んだ最初の行は文字の途中で切れうる（cp932 の 2 バイト目など）ので捨てる。
        data = data[data.find(b"\n") + 1:] if b"\n" in data else b""
    # .NET の Console は出力をリダイレクトすると OEM のコードページ（日本語の Windows は cp932）で書く。
    for encoding in ("utf-8", "cp932"):
        try:
            text = data.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    else:
        text = data.decode("utf-8", errors="replace")
    return text.splitlines()


def _file_name(token: str) -> str:
    return PureWindowsPath(token.rstrip("\\/.,;:)")).name or "…"


def _file_names_only(line: str) -> str:
    """行の中のパスらしい語をファイル名にし、長すぎる行を切る。"""
    line = _PATH_TOKEN.sub(lambda m: _file_name(m.group(0)), line).strip()
    return line if len(line) <= _MAX_LINE_CHARS else line[:_MAX_LINE_CHARS] + "…"


def _name_patterns(names) -> dict[str, re.Pattern]:
    """計測ファイルの名前ごとの検索（`S1.mzML` が `AS1.mzML` や `S11.mzML` に当たらない）。"""
    return {name: re.compile(r"(?<![\w.-])" + re.escape(name) + r"(?!\w)", re.IGNORECASE)
            for name in names}


def _names_in(line: str, patterns: dict[str, re.Pattern]) -> list[str]:
    return [name for name, pattern in patterns.items() if pattern.search(line)]


def diagnose_failure(log_path: Path, *, exit_code: int | None, raw_paths) -> dict:
    """完了しなかった Console の手掛かりを返す。何も見つからなければ空の dict。

    戻り値のキー（見つかったものだけ）: `cause`（`system_resources` / `broken_input` /
    `input_suspected`）・`message`・`suspect_files`（`[{"name", "reason"}]`。`reason` は
    `truncated` / `too_small` / `unreadable` / `named_in_error` / `last_named_before_error`）・
    `log_lines`（エラーの行。絶対パスはファイル名にしてある）。
    """
    raw_paths = [Path(p) for p in raw_paths or []]
    names = _name_patterns(sorted({p.name for p in raw_paths}))
    lines = _read_log_tail(Path(log_path))

    resource_hit = (exit_code in SYSTEM_RESOURCES_EXIT_CODES
                    or any(_RESOURCE_PATTERNS.search(line) for line in lines))
    error_indexes = [i for i, line in enumerate(lines) if _ERROR_PATTERN.search(line)]

    suspects: list[dict] = []
    seen: set[str] = set()

    def _add(name: str, reason: str) -> None:
        if name not in seen and len(suspects) < _MAX_SUSPECTS:
            seen.add(name)
            suspects.append({"name": name, "reason": reason})

    for problem in broken_open_format_files(raw_paths):
        _add(problem["name"], problem["reason"])
    broken = bool(suspects)
    for index in error_indexes:
        for name in _names_in(lines[index], names):
            _add(name, "named_in_error")
    if error_indexes and len(suspects) == 0:
        for line in reversed(lines[:error_indexes[0]]):
            found = _names_in(line, names)
            if found:
                _add(found[0], "last_named_before_error")
                break

    result: dict = {}
    if resource_hit:
        result["cause"] = "system_resources"
        result["message"] = (
            "資源不足で落ちた（終了コード 1450 = ERROR_NO_SYSTEM_RESOURCES か、ログのメモリ不足）。"
            "commit（物理メモリ＋ページファイル）の上限に当たった見込みが高い。ほかの重いプロセス"
            "（同時に走る別の Console を含む）を止めるかページファイルを広げてから、"
            "pipeline_resume(rerun_upstream=true) でやり直す。")
    elif broken:
        result["cause"] = "broken_input"
        result["message"] = ("途中で切れた・小さすぎる計測ファイルがある（suspect_files）。取り直すか"
                             "フォルダの外へ移してから、pipeline_run で受け付け直す。")
    elif suspects:
        result["cause"] = "input_suspected"
        result["message"] = ("ログのエラーの近くに名前の出た計測ファイルがある（suspect_files。見立てで"
                             "確定ではない）。そのファイルを MS-DIAL の GUI で開けるか確かめる。")
    if suspects:
        result["suspect_files"] = suspects
    error_lines = [_file_names_only(lines[i]) for i in error_indexes[-_MAX_LINES:]]
    if error_lines:
        result["log_lines"] = error_lines
    return result
