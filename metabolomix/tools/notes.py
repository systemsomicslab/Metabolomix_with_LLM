"""レポートと objective の読み書き（frontmatter 付き Markdown。stdlib だけ）。

MCP にも session にも依存しない純ロジック。`tools/reports.py` と `tools/objective.py` が
解析フォルダ配下 reports/ へ書くときに使う。

- frontmatter は `---` で囲んだ最小の YAML（スカラ・インラインリスト・1 段ネスト）。
  書き出し（`dump_frontmatter`）と読み取り（`parse_frontmatter`）はこのモジュールの中で閉じ、
  往復で値が変わらないことをテストで縛る。
- 書き込みは置き場ごとのロックファイルで直列にする（`write_note`）。
- 名前（`analysis_id`）は呼び出し側（LLM）由来なので、置き場の外へ出るものは `note_path`
  で弾く。
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import os
import re
import time


# --- frontmatter の読み取り ---

_UNESCAPES = {'"': '"', "\\": "\\", "n": "\n", "r": "\r"}


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        # dump_frontmatter が \\ \" 改行（\n）CR（\r）をエスケープする。1 回の走査で戻す
        # （順に置換すると、リテラルの「\」+「n」が改行に化ける）。
        return re.sub(r'\\(["\\nr])', lambda m: _UNESCAPES[m.group(1)], value[1:-1])
    if len(value) >= 2 and value[0] == value[-1] == "'":
        return value[1:-1]
    return value


def _split_inline(inner: str) -> list[str]:
    """インラインリストの中身を、引用符の中のカンマを区切りにせず分ける。"""
    items: list[str] = []
    current: list[str] = []
    quote: str | None = None
    escaped = False
    for ch in inner:
        if quote == '"' and escaped:
            escaped = False
        elif quote == '"' and ch == "\\":
            escaped = True
        elif quote is not None:
            if ch == quote:
                quote = None
        elif ch in "\"'" and not "".join(current).strip():
            # 引用符は項目の先頭でだけ開く（Parkinson's の ' は項目の中身）。
            quote = ch
        elif ch == ",":
            items.append("".join(current))
            current = []
            continue
        current.append(ch)
    items.append("".join(current))
    return [item for item in items if item.strip()]


def _scalar(value: str):
    value = value.strip()
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        return [_unquote(item) for item in _split_inline(inner)] if inner else []
    if value in ("True", "False"):
        return value == "True"
    return _unquote(value)


def _parse_meta(text: str) -> dict:
    """frontmatter の最小サブセット（スカラ・インラインリスト・1 段ネスト）を解釈する。"""
    result: dict = {}
    current_key: str | None = None
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip())
        line = raw.strip()
        if indent > 0 and current_key is not None and isinstance(result.get(current_key), dict):
            key, sep, val = line.partition(":")
            if sep:
                result[current_key][key.strip()] = _scalar(val)
            continue
        key, sep, val = line.partition(":")
        if not sep:
            continue
        key = key.strip()
        val = val.strip()
        if val == "":
            result[key] = {}
            current_key = key
        else:
            result[key] = _scalar(val)
            current_key = None
    return result


def parse_frontmatter(text: str) -> tuple[dict, str]:
    """``---`` で囲まれた frontmatter と本文に分割する。"""
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return {}, text
    end = None
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            end = i
            break
    if end is None:
        return {}, text
    meta = _parse_meta("\n".join(lines[1:end]))
    body = "\n".join(lines[end + 1:]).strip()
    return meta, body


# --- frontmatter の書き出し ---

def _dump_scalar(value) -> str:
    if isinstance(value, list):
        return "[" + ", ".join(_dump_scalar(v) for v in value) + "]"
    text = str(value)
    needs_quote = (
        text == ""
        or text.strip() != text
        or any(ch in text for ch in ":#[]{}\",'")
        or "\n" in text
        or "\r" in text
    )
    if needs_quote:
        escaped = (text.replace("\\", "\\\\").replace('"', '\\"')
                   .replace("\n", "\\n").replace("\r", "\\r"))
        return '"' + escaped + '"'
    return text


def dump_frontmatter(meta: dict) -> str:
    """dict を ``---`` 囲みの frontmatter 文字列にする（1 段ネストまで対応）。"""
    lines = ["---"]
    for key, value in meta.items():
        if value is None:  # 「None」という文字列として書かない（欄ごと書かない）
            continue
        if isinstance(value, dict):
            lines.append(f"{key}:")
            for inner_key, inner_value in value.items():
                if inner_value is not None:
                    lines.append(f"  {inner_key}: {_dump_scalar(inner_value)}")
        else:
            lines.append(f"{key}: {_dump_scalar(value)}")
    lines.append("---")
    return "\n".join(lines)


# --- 置き場のロックとノートの書き込み ---

LOCK_FILENAME = ".metabolomix.lock"
LOCK_TIMEOUT_SECONDS = 10.0
LOCK_POLL_SECONDS = 0.1


@contextmanager
def _directory_lock(directory: str | Path, timeout: float = LOCK_TIMEOUT_SECONDS):
    """置き場ごとの書き込みを、ロックファイル（原子的な作成）で直列にする。"""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    lock_path = directory / LOCK_FILENAME
    deadline = time.monotonic() + timeout
    fd: int | None = None
    while fd is None:
        try:
            fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, f"pid={os.getpid()}\n".encode("utf-8"))
        except FileExistsError:
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Timed out waiting for note lock: {lock_path}")
            time.sleep(LOCK_POLL_SECONDS)
    try:
        yield
    finally:
        if fd is not None:
            os.close(fd)
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass


def note_path(directory: str | Path, slug: str) -> Path:
    """``<directory>/<slug>.md`` を組み立てる。置き場の外に出る slug は拒否する。

    slug / analysis_id は呼び出し側（LLM）由来なので、素通しすると `../../..` や絶対パスで
    任意ファイルの読み書きに化ける。ここをパス組み立ての唯一の口にして封じる。
    書き込み側の make_slug とは別に、読み取り側にも効かせる必要がある。
    """
    directory = Path(directory)
    text = str(slug)
    if not text or text != text.strip() or set(text) <= {"."}:
        raise ValueError(f"不正なノート識別子です（空・空白・ドットのみ）: {slug!r}")
    if "/" in text or "\\" in text or "\x00" in text:
        raise ValueError(f"不正なノート識別子です（パス区切りを含む）: {slug!r}")
    path = directory / f"{text}.md"
    # ドライブレター相対（`C:name`）等、区切り文字を持たない脱出も残らず弾く。
    if path.resolve().parent != directory.resolve():
        raise ValueError(f"ノート置き場の外を指す識別子です: {slug!r}")
    return path


def write_note(directory: str | Path, slug: str, meta: dict, body: str) -> Path:
    """frontmatter 付きノートを書き出す（ディレクトリは必要なら作成）。

    slug の検証はロック取得（= mkdir）より**前**に行う。あとで検証すると、拒否した
    呼び出しが空の置き場だけを残す。
    """
    directory = Path(directory)
    path = note_path(directory, slug)
    with _directory_lock(directory):
        path.write_text(
            dump_frontmatter(meta) + "\n\n" + (body or "").strip() + "\n",
            encoding="utf-8",
        )
    return path


_SLUG_RE = re.compile(r"[^a-z0-9]+")


def make_slug(text: str, maxlen: int = 60) -> str:
    """任意文字列から ASCII の slug を作る。"""
    slug = _SLUG_RE.sub("-", (text or "").lower()).strip("-")
    return slug[:maxlen].strip("-") or "note"


# --- objective レコード（reports/<analysis_id>.objective.md） ---
# 小問は本文の "- Q1: ..." 形式（ラベル Q1, Q2, ...）。

_OBJ_SUBQ_RE = re.compile(r"-\s*Q(\d+)\s*[:：]\s*(.+)")
_SUBQ_HEADER = "## 派生する小問"
_EMERGENT_HEADER = "## 創発的に追記された小問"
_PLACEHOLDERS = {"（なし）", "（記載なし）", "（未確認。確認後に記入）"}

#: objective ノートのファイル名の接尾辞。レポート（`<make_slug(id)>.md`）と同じ
#: reports/ に置くので、同じ analysis_id でも互いを上書きしないよう名前を分ける
#: （make_slug は "." を残さないので、レポートの名前とは衝突しない）。
OBJECTIVE_SUFFIX = ".objective"


def objective_slug(analysis_id: str) -> str:
    """objective ノートの slug（`<analysis_id>.objective`）。"""
    return f"{analysis_id}{OBJECTIVE_SUFFIX}"


def parse_objective(path: str | Path) -> tuple[dict, list[tuple[str, str]], str]:
    """objective を (meta, [(label, text), ...], body) で返す。label は "Q1" など。"""
    text = Path(path).read_text(encoding="utf-8")
    meta, body = parse_frontmatter(text)
    subqs = [(f"Q{m.group(1)}", m.group(2).strip()) for m in _OBJ_SUBQ_RE.finditer(body)]
    return meta, subqs, body


def build_objective_body(
    confirmed_objective: str, sub_questions: list[str], evidence: list[str] | None = None
) -> str:
    lines = ["## 確定目的", confirmed_objective.strip() if confirmed_objective else "（未確認。確認後に記入）", ""]
    lines.append("## 推測の根拠（データ特徴）")
    lines.extend([f"- {e}" for e in evidence] if evidence else ["（記載なし）"])
    lines += ["", _SUBQ_HEADER]
    for i, question in enumerate(sub_questions, 1):
        lines.append(f"- Q{i}: {question}")
    lines += ["", _EMERGENT_HEADER, "（なし）"]
    return "\n".join(lines)


def write_objective(
    directory: str | Path,
    analysis_id: str,
    meta_fields: dict,
    sub_questions: list[str],
    evidence: list[str] | None = None,
) -> Path:
    """objective レコードを <directory>/<analysis_id>.objective.md として書き出す。"""
    meta = {"type": "objective", "analysis_id": analysis_id}
    meta.update(meta_fields)
    meta.setdefault("status", "active")
    meta["confirmed"] = bool(meta.get("confirmed_objective"))
    body = build_objective_body(meta.get("confirmed_objective", ""), sub_questions, evidence)
    return write_note(directory, objective_slug(analysis_id), meta, body)


def update_objective_meta(path: str | Path, updates: dict) -> Path:
    """frontmatter のみ更新し本文は温存する。"""
    path = Path(path)
    meta, body = parse_frontmatter(path.read_text(encoding="utf-8"))
    meta.update(updates)
    if "confirmed_objective" in updates:
        meta["confirmed"] = bool(updates["confirmed_objective"])
    return write_note(path.parent, path.stem, meta, body)


def _append_under_section(body: str, header: str, new_line: str) -> str:
    """指定セクション末尾に行を追記（プレースホルダは除去）。無ければ末尾に新設。"""
    lines = body.splitlines()
    out: list[str] = []
    inserted = False
    i = 0
    while i < len(lines):
        out.append(lines[i])
        if not inserted and lines[i].strip() == header:
            j = i + 1
            section: list[str] = []
            while j < len(lines) and not lines[j].startswith("## "):
                section.append(lines[j])
                j += 1
            cleaned = [s for s in section if s.strip() not in _PLACEHOLDERS]
            while cleaned and cleaned[-1].strip() == "":
                cleaned.pop()
            cleaned.append(new_line)
            cleaned.append("")
            out.extend(cleaned)
            inserted = True
            i = j
            continue
        i += 1
    if not inserted:
        out += ["", header, new_line, ""]
    return "\n".join(out)


def add_subquestions(path: str | Path, texts: list[str]) -> Path:
    """創発的な小問を採番（既存 Q の最大+1）して本文に追記する。"""
    path = Path(path)
    meta, subqs, body = parse_objective(path)
    next_n = (max(int(label[1:]) for label, _ in subqs) if subqs else 0) + 1
    new_body = body
    for text in texts:
        new_body = _append_under_section(new_body, _EMERGENT_HEADER, f"- Q{next_n}: {text}")
        next_n += 1
    return write_note(path.parent, path.stem, meta, new_body)
