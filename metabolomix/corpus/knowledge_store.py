"""蓄積ノート（knowledge / playbook / analyses）を読む純ロジック層。

MCP には依存しない。``server.py`` がこの関数群を ``@mcp.resource`` から呼ぶ。

各ノートは frontmatter 付き markdown。索引（INDEX）は実ファイルを持たず、
各ノートの frontmatter から **読まれた瞬間に動的生成** する（単一の真実の源＝
frontmatter、ドリフトゼロ）。関連ノートの取得は ``[[link]]`` グラフを
**構造予算付き**（max 1 hop / max 5 本体 / 約 15k トークン）で展開する。

設計上の分業:
    - 構造予算（ホップ数・本体数・トークン上限）= ここでサーバが機械強制する
    - 関連性ゲート（どの隣接を採用するか）= LLM が判断する（指示側）
      → 予算を超えた隣接は捨てず索引行に降格するので、存在は常に可視。
        被害は「本体を開かなかった」だけで回復可能（小さい側に倒す）。
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import os
import re
import time

# --- 構造予算（サーバが機械強制する上限） ---
MAX_HOP = 1            # 展開は直接リンク先までの1ホップ
MAX_BODIES = 5         # 返す本体数の上限（root を含む）
MAX_TOKENS = 15_000    # 返す総量の上限。arf2テーブル実績44k(窓22%)と競合しない水準
CHARS_PER_TOKEN = 4    # 粗いトークン見積り

_WIKILINK_RE = re.compile(r"\[\[([^\]]+)\]\]")
# frontmatter で明示される関連リンクのキー（本文の [[link]] と併せてグラフを成す）
_LINK_META_KEYS = ("related", "conflicts_with", "supersedes", "expected_biology")


# --- frontmatter パース ---
try:  # PyYAML があれば使う。無ければ管理下の最小サブセットを自前パース
    import yaml

    def _parse_meta(text: str) -> dict:
        data = yaml.safe_load(text) if text.strip() else {}
        return data if isinstance(data, dict) else {}

except ImportError:  # pragma: no cover - 環境依存のフォールバック

    def _parse_meta(text: str) -> dict:
        return _minimal_yaml(text)


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value


def _scalar(value: str):
    value = value.strip()
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_unquote(item) for item in inner.split(",")]
    return _unquote(value)


def _minimal_yaml(text: str) -> dict:
    """frontmatter の最小サブセット（スカラ・インラインリスト・1段ネスト）を解釈する。"""
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


# --- ノートモデル ---
@dataclass
class Note:
    slug: str
    path: Path
    meta: dict
    body: str

    @property
    def links(self) -> list[str]:
        """本文の ``[[link]]`` と frontmatter の関連キーから順序付きユニークなリンク先を返す。"""
        found: list[str] = []
        for match in _WIKILINK_RE.findall(self.body):
            found.append(match.strip())
        for key in _LINK_META_KEYS:
            value = self.meta.get(key)
            if isinstance(value, str):
                found.append(_strip_wikilink(value))
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, str):
                        found.append(_strip_wikilink(item))
        seen: set[str] = set()
        ordered: list[str] = []
        for slug in found:
            if slug and slug != self.slug and slug not in seen:
                seen.add(slug)
                ordered.append(slug)
        return ordered


def _strip_wikilink(value: str) -> str:
    return value.strip().strip("[]").strip()


def load_notes(directory: str | Path) -> dict[str, Note]:
    """ディレクトリ内の ``*.md`` を slug→Note で読み込む（ファイル名 stem が slug）。"""
    directory = Path(directory)
    notes: dict[str, Note] = {}
    if not directory.is_dir():
        return notes
    for path in sorted(directory.glob("*.md")):
        if path.name.upper() == "INDEX.MD":  # 実INDEXは持たない方針だが念のため除外
            continue
        meta, body = parse_frontmatter(path.read_text(encoding="utf-8"))
        notes[path.stem] = Note(path.stem, path, meta, body)
    return notes


def load_corpus(directories) -> dict[str, Note]:
    """複数ディレクトリを横断して slug→Note を統合する（リンク解決用）。

    playbook ノートは根拠の knowledge ノートを ``[[link]]`` で参照するため、
    リンク先は別コーパス（別ディレクトリ）に存在しうる。展開時はここで統合した
    全体から解決する。slug は全体で一意である前提（衝突時は先勝ち）。
    """
    corpus: dict[str, Note] = {}
    for directory in directories:
        for slug, note in load_notes(directory).items():
            corpus.setdefault(slug, note)
    return corpus


def _note_kind(note: Note) -> str:
    kind = note.meta.get("type")
    return kind if kind in ("knowledge", "playbook", "objective") else "knowledge"


# --- 索引（動的生成） ---
def _index_line(note: Note) -> str:
    kind = _note_kind(note)
    if kind == "playbook":
        return f"- {note.slug}: {note.meta.get('when_to_use', '(no when_to_use)')}"
    if kind == "objective":
        return f"- {note.slug}: {note.meta.get('comparison', '(objective)')}"
    desc = note.meta.get("description", "(no description)")
    strength = note.meta.get("claim_strength")
    suffix = f" [{strength}]" if strength else ""
    return f"- {note.slug}: {desc}{suffix}"


_INDEX_HEADER = {
    "knowledge": (
        "# knowledge index (論文由来の宣言的知識)\n"
        "各行は 1ノート。関連すると判断した slug だけを "
        "`lipidmix://knowledge/expand/<slug>` で本体取得すること。\n"
        "[strength] は claim_strength（established/suggested/speculative）。"
        "speculative を引用する際は必ずその旨を明示する。\n"
    ),
    "playbook": (
        "# playbook index (再利用可能な解析手順)\n"
        "各行は 1ノート。確定目的の小問に合致する slug だけを "
        "`lipidmix://playbook/expand/<slug>` で本体取得すること。\n"
    ),
}


def build_index(directory: str | Path, kind: str) -> str:
    """frontmatter から 1行索引を組み立てて返す（実INDEXファイルは作らない）。"""
    notes = load_notes(directory)
    header = _INDEX_HEADER.get(kind, f"# {kind} index\n")
    if not notes:
        return header + "\n(ノートはまだありません)\n"
    body = "\n".join(_index_line(note) for note in notes.values())
    return f"{header}\n{body}\n"


# --- グラフ展開（構造予算をサーバ強制） ---
def _render_body(note: Note) -> str:
    """本体を出典/手順メタ付きで描画する。各ノートは自身の type で描画する。

    LLM が引用根拠（knowledge の source）や使用ツール（playbook の tools）を
    その場で得られるようにする。
    """
    if _note_kind(note) == "playbook":
        tools = note.meta.get("tools", [])
        tools_str = ", ".join(tools) if isinstance(tools, list) else str(tools)
        meta_bits = [
            f"when_to_use: {note.meta.get('when_to_use', '')}",
            f"tools: {tools_str}",
        ]
    else:
        meta_bits = [
            f"source: {note.meta.get('source', '(MISSING — 出典なしノートは引用不可)')}",
            f"claim_strength: {note.meta.get('claim_strength', '?')}",
        ]
    meta_line = " | ".join(meta_bits)
    return f"## {note.slug}\n{meta_line}\n\n{note.body}\n"


def expand(slug: str, directories) -> str:
    """root ノート本体＋1ホップの隣接を、構造予算内で束ねて返す。

    ``directories`` で渡したコーパス全体からリンクを解決する（playbook→knowledge
    のクロスコーパス参照に対応）。隣接本体は予算（MAX_BODIES / MAX_TOKENS）内で
    貪欲に含め、溢れた分は捨てずに索引行へ降格する（存在を可視に保つ）。どれを
    実際に採用するかの関連性判断は受け手の LLM に委ねる。
    """
    corpus = load_corpus(directories)
    root = corpus.get(slug)
    if root is None:
        available = ", ".join(sorted(corpus)) or "(none)"
        return f"(not found: {slug})\n利用可能な slug: {available}\n"

    budget_chars = MAX_TOKENS * CHARS_PER_TOKEN
    sections: list[str] = []
    demoted: list[str] = []

    root_render = _render_body(root)
    sections.append(root_render)
    used_chars = len(root_render)
    bodies = 1

    for link in root.links:  # MAX_HOP = 1（推移展開しない）
        neighbor = corpus.get(link)
        if neighbor is None:
            demoted.append(f"- {link}: (missing — リンク先が存在しない)")
            continue
        candidate = _render_body(neighbor)
        if bodies < MAX_BODIES and used_chars + len(candidate) <= budget_chars:
            sections.append(candidate)
            used_chars += len(candidate)
            bodies += 1
        else:
            demoted.append(_index_line(neighbor))

    out = "\n".join(sections)
    if demoted:
        out += (
            "\n---\n# 予算外の隣接（索引行のみ。必要なら個別に expand すること）\n"
            + "\n".join(demoted)
            + "\n"
        )
    return out


# --- 保守支援: playbook の tools 参照を集める（実ツール名との整合チェック用） ---
def playbook_tool_references(directory: str | Path) -> dict[str, list[str]]:
    """各 playbook ノートが frontmatter で宣言する ``tools`` を slug→[tool] で返す。"""
    references: dict[str, list[str]] = {}
    for slug, note in load_notes(directory).items():
        tools = note.meta.get("tools", [])
        if isinstance(tools, str):
            tools = [tools]
        references[slug] = [t for t in tools if isinstance(t, str)]
    return references


# ======================================================================
# 文献探索（paper ingest）支援: frontmatter 書き出し / カバレッジ / _inbox
# ======================================================================

# --- frontmatter 書き出し（stage / promote 用） ---
def _dump_scalar(value) -> str:
    if isinstance(value, list):
        return "[" + ", ".join(_dump_scalar(v) for v in value) + "]"
    text = str(value)
    needs_quote = (
        text == ""
        or text.strip() != text
        or any(ch in text for ch in ':#[]{}",')
        or "\n" in text
    )
    if needs_quote:
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return text


def dump_frontmatter(meta: dict) -> str:
    """dict を ``---`` 囲みの frontmatter 文字列にする（1段ネストまで対応）。"""
    lines = ["---"]
    for key, value in meta.items():
        if isinstance(value, dict):
            lines.append(f"{key}:")
            for inner_key, inner_value in value.items():
                lines.append(f"  {inner_key}: {_dump_scalar(inner_value)}")
        else:
            lines.append(f"{key}: {_dump_scalar(value)}")
    lines.append("---")
    return "\n".join(lines)


LOCK_FILENAME = ".metabolomix.lock"
LOCK_TIMEOUT_SECONDS = 10.0
LOCK_POLL_SECONDS = 0.1


@contextmanager
def _directory_lock(directory: str | Path, timeout: float = LOCK_TIMEOUT_SECONDS):
    """Serialize note writes in one directory using an atomic lock file."""
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


@contextmanager
def _directory_locks(directories):
    """Acquire multiple directory locks in a stable order."""
    unique = sorted({str(Path(directory)) for directory in directories})
    stack = []
    try:
        for directory in unique:
            lock = _directory_lock(directory)
            lock.__enter__()
            stack.append(lock)
        yield
    finally:
        while stack:
            stack.pop().__exit__(None, None, None)


def note_path(directory: str | Path, slug: str) -> Path:
    """``<directory>/<slug>.md`` を組み立てる。置き場の外に出る slug は拒否する。

    slug / analysis_id は呼び出し側（LLM）由来であり、その材料には `_inbox` の抄録
    （paper_ingest が明示的に「非信頼データ」として扱うもの）が含まれる。素通しすると
    `../../..` や絶対パスで任意ファイルの読み書き・削除に化けるため、ここを唯一の
    パス組み立て口にして封じる。書き込み側の make_slug とは別に、読み取り/削除側にも
    効かせる必要がある。
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


def _write_note_unlocked(directory: str | Path, slug: str, meta: dict, body: str) -> Path:
    directory = Path(directory)
    path = note_path(directory, slug)
    directory.mkdir(parents=True, exist_ok=True)
    path.write_text(
        dump_frontmatter(meta) + "\n\n" + (body or "").strip() + "\n",
        encoding="utf-8",
    )
    return path


def write_note(directory: str | Path, slug: str, meta: dict, body: str) -> Path:
    """frontmatter 付きノートを書き出す（ディレクトリは必要なら作成）。

    slug の検証はロック取得より**前**に行う（promote / reject と同じ順序）。
    _directory_lock は mkdir を伴うので、あとで検証すると拒否した呼び出しが
    空の置き場だけを残す——ユーザが消した analyses/ が pytest のたびに
    復活していた経路がこれ。
    """
    directory = Path(directory)
    note_path(directory, slug)
    with _directory_lock(directory):
        return _write_note_unlocked(directory, slug, meta, body)


_SLUG_RE = re.compile(r"[^a-z0-9]+")


def make_slug(text: str, maxlen: int = 60) -> str:
    """任意文字列から ASCII の slug を作る（DOI/タイトル用）。"""
    slug = _SLUG_RE.sub("-", (text or "").lower()).strip("-")
    return slug[:maxlen].strip("-") or "note"


# --- 脂質クラス→検索語の語彙マップ（reference ノートで育てられる） ---
VOCAB_NOTE_SLUG = "lipid-class-search-vocabulary"
_DEFAULT_VOCAB: dict[str, list[str]] = {
    "pe p-": ["plasmalogen", "ethanolamine plasmalogen", "ether phospholipid"],
    "pc p-": ["plasmalogen", "choline plasmalogen", "ether phospholipid"],
    "pe o-": ["plasmanyl", "ether phospholipid"],
    "pc o-": ["plasmanyl", "ether phospholipid"],
    "tg": ["triacylglycerol", "triglyceride"],
    "dg": ["diacylglycerol"],
    "pc": ["phosphatidylcholine"],
    "pe": ["phosphatidylethanolamine"],
    "ps": ["phosphatidylserine"],
    "pi": ["phosphatidylinositol"],
    "pg": ["phosphatidylglycerol"],
    "fa": ["fatty acid"],
    "cer": ["ceramide"],
    "sm": ["sphingomyelin"],
    "lpc": ["lysophosphatidylcholine"],
    "che": ["cholesteryl ester", "cholesterol ester"],
}
_VOCAB_SEP_RE = re.compile(r"\s*(?:→|->|:)\s*")


def load_vocab(knowledge_dir: str | Path) -> dict[str, list[str]]:
    """既定語彙に reference ノート ``lipid-class-search-vocabulary`` の追記を重ねて返す。

    ノート本文は1行 ``KEY: syn1, syn2`` または ``KEY → syn1, syn2``（markdownの
    ``-``/``` ` ```/``*`` は許容）で記述する。
    """
    vocab = {k: list(v) for k, v in _DEFAULT_VOCAB.items()}
    note = load_notes(knowledge_dir).get(VOCAB_NOTE_SLUG)
    if note is None:
        return vocab
    for raw in note.body.splitlines():
        line = raw.strip().lstrip("-*").strip()
        if not line or line.startswith("#"):
            continue
        parts = _VOCAB_SEP_RE.split(line, maxsplit=1)
        if len(parts) != 2:
            continue
        key = parts[0].strip().strip("`").lower()
        syns = [s.strip().strip("`") for s in parts[1].split(",") if s.strip()]
        if key and syns:
            vocab[key] = syns
    return vocab


# --- カバレッジ判定（COVERED / WEAK / GAP の決定論ベースライン） ---
COVERAGE_COVERED_MIN = 0.15  # established 該当ノートがこの類似度以上なら COVERED
COVERAGE_WEAK_MIN = 0.06     # これ未満は該当なし＝GAP


def _normalize_for_match(text: str) -> str:
    return re.sub(r"\s+", "", (text or "").lower())


def _char_bigrams(text: str) -> set[str]:
    """言語非依存の文字bigram集合（日本語の無空白テキストにも有効）。"""
    norm = _normalize_for_match(text)
    if len(norm) < 2:
        return {norm} if norm else set()
    return {norm[i:i + 2] for i in range(len(norm) - 1)}


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def _expand_with_vocab(text: str, vocab: dict[str, list[str]]) -> str:
    low = (text or "").lower()
    extra: list[str] = []
    for key, syns in vocab.items():
        if key in low:
            extra.extend(syns)
    return f"{text} {' '.join(extra)}" if extra else text


def _note_match_text(note: Note) -> str:
    tags = note.meta.get("tags", [])
    tags_str = " ".join(tags) if isinstance(tags, list) else str(tags)
    return " ".join([note.meta.get("description", ""), tags_str, note.slug])


def coverage(
    subquestions: list[str],
    knowledge_dir: str | Path,
    vocab: dict[str, list[str]] | None = None,
) -> dict[str, dict]:
    """各小問 Qi を COVERED / WEAK / GAP に分類する決定論ベースライン。

    文字bigramのJaccard類似 + 語彙同義語展開で index と突合する。最終的な COVERED の
    真偽は呼び出し側（LLM）が該当ノートを ``expand`` して検証する前提。
    ``reference`` 型ノート（語彙マップ等）は突合対象から除外する。
    """
    if vocab is None:
        vocab = load_vocab(knowledge_dir)
    notes = {
        slug: note
        for slug, note in load_notes(knowledge_dir).items()
        if note.meta.get("type") != "reference"
    }
    note_grams = {
        slug: (_char_bigrams(_expand_with_vocab(_note_match_text(note), vocab)),
               note.meta.get("claim_strength", ""))
        for slug, note in notes.items()
    }

    result: dict[str, dict] = {}
    for question in subquestions:
        q_grams = _char_bigrams(_expand_with_vocab(question, vocab))
        scored = sorted(
            ((slug, round(_jaccard(q_grams, grams), 3), strength)
             for slug, (grams, strength) in note_grams.items()),
            key=lambda item: item[1],
            reverse=True,
        )
        best = scored[0] if scored else None
        if best is None or best[1] < COVERAGE_WEAK_MIN:
            state = "GAP"
        elif best[1] >= COVERAGE_COVERED_MIN and best[2] == "established":
            state = "COVERED"
        else:
            state = "WEAK"
        result[question] = {
            "state": state,
            "matches": [
                {"slug": s, "score": sc, "claim_strength": st}
                for s, sc, st in scored[:3] if sc >= COVERAGE_WEAK_MIN
            ],
        }
    return result


# --- _inbox（探索結果の隔離） ---
INBOX_DIRNAME = "_inbox"


def inbox_dir(knowledge_dir: str | Path) -> Path:
    return Path(knowledge_dir) / INBOX_DIRNAME


def load_inbox(knowledge_dir: str | Path) -> dict[str, Note]:
    return load_notes(inbox_dir(knowledge_dir))


def build_inbox_index(knowledge_dir: str | Path) -> str:
    """保留中ノートを analysis_id×Qi でグルーピングした1行索引で返す。"""
    notes = load_inbox(knowledge_dir)
    header = (
        "# knowledge inbox (pending review)\n"
        "speculative 隔離中。ingest_promote(slug, claim_strength) / ingest_reject(slug) で確定する。\n"
    )
    if not notes:
        return header + "\n(保留中のノートはありません)\n"

    def sort_key(item):
        _, note = item
        return (str(note.meta.get("found_for", "")), -float(note.meta.get("relevance_score", 0) or 0))

    lines: list[str] = []
    current_group = None
    for slug, note in sorted(notes.items(), key=sort_key):
        found_for = str(note.meta.get("found_for", "(unspecified)"))
        if found_for != current_group:
            lines.append(f"\n## {found_for}")
            current_group = found_for
        desc = note.meta.get("description") or note.meta.get("title") or slug
        score = note.meta.get("relevance_score", "?")
        query = note.meta.get("query", "")
        source = note.meta.get("source", "?")
        lines.append(f"- {slug} (score={score}) :: {desc}\n    query: {query} | source: {source}")
    return header + "\n".join(lines) + "\n"


def promote(slug: str, knowledge_dir: str | Path, claim_strength: str | None = None) -> Path:
    """``_inbox/<slug>.md`` を ``knowledge/<slug>.md`` へ昇格（status除去・強度設定可）。"""
    knowledge_dir = Path(knowledge_dir)
    inbox = inbox_dir(knowledge_dir)
    src = note_path(inbox, slug)
    with _directory_locks([knowledge_dir, inbox]):
        if not src.is_file():
            raise FileNotFoundError(f"inbox note not found: {slug}")
        meta, body = parse_frontmatter(src.read_text(encoding="utf-8"))
        meta.pop("status", None)
        meta.setdefault("type", "knowledge")
        if claim_strength:
            meta["claim_strength"] = claim_strength
        dest = _write_note_unlocked(knowledge_dir, slug, meta, body)
        src.unlink()
        return dest


def reject(slug: str, knowledge_dir: str | Path) -> bool:
    """``_inbox/<slug>.md`` を破棄する。存在しなければ False。"""
    inbox = inbox_dir(knowledge_dir)
    src = note_path(inbox, slug)
    with _directory_lock(inbox):
        if src.is_file():
            src.unlink()
            return True
        return False


# ======================================================================
# objective レコード（analyses/）のライフサイクル
# ======================================================================
# 小問は本文の "- Q1: ..." 形式（ラベル Q1,Q2... は found_for/search_log と統一）。
# 探索ログ(search_log)は frontmatter ではなく本文 "## 探索ログ（search_log）" 節に
# "- Q2 | date | query | hits=.. | promoted=.." 形式で持つ（クエリのカンマ等で
# frontmatter インラインリストが壊れるのを避けるため）。

_OBJ_SUBQ_RE = re.compile(r"-\s*Q(\d+)\s*[:：]\s*(.+)")
_OBJ_LOG_RE = re.compile(r"-\s*(Q\d+)\s*\|\s*(.+)")
_SUBQ_HEADER = "## 派生する小問（関連性ゲートの照合対象）"
_EMERGENT_HEADER = "## 創発的に追記された小問"
_SEARCHLOG_HEADER = "## 探索ログ（search_log）"
_PLACEHOLDERS = {"（なし）", "（記載なし）", "（未確認。確認後に記入）"}


def parse_objective(path: str | Path) -> tuple[dict, list[tuple[str, str]], str]:
    """objective を (meta, [(label, text), ...], body) で返す。label は "Q1" など。"""
    text = Path(path).read_text(encoding="utf-8")
    meta, body = parse_frontmatter(text)
    subqs = [(f"Q{m.group(1)}", m.group(2).strip()) for m in _OBJ_SUBQ_RE.finditer(body)]
    return meta, subqs, body


def searched_labels(path: str | Path) -> dict[str, str]:
    """探索ログ節から {label: ログ行} を返す（churn 防止の突合用）。"""
    _meta, body = parse_frontmatter(Path(path).read_text(encoding="utf-8"))
    result: dict[str, str] = {}
    for raw in body.splitlines():
        m = _OBJ_LOG_RE.match(raw.strip())
        if m:
            result.setdefault(m.group(1), raw.strip())
    return result


def build_objective_body(
    confirmed_objective: str, sub_questions: list[str], evidence: list[str] | None = None
) -> str:
    lines = ["## 確定目的", confirmed_objective.strip() if confirmed_objective else "（未確認。確認後に記入）", ""]
    lines.append("## 推測の根拠（データ特徴）")
    lines.extend([f"- {e}" for e in evidence] if evidence else ["（記載なし）"])
    lines += ["", _SUBQ_HEADER]
    for i, question in enumerate(sub_questions, 1):
        lines.append(f"- Q{i}: {question}")
    lines += ["", _EMERGENT_HEADER, "（なし）", "", _SEARCHLOG_HEADER, "（なし）"]
    return "\n".join(lines)


def write_objective(
    analyses_dir: str | Path,
    analysis_id: str,
    meta_fields: dict,
    sub_questions: list[str],
    evidence: list[str] | None = None,
) -> Path:
    """objective レコードを analyses/<analysis_id>.md として書き出す。"""
    meta = {"type": "objective", "analysis_id": analysis_id}
    meta.update(meta_fields)
    meta.setdefault("status", "active")
    meta["confirmed"] = bool(meta.get("confirmed_objective"))
    body = build_objective_body(meta.get("confirmed_objective", ""), sub_questions, evidence)
    return write_note(analyses_dir, analysis_id, meta, body)


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


def append_search_log(path: str | Path, label: str, date: str, query: str, hits, promoted) -> Path:
    """探索ログ節に1件追記する（churn 防止の記録）。"""
    path = Path(path)
    meta, body = parse_frontmatter(path.read_text(encoding="utf-8"))
    safe_query = str(query).replace("|", "/").replace("\n", " ").strip()
    line = f"- {label} | {date} | {safe_query} | hits={hits} | promoted={promoted}"
    new_body = _append_under_section(body, _SEARCHLOG_HEADER, line)
    return write_note(path.parent, path.stem, meta, new_body)


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
