"""消したツールの名前が、入口の文書とコードの文字列・コメントに残らないこと（Review Focus 5）。

B-2b で旧 ARF / ARF2 / PAI2 / EIC / 試料 / 文献のツールと `dataset_load` `dataset_handle_*`
`update_objective` `list_reports` が消えた。文書や `required_tools`・案内文・docstring に名前が
残ると、存在しないツールを呼ばせる。単語境界で数える（`_list_reports` や `dataset_loaded` は当たらない）。

対象外: `docs/superpowers/**`（設計・計画の履歴）と `docs/schema/**`（Key 番号表）。
"""
import io
import re
import tokenize
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REMOVED_TOOLS = """
arf2_annotate_identities arf2_parser arf_differential arf_exclude arf_export_differential
arf_export_evidence_bundle arf_list_classes arf_list_sample_roles arf_list_tags arf_parser
arf_pca_preprocessed arf_pca_species arf_plot_group_intensity arf_plot_species arf_plot_volcano
arf_preprocess dataset_handle_close dataset_handle_list dataset_handle_open dataset_load dcl_parser
eic_parser eic_plot_chromatograms eic_plot_compounds eic_search ingest_promote ingest_reject
ingest_review_queue ingest_stage knowledge_coverage list_reports log_search pai2_inspect_peak
pai2_parser paper_search plot_pca_loadings sample_search update_objective verify_peak_annotation
""".split()

WORD = re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(map(re.escape, REMOVED_TOOLS)) + r")(?![A-Za-z0-9_])")

#: 例外（ファイル相対パス -> 許す名前の集合）。理由を必ず書く。増やさず、文面を直すこと。
#: いまは無い。（このテスト自身は `tests/` にあり、走査の対象外。）
EXCEPTIONS: dict[str, set[str]] = {}

DOCS = ["README.md", "USAGE.md", "DEPLOY.md", "CLAUDE.md", "docs/cli.md"]


def _doc_files():
    files = [ROOT / name for name in DOCS]
    files += sorted((ROOT / "docs" / "workflow").rglob("*.md"))
    files += sorted((ROOT / "docs" / "output_format").rglob("*.md"))
    return [f for f in files if f.is_file()]


def _code_texts(path: Path):
    """`.py` の文字列・コメントの本文（docstring を含む）。名前そのもの（識別子）は数えない。"""
    source = path.read_text(encoding="utf-8")
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, IndentationError):
        yield source
        return
    for token in tokens:
        if token.type in (tokenize.STRING, tokenize.COMMENT):
            yield token.string
        elif getattr(tokenize, "FSTRING_MIDDLE", None) is not None and token.type == tokenize.FSTRING_MIDDLE:
            yield token.string


def _hits(path: Path, texts) -> set[str]:
    rel = path.relative_to(ROOT).as_posix()
    allowed = EXCEPTIONS.get(rel, set())
    found = set()
    for text in texts:
        found.update(WORD.findall(text))
    return found - allowed


def test_removed_tools_are_not_mentioned():
    offenders = {}
    for path in _doc_files():
        hits = _hits(path, [path.read_text(encoding="utf-8")])
        if hits:
            offenders[path.relative_to(ROOT).as_posix()] = sorted(hits)
    for path in sorted((ROOT / "metabolomix").rglob("*.py")):
        hits = _hits(path, _code_texts(path))
        if hits:
            offenders[path.relative_to(ROOT).as_posix()] = sorted(hits)
    assert not offenders, "消したツールの名前が残っている:\n" + "\n".join(
        f"  {path}: {', '.join(names)}" for path, names in offenders.items())


def test_word_boundary_has_no_false_positives():
    assert not WORD.search("_list_reports dataset_loaded my_dataset_load_x")
    assert WORD.search("`dataset_load` で読む")
    assert WORD.search("arf_parser()")
