"""ノートの読み書き（metabolomix/tools/notes.py）と objective レコードの検証。

文献と知識の蓄積の層を消したあとも、解析フォルダ配下 reports/ のレポートと objective は
frontmatter 付き Markdown で読み書きする。その純ロジックだけをここで縛る。
"""

import ast
import sys
import tempfile
import unittest
from pathlib import Path

from metabolomix.tools import notes

PACKAGE = Path(__file__).resolve().parent.parent / "metabolomix"


def test_notes_module_is_stdlib_only():
    tree = ast.parse(Path(notes.__file__).read_text(encoding="utf-8"))
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots |= {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            roots.add(node.module.split(".")[0])
    assert roots - {"__future__"} <= set(sys.stdlib_module_names)


def test_no_corpus_imports():
    """消した metabolomix.corpus を、どこも import していない。"""
    assert not (PACKAGE / "corpus").exists()
    offenders = []
    for path in PACKAGE.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""] + [
                    f"{node.module}.{alias.name}" for alias in node.names]
            if any(name.split(".")[:2] == ["metabolomix", "corpus"] for name in names):
                offenders.append(str(path))
    assert not offenders


class FrontmatterTests(unittest.TestCase):
    def test_parses_scalar_inline_list_and_nested_map(self):
        text = (
            "---\n"
            "type: report\n"
            'title: "クラス別に2群を見る"\n'
            "tools: [a_tool, b_tool]\n"
            "precondition:\n"
            "  polarity: any\n"
            "  group_structure: 2-group comparison\n"
            "---\n"
            "## 本文\n"
        )
        meta, body = notes.parse_frontmatter(text)
        self.assertEqual(meta["type"], "report")
        self.assertEqual(meta["title"], "クラス別に2群を見る")
        self.assertEqual(meta["tools"], ["a_tool", "b_tool"])
        self.assertEqual(meta["precondition"]["group_structure"], "2-group comparison")
        self.assertIn("## 本文", body)

    def test_no_frontmatter_returns_empty_meta(self):
        meta, body = notes.parse_frontmatter("# just a heading\ncontent")
        self.assertEqual(meta, {})
        self.assertIn("just a heading", body)

    def test_roundtrip_scalar_list_nested(self):
        meta = {
            "type": "report",
            "title": 'A title: with colon, and comma and "quote"',
            "groups": ["control, batch 1", "treated"],
            "empty": [],
            "flag": True,
            "precondition": {"polarity": "any", "group_structure": "2-group comparison"},
        }
        parsed, body = notes.parse_frontmatter(notes.dump_frontmatter(meta) + "\n\nbody")
        self.assertEqual(parsed, meta)
        self.assertEqual(body, "body")


class SlugTests(unittest.TestCase):
    def test_make_slug_is_ascii_and_stable(self):
        self.assertEqual(notes.make_slug("Group A vs B"), "group-a-vs-b")
        self.assertEqual(notes.make_slug("group-a-vs-b"), "group-a-vs-b")
        self.assertEqual(notes.make_slug("日本語だけ"), "note")
        self.assertLessEqual(len(notes.make_slug("x" * 200)), 60)


class SlugContainmentTests(unittest.TestCase):
    """analysis_id は呼び出し側（LLM）由来なので、置き場の外へ出られてはならない。"""

    def test_write_note_refuses_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            with self.assertRaises(ValueError):
                notes.write_note(directory, "../escaped", {"type": "report"}, "body")
            self.assertFalse((Path(tmp) / "escaped.md").exists())

    def test_write_note_refuses_backslash_and_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            for bad in ("..\\escaped", "", ".", "..", "sub/note"):
                with self.subTest(slug=bad), self.assertRaises(ValueError):
                    notes.write_note(directory, bad, {"type": "report"}, "body")

    def test_write_note_rejects_before_creating_directory(self):
        """拒否される slug では置き場そのものを作らない（消した置き場が復活しない）。"""
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            with self.assertRaises(ValueError):
                notes.write_note(directory, "../escaped", {"type": "report"}, "body")
            self.assertFalse(directory.exists())

    def test_normal_slug_still_writes_and_removes_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            path = notes.write_note(directory, "exp-1", {"type": "report"}, "body")
            self.assertEqual(path, directory / "exp-1.md")
            self.assertTrue(path.is_file())
            self.assertFalse((directory / notes.LOCK_FILENAME).exists())

    def test_directory_lock_times_out_when_already_locked(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / notes.LOCK_FILENAME).write_text("held", encoding="utf-8")
            with self.assertRaises(TimeoutError):
                with notes._directory_lock(directory, timeout=0.01):
                    pass


class ObjectiveStoreTests(unittest.TestCase):
    def _write(self, directory: Path) -> Path:
        return notes.write_objective(
            directory,
            "exp-1",
            {
                "dataset": "NEG / x",
                "polarity": "NEG",
                "groups": ["control", "treatment"],
                "comparison": "群間差",
                "biological_context": "",
                "inferred_objective": "推測",
                "confirmed_objective": "",
                "expected_biology": [],
            },
            ["どのクラスに差が出るか", "エーテル脂質は減るか"],
        )

    def test_write_and_parse_labeled_subquestions(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            self.assertEqual(path.name, "exp-1.objective.md")
            meta, subqs, _body = notes.parse_objective(path)
            self.assertEqual(meta["analysis_id"], "exp-1")
            self.assertEqual(meta["status"], "active")
            self.assertEqual([label for label, _ in subqs], ["Q1", "Q2"])
            self.assertIn("どのクラス", subqs[0][1])
            self.assertIs(meta["confirmed"], False)

    def test_update_meta_preserves_body(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            notes.update_objective_meta(
                path, {"confirmed_objective": "確定した目的", "biological_context": "RAW macrophage"})
            meta, subqs, _ = notes.parse_objective(path)
            self.assertEqual(meta["confirmed_objective"], "確定した目的")
            self.assertEqual(meta["biological_context"], "RAW macrophage")
            self.assertIs(meta["confirmed"], True)
            self.assertEqual(len(subqs), 2)  # 本文の小問は温存

    def test_add_subquestions_numbers_after_existing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            notes.add_subquestions(path, ["創発した問い"])
            _meta, subqs, _ = notes.parse_objective(path)
            self.assertEqual([label for label, _ in subqs], ["Q1", "Q2", "Q3"])
            self.assertIn("創発した問い", dict(subqs)["Q3"])

    def test_objective_slug_does_not_collide_with_report_slug(self):
        self.assertNotEqual(notes.objective_slug("exp-1"), notes.make_slug("exp-1"))


if __name__ == "__main__":
    unittest.main()
