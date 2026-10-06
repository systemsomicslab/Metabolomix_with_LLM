import tempfile
import unittest
from pathlib import Path

from metabolomix.corpus import knowledge_store
from metabolomix.corpus.knowledge_store import (
    Note,
    build_index,
    coverage,
    dump_frontmatter,
    expand,
    load_vocab,
    parse_frontmatter,
    playbook_tool_references,
    promote,
    reject,
    write_note as ks_write_note,
)


def write_note(directory: Path, slug: str, frontmatter: str, body: str = "") -> None:
    text = f"---\n{frontmatter.strip()}\n---\n{body}"
    (directory / f"{slug}.md").write_text(text, encoding="utf-8")


class FrontmatterTests(unittest.TestCase):
    def test_parses_scalar_inline_list_and_nested_map(self):
        text = (
            "---\n"
            "type: playbook\n"
            'when_to_use: "クラス別に2群を見る"\n'
            "tools: [arf2_parser, arf_parser]\n"
            "precondition:\n"
            "  polarity: any\n"
            "  group_structure: 2-group comparison\n"
            "---\n"
            "## 本文\n"
        )
        meta, body = parse_frontmatter(text)
        self.assertEqual(meta["type"], "playbook")
        self.assertEqual(meta["when_to_use"], "クラス別に2群を見る")
        self.assertEqual(meta["tools"], ["arf2_parser", "arf_parser"])
        self.assertEqual(meta["precondition"]["polarity"], "any")
        self.assertEqual(meta["precondition"]["group_structure"], "2-group comparison")
        self.assertIn("## 本文", body)

    def test_no_frontmatter_returns_empty_meta(self):
        meta, body = parse_frontmatter("# just a heading\ncontent")
        self.assertEqual(meta, {})
        self.assertIn("just a heading", body)


class NoteLinkTests(unittest.TestCase):
    def test_links_from_body_and_meta_dedup_and_drop_self(self):
        note = Note(
            slug="a",
            path=Path("a.md"),
            meta={"conflicts_with": ["b"], "expected_biology": ["c"]},
            body="see [[b]] and [[d]] and [[a]]",
        )
        self.assertEqual(note.links, ["b", "d", "c"])


class IndexTests(unittest.TestCase):
    def test_knowledge_index_shows_description_and_strength(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_note(
                directory,
                "note-x",
                'type: knowledge\ndescription: "X についての主張"\n'
                "claim_strength: speculative\nsource: \"Doe 2020\"",
            )
            index = build_index(directory, "knowledge")
            self.assertIn("- note-x: X についての主張 [speculative]", index)

    def test_playbook_index_shows_when_to_use(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_note(
                directory,
                "flow-y",
                'type: playbook\nwhen_to_use: "Y したいとき"\ntools: [arf_parser]',
            )
            index = build_index(directory, "playbook")
            self.assertIn("- flow-y: Y したいとき", index)

    def test_empty_directory_index_is_safe(self):
        with tempfile.TemporaryDirectory() as tmp:
            index = build_index(Path(tmp), "knowledge")
            self.assertIn("ノートはまだありません", index)


class ExpandTests(unittest.TestCase):
    def _make_linked_corpus(self, directory: Path) -> None:
        write_note(
            directory,
            "root",
            'type: knowledge\ndescription: "root"\nsource: "S"\nclaim_strength: established',
            body="links to [[child1]] and [[child2]]",
        )
        write_note(directory, "child1", 'type: knowledge\ndescription: "c1"\nsource: "S1"')
        write_note(directory, "child2", 'type: knowledge\ndescription: "c2"\nsource: "S2"')

    def test_expand_returns_root_and_one_hop_neighbors(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._make_linked_corpus(directory)
            out = expand("root", [directory])
            self.assertIn("## root", out)
            self.assertIn("## child1", out)
            self.assertIn("## child2", out)
            self.assertIn("source: S1", out)  # 出典が本体に同梱される

    def test_expand_resolves_links_across_directories(self):
        # playbook→knowledge のクロスコーパス参照を解決できること
        with tempfile.TemporaryDirectory() as tmp_pb, tempfile.TemporaryDirectory() as tmp_kn:
            pb, kn = Path(tmp_pb), Path(tmp_kn)
            write_note(
                pb,
                "flow",
                "type: playbook\nwhen_to_use: x\ntools: [arf_parser]",
                body="根拠は [[fact]] を参照",
            )
            write_note(kn, "fact", 'type: knowledge\ndescription: "事実"\nsource: "S"')
            out = expand("flow", [pb, kn])
            self.assertIn("## flow", out)
            self.assertIn("## fact", out)  # 別ディレクトリの本体が解決される
            self.assertIn("source: S", out)

    def test_missing_slug_lists_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._make_linked_corpus(directory)
            out = expand("nope", [directory])
            self.assertIn("not found: nope", out)
            self.assertIn("root", out)

    def test_missing_link_is_reported_not_crashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_note(
                directory,
                "root",
                'type: knowledge\ndescription: "root"\nsource: "S"',
                body="dangling [[ghost]]",
            )
            out = expand("root", [directory])
            self.assertIn("ghost", out)
            self.assertIn("missing", out)

    def test_budget_demotes_overflow_to_index_lines(self):
        # MAX_BODIES を小さくして、超過隣接が索引行に降格することを確認
        original = knowledge_store.MAX_BODIES
        knowledge_store.MAX_BODIES = 2  # root + 1 隣接のみ本体、残りは降格
        try:
            with tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                write_note(
                    directory,
                    "root",
                    'type: knowledge\ndescription: "root"\nsource: "S"',
                    body="[[c1]] [[c2]] [[c3]]",
                )
                for slug in ("c1", "c2", "c3"):
                    write_note(
                        directory, slug, f'type: knowledge\ndescription: "{slug}"\nsource: "S"'
                    )
                out = expand("root", [directory])
                # root + c1 は本体、c2/c3 は索引行に降格
                self.assertIn("## root", out)
                self.assertIn("## c1", out)
                self.assertIn("予算外の隣接", out)
                self.assertNotIn("## c2", out)
                self.assertIn("- c2:", out)
        finally:
            knowledge_store.MAX_BODIES = original


class PlaybookToolReferenceTests(unittest.TestCase):
    def test_collects_tools_per_slug(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_note(
                directory,
                "flow",
                "type: playbook\nwhen_to_use: x\ntools: [arf2_parser, arf_parser]",
            )
            refs = playbook_tool_references(directory)
            self.assertEqual(refs["flow"], ["arf2_parser", "arf_parser"])


class CoverageTests(unittest.TestCase):
    def _corpus(self, directory: Path) -> None:
        write_note(
            directory, "plasmalogen-oxidation",
            'type: knowledge\ndescription: "プラスマローゲンは酸化ストレスで選択的に減少する"\n'
            'claim_strength: established\nsource: "S"',
        )
        write_note(
            directory, "ceramide-apoptosis",
            'type: knowledge\ndescription: "セラミドはアポトーシスに関与する"\n'
            'claim_strength: speculative\nsource: "S"',
        )
        # reference ノートは突合対象外であることの確認用
        write_note(
            directory, "lipid-class-search-vocabulary",
            'type: reference\ndescription: "語彙マップ"',
        )

    def test_covered_weak_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._corpus(directory)
            result = coverage(
                [
                    "プラスマローゲンは酸化ストレスで減少するか",   # established 該当 -> COVERED
                    "セラミドはアポトーシスに関与するか",          # speculative 該当 -> WEAK
                    "解糖系のグルコース取り込み速度はどうか",       # 無関係 -> GAP
                ],
                directory,
            )
            states = {q: info["state"] for q, info in result.items()}
            self.assertEqual(states["プラスマローゲンは酸化ストレスで減少するか"], "COVERED")
            self.assertEqual(states["セラミドはアポトーシスに関与するか"], "WEAK")
            self.assertEqual(states["解糖系のグルコース取り込み速度はどうか"], "GAP")

    def test_reference_note_excluded_from_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._corpus(directory)
            result = coverage(["語彙マップ"], directory)
            matched = [m["slug"] for m in result["語彙マップ"]["matches"]]
            self.assertNotIn("lipid-class-search-vocabulary", matched)


class VocabTests(unittest.TestCase):
    def test_default_and_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            # 既定のみ
            vocab = load_vocab(directory)
            self.assertIn("plasmalogen", vocab["pe p-"])
            # reference ノートで追記/上書き
            write_note(
                directory, "lipid-class-search-vocabulary",
                'type: reference\ndescription: "v"',
                body="- HEX: hexosylceramide, glucosylceramide\nPE P-: my-synonym",
            )
            vocab = load_vocab(directory)
            self.assertIn("hexosylceramide", vocab["hex"])
            self.assertEqual(vocab["pe p-"], ["my-synonym"])


class PromoteRejectTests(unittest.TestCase):
    def _stage(self, knowledge_dir: Path, slug: str) -> None:
        ks_write_note(
            knowledge_store.inbox_dir(knowledge_dir), slug,
            {"type": "knowledge", "status": "pending", "claim_strength": "speculative",
             "title": "t", "source": "S", "found_for": "A/Q1"},
            "## Abstract\nbody",
        )

    def test_promote_moves_and_sets_strength(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._stage(directory, "paper-x")
            dest = promote("paper-x", directory, claim_strength="suggested")
            self.assertTrue(dest.is_file())
            self.assertEqual(dest.parent, directory)  # knowledge/ 直下へ
            self.assertFalse((knowledge_store.inbox_dir(directory) / "paper-x.md").exists())
            meta, _ = parse_frontmatter(dest.read_text(encoding="utf-8"))
            self.assertEqual(meta["claim_strength"], "suggested")
            self.assertNotIn("status", meta)  # status は除去

    def test_reject_deletes(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._stage(directory, "paper-y")
            self.assertTrue(reject("paper-y", directory))
            self.assertFalse((knowledge_store.inbox_dir(directory) / "paper-y.md").exists())
            self.assertFalse(reject("missing", directory))


class SlugContainmentTests(unittest.TestCase):
    """slug/analysis_id は呼び出し側（LLM）由来なので、ノート置き場の外へ出られてはならない。

    _inbox の抄録は非信頼データ（プロンプトインジェクション経路）であり、そこから
    導かれた識別子がそのままファイルパスになると任意ファイルの読み書き/削除に化ける。
    """

    def _outside_file(self, tmp: str) -> Path:
        """ノート置き場の1つ上に、消えては困る .md を置く。"""
        outside = Path(tmp) / "outside.md"
        outside.write_text("must survive", encoding="utf-8")
        return outside

    def test_reject_refuses_traversal_and_keeps_outside_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            outside = self._outside_file(tmp)
            knowledge_dir = Path(tmp) / "knowledge"
            knowledge_store.inbox_dir(knowledge_dir).mkdir(parents=True, exist_ok=True)
            with self.assertRaises(ValueError):
                reject("../../outside", knowledge_dir)
            self.assertTrue(outside.is_file())

    def test_promote_refuses_traversal_and_keeps_outside_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            outside = self._outside_file(tmp)
            knowledge_dir = Path(tmp) / "knowledge"
            knowledge_store.inbox_dir(knowledge_dir).mkdir(parents=True, exist_ok=True)
            with self.assertRaises(ValueError):
                promote("../../outside", knowledge_dir)
            self.assertTrue(outside.is_file())

    def test_write_note_refuses_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            with self.assertRaises(ValueError):
                ks_write_note(directory, "../escaped", {"type": "knowledge"}, "body")
            self.assertFalse((Path(tmp) / "escaped.md").exists())

    def test_write_note_refuses_backslash_and_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            for bad in ("..\\escaped", "", ".", "..", "sub/note"):
                with self.subTest(slug=bad), self.assertRaises(ValueError):
                    ks_write_note(directory, bad, {"type": "knowledge"}, "body")

    def test_write_note_rejects_before_creating_directory(self):
        """拒否される slug では置き場そのものを作らない。

        write_note は検証より先にディレクトリロック（= mkdir）を取っていたため、
        traversal を弾いたあとに空の置き場だけが残った。ユーザが消した analyses/ が
        pytest のたびに復活する経路がこれ。promote / reject は先に note_path で
        検証しており、write_note だけが例外だった。
        """
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            with self.assertRaises(ValueError):
                ks_write_note(directory, "../escaped", {"type": "knowledge"}, "body")
            self.assertFalse(directory.exists(), f"拒否したのに置き場が作られた: {directory}")

    def test_normal_slug_still_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "notes"
            path = ks_write_note(directory, "pe-p-vs-pe-o", {"type": "knowledge"}, "body")
            self.assertEqual(path, directory / "pe-p-vs-pe-o.md")
            self.assertTrue(path.is_file())


class NoteLockTests(unittest.TestCase):
    def test_write_note_removes_lock_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            ks_write_note(directory, "note", {"type": "knowledge"}, "body")
            self.assertTrue((directory / "note.md").is_file())
            self.assertFalse((directory / knowledge_store.LOCK_FILENAME).exists())

    def test_directory_lock_times_out_when_already_locked(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / knowledge_store.LOCK_FILENAME).write_text("held", encoding="utf-8")
            with self.assertRaises(TimeoutError):
                with knowledge_store._directory_lock(directory, timeout=0.01):
                    pass


class DumpFrontmatterTests(unittest.TestCase):
    def test_roundtrip_scalar_list_nested(self):
        meta = {
            "type": "knowledge",
            "title": "A title: with colon, and comma",
            "tags": ["a", "b"],
            "precondition": {"polarity": "any", "group_structure": "2-group comparison"},
        }
        text = dump_frontmatter(meta) + "\n\nbody"
        parsed, body = parse_frontmatter(text)
        self.assertEqual(parsed["type"], "knowledge")
        self.assertEqual(parsed["title"], "A title: with colon, and comma")
        self.assertEqual(parsed["tags"], ["a", "b"])
        self.assertEqual(parsed["precondition"]["group_structure"], "2-group comparison")


if __name__ == "__main__":
    unittest.main()
