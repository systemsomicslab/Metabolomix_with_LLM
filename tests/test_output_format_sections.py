"""output-format のトピック分割と、ツール単位のオンデマンド参照の検証。

トピックの正準は `mcp_core.OUTPUT_FORMAT_SECTIONS`。ここでは
  (1) 宣言した全トピックに実体ファイルがあること
  (2) B-2b のトピック集合（`dataset` があり、旧 `arf` `arf2` `pai2` が無い）
  (3) 意味論マーカーが失われていないこと
  (4) トピック resource が読めて既読管理されること
  (5) 未読トピックのツール出力にだけ誘導 1 行が付くこと
を固定する。
"""
import asyncio
import tempfile
import unittest
from pathlib import Path

from metabolomix.core import mcp_core
from metabolomix.core import session_state
from metabolomix.tools import resources as tools_resources
from tests.config_fixtures import config_env


B2_TOPICS = {"core", "dataset", "dcl", "eic", "identity", "mztab", "library", "curation",
             "pathway", "reanalysis"}

# 分割・統合で絶対に失ってはならない意味論。トピック -> そのトピックに必ず残る語。
REQUIRED_MARKERS = {
    "core": ["IsGapFilled", "34:1", "Unknown"],
    # 旧 arf / arf2 / pai2 の行の粒度・gap-fill・必須の注意（output-format-mapping.md が対応表）。
    "dataset": ["MasterAlignmentID", "IsGapFilled", "MasterPeakID", "build_pca_matrix",
                "Welch", "extract_arf2_data", "ChromatogramPeakFeature", "単一サンプル",
                "detected_mask", "feature_metadata"],
    "dcl": ["MSDec"],
    "eic": ["peak_top", "lipidmix.eic.v1", "lipidmix.eic.multi.v1"],
    "identity": ["GOSLIN", "MSI"],
    # SME（証拠）と SML（MS1 注釈）の区別は、この経路で最も誤読されやすい。
    "mztab": ["SMF_ID", "ms1_annotation", "mztab_sml", "gap-fill"],
}


def test_topics_are_b2_set():
    assert set(mcp_core.OUTPUT_FORMAT_SECTIONS) == B2_TOPICS
    for removed in ("arf", "arf2", "pai2"):
        assert removed not in mcp_core.OUTPUT_FORMAT_SECTIONS
        assert not (mcp_core.OUTPUT_FORMAT_DIR / f"{removed}.md").exists()


def test_topic_descriptions_name_only_surviving_tools():
    removed = ("dataset_handle_", "arf_", "arf2_", "pai2_", "dcl_parser", "eic_parser",
               "eic_search", "dataset_load")
    for topic, description in mcp_core.OUTPUT_FORMAT_SECTIONS.items():
        for name in removed:
            assert name not in description, f"{topic} の説明に消した名前: {name}"


class SectionFileTests(unittest.TestCase):
    def test_every_declared_section_has_a_file(self):
        self.assertTrue(mcp_core.OUTPUT_FORMAT_SECTIONS, "トピック宣言が空")
        for topic in mcp_core.OUTPUT_FORMAT_SECTIONS:
            with self.subTest(topic=topic):
                path = mcp_core.output_format_section_path(topic)
                self.assertTrue(path.is_file(), f"実体ファイルが無い: {path}")
                self.assertTrue(path.read_text(encoding="utf-8").strip())

    def test_every_section_has_a_one_line_description(self):
        for topic, description in mcp_core.OUTPUT_FORMAT_SECTIONS.items():
            with self.subTest(topic=topic):
                self.assertTrue(description.strip())
                self.assertNotIn("\n", description)

    def test_key_semantics_survive(self):
        for topic, markers in REQUIRED_MARKERS.items():
            text = mcp_core.output_format_section_path(topic).read_text(encoding="utf-8")
            for marker in markers:
                with self.subTest(topic=topic, marker=marker):
                    self.assertIn(marker, text)

    def test_core_indexes_every_other_topic(self):
        """核を読めば、どのトピックを追加で引けばよいか分かること。"""
        core = mcp_core.output_format_section_path("core").read_text(encoding="utf-8")
        for topic in mcp_core.OUTPUT_FORMAT_SECTIONS:
            if topic == "core":
                continue
            with self.subTest(topic=topic):
                self.assertIn(f"`{topic}`", core)


class TopicListsInLlmFacingTextTests(unittest.TestCase):
    """LLM に見せる案内文が、宣言済みトピックを 1 つも取りこぼさないこと。"""

    def _non_core_topics(self):
        return [t for t in mcp_core.OUTPUT_FORMAT_SECTIONS if t != "core"]

    def test_server_instructions_list_every_topic(self):
        for topic in self._non_core_topics():
            with self.subTest(topic=topic):
                self.assertIn(f"`{topic}`", mcp_core.MCP_INSTRUCTIONS)

    def test_section_template_description_lists_every_topic(self):
        templates = asyncio.run(mcp_core.mcp.list_resource_templates())
        by_uri = {str(t.uriTemplate): t for t in templates}
        description = by_uri["lipidmix://docs/output-format/{topic}"].description
        for topic in self._non_core_topics():
            with self.subTest(topic=topic):
                self.assertIn(topic, description)

    def test_static_resource_description_lists_every_topic(self):
        resources = asyncio.run(mcp_core.mcp.list_resources())
        by_uri = {str(r.uri): r for r in resources}
        description = by_uri["lipidmix://docs/output-format"].description
        for topic in self._non_core_topics():
            with self.subTest(topic=topic):
                self.assertIn(topic, description)


class SectionResourceTests(unittest.TestCase):
    def setUp(self):
        session_state.session = session_state.AnalysisSession()

    def test_section_resource_template_registered(self):
        templates = asyncio.run(mcp_core.mcp.list_resource_templates())
        uris = {str(t.uriTemplate) for t in templates}
        self.assertIn("lipidmix://docs/output-format/{topic}", uris)

    def test_reading_section_returns_body_and_marks_seen(self):
        text = tools_resources.output_format_section("eic")
        self.assertIn("peak_top", text)
        self.assertIn("eic", session_state.session.sections_seen)

    def test_unknown_topic_lists_valid_topics(self):
        with self.assertRaises(ValueError) as ctx:
            tools_resources.output_format_section("nope")
        message = str(ctx.exception)
        for topic in mcp_core.OUTPUT_FORMAT_SECTIONS:
            self.assertIn(topic, message)

    def test_removed_topics_are_unknown(self):
        for topic in ("arf", "arf2", "pai2"):
            with self.subTest(topic=topic), self.assertRaises(ValueError):
                tools_resources.output_format_section(topic)

    def test_unknown_topic_cannot_escape_the_directory(self):
        with self.assertRaises(ValueError):
            tools_resources.output_format_section("../../server")

    def test_core_resource_still_marks_output_format_seen(self):
        tools_resources.output_format_reference()
        self.assertTrue(session_state.session.output_format_seen)
        self.assertIn("core", session_state.session.sections_seen)


class SectionHintTests(unittest.TestCase):
    """ツール出力には、そのツールのトピックが未読のときだけ誘導 1 行が付く。"""

    def setUp(self):
        session_state.session = session_state.AnalysisSession()

    def test_hint_appended_for_unread_topic(self):
        out = session_state.session.maybe_prepend_caveat("BODY", topic="eic")
        self.assertIn("BODY", out)
        self.assertIn("lipidmix://docs/output-format/eic", out)

    def test_hint_suppressed_once_section_read(self):
        tools_resources.output_format_section("eic")
        out = session_state.session.maybe_prepend_caveat("BODY", topic="eic")
        self.assertNotIn("lipidmix://docs/output-format/eic", out)

    def test_hint_is_per_topic(self):
        tools_resources.output_format_section("eic")
        out = session_state.session.maybe_prepend_caveat("BODY", topic="dataset")
        self.assertIn("lipidmix://docs/output-format/dataset", out)

    def test_mode_off_disables_hint(self):
        with (tempfile.TemporaryDirectory() as tmp,
              config_env(Path(tmp), {"server.caveat_mode": "off"})):
            out = session_state.session.maybe_prepend_caveat("BODY", topic="eic")
        self.assertEqual(out, "BODY")

    def test_hint_mentions_lipid_grammar_only_for_lipid_assay(self):
        session_state.session.set_assay_kind("lipid")
        hint = session_state.session.section_hint("dataset")
        self.assertIn("脂質名文法", hint)

    def test_hint_drops_lipid_grammar_for_metabolite_assay(self):
        session_state.session.set_assay_kind("metabolite")
        hint = session_state.session.section_hint("dataset")
        self.assertIn("lipidmix://docs/output-format/dataset", hint)
        self.assertNotIn("脂質名文法", hint)


if __name__ == "__main__":
    unittest.main()
