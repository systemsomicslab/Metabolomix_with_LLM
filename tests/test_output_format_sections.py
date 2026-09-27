"""output-format のトピック分割と、ツール単位のオンデマンド参照の検証。

一枚岩の docs/output_format.md（約700行）を毎回全文 pull させると、解釈に不要な
節まで文脈を食う。パーサ/ツール単位のトピックへ割り、その出力を解釈する直前に
該当トピックだけを引かせるのが本改修の狙い。ここでは
  (1) 宣言した全トピックに実体ファイルがあること
  (2) 分割で意味論マーカーが失われていないこと
  (3) トピック resource が読めて既読管理されること
  (4) 未読トピックのツール出力にだけ誘導1行が付くこと
を固定する。
"""
import asyncio
import os
import unittest
from unittest import mock

from lipidmix.core import mcp_core
from lipidmix.core import session_state
from lipidmix.tools import resources as tools_resources


# 分割で絶対に失ってはならない意味論。トピック -> そのトピックに必ず残る語。
REQUIRED_MARKERS = {
    "core": ["IsGapFilled", "34:1", "Unknown"],
    "arf": ["MasterAlignmentID", "build_pca_matrix", "Welch"],
    "arf2": ["extract_arf2_data"],
    "pai2": ["ChromatogramPeakFeature", "単一サンプル"],
    "dcl": ["MSDec"],
    "eic": ["peak_top", "lipidmix.eic.v1", "lipidmix.eic.multi.v1"],
    "identity": ["GOSLIN", "MSI"],
    # SME（証拠）と SML（MS1 注釈）の区別は、この経路で最も誤読されやすい。
    "mztab": ["SMF_ID", "ms1_annotation", "mztab_sml", "gap-fill"],
}


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

    def test_key_semantics_survive_the_split(self):
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
                self.assertIn(topic, core)


class TopicListsInLlmFacingTextTests(unittest.TestCase):
    """LLM に見せる案内文が、宣言済みトピックを 1 つも取りこぼさないこと。

    正準は `OUTPUT_FORMAT_SECTIONS`。サーバ共通指示やリソース説明の一覧が
    それより短いと、案内に載らないトピック（過去には library・mztab）は
    LLM から存在しないのと同じになる。
    """

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

    def test_unknown_topic_cannot_escape_the_directory(self):
        with self.assertRaises(ValueError):
            tools_resources.output_format_section("../../server")

    def test_core_resource_still_marks_output_format_seen(self):
        tools_resources.output_format_reference()
        self.assertTrue(session_state.session.output_format_seen)
        self.assertIn("core", session_state.session.sections_seen)


class SectionHintTests(unittest.TestCase):
    """ツール出力には、そのツールのトピックが未読のときだけ誘導1行が付く。"""

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
        out = session_state.session.maybe_prepend_caveat("BODY", topic="arf")
        self.assertIn("lipidmix://docs/output-format/arf", out)

    def test_mode_off_disables_hint(self):
        with mock.patch.dict(os.environ, {"LIPIDMIX_CAVEAT_MODE": "off"}):
            out = session_state.session.maybe_prepend_caveat("BODY", topic="eic")
        self.assertEqual(out, "BODY")

    def test_topic_omitted_keeps_previous_behaviour(self):
        out = session_state.session.maybe_prepend_caveat("BODY")
        self.assertTrue(out.startswith(session_state.SEMANTICS_CAVEAT))
        self.assertIn("BODY", out)

    def test_hint_mentions_lipid_grammar_only_for_lipid_assay(self):
        session_state.session.set_assay_kind("lipid")
        hint = session_state.session.section_hint("arf")
        self.assertIn("脂質名文法", hint)

    def test_hint_drops_lipid_grammar_for_metabolite_assay(self):
        session_state.session.set_assay_kind("metabolite")
        hint = session_state.session.section_hint("arf")
        self.assertIn("lipidmix://docs/output-format/arf", hint)
        self.assertNotIn("脂質名文法", hint)


if __name__ == "__main__":
    unittest.main()
