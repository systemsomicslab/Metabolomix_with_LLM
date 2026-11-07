"""意味論 caveat ガード（output-format 未 pull 時に1回だけ前置）の単体テスト。

条件付きガードの4条件（未読時のみ / 1回のみ / リソース読了で解除 / MODE=off で無効）
と、フラグがデータ切替（reset_analysis_state / load_data）で消えないことを検証する。
"""
import json
import tempfile
import unittest
from pathlib import Path

from metabolomix.core import handles, session_state
from metabolomix.tools import handle_tools, objective
from tests.config_fixtures import config_env
from tests.reanalysis_fixtures import write_mztab


class TestSemanticsCaveatGuard(unittest.TestCase):
    def setUp(self):
        session_state.session = session_state.AnalysisSession()

    def test_defaults_flags_false(self):
        self.assertFalse(session_state.session.output_format_seen)
        self.assertFalse(session_state.session.caveat_emitted)

    def test_prepends_once_then_stops(self):
        s = session_state.session
        first = s.maybe_prepend_caveat("BODY-1")
        self.assertTrue(first.startswith(session_state.SEMANTICS_CAVEAT))
        self.assertIn("BODY-1", first)
        self.assertTrue(s.caveat_emitted)
        # 2回目は前置しない（本文そのまま）
        second = s.maybe_prepend_caveat("BODY-2")
        self.assertEqual(second, "BODY-2")

    def test_resource_read_suppresses_prepend(self):
        s = session_state.session
        s.output_format_seen = True  # output-format を fetch 済み相当
        out = s.maybe_prepend_caveat("BODY")
        self.assertEqual(out, "BODY")
        self.assertFalse(s.caveat_emitted)

    def test_mode_off_disables(self):
        s = session_state.session
        with (tempfile.TemporaryDirectory() as tmp,
              config_env(Path(tmp), {"server.caveat_mode": "off"})):
            out = s.maybe_prepend_caveat("BODY")
        self.assertEqual(out, "BODY")
        self.assertFalse(s.caveat_emitted)

    def test_mode_default_is_digest(self):
        s = session_state.session
        with tempfile.TemporaryDirectory() as tmp, config_env(Path(tmp), {}):  # キーが無い
            out = s.maybe_prepend_caveat("BODY")
        self.assertTrue(out.startswith(session_state.SEMANTICS_CAVEAT))

class TestResourceReadClearsGuard(unittest.TestCase):
    """output_format_reference() を実際に呼ぶと output_format_seen が立つこと。"""

    def setUp(self):
        session_state.session = session_state.AnalysisSession()

    def test_reading_resource_sets_seen(self):
        from metabolomix.tools import resources as tools_resources
        self.assertFalse(session_state.session.output_format_seen)
        text = tools_resources.output_format_reference()
        self.assertIn("オントロジー", text)
        self.assertTrue(session_state.session.output_format_seen)


class TestAssayKindDigest(unittest.TestCase):
    """意味論ダイジェストの種別行が assay_kind で切り替わること。

    共通部（SEMANTICS_CAVEAT）はどの種別でも前置の先頭にあり、その後ろに
    種別固有の1行だけが付く。行数は種別によらず一定。
    """

    def setUp(self):
        session_state.session = session_state.AnalysisSession()

    def test_default_kind_is_unknown(self):
        self.assertEqual(session_state.session.assay_kind, "unknown")

    def test_common_block_prefixes_every_kind(self):
        for kind in session_state.ASSAY_KINDS:
            digest = session_state.assay_digest(kind)
            self.assertTrue(digest.startswith(session_state.SEMANTICS_CAVEAT))
            # 共通部＋種別行1行。増やさない。
            extra = digest[len(session_state.SEMANTICS_CAVEAT):].strip().splitlines()
            self.assertEqual(len(extra), 1)

    def test_unknown_withholds_lipid_grammar(self):
        digest = session_state.assay_digest("unknown")
        self.assertNotIn("16:0/18:1", digest)
        self.assertIn("未確定", digest)

    def test_lipid_kind_carries_lipid_grammar(self):
        digest = session_state.assay_digest("lipid")
        self.assertIn("16:0/18:1", digest)

    def test_metabolite_kind_carries_candidate_rules(self):
        digest = session_state.assay_digest("metabolite")
        self.assertIn("アダクト", digest)
        self.assertNotIn("16:0/18:1", digest)

    def test_set_assay_kind_rejects_unknown_value(self):
        with self.assertRaises(ValueError):
            session_state.session.set_assay_kind("proteomics")

    def test_set_assay_kind_normalizes_case(self):
        self.assertEqual(session_state.session.set_assay_kind("LIPID"), "lipid")
        self.assertEqual(session_state.session.assay_kind, "lipid")

    def test_reemits_once_when_kind_becomes_known(self):
        s = session_state.session
        first = s.maybe_prepend_caveat("BODY-1")
        self.assertIn("未確定", first)
        # 種別が確定したら、確定後の規則を1回だけ届ける
        s.set_assay_kind("metabolite")
        second = s.maybe_prepend_caveat("BODY-2")
        self.assertIn("アダクト", second)
        # 同じ種別のままなら再注入しない
        third = s.maybe_prepend_caveat("BODY-3")
        self.assertEqual(third, "BODY-3")

    def test_resource_read_suppresses_reemission(self):
        s = session_state.session
        s.maybe_prepend_caveat("BODY-1")
        s.output_format_seen = True
        s.set_assay_kind("lipid")
        self.assertEqual(s.maybe_prepend_caveat("BODY-2"), "BODY-2")

class TestDigestTravelsWithLoadDataset(unittest.TestCase):
    """`load_dataset` の戻り値（JSON）の `caveat` に意味論ダイジェストが載る（B-2b Task 10）。

    旧パーサが消えて `maybe_prepend_caveat` を呼ぶ口が無くなったので、入口の `load_dataset` が
    JSON の最上位の `caveat` 文字列として運ぶ。1 回だけ・種別ごと・`caveat_mode=off` で止まる。
    """

    def setUp(self):
        session_state.session = session_state.AnalysisSession()
        handles.registry = handles.HandleRegistry()
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.mztab = write_mztab(Path(self._tmp.name) / "AlignResult-1.mzTab")

    def _load(self) -> dict:
        return json.loads(handle_tools.load_dataset(mztab_path=str(self.mztab), library="none"))

    def test_first_call_carries_the_digest_and_the_second_does_not(self):
        first = self._load()
        self.assertTrue(first["caveat"].startswith(session_state.SEMANTICS_CAVEAT))
        self.assertIn("未確定", first["caveat"])        # 種別はまだ unknown
        self.assertNotIn("caveat", self._load())

    def test_mode_off_withholds_the_digest(self):
        with config_env(Path(self._tmp.name), {"server.caveat_mode": "off"}):
            self.assertNotIn("caveat", self._load())

    def test_output_format_resource_read_withholds_the_digest(self):
        session_state.session.output_format_seen = True
        self.assertNotIn("caveat", self._load())

    def test_settling_the_kind_emits_the_digest_once_more(self):
        self._load()                                    # unknown 用が出る
        out = objective.record_objective(
            "exp-1", dataset="POS", polarity="positive", groups=["HC", "RA"],
            comparison="HC vs RA", sub_questions=["問い"], assay_kind="lipid")
        self.assertTrue(out.startswith(session_state.SEMANTICS_CAVEAT))
        self.assertIn("16:0/18:1", out)                 # 脂質の規則が届く
        self.assertIn("objective を記録", out)
        # 同じ種別のままなら、続く入口にも更新にも付かない。
        self.assertNotIn("caveat", self._load())
        again = objective.record_objective("exp-1", status="active", assay_kind="lipid")
        self.assertFalse(again.startswith(session_state.SEMANTICS_CAVEAT))

    def test_record_objective_without_assay_kind_adds_no_digest(self):
        self._load()
        out = objective.record_objective(
            "exp-2", dataset="POS", polarity="positive", groups=["HC", "RA"],
            comparison="HC vs RA", sub_questions=["問い"])
        self.assertFalse(out.startswith(session_state.SEMANTICS_CAVEAT))


if __name__ == "__main__":
    unittest.main()
