"""クライアント非依存のエラーエンベロープ。"""
import json
import unittest

from metabolomix.core import mcp_errors


class MissingStateEnvelopeTests(unittest.TestCase):
    def test_envelope_has_the_contracted_shape(self):
        raw = mcp_errors.missing_state(
            "preprocessed_matrix",
            ["arf_preprocess"],
            "前処理後の行列がありません。先に arf_preprocess を実行してください。",
        )
        payload = json.loads(raw)
        self.assertEqual(payload["error"]["code"], "missing_state")
        self.assertEqual(payload["error"]["state"], "preprocessed_matrix")
        self.assertEqual(payload["error"]["required_tools"], ["arf_preprocess"])
        self.assertIn("arf_preprocess", payload["error"]["message"])

    def test_message_survives_verbatim(self):
        """既存の日本語文面をそのまま運ぶ（Claude Desktop での見え方を変えない）。"""
        message = "先に arf_parser を実行してARFデータとタグファイルを読み込んでください。"
        payload = json.loads(mcp_errors.missing_state("arf_dataset", ["arf_parser"], message))
        self.assertEqual(payload["error"]["message"], message)

    def test_japanese_is_not_escaped(self):
        raw = mcp_errors.missing_state("arf_dataset", ["arf_parser"], "先に arf_parser を")
        self.assertIn("先に", raw)

    def test_alternatives_are_preserved_in_order(self):
        """required_tools は OR の代替候補。順序が意味を持つ。"""
        payload = json.loads(
            mcp_errors.missing_state(
                "eic_plot", ["eic_plot_chromatograms", "eic_plot_compounds"], "先に"
            )
        )
        self.assertEqual(
            payload["error"]["required_tools"],
            ["eic_plot_chromatograms", "eic_plot_compounds"],
        )

    def test_empty_required_tools_is_rejected(self):
        """復旧の手掛かりが無いエンベロープは契約違反。無言で出さない。"""
        with self.assertRaises(ValueError):
            mcp_errors.missing_state("arf_dataset", [], "先に arf_parser を")

    def test_empty_state_is_rejected(self):
        with self.assertRaises(ValueError):
            mcp_errors.missing_state("", ["arf_parser"], "先に arf_parser を")

    def test_empty_message_is_rejected(self):
        with self.assertRaises(ValueError):
            mcp_errors.missing_state("arf_dataset", ["arf_parser"], "")


class TestDatasetBadRequest(unittest.TestCase):
    def test_code_is_registered(self):
        from metabolomix.core.mcp_errors import MZTAB_ERROR_CODES
        self.assertIn("DATASET_BAD_REQUEST", MZTAB_ERROR_CODES)

    def test_envelope_carries_details(self):
        import json
        from metabolomix.core.mcp_errors import mztab_error
        parsed = json.loads(mztab_error(
            "DATASET_BAD_REQUEST", "群サイズ不足", {"n_a": 1, "n_b": 1}))
        self.assertEqual(parsed["error"]["code"], "DATASET_BAD_REQUEST")
        self.assertEqual(parsed["error"]["message"], "群サイズ不足")
        self.assertEqual(parsed["error"]["details"], {"n_a": 1, "n_b": 1})

    def test_is_not_missing_state(self):
        """引数エラーを missing_state と混同しないことを固定する。"""
        import json
        from metabolomix.core.mcp_errors import MISSING_STATE, mztab_error
        parsed = json.loads(mztab_error("DATASET_BAD_REQUEST", "x"))
        self.assertNotEqual(parsed["error"]["code"], MISSING_STATE)
        self.assertNotIn("required_tools", parsed["error"])


if __name__ == "__main__":
    unittest.main()
