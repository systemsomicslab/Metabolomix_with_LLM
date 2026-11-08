"""クライアント非依存のエラーエンベロープ。"""
import json
import unittest

from metabolomix.core import mcp_errors


class MissingStateEnvelopeTests(unittest.TestCase):
    def test_envelope_has_the_contracted_shape(self):
        raw = mcp_errors.missing_state(
            "preprocessed_matrix",
            ["dataset_preprocess"],
            "前処理後の行列がありません。先に dataset_preprocess を実行してください。",
        )
        payload = json.loads(raw)
        self.assertEqual(payload["error"]["code"], "missing_state")
        self.assertEqual(payload["error"]["state"], "preprocessed_matrix")
        self.assertEqual(payload["error"]["required_tools"], ["dataset_preprocess"])
        self.assertIn("dataset_preprocess", payload["error"]["message"])

    def test_message_survives_verbatim(self):
        """既存の日本語文面をそのまま運ぶ（Claude Desktop での見え方を変えない）。"""
        message = "先に load_dataset を実行してデータセットを読み込んでください。"
        payload = json.loads(mcp_errors.missing_state("dataset", ["load_dataset"], message))
        self.assertEqual(payload["error"]["message"], message)

    def test_japanese_is_not_escaped(self):
        raw = mcp_errors.missing_state("dataset", ["load_dataset"], "先に load_dataset を")
        self.assertIn("先に", raw)

    def test_alternatives_are_preserved_in_order(self):
        """required_tools は OR の代替候補。順序が意味を持つ。"""
        payload = json.loads(
            mcp_errors.missing_state(
                "analysis_matrix", ["dataset_preprocess", "dataset_build_matrix"], "先に"
            )
        )
        self.assertEqual(
            payload["error"]["required_tools"],
            ["dataset_preprocess", "dataset_build_matrix"],
        )

    def test_empty_required_tools_is_rejected(self):
        """復旧の手掛かりが無いエンベロープは契約違反。無言で出さない。"""
        with self.assertRaises(ValueError):
            mcp_errors.missing_state("dataset", [], "先に load_dataset を")

    def test_empty_state_is_rejected(self):
        with self.assertRaises(ValueError):
            mcp_errors.missing_state("", ["load_dataset"], "先に load_dataset を")

    def test_empty_message_is_rejected(self):
        with self.assertRaises(ValueError):
            mcp_errors.missing_state("dataset", ["load_dataset"], "")


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


class NeedsInputEnvelopeTests(unittest.TestCase):
    """利用者に尋ねてから呼び直す封筒。エラーではない（`error` キーを持たない）。"""

    OPTIONS = [{"value": "on", "label": "補正する"}, {"value": "off", "label": "補正しない"}]

    def test_needs_input_envelope_shape(self):
        raw = mcp_errors.needs_input(
            "preprocess.drift.mode", "内部標準で換算済みの値に QC のドリフト補正を重ねるか",
            self.OPTIONS, argument="drift_mode",
            retry_tools=["dataset_preprocess", "dataset_compare"])
        payload = json.loads(raw)
        self.assertEqual(set(payload), {"needs_input"})
        self.assertNotIn("error", payload)
        body = payload["needs_input"]
        self.assertEqual(body["field"], "preprocess.drift.mode")
        self.assertIn("ドリフト", body["reason"])
        self.assertEqual(body["options"], self.OPTIONS)
        self.assertEqual(body["argument"], "drift_mode")
        self.assertEqual(body["retry_tools"], ["dataset_preprocess", "dataset_compare"])
        self.assertIn("補正する", raw)   # 日本語をエスケープしない

    def test_argument_is_omitted_when_not_given(self):
        """spec の欄で答える経路（dataset_analysis_run の v2）は argument を持たない。"""
        body = json.loads(mcp_errors.needs_input(
            "preprocess.drift.mode", "理由", self.OPTIONS,
            retry_tools=["dataset_analysis_run"]))["needs_input"]
        self.assertNotIn("argument", body)
        self.assertEqual(body["field"], "preprocess.drift.mode")

    def test_required_parts_are_validated(self):
        ok = dict(field="f", reason="r", options=self.OPTIONS, retry_tools=["t"])
        for key, bad in (("field", ""), ("reason", ""), ("options", []),
                         ("options", [{"value": "on"}]), ("retry_tools", [])):
            with self.subTest(key=key, bad=bad):
                kwargs = {**ok, key: bad}
                with self.assertRaises(ValueError):
                    mcp_errors.needs_input(kwargs.pop("field"), kwargs.pop("reason"),
                                           kwargs.pop("options"), **kwargs)


if __name__ == "__main__":
    unittest.main()
