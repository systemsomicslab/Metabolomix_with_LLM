"""状態不足の箇所が missing_state エンベロープを返すことの契約テスト。

message は既存の日本語文面を保持する（エンベロープを解釈しないクライアントでも
LLM が読む内容が変わらないこと）。
"""
import json
import unittest

import server
from metabolomix.core import session_state


def envelope(raw):
    """本文がエンベロープならその error 部を返す。違えば None。"""
    try:
        payload = json.loads(raw)
    except (ValueError, TypeError):
        return None
    if not isinstance(payload, dict):
        return None
    error = payload.get("error")
    if isinstance(error, dict) and error.get("code") == "missing_state":
        return error
    return None


class MissingStateContractTests(unittest.TestCase):
    """空セッションで各ツールを呼び、エンベロープが返ることを確認する。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def assert_missing(self, raw, state, required_tools, message_contains):
        error = envelope(raw)
        self.assertIsNotNone(error, f"エンベロープでない: {raw[:200]}")
        self.assertEqual(error["state"], state)
        self.assertEqual(error["required_tools"], required_tools)
        self.assertIn(message_contains, error["message"])

    def assert_no_dataset(self, raw):
        self.assert_missing(raw, "dataset", ["load_dataset"], "load_dataset")

    def test_dataset_status(self):
        # 旧 arf_list_tags / arf_list_classes / arf_list_sample_roles（B-2b Task 8 で消えた）の後継
        self.assert_no_dataset(server.dataset_status())
        self.assert_no_dataset(server.dataset_status(view="samples"))
        self.assert_no_dataset(server.dataset_status(view="tags"))

    def test_dataset_exclude(self):
        self.assert_no_dataset(server.dataset_exclude(samples=["x"]))

    def test_dataset_preprocess(self):
        self.assert_no_dataset(server.dataset_preprocess())

    def test_dataset_pca(self):
        self.assert_no_dataset(server.dataset_pca())

    def test_dataset_compare(self):
        self.assert_no_dataset(server.dataset_compare(reference_group="a", test_group="b"))

    def test_dataset_feature_search(self):
        self.assert_no_dataset(server.dataset_feature_search(sort="intensity"))

    def test_save_figure_without_dataset(self):
        # 図の置き場はデータセット（ds.figures）。データセットが無ければ load_dataset を案内する。
        self.assert_no_dataset(server.save_figure("pca", analysis_id="x"))

    def test_dataset_plot_without_dataset(self):
        """どの kind でも、素の戻り値としてエンベロープを返す（旧 arf_plot_* / plot_pca_loadings /
        eic_plot_* の契約を引き継ぐ）。"""
        for kind, kwargs in [
            ("volcano", {}),
            ("pca", {}),
            ("loadings", {"output": "payload"}),
            ("group_intensity", {"items": ["PG"], "groups": ["x"], "output": "payload"}),
            ("species", {"items": ["PG"], "groups": ["x"], "output": "payload"}),
            ("species_pca", {"items": ["PG"], "groups": ["x"], "output": "payload"}),
            ("eic", {"feature_id": "1"}),
        ]:
            with self.subTest(kind=kind):
                self.assert_no_dataset(server.dataset_plot(kind=kind, **kwargs))


if __name__ == "__main__":
    unittest.main()
