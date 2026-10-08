"""状態不足15箇所が missing_state エンベロープを返すことの契約テスト。

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

    def test_arf_list_tags(self):
        self.assert_missing(
            server.arf_list_tags(), "arf_dataset", ["arf_parser", "load_dataset"], "タグファイル"
        )

    def test_arf_list_classes(self):
        self.assert_missing(
            server.arf_list_classes(), "arf_dataset", ["arf_parser", "load_dataset"], "ARFデータ"
        )

    def test_arf_list_sample_roles(self):
        self.assert_missing(
            server.arf_list_sample_roles(), "arf_dataset", ["arf_parser", "load_dataset"], "arf_parser"
        )

    def test_arf_exclude(self):
        self.assert_missing(
            server.arf_exclude(exclude_samples=["x"]), "arf_dataset", ["arf_parser", "load_dataset"], "arf_parser"
        )

    def test_arf_preprocess(self):
        self.assert_missing(
            server.arf_preprocess(), "arf_dataset", ["arf_parser", "load_dataset"], "arf_parser"
        )

    def test_arf_pca_preprocessed(self):
        self.assert_missing(
            server.arf_pca_preprocessed(),
            "preprocessed_matrix",
            ["arf_preprocess"],
            "前処理後の行列がありません",
        )

    def test_arf_differential(self):
        self.assert_missing(
            server.arf_differential(group_a="a", group_b="b"),
            "preprocessed_matrix",
            ["arf_preprocess"],
            "arf_preprocess",
        )

    def test_arf_plot_volcano(self):
        """他の12ツールと同じく、素の戻り値としてエンベロープを返す。

        以前は outputSchema 導出のため戻り値型が構造化に固定されており、失敗時に
        str を返せず例外の本文にエンベロープを載せていた。structured_output=False に
        した今はその制約が無いので、例外経由の特例を廃して契約を揃えてある。
        """
        self.assert_missing(
            server.arf_plot_volcano(),
            "differential_result",
            ["arf_differential"],
            "arf_differential",
        )

    def test_pai2_inspect_peak(self):
        self.assert_missing(
            server.pai2_inspect_peak(peak_name="PC 34:1"),
            "pai2_dataset",
            ["pai2_parser"],
            "pai2_parser",
        )

    def test_verify_peak_annotation(self):
        self.assert_missing(
            server.verify_peak_annotation(peak_name="PC 34:1"),
            "pai2_dataset",
            ["pai2_parser"],
            "pai2_parser",
        )

    def test_save_figure_pca(self):
        self.assert_missing(
            server.save_figure("pca", analysis_id="x"),
            "pca_result",
            # mzTab-M 経路（dataset_pca）も PCA の生産者なのでリプレイ候補に入る
            ["arf_parser", "arf_pca_preprocessed", "load_dataset", "dataset_pca"],
            "PCA",
        )

    def test_save_figure_volcano(self):
        self.assert_missing(
            server.save_figure("volcano", analysis_id="x"),
            "differential_result",
            # 図の保存は mzTab-M 経路（dataset_differential）からもできる。
            # arf_plot_volcano は ARF 専用なので候補は増やさない。
            ["arf_differential", "dataset_differential"],
            "arf_differential",
        )

    def test_arf_plot_group_intensity(self):
        """ARF 未読み込みの空セッションでは、素の戻り値としてエンベロープを返す。"""
        self.assert_missing(
            server.arf_plot_group_intensity(items=["PG"], groups=["x"], output="payload"),
            "arf_dataset",
            ["arf_parser", "load_dataset"],
            "load_dataset",
        )

    def test_arf_plot_species(self):
        self.assert_missing(
            server.arf_plot_species(items=["PG"], groups=["x"], output="payload"),
            "arf_dataset",
            ["arf_parser", "load_dataset"],
            "load_dataset",
        )

    def test_save_figure_group_intensity(self):
        self.assert_missing(
            server.save_figure("group_intensity", analysis_id="x"),
            "group_intensity_plot",
            ["arf_plot_group_intensity"],
            "arf_plot_group_intensity",
        )

    def test_save_figure_eic_offers_both_producers(self):
        """eic_plot は2つのツールのどちらでも作れる（OR の代替候補）。"""
        self.assert_missing(
            server.save_figure("eic", analysis_id="x"),
            "eic_plot",
            ["eic_plot_chromatograms", "eic_plot_compounds"],
            "eic_plot_chromatograms",
        )


class SampleSearchIsNotAMissingStateTests(unittest.TestCase):
    """sample_search は ARF 未ロードでも動くのが要件。エンベロープ化してはいけない。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_sample_search_failure_is_not_an_envelope(self):
        import tempfile

        with tempfile.TemporaryDirectory() as empty:
            raw = server.sample_search(directory=empty)
        self.assertIsNone(
            envelope(raw),
            "sample_search をエンベロープ化すると、ARF を必要としない検索のために "
            "arf_parser の自動リプレイが走ってしまう",
        )


if __name__ == "__main__":
    unittest.main()
