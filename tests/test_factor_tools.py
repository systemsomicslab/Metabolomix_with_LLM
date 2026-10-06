"""MCP ツールレベルの因子トークン挙動（sample_factors 純関数は test_sample_factors.py）。"""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import server
from metabolomix.core import session_state
from tests.test_server_class_filter import (
    fake_build_pca_matrix,
    fake_extract_peak_properties,
    fake_run_pca,
)


class DifferentialByNameTokenTests(unittest.TestCase):
    """Class ID は処置だけ、時点はサンプル名にしかない構成で時点を揃えた2群比較。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()
        names = [
            "20220902_RAW_ILG_6h_1_NEG", "20220902_RAW_ILG_6h_2_NEG",
            "20220902_RAW_ILG_6h_3_NEG",
            "20220901_RAW_control_6h_1_NEG", "20220901_RAW_control_6h_2_NEG",
            "20220901_RAW_control_6h_3_NEG",
            "20220902_RAW_ILG_0h_1_NEG", "20220901_RAW_control_0h_1_NEG",
        ]
        session_state.session.arf.feature_matrix = np.array([
            [50.0, 5.0], [52.0, 5.1], [48.0, 4.9],
            [10.0, 5.0], [11.0, 5.2], [9.5, 4.8],
            [30.0, 5.0], [30.0, 5.0],
        ])
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["Spot_0_height", "Spot_1_height"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("ILG" if "ILG" in n else "control"),
                "role": "sample", "batch": "d1"}
            for n in names
        }

    def test_time_matched_two_group_comparison(self):
        out = json.loads(server.arf_differential(group_a="ILG_6h", group_b="control_6h"))
        self.assertEqual(out["status"], "success")
        self.assertEqual(out["n_a"], 3)
        self.assertEqual(out["n_b"], 3)
        self.assertEqual(out["summary"]["n_significant"], 1)

    def test_payload_lists_resolved_samples(self):
        out = json.loads(server.arf_differential(group_a="ILG_6h", group_b="control_6h"))
        self.assertEqual(out["resolved_samples"]["group_a"], [
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_RAW_ILG_6h_2_NEG",
            "20220902_RAW_ILG_6h_3_NEG",
        ])
        self.assertTrue(any("比較サンプル" in c for c in out["caveats"]))

    def test_0h_samples_are_not_pooled_in(self):
        out = json.loads(server.arf_differential(group_a="ILG_6h", group_b="control_6h"))
        joined = " ".join(out["resolved_samples"]["group_a"] + out["resolved_samples"]["group_b"])
        self.assertNotIn("_0h_", joined)

    def test_treatment_only_spec_still_pools_all_timepoints(self):
        out = json.loads(server.arf_differential(group_a="ILG", group_b="control"))
        self.assertEqual(out["n_a"], 4)
        self.assertEqual(out["n_b"], 4)

    def test_blank_group_spec_returns_structured_error(self):
        """トークンが空集合になる群指定（空文字/空白/アンダースコアのみ）。

        expand_sample_specs はそういう spec を読み飛ばすため matches に現れない。
        _pool_group_labels が KeyError で MCP ツールの外まで抜けず、他の失敗経路と
        同じ構造化エラー JSON になることを固定する（arf_exclude 側と同じ形）。
        """
        for spec in ("", "   ", "___"):
            with self.subTest(spec=spec):
                out = json.loads(
                    server.arf_differential(group_a=spec, group_b="control_6h"))
                self.assertEqual(out["status"], "error")
                self.assertIn("因子トークン", out["message"])

    def test_same_class_id_groups_are_split_by_name_token(self):
        """同一 Class ID（ILG）内で時点だけが違う2群。サンプル集合は交差しない。"""
        out = json.loads(server.arf_differential(group_a="ILG_6h", group_b="ILG_0h"))
        self.assertEqual(out["status"], "success")
        self.assertEqual(out["n_a"], 3)
        self.assertEqual(out["n_b"], 1)
        a = set(out["resolved_samples"]["group_a"])
        b = set(out["resolved_samples"]["group_b"])
        self.assertEqual(a & b, set())
        self.assertEqual(b, {"20220902_RAW_ILG_0h_1_NEG"})

    def test_shared_class_id_is_disclosed_as_not_group_defining(self):
        """resolved_class_ids は両群とも ['ILG']。縮退比較と誤読されないよう開示する。"""
        out = json.loads(server.arf_differential(group_a="ILG_6h", group_b="ILG_0h"))
        self.assertEqual(out["resolved_class_ids"]["group_a"], ["ILG"])
        self.assertEqual(out["resolved_class_ids"]["group_b"], ["ILG"])
        self.assertTrue(
            any("Class ID では区別されず" in c for c in out["caveats"]),
            out["caveats"])

    def test_pool_caveat_fires_on_multi_sample_groups(self):
        """プール caveat は Class ID 数ではなくサンプル数で発火する。

        ILG_6h / control_6h はどちらも Class ID 1個・サンプル3件。旧条件
        （len(resolved[...]) > 1）では発火せず、プールしている事実が隠れていた。
        """
        out = json.loads(server.arf_differential(group_a="ILG_6h", group_b="control_6h"))
        pooled = [c for c in out["caveats"] if "プール群として解決" in c]
        self.assertEqual(len(pooled), 1, out["caveats"])
        self.assertIn("3 サンプル", pooled[0])


class _ParserFakeArfState:
    def __init__(self, names):
        self.current_file_path = None
        self.features = [{
            "MasterAlignmentID": 10,
            "AlignedPeakProperties": [[i, n, 100.0 + i] for i, n in enumerate(names)],
        }]
        self.filtered_features = None
        self.pca_result = None
        self.tag_index = {}
        self.class_index = None
        self.excluded_samples = set()
        self.excluded_spots = set()
        self.feature_matrix = None
        self.last_pca_plot = None
        self.last_differential = None
        self.sample_meta = {}
        self.preprocessing_recipe = {}
        self.pp_sample_names = None
        self.pp_feature_names = None

    def load_data(self, file_path, tag_directory=None):
        self.current_file_path = file_path
        return self.features


class _ParserFakeSession:
    """test_server_class_filter.FakeSession と同型だが、サンプル名を差し替えられる版。"""

    def __init__(self, names):
        self.arf = _ParserFakeArfState(names)

    def maybe_prepend_caveat(self, text, topic=None):
        return text


class ArfParserFactorTests(unittest.TestCase):
    def setUp(self):
        self.session = _ParserFakeSession([
            "20220901_RAW_control_0h_1_NEG",
            "20220901_RAW_control_6h_1_NEG",
            "20220902_RAW_ILG_0h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
            "20220901_QC_RAW_NEG_1",
        ])
        self.patches = [
            patch.object(session_state, "session", self.session),
            patch.object(server.arf_reader, "extract_peak_properties", fake_extract_peak_properties),
            patch.object(server.arf_reader, "build_pca_matrix", fake_build_pca_matrix),
            patch.object(server.arf_reader, "run_pca", fake_run_pca),
            patch.object(server.arf_reader, "get_pca_loading_features", return_value=[]),
        ]
        for item in self.patches:
            item.start()

    def tearDown(self):
        for item in reversed(self.patches):
            item.stop()

    def _run(self, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            arf_path = Path(tmp) / "test.arf"
            arf_path.touch()
            return server.arf_parser(str(arf_path), **kwargs)

    def test_group_factors_produce_cross_product_counts(self):
        result = self._run(group_factors=[["control", "ILG"], ["0h", "6h"]])
        self.assertIn("control|0h=1", result)
        self.assertIn("ILG|6h=1", result)
        self.assertIn("qc=1", result)

    def test_class_ids_filter_by_name_only_token(self):
        result = self._run(class_ids=["6h"])
        self.assertIn("Class IDフィルタ**: `6h`", result)
        rows = self.session.arf.filtered_features[0]["AlignedPeakProperties"]
        self.assertEqual(
            [row[1] for row in rows],
            ["20220901_RAW_control_6h_1_NEG", "20220902_RAW_ILG_6h_1_NEG"],
        )

    def test_qc_excluded_by_role_is_disclosed(self):
        result = self._run(class_ids=["raw"])
        self.assertIn("role により除外", result)
        self.assertIn("20220901_QC_RAW_NEG_1", result)

    def test_include_roles_brings_qc_back(self):
        self._run(class_ids=["raw"], include_roles=["sample", "qc"])
        rows = self.session.arf.filtered_features[0]["AlignedPeakProperties"]
        self.assertIn("20220901_QC_RAW_NEG_1", [row[1] for row in rows])

    def test_summary_discloses_selected_samples_not_just_class_ids(self):
        """class_ids=["6h"] は 6h のサンプルだけを選ぶ。Class ID だけを見せると
        「ILG と control の全サンプル」と誤読されるため、件数と実サンプル名を出す。
        """
        result = self._run(class_ids=["6h"])
        self.assertIn("選択サンプル**: 2 件", result)
        self.assertIn("20220901_RAW_control_6h_1_NEG", result)
        self.assertNotIn("展開先", result)

    def test_list_classes_vocabulary_is_not_narrowed_by_a_prior_filter(self):
        """絞り込む arf_parser の後でも、発見用の語彙はデータセット全体を見る。

        arf_list_classes は「何で絞れるか」の入口。フィルタ後の特徴量から語彙を
        作ると 0h/control が消え、ユーザーは「そんなサンプルは無い」と誤る。
        """
        self._run(class_ids=["6h"])
        payload = json.loads(server.arf_list_classes())
        tokens = payload["sample_token_vocabulary"]["tokens"]
        self.assertIn("0h", tokens)
        self.assertEqual(tokens["6h"]["samples"], 2)
        self.assertEqual(tokens["0h"]["samples"], 2)


class ArfPcaPreprocessedFactorTests(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        names = ["20220901_RAW_control_0h_1_NEG", "20220902_RAW_ILG_6h_1_NEG"]
        session_state.session.arf.feature_matrix = np.array([[1.0, 2.0], [3.0, 4.0]])
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["Spot_0_height", "Spot_1_height"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.features = []

    def test_group_factors_label_the_plot(self):
        with patch.object(server.arf_reader, "run_pca", fake_run_pca), \
             patch.object(server.arf_reader, "get_pca_loading_features", return_value=[]):
            result = server.arf_pca_preprocessed(
                group_factors=[["control", "ILG"], ["0h", "6h"]])
        self.assertIn("control|0h=1", result)
        self.assertIn("ILG|6h=1", result)


class ArfListClassesVocabularyTests(unittest.TestCase):
    def setUp(self):
        self.session = _ParserFakeSession([
            "20220901_RAW_control_6h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
            "20220901_QC_RAW_NEG_1",
        ])
        self.session.arf.current_file_path = "test.arf"
        self.patcher = patch.object(session_state, "session", self.session)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_reports_sample_token_vocabulary(self):
        payload = json.loads(server.arf_list_classes())
        tokens = payload["sample_token_vocabulary"]["tokens"]
        self.assertEqual(tokens["6h"]["samples"], 2)
        self.assertEqual(tokens["ilg"]["samples"], 1)

    def test_vocabulary_breaks_down_by_role(self):
        payload = json.loads(server.arf_list_classes())
        tokens = payload["sample_token_vocabulary"]["tokens"]
        self.assertEqual(tokens["raw"]["roles"], {"qc": 1, "sample": 2})

    def test_works_without_mddata(self):
        payload = json.loads(server.arf_list_classes())
        self.assertIsNone(payload["mddata_path"])
        self.assertEqual(payload["class_counts"], {})


def _exclude_row(file_id, name, height):
    """exclusions.roster が file_name を拾える現実的な生 row（len>25・data[18] 数値）。"""
    row = [0] * 26
    row[0] = file_id
    row[1] = name
    row[2] = file_id
    row[18] = float(height)
    return row


class ArfExcludeSpecTests(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        names = [
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_RAW_ILG_6h_2_NEG",
            "20220902_RAW_ILG_0h_1_NEG",
            "20220901_RAW_control_6h_1_NEG",
        ]
        session_state.session.arf.filtered_features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [_exclude_row(i, n, 10 + i) for i, n in enumerate(names)],
        }]

    def test_spec_excludes_all_matching_samples(self):
        out = json.loads(server.arf_exclude(exclude_samples=["ILG_6h"]))
        self.assertEqual(out["status"], "success")
        self.assertEqual(sorted(out["excluded_samples"]), [
            "20220902_RAW_ILG_6h_1_NEG", "20220902_RAW_ILG_6h_2_NEG",
        ])
        self.assertEqual(out["samples_after"], 2)

    def test_payload_discloses_spec_resolution(self):
        out = json.loads(server.arf_exclude(exclude_samples=["ILG_6h"]))
        self.assertEqual(out["resolved_samples"]["ILG_6h"], [
            "20220902_RAW_ILG_6h_1_NEG", "20220902_RAW_ILG_6h_2_NEG",
        ])

    def test_remove_mode_accepts_specs(self):
        server.arf_exclude(exclude_samples=["ILG_6h"])
        out = json.loads(server.arf_exclude(exclude_samples=["ILG_6h"], mode="remove"))
        self.assertEqual(out["excluded_samples"], [])
        self.assertEqual(session_state.session.arf.excluded_samples, set())

    def test_exact_name_still_wins(self):
        out = json.loads(server.arf_exclude(exclude_samples=["20220902_RAW_ILG_6h_1_NEG"]))
        self.assertEqual(out["excluded_samples"], ["20220902_RAW_ILG_6h_1_NEG"])

    def test_unmatched_spec_is_reported_not_raised(self):
        out = json.loads(server.arf_exclude(exclude_samples=["24h"]))
        self.assertEqual(out["status"], "success")
        self.assertIn("24h", out["unmatched_samples"])
        self.assertEqual(session_state.session.arf.excluded_samples, set())

    def test_empty_or_underscore_only_specs_are_unmatched_not_raised(self):
        """split_tokens が空集合を返す spec（空文字/空白/アンダースコアのみ）は
        expand_sample_specs 内で continue され matches に入らない。resolved 側が
        それを KeyError で落とさず unmatched_samples 行きにできることを確認する。
        """
        out = json.loads(server.arf_exclude(exclude_samples=["", "   ", "___"]))
        self.assertEqual(out["status"], "success")
        for spec in ("", "   ", "___"):
            self.assertIn(spec, out["unmatched_samples"])
        self.assertEqual(out["excluded_samples"], [])
        self.assertEqual(session_state.session.arf.excluded_samples, set())


class ArfExcludeQcRoleSpecTests(unittest.TestCase):
    """include_roles=None（tools_arf._resolve_exclude_specs）の意図: 除外指定は
    group 形成と異なり QC/blank ロールのサンプルも因子トークン spec で狙えてよい。
    """

    def setUp(self):
        session_state.session = server.AnalysisSession()
        names = [
            "20220902_RAW_ILG_6h_1_NEG",
            "20220901_QC_RAW_NEG_1",
        ]
        session_state.session.arf.filtered_features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [_exclude_row(i, n, 10 + i) for i, n in enumerate(names)],
        }]

    def test_factor_token_spec_excludes_qc_sample(self):
        out = json.loads(server.arf_exclude(exclude_samples=["QC"]))
        self.assertEqual(out["status"], "success")
        self.assertEqual(out["excluded_samples"], ["20220901_QC_RAW_NEG_1"])
        self.assertEqual(out["resolved_samples"]["QC"], ["20220901_QC_RAW_NEG_1"])


if __name__ == "__main__":
    unittest.main()
