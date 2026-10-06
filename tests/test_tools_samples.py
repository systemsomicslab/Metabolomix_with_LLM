import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from metabolomix.core import mcp_core
import server
from metabolomix.core import session_state


SAMPLES = [
    "20220901_RAW_control_6h_1_NEG",
    "20220902_RAW_ILG_6h_1_NEG",
    "20220902_RAW_ILG_6h_2_NEG",
    "20220901_QC_RAW_NEG_1",
]


def make_dir(tmp: str, timestamps=("202605151012",)) -> Path:
    """サンプルごとの .pai2 / .dcl を、指定バッチ分だけ作る。"""
    directory = Path(tmp)
    for stamp in timestamps:
        for name in SAMPLES:
            for suffix in (".pai2", ".dcl"):
                (directory / f"{name}_{stamp}{suffix}").touch()
    return directory


class SampleSearchFromDirectoryTests(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_vocabulary_mode_lists_all_samples_and_tokens(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(server.sample_search(directory=str(make_dir(tmp))))
        self.assertEqual(payload["status"], "success")
        self.assertEqual(payload["total_samples"], 4)
        self.assertEqual(payload["token_vocabulary"]["tokens"]["6h"]["samples"], 3)

    def test_spec_mode_returns_matching_samples_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(
                server.sample_search(specs=["ILG_6h"], directory=str(make_dir(tmp))))
        self.assertEqual(payload["matched"], 2)
        self.assertEqual(
            [s["name"] for s in payload["samples"]],
            ["20220902_RAW_ILG_6h_1_NEG", "20220902_RAW_ILG_6h_2_NEG"],
        )
        self.assertNotIn("token_vocabulary", payload)

    def test_returns_real_file_paths_per_extension(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = make_dir(tmp)
            payload = json.loads(
                server.sample_search(specs=["ILG_6h_1"], directory=str(directory)))
            files = payload["samples"][0]["files"]
            self.assertTrue(files[".pai2"].endswith("20220902_RAW_ILG_6h_1_NEG_202605151012.pai2"))
            self.assertTrue(Path(files[".dcl"]).is_file())

    def test_selects_latest_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = make_dir(tmp, timestamps=("202601010000", "202605151012"))
            payload = json.loads(
                server.sample_search(specs=["ILG_6h_1"], directory=str(directory)))
        self.assertIn("202605151012", payload["samples"][0]["files"][".pai2"])
        self.assertNotIn("202601010000", payload["samples"][0]["files"][".pai2"])

    def test_search_returns_all_roles_with_role_field(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(
                server.sample_search(specs=["raw"], directory=str(make_dir(tmp))))
        roles = {s["name"]: s["role"] for s in payload["samples"]}
        self.assertEqual(roles["20220901_QC_RAW_NEG_1"], "qc")
        self.assertEqual(payload["matched"], 4)

    def test_include_roles_narrows_the_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(server.sample_search(
                specs=["raw"], directory=str(make_dir(tmp)), include_roles=["sample"]))
        self.assertEqual(payload["matched"], 3)
        self.assertEqual(payload["excluded_by_role"], ["20220901_QC_RAW_NEG_1"])

    def test_vocabulary_mode_with_include_roles_drops_qc(self):
        # specs 省略（語彙モード）でも include_roles を適用し、QC/blank を
        # 結果とトークン語彙の両方から落とすことを検証する。
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(server.sample_search(
                directory=str(make_dir(tmp)), include_roles=["sample"]))
        self.assertEqual(payload["matched"], 3)
        self.assertEqual(payload["excluded_by_role"], ["20220901_QC_RAW_NEG_1"])
        self.assertNotIn(
            "20220901_QC_RAW_NEG_1",
            [s["name"] for s in payload["samples"]],
        )
        self.assertNotIn("qc", payload["token_vocabulary"]["tokens"].get("raw", {}).get("roles", {}))

    def test_vocabulary_mode_without_include_roles_keeps_qc(self):
        # include_roles 省略時は従来どおり全 role を含む（回帰防止）。
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(server.sample_search(directory=str(make_dir(tmp))))
        self.assertEqual(payload["matched"], 4)
        self.assertEqual(payload["excluded_by_role"], [])
        self.assertIn(
            "20220901_QC_RAW_NEG_1",
            [s["name"] for s in payload["samples"]],
        )

    def test_custom_extensions(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(server.sample_search(
                specs=["ILG_6h_1"], directory=str(make_dir(tmp)), extensions=[".dcl"]))
        self.assertEqual(list(payload["samples"][0]["files"]), [".dcl"])

    def test_zero_match_returns_error_with_vocabulary(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = json.loads(
                server.sample_search(specs=["24h"], directory=str(make_dir(tmp))))
        self.assertEqual(payload["status"], "error")
        self.assertIn("24h", payload["message"])
        self.assertIn("token_vocabulary", payload)

    def test_missing_directory_is_an_error(self):
        payload = json.loads(server.sample_search(directory="C:/no/such/dir"))
        self.assertEqual(payload["status"], "error")


class SampleSearchFromLoadedArfTests(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        session_state.session.arf.current_file_path = "loaded.arf"
        session_state.session.arf.features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [[i, n, 1.0] for i, n in enumerate(SAMPLES)],
        }]

    def test_uses_loaded_arf_sample_names(self):
        payload = json.loads(server.sample_search(specs=["ILG_6h"]))
        self.assertEqual(payload["source"], "loaded_arf")
        self.assertEqual(payload["matched"], 2)

    def test_file_id_comes_from_class_index_when_absent(self):
        payload = json.loads(server.sample_search(specs=["ILG_6h"]))
        self.assertIsNone(payload["samples"][0]["file_id"])

    def test_search_space_is_not_narrowed_by_a_prior_filter(self):
        """arf_parser の絞り込み（filtered_features）で検索対象が痩せないこと。

        sample_search は「どのサンプルが存在するか」の発見入口なので、直前の
        フィルタに関わらずデータセット全体を見る。
        """
        session_state.session.arf.filtered_features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [[1, "20220902_RAW_ILG_6h_1_NEG", 1.0]],
        }]
        payload = json.loads(server.sample_search())
        self.assertEqual(payload["total_samples"], len(SAMPLES))
        self.assertIn("control", payload["token_vocabulary"]["tokens"])

    def test_spec_outside_the_prior_filter_still_resolves(self):
        session_state.session.arf.filtered_features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [[1, "20220902_RAW_ILG_6h_1_NEG", 1.0]],
        }]
        payload = json.loads(server.sample_search(specs=["control"]))
        self.assertEqual(payload["status"], "success")
        self.assertEqual(
            [s["name"] for s in payload["samples"]], ["20220901_RAW_control_6h_1_NEG"])


if __name__ == "__main__":
    unittest.main()
