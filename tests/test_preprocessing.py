"""試料名・Class ID のトークンから役割（qc / blank / sample）を推定する（`preprocessing.detect_sample_roles`）。

v1 の前処理（正規化・足切り・ドリフト補正・補完）は D 段で消えた。前処理は統計エンジン
（`analysis/stats/preprocess.py`。`tests/test_stats_preprocess.py`）が持つ。
"""
import unittest

from metabolomix.analysis import preprocessing as pp


class TestDetectSampleRoles(unittest.TestCase):
    def test_qc_and_blank_detected_from_filename(self):
        names = [
            "20240311_QC_Cerebellum_ICR_NEG_1",
            "20240311_blank_Cerebellum_NEG_1",
            "20240311_Cerebellum_gf_AIN_1_NEG",
        ]
        roles = pp.detect_sample_roles(names)
        self.assertEqual(roles[names[0]], "qc")
        self.assertEqual(roles[names[1]], "blank")
        self.assertEqual(roles[names[2]], "sample")

    def test_role_from_class_id_token(self):
        names = ["s1"]
        roles = pp.detect_sample_roles(names, class_ids={"s1": "QC_pool"})
        self.assertEqual(roles["s1"], "qc")

    def test_config_overrides_tokens(self):
        names = ["ctrl_pooledqc_1"]
        roles = pp.detect_sample_roles(names, config={"qc_tokens": ["pooledqc"]})
        self.assertEqual(roles["ctrl_pooledqc_1"], "qc")
