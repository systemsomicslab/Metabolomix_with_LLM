import math
import unittest

from metabolomix.analysis import differential as diff


class TestVolcano(unittest.TestCase):
    def test_volcano_flags(self):
        results = [
            {"feature": "f0", "log2fc": 2.0, "p": 0.001, "q": 0.003, "mean_a": 4, "mean_b": 1},
            {"feature": "f1", "log2fc": -3.0, "p": 0.002, "q": 0.003, "mean_a": 1, "mean_b": 8},
            {"feature": "f2", "log2fc": 0.1, "p": 0.9, "q": 0.9, "mean_a": 1, "mean_b": 1},
        ]
        pts = diff.volcano_data(results, q_thr=0.05, log2fc_thr=1.0)
        flags = {p["feature"]: p["sig"] for p in pts}
        self.assertEqual(flags["f0"], "up")
        self.assertEqual(flags["f1"], "down")
        self.assertEqual(flags["f2"], "ns")

    def test_volcano_p_underflow_to_zero_is_drawn_at_the_float_floor(self):
        """p が 0.0 に潰れた（アンダーフロー）検定済みの行は NaN（検定不能）にせず、最小の正の float で
        -log10 を取る（R17）。保存した p は変えない。印の無い行と統計エンジンの行で同じ。"""
        import sys
        rows = [{"feature": "flagless", "log2fc": 3.0, "p": 0.0, "q": 0.0},
                {"feature": "engine", "log2fc": -3.0, "p": 0.0, "q": 0.0, "significant": True},
                {"feature": "neg", "log2fc": 3.0, "p": -1.0, "q": 0.0}]
        pts = {p["feature"]: p for p in diff.volcano_data(rows, q_thr=0.05, log2fc_thr=1.0)}
        floor = -math.log10(sys.float_info.min)
        self.assertAlmostEqual(floor, 307.6526555685888, places=9)
        self.assertEqual(pts["flagless"]["neg_log10_p"], floor)
        self.assertEqual(pts["flagless"]["sig"], "up")
        self.assertEqual(pts["engine"]["neg_log10_p"], floor)
        self.assertEqual(pts["engine"]["sig"], "down")
        self.assertTrue(math.isnan(pts["neg"]["neg_log10_p"]))          # 負の p は値として扱わない
        self.assertEqual(rows[0]["p"], 0.0)                              # 保存した p はそのまま

    def test_volcano_row_flag_wins_over_q_and_log2fc(self):
        """行が `significant` を持てば（統計エンジンの行。R16）、色はその印だけで決まる。"""
        rows = [{"feature": "judged_out", "log2fc": 4.0, "p": 1e-20, "q": 1e-18, "significant": False},
                {"feature": "kept", "log2fc": -4.0, "p": 1e-20, "q": 1e-18, "significant": True}]
        flags = {p["feature"]: p["sig"] for p in diff.volcano_data(rows, q_thr=0.05, log2fc_thr=1.0)}
        self.assertEqual(flags, {"judged_out": "ns", "kept": "down"})
