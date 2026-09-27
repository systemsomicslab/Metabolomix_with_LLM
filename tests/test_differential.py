import math
import unittest

import numpy as np

from lipidmix.analysis import differential as diff


class TestBhFdr(unittest.TestCase):
    def test_monotone_and_bounded(self):
        q = diff.bh_fdr([0.001, 0.01, 0.02, 0.5])
        self.assertTrue(all(0.0 <= v <= 1.0 for v in q))
        # smallest p gets smallest q
        self.assertLessEqual(q[0], q[3])

    def test_nan_preserved(self):
        q = diff.bh_fdr([0.01, float("nan"), 0.02])
        self.assertTrue(math.isnan(q[1]))
        self.assertFalse(math.isnan(q[0]))


class TestTwoGroup(unittest.TestCase):
    def test_clear_difference_significant(self):
        # feature 0 differs strongly between A and B; feature 1 does not.
        matrix = np.array([
            [10.0, 5.0],  # A
            [11.0, 5.1],  # A
            [ 9.5, 4.9],  # A
            [50.0, 5.0],  # B
            [52.0, 5.2],  # B
            [48.0, 4.8],  # B
        ])
        labels = ["A", "A", "A", "B", "B", "B"]
        res = diff.two_group_test(matrix, ["f0", "f1"], labels, "A", "B")
        by_feat = {r["feature"]: r for r in res}
        self.assertLess(by_feat["f0"]["p"], 0.05)
        self.assertGreater(by_feat["f1"]["p"], 0.05)
        # log2fc = log2(mean_b / mean_a)。f0 は A(~10) より B(50) が高いので正。
        self.assertGreater(by_feat["f0"]["log2fc"], 1.0)

    def test_missing_group_member_gives_nan(self):
        matrix = np.array([[10.0], [np.nan]])
        res = diff.two_group_test(matrix, ["f0"], ["A", "B"], "A", "B")
        self.assertTrue(math.isnan(res[0]["p"]))

    def test_log2fc_positive_when_group_b_higher(self):
        # 慣習: log2FC = log2(比較対象 / 基準)。group_a が基準、group_b が比較対象。
        # B(40) が A(10) より高いので正になる。log2((40+1)/(10+1)) ≒ 1.898
        matrix = np.array([[10.0], [10.0], [10.0], [40.0], [40.0], [40.0]])
        labels = ["A", "A", "A", "B", "B", "B"]
        res = diff.two_group_test(matrix, ["f0"], labels, "A", "B",
                                  log_transform=False)
        self.assertGreater(res[0]["log2fc"], 1.85)
        self.assertLess(res[0]["log2fc"], 1.95)

    def test_log2fc_sign_is_same_in_log_space(self):
        # log_transform=True の分岐でも向きが一致すること（2 分岐あるので両方を固定する）
        matrix = np.array([[10.0], [10.0], [10.0], [40.0], [40.0], [40.0]])
        labels = ["A", "A", "A", "B", "B", "B"]
        res = diff.two_group_test(matrix, ["f0"], labels, "A", "B",
                                  log_transform=True)
        self.assertGreater(res[0]["log2fc"], 1.85)
        self.assertLess(res[0]["log2fc"], 1.95)


    def test_t_has_the_same_sign_as_log2fc(self):
        # t は「group_b − group_a」。arf_differential の summary.top に log2fc と
        # 並んで出るので、向きが逆だと同じ行で符号が食い違う。両分岐を固定する。
        matrix = np.array([[10.0], [11.0], [9.0], [40.0], [42.0], [38.0]])
        labels = ["A", "A", "A", "B", "B", "B"]
        for log_transform in (False, True):
            res = diff.two_group_test(matrix, ["f0"], labels, "A", "B",
                                      log_transform=log_transform)
            self.assertGreater(res[0]["log2fc"], 0)
            self.assertGreater(res[0]["t"], 0, f"log_transform={log_transform}")


class TestAnova(unittest.TestCase):
    def test_three_group_difference(self):
        matrix = np.array([
            [1.0], [1.1], [0.9],     # G1
            [5.0], [5.1], [4.9],     # G2
            [9.0], [9.2], [8.8],     # G3
        ])
        labels = ["G1"] * 3 + ["G2"] * 3 + ["G3"] * 3
        res = diff.one_way_anova(matrix, ["f0"], labels)
        self.assertLess(res[0]["p"], 0.01)
        self.assertEqual(res[0]["n_groups"], 3)


class TestVolcanoAndConfound(unittest.TestCase):
    def test_add_fdr_and_volcano_flags(self):
        results = [
            {"feature": "f0", "log2fc": 2.0, "p": 0.001, "mean_a": 4, "mean_b": 1},
            {"feature": "f1", "log2fc": -3.0, "p": 0.002, "mean_a": 1, "mean_b": 8},
            {"feature": "f2", "log2fc": 0.1, "p": 0.9, "mean_a": 1, "mean_b": 1},
        ]
        withq = diff.add_fdr(results)
        self.assertIn("q", withq[0])
        pts = diff.volcano_data(withq, q_thr=0.05, log2fc_thr=1.0)
        flags = {p["feature"]: p["sig"] for p in pts}
        self.assertEqual(flags["f0"], "up")
        self.assertEqual(flags["f1"], "down")
        self.assertEqual(flags["f2"], "ns")

    def test_confounding_detected(self):
        groups = ["ctrl", "ctrl", "trt", "trt"]
        batches = ["d1", "d1", "d2", "d2"]  # group perfectly aligns with batch
        out = diff.check_confounding(groups, batches)
        self.assertTrue(out["confounded"])
        self.assertTrue(out["assessable"])

    def test_single_batch_is_not_assessable(self):
        # only one batch present -> cannot assess confounding, must NOT read as "clean".
        out = diff.check_confounding(["ctrl", "ctrl", "trt", "trt"],
                                     ["d1", "d1", "d1", "d1"])
        self.assertFalse(out["confounded"])
        self.assertFalse(out["assessable"])

    def test_no_batch_is_not_assessable(self):
        out = diff.check_confounding(["ctrl", "trt"], [None, None])
        self.assertFalse(out["assessable"])


class TestLogTransform(unittest.TestCase):
    def test_log_transform_uses_log_space_fold_change(self):
        # B is ~5x A on the linear scale for f0.
        matrix = np.array([
            [10.0], [11.0], [9.5],
            [50.0], [52.0], [48.0],
        ])
        labels = ["A", "A", "A", "B", "B", "B"]
        res = diff.two_group_test(matrix, ["f0"], labels, "A", "B", log_transform=True)
        # log2 空間の群平均差 ~ log2(50/10) ~ +2.3。生の平均は従来どおり返る。
        self.assertGreater(res[0]["log2fc"], 2.0)
        self.assertGreater(res[0]["mean_a"], 9.0)  # raw mean, not log mean
        self.assertLess(res[0]["p"], 0.05)

    def test_anova_log_transform_runs(self):
        matrix = np.array([[1.0], [1.1], [0.9], [50.0], [55.0], [45.0]])
        labels = ["G1", "G1", "G1", "G2", "G2", "G2"]
        res = diff.one_way_anova(matrix, ["f0"], labels, log_transform=True)
        self.assertLess(res[0]["p"], 0.05)


class TestFallbackDistributions(unittest.TestCase):
    """scipy 不在フォールバックが scipy と一致することを担保する（近似ではない）。"""

    def test_t_sf_matches_scipy(self):
        from scipy import stats
        for t, df in [(2.5, 4.0), (1.1, 12.0), (3.3, 2.0), (0.5, 30.0)]:
            expected = 2.0 * stats.t.sf(abs(t), df)
            self.assertAlmostEqual(diff._t_sf_two_sided(t, df), expected, places=8)

    def test_f_sf_matches_scipy_including_odd_df(self):
        from scipy import stats
        # odd df1 was the case the old chi2 series got wrong.
        for f, d1, d2 in [(4.0, 3, 8), (2.0, 2, 10), (7.5, 5, 6), (1.5, 1, 20)]:
            self.assertAlmostEqual(diff._f_sf(f, d1, d2), stats.f.sf(f, d1, d2), places=8)


if __name__ == "__main__":
    unittest.main()
