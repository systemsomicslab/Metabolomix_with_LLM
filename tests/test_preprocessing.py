import unittest
import numpy as np
from metabolomix.analysis import preprocessing as pp


class TestDetectQcStrata(unittest.TestCase):
    def test_detects_region_stratified_qc(self):
        # real brain naming: QC is stratified by tissue region
        names = ["20240311_QC_Cerebellum_ICR_NEG_1",
                 "20240311_QC_Hippocampus_ICR_NEG_2",
                 "20240311_Cerebellum_gf_NC_1_NEG"]
        roles = pp.detect_sample_roles(names)
        strata = pp.detect_qc_strata(names, roles)
        self.assertGreater(len(strata), 1)

    def test_uniform_qc_is_single_stratum(self):
        names = ["20220901_QC_1_NEG", "20220901_QC_2_NEG", "20220901_ctrl_1_NEG"]
        roles = pp.detect_sample_roles(names)
        strata = pp.detect_qc_strata(names, roles)
        self.assertLessEqual(len(strata), 1)


class TestPreprocessReportClarity(unittest.TestCase):
    def test_reports_total_features_removed(self):
        m = np.array([[1.0, 2.0, 3.0, 4.0], [1.5, 2.5, 3.5, 4.5]])
        names = ["s1", "s2"]
        roles = {"s1": "sample", "s2": "sample"}
        order = {"s1": 1, "s2": 2}
        _, _, rep = pp.preprocess(m, names, roles, order, {"impute": "none"})
        self.assertEqual(rep["features_removed_total"],
                         rep["features_before"] - rep["features_after"])

    def test_warns_when_all_features_removed(self):
        # QC has variance in every feature -> max_qc_rsd=0.0 removes them all
        m = np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [5.0, 6.0, 7.0], [6.0, 7.0, 8.0]])
        names = ["a_qc", "b_qc", "c", "d"]
        roles = pp.detect_sample_roles(names)
        order = {n: i for i, n in enumerate(names)}
        _, _, rep = pp.preprocess(m, names, roles, order,
                                  {"max_qc_rsd": 0.0, "impute": "none"})
        self.assertEqual(rep["features_after"], 0)
        self.assertTrue(any("残存" in c or "除去" in c for c in rep["caveats"]),
                        f"expected extreme-removal caveat, got {rep['caveats']}")


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


class TestNormalize(unittest.TestCase):
    def test_none_is_identity(self):
        m = np.array([[1.0, 2.0], [3.0, 4.0]])
        out, factors, report = pp.normalize(m, "none")
        np.testing.assert_array_equal(out, m)
        self.assertEqual(report["method"], "none")

    def test_tic_divides_by_row_sum(self):
        m = np.array([[1.0, 1.0], [2.0, 2.0]])
        out, factors, report = pp.normalize(m, "tic")
        # each row sums to 2 and 4 -> scaled so row sums equal the mean row sum (3)
        self.assertAlmostEqual(out[0, 0] / out[1, 0], (1 / 2) / (2 / 4))

    def test_median_uses_row_median(self):
        m = np.array([[2.0, 4.0], [10.0, 20.0]])
        out, factors, report = pp.normalize(m, "median")
        # raw row medians are 3.0 and 15.0; factors is normalized by their
        # median, but the ratio between rows must be preserved.
        self.assertAlmostEqual(factors[0] / factors[1], 3.0 / 15.0)

    def test_unknown_method_raises(self):
        with self.assertRaises(ValueError):
            pp.normalize(np.zeros((2, 2)), "bogus")

    def test_tic_factors_match_applied_divisor(self):
        m = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 9.0]])
        out, factors, report = pp.normalize(m, "tic")
        np.testing.assert_allclose(out, m / factors[:, None])

    def test_median_factors_match_applied_divisor(self):
        m = np.array([[2.0, 4.0], [10.0, 20.0], [1.0, 7.0]])
        out, factors, report = pp.normalize(m, "median")
        np.testing.assert_allclose(out, m / factors[:, None])

    def test_pqn_uses_qc_reference_when_roles_provided(self):
        m = np.array([
            [2.0, 4.0],
            [4.0, 8.0],
            [1.0, 2.0],
        ])
        sample_names = ["s1_qc", "s2_qc", "s3_sample"]
        roles = pp.detect_sample_roles(sample_names)
        out, factors, report = pp.normalize(
            m, "pqn", roles=roles, sample_names=sample_names
        )
        self.assertEqual(report["pqn_reference"], "qc_median")

    def test_pqn_uses_all_sample_reference_without_roles(self):
        m = np.array([[2.0, 4.0], [4.0, 8.0], [1.0, 2.0]])
        out, factors, report = pp.normalize(m, "pqn")
        self.assertEqual(report["pqn_reference"], "all_sample_median")

    def test_pqn_handles_zero_in_reference_row_without_inf_or_crash(self):
        # column 1's reference (median) is 0.0, which would otherwise cause
        # a division-by-zero when computing quotients.
        m = np.array([[2.0, 0.0], [4.0, 0.0], [1.0, 3.0]])
        out, factors, report = pp.normalize(m, "pqn")
        # the zero-reference column becomes NaN (not inf), other columns stay finite.
        self.assertTrue(np.all(np.isfinite(out[:, 0])))
        self.assertFalse(np.any(np.isinf(out)))

    def test_median_does_not_destroy_zero_median_samples(self):
        # A sparse sample whose row median is 0 (>50% features undetected=0) must
        # NOT be nuked to all-NaN; leave it unscaled and flag it in the report.
        m = np.array([
            [10.0, 12.0, 11.0, 9.0],  # dense: median 10.5
            [0.0, 0.0, 0.0, 5.0],     # sparse: median 0 -> would divide by zero
        ])
        out, factors, report = pp.normalize(m, "median")
        self.assertTrue(np.all(np.isfinite(out[1])),
                        "zero-median sample must stay finite, not become NaN")
        self.assertEqual(report["unscaled_samples"], 1)
        self.assertIn("caveat", report)

    def test_median_unscaled_sample_keeps_raw_values(self):
        m = np.array([[10.0, 12.0, 11.0, 9.0], [0.0, 0.0, 0.0, 5.0]])
        out, factors, report = pp.normalize(m, "median")
        # degenerate sample gets factor 1.0 -> values unchanged
        np.testing.assert_allclose(out[1], m[1])

    def test_preprocess_flags_normalization_dropped_samples(self):
        m = np.array([[10.0, 12.0, 11.0, 9.0], [0.0, 0.0, 0.0, 5.0]])
        names = ["s1", "s2"]
        roles = {"s1": "sample", "s2": "sample"}
        order = {"s1": 1, "s2": 2}
        _, _, report = pp.preprocess(
            m, names, roles, order, {"normalize": "median", "impute": "none"})
        self.assertTrue(any("正規化" in c for c in report["caveats"]),
                        f"expected normalization caveat, got {report['caveats']}")


class TestBlankFilter(unittest.TestCase):
    def test_background_feature_removed(self):
        # col0: sample >> blank (keep). col1: sample ~ blank (remove).
        m = np.array([
            [100.0, 10.0],  # sample
            [120.0, 11.0],  # sample
            [5.0, 9.0],     # blank
        ])
        roles = {"s1": "sample", "s2": "sample", "b1": "blank"}
        names = ["s1", "s2", "b1"]
        mask, report = pp.blank_filter(m, roles, names, min_fold=3.0)
        self.assertTrue(mask[0])
        self.assertFalse(mask[1])
        self.assertEqual(report["removed"], 1)

    def test_no_blank_keeps_all_with_caveat(self):
        m = np.array([[1.0, 2.0]])
        mask, report = pp.blank_filter(m, {"s1": "sample"}, ["s1"])
        self.assertTrue(mask.all())
        self.assertIn("caveat", report)
        # Task 11: caveatだけでなくstatusでも「未実施」を機械可読にする
        # （呼び出し側のpreprocess()がcaveatの文面をパースせずに済むように）。
        self.assertEqual(report["status"], "skipped")


class TestDriftCorrect(unittest.TestCase):
    def test_linear_drift_flattened(self):
        # feature drifts linearly with run order across QC; correction should
        # bring QC intensities close to a constant level.
        order = [1, 2, 3, 4, 5, 6]
        names = [f"n{i}" for i in order]
        roles = {n: ("qc" if i % 2 == 1 else "sample") for n, i in zip(names, order)}
        run_order = dict(zip(names, order))
        base = np.array([o * 10.0 for o in order])  # 10,20,30,40,50,60
        m = base.reshape(-1, 1)
        out, report = pp.qc_drift_correct(m, roles, names, run_order, min_qc=3)
        qc_vals = [out[i, 0] for i, n in enumerate(names) if roles[n] == "qc"]
        self.assertLess(np.std(qc_vals), np.std([10, 30, 50]))
        self.assertEqual(report["status"], "applied")

    def test_no_run_order_skips(self):
        m = np.array([[1.0], [2.0]])
        names = ["a", "b"]
        roles = {"a": "qc", "b": "sample"}
        out, report = pp.qc_drift_correct(m, roles, names, {"a": None, "b": None})
        np.testing.assert_array_equal(out, m)
        self.assertEqual(report["status"], "skipped")


class TestQcRsdFilter(unittest.TestCase):
    def test_high_rsd_feature_removed(self):
        # col0: stable in QC (keep). col1: highly variable in QC (remove).
        m = np.array([
            [100.0, 10.0],   # qc
            [101.0, 90.0],   # qc
            [ 99.0, 50.0],   # qc
            [ 50.0, 40.0],   # sample (ignored for RSD)
        ])
        names = ["q1", "q2", "q3", "s1"]
        roles = {"q1": "qc", "q2": "qc", "q3": "qc", "s1": "sample"}
        mask, report = pp.qc_rsd_filter(m, roles, names, max_rsd=0.30)
        self.assertTrue(mask[0])
        self.assertFalse(mask[1])

    def test_no_qc_keeps_all_with_caveat(self):
        m = np.array([[1.0, 2.0]])
        mask, report = pp.qc_rsd_filter(m, {"s1": "sample"}, ["s1"])
        self.assertTrue(mask.all())
        self.assertIn("caveat", report)
        # Task 11: blank_filterと対称にstatus="skipped"を持つ（同じ理由）。
        self.assertEqual(report["status"], "skipped")


class TestImpute(unittest.TestCase):
    def test_half_min_fills_per_feature(self):
        m = np.array([[10.0, np.nan], [20.0, 4.0], [30.0, 8.0]])
        out, report = pp.impute(m, "half_min")
        self.assertAlmostEqual(out[0, 1], 2.0)  # half of column-1 min (4) = 2
        self.assertFalse(np.isnan(out).any())

    def test_none_keeps_nan(self):
        m = np.array([[np.nan, 1.0]])
        out, report = pp.impute(m, "none")
        self.assertTrue(np.isnan(out[0, 0]))

    def test_all_nan_column_fills_zero(self):
        # column 1 is entirely NaN -> fallback fill value is 0.0
        m = np.array([[10.0, np.nan], [20.0, np.nan], [30.0, np.nan]])
        out_hm, _ = pp.impute(m, "half_min")
        self.assertTrue((out_hm[:, 1] == 0.0).all())
        self.assertFalse(np.isnan(out_hm).any())
        out_cm, _ = pp.impute(m, "column_mean")
        self.assertTrue((out_cm[:, 1] == 0.0).all())


class TestPreprocessOrchestration(unittest.TestCase):
    def test_recipe_applied_order_and_caveats(self):
        m = np.array([
            [100.0, 10.0, 5.0],  # sample
            [120.0, 90.0, 5.0],  # sample
            [110.0, 50.0, 5.0],  # qc
        ])
        names = ["s1", "s2", "q1"]
        roles = {"s1": "sample", "s2": "sample", "q1": "qc"}
        run_order = {"s1": 1, "s2": 2, "q1": 3}
        recipe = {"normalize": "median", "impute": "half_min",
                  "drift_correct": True, "max_qc_rsd": 0.30}
        out, kept_idx, report = pp.preprocess(m, names, roles, run_order, recipe)
        self.assertIn("normalize", report["recipe_applied"])
        # drift correction has only 1 QC -> skipped caveat present
        self.assertTrue(any("ドリフト" in c for c in report["caveats"]))
        self.assertEqual(out.shape[0], 3)


if __name__ == "__main__":
    unittest.main()


class TestQcInterspersion(unittest.TestCase):
    """QC-RLSC は QC が試料列に散在することを前提にする（kidney aging 実データで露見）。"""

    def _fixture(self, qc_orders, sample_orders):
        names = [f"QC_{i}" for i in range(len(qc_orders))] + \
                [f"S_{i}" for i in range(len(sample_orders))]
        roles = {n: ("qc" if n.startswith("QC") else "sample") for n in names}
        run_order = dict(zip(names, list(qc_orders) + list(sample_orders)))
        matrix = np.arange(1, len(names) * 3 + 1, dtype=float).reshape(len(names), 3)
        return matrix, names, roles, run_order

    def test_qc_appended_after_samples_is_skipped(self):
        # 試料 1-48 → Blank → QC 50-56。内挿が全て外挿になるので補正してはいけない。
        matrix, names, roles, run_order = self._fixture(range(50, 57), range(1, 49))
        out, rep = pp.qc_drift_correct(matrix, roles, names, run_order)
        self.assertEqual(rep["status"], "skipped")
        self.assertIn("挿入されていません", rep["caveat"])
        self.assertEqual(rep["qc_interspersion"]["covered"], 0)
        np.testing.assert_array_equal(out, matrix)

    def test_interspersed_qc_is_applied(self):
        matrix, names, roles, run_order = self._fixture([1, 5, 9, 13], [2, 3, 6, 7, 10, 11])
        _, rep = pp.qc_drift_correct(matrix, roles, names, run_order)
        self.assertEqual(rep["status"], "applied")
        self.assertEqual(rep["qc_interspersion"]["covered"], 6)
        self.assertNotIn("caveat", rep)

    def test_partial_coverage_is_caveated(self):
        matrix, names, roles, run_order = self._fixture([1, 2, 3, 4], [2, 3, 20, 21, 22])
        _, rep = pp.qc_drift_correct(matrix, roles, names, run_order)
        self.assertEqual(rep["status"], "applied")
        self.assertIn("外挿", rep["caveat"])


class TestDetectFailedQc(unittest.TestCase):
    """失敗 QC 注入を残すと qc_rsd_filter がほぼ全特徴を落とす（1345 -> 51 を実測）。"""

    def _fixture(self, qc_rows):
        names = [f"QC_{i}" for i in range(len(qc_rows))] + ["S_0", "S_1"]
        roles = {n: ("qc" if n.startswith("QC") else "sample") for n in names}
        matrix = np.array(list(qc_rows) + [[100.0, 100.0], [100.0, 100.0]])
        return matrix, names, roles

    def test_flags_low_intensity_qc(self):
        matrix, names, roles = self._fixture([
            [100.0, 100.0], [102.0, 98.0], [99.0, 101.0], [1.0, 0.5],
        ])
        rep = pp.detect_failed_qc(matrix, roles, names)
        self.assertEqual(rep["failed"], ["QC_3"])
        self.assertIn("arf_exclude", rep["caveat"])

    def test_healthy_qc_is_silent(self):
        matrix, names, roles = self._fixture([
            [100.0, 100.0], [102.0, 98.0], [99.0, 101.0], [98.0, 103.0],
        ])
        rep = pp.detect_failed_qc(matrix, roles, names)
        self.assertEqual(rep["failed"], [])
        self.assertNotIn("caveat", rep)

    def test_too_few_qc_is_not_judged(self):
        matrix, names, roles = self._fixture([[100.0, 100.0], [1.0, 1.0]])
        rep = pp.detect_failed_qc(matrix, roles, names)
        self.assertEqual(rep["failed"], [])


class TestDropSamplesByRole(unittest.TestCase):
    """ブランクは背景除去に使ったあと解析行列から外す（PCA/差次を支配させない）。

    QC は PCA に残す（クラスタの締まり具合の確認は PCA の主目的の一つ）ため、
    既定の除外対象は blank のみ。
    """

    def _matrix(self):
        return np.array([
            [10.0, 1.0],
            [20.0, 2.0],
            [15.0, 1.5],
            [0.1, 0.01],
        ])

    NAMES = ["s1", "s2", "qc1", "blank1"]
    ROLES = {"s1": "sample", "s2": "sample", "qc1": "qc", "blank1": "blank"}

    def test_drops_blank_rows_only_by_default(self):
        matrix, names, dropped = pp.drop_samples_by_role(
            self._matrix(), self.NAMES, self.ROLES,
        )
        self.assertEqual(names, ["s1", "s2", "qc1"])
        self.assertEqual(dropped, {"blank": ["blank1"]})
        self.assertEqual(matrix.shape, (3, 2))
        np.testing.assert_allclose(matrix[2], [15.0, 1.5])  # qc1 が残る

    def test_can_drop_multiple_roles(self):
        matrix, names, dropped = pp.drop_samples_by_role(
            self._matrix(), self.NAMES, self.ROLES, drop_roles=("blank", "qc"),
        )
        self.assertEqual(names, ["s1", "s2"])
        self.assertEqual(dropped, {"blank": ["blank1"], "qc": ["qc1"]})

    def test_no_matching_role_is_a_noop(self):
        names_only_samples = {"s1": "sample", "s2": "sample"}
        matrix, names, dropped = pp.drop_samples_by_role(
            np.array([[1.0], [2.0]]), ["s1", "s2"], names_only_samples,
        )
        self.assertEqual(names, ["s1", "s2"])
        self.assertEqual(dropped, {})
        self.assertEqual(matrix.shape, (2, 1))

    def test_dropping_every_row_yields_empty_matrix_not_error(self):
        matrix, names, dropped = pp.drop_samples_by_role(
            np.array([[1.0], [2.0]]), ["b1", "b2"], {"b1": "blank", "b2": "blank"},
        )
        self.assertEqual(names, [])
        self.assertEqual(dropped, {"blank": ["b1", "b2"]})
        self.assertEqual(matrix.shape[0], 0)


class TestRecipeHonesty(unittest.TestCase):
    """要求したステップが実際に適用されたかを、報告で区別する。

    `recipe_applied` に "drift_correct" が載るのに中身は status=skipped、という
    状態だと、「やった」と「できなかった」を呼び出し側が区別できない。
    """

    def _matrix(self):
        return np.array([[10.0, 20.0], [12.0, 21.0], [11.0, 19.0], [13.0, 22.0]])

    def test_requested_drift_correct_without_qc_is_reported_as_skipped(self):
        names = ["s1", "s2", "s3", "s4"]
        roles = {n: "sample" for n in names}
        run_order = {n: i + 1 for i, n in enumerate(names)}
        _, _, report = pp.preprocess(self._matrix(), names, roles, run_order,
                                  {"drift_correct": True})
        self.assertNotIn("drift_correct", report["recipe_applied"])
        self.assertIn("drift_correct", report["recipe_skipped"])

    def test_applied_steps_stay_in_recipe_applied(self):
        names = ["s1", "s2", "s3", "s4"]
        roles = {n: "sample" for n in names}
        run_order = {n: None for n in names}
        _, _, report = pp.preprocess(self._matrix(), names, roles, run_order,
                                  {"normalize": "tic"})
        self.assertIn("normalize", report["recipe_applied"])
        self.assertIn("impute", report["recipe_applied"])
        self.assertEqual(report["recipe_skipped"], [])

    def test_recipe_skipped_is_always_present(self):
        """キーの有無で分岐させると、呼び出し側が毎回 get() を書くことになる。"""
        names = ["s1", "s2", "s3", "s4"]
        roles = {n: "sample" for n in names}
        _, _, report = pp.preprocess(self._matrix(), names, roles,
                                  {n: None for n in names}, {})
        self.assertEqual(report["recipe_skipped"], [])

    def test_requested_blank_filter_without_blank_is_reported_as_skipped(self):
        """Task 11: blank_filterがcaveatだけ返す（statusが無い）と、ここで
        誤って recipe_applied に数えられていた（QCが無いバッチでblank_min_foldを
        要求したのに「適用した」ことになる）。"""
        names = ["s1", "s2", "s3", "s4"]
        roles = {n: "sample" for n in names}
        run_order = {n: i + 1 for i, n in enumerate(names)}
        _, _, report = pp.preprocess(self._matrix(), names, roles, run_order,
                                  {"blank_min_fold": 3.0})
        self.assertNotIn("blank_filter", report["recipe_applied"])
        self.assertIn("blank_filter", report["recipe_skipped"])

    def test_requested_qc_rsd_filter_without_qc_is_reported_as_skipped(self):
        names = ["s1", "s2", "s3", "s4"]
        roles = {n: "sample" for n in names}
        run_order = {n: i + 1 for i, n in enumerate(names)}
        _, _, report = pp.preprocess(self._matrix(), names, roles, run_order,
                                  {"max_qc_rsd": 0.30})
        self.assertNotIn("qc_rsd_filter", report["recipe_applied"])
        self.assertIn("qc_rsd_filter", report["recipe_skipped"])


class TestNoQcNoBlankBatch(unittest.TestCase):
    """QC もブランクも無いバッチでは、前処理の大半が原理的に効かない。

    それでも preprocess は success を返すので、`recipe_applied` を読まない限り
    「normalize と impute しかしていない」ことが分からない。実データの POS 60
    サンプルがこれで、60/60 が role=sample だった。
    """

    def _run(self, roles):
        names = list(roles)
        m = np.arange(len(names) * 3, dtype=float).reshape(len(names), 3) + 1.0
        return pp.preprocess(m, names, roles, {n: None for n in names}, {})[2]

    def test_caveat_when_neither_qc_nor_blank_present(self):
        report = self._run({f"s{i}": "sample" for i in range(4)})
        self.assertTrue(any("QC" in c and "ブランク" in c for c in report["caveats"]),
                        report["caveats"])

    def test_no_caveat_when_qc_present(self):
        roles = {f"s{i}": "sample" for i in range(3)}
        roles["qc_1"] = "qc"
        report = self._run(roles)
        self.assertFalse(any("いずれも適用できません" in c for c in report["caveats"]),
                         report["caveats"])

    def test_no_caveat_when_blank_present(self):
        roles = {f"s{i}": "sample" for i in range(3)}
        roles["blank_1"] = "blank"
        report = self._run(roles)
        self.assertFalse(any("いずれも適用できません" in c for c in report["caveats"]),
                         report["caveats"])

    def test_caveat_ignores_role_entry_not_in_sample_names(self):
        """rolesに残る評価対象外のエントリで、注意書きの判定を誤らせない。

        呼び出し側（dataset_analysis.build_dataset_pp_inputs）はinclude=falseの
        試料をmatrix/sample_namesから落とす一方、roles辞書は監査・表示用に
        全件分（除外分込み）のまま渡す設計にしている。present_rolesが
        `set(roles.values())`のようにroles辞書の全件を見ると、実際には
        渡されていない（=sample_namesに無い）"qc_stale"のroleを拾って
        しまい、実評価対象にQCが無いのに「QCがある」と誤認する。
        """
        names = [f"s{i}" for i in range(4)]
        roles = {n: "sample" for n in names}
        roles["qc_stale"] = "qc"  # sample_names には含まれない、除外済みQCの残骸
        m = np.arange(len(names) * 3, dtype=float).reshape(len(names), 3) + 1.0
        _, _, report = pp.preprocess(m, names, roles, {n: None for n in names}, {})
        self.assertTrue(any("QC" in c and "ブランク" in c for c in report["caveats"]),
                        report["caveats"])
