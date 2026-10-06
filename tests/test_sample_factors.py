import unittest

from metabolomix.msdial.classes import assign_sample_groups, filter_arf_by_class_ids
from metabolomix.msdial.tags import normalize_sample_name
from metabolomix.msdial.sample_factors import (
    SampleFacet,
    arf_sample_names,
    assign_factor_groups,
    build_sample_facets,
    expand_sample_specs,
    sample_tokens,
    split_tokens,
    token_vocabulary,
)


def name_class_index(name_to_class: dict[str, str]) -> dict:
    """discover_arf_class_index の戻りのうち、本モジュールが使う部分だけを作る。"""
    records = []
    by_file_name = {}
    for file_id, (name, class_id) in enumerate(name_to_class.items()):
        record = {"file_id": file_id, "file_name": name, "class_id": class_id}
        records.append(record)
        by_file_name[normalize_sample_name(name, strip_processing_timestamp=False)] = record
    return {"records": records, "by_file_id": {}, "by_file_name": by_file_name}


class SplitTokensTests(unittest.TestCase):
    def test_splits_on_underscore_and_casefolds(self):
        self.assertEqual(split_tokens("Cerebellum_gf_AIN"), frozenset({"cerebellum", "gf", "ain"}))

    def test_none_and_empty_yield_empty_set(self):
        self.assertEqual(split_tokens(None), frozenset())
        self.assertEqual(split_tokens("__"), frozenset())


class SampleTokensTests(unittest.TestCase):
    def test_strips_trailing_processing_timestamp(self):
        # 末尾12桁は MS-DIAL の処理タイムスタンプ（再処理ごとに変わる）で実験因子ではない。
        self.assertEqual(
            sample_tokens("20220902_RAW_ILG_6h_2_NEG_202605151012"),
            frozenset({"20220902", "raw", "ilg", "6h", "2", "neg"}),
        )

    def test_strips_known_measurement_suffix(self):
        self.assertEqual(sample_tokens("sample_A.wiff"), frozenset({"sample", "a"}))

    def test_merges_class_id_tokens(self):
        self.assertEqual(
            sample_tokens("s1", class_id="Cerebellum_gf_AIN"),
            frozenset({"s1", "cerebellum", "gf", "ain"}),
        )

    def test_merges_extra_label_tokens(self):
        self.assertEqual(sample_tokens("s1", extra="24M_GF"), frozenset({"s1", "24m", "gf"}))


class BuildSampleFacetsTests(unittest.TestCase):
    def test_multi_token_value_splits_into_two_tokens(self):
        # G_uralensis は2トークンに割れる。だから位置インデックスでの因子指定は使わない。
        name = "20220902_RAW_G_uralensis_6h_2_NEG"
        facets = build_sample_facets([name], name_class_index({name: "G"}))
        self.assertEqual(
            facets[name].tokens,
            frozenset({"20220902", "raw", "g", "uralensis", "6h", "2", "neg"}),
        )
        self.assertEqual(facets[name].class_id, "G")
        self.assertEqual(facets[name].file_id, 0)
        self.assertEqual(facets[name].role, "sample")

    def test_works_without_class_index(self):
        # .mddata が無いフォルダでもサンプル名トークンだけで成立する。
        facets = build_sample_facets(["20220901_RAW_LPS_6h_1_NEG"], None)
        facet = facets["20220901_RAW_LPS_6h_1_NEG"]
        self.assertIsNone(facet.class_id)
        self.assertIsNone(facet.file_id)
        self.assertIn("lps", facet.tokens)
        self.assertIn("6h", facet.tokens)

    def test_detects_qc_and_blank_roles(self):
        facets = build_sample_facets(["20240311_QC_Cerebellum_NEG_1", "Blank_01", "s1"], None)
        self.assertEqual(facets["20240311_QC_Cerebellum_NEG_1"].role, "qc")
        self.assertEqual(facets["Blank_01"].role, "blank")
        self.assertEqual(facets["s1"].role, "sample")

    def test_sample_meta_supplies_role_and_group_tokens(self):
        # arf_differential は class_index を持たず sample_meta["group"] だけを持つ経路。
        facets = build_sample_facets(
            ["s0", "qc1"],
            None,
            sample_meta={"s0": {"group": "24M_GF", "role": "sample"},
                         "qc1": {"group": "24M_GF", "role": "qc"}},
        )
        self.assertEqual(facets["s0"].tokens, frozenset({"s0", "24m", "gf"}))
        self.assertEqual(facets["qc1"].role, "qc")

    def test_preserves_input_order(self):
        facets = build_sample_facets(["b", "a", "c"], None)
        self.assertEqual(list(facets), ["b", "a", "c"])

    def test_is_frozen_dataclass(self):
        facet = build_sample_facets(["s1"], None)["s1"]
        self.assertIsInstance(facet, SampleFacet)
        with self.assertRaises(Exception):
            facet.name = "other"


class ArfSampleNamesTests(unittest.TestCase):
    def test_collects_names_in_first_appearance_order_without_duplicates(self):
        features = [
            {"AlignedPeakProperties": [[0, "sB", 1.0], [1, "sA", 2.0]]},
            {"AlignedPeakProperties": [[0, "sB", 3.0], [1, "sA", 4.0]]},
        ]
        self.assertEqual(arf_sample_names(features), ["sB", "sA"])

    def test_decodes_bytes_names_and_skips_malformed_rows(self):
        features = [{"AlignedPeakProperties": [[0, b"sA", 1.0], "not-a-row", [1]]}]
        self.assertEqual(arf_sample_names(features), ["sA"])

    def test_empty_input(self):
        self.assertEqual(arf_sample_names([]), [])
        self.assertEqual(arf_sample_names(None), [])


class ExpandSampleSpecsTests(unittest.TestCase):
    def setUp(self):
        self.facets = build_sample_facets([
            "20220901_RAW_control_6h_1_NEG",
            "20220901_RAW_LPS_6h_1_NEG",
            "20220902_RAW_ILG_0h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_RAW_ILG_6h_2_NEG",
        ], None)

    def test_single_token_matches_all_samples_with_token(self):
        matches, _ = expand_sample_specs(["ILG"], self.facets)
        self.assertEqual(matches["ILG"], [
            "20220902_RAW_ILG_0h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_RAW_ILG_6h_2_NEG",
        ])

    def test_multi_token_spec_is_and(self):
        matches, _ = expand_sample_specs(["ILG_6h"], self.facets)
        self.assertEqual(matches["ILG_6h"], [
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_RAW_ILG_6h_2_NEG",
        ])

    def test_token_absent_from_class_id_still_matches(self):
        # 6h は Class ID に無くサンプル名にしかない因子。これが本機能の眼目。
        matches, _ = expand_sample_specs(["6h"], self.facets)
        self.assertEqual(len(matches["6h"]), 4)

    def test_multiple_specs_are_independent(self):
        matches, _ = expand_sample_specs(["ILG_6h", "control_6h"], self.facets)
        self.assertEqual(len(matches["ILG_6h"]), 2)
        self.assertEqual(matches["control_6h"], ["20220901_RAW_control_6h_1_NEG"])

    def test_case_insensitive(self):
        matches, _ = expand_sample_specs(["ilg_6H"], self.facets)
        self.assertEqual(len(matches["ilg_6H"]), 2)

    def test_full_sample_name_matches_only_itself(self):
        matches, _ = expand_sample_specs(["20220902_RAW_ILG_6h_2_NEG"], self.facets)
        self.assertEqual(matches["20220902_RAW_ILG_6h_2_NEG"], ["20220902_RAW_ILG_6h_2_NEG"])

    def test_zero_match_raises_with_available_tokens(self):
        with self.assertRaises(ValueError) as ctx:
            expand_sample_specs(["24h"], self.facets)
        message = str(ctx.exception)
        self.assertIn("matched", message)
        self.assertIn("24h", message)
        self.assertIn("ilg", message)

    def test_blank_specs_are_skipped(self):
        matches, _ = expand_sample_specs(["ILG", "", "  "], self.facets)
        self.assertEqual(list(matches), ["ILG"])


class ExpandSampleSpecsRoleTests(unittest.TestCase):
    def setUp(self):
        self.facets = build_sample_facets([
            "20240311_Cerebellum_ICR_NEG_1",
            "20240311_Cerebellum_ICR_NEG_2",
            "20240311_QC_Cerebellum_ICR_NEG_1",
        ], None)

    def test_qc_is_excluded_by_default_and_reported(self):
        matches, excluded = expand_sample_specs(["cerebellum"], self.facets)
        self.assertEqual(len(matches["cerebellum"]), 2)
        self.assertEqual(excluded["cerebellum"], ["20240311_QC_Cerebellum_ICR_NEG_1"])

    def test_include_roles_can_bring_qc_back(self):
        matches, excluded = expand_sample_specs(
            ["cerebellum"], self.facets, include_roles=("sample", "qc"))
        self.assertEqual(len(matches["cerebellum"]), 3)
        self.assertEqual(excluded["cerebellum"], [])

    def test_include_roles_none_disables_role_filtering(self):
        matches, _ = expand_sample_specs(["cerebellum"], self.facets, include_roles=None)
        self.assertEqual(len(matches["cerebellum"]), 3)

    def test_all_hits_dropped_by_role_raises_and_names_include_roles(self):
        with self.assertRaises(ValueError) as ctx:
            expand_sample_specs(["qc"], self.facets)
        self.assertIn("include_roles", str(ctx.exception))


class AssignFactorGroupsTests(unittest.TestCase):
    def setUp(self):
        self.names = [
            "20220901_RAW_control_0h_1_NEG",
            "20220901_RAW_control_6h_1_NEG",
            "20220902_RAW_ILG_0h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_RAW_G_uralensis_6h_1_NEG",
        ]
        self.facets = build_sample_facets(self.names, name_class_index({
            "20220901_RAW_control_0h_1_NEG": "control",
            "20220901_RAW_control_6h_1_NEG": "control",
            "20220902_RAW_ILG_0h_1_NEG": "ILG",
            "20220902_RAW_ILG_6h_1_NEG": "ILG",
            "20220902_RAW_G_uralensis_6h_1_NEG": "G",
        }))

    def test_no_factors_falls_back_to_class_id(self):
        groups = assign_factor_groups(self.facets)
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "ILG")

    def test_single_axis_label_has_no_separator(self):
        # 1軸なら従来の group_levels と出力文字列が完全一致する（後方互換の要）。
        groups = assign_factor_groups(self.facets, group_levels=["0h", "6h"])
        self.assertEqual(groups["20220902_RAW_ILG_0h_1_NEG"], "0h")
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "6h")

    def test_two_axes_produce_cross_product_labels(self):
        groups = assign_factor_groups(self.facets, group_factors=[
            ["control", "ILG", "G_uralensis"],
            ["0h", "6h"],
        ])
        self.assertEqual(groups["20220901_RAW_control_0h_1_NEG"], "control|0h")
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "ILG|6h")
        self.assertEqual(groups["20220902_RAW_G_uralensis_6h_1_NEG"], "G_uralensis|6h")

    def test_custom_separator(self):
        groups = assign_factor_groups(
            self.facets, group_factors=[["ILG"], ["6h"]], sep=" / ")
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "ILG / 6h")

    def test_axis_without_hit_becomes_other(self):
        groups = assign_factor_groups(self.facets, group_factors=[["ILG"], ["24h"]])
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "ILG|other")

    def test_two_hits_within_one_axis_raises(self):
        with self.assertRaises(ValueError) as ctx:
            assign_factor_groups(self.facets, group_factors=[["ILG", "6h"]])
        self.assertIn("multiple", str(ctx.exception))

    def test_group_factors_wins_over_group_levels(self):
        groups = assign_factor_groups(
            self.facets, group_factors=[["0h"], ["ILG"]], group_levels=["control"])
        self.assertEqual(groups["20220902_RAW_ILG_0h_1_NEG"], "0h|ILG")

    def test_blank_axis_values_are_ignored(self):
        groups = assign_factor_groups(self.facets, group_factors=[["ILG", "", "  "]])
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "ILG")


class AssignFactorGroupsRoleTests(unittest.TestCase):
    def test_non_sample_roles_are_labeled_by_role_not_excluded(self):
        # QC は除外せず可視化する。QC の凝集は前処理品質の判断材料になるため。
        facets = build_sample_facets(
            ["20240311_ILG_6h_1_NEG", "20240311_QC_ILG_NEG_1"], None)
        groups = assign_factor_groups(facets, group_levels=["6h"])
        self.assertEqual(groups["20240311_ILG_6h_1_NEG"], "6h")
        self.assertEqual(groups["20240311_QC_ILG_NEG_1"], "qc")

    def test_role_labeling_does_not_apply_without_factors(self):
        facets = build_sample_facets(["20240311_QC_ILG_NEG_1"], None)
        self.assertIsNone(assign_factor_groups(facets)["20240311_QC_ILG_NEG_1"])


class TokenVocabularyTests(unittest.TestCase):
    def setUp(self):
        self.facets = build_sample_facets([
            "20220901_RAW_control_6h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
            "20220902_QC_RAW_NEG_1",
        ], None)
        self.vocab = token_vocabulary(self.facets)

    def test_counts_samples_per_token(self):
        self.assertEqual(self.vocab["tokens"]["raw"]["samples"], 3)
        self.assertEqual(self.vocab["tokens"]["6h"]["samples"], 2)
        self.assertEqual(self.vocab["tokens"]["ilg"]["samples"], 1)

    def test_breaks_down_by_role(self):
        self.assertEqual(self.vocab["tokens"]["raw"]["roles"], {"qc": 1, "sample": 2})
        self.assertEqual(self.vocab["tokens"]["ilg"]["roles"], {"sample": 1})

    def test_reports_positions_within_sample_name(self):
        self.assertEqual(self.vocab["tokens"]["raw"]["positions"], [1, 2])
        self.assertEqual(self.vocab["tokens"]["6h"]["positions"], [3])

    def test_by_position_lists_token_vocabulary(self):
        self.assertEqual(self.vocab["by_position"]["0"], ["20220901", "20220902"])
        self.assertEqual(set(self.vocab["by_position"]["2"]), {"control", "ilg", "raw"})

    def test_class_only_token_has_no_position(self):
        facets = build_sample_facets(["s1"], name_class_index({"s1": "Cerebellum_gf"}))
        vocab = token_vocabulary(facets)
        self.assertEqual(vocab["tokens"]["cerebellum"]["positions"], [])
        self.assertEqual(vocab["tokens"]["cerebellum"]["samples"], 1)

    def test_empty_facets(self):
        self.assertEqual(token_vocabulary({}), {"tokens": {}, "by_position": {}})


class FilterArfBySampleTokensTests(unittest.TestCase):
    def _features(self, names):
        return [{"MasterAlignmentID": 1,
                 "AlignedPeakProperties": [[i, n, 100.0 + i] for i, n in enumerate(names)]}]

    def setUp(self):
        self.names = [
            "20220901_RAW_control_0h_1_NEG",
            "20220901_RAW_control_6h_1_NEG",
            "20220902_RAW_ILG_0h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
        ]
        self.features = self._features(self.names)
        self.index = name_class_index({
            "20220901_RAW_control_0h_1_NEG": "control",
            "20220901_RAW_control_6h_1_NEG": "control",
            "20220902_RAW_ILG_0h_1_NEG": "ILG",
            "20220902_RAW_ILG_6h_1_NEG": "ILG",
        })

    def test_filters_by_token_absent_from_class_id(self):
        # 6h は Class ID に無い。これが通ることが本タスクの眼目。
        filtered, stats = filter_arf_by_class_ids(self.features, self.index, ["6h"])
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertEqual(kept, {"20220901_RAW_control_6h_1_NEG", "20220902_RAW_ILG_6h_1_NEG"})
        self.assertEqual(stats["after_sample_peaks"], 2)
        self.assertEqual(stats["before_sample_peaks"], 4)

    def test_multi_token_spec_and_or_across_specs(self):
        filtered, stats = filter_arf_by_class_ids(
            self.features, self.index, ["ILG_6h", "control_6h"])
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertEqual(kept, {"20220901_RAW_control_6h_1_NEG", "20220902_RAW_ILG_6h_1_NEG"})
        self.assertEqual(sorted(stats["matched_class_ids"]), ["ILG", "control"])
        self.assertEqual(len(stats["matched_samples"]), 2)

    def test_works_without_mddata(self):
        # class_index=None は従来 ValueError だった。名前トークンだけで成立させる。
        filtered, stats = filter_arf_by_class_ids(self.features, None, ["ILG"])
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertEqual(kept, {"20220902_RAW_ILG_0h_1_NEG", "20220902_RAW_ILG_6h_1_NEG"})
        self.assertEqual(stats["matched_class_ids"], [])
        self.assertEqual(stats["missing_samples"], 0)

    def test_qc_is_excluded_by_default_and_reported(self):
        names = self.names + ["20220901_QC_RAW_NEG_1"]
        filtered, stats = filter_arf_by_class_ids(self._features(names), None, ["raw"])
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertNotIn("20220901_QC_RAW_NEG_1", kept)
        self.assertEqual(stats["excluded_by_role"], ["20220901_QC_RAW_NEG_1"])

    def test_include_roles_can_bring_qc_back(self):
        names = self.names + ["20220901_QC_RAW_NEG_1"]
        filtered, stats = filter_arf_by_class_ids(
            self._features(names), None, ["raw"], include_roles=("sample", "qc"))
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertIn("20220901_QC_RAW_NEG_1", kept)
        self.assertEqual(stats["excluded_by_role"], [])


class FilterArfMissingClassPolicyTests(unittest.TestCase):
    """Class メタデータを解決できないサンプルの policy 分岐（統合トークン空間）。

    統合トークン空間では Class ID レコードが無いサンプルもサンプル名トークンで
    spec に一致するため、"exclude" が名前どおりに除外することを固定する。
    """

    def setUp(self):
        self.names = [
            "20220901_RAW_control_6h_1_NEG",
            "20220902_RAW_ILG_6h_1_NEG",
        ]
        self.features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [[i, n, 100.0 + i] for i, n in enumerate(self.names)],
        }]
        # ILG しか知らない class_index（control は Class ID レコードを持たない）。
        self.index = name_class_index({"20220902_RAW_ILG_6h_1_NEG": "ILG"})

    def test_error_policy_still_raises(self):
        with self.assertRaises(ValueError) as ctx:
            filter_arf_by_class_ids(self.features, self.index, ["6h"])
        self.assertIn("20220901_RAW_control_6h_1_NEG", str(ctx.exception))

    def test_exclude_policy_actually_drops_the_sample(self):
        filtered, stats = filter_arf_by_class_ids(
            self.features, self.index, ["6h"], missing_sample_policy="exclude")
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertEqual(kept, {"20220902_RAW_ILG_6h_1_NEG"})
        self.assertEqual(stats["matched_samples"], ["20220902_RAW_ILG_6h_1_NEG"])
        self.assertEqual(stats["after_sample_peaks"], 1)

    def test_exclude_policy_discloses_which_samples_were_dropped(self):
        _, stats = filter_arf_by_class_ids(
            self.features, self.index, ["6h"], missing_sample_policy="exclude")
        self.assertEqual(stats["missing_samples"], 1)
        self.assertEqual(
            stats["missing_samples_excluded"], ["20220901_RAW_control_6h_1_NEG"])


class AssignSampleGroupsFactorTests(unittest.TestCase):
    def test_group_factors_produce_cross_product(self):
        index = name_class_index({
            "20220901_RAW_control_0h_1_NEG": "control",
            "20220902_RAW_ILG_6h_1_NEG": "ILG",
        })
        groups = assign_sample_groups(
            ["20220901_RAW_control_0h_1_NEG", "20220902_RAW_ILG_6h_1_NEG"],
            index,
            group_factors=[["control", "ILG"], ["0h", "6h"]],
        )
        self.assertEqual(groups["20220901_RAW_control_0h_1_NEG"], "control|0h")
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "ILG|6h")

    def test_group_levels_work_without_class_index(self):
        groups = assign_sample_groups(
            ["20220902_RAW_ILG_6h_1_NEG"], None, group_levels=["6h"])
        self.assertEqual(groups["20220902_RAW_ILG_6h_1_NEG"], "6h")


if __name__ == "__main__":
    unittest.main()
