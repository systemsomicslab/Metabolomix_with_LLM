"""volcano プロットの構造化ペイロード（lipidmix.volcano.v1）のテスト。

間引きは決定的で、up/down は絶対に切らないことを固定する。ns の削減件数が
selection と caveats の両方に出ることも検証する（「間引かれた点＝有意でない」と
誤読されないための表示根拠）。
"""
import json
import math
import unittest

from metabolomix.plots import volcano as volcano_plot
from tests.config_fixtures import config_env


def _point(feature, log2fc, neg_log10_p, sig):
    return {"feature": feature, "log2fc": log2fc,
            "neg_log10_p": neg_log10_p, "sig": sig}


def _differential(volcano, **overrides):
    base = {"kind": "two_group", "a": "24M", "b": "9w", "n_a": 6, "n_b": 5,
            "q_threshold": 0.05, "log2fc_threshold": 1.0, "volcano": volcano}
    base.update(overrides)
    return base


class TestBuildVolcanoPlotPayload(unittest.TestCase):
    def test_payload_shape_and_metadata(self):
        payload = volcano_plot.build_volcano_plot_payload(_differential([
            _point("PC 34:1", 1.8, 3.4, "up"),
            _point("PE 36:2", -2.1, 4.0, "down"),
            _point("TG 52:3", 0.2, 0.5, "ns"),
        ]))
        self.assertEqual(payload["plot_schema"], "lipidmix.volcano.v1")
        self.assertEqual(payload["plot_type"], "scatter")
        self.assertEqual(payload["title"], "Volcano (24M vs 9w)")
        self.assertEqual(payload["comparison"],
                         {"group_a": "24M", "group_b": "9w", "n_a": 6, "n_b": 5})
        self.assertEqual(payload["axes"]["x"]["label"], "log2 fold change")
        self.assertEqual(payload["axes"]["y"]["label"], "-log10 p")
        self.assertEqual(payload["thresholds"], {"q": 0.05, "log2fc": 1.0})
        self.assertEqual(payload["render_hints"]["mode"], "markers")
        self.assertEqual(payload["render_hints"]["color_by"], "sig")
        self.assertEqual(payload["render_hints"]["guides"]["x"], [-1.0, 1.0])
        self.assertAlmostEqual(payload["render_hints"]["guides"]["y"][0], 1.3010, places=3)
        self.assertEqual(len(payload["points"]), 3)
        self.assertEqual(payload["selection"]["total"], 3)
        self.assertEqual(payload["selection"]["plotted"], 3)

    def test_custom_title_overrides_default(self):
        payload = volcano_plot.build_volcano_plot_payload(
            _differential([_point("PC 34:1", 1.8, 3.4, "up")]), title="肝 24M vs 9w",
        )
        self.assertEqual(payload["title"], "肝 24M vs 9w")

    def test_significant_points_are_never_thinned(self):
        volcano = [_point(f"sig-{i}", 2.0, 5.0, "up") for i in range(40)]
        volcano += [_point(f"sig-d-{i}", -2.0, 5.0, "down") for i in range(40)]
        volcano += [_point(f"ns-{i}", 0.1, 0.2, "ns") for i in range(500)]
        payload = volcano_plot.build_volcano_plot_payload(
            _differential(volcano), max_points=100,
        )
        kept = [p for p in payload["points"] if p["sig"] in ("up", "down")]
        self.assertEqual(len(kept), 80)
        self.assertEqual(payload["selection"]["significant_total"], 80)
        self.assertEqual(payload["selection"]["significant_plotted"], 80)
        self.assertEqual(payload["selection"]["ns_plotted"], 20)
        self.assertEqual(payload["selection"]["plotted"], 100)

    def test_ns_subsampling_is_deterministic(self):
        volcano = [_point(f"ns-{i}", 0.1, 0.2, "ns") for i in range(1000)]
        first = volcano_plot.build_volcano_plot_payload(
            _differential(volcano), max_points=50)
        second = volcano_plot.build_volcano_plot_payload(
            _differential(volcano), max_points=50)
        self.assertEqual([p["feature"] for p in first["points"]],
                         [p["feature"] for p in second["points"]])
        self.assertEqual(len(first["points"]), 50)

    def test_ns_thinning_is_reported_in_caveats(self):
        volcano = [_point(f"ns-{i}", 0.1, 0.2, "ns") for i in range(300)]
        payload = volcano_plot.build_volcano_plot_payload(
            _differential(volcano), max_points=50)
        self.assertTrue(any("間引" in note for note in payload["caveats"]))

    def test_nonfinite_points_are_dropped_and_explained(self):
        payload = volcano_plot.build_volcano_plot_payload(_differential([
            _point("ok", 1.8, 3.4, "up"),
            _point("no-fc", None, 3.4, "ns"),
            _point("nan-fc", math.nan, 3.4, "ns"),
            _point("nan-p", 1.0, math.nan, "ns"),
        ]))
        self.assertEqual([p["feature"] for p in payload["points"]], ["ok"])
        self.assertEqual(payload["selection"]["dropped_nonfinite"], 3)
        self.assertEqual(payload["selection"]["total"], 4)
        self.assertTrue(any("有限値でない" in note for note in payload["caveats"]))

    def test_significant_overflow_drops_all_ns_but_keeps_significant(self):
        volcano = [_point(f"sig-{i}", 2.0, 5.0, "up") for i in range(120)]
        volcano += [_point(f"ns-{i}", 0.1, 0.2, "ns") for i in range(30)]
        payload = volcano_plot.build_volcano_plot_payload(
            _differential(volcano), max_points=100)
        self.assertEqual(payload["selection"]["significant_plotted"], 120)
        self.assertEqual(payload["selection"]["ns_plotted"], 0)
        self.assertEqual(payload["selection"]["ns_total"], 30)
        self.assertEqual(payload["selection"]["plotted"], 120)
        self.assertTrue(any("全て省略" in note for note in payload["caveats"]))

    def test_q_and_p_axis_mismatch_is_always_caveated(self):
        payload = volcano_plot.build_volcano_plot_payload(
            _differential([_point("PC 34:1", 1.8, 3.4, "up")]))
        self.assertTrue(any("-log10(q" in note for note in payload["caveats"]))

    def test_unadjusted_confounded_status_is_carried_from_provenance(self):
        """spec 7.4(R16): allow_confounded=trueで継続した比較は、payloadのcaveatsにも
        その旨が要る——run_comparisonが`provenance.comparison.unadjusted_confounded`へ
        残す事実を、build_volcano_plot_payload自身の独立したcaveatsリストが
        これまで一度も読んでいなかった（Task17のcontroller ruling R16）。"""
        last = _differential(
            [_point("PC 34:1", 1.8, 3.4, "up")],
            provenance={"comparison": {"unadjusted_confounded": True}})
        payload = volcano_plot.build_volcano_plot_payload(last)
        self.assertTrue(any("未調整" in note for note in payload["caveats"]))

    def test_adjusted_comparisons_get_no_confounding_caveat(self):
        last = _differential(
            [_point("PC 34:1", 1.8, 3.4, "up")],
            provenance={"comparison": {"unadjusted_confounded": False}})
        payload = volcano_plot.build_volcano_plot_payload(last)
        self.assertFalse(any("未調整" in note for note in payload["caveats"]))

    def test_missing_thresholds_fall_back_to_defaults(self):
        last = {"kind": "two_group", "a": "A", "b": "B",
                "volcano": [_point("PC 34:1", 1.8, 3.4, "up")]}
        payload = volcano_plot.build_volcano_plot_payload(last)
        self.assertEqual(payload["thresholds"], {"q": 0.05, "log2fc": 1.0})
        self.assertIsNone(payload["comparison"]["n_a"])

    def test_rejects_invalid_max_points(self):
        last = _differential([_point("PC 34:1", 1.8, 3.4, "up")])
        for bad in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                volcano_plot.build_volcano_plot_payload(last, max_points=bad)


if __name__ == "__main__":
    unittest.main()
