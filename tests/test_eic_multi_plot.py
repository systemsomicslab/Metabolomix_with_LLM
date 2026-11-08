import json
import unittest

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from metabolomix.plots.eic import build_multi_compound_plot_payload, render_eic_plot
from tests.config_fixtures import config_env


def _candidate(spot_id, name, ontology, rt, mz):
    return {
        "spot_id": spot_id, "name": name, "ontology": ontology,
        "adduct": "[M+H]+", "rt": rt, "mz": mz, "height_average": 1.0,
    }


def _spot(spot_id, rt, mz, file_id, points, peak_top, max_intensity):
    samples = []
    if file_id is not None:
        samples.append({
            "file_id": file_id,
            "peak_left": points[0][0],
            "peak_top": peak_top,
            "peak_right": points[-1][0],
            "num_points": len(points),
            "mean_intensity": sum(p[1] for p in points) / len(points),
            "max_intensity": max_intensity,
            "chromatogram": [list(point) for point in points],
        })
    return {
        "spot_id": spot_id, "rt": rt, "ri": 0.0, "mz": mz, "drift": -1.0,
        "main_type": 0, "num_samples": len(samples),
        "selected_samples": len(samples), "selected_points": len(points),
        "samples": samples,
    }


class MultiCompoundPayloadTests(unittest.TestCase):
    def setUp(self):
        self.candidates = [
            _candidate(0, "PC(12:0/13:0)", "PC", 13.5, 636.4),
            _candidate(1, "Ceramide (d18:1/25:0)", "Cer", 22.1, 650.6),
        ]
        self.spots = [
            _spot(0, 13.5, 636.4, 7, [(13.4, 2.0), (13.5, 40.0), (13.6, 3.0)], 13.5, 40.0),
            _spot(1, 22.1, 650.6, 7, [(22.0, 5.0), (22.1, 90.0), (22.2, 4.0)], 22.1, 90.0),
        ]

    def _build(self, **overrides):
        candidates = overrides.pop("candidates", self.candidates)
        spots = overrides.pop("spots", self.spots)
        params = {
            "file_id": 7,
            "file_path": "alignment.EIC.aef",
            "arf2_path": "alignment.arf2",
            "queries": ["PC(12:0", "Ceramide"],
            "ontologies": [],
        }
        params.update(overrides)
        return build_multi_compound_plot_payload(candidates, spots, **params)

    def test_series_are_labelled_by_lipid_name_and_sorted_by_rt(self):
        payload = self._build()
        self.assertEqual(payload["plot_schema"], "lipidmix.eic.multi.v1")
        self.assertEqual(payload["plot_type"], "line")
        self.assertEqual(payload["sample"]["file_id"], 7)
        self.assertEqual(
            [item["label"] for item in payload["series"]],
            ["PC(12:0/13:0)", "Ceramide (d18:1/25:0)"],
        )
        self.assertEqual([item["spot_id"] for item in payload["series"]], [0, 1])
        self.assertEqual(payload["axes"]["x"]["unit"], "min")
        self.assertEqual(payload["axes"]["y"]["label"], "Intensity")

    def test_annotation_sits_on_the_apex_data_point(self):
        payload = self._build()
        annotation = payload["series"][1]["annotation"]
        self.assertEqual(annotation["text"], "Ceramide (d18:1/25:0) / 22.100")
        self.assertAlmostEqual(annotation["x"], 22.1, places=5)
        self.assertAlmostEqual(annotation["y"], 90.0, places=5)

    def test_annotation_y_follows_normalization(self):
        payload = self._build(normalize="per_trace_max")
        self.assertEqual(payload["axes"]["y"]["label"], "Relative intensity")
        self.assertAlmostEqual(payload["series"][1]["annotation"]["y"], 1.0, places=5)

    def test_rt_mismatch_is_dropped_with_reason(self):
        candidates = [self.candidates[0], _candidate(1, "Ceramide", "Cer", 30.0, 650.6)]
        payload = build_multi_compound_plot_payload(
            candidates, self.spots, file_id=7,
            file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
        )
        self.assertEqual([item["spot_id"] for item in payload["series"]], [0])
        self.assertEqual(
            payload["selection"]["dropped"],
            [{"spot_id": 1, "name": "Ceramide", "reason": "rt_mismatch"}],
        )
        self.assertTrue(any("rt_mismatch" in note for note in payload["caveats"]))

    def test_missing_spot_and_missing_file_id_are_distinguished(self):
        candidates = self.candidates + [_candidate(2, "TG(x)", "TG", 25.0, 800.0)]
        spots = self.spots + [_spot(2, 25.0, 800.0, None, [(25.0, 1.0)], 25.0, 1.0)]
        payload = build_multi_compound_plot_payload(
            candidates + [_candidate(9, "DG(y)", "DG", 9.0, 500.0)],
            spots, file_id=7,
            file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
        )
        reasons = {item["spot_id"]: item["reason"] for item in payload["selection"]["dropped"]}
        self.assertEqual(reasons[2], "file_id_absent")
        self.assertEqual(reasons[9], "spot_out_of_range")

    def test_top_n_keeps_strongest_and_records_the_rest(self):
        payload = self._build(top_n=1)
        self.assertEqual([item["spot_id"] for item in payload["series"]], [1])
        self.assertEqual(
            payload["selection"]["dropped"],
            [{"spot_id": 0, "name": "PC(12:0/13:0)", "reason": "below_top_n"}],
        )
        self.assertEqual(payload["selection"]["plotted"], 1)
        self.assertEqual(payload["selection"]["candidates"], 2)
        self.assertEqual(payload["selection"]["candidates_evaluated"], 2)
        self.assertEqual(payload["selection"]["top_n"], 1)

    def test_candidates_defaults_to_evaluated_count_when_total_matched_omitted(self):
        payload = self._build()
        self.assertEqual(payload["selection"]["candidates"], 2)
        self.assertEqual(payload["selection"]["candidates_evaluated"], 2)

    def test_candidates_reports_pre_prune_total_when_pruning_happened(self):
        # select_identity_candidates prunes to max_candidates before the builder
        # sees the list; total_matched carries the pre-prune query-match count so
        # `candidates` keeps its documented meaning even though only 2 of the 900
        # matches were actually read and verified.
        payload = self._build(total_matched=900)
        self.assertEqual(payload["selection"]["candidates"], 900)
        self.assertEqual(payload["selection"]["candidates_evaluated"], 2)
        self.assertGreater(
            payload["selection"]["candidates"],
            payload["selection"]["candidates_evaluated"],
        )

    def test_caveats_aggregate_below_top_n_into_one_line_with_a_count(self):
        n = 30
        candidates = [
            _candidate(i, f"TG(x{i})", "TG", 10.0 + i * 0.001, 800.0)
            for i in range(n)
        ]
        spots = [
            _spot(
                i, 10.0 + i * 0.001, 800.0, 7,
                [(9.9, 1.0), (10.0 + i * 0.001, float(n - i)), (10.1, 1.0)],
                10.0 + i * 0.001, float(n - i),
            )
            for i in range(n)
        ]
        payload = build_multi_compound_plot_payload(
            candidates, spots, file_id=7,
            file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
            top_n=5,
        )
        below_top_n_notes = [
            note for note in payload["caveats"] if "below_top_n" in note
        ]
        self.assertEqual(len(below_top_n_notes), 1)
        dropped_count = sum(
            1 for item in payload["selection"]["dropped"]
            if item["reason"] == "below_top_n"
        )
        self.assertEqual(dropped_count, n - 5)
        self.assertIn(str(n - 5), below_top_n_notes[0])
        self.assertIn("selection.dropped", below_top_n_notes[0])

    def test_named_reasons_cap_at_ten_with_an_overflow_note(self):
        mismatched = [
            _candidate(i, f"PC(x{i})", "PC", 999.0, 636.4) for i in range(12)
        ]
        mismatched_spots = [
            _spot(
                i, 13.5, 636.4, 7,
                [(13.4, 1.0), (13.5, 5.0), (13.6, 1.0)], 13.5, 5.0,
            )
            for i in range(12)
        ]
        candidates = mismatched + [_candidate(100, "PC(ok)", "PC", 13.5, 636.4)]
        spots = mismatched_spots + [
            _spot(
                100, 13.5, 636.4, 7,
                [(13.4, 1.0), (13.5, 5.0), (13.6, 1.0)], 13.5, 5.0,
            )
        ]
        payload = build_multi_compound_plot_payload(
            candidates, spots, file_id=7,
            file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
        )
        rt_notes = [
            note for note in payload["caveats"] if note.startswith("rt_mismatch")
        ]
        self.assertEqual(len(rt_notes), 1)
        note = rt_notes[0]
        self.assertIn("12", note)
        self.assertIn("ほか 2 件", note)
        named_spots = sum(f"spot_id={i}" in note for i in range(12))
        self.assertEqual(named_spots, 10)
        dropped_rt_mismatch = [
            item for item in payload["selection"]["dropped"]
            if item["reason"] == "rt_mismatch"
        ]
        self.assertEqual(len(dropped_rt_mismatch), 12)

    def test_unknown_name_falls_back_to_spot_and_mz_label(self):
        candidates = [_candidate(0, "Unknown", "", 13.5, 636.4)]
        payload = build_multi_compound_plot_payload(
            candidates, self.spots[:1], file_id=7,
            file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
        )
        self.assertEqual(payload["series"][0]["label"], "spot 0 (m/z 636.4000)")

    def test_no_verified_candidate_raises(self):
        candidates = [_candidate(0, "PC", "PC", 99.0, 636.4)]
        with self.assertRaisesRegex(ValueError, "rt_mismatch"):
            build_multi_compound_plot_payload(
                candidates, self.spots[:1], file_id=7,
                file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
            )

    def test_invalid_arguments_raise(self):
        with self.assertRaisesRegex(ValueError, "normalize"):
            self._build(normalize="zscore")
        with self.assertRaisesRegex(ValueError, "top_n"):
            self._build(top_n=0)


class MultiCompoundRenderTests(unittest.TestCase):
    def setUp(self):
        candidates = [
            _candidate(0, "PC(12:0/13:0)", "PC", 13.5, 636.4),
            _candidate(1, "Ceramide (d18:1/25:0)", "Cer", 22.1, 650.6),
        ]
        spots = [
            _spot(0, 13.5, 636.4, 7, [(13.4, 2.0), (13.5, 40.0), (13.6, 3.0)], 13.5, 40.0),
            _spot(1, 22.1, 650.6, 7, [(22.0, 5.0), (22.1, 90.0), (22.2, 4.0)], 22.1, 90.0),
        ]
        self.payload = build_multi_compound_plot_payload(
            candidates, spots, file_id=7,
            file_path="alignment.EIC.aef", arf2_path="alignment.arf2",
        )

    def test_renders_one_line_and_one_annotation_per_series(self):
        fig = render_eic_plot(self.payload)
        try:
            ax = fig.axes[0]
            self.assertEqual(len(ax.get_lines()), 2)
            texts = [annotation.get_text() for annotation in ax.texts]
            self.assertIn("Ceramide (d18:1/25:0) / 22.100", texts)
            self.assertEqual(ax.get_xlabel(), "RT (min)")
        finally:
            plt.close(fig)

    def test_annotations_can_be_switched_off(self):
        self.payload["render_hints"]["show_annotations"] = False
        fig = render_eic_plot(self.payload)
        try:
            self.assertEqual(list(fig.axes[0].texts), [])
        finally:
            plt.close(fig)

    def test_unknown_schema_raises(self):
        with self.assertRaisesRegex(ValueError, "schema"):
            render_eic_plot({"plot_schema": "lipidmix.eic.v99", "series": []})


class PlotSourceShapeRegressionTests(unittest.TestCase):
    """Important 1: widening PlotSource must not leak into eic_plot_chromatograms.

    ``PlotSource`` must stay a strict ``file``/``file_name`` TypedDict; the
    ``arf2_file`` field belongs only to ``eic_plot_compounds`` via a separate
    ``MultiPlotSource`` TypedDict.

    以前は FastMCP が publish する outputSchema でこれを縛っていたが、payload の
    二重送出を止めるため outputSchema を出さなくなった（structured_output=False）。
    契約の実体は TypedDict そのものなので、そちらを直接検査する。
    """

    def test_single_spot_source_shape_is_unchanged(self):
        from metabolomix.plots.eic import PlotSource

        self.assertEqual(set(PlotSource.__annotations__), {"file", "file_name"})
        self.assertEqual(set(PlotSource.__required_keys__), {"file", "file_name"})

    def test_multi_compound_source_shape_has_all_three_fields(self):
        from metabolomix.plots.eic import MultiPlotSource

        self.assertEqual(
            set(MultiPlotSource.__annotations__), {"file", "file_name", "arf2_file"},
        )
        self.assertEqual(
            set(MultiPlotSource.__required_keys__), {"file", "file_name", "arf2_file"},
        )


if __name__ == "__main__":
    unittest.main()
