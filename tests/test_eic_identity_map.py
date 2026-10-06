import unittest

from metabolomix.eic.identity_map import (
    MAX_CANDIDATES,
    select_identity_candidates,
    verify_spot_match,
)


def _record(alignment_id, name, ontology, rt, mz, height=0.0, adduct="[M+H]+"):
    return {
        "AlignmentID": alignment_id,
        "Name": name,
        "Ontology": ontology,
        "AdductType": adduct,
        "RT": rt,
        "MassCenter": mz,
        "HeightAverage": height,
    }


class SelectIdentityCandidatesTests(unittest.TestCase):
    def setUp(self):
        self.records = [
            _record(0, "PC(12:0/13:0)", "PC", 13.538, 636.4, height=10.0),
            _record(1, "TG(15:0/18:1/15:0) (d5)", "TG", 26.676, 800.7, height=30.0),
            _record(2, "Ceramide (d18:1/25:0)", "Cer", 22.148, 650.6, height=20.0),
            _record(3, "Unknown", "", 5.0, 300.0, height=5.0),
        ]

    def test_name_query_matches_case_insensitive_substring(self):
        candidates, caveats, total_matched = select_identity_candidates(
            self.records, names=["ceramide"],
        )
        self.assertEqual([c["spot_id"] for c in candidates], [2])
        self.assertEqual(candidates[0]["name"], "Ceramide (d18:1/25:0)")
        self.assertEqual(candidates[0]["ontology"], "Cer")
        self.assertAlmostEqual(candidates[0]["rt"], 22.148, places=3)
        self.assertEqual(caveats, [])
        self.assertEqual(total_matched, 1)

    def test_ontology_query_matches_exact_class_ignoring_case(self):
        candidates, _, total_matched = select_identity_candidates(
            self.records, ontologies=["tg"],
        )
        self.assertEqual([c["spot_id"] for c in candidates], [1])
        self.assertEqual(total_matched, 1)

    def test_name_and_ontology_results_are_a_deduplicated_union_sorted_by_spot_id(self):
        candidates, _, total_matched = select_identity_candidates(
            self.records, names=["TG(15:0"], ontologies=["TG", "PC"],
        )
        self.assertEqual([c["spot_id"] for c in candidates], [0, 1])
        self.assertEqual(total_matched, 2)

    def test_no_query_raises(self):
        with self.assertRaisesRegex(ValueError, "names"):
            select_identity_candidates(self.records, names=[], ontologies=None)

    def test_too_many_candidates_are_pruned_by_height_with_a_caveat(self):
        records = [
            _record(index, f"TG(x{index})", "TG", 10.0, 800.0, height=float(index))
            for index in range(5)
        ]
        candidates, caveats, total_matched = select_identity_candidates(
            records, ontologies=["TG"], max_candidates=2,
        )
        self.assertEqual([c["spot_id"] for c in candidates], [3, 4])
        self.assertEqual(len(caveats), 1)
        self.assertIn("HeightAverage", caveats[0])
        self.assertEqual(total_matched, 5)
        self.assertGreater(total_matched, len(candidates))

    def test_default_candidate_limit_is_three_hundred(self):
        self.assertEqual(MAX_CANDIDATES, 300)


class VerifySpotMatchTests(unittest.TestCase):
    def setUp(self):
        self.candidate = {
            "spot_id": 2, "name": "Ceramide (d18:1/25:0)", "ontology": "Cer",
            "adduct": "[M+H]+", "rt": 22.148, "mz": 650.6, "height_average": 1.0,
        }

    def test_matching_spot_returns_none(self):
        spot = {"spot_id": 2, "rt": 22.150, "mz": 650.605}
        self.assertIsNone(verify_spot_match(self.candidate, spot))

    def test_rt_outside_tolerance_returns_rt_mismatch(self):
        spot = {"spot_id": 2, "rt": 22.5, "mz": 650.6}
        self.assertEqual(verify_spot_match(self.candidate, spot), "rt_mismatch")

    def test_mz_outside_tolerance_returns_mz_mismatch(self):
        spot = {"spot_id": 2, "rt": 22.148, "mz": 651.0}
        self.assertEqual(verify_spot_match(self.candidate, spot), "mz_mismatch")


if __name__ == "__main__":
    unittest.main()
