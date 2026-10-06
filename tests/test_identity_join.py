"""spot id と .arf2 同定情報の結合(純関数)。"""
import unittest

from metabolomix.arf import identity_join


def _catalog():
    return {
        1: {"MasterAlignmentID": 1, "Name": "PC 34:1", "Ontology": "PC",
            "InChIKey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "SMILES": "CCO",
            "MassCenter": 760.5851, "RT": 12.34},
        2: {"MasterAlignmentID": 2, "Name": "Unknown", "Ontology": "",
            "InChIKey": "", "SMILES": "", "MassCenter": 100.0, "RT": 1.0},
    }


class TestSpotIdOf(unittest.TestCase):
    def test_extracts_id(self):
        self.assertEqual(identity_join.spot_id_of("Spot_474_height"), 474)

    def test_returns_none_for_other_shapes(self):
        self.assertIsNone(identity_join.spot_id_of("f0"))


class TestJoinIdentity(unittest.TestCase):
    def test_keeps_only_rows_with_inchikey(self):
        results = [
            {"feature": "Spot_1_height", "mean_a": 10.0, "mean_b": 40.0,
             "log2fc": 1.9, "p": 0.001, "q": 0.01},
            {"feature": "Spot_2_height", "mean_a": 5.0, "mean_b": 5.0,
             "log2fc": 0.0, "p": 0.9, "q": 0.95},
        ]
        rows, report = identity_join.join_identity(results, _catalog())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["spot_id"], 1)
        self.assertEqual(rows[0]["inchikey"], "AAAAAAAAAAAAAA-BBBBBBBBBB-C")
        self.assertEqual(rows[0]["ontology"], "PC")
        self.assertAlmostEqual(rows[0]["mz"], 760.5851)
        self.assertEqual(report["n_features_total"], 2)
        self.assertEqual(report["n_with_inchikey"], 1)
        self.assertEqual(report["n_unannotated"], 1)

    def test_unknown_spot_id_counts_as_unannotated(self):
        results = [{"feature": "Spot_999_height", "mean_a": 1.0, "mean_b": 1.0,
                    "log2fc": 0.0, "p": 0.5, "q": 0.6}]
        rows, report = identity_join.join_identity(results, _catalog())
        self.assertEqual(rows, [])
        self.assertEqual(report["n_unannotated"], 1)

    def test_non_spot_feature_name_counts_as_unannotated(self):
        results = [{"feature": "f0", "mean_a": 1.0, "mean_b": 1.0,
                    "log2fc": 0.0, "p": 0.5, "q": 0.6}]
        rows, report = identity_join.join_identity(results, _catalog())
        self.assertEqual(rows, [])
        self.assertEqual(report["n_unannotated"], 1)
