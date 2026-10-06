import unittest
from metabolomix.arf import exclusions


def _row(file_id, name, height):
    """現実的な AlignmentChromPeakFeature 生 row（len>25・data[18] 数値）so that
    arf_reader._convert_to_alignment_feature exercises its real path and derives
    file_name from the string in data[:10]."""
    row = [0] * 26
    row[0] = file_id
    row[1] = name
    row[2] = file_id
    row[18] = float(height)
    return row


def _spot(master_id, rows):
    return {"MasterAlignmentID": master_id,
            "AlignedPeakProperties": [list(r) for r in rows]}


def _fixture():
    return [
        _spot(1, [_row(0, "sA", 10), _row(1, "sB", 20), _row(2, "sC", 30)]),
        _spot(2, [_row(0, "sA", 11), _row(1, "sB", 21), _row(2, "sC", 31)]),
    ]


class TestPruneSpots(unittest.TestCase):
    def test_empty_exclusions_returns_input_unchanged(self):
        spots = _fixture()
        out = exclusions.prune_spots(spots, set(), set())
        self.assertEqual(out, spots)

    def test_exclude_sample_drops_entry_from_all_spots(self):
        spots = _fixture()
        out = exclusions.prune_spots(spots, {"sB"}, set())
        for spot in out:
            names = {row[1] for row in spot["AlignedPeakProperties"]}
            self.assertNotIn("sB", names)
            self.assertEqual(names, {"sA", "sC"})

    def test_exclude_spot_drops_whole_spot(self):
        spots = _fixture()
        out = exclusions.prune_spots(spots, set(), {1})
        ids = {s["MasterAlignmentID"] for s in out}
        self.assertEqual(ids, {2})

    def test_exclude_both(self):
        spots = _fixture()
        out = exclusions.prune_spots(spots, {"sA"}, {2})
        self.assertEqual([s["MasterAlignmentID"] for s in out], [1])
        names = {row[1] for row in out[0]["AlignedPeakProperties"]}
        self.assertEqual(names, {"sB", "sC"})

    def test_non_destructive(self):
        spots = _fixture()
        exclusions.prune_spots(spots, {"sB"}, {1})
        # 元データは不変
        self.assertEqual(len(spots), 2)
        self.assertEqual(len(spots[0]["AlignedPeakProperties"]), 3)

    def test_unknown_names_are_noop(self):
        spots = _fixture()
        out = exclusions.prune_spots(spots, {"zzz"}, {999})
        self.assertEqual(len(out), 2)
        self.assertEqual(len(out[0]["AlignedPeakProperties"]), 3)


class TestRoster(unittest.TestCase):
    def test_roster_returns_names_and_ids(self):
        names, ids = exclusions.roster(_fixture())
        self.assertEqual(names, {"sA", "sB", "sC"})
        self.assertEqual(ids, {1, 2})


if __name__ == "__main__":
    unittest.main()
