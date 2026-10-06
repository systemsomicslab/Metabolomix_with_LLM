import json
import unittest

import server
from metabolomix.core import session_state


def _row(file_id, name, height):
    """現実的な AlignmentChromPeakFeature 生 row（len>25・data[18] 数値）。
    exclusions.roster/_convert_to_alignment_feature が実パスで file_name を拾える。"""
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


class TestArfExclude(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        session_state.session.arf.filtered_features = _fixture()

    def test_requires_data(self):
        session_state.session.arf.filtered_features = None
        out = json.loads(server.arf_exclude(exclude_samples=["sA"]))
        self.assertEqual(out["error"]["code"], "missing_state")

    def test_add_sample_updates_set_and_counts(self):
        out = json.loads(server.arf_exclude(exclude_samples=["sB"]))
        self.assertEqual(out["status"], "success")
        self.assertIn("sB", out["excluded_samples"])
        self.assertEqual(out["samples_before"], 3)
        self.assertEqual(out["samples_after"], 2)
        self.assertEqual(session_state.session.arf.excluded_samples, {"sB"})

    def test_add_spot_updates_set_and_counts(self):
        out = json.loads(server.arf_exclude(exclude_spots=[1]))
        self.assertEqual(out["spots_before"], 2)
        self.assertEqual(out["spots_after"], 1)
        self.assertEqual(session_state.session.arf.excluded_spots, {1})

    def test_remove_re_includes(self):
        server.arf_exclude(exclude_samples=["sB"])
        out = json.loads(server.arf_exclude(exclude_samples=["sB"], mode="remove"))
        self.assertNotIn("sB", out["excluded_samples"])
        self.assertEqual(session_state.session.arf.excluded_samples, set())

    def test_clear_empties_all(self):
        server.arf_exclude(exclude_samples=["sB"], exclude_spots=[1])
        out = json.loads(server.arf_exclude(mode="clear"))
        self.assertEqual(out["excluded_samples"], [])
        self.assertEqual(out["excluded_spots"], [])
        self.assertEqual(session_state.session.arf.excluded_samples, set())
        self.assertEqual(session_state.session.arf.excluded_spots, set())

    def test_list_reports_without_change(self):
        server.arf_exclude(exclude_samples=["sB"])
        out = json.loads(server.arf_exclude(mode="list"))
        self.assertEqual(out["mode"], "list")
        self.assertEqual(out["excluded_samples"], ["sB"])
        self.assertEqual(session_state.session.arf.excluded_samples, {"sB"})

    def test_unmatched_are_reported_and_not_added(self):
        out = json.loads(server.arf_exclude(exclude_samples=["zzz"], exclude_spots=[999]))
        self.assertIn("zzz", out["unmatched_samples"])
        self.assertIn(999, out["unmatched_spots"])
        self.assertEqual(session_state.session.arf.excluded_samples, set())
        self.assertEqual(session_state.session.arf.excluded_spots, set())
        self.assertTrue(out["caveats"])


class TestSampleRolesExcludedFlag(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        session_state.session.arf.filtered_features = _fixture()
        session_state.session.arf.class_index = None

    def test_roles_mark_excluded_samples(self):
        session_state.session.arf.excluded_samples.add("sB")
        rows = _parse_tsv(server.arf_list_sample_roles())
        self.assertEqual(rows["sB"]["excluded"], "True")
        self.assertEqual(rows["sA"]["excluded"], "False")


def _parse_tsv(text: str) -> dict[str, dict[str, str]]:
    """`# ` で始まるヘッダ行を捨て、TSV 表を {sample: 行} に読み直す。"""
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    columns = lines[0].split("\t")
    return {
        dict(zip(columns, line.split("\t")))["sample"]: dict(zip(columns, line.split("\t")))
        for line in lines[1:]
    }


if __name__ == "__main__":
    unittest.main()
