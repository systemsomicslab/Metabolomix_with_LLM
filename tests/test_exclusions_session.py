import io
import unittest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state


class TestExclusionSessionState(unittest.TestCase):
    def setUp(self):
        session_state.session = session_state.AnalysisSession()

    def test_defaults_are_empty_sets(self):
        self.assertEqual(session_state.session.arf.excluded_samples, set())
        self.assertEqual(session_state.session.arf.excluded_spots, set())

    def test_new_file_load_resets_exclusions(self):
        s = session_state.session.arf
        s.excluded_samples.add("sA")
        s.excluded_spots.add(1)
        # 新ファイル読込（deserialize を差し替えてディスク非依存に）
        orig = arf_reader.deserialize
        arf_reader.deserialize = lambda buf: []
        self.addCleanup(setattr, arf_reader, "deserialize", orig)
        # discover_* は .arf 経路で呼ばれるので空データで無害に通す
        try:
            s.load_data("dummy_path.arf")
        except Exception:
            pass  # discover 系がファイル不在で落ちても、リセットは load_data 冒頭〜features 設定後
        self.assertEqual(s.excluded_samples, set())
        self.assertEqual(s.excluded_spots, set())


if __name__ == "__main__":
    unittest.main()
