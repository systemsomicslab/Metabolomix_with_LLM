"""群別強度プロット（metabolomix/plots/group_intensity.py ほか）の純関数テスト。"""
from __future__ import annotations

from metabolomix.msdial import analysis_params


def test_param_file_min_peak_height_is_read(tmp_path):
    p = tmp_path / "Dataset_x_param_202610061200.txt"
    p.write_text("Ion mode: Negative\nMinimum peak height: 1000\n", encoding="utf-8")
    assert analysis_params.read_analysis_params(p)["min_peak_height"] == 1000.0


def test_param_file_without_min_peak_height_gives_none(tmp_path):
    p = tmp_path / "Dataset_x_param_202610061200.txt"
    p.write_text("Ion mode: Negative\n", encoding="utf-8")
    assert analysis_params.read_analysis_params(p)["min_peak_height"] is None
