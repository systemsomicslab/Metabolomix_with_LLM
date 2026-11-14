"""質量精度から勧める許容幅（recommended_tolerance。A3-T12）。"""
import pytest

from metabolomix.analysis.mass_accuracy import (MIN_N_FOR_TOLERANCE, band_summary,
                                                recommended_tolerance)


def test_needs_ten_identifications():
    out = recommended_tolerance([1.0] * (MIN_N_FOR_TOLERANCE - 1))
    assert out["mz_tol_ppm"] is None and out["reason"] == "too_few"
    assert out["n"] == 9 and out["min_n"] == 10
    assert recommended_tolerance([])["mz_tol_ppm"] is None


def test_width_is_abs_median_plus_three_raw_mad_rounded_up():
    values = [-3.0, -2.0, -2.0, -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 4.0]
    out = recommended_tolerance(values)
    # 中央値 -1.0、|偏差| = 2,1,1,0,0,0,1,1,2,5 → MAD 1.0（尺度の補正なし）→ 1 + 3 = 4
    assert out["median_ppm"] == pytest.approx(-1.0) and out["mad_ppm"] == pytest.approx(1.0)
    assert out["mz_tol_ppm"] == 4.0 and out["basis"] == "ceil(|median| + 3*MAD)"
    assert "reason" not in out
    shifted = recommended_tolerance([v + 0.3 for v in values])   # 0.7 + 3 = 3.7 → 4
    assert shifted["mz_tol_ppm"] == 4.0


def test_exact_integer_is_not_rounded_up_by_float_noise():
    # 中央値 0.1 と MAD 0.3 → 0.1 + 0.9 は浮動小数で 1.0000000000000002 になりうる。1 のまま
    values = [0.1 - 0.3] * 5 + [0.1] + [0.1 + 0.3] * 5
    assert recommended_tolerance(values)["mz_tol_ppm"] == 1.0


def test_never_below_one_ppm():
    assert recommended_tolerance([0.0] * 12)["mz_tol_ppm"] == 1.0


def test_band_summary_carries_the_recommendation_per_band_and_overall():
    rows = [{"theoretical_mz": 150.0, "ppm": -8.0}] * 10 + [{"theoretical_mz": 400.0, "ppm": 2.0}] * 3
    out = band_summary(rows, [200.0])
    low, high = out["bands"]
    assert low["recommended_mz_tol_ppm"] == 8.0
    assert high["recommended_mz_tol_ppm"] is None              # 3 件しかない
    assert out["recommended_tolerance"]["n"] == 13
    assert out["recommended_tolerance"]["mz_tol_ppm"] == 8.0   # 中央値 -8、MAD 0
