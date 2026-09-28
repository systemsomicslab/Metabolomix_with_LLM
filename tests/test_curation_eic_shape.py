import math

from lipidmix.curation import eic_shape

TH = {"eic_min_points": 5, "eic_min_r2": 0.8, "eic_max_maxima": 2,
      "eic_pass_frac": 0.5, "eic_borderline_frac": 0.2, "eic_rt_scatter_sd": 0.1}


def gaussian(center=10.0, sigma=0.05, height=1000.0, n=41, step=0.01):
    start = center - (n // 2) * step
    return [[start + i * step, height * math.exp(-((start + i * step - center) ** 2) / (2 * sigma ** 2))]
            for i in range(n)]


def test_clean_gaussian_is_good():
    shape = eic_shape.sample_shape(gaussian(), left=9.85, top=10.0, right=10.15)
    assert shape["apex_in_window"] is True
    assert shape["gauss_r2"] > 0.95
    assert shape["n_maxima"] == 1
    assert eic_shape.sample_is_good(shape, TH)


def test_jagged_trace_has_many_maxima():
    trace = gaussian()
    for i in range(1, len(trace) - 1, 3):
        trace[i][1] *= 1.6
    shape = eic_shape.sample_shape(trace, left=9.85, top=10.0, right=10.15)
    assert shape["n_maxima"] > 2
    assert not eic_shape.sample_is_good(shape, TH)


def test_apex_at_the_window_edge_is_not_in_window():
    shape = eic_shape.sample_shape(gaussian(center=10.15), left=9.85, top=10.15, right=10.15)
    assert shape["apex_in_window"] is False


def test_empty_window_is_reported_not_raised():
    shape = eic_shape.sample_shape(gaussian(), left=20.0, top=20.1, right=20.2)
    assert shape["n_points"] == 0
    assert shape["gauss_r2"] is None


def test_spot_band_uses_only_detected_samples():
    good = {"file_id": 0, "detected": True, "chromatogram": gaussian(),
            "left": 9.85, "top": 10.0, "right": 10.15}
    flat = {"file_id": 1, "detected": False, "chromatogram": [[9.9, 1.0], [10.0, 1.0], [10.1, 1.0]],
            "left": 9.85, "top": 10.0, "right": 10.15}
    result = eic_shape.spot_shape([good, flat], TH)
    assert result["n_detected"] == 1
    assert result["good_fraction"] == 1.0
    assert result["band"] == "PASS"


def test_spot_without_detected_samples_is_unknown():
    flat = {"file_id": 1, "detected": False, "chromatogram": gaussian(),
            "left": 9.85, "top": 10.0, "right": 10.15}
    assert eic_shape.spot_shape([flat], TH)["band"] == "UNKNOWN"
