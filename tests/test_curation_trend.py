from lipidmix.curation import trend

TH = {"trend_min_points": 5, "trend_outlier_z": 3.0, "trend_min_r2": 0.7}


def test_composition_from_species_and_molecular_species():
    assert trend.composition("PC 34:1") == (34, 1)
    assert trend.composition("PC 16:0_18:1") == (34, 1)
    assert trend.composition("low score: PE O-38:5") == (38, 5)
    assert trend.composition("Cer 42:1;O2|Cer 18:1;O2/24:0") == (42, 1)


def test_composition_is_none_for_non_lipid_names():
    assert trend.composition("RIKEN N-VS1 ID-45 from Mouse") is None
    assert trend.composition("") is None
    assert trend.composition(None) is None


def _pc(spot_id, carbon, db, rt=None):
    rt = 2.0 + 0.5 * carbon - 0.8 * db if rt is None else rt
    return {"spot_id": spot_id, "ontology": "PC", "rt": rt, "mz": 400 + 14 * carbon - 2 * db,
            "carbon": carbon, "db": db}


def test_additive_model_flags_the_one_off_point():
    points = [_pc(i, c, d) for i, (c, d) in enumerate(
        [(32, 0), (34, 0), (36, 0), (34, 1), (36, 1), (38, 1), (36, 2), (38, 4)])]
    points.append(_pc(99, 34, 2, rt=30.0))
    result = trend.fit_trends(points, TH)
    assert result["classes"]["PC"]["n"] == 9
    assert result["spots"][99]["outlier"] is True
    assert result["spots"][99]["reliable"] is True
    assert result["spots"][0]["outlier"] is False


def test_class_with_too_few_points_is_not_fitted():
    result = trend.fit_trends([_pc(0, 34, 1), _pc(1, 36, 1)], TH)
    assert "PC" not in result["classes"]
    assert result["spots"] == {}


def test_groups_report_per_unsaturation_fit_quality():
    points = [_pc(i, c, 1) for i, c in enumerate([32, 34, 36, 38, 40])]
    group = trend.fit_trends(points, TH)["groups"]["PC"][1]
    assert group["n"] == 5
    assert group["r2"] > 0.99


def test_linear_inliers_without_outlier_flag_nothing():
    points = [_pc(i, c, d) for i, (c, d) in enumerate(
        [(32, 0), (34, 0), (36, 0), (34, 1), (36, 1), (38, 1), (36, 2), (38, 4)])]
    result = trend.fit_trends(points, TH)
    assert all(not v["outlier"] for v in result["spots"].values())
