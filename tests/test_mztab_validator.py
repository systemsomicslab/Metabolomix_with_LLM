# tests/test_mztab_validator.py
from metabolomix.mztab.validator import validate_mztab, detect_quantification_measure

_BASE = {
    "metadata": {
        "mzTab-version": "2.0.0-M",
        "mzTab-mode": "Complete",
        "mzTab-type": "Quantification",
        "assay[1]-ms_run_ref": "ms_run[1]",
    },
    "sections": {
        "SMF": {"header": ["SMF_ID", "abundance_assay[1]"], "rows": [{"SMF_ID": "1", "abundance_assay[1]": "100.0"}], "warnings": []},
    },
    "warnings": [],
}


def _make(overrides=None, sections_overrides=None):
    import copy
    r = copy.deepcopy(_BASE)
    if overrides:
        r["metadata"].update(overrides)
    if sections_overrides:
        r["sections"].update(sections_overrides)
    return r


def test_valid_passes():
    result = validate_mztab(_BASE)
    assert result["ok"] is True
    assert result["errors"] == []


def test_wrong_version_fails():
    result = validate_mztab(_make({"mzTab-version": "2.1.0-M"}))
    assert result["ok"] is False
    assert any("2.0.0-M" in e for e in result["errors"])


def test_missing_smf_section_fails():
    import copy
    r = copy.deepcopy(_BASE)
    r["sections"].pop("SMF")
    result = validate_mztab(r)
    assert result["ok"] is False
    assert any("SMF" in e for e in result["errors"])


def test_missing_assay_fails():
    result = validate_mztab(_make({"assay[1]-ms_run_ref": None}))
    assert result["ok"] is False


def test_sme_trailing_warning_propagated():
    r = _make(sections_overrides={
        "SME": {"header": ["SME_ID"], "rows": [], "warnings": ["L10: trailing empty column removed"]}
    })
    result = validate_mztab(r)
    assert any("trailing" in w.lower() for w in result["warnings"])


def test_detect_height_from_filename():
    measure, conf = detect_quantification_measure(
        _BASE, "Height_AlignmentResult_2026.mzTab"
    )
    assert measure == "peak_height"
    assert conf in ("verified", "inferred")


def test_detect_area_from_filename():
    measure, conf = detect_quantification_measure(
        _BASE, "Area_AlignmentResult_2026.mzTab"
    )
    assert measure == "peak_area_above_zero"
    assert conf in ("verified", "inferred")


def test_detect_unknown_from_generic_filename():
    measure, conf = detect_quantification_measure(_BASE, "result.mzTab")
    assert measure is None or conf == "unknown"
