"""`metabolomix.console.method_keys`（Console が読むキーの照合表と上書きの正規化）を検証する。"""
from __future__ import annotations

import pytest

from metabolomix.console import method_file as method_file_mod
from metabolomix.console import method_keys as mk


def _normalize(overrides, *, from_method=False):
    return mk.normalize_overrides(overrides, alignment_reference_from_method=from_method)


def _error(overrides, *, from_method=False) -> mk.MethodOverrideError:
    with pytest.raises(mk.MethodOverrideError) as info:
        _normalize(overrides, from_method=from_method)
    return info.value


def test_values_are_normalized_to_what_console_parses():
    """書き方が違っても同じ文字列になる（要求の内容 hash が揃う）。"""
    out = _normalize({
        "minimum peak height": 50000.0,
        "  N   percent detected in one group ": "50",
        "Force insert peaks in gap filling": "FALSE",
        "Retention time tolerance for alignment": "2e-1",
        "smoothing method": "linearweightedmovingaverage",
        "Solvent type": "hcoonh4",
    })
    assert out == {
        "Force insert peaks in gap filling": "False",
        "Minimum peak height": "50000",
        "N percent detected in one group": "50",
        "Retention time tolerance for alignment": "0.2",
        "Smoothing method": "LinearWeightedMovingAverage",
        "Solvent type": "HCOONH4",
    }
    assert list(out) == sorted(out)


@pytest.mark.parametrize("value", ["1e5", "50000.5", 50000.5, True, "", None, [1]])
def test_int_keys_refuse_values_that_int_tryparse_drops(value):
    """`Minimum peak height` は Console で int。通らない値は黙って既定に戻るので断る。"""
    err = _error({"Minimum peak height": value})
    assert err.details["key"] == "Minimum peak height"


def test_unknown_key_is_refused_with_close_spellings():
    err = _error({"Minimum peek height": 10000})
    assert "Minimum peak height" in err.details["suggestions"]
    assert err.details["upstream_commit"] == mk.UPSTREAM_COMMIT


def test_keys_console_cannot_read_are_refused_with_the_reason():
    """上流の case の綴りの誤りで、Console は `Sigma window value` を読まない。"""
    err = _error({"Sigma window value": 0.5})
    assert "読みません" in err.message


@pytest.mark.parametrize("key", ["Ion mode", "Msp file path", "Lbm file path", "Text DB file path",
                                 "MSP annotator settings file path"])
def test_keys_owned_by_other_request_fields_are_refused(key):
    err = _error({key: "x"})
    assert "use" in err.details


def test_alignment_reference_only_with_alignment_reference_method():
    assert _error({"Alignment reference file ID": 3}).details["use"] == "alignment_reference"
    assert _normalize({"Alignment reference file ID": 3}, from_method=True) == {
        "Alignment reference file ID": "3"}


def test_aliases_of_one_console_value_cannot_both_be_given():
    err = _error({"Mass accuracy": 0.01, "MS1 tolerance for centroid": 0.005})
    assert err.details["other"] == "Mass accuracy"
    assert mk.alias_keys("mass accuracy") == ("MS1 tolerance for centroid",)


def test_duplicate_key_with_different_case_is_refused():
    assert "2 度" in _error({"Minimum peak height": 1, "MINIMUM PEAK HEIGHT": 2}).message


def test_values_must_be_single_line_ascii():
    """改行を許すと、値の中から別のキーの行を差し込める。"""
    _error({"Searched adduct ions": "[M+H]+\nIon mode: Negative"})
    _error({"Searched adduct ions": "[M+H]+，[M+Na]+"})


def test_enum_values_are_checked():
    err = _error({"Blank filtering": "SampleAveOverBlankAve"})
    assert err.details["options"] == ["SampleMaxOverBlankAve"]


def test_non_dict_is_refused():
    _error(["Minimum peak height", 1])


def test_table_covers_the_keys_this_repository_writes():
    """計画とテンプレートが書くキーは、どれも Console が読む表にある（綴りの誤りの検出）。"""
    written = [method_file_mod.ION_MODE_KEY, method_file_mod.ADDUCT_KEY, method_file_mod.LBM_KEY,
               method_file_mod.MSP_KEY, method_file_mod.TEXT_DB_KEY,
               method_file_mod.RT_REFERENCE_KEY, method_file_mod.ALIGNMENT_REFERENCE_KEY,
               method_file_mod.LBM_RT_SCORING_KEY, method_file_mod.LBM_RT_FILTERING_KEY,
               *method_file_mod.CONSOLE_PATH_KEYS, *method_file_mod.SETTINGS_PATH_KEYS,
               *method_file_mod.KEY_PARAM_KEYS]
    missing = [key for key in written if key.lower() not in mk.CONSOLE_KEYS]
    assert missing == []


def test_every_enum_key_lists_its_values_and_reserved_keys_exist():
    for key, row in mk.CONSOLE_KEYS.items():
        if row.kind == "enum":
            assert key in mk._ENUM_VALUES, key
    assert set(mk.RESERVED_KEYS) <= set(mk.CONSOLE_KEYS)
    assert mk.ALIGNMENT_REFERENCE_KEY in mk.CONSOLE_KEYS
    assert not set(mk.UNREADABLE_KEYS) & set(mk.CONSOLE_KEYS)
