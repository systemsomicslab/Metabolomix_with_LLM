"""質量の窓の中心のずれ（ppm_offset）。帯の選び方・不正な形・新しい知見の候補の窓・rule_hash。"""
import pytest

from metabolomix.analysis.mass_accuracy import (normalize_novel_spec, novel_candidates,
                                                normalize_ppm_offset, offset_for, rule_hash)
from metabolomix.core.atomic_io import DomainError

#: ppm_offset を足す前（v8）の実装で計算した、最小の spec の rule_hash
_V8_HASH_OF_BASE = "37c0a80cd04a5a34eb32e842d128d880e1a737c641619f05c1fe770551a009c1"


class _Ds:
    def __init__(self, ppm, theo):
        self.feature_ids = ["1"]
        self.feature_metadata = {"1": {"rt": 2.0}}
        exp = theo * (1 + ppm / 1e6)
        self.feature_candidates = {"1": [{"rank": 1, "exp_mass_to_charge": exp,
                                          "theoretical_mass_to_charge": theo,
                                          "chemical_name": "x", "inchikey": "AAAAAAAAAAAAAA-B",
                                          "adduct": "[M+H]+"}]}


def _novel_ds(ppm, theo):
    stats = {"1": {"log2fc_test_over_ref": 1.0, "p": 0.001, "q": 0.001}}
    return _Ds(ppm, theo), stats


def test_offset_band_edge_goes_to_next_band():
    spec = normalize_ppm_offset({"default": 1.0, "bands": [{"max_mz": 200, "offset": -8.0}]},
                                code="NOVEL_SPEC_INVALID")
    assert offset_for(199.999, spec) == -8.0 and offset_for(200.0, spec) == 1.0
    assert offset_for(150.0, None) == 0.0


@pytest.mark.parametrize("bad", [{"bands": [{"max_mz": 200}]}, {"default": float("nan")},
                                 {"default": 0, "extra": 1}, "x"])
def test_bad_offset_names_the_field(bad):
    with pytest.raises(DomainError) as info:
        normalize_ppm_offset(bad, code="NOVEL_SPEC_INVALID")
    assert info.value.details["field"] == "ppm_offset"


def test_none_stays_none_and_bands_are_sorted():
    assert normalize_ppm_offset(None, code="X") is None
    got = normalize_ppm_offset({"bands": [{"max_mz": 400, "offset": 1}, {"max_mz": 200, "offset": 2}]},
                               code="X")
    assert [b["max_mz"] for b in got["bands"]] == [200.0, 400.0] and got["default"] == 0.0


def test_novel_candidate_inside_shifted_window_only():
    ds, stats = _novel_ds(ppm=-9.0, theo=150.0)                # 同定の誤差 -9 ppm
    base = {"statistics": {"result_id": "r"}, "ppm": {"default": 5.0}}
    plain = normalize_novel_spec(base)
    shifted = normalize_novel_spec({**base, "ppm_offset": {"bands": [{"max_mz": 200, "offset": -8.3}]}})
    assert novel_candidates(ds, plain, stats, group_ref="A", group_test="B") == []
    [row] = novel_candidates(ds, shifted, stats, group_ref="A", group_test="B")
    assert row["ppm"] == pytest.approx(-9.0)


def test_rule_hash_unchanged_without_offset_and_changed_with_it():
    base = {"statistics": {"result_id": "r"}}
    assert rule_hash(normalize_novel_spec(base)) == _V8_HASH_OF_BASE
    assert "ppm_offset" not in normalize_novel_spec(base)
    assert rule_hash(normalize_novel_spec({**base, "ppm_offset": {"default": -3}})) != _V8_HASH_OF_BASE
