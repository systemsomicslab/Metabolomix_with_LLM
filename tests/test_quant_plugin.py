"""規則プラグインの読み込み（偽のプラグインで）。"""
from __future__ import annotations

import sys

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.quant import plugin
from tests.quant_fake_rules import write_fake_rules

REF = ("TargetClass,StandardName,Coefficient,Unit\nPC,IS-A,2.0,nmol\n"
       "Any others,IS-B,1.0,nmol\n")
ADDUCTS = "Class,Adduct,Ion mode,IsSelected\nPC,[M+H]+,Positive,True\n"


def _point(tmp_path, monkeypatch, folder):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text(f"[quant]\nrules = '{folder}'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))


def _feat(i, name, cls, heights, detected=None):
    return {"feature_id": i, "name": name, "lipid_class": cls, "mz": 1.0, "rt": 1.0,
            "adduct": "[M+H]+", "heights": heights, "detected": detected or [True] * len(heights)}


def _sel(fid, name, role="target"):
    return {"feature_id": fid, "name": name, "comment": "", "msms_assigned": "True",
            "ontology": "PC", "adduct": "[M+H]+", "mode": "Positive", "role": role}


def test_rules_unconfigured():
    with pytest.raises(DomainError) as e:
        plugin.load_rules()
    assert e.value.code == "quant_rules_unavailable"
    assert "[quant]" in e.value.message and "rules" in e.value.message


def test_rules_folder_missing(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, tmp_path / "nope")
    with pytest.raises(DomainError) as e:
        plugin.load_rules()
    assert e.value.code == "quant_rules_unavailable"


def test_rules_import_failure_leaves_no_modules(tmp_path, monkeypatch):
    folder = write_fake_rules(tmp_path, extra="\nraise RuntimeError('boom')\n")
    _point(tmp_path, monkeypatch, folder)
    before = set(sys.modules)
    with pytest.raises(DomainError) as e:
        plugin.load_rules()
    assert e.value.code == "quant_rules_unavailable" and "boom" in e.value.message
    assert not [k for k in set(sys.modules) - before if "quant_rules" in k]


def test_rules_loaded_from_fake(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, write_fake_rules(tmp_path))
    p = plugin.load_rules()
    assert p.rules_version == "fake-rules-1" and p.label == "msdial_rules"
    assert "msdial_rules" not in sys.modules
    assert not any(k.startswith("msdial_rules.") for k in sys.modules)
    assert plugin.load_rules() is p  # 同じ指紋は使い回す
    assert plugin.rules_settings(p) == {"rules_version": "fake-rules-1", "digest": p.digest,
                                        "api_version": 1}
    assert f"_lipidmix_quant_rules_{p.digest[:8]}" in sys.modules


def test_incompatible_api_version(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, write_fake_rules(tmp_path, api_version=2))
    with pytest.raises(DomainError) as e:
        plugin.load_rules()
    assert e.value.code == "quant_rules_incompatible"


def test_digest_changes_with_code_and_reloads(tmp_path, monkeypatch):
    folder = write_fake_rules(tmp_path)
    _point(tmp_path, monkeypatch, folder)
    first = plugin.load_rules()
    api = folder / "api.py"
    api.write_text(api.read_text(encoding="utf-8") + "\n#\n", encoding="utf-8")
    second = plugin.load_rules()
    assert second.digest != first.digest and second.api is not first.api
    (folder / "__pycache__").mkdir()
    (folder / "__pycache__" / "x.py").write_text("junk", encoding="utf-8")
    assert plugin.load_rules().digest == second.digest


def test_rejected_rule_is_wrapped(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, write_fake_rules(tmp_path))
    p = plugin.load_rules()
    with pytest.raises(DomainError) as e:
        plugin.call(p, "normalize", features=[], n_samples=1, reference_text="bad")
    assert e.value.code == "QUANT_RULES_REJECTED"
    assert e.value.details["rule_code"] == "reference_invalid"


def test_fake_normalize(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, write_fake_rules(tmp_path))
    p = plugin.load_rules()
    feats = [_feat(1, "IS-A", "PC", [10.0, 0.0]), _feat(2, "IS-B", "PE", [5.0, 5.0]),
             _feat(3, "PC 34:1", "PC", [20.0, 20.0]), _feat(4, "PE 36:2", "PE", [10.0, None]),
             _feat(5, "TG 50:0", "TG", [1.0, 1.0])]
    out = plugin.call(p, "normalize", features=feats, n_samples=2, reference_text=REF)
    assert out["values"][2] == [4.0, None]
    assert out["statuses"][2] == ["ok", "standard_invalid_height"]
    assert out["values"][3] == [2.0, None] and out["statuses"][3][1] == "missing"
    assert sorted(out["any_others_classes"]) == ["PE", "TG"]
    assert out["roles"][:2] == ["internal_standard"] * 2
    out = plugin.call(p, "normalize", features=feats[2:], n_samples=2, reference_text=REF)
    assert out["assignments"]["PC"]["status"] == "standard_missing"
    twin = feats + [_feat(6, "IS-A", "PC", [4.0, 4.0])]
    out = plugin.call(p, "normalize", features=twin, n_samples=2, reference_text=REF)
    assert out["assignments"]["PC"]["status"] == "standard_ambiguous"
    out = plugin.call(p, "normalize", features=twin, n_samples=2, reference_text=REF,
                      bindings={"PC": 6})
    assert out["values"][2] == [10.0, 10.0]
    std = ("TargetClass\tStandardName\tStandardAdduct\nPC\tIS-A\t[M+H]+\n"
           "Any others\tIS-B\t[M+H]+\n")
    text = plugin.call(p, "coefficient_reference_text", reference_text=std,
                       concentrations={"IS-A": 2.0, "IS-B": 1.0}, solvent_volume=10.0,
                       cells_millions=4.0)
    lines = text.splitlines()
    assert lines[0] == "TargetClass,StandardName,Coefficient,Unit,StandardAdduct"
    assert lines[1] == "PC,IS-A,5.0,pmol / 10^6 cells,[M+H]+"
    assert lines[2] == "Any others,IS-B,2.5,pmol / 10^6 cells,[M+H]+"
    out = plugin.call(p, "normalize", features=feats, n_samples=2, reference_text=text)
    assert out["values"][2] == [10.0, None]  # 20/10*5; second sample IS height 0
    assert out["values"][3][0] == 5.0  # 10/5*2.5 via Any others
    with pytest.raises(DomainError) as e:
        plugin.call(p, "coefficient_reference_text", reference_text=std,
                    concentrations={"IS-A": 2.0}, solvent_volume=10.0, cells_millions=4.0)
    assert e.value.details["rule_code"] == "concentration_missing"


def test_fake_match_and_select(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, write_fake_rules(tmp_path))
    p = plugin.load_rules()
    m = plugin.call(p, "match_samples",
                    names_by_mode={"pos": ["s1_pos", "s2_pos"], "neg": ["s1_neg", "s3"]})
    assert m["ids_by_mode"] == {"pos": ["s1"], "neg": ["s1"]}
    assert m["mismatched"]["pos"] == ["s2"]
    sel = plugin.call(p, "select", adduct_text=ADDUCTS,
                      features=[_sel("a", "PC 34:1"), _sel("b", "Unknown"),
                                _sel("c", "IS-A", "internal_standard")])
    assert [s["action"] for s in sel] == ["keep", "exclude", "exclude"]


def test_fake_match_samples_with_mapping(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, write_fake_rules(tmp_path))
    p = plugin.load_rules()
    rows = [{"mode": "pos", "source_sample": "a_pos", "sample_id": "S1", "group": "g", "sample_type": "x"},
            {"mode": "pos", "source_sample": "b_pos", "sample_id": "S2", "group": "g", "sample_type": "x"},
            {"mode": "neg", "source_sample": "zz", "sample_id": "S1", "group": "g", "sample_type": "x"},
            {"mode": "neg", "source_sample": "yy", "sample_id": "S3", "group": "g", "sample_type": "x"}]
    names = {"pos": ["b_pos", "a_pos"], "neg": ["zz", "yy"]}
    m = plugin.call(p, "match_samples", names_by_mode=names, mapping_rows=rows)
    assert m["rule"] == "mapping_csv"
    assert m["ids_by_mode"] == {"pos": ["S1"], "neg": ["S1"]}
    assert m["mismatched"] == {"pos": ["S2"], "neg": ["S3"]}
    assert m["rows_by_mode"]["pos"]["a_pos"]["sample_id"] == "S1"
    # the mapping resolves a pair the suffix rule cannot
    assert plugin.call(p, "match_samples", names_by_mode=names)["ids_by_mode"]["neg"] == []
    with pytest.raises(DomainError) as e:
        plugin.call(p, "match_samples", names_by_mode={"pos": ["c_pos"]}, mapping_rows=rows)
    assert e.value.details["rule_code"] == "mapping_invalid"
