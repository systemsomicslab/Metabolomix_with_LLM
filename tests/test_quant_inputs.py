"""定量の入力表（参照表・付加体の採用表）。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from metabolomix.core import user_config
from metabolomix.core.atomic_io import DomainError
from metabolomix.quant import inputs

CSV = "TargetClass,StandardName,Coefficient,Unit\nPC,IS-A,2.0,nmol\n"


@pytest.fixture
def ref(tmp_path):
    path = tmp_path / "nested" / "ref-table.csv"
    path.parent.mkdir()
    path.write_bytes(CSV.encode())
    return path


def test_reference_is_required():
    with pytest.raises(DomainError) as e:
        inputs.reference_input(None)
    assert e.value.code == "QUANT_BAD_REQUEST"


def test_missing_reference_file(tmp_path):
    with pytest.raises(DomainError) as e:
        inputs.reference_input(str(tmp_path / "none.csv"))
    assert e.value.code == "QUANT_REFERENCE_NOT_FOUND"


def test_concentrations_need_volume_and_cells(ref):
    with pytest.raises(DomainError) as e:
        inputs.reference_input(str(ref), concentrations={"IS-A": 1.0}, solvent_volume=100)
    assert e.value.code == "QUANT_BAD_REQUEST"
    got = inputs.reference_input(str(ref), concentrations={"IS-A": 1}, solvent_volume=100,
                                 cells_millions=2.0, unit="nmol")
    assert got.concentrations == {"IS-A": 1.0}


def test_concentration_method_excludes_cell_inputs(ref):
    for extra in ({"concentrations": {"IS-A": 1.0}, "solvent_volume": 1, "cells_millions": 1},
                  {"solvent_volume": 1}):
        with pytest.raises(DomainError) as e:
            inputs.reference_input(str(ref), method="standard_concentration", **extra)
        assert e.value.code == "QUANT_BAD_REQUEST"
    got = inputs.reference_input(str(ref), method="standard_concentration")
    assert got.method == "standard_concentration"


def test_sha_and_label(ref):
    got = inputs.reference_input(str(ref))
    assert got.sha256 == hashlib.sha256(CSV.encode()).hexdigest()
    assert got.label == "ref-table.csv"
    assert "nested" not in got.label
    assert got.text == CSV


def test_adduct_table_from_config(tmp_path, monkeypatch):
    assert inputs.adduct_table() is None
    table = tmp_path / "adducts.csv"
    table.write_text("Class,Adduct,Ion mode,IsSelected\nPC,[M+H]+,Positive,True\n",
                     encoding="utf-8")
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text(f"[quant]\nadducts = '{table}'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    got = inputs.adduct_table()
    assert got.label == "adducts.csv" and "IsSelected" in got.text


def test_adduct_table_pointing_nowhere(tmp_path, monkeypatch):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[quant]\nadducts = 'gone.csv'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    with pytest.raises(DomainError) as e:
        inputs.adduct_table()
    assert e.value.code == "QUANT_REFERENCE_NOT_FOUND"


def test_quant_settings_are_config_file_keys(tmp_path, monkeypatch):
    # B-2 で設定は設定ファイルだけになった（環境変数は無い）。定量の 2 キーも普通の path の設定。
    assert user_config.SETTINGS["quant.rules"] == "path"
    assert user_config.SETTINGS["quant.adducts"] == "path"
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[quant]\nrules = 'x'\nadducts = 'y'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    s = user_config.get_setting("quant.rules")
    assert s.config_file == str(cfg) and Path(s.value) == (tmp_path / "x").resolve()
    assert "unknown_keys" not in user_config.describe_missing("msdial.exe")
    assert "環境変数" not in user_config.missing_hint("quant.rules", "規則")


def test_quant_error_matches_mztab_envelope():
    from metabolomix.core.mcp_errors import mztab_error
    assert inputs.quant_error("X", "m", {"a": 1}) == mztab_error("X", "m", {"a": 1})
    assert json.loads(inputs.quant_error("X", "m"))["error"] == {"code": "X", "message": "m"}



def test_non_utf8_tables_are_bad_request(tmp_path, monkeypatch):
    """cp932 など UTF-8 でない表は、生の UnicodeDecodeError でなく名前の末尾だけの封筒にする。"""
    path = tmp_path / "deep" / "ref-sjis.csv"
    path.parent.mkdir()
    path.write_bytes("TargetClass,StandardName,Coefficient,Unit\nPC,内部標準,2.0,nmol\n"
                     .encode("cp932"))
    with pytest.raises(DomainError) as e:
        inputs.reference_input(str(path))
    assert e.value.code == "QUANT_BAD_REQUEST"
    assert "ref-sjis.csv" in e.value.message and str(path.parent) not in e.value.message
    assert str(path.parent) not in json.dumps(e.value.details or {})
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text(f"[quant]\nadducts = '{path}'\n", encoding="utf-8")
    monkeypatch.setenv(user_config.CONFIG_ENV, str(cfg))
    with pytest.raises(DomainError) as e:
        inputs.adduct_table()
    assert e.value.code == "QUANT_BAD_REQUEST"
