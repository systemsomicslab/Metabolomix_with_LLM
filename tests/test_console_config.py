"""Console 経路が設定ファイルから実行体・LBM を引く（spec 2026-10-06 §6・§7）。"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from lipidmix.console.runner import MsdialExeNotFoundError, get_exe_path

BS = "\\"


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    path = tmp_path / "conf" / "lipidmix.local.toml"
    path.parent.mkdir()
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(path))
    return path


def _method(tmp_path: Path) -> Path:
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    return method


def _plan(tmp_path: Path) -> dict:
    from lipidmix.tools.console_tools import console_plan
    return json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(_method(tmp_path)),
                                   polarity="negative", measure="peak_height"))


# ---------- get_exe_path ----------

def test_get_exe_path_reads_the_config_file(cfg, tmp_path):
    exe = tmp_path / "MSDIALCUI.exe"
    exe.touch()
    cfg.write_text(f"[msdial]\nexe = '{exe}'\n", encoding="utf-8")
    assert get_exe_path() == str(exe)


def test_config_exe_pointing_at_a_missing_file_names_the_config_file(cfg, tmp_path):
    """設定ファイルの typo は NOT_FOUND で、出どころ（設定ファイルの [msdial] exe）を名指しする。"""
    exe = tmp_path / "no_such" / "MSDIALCUI.exe"
    cfg.write_text(f"[msdial]\nexe = '{exe}'\n", encoding="utf-8")
    with pytest.raises(MsdialExeNotFoundError) as exc:
        get_exe_path()
    assert exc.value.code == "MSDIAL_EXE_NOT_FOUND"
    assert "lipidmix.local.toml の [msdial] exe が指すファイルがありません" in str(exc.value)
    assert exc.value.details["source"] == "config_file"
    assert exc.value.details["setting"] == "msdial.exe"


def test_env_exe_is_not_checked_for_existence(cfg, monkeypatch):
    """環境変数の値は PATH 解決される素の名前でありうるので、存在確認しない。"""
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    assert get_exe_path() == "fake.exe"


def test_unset_exe_names_the_config_file_and_the_key(cfg):
    with pytest.raises(MsdialExeNotFoundError) as exc:
        get_exe_path()
    assert exc.value.code == "MSDIAL_EXE_NOT_FOUND"
    assert "lipidmix.example.toml" in str(exc.value)
    assert "[msdial] exe" in str(exc.value)
    assert exc.value.details["config_file"] == str(cfg)
    assert exc.value.details["config_file_exists"] is False


def test_broken_config_file_is_config_invalid(cfg):
    cfg.write_text("[msdial\n", encoding="utf-8")
    with pytest.raises(MsdialExeNotFoundError) as exc:
        get_exe_path()
    assert exc.value.code == "CONFIG_INVALID"
    assert exc.value.details["config_file"] == str(cfg)


# ---------- console_plan の封筒 ----------

def test_console_plan_envelope_carries_the_config_details(tmp_path, cfg):
    err = _plan(tmp_path)["error"]
    assert err["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert err["details"]["human_action_required"] is True
    assert err["details"]["setting"] == "msdial.exe"
    assert err["details"]["restart_required"] is False
    assert any("[msdial] exe" in step for step in err["details"]["how_to_set"])


def test_console_plan_with_a_missing_config_exe_is_not_found_not_not_console(tmp_path, cfg):
    exe = tmp_path / "no_such" / "MSDIALCUI.exe"
    cfg.write_text(f"[msdial]\nexe = '{exe}'\n", encoding="utf-8")
    err = _plan(tmp_path)["error"]
    assert err["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert "lipidmix.local.toml の [msdial] exe" in err["message"]
    assert err["details"]["setting"] == "msdial.exe"
    assert err["details"]["human_action_required"] is True


def test_console_plan_reports_config_invalid_with_the_quote_hint(tmp_path, cfg):
    cfg.write_text('[msdial]\nexe = "C:' + BS + 'MS-DIAL' + BS + 'x.exe"\n', encoding="utf-8")
    err = _plan(tmp_path)["error"]
    assert err["code"] == "CONFIG_INVALID"
    assert "単一引用符" in err["message"]


# ---------- console_method_template ----------

def _neg_param(path: Path) -> Path:
    path.write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\nLbm file path: \n"
        "Searched adduct ions: [M-H]-\n", encoding="ascii")
    return path


def test_template_uses_the_exe_from_the_config_file_to_find_the_lbm(tmp_path, cfg):
    """exe フォルダの *.lbm2 を 1 件だけ自動採用する GUI 流の推定が、設定ファイルの exe でも効く。"""
    from lipidmix.tools.console_tools import console_method_template
    app = tmp_path / "app"
    app.mkdir()
    (app / "lib.lbm2").touch()
    (app / "MSDIALCUI.exe").touch()
    cfg.write_text(f"[msdial]\nexe = '{app / 'MSDIALCUI.exe'}'\n", encoding="utf-8")
    parsed = json.loads(console_method_template(
        out_path=str(tmp_path / "param_POS.txt"), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["status"] == "written"
    assert parsed["lbm"]["source"] == "exe_dir"


def test_template_reports_config_invalid(tmp_path, cfg):
    from lipidmix.tools.console_tools import console_method_template
    cfg.write_text("[msdial\n", encoding="utf-8")
    parsed = json.loads(console_method_template(
        out_path=str(tmp_path / "param_POS.txt"), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["error"]["code"] == "CONFIG_INVALID"


# ---------- pipeline の受付 ----------

def test_pipeline_intake_resolves_the_exe_from_the_config_file(cfg, tmp_path):
    from lipidmix.pipeline import service
    exe = tmp_path / "MSDIALCUI.exe"
    exe.touch()
    cfg.write_text(f"[msdial]\nexe = '{exe}'\n", encoding="utf-8")
    assert service._resolve_exe_path() == exe


def test_pipeline_intake_wraps_config_invalid(cfg):
    from lipidmix.core.atomic_io import DomainError
    from lipidmix.pipeline import service
    cfg.write_text("[msdial\n", encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        service._resolve_exe_path()
    assert exc.value.code == "CONFIG_INVALID"


def test_pipeline_intake_keeps_the_not_found_code_with_details(cfg):
    from lipidmix.core.atomic_io import DomainError
    from lipidmix.pipeline import service
    with pytest.raises(DomainError) as exc:
        service._resolve_exe_path()
    assert exc.value.code == "MSDIAL_EXE_NOT_FOUND"
    assert exc.value.details["setting"] == "msdial.exe"


# ---------- LBM（[msdial] lbm） ----------

def test_template_uses_lbm_from_the_config_file_while_exe_comes_from_env(tmp_path, cfg, monkeypatch):
    """出どころが混ざる構成: exe は環境変数、LBM は設定ファイル。"""
    from lipidmix.tools.console_tools import console_method_template
    app = tmp_path / "app"
    app.mkdir()
    (app / "bundled.lbm2").touch()  # exe フォルダの LBM より設定ファイルの LBM が先
    lib = tmp_path / "libs" / "chosen.lbm2"
    lib.parent.mkdir()
    lib.touch()
    monkeypatch.setenv("MSDIAL_EXE", str(app / "MSDIALCUI.exe"))
    cfg.write_text(f"[msdial]\nlbm = '{lib}'\n", encoding="utf-8")
    out = tmp_path / "param_POS.txt"
    parsed = json.loads(console_method_template(
        out_path=str(out), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["lbm"]["source"] == "config_file"
    assert "chosen.lbm2" in out.read_text(encoding="ascii")


def test_template_reports_config_invalid_from_the_lbm_lookup(tmp_path, cfg, monkeypatch):
    """exe を環境変数で渡していても、LBM を引くときに壊れた設定ファイルに当たれば止める。"""
    from lipidmix.tools.console_tools import console_method_template
    app = tmp_path / "app"
    app.mkdir()
    monkeypatch.setenv("MSDIAL_EXE", str(app / "MSDIALCUI.exe"))
    cfg.write_text("[msdial\n", encoding="utf-8")
    parsed = json.loads(console_method_template(
        out_path=str(tmp_path / "param_POS.txt"), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["error"]["code"] == "CONFIG_INVALID"


def test_pipeline_pins_the_lbm_from_the_config_file(tmp_path, cfg):
    from lipidmix.pipeline import inputs
    lib = tmp_path / "chosen.lbm2"
    lib.write_bytes(b"lbm")
    cfg.write_text(f"[msdial]\nlbm = '{lib}'\n", encoding="utf-8")
    pinned = inputs._resolve_lbm_pinned({}, tmp_path / "param.txt",
                                        str(tmp_path / "app" / "MSDIALCUI.exe"), {}, tmp_path)
    assert pinned["source"] == "config_file"
    assert Path(pinned["path"]) == lib.resolve()


def test_pipeline_lbm_lookup_wraps_config_invalid(tmp_path, cfg):
    from lipidmix.core.atomic_io import DomainError
    from lipidmix.pipeline import inputs
    cfg.write_text("[msdial\n", encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        inputs._resolve_lbm_pinned({}, tmp_path / "param.txt",
                                   str(tmp_path / "app" / "MSDIALCUI.exe"), {}, tmp_path)
    assert exc.value.code == "CONFIG_INVALID"
