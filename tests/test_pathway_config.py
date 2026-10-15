import json

import pytest

from metabolomix.core import mcp_errors, user_config


def test_mslipidmapper_settings_read_env_first(monkeypatch, tmp_path):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[mslipidmapper]\nrscript = 'R/bin/Rscript.exe'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    got = user_config.get_setting("mslipidmapper.rscript")
    assert got.source == "config_file"
    assert got.value == str(tmp_path / "R" / "bin" / "Rscript.exe")
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", "X:/Rscript.exe")
    assert user_config.get_setting("mslipidmapper.rscript").value == "X:/Rscript.exe"


@pytest.mark.parametrize("key,env", [
    ("mslipidmapper.rscript", "MSLIPIDMAPPER_RSCRIPT"),
    ("mslipidmapper.browser", "MSLIPIDMAPPER_BROWSER"),
    ("mslipidmapper.api_url", "MSLIPIDMAPPER_API_URL"),
    ("mslipidmapper.api_token", "MSLIPIDMAPPER_API_TOKEN"),
])
def test_mslipidmapper_setting_env_names(key, env):
    assert user_config.SETTINGS[key] == env


def test_api_url_in_config_file_is_not_resolved_as_a_path(monkeypatch, tmp_path):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[mslipidmapper]\napi_url = 'http://10.0.0.1:8000/'\napi_token = ' t0k '\n",
                   encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    assert user_config.get_setting("mslipidmapper.api_url").value == "http://10.0.0.1:8000/"
    assert user_config.get_setting("mslipidmapper.api_token").value == "t0k"


def test_empty_api_url_is_unset(monkeypatch, tmp_path):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[mslipidmapper]\napi_url = ''\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    assert user_config.get_setting("mslipidmapper.api_url") is None


def test_pathway_error_envelope():
    body = json.loads(mcp_errors.pathway_error("MSLM_TIMEOUT", "遅い", {"timeout_s": 300}))
    assert body == {"error": {"code": "MSLM_TIMEOUT", "message": "遅い",
                              "details": {"timeout_s": 300}}}


def test_pathway_error_rejects_unknown_code():
    with pytest.raises(ValueError):
        mcp_errors.pathway_error("NOPE", "x")


def test_pathway_error_codes_are_the_spec_set():
    assert mcp_errors.PATHWAY_ERROR_CODES == frozenset({
        "MSLM_NOT_CONFIGURED", "MSLM_R_PACKAGES_MISSING", "MSLM_RUN_FAILED",
        "MSLM_TIMEOUT", "MSLM_OUTPUT_INVALID", "PATHWAY_RENDER_FAILED",
        "PATHWAY_BAD_REQUEST", "MSLM_REMOTE_UNAVAILABLE"})
