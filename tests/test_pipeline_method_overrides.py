"""要求の `method_overrides`（実効メソッドの上書き）と、実効メソッドの主なキーを見せる経路を検証する。

要求の検査と正規化（`request`）→ 計画（`inputs.inspect_inputs`）→ 実効メソッドの書き出し
（`stage_inputs`）→ 受領・状態・読み込みが返す `unit_method_params` の順に見る。Console は動かさない。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from metabolomix.console import method_file as method_file_mod
from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import inputs as inputs_mod
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import service
from tests.config_fixtures import use_config
from tests.pipeline_fixtures import make_source


def _root(tmp_path, name="source"):
    root = tmp_path / name
    root.mkdir()
    return root


def _code(fn) -> DomainError:
    with pytest.raises(DomainError) as info:
        fn()
    return info.value


def _prepare_source(tmp_path, monkeypatch):
    source = make_source(tmp_path / "source")
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["exe"]})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    return source


# ---------- 要求 ----------

def test_single_polarity_request_normalizes_the_overrides(tmp_path):
    req = request_mod.resolve_request(_root(tmp_path), {
        "polarity": "positive",
        "method_overrides": {"minimum peak height": 50000.0, "N percent detected in one group": 50}})
    assert req["units"][0]["method_overrides"] == {
        "Minimum peak height": "50000", "N percent detected in one group": "50"}
    assert req["value_sources"]["units"][0]["method_overrides"] == "explicit"


def test_two_polarity_units_carry_their_own_overrides(tmp_path):
    root = _root(tmp_path)
    req = request_mod.resolve_request(root, {"units": [
        {"polarity": "positive", "raw_dir": "pos",
         "method_overrides": {"Minimum peak height": 30000}},
        {"polarity": "negative", "raw_dir": "neg"}]})
    assert req["units"][0]["method_overrides"] == {"Minimum peak height": "30000"}
    assert "method_overrides" not in req["units"][1]
    assert request_mod.unit_request(req, 1, root)["method_overrides"] == {}


def test_overrides_from_the_request_file_layer_with_explicit_keys(tmp_path):
    """1 極性の形はキーごとに重なる: ファイルの要求に上書きだけを明示で足せる。"""
    root = _root(tmp_path)
    (root / "analysis-request.json").write_text(
        json.dumps({"polarity": "negative", "library_mode": "msp_only"}), encoding="utf-8")
    req = request_mod.resolve_request(root, {"method_overrides": {"Minimum peak height": 20000}})
    assert req["library_mode"] == "msp_only"
    assert req["units"][0]["polarity"] == "negative"
    assert req["units"][0]["method_overrides"] == {"Minimum peak height": "20000"}


def test_empty_overrides_do_not_change_the_fingerprint(tmp_path):
    """`{}` は省いたのと同じ（既存の run の内容 hash を動かさない）。"""
    root = _root(tmp_path)
    plain = request_mod.resolve_request(root, {"polarity": "positive"})
    empty = request_mod.resolve_request(root, {"polarity": "positive", "method_overrides": {}})
    assert "method_overrides" not in empty["units"][0]
    assert request_mod.request_fingerprint(plain) == request_mod.request_fingerprint(empty)


def test_overrides_change_the_fingerprint_but_spelling_does_not(tmp_path):
    root = _root(tmp_path)
    plain = request_mod.resolve_request(root, {"polarity": "positive"})
    a = request_mod.resolve_request(root, {"polarity": "positive",
                                           "method_overrides": {"Minimum peak height": 50000}})
    b = request_mod.resolve_request(root, {"polarity": "positive",
                                           "method_overrides": {"minimum peak height": "50000"}})
    assert request_mod.request_fingerprint(a) == request_mod.request_fingerprint(b)
    assert request_mod.request_fingerprint(a) != request_mod.request_fingerprint(plain)


def test_bad_overrides_are_refused_at_the_request(tmp_path):
    root = _root(tmp_path)
    err = _code(lambda: request_mod.resolve_request(
        root, {"polarity": "positive", "method_overrides": {"Minimum peak height": "1e5"}}))
    assert err.code == "PIPELINE_REQUEST_INVALID"
    assert err.details["key"] == "Minimum peak height"
    err = _code(lambda: request_mod.resolve_request(
        root, {"polarity": "positive", "method_overrides": {"Minimum peek height": 1}}))
    assert "Minimum peak height" in err.details["suggestions"]
    err = _code(lambda: request_mod.resolve_request(
        root, {"polarity": "positive", "method_overrides": None}))
    assert err.code == "PIPELINE_REQUEST_INVALID"


def test_alignment_reference_override_follows_alignment_reference(tmp_path):
    root = _root(tmp_path)
    overrides = {"Alignment reference file ID": 2}
    assert _code(lambda: request_mod.resolve_request(
        root, {"polarity": "positive", "method_overrides": overrides})).details["use"] \
        == "alignment_reference"
    req = request_mod.resolve_request(root, {"polarity": "positive", "alignment_reference": "method",
                                             "method_overrides": overrides})
    assert req["units"][0]["method_overrides"] == {"Alignment reference file ID": "2"}


def test_resume_cannot_change_the_overrides(tmp_path):
    """上書きは上流の条件。変えるなら新しい run。"""
    req = request_mod.resolve_request(_root(tmp_path), {"polarity": "positive"})
    err = _code(lambda: request_mod.merge_updates(
        req, {"method_overrides": {"Minimum peak height": 1}}))
    assert err.code == "NEW_PIPELINE_REQUIRED"
    err = _code(lambda: request_mod.merge_updates(
        req, {"units": {"positive": {"method_overrides": {"Minimum peak height": 1}}}}))
    assert err.code == "NEW_PIPELINE_REQUIRED"


def test_stored_requests_without_overrides_still_validate(tmp_path):
    """このキーを足す前に保存した要求（単位に `method_overrides` が無い）も正規形として読める。"""
    req = request_mod.resolve_request(_root(tmp_path), {"polarity": "positive"})
    assert request_mod.validate_request(req, internal=True) is req


# ---------- 計画と実効メソッド ----------

def _plan(source, explicit):
    root = Path(source["root"])
    flat = request_mod.unit_request(request_mod.resolve_request(root, explicit), 0, root)
    return inputs_mod.inspect_inputs(root, flat, exe_path=source["exe"])


def test_overrides_reach_the_effective_method_and_the_original_is_untouched(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    method = source["method"]
    method.write_text(method.read_text(encoding="ascii")
                      + "Minimum peak height: 1000\nMass accuracy: 0.01\n",
                      encoding="ascii", newline="\n")
    original = method.read_bytes()

    plan = _plan(source, {"method_overrides": {
        "Minimum peak height": 50000, "MS1 tolerance for centroid": 0.005,
        "N percent detected in one group": 50}})

    # 原本の別名（Mass accuracy）も同じ値で差し替える。宣言の無かったキーは from が null。
    assert plan["method"]["user_overrides"] == {
        "MS1 tolerance for centroid": {"from": None, "to": "0.005"},
        "Mass accuracy": {"from": "0.01", "to": "0.005"},
        "Minimum peak height": {"from": "1000", "to": "50000"},
        "N percent detected in one group": {"from": None, "to": "50"},
    }
    assert plan["method"]["key_params"]["Minimum peak height"] == "50000"
    assert "Minimum peak height" not in plan["method"]["key_params_undeclared"]

    snapshot = inputs_mod.stage_inputs(plan, tmp_path / "unit")
    keys = method_file_mod.read_method_keys(tmp_path / "unit" / "inputs" / "effective-method.txt")
    assert keys["minimum peak height"] == "50000"
    assert keys["mass accuracy"] == "0.005"
    assert keys["ms1 tolerance for centroid"] == "0.005"
    assert keys["n percent detected in one group"] == "50"
    text = (tmp_path / "unit" / "inputs" / "effective-method.txt").read_text(encoding="ascii")
    assert text.count("Minimum peak height") == 1   # 原本の行を差し替える（後ろに足さない）
    assert method.read_bytes() == original

    params = receipt_mod.unit_method_params(tmp_path / "unit", snapshot)
    assert params["key_params"]["Minimum peak height"] == "50000"
    assert params["overrides"]["Minimum peak height"] == {"from": "1000", "to": "50000"}


def test_overrides_that_clash_with_what_the_plan_writes_are_refused(tmp_path, monkeypatch):
    """msp_then_lbm は LBM の照合で RT を使わない行を自分で書く。上書きと重なれば断る。"""
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, lbm=True)
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["exe"]})
    write_lab_msp_config(tmp_path, monkeypatch, lbm=True)

    err = _code(lambda: _plan(source, {"library_mode": "msp_then_lbm", "method_overrides": {
        "Use retention information for LBM-based annotation scoring": True}}))
    assert err.code == "PIPELINE_REQUEST_INVALID"
    assert err.details["library_mode"] == "msp_then_lbm"


def test_plan_reports_overrides_and_a_new_run_is_made_for_them(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    calls = []

    def fake_launch_detached(command, *, cwd, log_path):
        calls.append(command)
        return {"launched": True, "pid": 1000 + len(calls)}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.1)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    planned = service.plan_pipeline(source["root"], {
        "method_overrides": {"Minimum peak height": 50000}})
    method = planned["resolved"]["units"]["negative"]["method"]
    assert method["user_overrides"] == {"Minimum peak height": {"from": None, "to": "50000"}}
    assert method["key_params"]["Minimum peak height"] == "50000"

    first = service.start_pipeline(source["root"])
    second = service.start_pipeline(source["root"], {
        "method_overrides": {"Minimum peak height": 50000}})
    assert first["pipeline_path"] != second["pipeline_path"]
    # 計画だけの run（上と同じ要求）は起動のときに使い回す。
    assert second["pipeline_path"] == planned["pipeline_path"]


def test_unit_method_params_is_none_before_the_effective_method_exists(tmp_path):
    assert receipt_mod.unit_method_params(tmp_path, {"method": {}}) is None
    assert receipt_mod.unit_method_params(
        tmp_path, {"method": {"effective_relative_path": "inputs/effective-method.txt"}}) is None
