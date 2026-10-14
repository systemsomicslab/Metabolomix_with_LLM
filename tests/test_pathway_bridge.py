# tests/test_pathway_bridge.py
import json
import sys
import textwrap

import pytest

from metabolomix.pathway import bridge

pytestmark = pytest.mark.skipif(sys.platform != "win32", reason="start_owned_process は Windows 専用")


def _fake(tmp_path, body: str):
    """偽ブリッジ: Python スクリプト。argv[1] が renv プロジェクト、argv[2] が request.json。"""
    script = tmp_path / "fake_bridge.py"
    script.write_text(textwrap.dedent("""
        import json, sys, time, pathlib
        assert pathlib.Path(sys.argv[1]).name == "r"
        req = json.loads(pathlib.Path(sys.argv[2]).read_text(encoding="utf-8"))
        out = pathlib.Path(req["out_dir"])
    """) + textwrap.dedent(body), encoding="utf-8")
    return script


@pytest.fixture
def configured(monkeypatch, tmp_path):
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", sys.executable)
    renv = tmp_path / "r"
    (renv / "renv").mkdir(parents=True)
    (renv / "renv" / "activate.R").write_text("# fake", encoding="utf-8")
    monkeypatch.setattr(bridge, "RENV_PROJECT", renv)
    return tmp_path


def test_missing_renv_environment(configured, monkeypatch):
    monkeypatch.setattr(bridge, "RENV_PROJECT", configured / "no_r")
    run_dir = configured / "run"
    run_dir.mkdir()
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_R_PACKAGES_MISSING"
    assert "r/setup.R" in exc.value.message


def test_not_configured(monkeypatch):
    monkeypatch.delenv("MSLIPIDMAPPER_RSCRIPT", raising=False)
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.resolve_rscript()
    assert exc.value.code == "MSLM_NOT_CONFIGURED"


def test_configured_but_missing_file(monkeypatch, tmp_path):
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", str(tmp_path / "nope.exe"))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.resolve_rscript()
    assert exc.value.code == "MSLM_NOT_CONFIGURED"


def test_ok_result_is_returned(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        (out / "bridge-result.json").write_text(json.dumps(
            {"status": "ok", "code": None, "message": "", "outputs": {"remodeling": "x.json"},
             "mapped_nodes": {"remodeling": 3}, "package_version": "0.1"}), encoding="utf-8")
    """))
    res = bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert res["mapped_nodes"] == {"remodeling": 3}
    assert json.loads((run_dir / "request.json").read_text(encoding="utf-8"))["mode"] == "measure"


def test_error_result_carries_code(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        (out / "bridge-result.json").write_text(json.dumps(
            {"status": "error", "code": "MSLM_R_PACKAGES_MISSING",
             "message": "MSLipidMapper が入っていません"}), encoding="utf-8")
        sys.exit(2)
    """))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_R_PACKAGES_MISSING"


def test_nonzero_without_result_is_run_failed_with_log_tail(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        for i in range(40):
            print(f"line {i}", file=sys.stderr)
        sys.exit(1)
    """))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_RUN_FAILED"
    tail = exc.value.details["log_tail"]
    assert tail.splitlines()[-1] == "line 39" and len(tail.splitlines()) == 20
    assert (run_dir / "bridge.log").is_file()


def test_timeout_kills_tree(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        import subprocess
        subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
        time.sleep(60)
    """))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=2)
    assert exc.value.code == "MSLM_TIMEOUT"


def test_zero_exit_without_result_is_output_invalid(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, "pass\n"))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_OUTPUT_INVALID"
