# tests/test_pathway_bridge.py
import hashlib
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
    _restore_library(renv)
    monkeypatch.setattr(bridge, "RENV_PROJECT", renv)
    return tmp_path


def _restore_library(renv):
    """`Rscript r/setup.R` の後の状態を模す。

    renv::restore の配置（library/<platform>/R-x.y/<arch>/<pkg>/DESCRIPTION）に加え、
    setup.R が復元に使った renv.lock の md5 を目印として残す。
    """
    pkg = renv / "renv" / "library" / "windows" / "R-4.6" / "x86_64-w64-mingw32" / "MSLipidMapper"
    pkg.mkdir(parents=True)
    (pkg / "DESCRIPTION").write_text("Package: MSLipidMapper\n", encoding="utf-8")
    lock = renv / "renv.lock"
    lock.write_text('{"R": {"Version": "4.6.0"}}\n', encoding="utf-8")
    (renv / bridge.RESTORE_STAMP).write_text(
        hashlib.md5(lock.read_bytes()).hexdigest() + "\n", encoding="utf-8")


def test_library_without_stamp_is_stale(configured, monkeypatch):
    """目印の無いライブラリ（目印導入前の復元・手で入れたもの）は lock との一致を保証できない。"""
    (bridge.RENV_PROJECT / bridge.RESTORE_STAMP).unlink()
    run_dir = configured / "run"
    run_dir.mkdir()
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_R_PACKAGES_MISSING"
    assert "r/setup.R" in exc.value.message
    assert exc.value.details["reason"] == "stale"


def test_lock_changed_since_restore_is_stale(configured, monkeypatch):
    """更新で renv.lock が変わったのに再復元していない＝古い版のまま黙って動かさない。"""
    (bridge.RENV_PROJECT / "renv.lock").write_text('{"R": {"Version": "4.6.1"}}\n', encoding="utf-8")
    run_dir = configured / "run"
    run_dir.mkdir()
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_R_PACKAGES_MISSING"
    assert exc.value.details["reason"] == "stale"
    assert exc.value.details["setup"] == "Rscript r/setup.R"


def test_browser_setting_wins(monkeypatch, tmp_path):
    exe = tmp_path / "chosen.exe"
    exe.write_text("", encoding="utf-8")
    other = tmp_path / "msedge.exe"
    other.write_text("", encoding="utf-8")
    monkeypatch.setenv("MSLIPIDMAPPER_BROWSER", str(exe))
    monkeypatch.setattr(bridge, "BROWSER_CANDIDATES", (str(other),))
    assert bridge.browser_path() == str(exe)


def test_browser_is_detected_from_standard_locations(monkeypatch, tmp_path):
    """Edge は通常 Program Files (x86) にあり、上流の自動検出はそこを見ない。こちらで探して渡す。"""
    edge = tmp_path / "x86" / "msedge.exe"
    edge.parent.mkdir()
    edge.write_text("", encoding="utf-8")
    monkeypatch.delenv("MSLIPIDMAPPER_BROWSER", raising=False)
    monkeypatch.setattr(bridge, "BROWSER_CANDIDATES",
                        (str(tmp_path / "missing" / "chrome.exe"), str(edge)))
    assert bridge.browser_path() == str(edge)


def test_no_browser_found_leaves_detection_to_upstream(monkeypatch, tmp_path):
    monkeypatch.delenv("MSLIPIDMAPPER_BROWSER", raising=False)
    monkeypatch.setattr(bridge, "BROWSER_CANDIDATES", (str(tmp_path / "nope.exe"),))
    assert bridge.browser_path() is None


def test_browser_candidates_cover_edge_x86():
    assert any("Program Files (x86)" in c and c.endswith("msedge.exe")
               for c in bridge.BROWSER_CANDIDATES)


def test_missing_renv_environment(configured, monkeypatch):
    # 実際の未復元状態: activate.R（追跡物）はあるが library が無い
    renv = configured / "fresh" / "r"
    (renv / "renv").mkdir(parents=True)
    (renv / "renv" / "activate.R").write_text("# fake", encoding="utf-8")
    monkeypatch.setattr(bridge, "RENV_PROJECT", renv)
    run_dir = configured / "run"
    run_dir.mkdir()
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_R_PACKAGES_MISSING"
    assert "r/setup.R" in exc.value.message
    assert exc.value.details["setup"] == "Rscript r/setup.R"
    assert exc.value.details["renv_project"] == str(renv)
    assert not (run_dir / "request.json").exists()


def test_restored_library_proceeds(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        (out / "bridge-result.json").write_text(json.dumps(
            {"status": "ok", "outputs": {}, "mapped_nodes": {}}), encoding="utf-8")
    """))
    res = bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert res["status"] == "ok"


def test_spawn_failure_is_run_failed(configured, monkeypatch):
    from metabolomix.core.atomic_io import DomainError

    def boom(*a, **k):
        raise DomainError("UNSUPPORTED_PLATFORM", "Windows 専用です")
    monkeypatch.setattr(bridge, "start_owned_process", boom)
    run_dir = configured / "run"
    run_dir.mkdir()
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_RUN_FAILED"
    assert exc.value.details["cause_code"] == "UNSUPPORTED_PLATFORM"
    assert exc.value.details["cause_message"] == "Windows 専用です"


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
