# tests/test_version_update.py
"""origin/main からの遅れ検知と、その通知の載せ方を固定する。

自動 pull はしない。走っている MCP サーバの裏でコードが変わると、import 済み
モジュールだけが古いまま残る混在状態になるため（docs/HISTRY.md 2026-09-04(4)）。
ここで縛るのは「黙るべきときに黙るか」が中心——通知機能が解析の邪魔をしては
いけないので、判定できない場合はすべて通知なしに倒す。
"""
import json as _json

import pytest

from metabolomix.core import version


@pytest.fixture(autouse=True)
def reset_version_state():
    """モジュール大域のキャッシュ・fetch 完了フラグをテスト間で持ち越さない。"""
    def _clear():
        version._update_cache = version._UNSET
        version._fetch_ok.clear()
    _clear()
    yield
    _clear()


def _fake_git(responses):
    """`_git` の差し替え。args の先頭サブコマンドで引く。"""
    def _run(*args, timeout=None):
        return responses.get(args[0])
    return _run


# ---------- 判定 ----------

def test_silent_until_fetch_succeeds(monkeypatch):
    """fetch が終わっていない間は通知しない（古い ref で誤った断定をしない）。"""
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "main", "rev-list": "7"}))
    assert version.update_status() is None


def test_reports_behind_count_on_main(monkeypatch):
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "main", "rev-list": "7"}))
    version._fetch_ok.set()
    status = version.update_status()
    assert status["behind"] == 7
    assert "git pull" in status["message"]
    assert "再起動" in status["message"]


def test_silent_when_up_to_date(monkeypatch):
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "main", "rev-list": "0"}))
    version._fetch_ok.set()
    assert version.update_status() is None


def test_silent_on_feature_branch(monkeypatch):
    """開発者は main 以外で作業する。そこで毎回鳴らせない。"""
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "feature/x", "rev-list": "12"}))
    version._fetch_ok.set()
    assert version.update_status() is None


def test_silent_on_detached_head(monkeypatch):
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "HEAD", "rev-list": "3"}))
    version._fetch_ok.set()
    assert version.update_status() is None


def test_silent_when_git_unavailable(monkeypatch):
    """git が無い・リポジトリでない場合は `_git` が None を返す。"""
    monkeypatch.setattr(version, "_git", _fake_git({}))
    version._fetch_ok.set()
    assert version.update_status() is None


def test_silent_when_count_is_not_a_number(monkeypatch):
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "main", "rev-list": "??"}))
    version._fetch_ok.set()
    assert version.update_status() is None


def test_status_is_computed_once(monkeypatch):
    """プロセス生存中は固定（起動時のコードが何かを示す値なので）。"""
    calls = []

    def _run(*args, timeout=None):
        calls.append(args[0])
        return {"rev-parse": "main", "rev-list": "2"}[args[0]]

    monkeypatch.setattr(version, "_git", _run)
    version._fetch_ok.set()
    assert version.update_status()["behind"] == 2
    calls.clear()
    assert version.update_status()["behind"] == 2
    assert calls == []


# ---------- fetch の起動 ----------

def test_fetch_failure_leaves_checker_silent(monkeypatch):
    """オフライン・認証失敗は例外にせず、通知なしとして扱う。"""
    monkeypatch.setattr(version, "_git", _fake_git({
        "rev-parse": "main", "rev-list": "5"}))  # fetch だけ None
    version._run_update_fetch()
    assert not version._fetch_ok.is_set()
    assert version.update_status() is None


def test_fetch_success_enables_reporting(monkeypatch):
    monkeypatch.setattr(version, "_git", _fake_git({
        "fetch": "", "rev-parse": "main", "rev-list": "5"}))
    version._run_update_fetch()
    assert version._fetch_ok.is_set()
    assert version.update_status()["behind"] == 5


def test_start_update_check_runs_once(monkeypatch):
    started = []
    monkeypatch.setattr(version, "_check_started", False)
    monkeypatch.setattr(version.threading, "Thread",
                        lambda **kw: type("T", (), {"start": lambda self: started.append(kw)})())
    version.start_update_check()
    version.start_update_check()
    assert len(started) == 1
    assert started[0]["daemon"] is True


# ---------- 通知の載せ方 ----------

def test_dataset_status_carries_notice(monkeypatch):
    import numpy as np
    from metabolomix.core import session_state
    from metabolomix.mztab.dataset_state import DatasetState
    from metabolomix.tools.mztab_tools import dataset_status

    ds = DatasetState()
    ds.feature_matrix = np.ones((2, 2))
    ds.sample_names = ["s1", "s2"]
    ds.feature_ids = ["1", "2"]
    session_state.session = session_state.AnalysisSession()
    session_state.session.dataset = ds

    monkeypatch.setattr(version, "update_status", lambda: None)
    assert "update_available" not in _json.loads(dataset_status())

    monkeypatch.setattr(version, "update_status",
                        lambda: {"behind": 4, "message": "更新があります"})
    parsed = _json.loads(dataset_status())
    assert parsed["update_available"]["behind"] == 4


def test_console_status_carries_notice(tmp_path, monkeypatch):
    from metabolomix.console.job_manager import create_job
    from metabolomix.tools.console_tools import console_status

    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="negative", measure="peak_height")

    monkeypatch.setattr(version, "update_status", lambda: None)
    assert "update_available" not in _json.loads(console_status(str(job_path)))

    monkeypatch.setattr(version, "update_status",
                        lambda: {"behind": 1, "message": "更新があります"})
    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["update_available"]["behind"] == 1


def test_load_dataset_shows_notice(tmp_path, monkeypatch):
    """入口ツールは Markdown を返すので、キーではなく 1 行のブロックで出す。"""
    from metabolomix.core import mcp_core, session_state
    from metabolomix.tools.dataset import load_dataset

    session_state.session = session_state.AnalysisSession()
    original = mcp_core.DATA_DIR
    try:
        monkeypatch.setattr(version, "update_status", lambda: None)
        assert "更新があります" not in load_dataset(str(tmp_path))

        session_state.session = session_state.AnalysisSession()
        monkeypatch.setattr(version, "update_status",
                            lambda: {"behind": 2, "message": "更新があります"})
        assert "更新があります" in load_dataset(str(tmp_path))
    finally:
        mcp_core.DATA_DIR = original
