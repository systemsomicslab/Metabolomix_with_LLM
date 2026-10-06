# tests/test_server_update.py
"""配布先のクローンを更新する `server_update` の振る舞いを固定する。

このツールは**利用者の作業ツリーを書き換える**唯一のツールなので、縛るのは
「やらないこと」が中心になる。壊し方は 2 つあって、どちらも黙って起きる:

1. 手元の未コミット変更を pull で巻き込む（メンバーが data/ の外に置いた
   設定やメモが消える）。
2. 依存の更新に失敗したまま「更新しました」と言う（次回起動で ImportError）。

プロセス操作はしない。再起動は利用者が手で行う（DEPLOY.md）。走っている
サーバの裏でコードを差し替えても import 済みモジュールは古いままなので、
「再起動してください」を戻り値から落とさないことも回帰対象にする。
"""
import json as _json

import pytest

from metabolomix.core import version


@pytest.fixture(autouse=True)
def reset_version_state():
    def _clear():
        version._update_cache = version._UNSET
        version._fetch_ok.clear()
    _clear()
    yield
    _clear()


class FakeGit:
    """`_git` の差し替え。呼ばれた引数列を記録し、前方一致で応答を引く。

    `test_version_update.py` の `_fake_git` はサブコマンド名だけで引くが、
    更新経路は `rev-parse --abbrev-ref HEAD` と `rev-parse --short HEAD` の
    ように同じサブコマンドを別の意図で 2 回呼ぶため、引数列で引く。
    """

    def __init__(self, responses):
        self.responses = responses
        self.calls = []

    def __call__(self, *args, timeout=None):
        self.calls.append(args)
        for prefix, value in self.responses.items():
            if args[:len(prefix)] == prefix:
                return value
        return None

    def ran(self, subcommand):
        return any(call[0] == subcommand for call in self.calls)


#: 遅れが 1 コミットあり、requirements.txt は変わっていない正常系。
CLEAN_BEHIND_ONE = {
    ("rev-parse", "--abbrev-ref", "HEAD"): "main",
    ("status", "--porcelain"): "",
    ("fetch",): "",
    ("rev-list", "--count"): "1",
    ("diff", "--name-only"): "metabolomix/arf/reader.py\nUSAGE.md",
    ("merge", "--ff-only"): "Updating 4253408..abc1234",
    ("rev-parse", "--short", "HEAD"): "abc1234",
}


def _with(**overrides):
    """CLEAN_BEHIND_ONE を 1 項目だけ差し替えた応答表を作る。"""
    responses = dict(CLEAN_BEHIND_ONE)
    for key, value in overrides.items():
        responses[tuple(key.split("__"))] = value
    return responses


@pytest.fixture
def no_pip(monkeypatch):
    """pip を呼ばせない。呼ばれたら記録だけして成功を返す。"""
    calls = []
    monkeypatch.setattr(version, "_pip_install",
                        lambda: (calls.append(True), (True, ""))[1])
    return calls


# ---------- 拒否（手元を壊さないための門） ----------

def test_refuses_when_worktree_is_dirty(monkeypatch, no_pip):
    """未コミット変更があるときは pull しない。消えたら復元できないため。"""
    git = FakeGit(_with(**{"status__--porcelain": " M metabolomix/arf/reader.py"}))
    monkeypatch.setattr(version, "_git", git)

    result = version.apply_update()

    assert result["status"] == "refused"
    assert result["reason"] == "dirty_worktree"
    assert not git.ran("merge")


def test_refuses_when_not_on_main(monkeypatch, no_pip):
    """開発者の feature ブランチを main へ早送りしない。"""
    git = FakeGit(_with(**{"rev-parse__--abbrev-ref__HEAD": "feat/x"}))
    monkeypatch.setattr(version, "_git", git)

    result = version.apply_update()

    assert result["status"] == "refused"
    assert result["reason"] == "not_on_main"
    assert not git.ran("merge")


def test_refuses_when_fast_forward_is_impossible(monkeypatch, no_pip):
    """ローカルに独自コミットがあると ff-only が失敗する。merge も rebase もしない。"""
    git = FakeGit(_with(**{"merge__--ff-only": None}))
    monkeypatch.setattr(version, "_git", git)

    result = version.apply_update()

    assert result["status"] == "refused"
    assert result["reason"] == "not_fast_forward"
    assert no_pip == []


def test_refuses_when_fetch_fails(monkeypatch, no_pip):
    """オフライン・認証失敗。古い ref のまま「最新です」と言わない。"""
    git = FakeGit(_with(fetch=None))
    monkeypatch.setattr(version, "_git", git)

    result = version.apply_update()

    assert result["status"] == "refused"
    assert result["reason"] == "fetch_failed"
    assert not git.ran("merge")


def test_refuses_when_git_unavailable(monkeypatch, no_pip):
    """git が無い・リポジトリでない（zip 配布など）。"""
    monkeypatch.setattr(version, "_git", FakeGit({}))

    result = version.apply_update()

    assert result["status"] == "refused"
    assert result["reason"] == "git_unavailable"


# ---------- 何もしない正常系 ----------

def test_reports_up_to_date_without_pulling(monkeypatch, no_pip):
    git = FakeGit(_with(**{"rev-list__--count": "0"}))
    monkeypatch.setattr(version, "_git", git)

    result = version.apply_update()

    assert result["status"] == "up_to_date"
    assert not git.ran("merge")
    assert no_pip == []


def test_up_to_date_wins_over_dirty_worktree(monkeypatch, no_pip):
    """遅れていないなら手元が汚れていても拒否ではない（何もしないので安全）。"""
    git = FakeGit(_with(**{"rev-list__--count": "0",
                           "status__--porcelain": " M check.py"}))
    monkeypatch.setattr(version, "_git", git)

    assert version.apply_update()["status"] == "up_to_date"


# ---------- 更新の適用 ----------

def test_pulls_and_asks_for_restart(monkeypatch, no_pip):
    """pull だけでは古いプロセスが動き続ける。再起動の指示を落とさない。"""
    git = FakeGit(CLEAN_BEHIND_ONE)
    monkeypatch.setattr(version, "_git", git)

    result = version.apply_update()

    assert result["status"] == "updated"
    assert result["behind"] == 1
    assert result["head"] == "abc1234"
    assert "再起動" in result["message"]


# ---------- 依存の更新 ----------

def test_installs_dependencies_when_requirements_changed(monkeypatch):
    """git pull は requirements.txt の中身を反映するだけで、入れてはくれない。"""
    monkeypatch.setattr(version, "_git", FakeGit(_with(
        **{"diff__--name-only": "requirements.txt\nserver.py"})))
    calls = []
    monkeypatch.setattr(version, "_pip_install",
                        lambda: (calls.append(True), (True, ""))[1])

    result = version.apply_update()

    assert result["dependencies"] == "installed"
    assert len(calls) == 1


def test_skips_dependency_install_when_requirements_unchanged(monkeypatch, no_pip):
    monkeypatch.setattr(version, "_git", FakeGit(CLEAN_BEHIND_ONE))

    result = version.apply_update()

    assert result["dependencies"] == "unchanged"
    assert no_pip == []


def test_reports_dependency_failure_without_claiming_success(monkeypatch):
    """pip が落ちたら、コードだけ新しく依存が古い状態になる。黙って成功にしない。"""
    monkeypatch.setattr(version, "_git", FakeGit(_with(
        **{"diff__--name-only": "requirements.txt"})))
    monkeypatch.setattr(version, "_pip_install",
                        lambda: (False, "ERROR: No matching distribution"))

    result = version.apply_update()

    assert result["status"] == "updated_with_warning"
    assert result["dependencies"] == "failed"
    assert "No matching distribution" in result["message"]


# ---------- 通知との接続 ----------

def test_notice_names_the_tool(monkeypatch):
    """`/` スラッシュコマンドは local stdio では使えない（claude-code#82045）。
    利用者がツールに辿り着く唯一の導線は通知本文なので、ツール名を書く。"""
    monkeypatch.setattr(version, "_git", FakeGit(_with(
        **{"rev-list__--count": "3"})))
    version._fetch_ok.set()

    assert "server_update" in version.update_status()["message"]


# ---------- MCP 公開面 ----------

def test_tool_returns_compact_json(monkeypatch):
    from metabolomix.tools import maintenance

    monkeypatch.setattr(version, "apply_update",
                        lambda: {"status": "up_to_date", "head": "abc1234"})

    payload = maintenance.server_update()

    assert _json.loads(payload)["status"] == "up_to_date"
    assert "\n" not in payload  # json_payload（indent なし）で返している


def test_tool_reports_running_version(monkeypatch, no_pip):
    """どの版が動いていたかが分からないと、更新前後の比較ができない。"""
    from metabolomix.tools import maintenance

    monkeypatch.setattr(version, "_git", FakeGit(CLEAN_BEHIND_ONE))

    parsed = _json.loads(maintenance.server_update())

    assert parsed["running_version"] == version.server_version()
