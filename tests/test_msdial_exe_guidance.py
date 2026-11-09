"""Console の実行体（`[msdial] exe`）の未設定は人間にしか直せない。封筒はそれを明示し、手順を渡す。

MCP クライアントは設定ファイル（lipidmix.local.toml）を編集できない。ここで LLM に
「設定してください」とだけ返すと、LLM は設定を試みて失敗するか、黙って諦める。
誰の作業かを型で示す。設定ファイルは呼ばれるたびに読むので、再起動は要らない。
"""
from __future__ import annotations

import json as _json

from metabolomix.console.runner import msdial_exe_candidates


def test_finds_console_executables_under_search_roots(tmp_path):
    app = tmp_path / "MSDIAL.v5.5"
    app.mkdir()
    (app / "MSDIALCUI.exe").touch()
    (app / "MSDIAL.exe").touch()  # GUI は候補にしない
    found = msdial_exe_candidates([tmp_path])
    assert [p.rsplit("\\", 1)[-1].rsplit("/", 1)[-1] for p in found] == ["MSDIALCUI.exe"]


def test_returns_empty_for_missing_roots(tmp_path):
    assert msdial_exe_candidates([tmp_path / "nope"]) == []


def test_does_not_descend_indefinitely(tmp_path):
    """ホーム配下を無制限に歩くと計画が数十秒止まる。"""
    deep = tmp_path / "a" / "b" / "c" / "d" / "e"
    deep.mkdir(parents=True)
    (deep / "MSDIALCUI.exe").touch()
    assert msdial_exe_candidates([tmp_path], max_depth=2) == []


def _unresolved_exe_envelope() -> dict:
    """実行体が未設定のとき、pipeline の受付（pipeline_plan / pipeline_run）が返すエラー。"""
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import service
    try:
        service._resolve_exe_path()
    except DomainError as exc:
        return {"error": {"code": exc.code, "message": exc.message, "details": exc.details}}
    raise AssertionError("実行体が解けてしまった（テスト環境に設定がある）")


def test_exe_envelope_says_a_human_must_act():
    err = _unresolved_exe_envelope()["error"]
    assert err["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert err["details"]["human_action_required"] is True
    assert "restart" in _json.dumps(err["details"]).lower() or            "再起動" in _json.dumps(err["details"], ensure_ascii=False)


def test_exe_envelope_offers_a_copyable_command():
    assert _unresolved_exe_envelope()["error"]["details"]["how_to_set"]
