"""MSDIAL_EXE 未設定は人間にしか直せない。封筒はそれを明示し、手順を渡す。

MCP クライアントは環境変数を設定できないし、設定しても**起動中のサーバには
反映されない**。ここで LLM に「設定してください」とだけ返すと、LLM は設定を
試みて失敗するか、黙って諦める。誰の作業かを型で示す。
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


def test_console_plan_envelope_says_a_human_must_act(tmp_path, monkeypatch):
    monkeypatch.delenv("MSDIAL_EXE", raising=False)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    err = parsed["error"]
    assert err["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert err["details"]["human_action_required"] is True
    assert "restart" in _json.dumps(err["details"]).lower() or \
           "再起動" in _json.dumps(err["details"], ensure_ascii=False)


def test_console_plan_envelope_offers_a_copyable_command(tmp_path, monkeypatch):
    monkeypatch.delenv("MSDIAL_EXE", raising=False)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["details"]["how_to_set"]
