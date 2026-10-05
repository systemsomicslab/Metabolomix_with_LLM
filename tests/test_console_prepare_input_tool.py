"""console_prepare_input（MCP ツール層）と、混在エラーが自分の直し方を指すこと。

「正しいエラーメッセージ」と「実行可能な指示」は別物。MIXED_RAW_FORMATS の文面は
正しかったが、MCP クライアントにはフォルダを作る手段が無かった。封筒が
required_tools で自分の直し方を指し、その先が実在するツールであることを縛る。
"""
from __future__ import annotations

import json as _json
from pathlib import Path


def _sciex(root: Path, stem: str):
    (root / f"{stem}.wiff").write_text("p", encoding="ascii")
    (root / f"{stem}.wiff.scan").write_text("s", encoding="ascii")
    (root / f"{stem}.wiff2").write_text("n", encoding="ascii")


def test_prepare_input_creates_single_format_folder(tmp_path):
    from lipidmix.tools.console_tools import console_prepare_input
    src = tmp_path / "raw"
    src.mkdir()
    _sciex(src, "a")
    _sciex(src, "b")
    parsed = _json.loads(console_prepare_input(dataset_root=str(src),
                                               keep_extension="wiff",
                                               out_dir=str(tmp_path / "out")))
    assert parsed["status"] == "prepared"
    assert parsed["primary"] == 2
    assert parsed["formats"] == {"wiff": 2}


def test_prepare_input_defaults_out_dir_beside_source(tmp_path):
    """置き場所を毎回考えさせない。既定は <元フォルダ名>_<拡張子> の兄弟。"""
    from lipidmix.tools.console_tools import console_prepare_input
    src = tmp_path / "POS"
    src.mkdir()
    _sciex(src, "a")
    parsed = _json.loads(console_prepare_input(dataset_root=str(src),
                                               keep_extension="wiff"))
    assert Path(parsed["out_dir"]).parent == tmp_path
    assert Path(parsed["out_dir"]).name == "POS_wiff"


def test_prepare_input_returns_error_envelope_not_exception(tmp_path):
    from lipidmix.tools.console_tools import console_prepare_input
    src = tmp_path / "raw"
    src.mkdir()
    _sciex(src, "a")
    parsed = _json.loads(console_prepare_input(dataset_root=str(src),
                                               keep_extension="raw"))
    assert parsed["error"]["code"] == "INPUT_PREP_FAILED"


def _console_env(tmp_path: Path, monkeypatch) -> Path:
    """Console 実行体と LBM を偽装し、Negative のメソッドファイルを返す。"""
    from lipidmix.core import session_state
    session_state.session = session_state.AnalysisSession()
    app = tmp_path / "app"
    app.mkdir()
    (app / "x.lbm2").touch()
    monkeypatch.setenv("MSDIAL_EXE", str(app / "MSDIALCUI.exe"))
    monkeypatch.delenv("MSDIAL_LBM", raising=False)
    monkeypatch.setattr("lipidmix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    return method


def test_mixed_raw_formats_offers_choices_and_retries_plan(tmp_path, monkeypatch):
    """混在は止めて終わりにせず、選べる形式を示して console_plan 自身へ戻す。"""
    method = _console_env(tmp_path, monkeypatch)
    src = tmp_path / "raw"
    src.mkdir()
    _sciex(src, "a")
    _sciex(src, "b")
    from lipidmix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(src), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "MIXED_RAW_FORMATS"
    assert parsed["error"]["required_tools"] == ["console_plan"]
    assert parsed["error"]["details"]["choices"] == [{"keep_extension": "wiff", "count": 2},
                                          {"keep_extension": "wiff2", "count": 2}]
    # 選ぶ前に何も作らない（呼び直しで別の形式を選べるように）。
    assert not (tmp_path / "raw_wiff").exists()
    assert not (tmp_path / "raw_wiff2").exists()


def test_keep_extension_on_mixed_folder_plans_on_prepared_sibling(tmp_path, monkeypatch):
    method = _console_env(tmp_path, monkeypatch)
    src = tmp_path / "raw"
    src.mkdir()
    _sciex(src, "a")
    _sciex(src, "b")
    before = sorted(p.name for p in src.iterdir())
    from lipidmix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(src), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      keep_extension="wiff2"))
    assert parsed["status"] == "planned"
    assert parsed["input_count"] == 2
    assert parsed["keep_extension"] == "wiff2"
    prepared = tmp_path / "raw_wiff2"
    assert parsed["prepared_input"]["out_dir"] == str(prepared)
    assert parsed["prepared_input"]["source"] == str(src)
    job = _json.loads(Path(parsed["job_path"]).read_text(encoding="utf-8"))
    assert job["source"]["dataset_root"] == str(prepared)
    assert sorted(p.suffix for p in prepared.iterdir() if p.is_file()) == [".wiff2", ".wiff2"]
    assert sorted(p.name for p in src.iterdir()) == before


def test_keep_extension_absent_from_folder_is_rejected(tmp_path, monkeypatch):
    method = _console_env(tmp_path, monkeypatch)
    src = tmp_path / "raw"
    src.mkdir()
    _sciex(src, "a")
    from lipidmix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(src), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      keep_extension="mzml"))
    assert parsed["error"]["code"] == "MIXED_RAW_FORMATS"
    assert parsed["error"]["details"]["formats"] == {"wiff": 1, "wiff2": 1}
    assert not (tmp_path / "raw_mzml").exists()


def test_keep_extension_on_single_format_folder_plans_in_place(tmp_path, monkeypatch):
    method = _console_env(tmp_path, monkeypatch)
    src = tmp_path / "raw"
    src.mkdir()
    (src / "a.mzML").write_text("m", encoding="ascii")
    from lipidmix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(src), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      keep_extension=".mzML"))
    assert parsed["status"] == "planned"
    assert "prepared_input" not in parsed
    job = _json.loads(Path(parsed["job_path"]).read_text(encoding="utf-8"))
    assert job["source"]["dataset_root"] == str(src)
    assert not (tmp_path / "raw_mzml").exists()


def test_prepared_folder_then_plans_cleanly(tmp_path, monkeypatch):
    """作った先がそのまま console_plan を通ること（この2つが噛み合わないと意味がない）。"""
    method = _console_env(tmp_path, monkeypatch)
    src = tmp_path / "raw"
    src.mkdir()
    _sciex(src, "a")

    from lipidmix.tools.console_tools import console_plan, console_prepare_input
    prepared = _json.loads(console_prepare_input(dataset_root=str(src),
                                                 keep_extension="wiff"))
    parsed = _json.loads(console_plan(dataset_root=prepared["out_dir"],
                                      method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["status"] == "planned"
    assert parsed["input_count"] == 1
