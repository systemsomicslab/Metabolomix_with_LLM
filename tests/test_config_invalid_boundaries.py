"""設定ファイル（lipidmix.local.toml）が読めないときの、ツールの境界のふるまい。

置き場所・サーバの振る舞いを設定ファイルへ移したので（計画 B-2a の決めたこと 2）、壊れた
設定ファイルがパーサ・記録・ハンドル・ライブラリの各ツールに届く。どのツールも生の例外で
落ちず、`CONFIG_INVALID` の封筒を返すか（記録・ハンドルの一覧・ライブラリ）、既定で続けて
問題を 1 行添える（意味論ダイジェスト）。設定ファイルの典型的な壊れ方として、移した旧キー
`[msdial] lbm` が残っている形を使う。
"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import handles, mcp_core, session_state, user_config
from tests.config_fixtures import use_config
from tests.curation_fixtures import write_alignment_set


@pytest.fixture
def broken_config(tmp_path, monkeypatch):
    """旧 `[msdial] lbm` が残った設定ファイル（どのキーを読んでも CONFIG_INVALID）。"""
    path = tmp_path / "conf" / "lipidmix.local.toml"
    path.parent.mkdir()
    path.write_text(f"[msdial]\nlbm = '{tmp_path / 'old.lbm2'}'\n", encoding="utf-8")
    monkeypatch.setenv(user_config.CONFIG_ENV, str(path))
    user_config._cache.clear()
    with pytest.raises(user_config.ConfigInvalidError):
        user_config.get_setting("server.caveat_mode")
    return path


@pytest.fixture
def fresh_session(monkeypatch, tmp_path):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path / "data")
    return session_state.session


def _error(out: str) -> dict:
    return json.loads(out)["error"]


# ---------- 意味論ダイジェスト: 落ちずに digest で続け、1 行添える ----------

def test_arf2_parser_survives_a_broken_config(tmp_path, broken_config, fresh_session):
    from metabolomix.arf2.tools import arf2_parser

    paths = write_alignment_set(tmp_path / "neg")
    out = arf2_parser(file_path=str(paths["arf2"]))
    assert out.startswith(session_state.SEMANTICS_CAVEAT)          # 既定の digest で続ける
    assert "CONFIG_INVALID" in out and "[server] caveat_mode" in out
    assert "lipidmix.local.toml" in out
    assert "[library] lbm" in out                                  # 何が問題かを言う
    assert str(tmp_path / "conf") not in out                       # 設定ファイルの置き場所は出さない
    assert out.count("CONFIG_INVALID") == 1


def test_load_dataset_survives_a_broken_config(tmp_path, broken_config, fresh_session):
    import server

    paths = write_alignment_set(tmp_path / "neg")
    out = server.load_dataset(directory=str(paths["arf2"].parent))
    assert "データセット読み込み" in out
    assert session_state.SEMANTICS_CAVEAT in out
    assert "CONFIG_INVALID" in out
    assert mcp_core.DATA_DIR == paths["arf2"].parent
    assert session_state.session.arf.features                      # arf の段まで進んでいる


def test_caveat_off_still_works_and_the_note_is_absent(tmp_path, monkeypatch, fresh_session):
    use_config(monkeypatch, tmp_path, {"server.caveat_mode": "off"})
    assert fresh_session.maybe_prepend_caveat("BODY", topic="arf") == "BODY"


# ---------- 記録（レポート・objective）: CONFIG_INVALID の封筒 ----------

def test_report_tools_return_the_envelope(broken_config, fresh_session):
    from metabolomix.tools import reports

    for out in (reports.write_report("a-1", "DS", "本文"),
                reports.read_report("a-1"),
                reports.list_reports()):
        err = _error(out)
        assert err["code"] == "CONFIG_INVALID"
        assert "[library] lbm" in err["message"]
        assert err["details"]["config_file"] == str(broken_config)
    assert not (mcp_core.DATA_DIR / "reports").exists()            # 黙って解析フォルダへ書かない


def test_objective_tools_return_the_envelope_without_half_applied_state(broken_config, fresh_session):
    from metabolomix.tools import objective

    out = objective.record_objective(
        analysis_id="a-1", dataset="DS", polarity="NEG", groups=["a", "b"], comparison="a vs b",
        sub_questions=["Q"], assay_kind="lipid")
    assert _error(out)["code"] == "CONFIG_INVALID"
    assert fresh_session.assay_kind == "unknown"                   # 種別だけが切り替わらない
    for out in (objective.update_objective("a-1", status="done"),
                objective.log_search("a-1", "Q1", "q", 0),
                objective.knowledge_coverage("a-1")):
        assert _error(out)["code"] == "CONFIG_INVALID"


def test_record_objective_rejects_an_unknown_assay_kind_before_touching_anything(fresh_session):
    from metabolomix.tools import objective

    out = objective.record_objective(
        analysis_id="a-1", dataset="DS", polarity="NEG", groups=["a", "b"], comparison="a vs b",
        sub_questions=["Q"], assay_kind="protein")
    assert "assay_kind は lipid/metabolite/unknown のいずれかです" in out
    assert fresh_session.assay_kind == "unknown"
    assert not (mcp_core.DATA_DIR / "reports").exists()


def test_pathway_plot_returns_the_envelope(broken_config, fresh_session):
    from metabolomix.pathway import tools as pathway_tools

    fresh_session.pathway.last_run_id = "pw_20260101_000000_abcdef"
    assert _error(pathway_tools.pathway_plot())["code"] == "CONFIG_INVALID"


# ---------- ハンドル: 一覧は封筒、上限の読めない要約は None ----------

def test_dataset_handle_list_returns_the_envelope(broken_config):
    from metabolomix.tools import handle_tools

    err = _error(handle_tools.dataset_handle_list())
    assert err["code"] == "CONFIG_INVALID"
    assert err["details"]["config_file"] == str(broken_config)


def test_dataset_handle_open_refuses_with_the_envelope(tmp_path, broken_config):
    from metabolomix.tools import handle_tools

    paths = write_alignment_set(tmp_path / "neg")
    err = _error(handle_tools.dataset_handle_open(directory=str(paths["arf2"].parent)))
    assert err["code"] == "CONFIG_INVALID"
    assert handles.registry.list() == []


def test_max_handles_or_none_hides_a_broken_config(broken_config):
    from metabolomix.tools import handle_tools

    assert handle_tools._max_handles_or_none() is None


# ---------- ライブラリ: CONFIG_INVALID を言う ----------

_MSP = """NAME: GABA
PRECURSORMZ: 104.0706
PRECURSORTYPE: [M+H]+
IONMODE: Positive
Num Peaks: 1
87.04 100
"""


def test_library_load_reports_config_invalid(tmp_path, broken_config, fresh_session):
    from metabolomix.library import tools as library_tools

    msp = tmp_path / "lib.msp"
    msp.write_text(_MSP, encoding="utf-8")
    out = json.loads(library_tools.library_load(file_path=str(msp)))
    assert out["status"] == "error"
    assert out["code"] == "CONFIG_INVALID"
    assert out["details"]["config_file"] == str(broken_config)


def test_remembered_library_note_says_config_invalid(tmp_path, monkeypatch):
    """記憶したライブラリを開く段で設定ファイルが読めなくなっていたら、ハンドルは開き理由に CONFIG_INVALID。"""
    from metabolomix.library import store
    from metabolomix.tools import handle_tools
    from tests.dbs_fixture import record, write_dbs
    from tests.msdial_folder_fixtures import write_msdial_folder

    folder = tmp_path / "gui"
    write_msdial_folder(folder, spots=[{"name": "GABA"}], samples=[
        {"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(3)])
    write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record()], ["k"])])
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})

    def broken_cache_dir():
        raise user_config.ConfigInvalidError("lipidmix.local.toml の TOML 構文エラー: x",
                                             config_file=str(tmp_path / "lipidmix.local.toml"))

    monkeypatch.setattr(store, "cache_dir", broken_cache_dir)
    out = json.loads(handle_tools.dataset_handle_open(directory=str(folder), library="auto"))
    assert out["dataset_handle"].startswith("dh_")
    assert "CONFIG_INVALID" in out["library_note"]
    assert "ConfigInvalidError" not in out["library_note"]
