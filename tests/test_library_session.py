"""ライブラリのセッションスロットとパス解決。"""
from metabolomix.core import mcp_core, path_resolvers, session_state


def test_the_new_slot_does_not_disturb_the_others():
    s = session_state.AnalysisSession()
    assert s.library.store is None
    assert s.library.source_path is None
    assert s.library.last_match is None
    # パーサ別のスロットは無い（B-2b で消した）。残りのスロットが消えていないこと。
    for name in ("curation", "pathway"):
        assert getattr(s, name) is not None
    assert s.dataset is None


def test_a_dbs_is_preferred_over_an_msp(tmp_path, monkeypatch):
    (tmp_path / "other.msp").write_text("NAME: x\n", encoding="utf-8")
    (tmp_path / "P_Loaded.msp2.dbs").write_bytes(b"PK\x03\x04")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)   # 実物は pathlib.Path
    assert path_resolvers.resolve_library_path().endswith("P_Loaded.msp2.dbs")


def test_an_msp_is_used_when_no_dbs_exists(tmp_path, monkeypatch):
    (tmp_path / "lib.msp").write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)   # 実物は pathlib.Path
    assert path_resolvers.resolve_library_path().endswith("lib.msp")


def test_an_explicit_path_wins(tmp_path, monkeypatch):
    (tmp_path / "P_Loaded.msp2.dbs").write_bytes(b"PK\x03\x04")
    explicit = tmp_path / "explicit.msp"
    explicit.write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)   # 実物は pathlib.Path
    assert path_resolvers.resolve_library_path(str(explicit)) == str(explicit)


def test_nothing_found_returns_none(tmp_path, monkeypatch):
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)   # 実物は pathlib.Path
    assert path_resolvers.resolve_library_path() is None


def test_bare_msp2_and_lbm2_are_never_candidates(tmp_path, monkeypatch):
    # .msp2（ASCII .msp 指定時のみ実体が入り、脂質経路では0バイトのことがある）と
    # .lbm2（脂質専用 in-silico ライブラリ、本機能の入口としては出さない）は、
    # 他に何も無くても候補になってはいけない。
    (tmp_path / "Loaded.msp2").write_bytes(b"")
    (tmp_path / "Loaded.lbm2").write_bytes(b"PK\x03\x04")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)   # 実物は pathlib.Path
    assert path_resolvers.resolve_library_path() is None


# --------------------------------------------------------------------------
# 研究室ライブラリは pos / neg の 2 ファイルで、lipidmix.local.toml の
# [library] msp_positive / msp_negative で置き場所を指す。一度に読むのは 1 つだけ。複数の候補が
# あるときは更新日時で黙って 1 つ選ばず、どれにするかを呼び出し側に問い返す
# （極性違いのライブラリで照合しても、候補が少し減るだけで誤りに気づけない）。
# --------------------------------------------------------------------------
import pytest  # noqa: E402

from metabolomix.core.path_resolvers import LibraryPathError  # noqa: E402
from tests.config_fixtures import use_config  # noqa: E402


def _libs(tmp_path):
    pos = tmp_path / "libs" / "lab_pos.msp"
    neg = tmp_path / "libs" / "lab_neg.msp"
    pos.parent.mkdir()
    pos.write_text("NAME: p\n", encoding="utf-8")
    neg.write_text("NAME: n\n", encoding="utf-8")
    return pos, neg


def test_the_setting_for_the_requested_polarity_is_used(tmp_path, monkeypatch):
    pos, neg = _libs(tmp_path)
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path / "data")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(pos),
                                       "library.msp_negative": str(neg)})
    assert path_resolvers.resolve_library_path(ion_mode="negative") == str(neg)
    assert path_resolvers.resolve_library_path(ion_mode="Positive") == str(pos)


def test_the_requested_polarity_beats_a_dbs_in_the_data_dir(tmp_path, monkeypatch):
    pos, _ = _libs(tmp_path)
    (tmp_path / "P_Loaded.msp2.dbs").write_bytes(b"PK\x03\x04")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(pos)})
    assert path_resolvers.resolve_library_path(ion_mode="positive") == str(pos)


def test_an_explicit_path_beats_the_settings(tmp_path, monkeypatch):
    pos, neg = _libs(tmp_path)
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(pos)})
    assert path_resolvers.resolve_library_path(str(neg), ion_mode="positive") == str(neg)


def test_a_setting_pointing_nowhere_is_an_error_not_a_silent_fallback(tmp_path, monkeypatch):
    (tmp_path / "other.msp").write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    use_config(monkeypatch, tmp_path,
               {"library.msp_negative": str(tmp_path / "gone" / "secret_neg.msp")})
    with pytest.raises(LibraryPathError) as exc:
        path_resolvers.resolve_library_path(ion_mode="negative")
    assert exc.value.code == "MSP_ENV_NOT_FOUND"
    assert "[library] msp_negative" in exc.value.message
    assert "gone" not in exc.value.message   # 置き場所（ディレクトリ）は戻り値に出さない


def test_an_unset_polarity_setting_falls_back_to_the_data_dir(tmp_path, monkeypatch):
    (tmp_path / "lib.msp").write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    assert path_resolvers.resolve_library_path(ion_mode="negative").endswith("lib.msp")


def test_both_settings_without_a_polarity_is_ambiguous(tmp_path, monkeypatch):
    pos, neg = _libs(tmp_path)
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path / "data")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(pos),
                                       "library.msp_negative": str(neg)})
    with pytest.raises(LibraryPathError) as exc:
        path_resolvers.resolve_library_path()
    assert exc.value.code == "MSP_AMBIGUOUS"
    assert "ion_mode" in exc.value.message


def test_a_single_setting_is_used_without_a_polarity(tmp_path, monkeypatch):
    _, neg = _libs(tmp_path)
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path / "data")
    use_config(monkeypatch, tmp_path, {"library.msp_negative": str(neg)})
    assert path_resolvers.resolve_library_path() == str(neg)


def test_a_dbs_in_the_data_dir_still_wins_without_a_polarity(tmp_path, monkeypatch):
    """`*_Loaded.msp2.dbs` はその run が実際に使った絞り込み済みの参照（spec 2026-09-19）。"""
    pos, neg = _libs(tmp_path)
    (tmp_path / "P_Loaded.msp2.dbs").write_bytes(b"PK\x03\x04")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(pos),
                                       "library.msp_negative": str(neg)})
    assert path_resolvers.resolve_library_path().endswith("P_Loaded.msp2.dbs")


def test_several_msp_in_the_data_dir_are_ambiguous_rather_than_picked_by_mtime(tmp_path, monkeypatch):
    (tmp_path / "a_pos.msp").write_text("NAME: x\n", encoding="utf-8")
    (tmp_path / "b_neg.msp").write_text("NAME: y\n", encoding="utf-8")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    with pytest.raises(LibraryPathError) as exc:
        path_resolvers.resolve_library_path()
    assert exc.value.code == "MSP_AMBIGUOUS"
    assert "a_pos.msp" in exc.value.message and "b_neg.msp" in exc.value.message


def test_an_explicit_path_that_does_not_exist_is_an_error(tmp_path, monkeypatch):
    """以前は data ディレクトリの探索へ黙って落ちていた（別のライブラリを掴む）。"""
    (tmp_path / "lib.msp").write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    with pytest.raises(LibraryPathError) as exc:
        path_resolvers.resolve_library_path(str(tmp_path / "typo.msp"))
    assert exc.value.code == "LIBRARY_NOT_FOUND"


def test_an_unknown_polarity_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    with pytest.raises(LibraryPathError) as exc:
        path_resolvers.resolve_library_path(ion_mode="both")
    assert exc.value.code == "INVALID_ION_MODE"
