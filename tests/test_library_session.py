"""ライブラリのセッションスロットとパスの確認（探索は無い。B-2b）。"""
import pytest

from metabolomix.core import path_resolvers, session_state
from metabolomix.core.path_resolvers import LibraryPathError


def test_the_new_slot_does_not_disturb_the_others():
    s = session_state.AnalysisSession()
    assert s.library.store is None
    assert s.library.source_path is None
    assert s.library.last_match is None
    # パーサ別のスロットは無い（B-2b で消した）。残りのスロットが消えていないこと。
    for name in ("curation", "pathway"):
        assert getattr(s, name) is not None
    assert s.dataset is None


def test_an_explicit_path_is_returned_as_is(tmp_path):
    explicit = tmp_path / "explicit.msp"
    explicit.write_text("NAME: x\n", encoding="utf-8")
    assert path_resolvers.resolve_library_path(str(explicit)) == str(explicit)


def test_no_path_means_none_even_if_libraries_are_around(tmp_path, monkeypatch):
    """`.dbs` / `.msp` が手元にあっても探さない（旧版は `*_Loaded.msp2.dbs` → `*.msp` と探した）。"""
    (tmp_path / "P_Loaded.msp2.dbs").write_bytes(b"PK\x03\x04")
    (tmp_path / "lib.msp").write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    assert path_resolvers.resolve_library_path() is None


def test_an_explicit_path_that_does_not_exist_is_an_error(tmp_path, monkeypatch):
    """以前は探索へ黙って落ちていた（別のライブラリを掴む）。"""
    (tmp_path / "lib.msp").write_text("NAME: x\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    with pytest.raises(LibraryPathError) as exc:
        path_resolvers.resolve_library_path(str(tmp_path / "typo.msp"))
    assert exc.value.code == "LIBRARY_NOT_FOUND"


def test_library_load_hides_file_name(tmp_path, monkeypatch):
    import json

    from metabolomix.library import tools as library_tools
    from tests.config_fixtures import use_config
    from tests.dbs_fixture import record, write_dbs

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    session_state.session = session_state.AnalysisSession()
    lib = tmp_path / "Zebra_SecretLab_VS21.msp2.dbs"
    write_dbs(lib, [("lib", [record(mz=500.0, ion_mode=0)], ["k"])])
    shown = json.loads(library_tools.library_load(file_path=str(lib)))
    assert shown["file"] == lib.name  # 既定は今までどおり
    raw = library_tools.library_load(file_path=str(lib), library_names=False)
    out = json.loads(raw)
    assert out["status"] == "success" and out["source_sha256"] and out["record_count"] == 1
    assert "Zebra" not in raw and "SecretLab" not in raw and "file" not in out
    bad = json.loads(library_tools.library_load(file_path=str(lib), library_names="x"))
    assert bad["status"] == "error"
    session_state.session = session_state.AnalysisSession()
