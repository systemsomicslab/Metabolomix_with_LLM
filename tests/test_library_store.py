"""参照ライブラリの永続 store。"""
import textwrap

import pytest

from lipidmix.library import store

_MSP = textwrap.dedent("""\
    NAME: A
    PRECURSORMZ: 100.0
    IONMODE: Positive
    Num Peaks: 2
    50.0 999
    80.0 500

    NAME: B
    PRECURSORMZ: 100.005
    IONMODE: Positive
    Num Peaks: 1
    50.0 999

    NAME: C
    PRECURSORMZ: 100.0
    IONMODE: Negative
    Num Peaks: 1
    50.0 999

    NAME: D
    PRECURSORMZ: 200.0
    IONMODE: Positive
    Num Peaks: 1
    50.0 999
""")


@pytest.fixture()
def library(tmp_path, monkeypatch):
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "lib.msp"
    path.write_text(_MSP, encoding="utf-8")
    return path


def test_the_store_is_built_once_and_reopened(library):
    s1 = store.open_store(library)
    built_at = store.store_path_for(library).stat().st_mtime_ns
    assert s1.record_count == 4
    s1.close()

    s2 = store.open_store(library)
    assert store.store_path_for(library).stat().st_mtime_ns == built_at  # 作り直していない
    assert s2.record_count == 4
    s2.close()


def test_a_changed_source_gets_a_different_store(library, tmp_path):
    first = store.store_path_for(library)
    library.write_text(_MSP + "\nNAME: E\nPRECURSORMZ: 300.0\nNum Peaks: 0\n", encoding="utf-8")
    assert store.store_path_for(library) != first


def test_candidates_are_filtered_by_mz_window_and_polarity(library):
    s = store.open_store(library)
    try:
        hits = s.candidates(100.0, mz_tol=0.01, ion_mode="positive")
        assert sorted(r["name"] for r in hits) == ["A", "B"]   # C は極性違い、D は窓の外
        assert hits[0]["spectrum"] == [[50.0, 999.0], [80.0, 500.0]]
    finally:
        s.close()


def test_omitting_the_ion_mode_keeps_both_polarities(library):
    s = store.open_store(library)
    try:
        assert sorted(r["name"] for r in s.candidates(100.0, mz_tol=0.01)) == ["A", "B", "C"]
    finally:
        s.close()


def test_the_summary_reports_what_was_loaded(library):
    s = store.open_store(library)
    try:
        summary = s.summary()
        assert summary["record_count"] == 4
        assert summary["ion_modes"] == {"positive": 3, "negative": 1}
        assert summary["source_sha256"][:2].isalnum()
        assert summary["search_params"] is None       # .msp には同梱されない
    finally:
        s.close()


def test_library_id_and_record_index_round_trip(library):
    """record.RECORD_FIELDS は 13 フィールドが「採点エンジンの入力契約」（spec）。
    library_id / record_index を落とすと、どの元レコードが根拠かを辿れなくなる。"""
    s = store.open_store(library)
    try:
        hits = {r["name"]: r for r in s.candidates(100.0, mz_tol=0.01)}
        assert sorted(hits) == ["A", "B", "C"]
        # .msp は library_id を持たない（単一ファイル内なので区別が要らない）
        assert all(r["library_id"] is None for r in hits.values())
        # record_index は .msp 内の出現順の通し番号（A, B, C, D の順）
        assert hits["A"]["record_index"] == 0
        assert hits["B"]["record_index"] == 1
        assert hits["C"]["record_index"] == 2
    finally:
        s.close()


_MSP_WITH_RT = textwrap.dedent("""\
    NAME: RT_KNOWN
    PRECURSORMZ: 150.0
    IONMODE: Positive
    RETENTIONTIME: 5.0
    Num Peaks: 1
    50.0 999

    NAME: RT_UNKNOWN
    PRECURSORMZ: 150.0
    IONMODE: Positive
    Num Peaks: 1
    50.0 999

    NAME: RT_FAR
    PRECURSORMZ: 150.0
    IONMODE: Positive
    RETENTIONTIME: 50.0
    Num Peaks: 1
    50.0 999
""")


_MSP_CLASSES = textwrap.dedent("""\
    NAME: A1
    PRECURSORMZ: 100.0
    IONMODE: Positive
    COMPOUNDCLASS: LPE
    Num Peaks: 1
    50.0 999

    NAME: A2
    PRECURSORMZ: 100.0
    IONMODE: Positive
    COMPOUNDCLASS: LPE
    Num Peaks: 1
    50.0 999

    NAME: B1
    PRECURSORMZ: 100.0
    IONMODE: Positive
    COMPOUNDCLASS: PC
    Num Peaks: 1
    50.0 999

    NAME: C1
    PRECURSORMZ: 100.0
    IONMODE: Positive
    Num Peaks: 1
    50.0 999
""")


def test_compound_class_counts_reports_top_classes_descending(tmp_path, monkeypatch):
    """レコードに compound_class 未設定（None）が混ざっても集計から除外される。"""
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "classes.msp"
    path.write_text(_MSP_CLASSES, encoding="utf-8")

    s = store.open_store(path)
    try:
        assert s.compound_class_counts(top_n=10) == [
            {"name": "LPE", "count": 2},
            {"name": "PC", "count": 1},
        ]
        assert s.compound_class_counts(top_n=1) == [{"name": "LPE", "count": 2}]
    finally:
        s.close()


def test_compound_class_counts_does_not_rescan_the_source_on_cache_hit(tmp_path, monkeypatch):
    """`library_load` を 2 回呼んでも、2 回目はキャッシュを開くだけで元ファイルを
    再走査してはいけない（SQLite store を作った目的そのものが壊れるため）。"""
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "classes.msp"
    path.write_text(_MSP_CLASSES, encoding="utf-8")

    s1 = store.open_store(path)
    s1.close()

    def _boom(*_args, **_kwargs):
        raise AssertionError("iter_records が呼ばれた（cache hit のはずが元ファイルを再走査した）")

    monkeypatch.setattr(store.msp_reader, "iter_records", _boom)
    monkeypatch.setattr(store.dbs_reader, "iter_records", _boom)

    s2 = store.open_store(path)  # rebuild=False: キャッシュを開くだけのはず
    try:
        assert s2.compound_class_counts(top_n=10) == [
            {"name": "LPE", "count": 2},
            {"name": "PC", "count": 1},
        ]
    finally:
        s2.close()


def test_rt_filter_keeps_null_rt_but_drops_rt_outside_tolerance(tmp_path, monkeypatch):
    """brief 名指しの要件: `.msp` の RT は別の LC 条件で測られていることがあるので、
    RT を持たないレコードを rt 指定クエリで落としてはいけない
    （`AND (rt IS NULL OR ABS(rt - ?) <= ?)`）。RT が窓の外のレコードは落ちることも
    同じテストで押さえ、`rt IS NULL` だけが特別扱いされることを示す。"""
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "rt_lib.msp"
    path.write_text(_MSP_WITH_RT, encoding="utf-8")

    s = store.open_store(path)
    try:
        hits = s.candidates(150.0, mz_tol=0.01, rt=5.0, rt_tol=0.5)
        assert sorted(r["name"] for r in hits) == ["RT_KNOWN", "RT_UNKNOWN"]  # RT_FAR は窓の外
    finally:
        s.close()


_MSP_MISSING_PRECURSOR = textwrap.dedent("""\
    NAME: NO_PRECURSOR
    IONMODE: Positive
    Num Peaks: 1
    50.0 999

    NAME: HAS_PRECURSOR
    PRECURSORMZ: 100.0
    IONMODE: Positive
    Num Peaks: 1
    50.0 999
""")


def test_a_record_without_precursor_mz_is_skipped_not_fatal(tmp_path, monkeypatch):
    """Important 3 の再発防止: `record.precursor_mz REAL NOT NULL`（store のスキーマ）
    に対し、reader は PRECURSORMZ 欠損を契約どおり `None` に潰すだけで例外にしない
    （`msp.py` の `_read_record`）。以前は 1 レコードの欠損で
    `IntegrityError: NOT NULL constraint failed: record.precursor_mz` となり
    ライブラリ全体が構築できなかった。欠損レコードだけ読み飛ばし、件数を meta に残す。"""
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "missing_precursor.msp"
    path.write_text(_MSP_MISSING_PRECURSOR, encoding="utf-8")

    s = store.open_store(path)
    try:
        assert s.record_count == 1  # NO_PRECURSOR は読み飛ばされる
        assert s.summary()["skipped_no_precursor_mz"] == 1
        hits = s.candidates(100.0, mz_tol=0.01)
        assert [r["name"] for r in hits] == ["HAS_PRECURSOR"]
    finally:
        s.close()


def test_skipped_no_precursor_mz_defaults_to_zero_for_a_clean_library(library):
    s = store.open_store(library)
    try:
        assert s.summary()["skipped_no_precursor_mz"] == 0
    finally:
        s.close()


def test_candidates_rejects_an_rt_without_a_tolerance(library):
    """Minor 11 の再発防止: `candidates(rt=..., rt_tol=None)` は SQL の
    `ABS(rt - ?) <= NULL` が NULL（偽）になるため `rt IS NULL` の行しか
    返さない黙った縮退を起こす。シグネチャ上 `rt_tol` は省略できるように
    見えるので、呼び出し側の取り違えを実行時に検出する。"""
    s = store.open_store(library)
    try:
        with pytest.raises(ValueError):
            s.candidates(100.0, mz_tol=0.01, rt=5.0, rt_tol=None)
    finally:
        s.close()


def test_candidates_ion_mode_comparison_is_case_insensitive(library):
    """Critical 1 の再発防止: store の `ion_mode` 列は小文字（`record.ION_MODES`）
    だが、呼び出し側（`peak_verification._spectral_match_for_feature`）は
    PAI2 由来の `IonMode.Positive.name`（`"Positive"`、大文字始まり）をそのまま
    渡す。以前は SQLite の既定 BINARY 照合で一致せず常に 0 件になっていた。"""
    s = store.open_store(library)
    try:
        hits = s.candidates(100.0, mz_tol=0.01, ion_mode="Positive")
        assert sorted(r["name"] for r in hits) == ["A", "B"]
        hits_upper = s.candidates(100.0, mz_tol=0.01, ion_mode="POSITIVE")
        assert sorted(r["name"] for r in hits_upper) == ["A", "B"]
    finally:
        s.close()


# --------------------------------------------------------------------------
# 大容量 `.msp`（研究室ライブラリ、pos ≈ 1.2 GB）。
# --------------------------------------------------------------------------
_MSP_NO_ION_MODE = textwrap.dedent("""\
    NAME: NO_MODE
    PRECURSORMZ: 100.0
    Num Peaks: 1
    50.0 999

    NAME: POS
    PRECURSORMZ: 100.0
    IONMODE: Positive
    Num Peaks: 1
    50.0 999

    NAME: NEG
    PRECURSORMZ: 100.0
    IONMODE: Negative
    Num Peaks: 1
    50.0 999
""")


def test_records_without_an_ion_mode_survive_the_polarity_filter(tmp_path, monkeypatch):
    """極性ごとに分かれた `.msp` は IONMODE 欄を持たないことがある。NULL を
    `ion_mode = ?` で弾くと候補が黙って 0 件になる。極性不明は候補に残す。"""
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "lib.msp"
    path.write_text(_MSP_NO_ION_MODE, encoding="utf-8")
    s = store.open_store(path)
    try:
        hits = s.candidates(100.0, mz_tol=0.01, ion_mode="positive")
        assert sorted(r["name"] for r in hits) == ["NO_MODE", "POS"]
        assert s.summary()["records_without_ion_mode"] == 1
    finally:
        s.close()


def test_non_utf8_lines_are_recorded_in_the_summary(tmp_path, monkeypatch):
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "lib.msp"
    path.write_bytes("NAME: グルタミン酸\nPRECURSORMZ: 148.06\nNum Peaks: 0\n".encode("cp932"))
    s = store.open_store(path)
    try:
        assert s.summary()["non_utf8_lines"] == 1
        assert s.candidates(148.06, mz_tol=0.01)[0]["name"] == "グルタミン酸"
    finally:
        s.close()


def _count_hashing(monkeypatch):
    calls = []
    real = store.source_sha256

    def spy(path):
        calls.append(path)
        return real(path)

    monkeypatch.setattr(store, "source_sha256", spy)
    return calls


def test_reopening_an_unchanged_source_does_not_rehash_it(library, monkeypatch):
    """1.2 GB の sha256 を `library_load` のたびに払わない。サイズと更新時刻が
    前回と同じなら、前回の sha256 を使い回す。"""
    calls = _count_hashing(monkeypatch)
    store.open_store(library).close()
    assert len(calls) == 1
    s = store.open_store(library)
    try:
        assert len(calls) == 1
        assert s.record_count == 4
    finally:
        s.close()


def test_a_modified_source_is_rehashed_and_rebuilt(library, monkeypatch):
    calls = _count_hashing(monkeypatch)
    store.open_store(library).close()
    library.write_text(_MSP + "\nNAME: E\nPRECURSORMZ: 300.0\nNum Peaks: 0\n", encoding="utf-8")
    s = store.open_store(library)
    try:
        assert len(calls) == 2
        assert s.record_count == 5
    finally:
        s.close()


def test_rebuild_rehashes_even_when_the_source_looks_unchanged(library, monkeypatch):
    """`rebuild=True` は「元ファイルが壊れている疑い」のときに使う。記憶した
    sha256 を信じない。"""
    calls = _count_hashing(monkeypatch)
    store.open_store(library).close()
    store.open_store(library, rebuild=True).close()
    assert len(calls) == 2


def test_a_corrupt_digest_index_is_ignored(library, monkeypatch):
    store.open_store(library).close()
    (store.cache_dir() / store.DIGEST_INDEX_NAME).write_text("{not json", encoding="utf-8")
    s = store.open_store(library)
    try:
        assert s.record_count == 4
    finally:
        s.close()


def test_the_cli_prebuilds_the_store_for_the_env_var_library(tmp_path, monkeypatch, capsys):
    """MCP クライアントのタイムアウトを避けるため、初回の構築は手元で済ませられる。"""
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    lib = tmp_path / "outside" / "lab_neg.msp"
    lib.parent.mkdir()
    lib.write_text(_MSP, encoding="utf-8")
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(lib))

    assert store.main(["--ion-mode", "negative"]) == 0

    out = capsys.readouterr().out
    assert "lab_neg.msp" in out and "4" in out
    assert "outside" not in out                     # 置き場所は出さない
    assert store.store_path_for(lib).exists()


def test_the_cli_reports_resolution_errors_without_a_traceback(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(tmp_path / "gone.msp"))
    assert store.main(["--ion-mode", "negative"]) == 2
    assert "MSP_ENV_NOT_FOUND" in capsys.readouterr().err


def test_record_by_scan_id_uses_the_mz_index_and_the_scan_id(library):
    s = store.open_store(library)
    try:
        hit = s.record_by_scan_id(1, library_id=None, precursor_mz=100.004, mz_tol=0.05)
        assert hit["name"] == "B"
        assert hit["record_index"] == 1
        assert s.record_by_scan_id(3, library_id=None, precursor_mz=100.0, mz_tol=0.05) is None
    finally:
        s.close()


def test_record_by_scan_id_filters_by_library_id_when_the_store_has_one(library):
    s = store.open_store(library)
    try:
        # .msp 由来は library_id が NULL。NULL は「区別不要」として残す。
        assert s.record_by_scan_id(0, library_id="Msp1_lib", precursor_mz=100.0,
                                   mz_tol=0.05)["name"] == "A"
    finally:
        s.close()


def test_record_by_scan_id_rejects_a_row_whose_library_id_differs(library):
    s = store.open_store(library)
    try:
        # .dbs 由来の store を模して、行に非 NULL の library_id を入れる。
        s._conn.execute("UPDATE record SET library_id = 'LibA' WHERE record_index = 0")
        assert s.record_by_scan_id(0, library_id="LibB", precursor_mz=100.0, mz_tol=0.05) is None
        assert s.record_by_scan_id(0, library_id="LibA", precursor_mz=100.0,
                                   mz_tol=0.05)["name"] == "A"
    finally:
        s.close()


def test_library_id_from_annotator_strips_the_trailing_counter():
    assert store.library_id_from_annotator("Msp20260116160945_NCDK_dev_1") == "Msp20260116160945_NCDK_dev"
    assert store.library_id_from_annotator("plain") == "plain"
    assert store.library_id_from_annotator(None) is None
