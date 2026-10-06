"""キュレーションが、MS-DIAL Console 由来の照合結果から参照スペクトルを引けること。

Console の照合結果の AnnotatorID は `.dbs` のライブラリ名と一致しない（LBM なら
ファイルのパス）。ライブラリ名を AnnotatorID から文字列で推測していたため、Console で
作ったアラインメントでは参照が 1 件も引けず、ミラー図も再採点も出なかった。
"""
from __future__ import annotations

from metabolomix.library import store


def _console_store(tmp_path, monkeypatch):
    from tests.dbs_fixture import record, write_dbs
    monkeypatch.setenv(store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = write_dbs(tmp_path / "Project_Loaded.msp2.dbs", [
        ("LbmDB", [record(name="PC 34:1", mz=760.585, inchikey="")],
         ["C:/lib/NCDK_conventional.lbm2"]),
    ])
    return store.open_store(path)


_MATCH = {"library_id": 0, "annotator_id": "C:/lib/NCDK_conventional.lbm2",
          "inchikey": "", "name": "PC 34:1", "total_score": 0.9}


def test_review_evidence_resolves_a_console_lbm_reference(tmp_path, monkeypatch):
    from metabolomix.curation import evidence
    s = _console_store(tmp_path, monkeypatch)
    try:
        reference = evidence._reference(s, _MATCH, 760.586)
        assert reference is not None and reference["name"] == "PC 34:1"
        assert reference["spectrum"]
    finally:
        s.close()


def test_suggestion_resolves_a_console_lbm_reference(tmp_path, monkeypatch):
    from metabolomix.curation import candidates
    s = _console_store(tmp_path, monkeypatch)
    try:
        record = candidates._resolve_match(s, _MATCH, 760.586)
        assert record is not None and record["name"] == "PC 34:1"
    finally:
        s.close()
