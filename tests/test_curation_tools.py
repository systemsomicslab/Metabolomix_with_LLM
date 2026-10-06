import json

import pytest

from metabolomix.core import mcp_core, session_state
from metabolomix.curation import flags, suggest
from metabolomix.library import store as library_store
from metabolomix.tools import curation_tools
from tests.curation_fixtures import write_alignment_set, write_suggest_set


@pytest.fixture()
def ready(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg")
    monkeypatch.setattr(mcp_core, "DATA_DIR", paths["arf2"].parent)
    session_state.session.__init__()
    session_state.session.library.store = library_store.open_store(paths["msp"])
    yield paths
    session_state.session.library.store.close()
    session_state.session.__init__()


def test_review_requires_a_loaded_library(tmp_path, monkeypatch):
    paths = write_alignment_set(tmp_path / "neg")
    monkeypatch.setattr(mcp_core, "DATA_DIR", paths["arf2"].parent)
    session_state.session.__init__()
    body = json.loads(curation_tools.curation_review())
    assert body["error"]["code"] == "missing_state"
    assert body["error"]["required_tools"] == ["library_load"]


def test_review_returns_summary_and_writes_the_viewer(ready):
    body = json.loads(curation_tools.curation_review())
    assert body["counts"]["ok"] == 1
    assert body["table"].splitlines()[0].startswith("spot_id\tname")
    assert body["html_path"].endswith(".html")
    assert "EIC" not in body and "spots" not in body       # 座標を LLM に返さない
    assert session_state.session.curation.last_review_id == body["review_id"]


def test_submit_via_pasted_text_records_flags(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    text = flags.build_submission_text(review_id, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    body = json.loads(curation_tools.curation_submit(submission_text="送ります\n" + text))
    assert body["recorded"] == 1
    listing = json.loads(curation_tools.curation_flags())
    assert listing["table"].splitlines()[1].split("\t")[:2] == ["1", "wrong"]


def test_submit_rejects_spots_outside_the_review_without_writing(ready):
    review_id = json.loads(curation_tools.curation_review(name_contains="34:1"))["review_id"]
    body = json.loads(curation_tools.curation_submit(
        review_id=review_id, flags=[{"spot_id": 0, "flag": "wrong"}, {"spot_id": 1, "flag": "wrong"}]))
    assert body["status"] == "error"
    assert json.loads(curation_tools.curation_flags())["n_flags"] == 0


def test_submit_rejects_unknown_review(ready):
    body = json.loads(curation_tools.curation_submit(review_id="cr-nope", flags=[{"spot_id": 0, "flag": "wrong"}]))
    assert body["status"] == "error"


def test_submit_rejects_a_submission_naming_the_same_spot_twice(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    body = json.loads(curation_tools.curation_submit(
        review_id=review_id,
        flags=[{"spot_id": 0, "flag": "wrong"}, {"spot_id": 0, "flag": "suspect"}]))
    assert body["status"] == "error"
    assert json.loads(curation_tools.curation_flags())["n_flags"] == 0


def test_submit_rejects_when_the_alignment_changed_since_the_review(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    with open(ready["arf2"], "ab") as handle:
        handle.write(b"\x00")          # .arf2 の sha256 をレビュー時から変える
    body = json.loads(curation_tools.curation_submit(
        review_id=review_id, flags=[{"spot_id": 0, "flag": "wrong"}]))
    assert body["status"] == "error"
    assert json.loads(curation_tools.curation_flags())["n_flags"] == 0


def test_view_data_pages_the_saved_review(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    body = json.loads(curation_tools.curation_view_data(review_id, 0))
    assert body["n_pages"] == 1
    assert {s["spot_id"] for s in body["spots"]} == {0, 1}
    assert body["review"]["review_id"] == review_id
    assert "spots" not in body["review"]


def test_too_many_spots_is_refused(ready, monkeypatch):
    monkeypatch.setattr(curation_tools.evidence, "MAX_SPOTS", 1)
    body = json.loads(curation_tools.curation_review())
    assert body["status"] == "error"
    assert "ontology" in body["message"]


import asyncio

import server


def test_review_tool_advertises_the_viewer_resource():
    tools = {t.name: t for t in asyncio.run(server.mcp.list_tools())}
    assert tools["curation_review"].meta["ui"]["resourceUri"] == curation_tools.VIEWER_URI
    assert tools["curation_view_data"].meta["ui"]["visibility"] == ["app"]


def test_viewer_resource_is_html_for_mcp_apps():
    resources = {str(r.uri): r for r in asyncio.run(server.mcp.list_resources())}
    assert resources[curation_tools.VIEWER_URI].mimeType == "text/html;profile=mcp-app"


# ---------- I2: TSV の行数上限と総数 ----------

def test_review_table_is_capped_and_reports_the_total(ready):
    body = json.loads(curation_tools.curation_review(max_rows=0))
    assert body["n_table_rows_total"] == 1
    assert body["n_table_rows_shown"] == 0
    assert len(body["table"].splitlines()) == 1          # 見出しだけ
    assert "html_path" in body["table_note"]


def test_review_table_with_max_rows_one_returns_one_row(ready):
    flag_dir = flags.curation_dir(ready["arf2"])
    flags.FlagStore(flag_dir).append([{"spot_id": 0, "flag": "suspect"}],
                                     alignment=flags.alignment_key(ready["arf2"]),
                                     review_id="r", source="user")
    body = json.loads(curation_tools.curation_review(max_rows=1))
    assert body["n_table_rows_total"] == 2
    assert body["n_table_rows_shown"] == 1
    assert len(body["table"].splitlines()) == 2


# ---------- I4: 再起動・データフォルダ変更の後でも貼った文で送れる ----------

def test_submit_after_a_restart_uses_the_arf2_path_in_the_text(ready, tmp_path, monkeypatch):
    first = json.loads(curation_tools.curation_review())
    saved = json.loads((flags.curation_dir(ready["arf2"])
                        / f"review-{first['review_id']}.json").read_text(encoding="utf-8"))
    text = flags.build_submission_text(first["review_id"],
                                       [{"spot_id": 1, "flag": "wrong", "note": ""}],
                                       saved["arf2_path"])
    library = session_state.session.library.store
    session_state.session.__init__()                       # サーバ再起動を模す
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.setattr(mcp_core, "DATA_DIR", elsewhere)
    try:
        body = json.loads(curation_tools.curation_submit(submission_text=text))
    finally:
        session_state.session.library.store = library
    assert body["status"] == "ok", body
    assert body["recorded"] == 1


def test_submit_after_a_restart_accepts_file_path(ready, tmp_path, monkeypatch):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    library = session_state.session.library.store
    session_state.session.__init__()
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path / "nowhere")
    try:
        missing = json.loads(curation_tools.curation_submit(
            review_id=review_id, flags=[{"spot_id": 1, "flag": "wrong"}]))
        body = json.loads(curation_tools.curation_submit(
            review_id=review_id, flags=[{"spot_id": 1, "flag": "wrong"}],
            file_path=str(ready["arf2"])))
    finally:
        session_state.session.library.store = library
    assert missing["status"] == "error"
    assert "file_path" in missing["message"]
    assert "curation_review をやり直" not in missing["message"]
    assert body["status"] == "ok", body


# ---------- M1: review_id の形式を path に使う前に検める ----------

def test_submit_and_view_data_reject_a_malformed_review_id(ready):
    for body in (json.loads(curation_tools.curation_submit(
                     review_id="../../x", flags=[{"spot_id": 0, "flag": "wrong"}])),
                 json.loads(curation_tools.curation_view_data("../../x"))):
        assert body["status"] == "error"
        assert "review_id" in body["message"]


# ---------- M3: 知らない file_ids は先に弾く ----------

def test_review_rejects_unknown_file_ids(ready):
    body = json.loads(curation_tools.curation_review(file_ids=[0, 77]))
    assert body["status"] == "error"
    assert body["missing_file_ids"] == [77]


# ---------- I6: 壊れた flags.jsonl ----------

def _break_flags_file(arf2):
    store = flags.FlagStore(flags.curation_dir(arf2))
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=flags.alignment_key(arf2),
                 review_id="r", source="user")
    with open(store.path, "a", encoding="utf-8") as handle:
        handle.write('{"spot_id": 1, "fl')
    return store.path


def test_malformed_flags_file_is_reported_by_every_curation_tool(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    path = _break_flags_file(ready["arf2"])
    before = path.read_text(encoding="utf-8")
    for body in (json.loads(curation_tools.curation_flags()),
                 json.loads(curation_tools.curation_review()),
                 json.loads(curation_tools.curation_submit(
                     review_id=review_id, flags=[{"spot_id": 0, "flag": "wrong"}]))):
        assert body["status"] == "error", body
        assert body["flags_file"] == str(path)
        assert body["line"] == 2
    assert path.read_text(encoding="utf-8") == before       # 壊れたまま追記しない


# ---------- I7 / M10: warnings ----------

def test_review_warns_about_orphaned_flags(ready):
    flags.FlagStore(flags.curation_dir(ready["arf2"])).append(
        [{"spot_id": 1, "flag": "wrong"}],
        alignment={"alignment_file": ready["arf2"].name, "alignment_sha256": "0" * 64},
        review_id="r", source="user")
    body = json.loads(curation_tools.curation_review())
    assert any("1 件" in w and "以前の版" in w for w in body["warnings"]), body["warnings"]


def test_review_warns_about_missing_siblings(ready):
    folder = ready["arf2"].parent
    (folder / "AlignmentResult_x.EIC.aef").unlink()
    (folder / "AlignmentResult_x.dcl").unlink()
    body = json.loads(curation_tools.curation_review())
    assert any("dcl" in w and "eic" in w for w in body["warnings"]), body["warnings"]


# ---------- M7 ----------

def test_review_trend_summary_has_outlier_counts(ready, monkeypatch):
    real = curation_tools.review.trend.fit_trends

    def fake(points, th):
        result = real(points, th)
        result["classes"]["PC"] = {"n": 5, "r2": 0.9, "coef": [0.0], "scale": 0.02, "n_outliers": 2}
        return result

    monkeypatch.setattr(curation_tools.review.trend, "fit_trends", fake)
    body = json.loads(curation_tools.curation_review())
    assert body["trend"]["PC"] == {"n": 5, "r2": 0.9, "n_outliers": 2}


# --- 送信で `_tags.xml` の Misannotation を反映（ユーザー決定 2026-09-29）:
# 「間違い」→ 付ける、取消（clear）→ 外す、「疑わしい」→ 触らない。flags.jsonl が正本。---

from metabolomix.msdial import tags as msdial_tags


def _submit(review_id, entries):
    text = flags.build_submission_text(review_id, entries)
    return json.loads(curation_tools.curation_submit(submission_text=text))


def test_submitting_wrong_sets_misannotation_in_the_alignment_tags(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    body = _submit(review_id, [{"spot_id": 1, "flag": "wrong", "note": "自動: X"},
                               {"spot_id": 0, "flag": "suspect", "note": ""}])
    tag_path = msdial_tags.alignment_tag_path(ready["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {1: frozenset({3})}
    assert body["tags_xml"]["added"] == [1] and body["tags_xml"]["removed"] == []
    assert body["tags_xml"]["path"] == str(tag_path)
    assert "MS-DIAL" in body["tags_xml"]["note"]


def test_clearing_removes_misannotation(ready):
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    _submit(review_id, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    body = _submit(review_id, [{"spot_id": 1, "flag": "clear", "note": ""}])
    assert msdial_tags.parse_tag_file(msdial_tags.alignment_tag_path(ready["arf2"]))["peaks"] == {}
    assert body["tags_xml"]["removed"] == [1]


def test_existing_tags_file_is_backed_up_before_writing(ready):
    tag_path = msdial_tags.alignment_tag_path(ready["arf2"])
    msdial_tags.update_alignment_tag(tag_path, tag_id=1, add=[0], remove=[])     # Confirmed を人が付けていた
    before = tag_path.read_bytes()
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    body = _submit(review_id, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    backup = body["tags_xml"]["backup"]
    assert backup and open(backup, "rb").read() == before
    assert flags.curation_dir(ready["arf2"]) / "tags-backup" in [p for p in __import__("pathlib").Path(backup).parents]
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {0: frozenset({1}), 1: frozenset({3})}


def test_a_failed_tags_write_keeps_the_recorded_flags(ready):
    msdial_tags.alignment_tag_path(ready["arf2"]).mkdir()     # ファイルの位置にフォルダ＝書けない
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    body = _submit(review_id, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    assert body["status"] == "ok" and body["recorded"] == 1
    assert "error" in body["tags_xml"]
    listing = json.loads(curation_tools.curation_flags())
    assert listing["table"].splitlines()[1].split("\t")[:2] == ["1", "wrong"]


def test_a_corrupt_tags_file_is_reported_and_left_untouched(ready):
    tag_path = msdial_tags.alignment_tag_path(ready["arf2"])
    tag_path.write_text("<PeakSpotTags><Peaks>", encoding="utf-8")      # 途中で切れた XML
    review_id = json.loads(curation_tools.curation_review())["review_id"]
    body = _submit(review_id, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    assert body["recorded"] == 1 and "error" in body["tags_xml"]
    assert tag_path.read_text(encoding="utf-8") == "<PeakSpotTags><Peaks>"


def test_submit_is_annotated_as_destructive_because_it_rewrites_the_tags_file():
    import asyncio
    import server
    tools = {t.name: t for t in asyncio.run(server.mcp.list_tools())}
    annotations = tools["curation_submit"].annotations
    assert annotations.destructiveHint is True and annotations.readOnlyHint is False
    assert tools["curation_review"].annotations.destructiveHint is False


@pytest.fixture()
def suggest_env(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_suggest_set(tmp_path / "neg")
    monkeypatch.setattr(mcp_core, "DATA_DIR", paths["arf2"].parent)
    session_state.session.__init__()
    s = library_store.open_store(paths["msp"])
    session_state.session.library.store = s
    yield paths
    s.close()
    session_state.session.__init__()


def test_suggest_requires_a_review(suggest_env):
    out = json.loads(curation_tools.curation_suggest(file_path=str(suggest_env["arf2"])))
    assert out["error"]["code"] == "missing_state" and out["error"]["required_tools"] == ["curation_review"]


def test_suggest_requires_the_library(suggest_env):
    session_state.session.library.store = None
    out = json.loads(curation_tools.curation_suggest(file_path=str(suggest_env["arf2"])))
    assert out["error"]["required_tools"] == ["library_load"]


def test_suggest_then_submit_assign_and_redundant(suggest_env):
    arf2 = str(suggest_env["arf2"])
    json.loads(curation_tools.curation_review(file_path=arf2))
    out = json.loads(curation_tools.curation_suggest(file_path=arf2))
    sid = out["suggestion_id"]
    assert out["counts"]["targets"]["unannotated"] == 3
    assert out["table"].splitlines()[0].split("\t") == suggest.TSV_COLUMNS
    assert out["html_path"].endswith(f"suggest-{sid}.html")
    assert "spots" not in out                                     # 座標・EIC は戻り値に入れない
    assert "path" in out["library"]                               # ライブラリの出所（sha256 だけでは分からない）
    text = 'CURATION_SUBMIT ' + json.dumps({"review_id": sid, "arf2_path": arf2, "flags": [
        {"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "sum", "note": ""},
        {"spot_id": 3, "flag": "redundant", "candidate": "R1", "note": ""}]})
    submitted = json.loads(curation_tools.curation_submit(submission_text=text))
    assert submitted["status"] == "ok" and submitted["recorded"] == 2
    assert submitted["n_assign"] == 1 and submitted["n_redundant"] == 1
    assert submitted["tags_xml"]["added"] == [] and submitted["tags_xml"]["removed"] == []   # _tags.xml は不変
    table = json.loads(curation_tools.curation_flags(file_path=arf2))["table"].splitlines()
    assert table[0] == "spot_id\tflag\tname\tof\tnote\tsource\tts"
    assert table[1].split("\t")[:4] == ["2", "assign", "PE 36:2", ""]
    assert table[2].split("\t")[:4] == ["3", "redundant", "", "0"]


def test_submit_rejects_a_forged_candidate(suggest_env):
    arf2 = str(suggest_env["arf2"])
    json.loads(curation_tools.curation_review(file_path=arf2))
    sid = json.loads(curation_tools.curation_suggest(file_path=arf2))["suggestion_id"]
    out = json.loads(curation_tools.curation_submit(review_id=sid, file_path=arf2, flags=[
        {"spot_id": 2, "flag": "assign", "candidate": "L7", "level": "sum"}]))
    assert out["status"] == "error"
    flags_path = suggest_env["arf2"].parent / "curation" / "flags.jsonl"
    assert not flags_path.exists() or '"assign"' not in flags_path.read_text(encoding="utf-8")


def test_suggest_survives_a_corrupt_review_file(suggest_env):
    arf2 = str(suggest_env["arf2"])
    json.loads(curation_tools.curation_review(file_path=arf2))
    (flags.curation_dir(suggest_env["arf2"]) / "review-cr-99999999-999999-ffff.json").write_text(
        "{not json", encoding="utf-8")
    out = json.loads(curation_tools.curation_suggest(file_path=arf2))
    assert "suggestion_id" in out and "error" not in out
