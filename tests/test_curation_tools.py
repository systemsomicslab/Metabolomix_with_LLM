import json

import pytest

from lipidmix.core import mcp_core, session_state
from lipidmix.curation import flags
from lipidmix.library import store as library_store
from lipidmix.tools import curation_tools
from tests.curation_fixtures import write_alignment_set


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
