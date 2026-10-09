import http.client
import json
import threading

import pytest

from metabolomix.curation import flags, submission, submit_server
from metabolomix.library import store as library_store
from metabolomix.msdial import tags as msdial_tags
from metabolomix.core import mcp_core, session_state
from metabolomix.tools import curation_tools
from tests.curation_fixtures import build_saved_review, write_suggest_set


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now

    def advance(self, minutes: float):
        self.now += minutes * 60


@pytest.fixture()
def clock(monkeypatch):
    fake = FakeClock()
    monkeypatch.setattr(submit_server, "_clock", fake)
    return fake


def _post(port, path, body=None, *, raw=None, headers=None, method="POST"):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    try:
        data = raw if raw is not None else (None if body is None else json.dumps(body).encode("utf-8"))
        conn.request(method, path, body=data,
                     headers={"Content-Type": "text/plain;charset=UTF-8", **(headers or {})})
        res = conn.getresponse()
        text = res.read().decode("utf-8")
        return res.status, (json.loads(text) if text else None)
    finally:
        conn.close()


def test_register_starts_the_endpoint_and_later_registrations_share_it():
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    b = submit_server.register("suggest", "cs-20260101-000000-bbbb", "C:/x.arf2")
    assert a["port"] == b["port"] == submit_server.current_port()
    assert a["token"] != b["token"] and len(a["token"]) >= 40
    assert a["idle_timeout_min"] == 30
    status, body = _post(a["port"], "/v1/ping", {"token": b["token"]})
    assert status == 200
    assert body["status"] == "ok" and body["kind"] == "suggest"
    assert body["review_id"] == "cs-20260101-000000-bbbb" and body["expires_at"]


def test_the_last_unregister_stops_the_endpoint():
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    b = submit_server.register("review", "cr-20260101-000000-bbbb", "C:/x.arf2")
    submit_server.unregister(a["token"])
    assert submit_server.is_running()
    submit_server.unregister(b["token"])
    assert not submit_server.is_running() and submit_server.current_port() is None
    with pytest.raises(OSError):
        _post(a["port"], "/v1/ping", {"token": b["token"]})


def test_register_rejects_an_unknown_kind():
    with pytest.raises(ValueError):
        submit_server.register("other", "cr-20260101-000000-aaaa", "C:/x.arf2")


def test_ping_extends_the_deadline_and_idle_registrations_are_swept(clock):
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    clock.advance(29)
    assert _post(a["port"], "/v1/ping", {"token": a["token"]})[0] == 200     # heartbeat で延びる
    clock.advance(29)
    submit_server.sweep()
    assert submit_server.is_running()
    clock.advance(2)
    submit_server.sweep()
    assert not submit_server.is_running()


def test_an_expired_token_is_refused_before_the_sweep(clock):
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    clock.advance(31)
    assert _post(a["port"], "/v1/ping", {"token": a["token"]})[0] == 401


@pytest.mark.parametrize("raw", [b"[1, 2]", b'"token"', b'{"token": 5}', '{"token": "ｔｏｋｅｎ"}'.encode("utf-8")])
def test_non_object_bodies_and_bad_tokens_are_400_or_401(raw):
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    status, body = _post(a["port"], "/v1/ping", raw=raw)
    assert status in (400, 401) and body["status"] == "error"
    assert _post(a["port"], "/v1/ping", {"token": a["token"]})[0] == 200    # 受け口は生きている


def test_options_preflight_allows_private_network_access():
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    conn = http.client.HTTPConnection("127.0.0.1", a["port"], timeout=10)
    try:
        conn.request("OPTIONS", "/v1/submit", headers={"Access-Control-Request-Private-Network": "true"})
        res = conn.getresponse(); res.read()
    finally:
        conn.close()
    assert res.status == 204
    assert res.getheader("Access-Control-Allow-Origin") == "*"
    assert res.getheader("Access-Control-Allow-Private-Network") == "true"


@pytest.fixture()
def saved(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    built = build_saved_review(tmp_path / "neg")
    r = built["review"]
    built["endpoint"] = submit_server.register("review", r["review_id"], r["arf2_path"])
    return built


def _submit(saved, flag_entries, **overrides):
    ep, r = saved["endpoint"], saved["review"]
    body = {"token": ep["token"], "review_id": r["review_id"], "flags": flag_entries, **overrides}
    return _post(ep["port"], "/v1/submit", body)


def _n_rows(saved):
    return len(flags.FlagStore(flags.curation_dir(saved["review"]["arf2_path"])).rows())


def test_submit_records_as_the_user_and_updates_the_tags_file(saved):
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong", "note": "x"}])
    assert status == 200 and body["recorded"] == 1 and body["tags_xml"]["added"] == [1]
    rows = flags.FlagStore(flags.curation_dir(saved["review"]["arf2_path"])).rows()
    assert rows[-1]["source"] == "user"
    tag_path = msdial_tags.alignment_tag_path(saved["paths"]["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {1: frozenset({3})}


def test_submit_can_be_sent_in_several_parts(saved):
    assert _submit(saved, [{"spot_id": 1, "flag": "wrong", "note": ""}])[0] == 200
    assert _submit(saved, [{"spot_id": 0, "flag": "confirmed", "note": ""}])[0] == 200
    assert _n_rows(saved) == 2


def test_a_mismatched_review_id_is_403_and_writes_nothing(saved):
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}], review_id="cr-20000101-000000-abcd")
    assert status == 403 and body["status"] == "error"
    assert _n_rows(saved) == 0


def test_a_foreign_host_header_is_403(saved):
    ep = saved["endpoint"]
    status, _ = _post(ep["port"], "/v1/ping", {"token": ep["token"]},
                      headers={"Host": f"evil.example:{ep['port']}"})
    assert status == 403


def test_an_oversized_body_is_413(saved, monkeypatch):
    monkeypatch.setattr(submit_server, "MAX_BODY_BYTES", 64)
    status, _ = _submit(saved, [{"spot_id": 1, "flag": "wrong", "note": "x" * 100}])
    assert status == 413
    assert _n_rows(saved) == 0


def test_non_json_is_400_unknown_path_404_and_get_405(saved):
    ep = saved["endpoint"]
    assert _post(ep["port"], "/v1/submit", raw=b"not json")[0] == 400
    assert _post(ep["port"], "/v1/nope", {"token": ep["token"]})[0] == 404
    assert _post(ep["port"], "/v1/ping", method="GET")[0] == 405


def test_a_negative_content_length_is_400_without_blocking(saved):
    ep = saved["endpoint"]
    status, body = _post(ep["port"], "/v1/submit", raw=b"{}", headers={"Content-Length": "-1"})
    assert status == 400 and body["status"] == "error"
    assert _post(ep["port"], "/v1/ping", {"token": ep["token"]})[0] == 200


def test_an_invalid_flag_is_400_and_writes_nothing(saved):
    status, body = _submit(saved, [{"spot_id": 1, "flag": "nope"}])
    assert status == 400 and "flag" in body["message"]
    assert _n_rows(saved) == 0


def test_an_alignment_change_is_409(saved):
    with open(saved["paths"]["arf2"], "ab") as handle:
        handle.write(bytes([0]))
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert status == 409 and "curation_review" in body["message"]


def test_a_deleted_review_is_410(saved):
    r = saved["review"]
    (flags.curation_dir(r["arf2_path"]) / f"review-{r['review_id']}.json").unlink()
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert status == 410 and body["status"] == "error"


def test_an_unexpected_error_is_500_and_the_endpoint_survives(saved, monkeypatch):
    def boom(*args, **kwargs):
        raise RuntimeError("boom")
    monkeypatch.setattr(submission, "submit_flags", boom)
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert status == 500 and body["status"] == "error"
    ep = saved["endpoint"]
    assert _post(ep["port"], "/v1/ping", {"token": ep["token"]})[0] == 200


def test_finish_unregisters_the_token(saved):
    other = submit_server.register("review", "cr-20260101-000000-zzzz", "C:/x.arf2")   # 待受けを残す
    ep = saved["endpoint"]
    status, body = _post(ep["port"], "/v1/finish", {"token": ep["token"]})
    assert status == 200 and body == {"status": "ok"}
    assert _post(ep["port"], "/v1/ping", {"token": ep["token"]})[0] == 401
    assert _post(ep["port"], "/v1/ping", {"token": other["token"]})[0] == 200


def test_concurrent_submits_are_both_recorded(saved):
    results = []
    threads = [threading.Thread(target=lambda sid=sid: results.append(
        _submit(saved, [{"spot_id": sid, "flag": "wrong", "note": ""}])[0])) for sid in (0, 1)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=60)
    assert results == [200, 200]
    assert _n_rows(saved) == 2
    tag_path = msdial_tags.alignment_tag_path(saved["paths"]["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {0: frozenset({3}), 1: frozenset({3})}


def test_requests_write_nothing_to_stdout(saved, capfd):
    capfd.readouterr()
    ep = saved["endpoint"]
    _post(ep["port"], "/v1/ping", {"token": ep["token"]})
    _post(ep["port"], "/v1/nope", {"token": ep["token"]})
    _submit(saved, [{"spot_id": 1, "flag": "nope"}])
    _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert capfd.readouterr().out == ""


def test_a_suggestion_submit_over_http_records_assign_and_redundant(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_suggest_set(tmp_path / "neg")
    monkeypatch.setattr(mcp_core, "DATA_DIR", paths["arf2"].parent)
    session_state.session.__init__()
    lib = library_store.open_store(paths["msp"])
    session_state.session.library.store = lib
    try:
        arf2 = str(paths["arf2"])
        json.loads(curation_tools.curation_review(file_path=arf2))
        sid = json.loads(curation_tools.curation_suggest(file_path=arf2))["suggestion_id"]
        ep = submit_server.register("suggest", sid, arf2)
        status, body = _post(ep["port"], "/v1/submit", {"token": ep["token"], "review_id": sid, "flags": [
            {"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "sum", "note": ""},
            {"spot_id": 3, "flag": "redundant", "candidate": "R1", "note": ""}]})
        assert status == 200, body
        assert body["recorded"] == 2 and body["n_assign"] == 1 and body["n_redundant"] == 1
        assert body["tags_xml"]["added"] == [] and body["tags_xml"]["removed"] == []   # _tags.xml は不変
    finally:
        lib.close()
        session_state.session.__init__()
