import http.client
import json

import pytest

from metabolomix.curation import submit_server


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
