# tests/test_pathway_remote.py
"""研究室サーバ（MSLipidMapper REST API）経由の実行と、ローカルの R へのフォールバック。

サーバは httpx.MockTransport の偽物。ローカルの R は bridge._run_local を差し替える
（実物の Rscript の契約は test_pathway_bridge.py が縛る）。
"""
import json

import httpx
import pytest

from metabolomix.pathway import bridge, remote, tools

from tests.test_pathway_jsonld import _graph
from tests.test_pathway_tools import lipid_arf  # noqa: F401  （fixture）

URL = "http://lab.example:8000"


class FakeServer:
    """上流 R/rest_api.R の振る舞いを模す最小の偽サーバ。"""

    def __init__(self, *, health=200, submit=202, final="succeeded", polls_before_done=1,
                 result=None, artifacts=None, missing_artifacts=()):
        self.health, self.submit, self.final = health, submit, final
        self.polls_before_done = polls_before_done
        self.result = result or {}
        self.artifacts = artifacts or {}
        self.missing_artifacts = set(missing_artifacts)
        self.requests: list[httpx.Request] = []
        self.config = None
        self.input_name = None

    def __call__(self, request: httpx.Request) -> httpx.Response:
        self.requests.append(request)
        path = request.url.path
        if path == "/healthz":
            return httpx.Response(self.health, json={"ok": True} if self.health == 200 else {"error": "x"})
        if path == "/v1/jobs" and request.method == "POST":
            body = request.read()
            # multipart の config パートから JSON（＝YAML）を取り出す
            part = body.split(b'name="config"', 1)[1].split(b"\r\n\r\n", 1)[1]
            self.config = json.loads(part.split(b"\r\n--", 1)[0])
            self.input_name = body.split(b'name="input"; filename="', 1)[1].split(b'"', 1)[0].decode()
            if self.submit != 202:
                return httpx.Response(self.submit, json={"error": "Unknown config section."})
            return httpx.Response(202, json={"job_id": "job_1", "status": "queued"})
        if path == "/v1/jobs/job_1":
            if self.polls_before_done > 0:
                self.polls_before_done -= 1
                return httpx.Response(200, json={"job_id": "job_1", "status": "running"})
            if self.final == "failed":
                return httpx.Response(200, json={"job_id": "job_1", "status": "failed",
                                                 "error": "Analysis failed; consult the server worker.log."})
            if self.final == "running":
                return httpx.Response(200, json={"job_id": "job_1", "status": "running"})
            produced = [n for n in self.artifacts if n not in self.missing_artifacts]
            return httpx.Response(200, json={"job_id": "job_1", "status": "succeeded",
                                             "result": {**self.result, "artifacts": produced}})
        if path.startswith("/v1/jobs/job_1/artifacts/"):
            name = path.rsplit("/", 1)[1]
            if name not in self.artifacts:
                return httpx.Response(404, json={"error": "Artifact not found."})
            return httpx.Response(200, content=self.artifacts[name])
        return httpx.Response(404, json={"error": "404 - Resource Not Found"})


@pytest.fixture
def no_wait(monkeypatch):
    monkeypatch.setattr(remote, "_sleep", lambda s: None)


@pytest.fixture
def serve(monkeypatch, no_wait):
    def install(server):
        monkeypatch.setenv("MSLIPIDMAPPER_API_URL", URL + "/")
        monkeypatch.setattr(remote, "_transport", httpx.MockTransport(server))
        return server
    return install


@pytest.fixture
def local_calls(monkeypatch):
    """ローカルの R を差し替える。calls に request を積み、成功を返す。"""
    calls = []

    def run(request, run_dir, *, timeout_s=300.0):
        calls.append(request)
        return {"status": "ok", "outputs": {"remodeling": "local.json"}, "mapped_nodes": {},
                "package_version": "0.9", "package_sha": "abc"}
    monkeypatch.setattr(bridge, "_run_local", run)
    return calls


def _measure_request(tmp_path, chains=()):
    mztab = tmp_path / "input.mztab"
    mztab.write_text("MTD\tmzTab-version\t2.0.0-M\n", encoding="utf-8")
    return {"mode": "measure", "mztab": str(mztab), "out_dir": str(tmp_path),
            "pathways": ["remodeling", "global"],
            "networks": {p: bridge.PATHWAY_NETWORKS[p] for p in ("remodeling", "global")},
            "acyl_chains": list(chains), "aggregation": "mean",
            "normalization": "normalize=pqn", "unit": "arbitrary", "dataset_id": "pw_x"}


def _render_request(tmp_path):
    req = _measure_request(tmp_path)
    req.update(mode="render", pathways=["ceramide"], networks={"ceramide": "x"},
               browser=None, plot_type="box", x_order=["P", "C"], timeout_seconds=90)
    return req


def _measure_server(**kw):
    return FakeServer(
        result={"measurements": {"pathways": {"remodeling": {"measured_nodes": 3},
                                              "ceramide": {"measured_nodes": 0},
                                              "global": {"measured_nodes": 5}}}},
        artifacts={f"{p}-measurements.json": f'{{"p":"{p}"}}'.encode()
                   for p in ("remodeling", "ceramide", "global")}, **kw)


def test_measure_runs_on_the_server(serve, local_calls, tmp_path):
    server = serve(_measure_server())
    res = bridge.run_bridge(_measure_request(tmp_path, chains=["18:1"]), tmp_path)

    assert local_calls == []
    assert res["backend"] == {"kind": "remote", "api_url": URL, "job_id": "job_1"}
    assert "remote_fallback" not in res
    assert res["outputs"] == {p: str(tmp_path / f"{p}-measurements.json")
                              for p in ("remodeling", "global")}
    assert json.loads((tmp_path / "global-measurements.json").read_text()) == {"p": "global"}
    assert not (tmp_path / "ceramide-measurements.json").exists()  # 頼んでいない分は落とさない
    assert res["mapped_nodes"] == {"remodeling": 3, "global": 5}
    assert server.input_name == "input.mztab"
    # 前処理の記述文字列は API が受け付けないので none、鎖は include_subsets と組で送る
    assert server.config == {
        "input": {"format": "mztab-m"},
        "analysis": {"normalization": "none", "aggregation": "mean", "group_column": "class"},
        "measurements": {"unit": "arbitrary", "include_subsets": True, "acyl_chains": ["18:1"]},
        "output": {"format": ["jsonld"]}}
    record = json.loads((tmp_path / remote.JOB_RECORD).read_text(encoding="utf-8"))
    assert record["job_id"] == "job_1" and record["status"]["status"] == "succeeded"


def test_measure_without_chains_does_not_ask_for_subsets(serve, local_calls, tmp_path):
    server = serve(_measure_server())
    bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert server.config["measurements"] == {"unit": "arbitrary"}


def test_render_runs_on_the_server(serve, local_calls, tmp_path):
    server = serve(FakeServer(result={"mapped_nodes": {"ceramide": 4}},
                              artifacts={"ceramide.pdf": b"%PDF-1.4"}))
    res = bridge.run_bridge(_render_request(tmp_path), tmp_path)
    assert res["outputs"] == {"ceramide": str(tmp_path / "ceramide.pdf")}
    assert (tmp_path / "ceramide.pdf").read_bytes() == b"%PDF-1.4"
    assert res["mapped_nodes"] == {"ceramide": 4}
    assert server.config["plot"] == {"type": "box", "x_order": ["P", "C"]}
    assert server.config["pathway"] == {"timeout_seconds": 90}
    assert server.config["output"]["format"] == ["pdf"]


def test_token_is_sent_but_not_recorded(serve, local_calls, tmp_path, monkeypatch):
    monkeypatch.setenv("MSLIPIDMAPPER_API_TOKEN", "s3cret")
    server = serve(_measure_server())
    bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert all(r.headers["Authorization"] == "Bearer s3cret" for r in server.requests)
    assert "s3cret" not in (tmp_path / remote.JOB_RECORD).read_text(encoding="utf-8")


def test_without_api_url_the_server_is_never_contacted(local_calls, tmp_path, monkeypatch):
    def refuse(request):
        raise AssertionError("サーバに触れてはいけない")
    monkeypatch.setattr(remote, "_transport", httpx.MockTransport(refuse))
    res = bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert len(local_calls) == 1
    assert res["backend"] == {"kind": "local_r"} and "remote_fallback" not in res


@pytest.mark.parametrize("server,code", [
    (FakeServer(health=404), "MSLM_REMOTE_UNAVAILABLE"),
    (FakeServer(health=503), "MSLM_REMOTE_UNAVAILABLE"),
    (FakeServer(submit=404), "MSLM_REMOTE_UNAVAILABLE"),
    (FakeServer(submit=400), "MSLM_RUN_FAILED"),
    (FakeServer(final="failed"), "MSLM_RUN_FAILED"),
    (_measure_server(missing_artifacts={"global-measurements.json"}), "MSLM_OUTPUT_INVALID"),
])
def test_server_failure_falls_back_to_local_r(serve, local_calls, tmp_path, server, code):
    serve(server)
    res = bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert len(local_calls) == 1
    assert res["backend"] == {"kind": "local_r"}
    assert res["remote_fallback"]["code"] == code
    assert res["remote_fallback"]["api_url"] == URL
    assert res["package_version"] == "0.9"


def test_unreachable_server_falls_back_to_local_r(local_calls, tmp_path, monkeypatch, no_wait):
    def down(request):
        raise httpx.ConnectError("refused", request=request)
    monkeypatch.setenv("MSLIPIDMAPPER_API_URL", URL)
    monkeypatch.setattr(remote, "_transport", httpx.MockTransport(down))
    res = bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert len(local_calls) == 1
    assert res["remote_fallback"]["code"] == "MSLM_REMOTE_UNAVAILABLE"


def test_server_timeout_falls_back(serve, local_calls, tmp_path, monkeypatch):
    clock = iter(range(0, 10_000, 100))
    monkeypatch.setattr(remote, "_monotonic", lambda: next(clock))
    serve(FakeServer(final="running"))
    res = bridge.run_bridge(_measure_request(tmp_path), tmp_path, timeout_s=250)
    assert res["remote_fallback"]["code"] == "MSLM_TIMEOUT"


def test_server_down_and_no_local_r_raises_remote_error_with_local_reason(serve, tmp_path):
    serve(FakeServer(health=404))
    with pytest.raises(bridge.BridgeError) as exc:  # MSLIPIDMAPPER_RSCRIPT は conftest が消している
        bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert exc.value.code == "MSLM_REMOTE_UNAVAILABLE"
    assert exc.value.details["remote"]["api_url"] == URL
    assert exc.value.details["local"]["code"] == "MSLM_NOT_CONFIGURED"
    assert "ローカルの R も使えません" in exc.value.message


def test_failed_render_job_without_local_r_is_render_failed(serve, tmp_path):
    serve(FakeServer(final="failed"))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge(_render_request(tmp_path), tmp_path)
    assert exc.value.code == "PATHWAY_RENDER_FAILED"


def test_local_failure_after_fallback_keeps_both_reasons(serve, tmp_path, monkeypatch):
    serve(FakeServer(health=404))

    def broken(request, run_dir, *, timeout_s=300.0):
        raise bridge.BridgeError("MSLM_RUN_FAILED", "R が落ちた", {"log": "bridge.log"})
    monkeypatch.setattr(bridge, "_run_local", broken)
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge(_measure_request(tmp_path), tmp_path)
    assert exc.value.code == "MSLM_RUN_FAILED"
    assert exc.value.details["log"] == "bridge.log"
    assert exc.value.details["remote_fallback"]["code"] == "MSLM_REMOTE_UNAVAILABLE"


def test_pathway_project_through_the_server(lipid_arf, serve, tmp_path):  # noqa: F811
    doc = _graph()
    for x in doc["@graph"]:
        if x["@type"] == "mslm:Sample":
            x["sampleGroup"] = "C" if x["identifier"] == "S1" else "P"
    body = json.dumps(doc).encode()
    serve(FakeServer(result={"measurements": {"pathways": {"remodeling": {"measured_nodes": 2}}}},
                     artifacts={f"{p}-measurements.json": body
                                for p in ("remodeling", "ceramide", "global")}))
    out = json.loads(tools.pathway_project(source="arf", group_a="C", group_b="P",
                                           pathways=["remodeling"]))
    assert "error" not in out, out
    assert out["backend"] == f"研究室サーバ {URL}（job job_1）"
    manifest = json.loads((tmp_path / "reports" / "pathway" / out["run_id"] / "manifest.json")
                          .read_text(encoding="utf-8"))
    assert manifest["backend"] == {"kind": "remote", "api_url": URL, "job_id": "job_1"}
    assert manifest["remote_fallback"] is None
