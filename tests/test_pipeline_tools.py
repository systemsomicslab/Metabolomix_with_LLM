"""pipeline系5件のMCPツール（Task18: 薄いMCP層。実体はmetabolomix.pipeline.service/
metabolomix.pipeline.recovery）。

ここでは各ツールが (1) 対応するservice/recovery関数へ引数をそのまま渡し、
(2) 戻り値をjson_payloadで返し、(3) DomainErrorをconsole_errorへ変換し、
(4) 戻り値が行列・スコア・loadings・volcano点列を一切含まないcompactな
発送receiptであることだけを検証する。service/recovery自体の振る舞い
（受付冪等性・handler配線・resume判定等）はtest_pipeline_engine.py /
test_pipeline_recovery.py / test_pipeline_report.py 等が別に検証する。
"""
from __future__ import annotations

import json

from metabolomix.core.atomic_io import DomainError


def test_run_returns_only_compact_dispatch_receipt(monkeypatch):
    monkeypatch.setattr("metabolomix.pipeline.service.start_pipeline",
        lambda *a, **kw: {"status": "running", "pipeline_id": "p1",
                         "pipeline_path": "C:/fake/pipeline-run.json"})
    from metabolomix.tools.pipeline_tools import pipeline_run
    payload = json.loads(pipeline_run("C:/fake/source"))
    assert payload["pipeline_id"] == "p1"
    assert not {"matrix", "scores", "loadings", "volcano"} & payload.keys()


def test_run_forwards_dataset_root_as_path_and_request_args(monkeypatch):
    """`dataset_root: str` は `Path` へ変換して渡す。request/request_idはそのまま。"""
    from pathlib import Path

    captured = {}

    def fake_start(dataset_root, request, request_id):
        captured["dataset_root"] = dataset_root
        captured["request"] = request
        captured["request_id"] = request_id
        return {"status": "planned", "pipeline_id": "p2", "pipeline_path": "C:/fake/p"}

    monkeypatch.setattr("metabolomix.pipeline.service.start_pipeline", fake_start)
    from metabolomix.tools.pipeline_tools import pipeline_run
    pipeline_run("C:/fake/source", {"polarity": "negative"}, "req-1")

    assert captured["dataset_root"] == Path("C:/fake/source")
    assert captured["request"] == {"polarity": "negative"}
    assert captured["request_id"] == "req-1"


def test_run_converts_domain_error_to_console_error_envelope(monkeypatch):
    def raise_domain_error(*a, **kw):
        raise DomainError("MSDIAL_EXE_NOT_FOUND", "execが見つかりません", {"hint": "x"})

    monkeypatch.setattr("metabolomix.pipeline.service.start_pipeline", raise_domain_error)
    from metabolomix.tools.pipeline_tools import pipeline_run
    payload = json.loads(pipeline_run("C:/fake/source"))
    assert payload["error"]["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert payload["error"]["details"] == {"hint": "x"}


def test_plan_forwards_to_plan_pipeline_and_does_not_launch(monkeypatch):
    """`pipeline_plan` は `plan_pipeline`（起動しない受付専用関数）だけを呼ぶ。

    `start_pipeline`を誤って呼んでいれば、こちらをpatchしていないため
    実サービス層が動いてしまい、フォルダ不在などで例外になって落ちる。
    """
    from pathlib import Path

    monkeypatch.setattr("metabolomix.pipeline.service.plan_pipeline",
        lambda dataset_root, request, request_id: {
            "status": "planned", "pipeline_id": "p3",
            "pipeline_path": str(Path(dataset_root) / "pipeline-run.json"),
            "launched": False})
    from metabolomix.tools.pipeline_tools import pipeline_plan
    payload = json.loads(pipeline_plan("C:/fake/source2"))
    assert payload["launched"] is False
    assert payload["pipeline_id"] == "p3"


def test_status_forwards_include_details_and_is_read_only_receipt(monkeypatch):
    from pathlib import Path

    captured = {}

    def fake_read_status(pipeline_path, *, include_details):
        captured["pipeline_path"] = pipeline_path
        captured["include_details"] = include_details
        return {"status": "running", "pipeline_id": "p4", "observed_health": "ok"}

    monkeypatch.setattr("metabolomix.pipeline.recovery.read_status", fake_read_status)
    from metabolomix.tools.pipeline_tools import pipeline_status
    payload = json.loads(pipeline_status("C:/fake/pipeline-run.json", include_details=True))

    assert captured["pipeline_path"] == Path("C:/fake/pipeline-run.json")
    assert captured["include_details"] is True
    assert payload["status"] == "running"


def test_resume_forwards_updates_request_id_and_rerun_upstream(monkeypatch):
    captured = {}

    def fake_resume(pipeline_path, updates, request_id, rerun_upstream):
        captured["updates"] = updates
        captured["request_id"] = request_id
        captured["rerun_upstream"] = rerun_upstream
        return {"status": "planned", "pipeline_path": str(pipeline_path)}

    monkeypatch.setattr("metabolomix.pipeline.service.resume_pipeline", fake_resume)
    from metabolomix.tools.pipeline_tools import pipeline_resume
    pipeline_resume("C:/fake/p", updates={"target": "differential"},
                    request_id="req-2", rerun_upstream=True)

    assert captured["updates"] == {"target": "differential"}
    assert captured["request_id"] == "req-2"
    assert captured["rerun_upstream"] is True


def test_cancel_forwards_pipeline_path_and_is_idempotent_by_annotation(monkeypatch):
    from pathlib import Path

    captured = {}

    def fake_request_cancel(pipeline_path):
        captured["pipeline_path"] = pipeline_path
        return {"pipeline_id": "p5", "accepted": True, "status": "cancelled"}

    monkeypatch.setattr("metabolomix.pipeline.recovery.request_cancel", fake_request_cancel)
    from metabolomix.tools.pipeline_tools import pipeline_cancel
    payload = json.loads(pipeline_cancel("C:/fake/p"))

    assert captured["pipeline_path"] == Path("C:/fake/p")
    assert payload["accepted"] is True


def test_pipeline_status_is_read_only_others_are_not():
    """brief拘束: pipeline_statusだけreadOnlyHint=true、他4件はfalse。

    `@mcp.tool`はannotationsを関数自体には残さずFastMCPのTool登録側へ持つため、
    `tests/test_tool_annotations.py`と同じ経路（`server.mcp.list_tools()`）で読む。
    """
    import asyncio

    import server

    tools = {t.name: t.annotations for t in asyncio.run(server.mcp.list_tools())}
    assert tools["pipeline_status"].readOnlyHint is True
    for name in ("pipeline_plan", "pipeline_run", "pipeline_resume", "pipeline_cancel"):
        assert tools[name].readOnlyHint is False, name
