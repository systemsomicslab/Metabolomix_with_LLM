"""MSLipidMapper の REST API（研究室サーバ）でジョブを回し、成果物を run ディレクトリへ落とす。

API は上流の R/rest_api.R（`GET /healthz` → `POST /v1/jobs`（input と config の multipart）→
`GET /v1/jobs/<id>` を succeeded / failed まで → `GET /v1/jobs/<id>/artifacts/<name>`）。
config は YAML だが、JSON は YAML の部分集合なので JSON で送る（YAML の依存を持ち込まない）。
上流のサーバはジョブごとに 3 つのパスウェイを全部計算する（選べない）。

ここはどこにフォールバックするかを決めない（bridge.run_bridge が決める）。失敗はすべて
RemoteError で上げる。再試行しない。サーバ側のジョブは取り消せない（API が無い）。
session_state / mcp_core を import しない。
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import httpx

#: 生存確認と接続の待ち時間。届かないサーバで待たされず、すぐローカルの R へ回すため短く。
CONNECT_TIMEOUT_S = 5.0
#: 1 回の HTTP 往復（アップロード・成果物の取得を含む）の待ち時間。
REQUEST_TIMEOUT_S = 60.0
POLL_INTERVAL_S = 2.0
#: run ディレクトリに残す送受信の記録（bridge.log に相当する出どころ）。
JOB_RECORD = "remote-job.json"

#: テストが httpx.MockTransport を差し込む口。
_transport: httpx.BaseTransport | None = None
_sleep = time.sleep
_monotonic = time.monotonic


class RemoteError(Exception):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details or {}


def _error_text(response: httpx.Response) -> str:
    try:
        body = response.json()
    except ValueError:
        return response.text[:200]
    return str(body.get("error") or body) if isinstance(body, dict) else str(body)[:200]


def _client(base_url: str, token: str | None) -> httpx.Client:
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    return httpx.Client(base_url=base_url.rstrip("/"), headers=headers, transport=_transport,
                        timeout=httpx.Timeout(REQUEST_TIMEOUT_S, connect=CONNECT_TIMEOUT_S))


def _check_health(client: httpx.Client, details: dict) -> None:
    try:
        response = client.get("/healthz", timeout=CONNECT_TIMEOUT_S)
    except httpx.HTTPError as exc:
        raise RemoteError("MSLM_REMOTE_UNAVAILABLE",
                          f"MSLipidMapper サーバに接続できません: {type(exc).__name__}",
                          details) from exc
    try:
        ok = response.status_code == 200 and response.json().get("ok") is True
    except (ValueError, AttributeError):
        ok = False
    if not ok:
        raise RemoteError("MSLM_REMOTE_UNAVAILABLE",
                          f"MSLipidMapper サーバの生存確認が通りません（HTTP {response.status_code}）。",
                          {**details, "http_status": response.status_code})


def _submit(client: httpx.Client, mztab: Path, config: dict, details: dict) -> str:
    files = {"input": (mztab.name, mztab.read_bytes(), "text/plain"),
             "config": ("config.yml", json.dumps(config, ensure_ascii=False).encode("utf-8"),
                        "application/x-yaml")}
    try:
        response = client.post("/v1/jobs", files=files)
    except httpx.HTTPError as exc:
        raise RemoteError("MSLM_REMOTE_UNAVAILABLE",
                          f"ジョブを投入できません: {type(exc).__name__}", details) from exc
    if response.status_code != 202:
        # 404 / 5xx はサーバ（API）側の不在・不調、それ以外の 4xx は要求の拒否
        code = ("MSLM_REMOTE_UNAVAILABLE" if response.status_code == 404 or response.status_code >= 500
                else "MSLM_RUN_FAILED")
        raise RemoteError(code,
                          f"ジョブの投入が拒否されました（HTTP {response.status_code}）: "
                          f"{_error_text(response)}",
                          {**details, "http_status": response.status_code})
    try:
        return str(response.json()["job_id"])
    except (ValueError, KeyError, TypeError) as exc:
        raise RemoteError("MSLM_OUTPUT_INVALID", "ジョブ投入の応答に job_id がありません。",
                          details) from exc


def _wait(client: httpx.Client, job_id: str, timeout_s: float, details: dict) -> dict:
    deadline = _monotonic() + timeout_s
    while True:
        try:
            response = client.get(f"/v1/jobs/{job_id}")
            status = response.json() if response.status_code == 200 else None
        except (httpx.HTTPError, ValueError) as exc:
            raise RemoteError("MSLM_REMOTE_UNAVAILABLE",
                              f"ジョブの状態を取得できません: {type(exc).__name__}", details) from exc
        if not isinstance(status, dict):
            raise RemoteError("MSLM_REMOTE_UNAVAILABLE",
                              f"ジョブの状態を取得できません（HTTP {response.status_code}）: "
                              f"{_error_text(response)}",
                              {**details, "http_status": response.status_code})
        if status.get("status") == "succeeded":
            return status
        if status.get("status") == "failed":
            raise RemoteError("MSLM_RUN_FAILED",
                              f"サーバ側の MSLipidMapper が失敗しました: {status.get('error') or ''}"
                              f"（詳細はサーバの worker.log）",
                              {**details, "remote_status": status})
        if _monotonic() >= deadline:
            raise RemoteError("MSLM_TIMEOUT",
                              f"サーバのジョブが {timeout_s:.0f} 秒以内に終わりませんでした"
                              f"（サーバ側のジョブは取り消していません）。",
                              {**details, "timeout_s": timeout_s})
        _sleep(POLL_INTERVAL_S)


def _download(client: httpx.Client, job_id: str, name: str, dest: Path, details: dict) -> None:
    try:
        response = client.get(f"/v1/jobs/{job_id}/artifacts/{name}")
    except httpx.HTTPError as exc:
        raise RemoteError("MSLM_REMOTE_UNAVAILABLE",
                          f"成果物 {name} を取得できません: {type(exc).__name__}", details) from exc
    if response.status_code != 200:
        raise RemoteError("MSLM_OUTPUT_INVALID",
                          f"成果物 {name} を取得できません（HTTP {response.status_code}）。",
                          {**details, "artifact": name, "http_status": response.status_code})
    dest.write_bytes(response.content)


def run_job(base_url: str, token: str | None, mztab: Path, config: dict,
            artifacts: list[str], run_dir: Path, *, timeout_s: float) -> dict:
    """ジョブを 1 つ回し、`artifacts`（成果物名）を run_dir に同名で保存する。

    戻り値はサーバの状態応答の `result`（`job_id` を足したもの）。
    送受信の記録は run_dir/remote-job.json（トークンは書かない）。
    """
    run_dir = Path(run_dir)
    details = {"api_url": base_url, "record": str(run_dir / JOB_RECORD)}
    record: dict = {"api_url": base_url, "config": config, "artifacts": artifacts}

    def save_record():
        (run_dir / JOB_RECORD).write_text(json.dumps(record, ensure_ascii=False, indent=1),
                                          encoding="utf-8")

    try:
        with _client(base_url, token) as client:
            _check_health(client, details)
            job_id = _submit(client, Path(mztab), config, details)
            record["job_id"] = details["job_id"] = job_id
            status = _wait(client, job_id, timeout_s, details)
            record["status"] = status
            result = status.get("result") if isinstance(status.get("result"), dict) else {}
            produced = result.get("artifacts")
            absent = [a for a in artifacts if isinstance(produced, list) and a not in produced]
            if absent:
                raise RemoteError("MSLM_OUTPUT_INVALID",
                                  f"サーバのジョブに成果物がありません: {', '.join(absent)}",
                                  {**details, "artifacts": produced})
            for name in artifacts:
                _download(client, job_id, name, run_dir / name, details)
    except RemoteError as exc:
        record["error"] = {"code": exc.code, "message": exc.message}
        save_record()
        raise
    save_record()
    return {**result, "job_id": job_id}
