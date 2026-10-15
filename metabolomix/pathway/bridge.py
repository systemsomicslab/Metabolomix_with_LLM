# metabolomix/pathway/bridge.py
"""MSLipidMapper を呼ぶ。研究室サーバの REST API が設定されていればそちらを先に使い、
使えなければローカルの R（Rscript サブプロセス）へ回す。

入口は run_bridge だけ（request の dict → bridge-result.json と同じ形の dict）。
ローカルの契約は request.json → bridge-result.json（R 側は mslm_bridge.R）。Console と同じく
自前の Job Object で起動し、timeout ではプロセスツリーごと終了する。リモートは remote.py。
同じ経路での再試行はしない。session_state / mcp_core を import しない。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from metabolomix.core import user_config
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.process_control import start_owned_process
from metabolomix.pathway import remote

BRIDGE_SCRIPT = Path(__file__).with_name("mslm_bridge.R")
#: R の仮想環境（renv）。パッケージはここにだけ入る（spec §2）。
RENV_PROJECT = user_config.REPO_ROOT / "r"
PATHWAY_NETWORKS = {"remodeling": "remodeling_lipidonly.cyjs",
                    "ceramide": "ceramidepathway_lipidonly.cyjs",
                    "global": "global.cyjs"}
_LOG_TAIL_LINES = 20
#: 復元済みライブラリの目印。renv の既定配置は library/<platform>/R-x.y/<arch>/<pkg>、
#: 古い renv は library/R-x.y/<arch>/<pkg>。activate.R は追跡物なので目印にならない。
_RESTORED_MARKERS = ("renv/library/*/*/*/MSLipidMapper/DESCRIPTION",
                     "renv/library/*/*/MSLipidMapper/DESCRIPTION")
#: `r/setup.R` が復元に使った renv.lock の md5（RENV_PROJECT 相対。library/ ごと git の追跡外）。
#: 今の renv.lock と食い違えば、更新で lock が変わったのに再復元していない＝古い版のまま。
RESTORE_STAMP = "renv/library/restored-lock.md5"
SETUP_COMMAND = "Rscript r/setup.R"
#: 描画に使うヘッドレス Chromium の標準の置き場所（設定が無いときに順に探す）。
#: Edge は通常 Program Files (x86) にあり、上流の自動検出はそこを見ない。
BROWSER_CANDIDATES = (
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
)


class BridgeError(Exception):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details or {}


def resolve_rscript() -> str:
    key = "mslipidmapper.rscript"
    setting = user_config.get_setting(key)
    if setting is None:
        raise BridgeError("MSLM_NOT_CONFIGURED",
                          user_config.missing_hint(key, "Rscript（MSLipidMapper 入りの R）"),
                          user_config.describe_missing(key))
    if not Path(setting.value).is_file():
        raise BridgeError("MSLM_NOT_CONFIGURED",
                          f"{user_config.setting_label(setting)} が指す Rscript がありません。",
                          user_config.describe_missing(key, setting))
    return setting.value


def browser_path() -> str | None:
    """設定（環境変数 → 設定ファイル）を優先し、無ければ標準の置き場所から探す。

    どれも無ければ None（上流の自動検出に任せる。見つからなければ描画が
    PATHWAY_RENDER_FAILED になる）。
    """
    setting = user_config.get_setting("mslipidmapper.browser")
    if setting:
        return setting.value
    return next((c for c in BROWSER_CANDIDATES if Path(c).is_file()), None)


def _log_tail(log_path: Path) -> str:
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-_LOG_TAIL_LINES:])


def _environment_restored() -> bool:
    return any(next(RENV_PROJECT.glob(pattern), None) is not None
               for pattern in _RESTORED_MARKERS)


def _environment_matches_lock() -> bool:
    """setup.R の目印が今の renv.lock と一致するか（目印が無い・lock が読めないときは False）。"""
    try:
        stamp = (RENV_PROJECT / RESTORE_STAMP).read_text(encoding="utf-8").strip().lower()
        current = hashlib.md5((RENV_PROJECT / "renv.lock").read_bytes()).hexdigest()
    except OSError:
        return False
    return stamp == current


def _check_environment() -> None:
    details = {"setup": SETUP_COMMAND, "renv_project": str(RENV_PROJECT)}
    if not _environment_restored():
        raise BridgeError("MSLM_R_PACKAGES_MISSING",
                          f"R の仮想環境が未復元です（MSLipidMapper が r/renv/library にありません）。"
                          f"リポジトリ直下で `{SETUP_COMMAND}` を一度実行してください。",
                          {**details, "reason": "not_restored"})
    if not _environment_matches_lock():
        raise BridgeError("MSLM_R_PACKAGES_MISSING",
                          f"R の仮想環境が r/renv.lock と一致しません（更新で lock が変わった、"
                          f"または目印の導入前に復元した環境です）。リポジトリ直下で "
                          f"`{SETUP_COMMAND}` を再実行してください。",
                          {**details, "reason": "stale"})


#: ローカルの R が「使えない」ことを表すコード（「動かして失敗した」とは区別する）。
_LOCAL_UNAVAILABLE = frozenset({"MSLM_NOT_CONFIGURED", "MSLM_R_PACKAGES_MISSING"})


def api_url() -> str | None:
    """[mslipidmapper] api_url（環境変数 MSLIPIDMAPPER_API_URL）。未設定なら None。"""
    setting = user_config.get_setting("mslipidmapper.api_url")
    return setting.value.rstrip("/") if setting else None


def _api_token() -> str | None:
    setting = user_config.get_setting("mslipidmapper.api_token")
    return setting.value if setting else None


def _remote_config(request: dict) -> tuple[dict, dict[str, str]]:
    """request（ローカルと同じ形）→ (API の config, パスウェイ → 成果物名)。

    値はこちらの前処理後行列なので analysis.normalization は none（API は前処理の
    記述文字列を受け付けない。正規化の出どころは manifest.json の normalization）。
    """
    base = {"input": {"format": "mztab-m"},
            "analysis": {"normalization": "none", "aggregation": request["aggregation"],
                         "group_column": "class"}}
    if request["mode"] == "measure":
        chains = list(request.get("acyl_chains") or [])
        measurements = {"unit": request.get("unit") or "arbitrary"}
        if chains:  # 上流は acyl_chains に include_subsets: true を要求する
            measurements.update(include_subsets=True, acyl_chains=chains)
        config = {**base, "measurements": measurements, "output": {"format": ["jsonld"]}}
        names = {p: f"{p}-measurements.json" for p in request["pathways"]}
    elif request["mode"] == "render":
        config = {**base, "pathway": {"timeout_seconds": request["timeout_seconds"]},
                  "plot": {"type": request["plot_type"], "x_order": list(request["x_order"])},
                  "output": {"format": ["pdf"], "paper_size": "A4", "orientation": "landscape"}}
        names = {p: f"{p}.pdf" for p in request["pathways"]}
    else:
        raise BridgeError("PATHWAY_BAD_REQUEST", f"未知の mode: {request['mode']}")
    return config, names


def _mapped_nodes(request: dict, result: dict) -> dict:
    if request["mode"] == "measure":
        per = ((result.get("measurements") or {}).get("pathways")) or {}
        return {p: (per.get(p) or {}).get("measured_nodes") for p in request["pathways"]}
    per = result.get("mapped_nodes") or {}
    return {p: per.get(p) for p in request["pathways"]}


def _run_remote(base_url: str, request: dict, run_dir: Path, timeout_s: float) -> dict:
    config, names = _remote_config(request)
    try:
        result = remote.run_job(base_url, _api_token(), Path(request["mztab"]), config,
                                list(names.values()), Path(run_dir), timeout_s=timeout_s)
    except remote.RemoteError as exc:
        # 描画ジョブの失敗はローカルと同じく描画だけの失敗（pathway_project の結果は有効）
        code = ("PATHWAY_RENDER_FAILED" if request["mode"] == "render" and exc.code == "MSLM_RUN_FAILED"
                else exc.code)
        raise BridgeError(code, exc.message, exc.details) from exc
    # サーバは版を返さないので package_version / package_sha は不明（None）
    return {"status": "ok", "code": None, "message": "",
            "package_version": None, "package_sha": None,
            "outputs": {p: str(Path(run_dir) / n) for p, n in names.items()},
            "mapped_nodes": _mapped_nodes(request, result),
            "backend": {"kind": "remote", "api_url": base_url, "job_id": result.get("job_id")}}


def run_bridge(request: dict, run_dir: Path, *, timeout_s: float = 300.0) -> dict:
    """MSLipidMapper を 1 回呼ぶ。戻り値の `backend` がどちらで動いたかを示す。

    api_url が設定されていればまずサーバへ。サーバが使えない・失敗したとき
    （未到達・404・ジョブ失敗・時間切れ・成果物の欠落）は、ローカルの R へ回し、
    その理由を `remote_fallback` に残す。ローカルの R も使えなければ（未設定・
    仮想環境が未復元）、サーバ側のエラーを上げる（details.local にローカルの理由）。
    """
    url = api_url()
    if url is None:
        return {**_run_local(request, run_dir, timeout_s=timeout_s), "backend": {"kind": "local_r"}}
    try:
        return _run_remote(url, request, run_dir, timeout_s)
    except BridgeError as remote_exc:
        fallback = {"api_url": url, "code": remote_exc.code, "message": remote_exc.message}
    try:
        result = _run_local(request, run_dir, timeout_s=timeout_s)
    except BridgeError as local_exc:
        if local_exc.code in _LOCAL_UNAVAILABLE:
            raise BridgeError(fallback["code"],
                              f"{fallback['message']} ローカルの R も使えません: {local_exc.message}",
                              {"remote": fallback,
                               "local": {"code": local_exc.code, "message": local_exc.message,
                                         **local_exc.details}}) from local_exc
        raise BridgeError(local_exc.code, local_exc.message,
                          {**local_exc.details, "remote_fallback": fallback}) from local_exc
    return {**result, "backend": {"kind": "local_r"}, "remote_fallback": fallback}


def _run_local(request: dict, run_dir: Path, *, timeout_s: float = 300.0) -> dict:
    rscript = resolve_rscript()
    _check_environment()
    run_dir = Path(run_dir)
    request_path = run_dir / "request.json"
    result_path = run_dir / "bridge-result.json"
    log_path = run_dir / "bridge.log"
    result_path.unlink(missing_ok=True)
    request_path.write_text(json.dumps(request, ensure_ascii=False), encoding="utf-8")

    try:
        proc = start_owned_process([rscript, str(BRIDGE_SCRIPT), str(RENV_PROJECT), str(request_path)],
                                   cwd=run_dir, log_path=log_path)
    except DomainError as exc:
        raise BridgeError("MSLM_RUN_FAILED", f"Rscript を起動できません: {exc.message}",
                          {"cause_code": exc.code, "cause_message": exc.message,
                           "log": str(log_path)}) from exc
    try:
        try:
            exit_code = proc.wait(timeout=timeout_s)
        except DomainError as exc:
            if exc.code != "PROCESS_WAIT_TIMEOUT":
                raise
            proc.terminate_tree()
            raise BridgeError("MSLM_TIMEOUT",
                              f"MSLipidMapper が {timeout_s:.0f} 秒以内に終わりませんでした。",
                              {"timeout_s": timeout_s, "log": str(log_path)}) from exc
    finally:
        proc.close()

    result = None
    if result_path.is_file():
        try:
            result = json.loads(result_path.read_text(encoding="utf-8"))
        except ValueError:
            result = None
    if result is not None and result.get("status") == "error":
        raise BridgeError(result.get("code") or "MSLM_RUN_FAILED",
                          result.get("message") or "MSLipidMapper が失敗しました。",
                          {"log_tail": _log_tail(log_path), "log": str(log_path)})
    if exit_code != 0:
        raise BridgeError("MSLM_RUN_FAILED",
                          f"Rscript が終了コード {exit_code} で終わりました。",
                          {"exit_code": exit_code, "log_tail": _log_tail(log_path),
                           "log": str(log_path)})
    if result is None:
        raise BridgeError("MSLM_OUTPUT_INVALID",
                          "Rscript は正常終了しましたが bridge-result.json がありません。",
                          {"log_tail": _log_tail(log_path), "log": str(log_path)})
    return result
