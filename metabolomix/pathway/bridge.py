# metabolomix/pathway/bridge.py
"""MSLipidMapper（R パッケージ）を Rscript サブプロセスで呼ぶ。

契約は request.json → bridge-result.json（R 側は mslm_bridge.R）。Console と同じく
自前の Job Object で起動し、timeout ではプロセスツリーごと終了する。再試行しない。
session_state / mcp_core を import しない。
"""
from __future__ import annotations

import json
from pathlib import Path

from metabolomix.core import user_config
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.process_control import start_owned_process

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
SETUP_COMMAND = "Rscript r/setup.R"


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
    setting = user_config.get_setting("mslipidmapper.browser")
    return setting.value if setting else None


def _log_tail(log_path: Path) -> str:
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-_LOG_TAIL_LINES:])


def _environment_restored() -> bool:
    return any(next(RENV_PROJECT.glob(pattern), None) is not None
               for pattern in _RESTORED_MARKERS)


def run_bridge(request: dict, run_dir: Path, *, timeout_s: float = 300.0) -> dict:
    rscript = resolve_rscript()
    if not _environment_restored():
        raise BridgeError("MSLM_R_PACKAGES_MISSING",
                          f"R の仮想環境が未復元です（MSLipidMapper が r/renv/library にありません）。"
                          f"リポジトリ直下で `{SETUP_COMMAND}` を一度実行してください。",
                          {"setup": SETUP_COMMAND, "renv_project": str(RENV_PROJECT)})
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
