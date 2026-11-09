"""metabolomix.pipeline を実プロセスとして検証するための試験専用workerエントリ。

production worker（`metabolomix/pipeline/worker.py`）は`--pipeline`しか受け付けない
——試験用の環境変数も任意import名も持たせない（それ自体がコード実行の入口に
なる）。したがって「実プロセスとして起動されて初めて確認できる」振る舞いを
確かめるには、本番からは決して呼ばれない別のエントリが要る。それがこのファイル
（MCPからは到達できない試験専用の経路）。
**importせずスクリプトパスとして`python tests/pipeline_worker_harness.py`で起動する**
（`tests/fixtures/fake_console.py` と同じ流儀）。

2つのモードがある。

``--scenario`` 無し
    すべてのhandlerを「呼ばれたら成功する」合成版（`synthetic_handlers`）に差し替えて
    `run_engine`を回す。engine自体（stage順序・owner lock・冪等性）だけを見るためのモード。

``--scenario <name>``（Task 19）
    **本番の`metabolomix.pipeline.service.build_handlers()`をそのまま使う。**
    差し替えるのは`metabolomix.console.runner.build_msdial_cmd`だけ——起動する
    Consoleのコマンドラインを`tests/fixtures/fake_console.py`へ向ける。
    配置・mzTab-M と解析フォルダの検証・受領はすべて本物が走る。`supervise`
    （起動・監視・停止・収集・完了判定）も本物。

使い方::

    python tests/pipeline_worker_harness.py --pipeline <pipeline_root>
        [--sleep-stage <stage_id> --sleep-seconds <n>]
        [--scenario success|nonzero|hang|invalid|missing_sample|unmatched_arf [--no-inchikey]]
        [--counter <path>]
        [--launch-detached --launch-info <path>]

`--sleep-stage` はどちらのモードでも効く（合成版は自分で待ち、実handlerは
呼出し前に待つ薄いラッパを被せる）。owner lockが実際に保持され続けること・
工程の途中でworkerを失ったときの挙動・stage境界での取消を、実時間で確かめる
ための唯一の仕掛け。

`--launch-detached` は、このプロセス自身が**起動役**になり、同じ引数から
`--launch-detached`だけを取り除いたコマンドを`launch_detached`で切り離して
起こし、起動情報を`--launch-info`へ書いて即座に終了する。
「MCP相当の親プロセスが消えても解析が続く」を、起動役が本当に居なくなった
状態で観測するために使う。
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from metabolomix.core.atomic_io import DomainError  # noqa: E402
from metabolomix.pipeline.engine import run_engine  # noqa: E402

#: 偽Consoleスクリプト（importせずスクリプトパスとして渡す）。
_FAKE_CONSOLE = Path(__file__).resolve().parent / "fixtures" / "fake_console.py"

#: build_stagesが組み立てうるhandlerキーの全量（stage_idではなくhandlerキー単位）。
HANDLER_KEYS = ("prepare_input", "upstream", "validate_outputs", "receipt")


def synthetic_handlers(calls: list | None = None, *, fail_stages=(), unmatched_units=(),
                       sleep_stage: str | None = None, sleep_seconds: float = 0.0) -> dict:
    """呼ばれたら成功する合成handler一式（engine自体＝工程順・owner lock・単位の扱いを見る）。

    `fail_stages` に挙げた stage_id は failed を返す。`unmatched_units` の極性は
    validate_outputs が解析フォルダの照合に失敗した形（`alignment: None`）を返す。
    `receipt` は本物の `receipt.write_receipt` で受領を書く。
    """
    fail_stages = set(fail_stages)
    unmatched_units = set(unmatched_units)

    def make(name):
        def handler(context: dict) -> dict:
            stage_id = context["stage_id"]
            if calls is not None:
                calls.append(stage_id)
            if sleep_stage is not None and stage_id == sleep_stage and sleep_seconds > 0:
                marker = stage_marker_path(context["pipeline_root"], sleep_stage)
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.write_text("entered", encoding="ascii")
                time.sleep(sleep_seconds)
            if stage_id in fail_stages:
                return {"status": "failed", "result_refs": [], "warnings": [],
                        "error": {"code": "SYNTHETIC_FAILURE", "message": stage_id,
                                  "details": {}}}
            unit = context.get("unit")
            if name == "upstream":
                return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None,
                        "unit_updates": {"upstream": {
                            "console_job_path": None, "execution_id": f"exec-fake-{unit}",
                            "verification": {"status": "completed"}}}}
            if name == "validate_outputs":
                alignment = (None if unit in unmatched_units
                             else {"arf": f"{unit}_PeakProperties.arf", "arf2": f"{unit}.arf2"})
                return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None,
                        "unit_updates": {"outputs": {"mztab": None, "alignment": alignment,
                                                     "gui_project": None}}}
            if name == "receipt":
                from metabolomix.pipeline import receipt as receipt_mod
                from metabolomix.pipeline import store
                record = store.load_run(context["pipeline_root"])
                _receipt, ref = receipt_mod.write_receipt(record, context["pipeline_root"])
                return {"status": "succeeded", "result_refs": [ref], "warnings": [],
                        "error": None}
            return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None}
        return handler

    return {key: make(key) for key in HANDLER_KEYS}


def _build_handlers(sleep_stage: str | None, sleep_seconds: float) -> dict:
    return synthetic_handlers(sleep_stage=sleep_stage, sleep_seconds=sleep_seconds)


# ---------- Task 19: 実handler ＋ 偽Consoleコマンド ----------

def _install_fake_console(scenario: str, counter: str | None, no_inchikey: bool) -> None:
    """Consoleの**コマンドライン組み立てだけ**を偽Consoleへ向ける。

    差し替えるのは「何を起動するか」だけで、起動・監視・停止・収集・完了判定は
    本物の`metabolomix.console.execution.supervise`を通る。`supervise`の`command`
    引数はMCPの公開引数にできない（任意コマンドの実行口になる）ため、唯一の
    組み立て場所である`build_msdial_cmd`をこのworkerプロセス内で差し替える
    ——`tests/pipeline_fixtures.py::use_fake_console`と同じ理由・同じ場所で、
    `-i`/`-o`/`-m`/`-p`は実Consoleと同じ意味のまま渡す。
    """
    from metabolomix.console import runner as console_runner

    def build_msdial_cmd(exe, dataset_root, msdial_out_dir, method_file,
                         save_project=False):
        command = [sys.executable, str(_FAKE_CONSOLE), "--scenario", scenario,
                   "-i", str(dataset_root), "-o", str(msdial_out_dir),
                   "-m", str(method_file)]
        if save_project:
            command.append("-p")
        if counter:
            command += ["--counter", str(counter)]
        if no_inchikey:
            command.append("--no-inchikey")
        return command

    console_runner.build_msdial_cmd = build_msdial_cmd


def stage_marker_path(pipeline_root, stage_id: str) -> Path:
    """`--sleep-stage`のstageへ入ったことを知らせる小さな目印ファイル。

    「workerが今このstageに居る」を、状態ファイルをポーリングせずに知るための
    唯一の合図（Windowsでは`pipeline-run.json`を読むだけでworkerの`os.replace`を
    壊しうるため、待ち合わせに状態ファイルを使わない）。取消や強制終了を
    「確実に工程の途中で」当てるのに使う。
    """
    safe = stage_id.replace(":", "_")
    return Path(pipeline_root) / "control" / f"harness-entered-{safe}"


def _with_sleep(handlers: dict, sleep_stage: str | None, sleep_seconds: float) -> dict:
    """指定stageのhandlerだけ、目印を置いてから呼出し**前**に待つラッパを被せる。

    workerを工程の途中で失わせる／stage境界の取消を確実に踏ませるために、
    実時間の窓を作るためだけのもの。handlerの戻り値には一切触れない。
    """
    if not sleep_stage or sleep_seconds <= 0:
        return handlers

    def wrap(func):
        def handler(context: dict) -> dict:
            if context.get("stage_id") == sleep_stage:
                marker = stage_marker_path(context["pipeline_root"], sleep_stage)
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.write_text("entered", encoding="ascii")
                time.sleep(sleep_seconds)
            return func(context)
        return handler

    return {key: wrap(func) for key, func in handlers.items()}


def _real_handlers(args) -> dict:
    from metabolomix.pipeline.service import build_handlers

    _install_fake_console(args.scenario, args.counter, args.no_inchikey)
    return _with_sleep(build_handlers(), args.sleep_stage, args.sleep_seconds)


# ---------- 切り離し起動役 ----------

def _relaunch_detached(args, argv: list) -> int:
    """自分と同じコマンドから`--launch-detached`だけを外して切り離し起動する。

    起動情報（pid・identity）を`--launch-info`へ書いてから即座に終了する
    ——このプロセスが消えたあともworkerが進み続けることを、呼び出し側が
    「起動役の終了」を待ってから観測できるようにするため。
    """
    from metabolomix.core.process_control import launch_detached

    skip_next = False
    forwarded: list = []
    for token in argv:
        if skip_next:
            skip_next = False
            continue
        if token == "--launch-detached":
            continue
        if token == "--launch-info":
            skip_next = True
            continue
        forwarded.append(token)

    command = [sys.executable, str(Path(__file__).resolve()), *forwarded]
    log_path = Path(args.launch_info).with_suffix(".detached.log")
    info = launch_detached(command, cwd=_REPO_ROOT, log_path=log_path)
    Path(args.launch_info).write_text(json.dumps(info), encoding="utf-8")
    return 0


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="pipeline_worker_harness",
        description="run_engineを実プロセスとして起動する（試験専用）。")
    parser.add_argument("--pipeline", required=True)
    parser.add_argument("--sleep-stage", dest="sleep_stage", default=None)
    parser.add_argument("--sleep-seconds", dest="sleep_seconds", type=float, default=0.0)
    parser.add_argument("--scenario", dest="scenario", default=None,
                        help="指定すると本番のbuild_handlers()を使い、Consoleだけ偽物にする")
    parser.add_argument("--counter", dest="counter", default=None)
    parser.add_argument("--no-inchikey", dest="no_inchikey", action="store_true")
    parser.add_argument("--launch-detached", dest="launch_detached", action="store_true")
    parser.add_argument("--launch-info", dest="launch_info", default=None)
    args = parser.parse_args(argv)

    if args.launch_detached:
        if not args.launch_info:
            parser.error("--launch-detached には --launch-info が必要です")
        return _relaunch_detached(args, argv)

    if args.scenario:
        handlers = _real_handlers(args)
    else:
        handlers = _build_handlers(args.sleep_stage, args.sleep_seconds)

    try:
        result = run_engine(Path(args.pipeline), handlers)
    except DomainError as exc:
        print(json.dumps({"error_code": exc.code, "message": exc.message}))
        return 1
    print(json.dumps({"status": result.get("status")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
