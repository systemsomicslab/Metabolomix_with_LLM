"""生データフォルダpipelineの独立ワーカー（CLIエントリ）。

**sessionを一切importしない。** `metabolomix.core.session_state` / `metabolomix.core.mcp_core` /
`metabolomix.tools.*` はここからimportしてはいけない（`tests/test_pipeline_engine.py` が
ASTで検査する）。解析の進行はすべて `pipeline-run.json`（`metabolomix.pipeline.store`）と
このプロセスだけが持つ `runtime` に住み、MCPサーバの共有sessionへは一切触れない。

実際の工程handler一式（`build_handlers`）はTask18の `metabolomix.pipeline.service` が
実装する。ここでは意図的に**遅延import**にする——モジュール読込時点（`python -m` の
import、あるいは本モジュールをASTだけ検査する試験）で `service.py` の不在により
落ちないようにするため。Task18が実装するまで、このモジュールはpublicなMCPツールへ
登録しない（brief拘束）。

CLI:

    python -m metabolomix.pipeline.worker --pipeline <pipeline_root>

試験専用の注入口は持たない。任意のimport名・コマンド・環境変数を受け付ける
入力経路は一切作らない——それ自体がコード実行の入口になる（brief拘束）。
engineを注入handlerで検証する試験は、この production worker とは別の
`tests/pipeline_worker_harness.py` を使う（`metabolomix.console.worker` の
`--console-arg` と同じく、MCPからは到達できない試験専用の経路）。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

__all__ = ["main", "run_worker"]


def run_worker(pipeline_path: Path) -> dict:
    """本番handler一式（Task18 `build_handlers`）でrun_engineを1回分進める。

    importをここへ遅延させているのは、`metabolomix.pipeline.service` が無くても
    本モジュール自体のimport（や、それを静的にASTだけ検査する試験）を
    壊さないため。
    """
    from metabolomix.pipeline.engine import run_engine
    from metabolomix.pipeline.service import build_handlers  # Task18が実装する

    handlers = build_handlers()
    return run_engine(Path(pipeline_path), handlers)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="metabolomix.pipeline.worker",
        description="pipeline-run.json を1回分進める独立ワーカー（sessionを持たない）。")
    parser.add_argument("--pipeline", required=True,
                        help="pipeline_root（pipeline-run.json のあるディレクトリ）")
    args = parser.parse_args(argv)

    result = run_worker(Path(args.pipeline))
    print(f"status={result.get('status')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
