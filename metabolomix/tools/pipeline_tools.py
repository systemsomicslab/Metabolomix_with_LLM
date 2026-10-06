"""生データフォルダpipelineのMCPツール（Task18: spec §10.1の5件）。

実体はすべて `metabolomix.pipeline.service` / `metabolomix.pipeline.recovery` にある薄い
ファサード。ここでの責務は3つだけ:

1. MCP引数（`str`/`dict`）を各関数の期待する型（`Path`/`dict`）へ変換して渡す。
2. 戻り値（compactな発送receipt。行列・スコア・loadings・volcano点列は含まない）を
   `json_payload` で返す。
3. `DomainError` を `console_error` の機械可読エンベロープへ変換する。

サービス層のimportをすべて関数本体へ遅延させているのは、他の `*_tools.py`
（`console_tools.py` 等）と同じ流儀——モジュール読込コストを呼び出し時まで
遅らせるためで、循環import防止の意味は無い（サービス層はsession/mcp_coreを
importしないため）。
"""
from __future__ import annotations

from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import console_error
from metabolomix.core.serialization import json_payload

__all__ = [
    "pipeline_plan",
    "pipeline_run",
    "pipeline_status",
    "pipeline_resume",
    "pipeline_cancel",
]

#: pipeline_plan/run/resume/cancelに共通するannotations。
#: readOnlyHint=False（ファイルへ書く）、destructiveHint=False（rawを消さない）、
#: idempotentHint=True（同一request_idの再送・二重cancelは同じ結果に落ち着く。
#: spec §9.3「request_idの再送は同一内容なら同じ結果」「二重取消は冪等」）。
_LOCAL_WRITE = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_plan(dataset_root: str, request: dict | None = None,
                  request_id: str | None = None) -> str:
    """入力検査・不足情報・固定要求の保存のみを行う。Consoleは起動しない。

    Parameters
    ----------
    dataset_root:
        生データフォルダのパス。MS-DIAL が読む 12 形式
        （.abf / .ibf / .cdf / .mzml / .wiff / .raw / .d / .wiff2 /
        .qgd / .lcd / .lrp / .imzml）が入っているフォルダを指します。
        Agilent・Bruker の .d と Waters の .raw は**フォルダ**が 1 検体です。
    request:
        `pipeline-request.v1` のトップレベル項目（省略可、spec §10.1）。
        `target` / `polarity` / `method_file` / `lbm_file` / `sample_manifest` /
        `preprocess` / `comparisons` 等。相対パスは `dataset_root` 基準。
    request_id:
        冪等性キー。同一内容の再送は同じ結果を返し、別内容は `IDEMPOTENCY_CONFLICT`。

    `dataset_root` 直下に `analysis-request.json`（同じトップレベル項目を持つ
    JSON）と `sample-manifest.tsv` があれば既定名として読む。要求値の優先順位は
    **`request` で明示した値 > `analysis-request.json` > 既定値**（spec §7.1)で、
    `preprocess` だけは子キー単位で重なる。どの値がどこから来たかは保存された
    要求の `value_sources`（`explicit` / `request_file` / `default`）で確認できる。
    置いてあるファイルが壊れている・未知キーを持つ場合は黙って無視せず
    `PIPELINE_REQUEST_INVALID` にする。

    receiptの `resolved` に、解決済みの実行条件——`method`（`source_path` /
    `sha256`）・`lbm`・`polarity`（`value` とその `source`）——が入る。何も
    起動しないうちにこの3項目を確認するのがこのツールの用途。

    起動せずに `needs_input` を保存することがある（不正なsample_manifest等、
    spec §9.2「既知の不正入力はConsole起動前に拒否する」）。条件だけ確認したい
    利用者はこちらを使い、実際に進めたいときは `pipeline_run` を使う。
    """
    from metabolomix.pipeline.service import plan_pipeline
    try:
        return json_payload(plan_pipeline(Path(dataset_root), request, request_id))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_run(dataset_root: str, request: dict | None = None,
                 request_id: str | None = None) -> str:
    """計画と起動を一括実行する。生データから MS-DIAL Console を回す補助経路の入口。

    主経路は MS-DIAL GUI で解析済みの出力を `load_dataset` で読むこと。生データ
    だけのフォルダでは、まず GUI での解析を勧め、Console での自動解析は選択肢
    として示す。このツールはユーザーがそれを選んだときにだけ呼ぶ。起動前に解決
    結果（method file・LBM・極性）を見たいときは `pipeline_plan`。引数は
    `pipeline_plan` と同じ。

    短時間で `pipeline_path` を返す——起動受理は、workerのidentity保存と
    起動handshakeを確認した後に返すが、応答待ちには短い上限があり、上限に
    達しても起動失敗と決め付けず再起動もしない（`launch.handshake` が
    `"not_confirmed"` になるだけで、`pipeline_status` で追跡できる）。

    `launch` キーは実際にworkerを起動した（`launched=true`）ときだけ含まれる。
    不正なsample_manifest等を検出した場合や、`find_or_create_run` が既存run
    （活動中／completed／失敗・取消・部分完了済み）を再利用した場合は
    `launched=false` のみを返し、`launch` キー自体を持たない（起動していない
    以上、handshake結果も存在しない）。既存runを本当に進めたい場合は
    `pipeline_resume` を使う。
    """
    from metabolomix.pipeline.service import start_pipeline
    try:
        return json_payload(start_pipeline(Path(dataset_root), request, request_id))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def pipeline_status(pipeline_path: str, include_details: bool = False) -> str:
    """状態・不足情報・成果物参照を読む。最終化や起動は一切しない（読取専用）。

    Parameters
    ----------
    pipeline_path:
        `pipeline_run` / `pipeline_plan` が返した `pipeline_path`（pipeline-run.json
        そのもの、またはその親ディレクトリのどちらでもよい）。
    include_details:
        True なら永続レコード全体（`record`）も返す（既定False。通常は
        `status` / `stage_statuses` / `needs_input` / `warnings` で足りる）。

    `status="running"` のときだけ、記録済みworkerの生存確認 `observed_health`
    を追加で行う。workerが消えていれば `observed_health="worker_missing"` と
    `recovery_hint`（`pipeline_resume` を促す）を返すが、statusフィールド自体は
    書き換えない——監視・成果物確定をポーリング呼出しに依存させない（spec）。
    """
    from metabolomix.pipeline.recovery import read_status
    try:
        return json_payload(read_status(Path(pipeline_path), include_details=include_details))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_resume(pipeline_path: str, updates: dict | None = None,
                    request_id: str | None = None, rerun_upstream: bool = False) -> str:
    """入力訂正・下流revision作成・停止工程からの再開。

    Parameters
    ----------
    pipeline_path:
        対象runの `pipeline_path`（`pipeline_status` と同じく、`pipeline-run.json`
        そのもの、またはその親ディレクトリ＝pipeline_rootのどちらでもよい）。
    updates:
        変更したい項目のみ（`target` / `sample_manifest` / `preprocess` /
        `comparisons` に限定。spec §10.1「resumeで変更可能なのは...に限定する」）。
    request_id:
        冪等性キー。同一内容の再送は新しいrevisionを作らず直前の結果を返す。
    rerun_upstream:
        True を明示しない限り、Console実行はやり直さない（自動再試行はしない）。
        監視を失ったまま稼働中のConsoleがある場合は `EXECUTION_UNRESOLVED` で
        待つ（終了を確認できないものを勝手に完了扱いしない）。
    """
    from metabolomix.pipeline.service import resume_pipeline
    try:
        return json_payload(
            resume_pipeline(Path(pipeline_path), updates, request_id, rerun_upstream))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_cancel(pipeline_path: str) -> str:
    """取消要求を保存する。二重取消は冪等。

    **rawデータを削除しない。** ここで保存するのは協調的な取消フラグだけで、
    上流Consoleの停止は監視経路が、下流の停止はstage境界が行う。受理
    （このツール呼び出しが成功したこと）と、実際にworkerが停止を確認した
    こと（`cancelled` への確定）は別状態——確定は `pipeline_status` で確認する。
    """
    from metabolomix.pipeline.recovery import request_cancel
    try:
        return json_payload(request_cancel(Path(pipeline_path)))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)
