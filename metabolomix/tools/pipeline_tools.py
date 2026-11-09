"""生データフォルダpipelineのMCPツール（5件。pipeline は上流＝Console と解析フォルダの受領まで）。

実体はすべて `metabolomix.pipeline.service` / `metabolomix.pipeline.recovery` にある薄い
ファサード。ここでの責務は3つだけ:

1. MCP引数（`str`/`dict`）を各関数の期待する型（`Path`/`dict`）へ変換して渡す。
2. 戻り値（compactな発送の要約。行列・スコアは含まない）を
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
        `pipeline-request.v3`（省略可。`schema` 省略・`pipeline-request.v1` も受けて v3 で
        保存する。`pipeline-request.v2` は `PIPELINE_REQUEST_INVALID`）。共有の項目は
        `library_mode` / `measure` / `keep_extension` / `timeout_s` / `save_project` /
        `output_root` / `alignment_reference`。1 極性の形は `polarity` / `method_file` /
        `lbm_file` / `msp_file` / `sample_manifest` をトップに置く。2 極性の形は
        `units: [{"polarity", "raw_dir", "method_file"?, "lbm_file"?, "msp_file"?,
        "sample_manifest"?, "keep_extension"?}]`（1〜2 件。`polarity` と `raw_dir`
        （`dataset_root` 相対）は必須で、重ならない。単位の `sample_manifest` の既定は
        `<raw_dir>/sample-manifest.tsv`）。`units` とトップの単位の項目は混ぜられない。
        下流の項目（`target` / `preprocess` / `comparisons`）は `PIPELINE_DOWNSTREAM_REMOVED`
        （`details.replacements`: `load_dataset` → `dataset_preprocess` / `dataset_compare` /
        `dataset_pca`）。相対パスは `dataset_root` 基準。`library_mode="msp_only"` は LBM を使わず研究室の参照
        ライブラリ（.msp）だけで同定する（`msp_file` か極性の設定 `[library]
        msp_positive` / `msp_negative`。メソッドの Text DB・注釈器設定表の宣言は空にする）。
        `library_mode="msp_then_lbm"`（脂質の run 用）は研究室 MSP（priority 2 の注釈器設定表）を
        優先し、`[library] lbm` の LBM で補う（LBM の照合で RT を使わない。LBM は
        メソッドの宣言・実行体のフォルダを見ない。無ければ `LBM_NOT_CONFIGURED`）。
        `lbm_file` との同時指定（msp_only / msp_then_lbm）、`auto` での `msp_file` は
        `PIPELINE_REQUEST_INVALID`。
        `alignment_reference`（既定 `"median_injection"`）は MS-DIAL の
        `Alignment reference file ID` を、全試料（`include` を問わない）の注入順の中央
        （偶数なら小さい方）の試料にする。注入順は実験情報シートの `injection_order` から
        取り、1 件でも欠けるときはファイル名順の中央（戻り値の
        `injection_order_source` が `"filename"`）。`"method"` はメソッドの宣言をそのまま使う。
        resume では変えられない（`NEW_PIPELINE_REQUIRED`）。
    request_id:
        冪等性キー。同一内容の再送は同じ結果を返し、別内容は `IDEMPOTENCY_CONFLICT`。

    `dataset_root` 直下に `analysis-request.json`（同じ項目を持つ JSON）と
    `sample-manifest.tsv` があれば既定名として読む。要求値の優先順位は
    **`request` で明示した値 > `analysis-request.json` > 既定値**（spec §7.1)
    （`units` は層ごとに丸ごと）。どの値がどこから来たかは保存された要求の
    `value_sources`（`explicit` / `request_file` / `default`）で確認できる。置いてある
    ファイルが壊れている・未知キーを持つ場合は黙って無視せず `PIPELINE_REQUEST_INVALID`
    にする。ただしファイルの下流の項目（`target` / `preprocess` / `comparisons`）は
    読み飛ばして警告 `REQUEST_FILE_DOWNSTREAM_IGNORED` を返す（Console は回す）。

    戻り値の `resolved.units.<極性>` に、単位ごとの解決済みの実行条件——`method`（`source_path` /
    `sha256`）・`lbm`（msp_then_lbm ではファイル名・sha256・出どころだけ）・`polarity`
    （`value` とその `source`）・`library_mode`・`msp`（ファイル名・sha256・出どころ、
    msp_then_lbm では `priority` も。置き場所は載せない）・`alignment_reference`
    （`file_id`・`sample_file`・`injection_order_source`・`method_value`。v8 より前の計画から
    再利用した run では null）・msp_only / msp_then_lbm では実効の足切り `dot_product_cutoffs`
    （`weighted`・`simple`・`reverse`・`scale: "square_root"`・種類ごとの `source`）——が入る。何も
    起動しないうちにこれらを確認するのがこのツールの用途。

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
    """計画と起動を一括実行する。生データから MS-DIAL Console を回す、対等な 2 つ目の入口。

    結果は極性ごとの解析フォルダ（`<run>/<極性>/input/`。MS-DIAL の `.arf` 等）と受領
    （`results/receipt.json`）で、続きは `load_dataset(pipeline_path=…, polarity=…)` から
    （前処理・比較・PCA は `dataset_preprocess` / `dataset_compare` / `dataset_pca`）。
    先にフォルダに MS-DIAL を動かした形跡（.arf 等の出力）があるかを確かめ、
    あれば `load_dataset` で進む。形跡が無いときだけ、ユーザーに Console を回してよいか
    許可を取り、同意を得てからこのツールを呼ぶ。起動前に解決結果（method file・LBM・
    極性）を見たいときは `pipeline_plan`。引数は `pipeline_plan` と同じ。2 極性の要求は
    要求の順に 1 極性ずつ回し、片方が失敗してももう片方を続ける（run は `partial`）。

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
        `status` / `units` / `stage_statuses` / `needs_input` / `warnings` で足りる）。

    `units` は極性ごとの状態（`pending` / `running` / `validated` / `failed`）。受領を
    書いた後は `receipt`（検証を通った極性の解析フォルダの場所・照合した `.arf`・
    次に呼ぶ `load_dataset(...)`）も返す。旧形式（pipeline-run.v1）の run も読め、
    `legacy: true` と注意 1 行が付く（再開・取消はできない）。

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
    """試料の対応表の訂正・停止した工程からの再開。

    Parameters
    ----------
    pipeline_path:
        対象runの `pipeline_path`（`pipeline_status` と同じく、`pipeline-run.json`
        そのもの、またはその親ディレクトリ＝pipeline_rootのどちらでもよい）。
    updates:
        `sample_manifest`（1 極性の run）か `units: {<極性>: {"sample_manifest": …}}` だけ。
        変えると受領だけを書き直す。下流の項目（`target` / `preprocess` / `comparisons`）は
        `PIPELINE_DOWNSTREAM_REMOVED`、上流の条件（メソッド・ライブラリ・極性…）は
        `NEW_PIPELINE_REQUIRED`。旧形式（pipeline-run.v1）の run は `PIPELINE_RUN_LEGACY`
        （読むなら `load_dataset(pipeline_path=…)`、Console をやり直すなら新しい `pipeline_run`）。
    request_id:
        冪等性キー。同一内容の再送は新しいrevisionを作らず直前の結果を返す。
    rerun_upstream:
        True を明示しない限り、Console実行はやり直さない（自動再試行はしない）。
        True は upstream が成功していない極性だけをやり直す（新しい attempt）。
        省くと、失敗した極性はそのままにして他の極性を進める（`not_retried_units`）。
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
    Consoleの停止は監視経路が、工程の停止はstage境界が行う。旧形式（pipeline-run.v1）の
    run は `PIPELINE_RUN_LEGACY`。受理
    （このツール呼び出しが成功したこと）と、実際にworkerが停止を確認した
    こと（`cancelled` への確定）は別状態——確定は `pipeline_status` で確認する。
    """
    from metabolomix.pipeline.recovery import request_cancel
    try:
        return json_payload(request_cancel(Path(pipeline_path)))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)
