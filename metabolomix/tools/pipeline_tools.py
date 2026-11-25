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
    """入力を検査し、不足情報と固定した要求を保存する。Console は起動しない（起動は `pipeline_run`）。

    Parameters
    ----------
    dataset_root:
        生データフォルダ。MS-DIAL が読む 12 形式（.abf / .ibf / .cdf / .mzml / .wiff / .raw / .d /
        .wiff2 / .qgd / .lcd / .lrp / .imzml）。Agilent・Bruker の .d と Waters の .raw はフォルダが 1 検体。
    request:
        `pipeline-request.v3`（省略可。`schema` 省略と `pipeline-request.v1` は v3 で保存、
        `pipeline-request.v2` は `PIPELINE_REQUEST_INVALID`）。
        共有の項目: `library_mode` / `measure` / `keep_extension` / `timeout_s` / `save_project` /
        `output_root` / `alignment_reference`。
        1 極性: `polarity` / `method_file` / `lbm_file` / `msp_file` / `sample_manifest` / `method_overrides`
        をトップに置く。
        2 極性: `units: [{"polarity", "raw_dir", "method_file"?, "lbm_file"?, "msp_file"?,
        "sample_manifest"?, "keep_extension"?, "method_overrides"?}]`（1〜2 件。`polarity` と `raw_dir`（`dataset_root` 相対）は
        必須で重ならない。`sample_manifest` の既定は `<raw_dir>/sample-manifest.tsv`）。`units` とトップの
        単位の項目は混ぜられない。相対パスは `dataset_root` 基準。
        知らない項目は `PIPELINE_REQUEST_INVALID`（`details.unknown_keys`）。消えた下流の項目
        （`target` / `preprocess` / `comparisons`）は文面と `details.replacements`
        （`load_dataset` → `dataset_preprocess` / `dataset_compare` / `dataset_pca`）で案内する。
        `library_mode="msp_only"`: LBM を使わず研究室の .msp だけで同定する（`msp_file` か `[library]
        msp_positive` / `msp_negative`。メソッドの Text DB・注釈器設定表の宣言は空にする）。
        `library_mode="msp_then_lbm"`（脂質の run 用）: 研究室 MSP（注釈器設定表の priority 2）を優先し、
        `[library] lbm` の LBM で補う（LBM の照合は RT を使わず、メソッドの宣言・実行体のフォルダを見ない。
        無ければ `LBM_NOT_CONFIGURED`）。
        `lbm_file` と msp_only / msp_then_lbm の併用、`auto` での `msp_file` は `PIPELINE_REQUEST_INVALID`。
        `alignment_reference`（既定 `"median_injection"`）: MS-DIAL の `Alignment reference file ID` を、全試料
        （`include` を問わない）の注入順の中央（偶数なら小さい方）にする。注入順は実験情報シートの
        `injection_order`、1 件でも欠ければファイル名順（`injection_order_source: "filename"`）。`"method"` は
        メソッドの宣言のまま。resume では変えられない（`NEW_PIPELINE_REQUIRED`）。
        `method_overrides`（単位ごと）: `{Console のキー: 値}`。メソッドの原本は変えず、実効メソッドの値を
        上書きする（例 `{"Minimum peak height": 50000, "N percent detected in one group": 50}`）。同じ値の
        受け先を持つ別名が原本にあれば、それも同じ値にする。Console が読むキーと型（MsdialWorkbench afd5f95 の
        ConfigParser）で検査し、知らないキー（`details.suggestions` に近い綴り）・現用の Console のビルドが読まない
        キー（`Number of threads`・`Lbm annotator priority`・`Detailed alignment provenance`・`Annotation candidates`
        とその別名）・型に合わない値（`Minimum peak
        height` は整数。`1e5` は不可）・ほかの欄が受け持つキー（`Ion mode`・ライブラリやファイルのパス。
        `Alignment reference file ID` は `alignment_reference="method"` のときだけ）・計画が自分で書くキーと
        重なるものは `PIPELINE_REQUEST_INVALID`。上書きが変われば別の run になり、resume では変えられない。
    request_id:
        冪等性キー。同じ内容の再送は同じ結果、別内容は `IDEMPOTENCY_CONFLICT`。

    `dataset_root` 直下の `analysis-request.json`（同じ項目の JSON）と `sample-manifest.tsv` を既定で読む。
    優先順は `request` > `analysis-request.json` > 既定値（`units` は層ごとに丸ごと）。出どころは保存した
    要求の `value_sources`（`explicit` / `request_file` / `default`）。ファイルが壊れている・未知キーを持つなら
    `PIPELINE_REQUEST_INVALID`。ただしファイルの下流の項目は読み飛ばして警告 `REQUEST_FILE_DOWNSTREAM_IGNORED`
    （Console は回す）、明示の `units` に負けたファイルのトップの単位のキーは警告 `REQUEST_FILE_UNIT_KEYS_IGNORED`。

    戻り値の `resolved.units.<極性>`: `method`（`source_path` / `sha256` / 結果を大きく変えるキーの実効値
    `key_params` と宣言の無い（Console の既定で走る）キー `key_params_undeclared` / 上書きしたキーの
    `user_overrides`（`{"from": 原本の値か null, "to": 値}`））・`lbm`（msp_then_lbm ではファイル名・
    sha256・出どころだけ）・`polarity`（`value` / `source`）・`library_mode`・`msp`（ファイル名・sha256・出どころ、
    msp_then_lbm では `priority`。置き場所は載せない）・`alignment_reference`（`file_id`・`sample_file`・
    `injection_order_source`・`method_value`。この欄の無い古い計画を再利用した run では null）・msp_only /
    msp_then_lbm では実効の足切り `dot_product_cutoffs`（`weighted`・`simple`・`reverse`・`scale: "square_root"`・
    種類ごとの `source`）。

    既知の不正入力（不正な sample_manifest など）は起動せずに `needs_input` を保存する。開いた形式（mzML・imzML・
    mzXML）の計測ファイルが閉じタグで終わらない・32 KiB 未満（imzML は閉じタグだけ。スペクトルは隣の .ibd）・
    mzML の見出しに中身の無い自己終了の `referenceableParamGroup` など（MS-DIAL がスペクトルを 1 本も読まず、Console は
    exit 0 のままピーク 0 で終わる）なら、run を作る前に `PIPELINE_INPUT_INVALID`（`details.files` に
    `{"name", "size", "reason"}`。見出しの問題は `reason: "empty_header_element"` と `element`）。
    """
    from metabolomix.pipeline.service import plan_pipeline
    try:
        return json_payload(plan_pipeline(Path(dataset_root), request, request_id))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_run(dataset_root: str, request: dict | None = None,
                 request_id: str | None = None) -> str:
    """計画と起動を一度に行う。生データから MS-DIAL Console を回す（`load_dataset` と対等な入口）。
    引数は `pipeline_plan` と同じ（起動前に解決結果を見るのは `pipeline_plan`）。

    出力は極性ごとの解析フォルダ（`<run>/<極性>/input/`。`.arf` 等）と受領（`results/receipt.json`）。続きは
    `load_dataset(pipeline_path=…, polarity=…)`。MS-DIAL の出力が既にあるフォルダは Console を回さず
    `load_dataset` で開ける。2 極性は要求の順に 1 極性ずつ回し、片方が失敗しても続ける（run は `partial`）。
    実験情報シートがあれば Console の `-i` に試料の表（解析フォルダの `console-input.csv`。群・役割・注入順）を
    渡し、`N percent detected in one group` を群ごとに効かせる。無ければフォルダ（各ファイルが別の群）。
    どちらかは受領の `console.input` と `processing` / `method_params` の `console_input`。

    worker の identity の保存と起動 handshake を確かめて `pipeline_path` を返す。待ちには短い上限があり、
    達しても失敗扱い・再起動はせず `launch.handshake: "not_confirmed"`（`pipeline_status` で追う）。
    `launch` は実際に起動した（`launched=true`）ときだけ付く。不正な sample_manifest など、または既存の run
    （活動中・completed・失敗・取消・部分完了）を再利用したときは `launched=false` だけ（既存の run を進めるのは
    `pipeline_resume`）。

    `launch.worker_in_job: true` は worker が Windows の Job（何段目でも）に残っていること。中継が直近の Job を
    抜けたのに worker が Job にいて、その Job が閉じたら止めるか分からないなら、Job が閉じると worker と Console
    も止まりうるので警告 `WORKER_IN_JOB`（`details.breakaway` / `relay_breakaway`。止まったら `pipeline_status` →
    `pipeline_resume`）。中継の直近の Job が閉じても止めないもの（`launch.relay_breakaway: "inherited_job"`）なら
    警告せず、外側の Job は確かめていないことを `launch.outer_jobs: "unchecked"` と `launch.outer_jobs_note` に
    書く（run の記録の `worker` にも同じ注記）。

    Console は同じ設定ファイルを読むサーバ・worker の間で同時に 1 本（run を跨ぐ。待ちは `pipeline_status` の
    `waiting`）。起動前の commit の空きが `[msdial] min_free_commit_gb`（既定 8 GB）を下回れば警告
    `LOW_COMMIT_MEMORY`（止めない）。Console が完了しなければ `MSDIAL_EXECUTION_FAILED` に `details.cause`
    （`system_resources`: 終了コード 1450 などの資源不足 / `broken_input` / `input_suspected`）と
    `details.diagnosis`（疑わしいファイル・ログのエラーの行）が付く（`pipeline_status` の `errors` でも読める）。
    """
    from metabolomix.pipeline.service import start_pipeline
    try:
        return json_payload(start_pipeline(Path(dataset_root), request, request_id))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def pipeline_status(pipeline_path: str, include_details: bool = False) -> str:
    """run の状態・不足情報・成果物の参照を読む（読み取り専用。最終化も起動もしない）。

    Parameters
    ----------
    pipeline_path:
        `pipeline_run` / `pipeline_plan` が返した `pipeline_path`（`pipeline-run.json` かその親フォルダ）。
    include_details:
        True なら永続レコード全体（`record`）も返す（既定 False。普段は `status` / `units` /
        `stage_statuses` / `errors` / `needs_input` / `warnings` で足りる）。

    `units` は極性ごとの状態（`pending` / `running` / `validated` / `failed`）。受領を書いた後は `receipt`
    （検証を通った極性の解析フォルダ・照合した `.arf`・次に呼ぶ `load_dataset(...)`・`method_params`: 極性ごとの
    実効メソッドの `key_params` / `undeclared` / `overrides`）も返す。値を変えて処理し直すのは、要求に
    `method_overrides` を足した `pipeline_run`（新しい run）。旧形式
    （pipeline-run.v1）の run も読め、`legacy: true` と注意 1 行が付く（再開・取消はできない）。

    `status="running"` のときだけ worker の生存を確かめ、居なければ `observed_health="worker_missing"` と
    `recovery_hint`（`pipeline_resume` を促す）を返す（`status` は書き換えない）。この機械のほかの run の Console
    を待つ工程は `waiting`（工程ごとの `reason: "console_slot"`・`since`・相手の `holder`）。待っていた worker が
    居なければ `waiting` ではなく `interrupted_while_waiting`（中断時の待ちの記録）に出し、Console をまだ起こして
    いない極性を `recovery_hint.console_not_started_units` に挙げる（`pipeline_resume` は `rerun_upstream` なしで
    その極性の最初の Console から続ける）。worker が、閉じたら止めるか分からない Windows の Job に残っていれば
    `warnings` に `WORKER_IN_JOB`。

    失敗・取消で終わった工程は `errors`（`{<stage_id>: {"code", "message", "cause"?, "suspect_files"?}}`）。
    `cause` は Console が完了しなかった原因（`system_resources` / `broken_input` / `input_suspected`）、
    `suspect_files` は `[{"name", "reason"}]`（配置の工程で断った切れたファイルも）。ログの行
    （`details.diagnosis.log_lines`）は `include_details=True` の `record` の工程の `error` にある。
    """
    from metabolomix.pipeline.recovery import read_status
    try:
        return json_payload(read_status(Path(pipeline_path), include_details=include_details))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_resume(pipeline_path: str, updates: dict | None = None,
                    request_id: str | None = None, rerun_upstream: bool = False) -> str:
    """試料の対応表を訂正する・止まった工程から再開する。

    Parameters
    ----------
    pipeline_path:
        対象 run の `pipeline_path`（`pipeline-run.json` かその親フォルダ）。
    updates:
        `sample_manifest`（1 極性の run）か `units: {<極性>: {"sample_manifest": …}}` だけ。変えると受領だけを
        書き直す。成功済みの極性の群別絞り込みも訂正後のシートで再計算する場合は、訂正後のシートを指定し
        request_id を省いて pipeline_run を呼ぶ。知らない項目は `PIPELINE_REQUEST_INVALID`（消えた下流の項目 `target` / `preprocess` /
        `comparisons` なら `details.replacements`）、上流の条件（メソッド・ライブラリ・極性…）は
        `NEW_PIPELINE_REQUIRED`。旧形式（pipeline-run.v1）の run は `PIPELINE_RUN_LEGACY`（読むなら
        `load_dataset(pipeline_path=…)`、Console をやり直すなら新しい `pipeline_run`）。
    request_id:
        冪等性キー。同じ内容の再送は新しい revision を作らず、1 回目と同じ内訳（`console_not_started_units` /
        `not_retried_units` / `held_units`）を `reused: true` で返す。
    rerun_upstream:
        既定 False は Console をやり直さない（自動再試行もしない）。失敗した極性はそのままで他の極性を進め
        （`not_retried_units`）、まだ回していない極性の Console は起こす（1 回目）。
        True は upstream が成功していない極性だけを新しい attempt でやり直す。生きている worker が回している
        極性の Console は捨てず `held_units` に並べる。
        worker が居なくなって `running` のままの極性は Console の状態が分からないので `UPSTREAM_RERUN_REQUIRED`
        （True を求める）。ただし Console の順番待ちの間に worker が消え、Console をまだ起こしていない極性（今の
        attempt の `analysis-job.json` が無い。`pipeline_status` の `recovery_hint.console_not_started_units`）は、
        False でもその極性の最初の Console から続け、戻り値の `console_not_started_units` に並べる。
        監視を失ったまま稼働中の Console があれば `EXECUTION_UNRESOLVED`（終了を確かめるまで完了扱いしない）。

    起こした worker の Job の扱いは `pipeline_run` と同じ（警告 `WORKER_IN_JOB`、または `inherited_job` なら
    `launch.outer_jobs: "unchecked"` と `launch.outer_jobs_note`）。
    """
    from metabolomix.pipeline.service import resume_pipeline
    try:
        return json_payload(
            resume_pipeline(Path(pipeline_path), updates, request_id, rerun_upstream))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False)
def pipeline_cancel(pipeline_path: str) -> str:
    """取消を要求する（協調的な取消フラグを保存するだけ。二重取消は冪等）。raw データは消さない。

    Console は監視経路が、工程は stage の境目で止める。受理（この呼び出しの成功）と停止の確定（`cancelled`）は
    別で、確定は `pipeline_status` で見る。旧形式（pipeline-run.v1）の run は `PIPELINE_RUN_LEGACY`。
    """
    from metabolomix.pipeline.recovery import request_cancel
    try:
        return json_payload(request_cancel(Path(pipeline_path)))
    except DomainError as exc:
        return console_error(exc.code, exc.message, exc.details or None)
