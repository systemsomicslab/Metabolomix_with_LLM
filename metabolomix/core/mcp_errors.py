"""クライアント非依存のエラーエンベロープ。

MCP クライアントが「どのツールを先に呼べばよいか」を、サーバ固有の日本語文面を
解釈せずに機械的に判断できるようにする。汎用クライアント（Use-LLLM 等）は
本文を JSON として読み、`error.code` が `missing_state` なら
`error.required_tools` を再実行して状態を復元できる。

なぜ本文の JSON なのか: 使用中の MCP SDK では `isError=true` を作る経路
（lowlevel/server.py の _make_error_result）がテキストのみを返して
structuredContent を捨てる。さらに FastMCP は戻り値アノテーションから
outputSchema を導出するため、`-> str` のツールが失敗時だけ dict を返すこともできない。
機械可読なエラーはテキストに載せるしかない。

もう 1 つの封筒 `needs_input`（`{"needs_input": {field, reason, options, argument?, retry_tools}}`）は
**エラーではない**（`error` キーを持たない）。既定を決めずに利用者へ尋ねる問い（例: 内部標準で
換算済みの値に QC のドリフト補正を重ねるか）で、クライアントは `options` を利用者に示し、選ばれた
`value` を `argument`（無ければ spec の `field`）に入れて `retry_tools` のどれかを呼び直す。

このモジュールは依存グラフの leaf（stdlib と同じく leaf の
metabolomix.core.serialization と metabolomix.core.atomic_io のみ）。tools_* / server を import しない。
"""
import functools

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.serialization import json_payload

MISSING_STATE = "missing_state"


def missing_state(state: str, required_tools: list[str], message: str) -> str:
    """前提となるセッション状態が無いことを、契約どおりのエンベロープで返す。

    引数:
        state: 欠けている状態の識別子（例 "preprocessed_matrix"）。クライアントは
            これを不透明な文字列として扱い、解釈しない。開示とログのためにある。
        required_tools: その状態を作れるツール名の**代替候補（OR）**。サーバ名は
            含めない（クライアント側の名前空間はクライアントが決める）。先頭ほど優先。
        message: LLM と人間が読む説明。既存の日本語文面をそのまま渡すこと。
            エンベロープを解釈しないクライアントではこれだけが見える。

    3つとも欠かせない。復旧の手掛かりが無いエンベロープを無言で出すと、
    クライアントは「状態不足だが何もできない」状態に陥り、原因究明を誤らせる。
    """
    if not state:
        raise ValueError("state は必須です（欠けている状態の識別子）。")
    if not required_tools:
        raise ValueError("required_tools は必須です（状態を作れるツールの候補）。")
    if not message:
        raise ValueError("message は必須です（LLM と人間が読む説明）。")
    return json_payload({
            "error": {
                "code": MISSING_STATE,
                "state": state,
                "required_tools": list(required_tools),
                "message": message,
            }
        })


NEEDS_INPUT = "needs_input"


def needs_input(field: str, reason: str, options: list[dict], *, argument: str | None = None,
                retry_tools: list[str]) -> str:
    """利用者に尋ねてから呼び直す封筒を返す（エラーではない）。

    引数:
        field: 答えが入る設定の欄（spec の路。例 "preprocess.drift.mode"）。
        reason: なぜ尋ねるのか（LLM と人間が読む説明）。
        options: 選択肢 `[{"value", "label"}, ...]`。`value` を引数に入れて呼び直す。
        argument: 呼び直すツールの引数名（例 "drift_mode"）。spec で答える経路
            （`dataset_analysis_run` の v2）では渡さず、封筒にも載せない（`field` が答えの場所）。
        retry_tools: 答えを入れて呼び直せるツール名の候補（OR）。

    `missing_state` と同じく、復旧の手掛かりが欠けた封筒は出さない（ValueError）。
    """
    if not field:
        raise ValueError("field は必須です（答えが入る設定の欄）。")
    if not reason:
        raise ValueError("reason は必須です（尋ねる理由）。")
    if not options or any(not isinstance(o, dict) or "value" not in o or "label" not in o
                          for o in options):
        raise ValueError("options は value と label を持つ選択肢の 1 つ以上のリストです。")
    if not retry_tools:
        raise ValueError("retry_tools は必須です（答えを入れて呼び直せるツールの候補）。")
    body: dict = {"field": field, "reason": reason, "options": [dict(o) for o in options]}
    if argument:
        body["argument"] = argument
    body["retry_tools"] = list(retry_tools)
    return json_payload({NEEDS_INPUT: body})


#: `argument`（drift_mode）を受けて呼び直せるツール（`needs_input` の retry_tools の意味）と、その引数名。
#: 受けるのは dataset_preprocess だけ。比較・PCA・パスウェイは返った matrix_id を渡して呼び直す。
DRIFT_ARGUMENT = "drift_mode"
DRIFT_RETRY_TOOLS = ["dataset_preprocess"]


def default_matrix_needs_input(needs, tool: str) -> str:
    """既定の前処理が利用者に尋ねたいとき（`needs` は前処理の `NeedsInput`）の封筒。

    `tool` は drift_mode の引数を持たないツール（比較・PCA・パスウェイ）の名前で、理由に呼び直し先
    （`dataset_preprocess(drift_mode=...)` で行列を作り、その matrix_id を渡す）を書く。
    """
    return needs_input(
        needs.field,
        f"{needs.reason} {tool} には drift_mode 引数がありません。"
        f"dataset_preprocess(drift_mode=...) で行列を作り、その matrix_id を "
        f"{tool} に渡して呼び直してください。",
        needs.options, argument=DRIFT_ARGUMENT, retry_tools=list(DRIFT_RETRY_TOOLS))


CONSOLE_ERROR_CODES = frozenset({
    "MSDIAL_EXE_NOT_FOUND",
    # 設定ファイル（lipidmix.local.toml）が読めない。
    "CONFIG_INVALID",
    "MSDIAL_TIMEOUT",
    "MSDIAL_NONZERO_EXIT",
    "NO_JOB_OUTPUT",
    "JOB_NOT_FOUND",
    "JOB_NOT_PLANNED",
    "DATASET_ROOT_IN_REPO",
    "METHOD_FILE_NOT_FOUND",
    # MS-DIAL Console が読む ASCII の key/value テキストではない場合。
    "METHOD_FILE_NOT_TEXT",
    # [msdial] exe が Console ではなく GUI を指している場合。
    "MSDIAL_EXE_NOT_CONSOLE",
    # Console の起動検査（--help）が時間切れを繰り返した場合。Console でないとは
    # 判定できていない（冷間起動の遅れでも起きる）ので NOT_CONSOLE と分ける。
    "MSDIAL_EXE_PROBE_TIMEOUT",
    # MS-DIAL が対象とする計測拡張子が混在する、または 0 種類の場合。
    "MIXED_RAW_FORMATS",
    # method_file 省略時に、使えるパラメータファイルが見つからなかった場合。
    "METHOD_FILE_NOT_GIVEN",
    # メソッドファイルの Ion mode と宣言 polarity が食い違う場合。
    "METHOD_FILE_POLARITY_MISMATCH",
    # 脂質ライブラリ（.lbm2）を解決できない／候補が 1 件に絞れない場合。
    # 空のまま実行すると MS-DIAL は警告なしで同定 0 件のまま完走する。
    "LBM_NOT_FOUND",
    "LBM_AMBIGUOUS",
    # pipeline の library_mode="msp_then_lbm" で `[library] lbm` が未設定の場合（v8）。
    "LBM_NOT_CONFIGURED",
    # console_prepare_input が単一フォーマットの入力フォルダを作れなかった場合。
    "INPUT_PREP_FAILED",
    # MS-DIAL 実行自体は成功したが、その後の生成物収集・ハッシュ計算・
    # analysis-job.json への保存で失敗した場合。NO_JOB_OUTPUT（生成物が
    # そもそも無い）とは別のエラー: こちらは生成物はあるが確定処理に失敗した状態。
    "JOB_POST_RUN_FAILED",
})


def console_error(code: str, message: str, details: dict | None = None,
                  required_tools: list[str] | None = None) -> str:
    """Console 実行層固有のエラーを機械可読エンベロープで返す。

    `required_tools` は missing_state と同じ意味の**代替候補（OR）**で、
    「このエラーを解消できるツール」を指す。文面がいくら正しくても、
    MCP クライアントに実行できない手順（フォルダを作る等）を指示していては
    復旧できない。封筒が自分の直し方を機械可読に指すためにある。
    """
    if not code or not message:
        raise ValueError("code と message は必須です。")
    payload: dict = {"code": code, "message": message}
    if details is not None:
        payload["details"] = details
    if required_tools:
        payload["required_tools"] = list(required_tools)
    return json_payload({"error": payload})


MZTAB_ERROR_CODES = frozenset({
    "MZTAB_NOT_FOUND",
    "MZTAB_STRUCTURE_INVALID",
    "QUANTIFICATION_CONFLICT",
    "AMBIGUOUS_PRIMARY_MZTAB",
    "UNSUPPORTED_AREA_CONSOLE",
    "SAMPLE_DESIGN_MISSING",
    "COMPANION_ARTIFACT_MISSING",
    "ARTIFACT_HASH_MISMATCH",
    "POLARITY_MISMATCH",
    # DatasetState 解析層の引数エラー。missing_state と違い、別のツールを先に
    # 呼んでも直らない（引数を直して呼び直すしかない）。
    "DATASET_BAD_REQUEST",
    # v2 の matrix recipe が schema に合わない。DATASET_BAD_REQUEST と分けるのは、
    # 直す場所が recipe の中だと一目で分かるようにするため。
    "MATRIX_RECIPE_INVALID",
    # キュレーションのフラグ記録（curation/flags.jsonl）に読めない行がある。
    # wrong を黙って落とさないよう、差次的エクスポートを書き出さずに止める。
    "CURATION_FLAGS_INVALID",
    # dataset_feature_search の引数エラー（どれも引数を直して呼び直すしかない）。
    # 目標の形・値が不正（details に target_id と field）。
    "FEATURE_TARGET_INVALID",
    # 極性が決まらず、目標にも polarity・adduct_tiers・(expected_mz と adduct) が無い。
    "POLARITY_UNKNOWN",
    # 窓検索の条件が 1 つも無い（全特徴を黙って列挙しない）。
    "FEATURE_SEARCH_EMPTY_QUERY",
    # 窓の下限が上限より大きい。
    "FEATURE_SEARCH_BAD_RANGE",
    # `.arf` が大きすぎて読めない（伸長後の大きさの見積もりが空きメモリを超える、または読む途中で
    # メモリが尽きた）。details に path_name / estimated_bytes / available_bytes / remedy。
    "ARF_TOO_LARGE",
})


def mztab_error(code: str, message: str, details: dict | None = None) -> str:
    """mzTab-M 処理固有のエラーを機械可読エンベロープで返す。

    code は MZTAB_ERROR_CODES の値を使う。missing_state とは用途が異なる:
    こちらはバリデーション失敗・契約違反のような確定エラーで、
    「前提状態が無い」という復旧可能な状態不足とは意味が異なる。
    """
    if not code or not message:
        raise ValueError("code と message は必須です。")
    payload: dict = {"code": code, "message": message}
    if details is not None:
        payload["details"] = details
    return json_payload({"error": payload})


PATHWAY_ERROR_CODES = frozenset({
    # [mslipidmapper] rscript が未設定、または指す先が無い（mslipidmapper_api_url も無いか使えないとき）。
    "MSLM_NOT_CONFIGURED",
    # MSLipidMapper（R パッケージ）か、その公開関数が無い。自動インストールはしない。
    "MSLM_R_PACKAGES_MISSING",
    # [mslipidmapper] mslipidmapper_api_url のサーバに届かない・生存確認が通らない・API が 404/5xx を返した
    # （ローカルの R も使えなかったときだけ表に出る。使えればそちらへ回る）。
    "MSLM_REMOTE_UNAVAILABLE",
    # Rscript が非ゼロ終了した（stderr 末尾を details に、全文は run の bridge.log）。
    "MSLM_RUN_FAILED",
    # Rscript が制限時間内に終わらなかった（プロセスツリーごと終了済み）。
    "MSLM_TIMEOUT",
    # JSON-LD が無い・壊れている・schemaVersion が想定外（上流の仕様変更の検知）。
    "MSLM_OUTPUT_INVALID",
    # PDF 描画だけが失敗した（ブラウザ不在・unpkg 不達）。pathway_project の結果は有効。
    "PATHWAY_RENDER_FAILED",
    # 引数を直して呼び直すしかない誤り（未知のパスウェイ名・群が小さい・鎖がデータに無い）。
    "PATHWAY_BAD_REQUEST",
})


def pathway_error(code: str, message: str, details: dict | None = None) -> str:
    """パスウェイ投影（MSLipidMapper 連携）固有の確定エラーを返す。"""
    if code not in PATHWAY_ERROR_CODES:
        raise ValueError(f"未知のパスウェイエラーコード: {code!r}")
    if not message:
        raise ValueError("message は必須です。")
    payload: dict = {"code": code, "message": message}
    if details is not None:
        payload["details"] = details
    return json_payload({"error": payload})


def config_invalid_as_envelope(fn):
    """記録の置き場所を決めるときの `DomainError` を、ツールの封筒にして返す。

    記録の置き場所（`[paths] reports_dir` と、データセットの解析フォルダ配下 `reports/`）を
    読むツール（レポート・objective・図の保存・パスウェイ）が対象。
    - `CONFIG_INVALID`（設定ファイルが読めない）: そのままの封筒。以前は例外のまま落ちていた。
    - `NO_REPORT_DIR`（設定の置き場所が無く、データセットも無い）: `missing_state("dataset",
      ["load_dataset"])`。
    ほかのコードの DomainError はそのまま送出する。`@mcp.tool` の内側に付ける
    （`functools.wraps` が引数の形と docstring を保つ）。
    """
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except DomainError as exc:
            if exc.code == "NO_REPORT_DIR":
                return missing_state("dataset", ["load_dataset"], exc.message)
            if exc.code != "CONFIG_INVALID":
                raise
            return json_payload({"error": {"code": exc.code, "message": exc.message,
                                           "details": exc.details}})
    return wrapper
