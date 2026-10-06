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

このモジュールは依存グラフの leaf（stdlib と同じく leaf の
metabolomix.core.serialization のみ）。tools_* / server を import しない。
"""
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
    # MSDIAL_EXE が Console ではなく GUI を指している場合。
    "MSDIAL_EXE_NOT_CONSOLE",
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
