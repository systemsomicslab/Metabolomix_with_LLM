"""pipeline-request.v1 の厳密な検証・解決・更新契約（spec §10.1）。

このモジュールは以降のTask（13〜18）がrequestを読み書きする唯一の入口になる。
守るべき区別は次の2つ:

1. **未指定 と 明示null の区別**。例えば ``sample_manifest`` を省略すると
   「既定シート名を探索する」(既定探索)ことを意味し、明示的に ``null`` を
   渡すと「探索せず自動一覧生成へ切り替える」(明示解除)ことを意味する——
   両者ともPython上の値は``None``だが、``value_sources`` の出所
   （``"default"`` か ``"explicit"``/``"explicit_update"``）で区別する。
   同様に ``preprocess.blank_min_fold`` / ``max_qc_rsd`` は明示 ``null`` で
   無効化でき、``normalize=none`` ``drift_correct=false``
   ``min_detection_rate=0.0`` も明示無効化として ``auto`` と区別される
   （値そのものは既定と一致しうるので、ここでも出所が唯一の手がかり）。

2. **要求済みの上流条件 と resumeで訂正可能な下流条件の区別**。
   ``UPDATABLE`` に挙げた4キー（target・sample_manifest・preprocess・
   comparisons）だけが ``merge_updates`` で変更でき、それ以外
   （method_file・lbm_file・polarity・measure・keep_extension・timeout_s・
   save_project・output_root、および内部専用キーeffective_target・
   value_sources）を変えようとすると ``NEW_PIPELINE_REQUIRED`` になる——
   上流の実行条件が変わるなら別pipelineを作るべきで、resumeを
   「同じ入力を繰り返し送るだけの無限ループ」にしないための境界線。

解決順序は「明示値 > (前revisionを読み戻した)既存値 > 既定値」の3段。
``resolve_request`` が最初の revision を「既定値」から作り、
``merge_updates`` がその戻り値（＝保存され読み戻された「既存値」）へ
新しい明示値を重ねる、という形でこの3段が表現される。

初回解決（``resolve_request``）の「既定値」の手前には、もう1段
**元フォルダ直下の ``analysis-request.json``**（spec §7.1「要求値の優先順位は
MCPで明示した値 > analysis-request.json > 本specの既定」）が入る。ファイルの
値は ``value_sources`` で ``"request_file"`` として区別され、MCPで明示した値が
あればそちらが勝つ。ファイルの不正（JSONとして壊れている・オブジェクトでない・
未知キー・不許可のnull）は明示値と同じ土俵で ``PIPELINE_REQUEST_INVALID``
にする——「置いてあるのに黙って無視された」という一番わかりにくい失敗を
作らない。

不正JSON・未知キー・不正値はこの場でDomainError（コード
``PIPELINE_REQUEST_INVALID``）にする。上流情報が単に「まだ無い」ケースを
``missing_state`` 風の封筒に化けさせて同じ呼び出しを無限反復させるのは
後続タスク（pipeline_plan/pipeline_run, Task 13/14）の責務であり、ここでは
やらない。
"""
from __future__ import annotations

import copy
import json
import math
import re
from pathlib import Path

from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.pipeline import request_v2

__all__ = [
    "REQUEST_FILE_NAME",
    "SCHEMA",
    "UPDATABLE",
    "read_request_file",
    "merge_updates",
    "request_fingerprint",
    "resolve_request",
    "validate_request",
]

SCHEMA = "pipeline-request.v1"

#: 元フォルダ直下から探す任意の要求ファイル（spec §7.1）。存在すれば
#: 「MCPで明示した値 > analysis-request.json > 本specの既定」の中段になる。
REQUEST_FILE_NAME = "analysis-request.json"

#: resumeで変更可能なトップレベルキー（spec §10.1: 「resumeで変更可能なのは
#: target、sample_manifest、preprocess、comparisonsに限定する」）。
UPDATABLE = {"target", "sample_manifest", "preprocess", "comparisons"}

#: request/updatesが受け付けるトップレベルキー（spec §10.1本文の列挙そのまま）。
_TOP_LEVEL_KEYS = frozenset({
    "schema", "target", "method_file", "lbm_file", "polarity", "measure",
    "keep_extension", "timeout_s", "save_project", "output_root",
    "sample_manifest", "preprocess", "comparisons",
})

#: 内部でだけ使う2キー。外部入力（explicit/updates）には現れてはいけない。
_INTERNAL_KEYS = frozenset({"effective_target", "value_sources"})

_TARGET_VALUES = frozenset({"auto", "exploratory", "differential"})
_POLARITY_VALUES = frozenset({"positive", "negative"})
_MEASURE_VALUES = frozenset({"peak_height"})  # 初期版はpeak_heightのみ(spec §10.1)
_NORMALIZE_VALUES = frozenset({"auto", "none", "tic", "median", "pqn"})
_IMPUTE_VALUES = frozenset({"none", "half_min", "knn", "column_mean"})
#: 自動前処理は現行版ではconservative-v1のみ（common-context.md）。
_POLICY_VALUES = frozenset({"conservative-v1"})

_PREPROCESS_KEYS = frozenset({
    "policy", "normalize", "blank_min_fold", "drift_correct", "max_qc_rsd",
    "impute", "min_detection_rate",
})

_COMPARISON_KEYS = frozenset({
    "comparison_id", "reference_group", "test_group", "q_threshold",
    "log2fc_threshold", "log_transform", "allow_confounded",
})
_COMPARISON_REQUIRED_KEYS = frozenset({"comparison_id", "reference_group", "test_group"})

#: comparison_idは結果ディレクトリ名やstage id(`export:<comparison_id>`等)に
#: 使うため、パス区切り・"."を含む脱出パターンを一切許さない安全な文字だけに絞る。
_SAFE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")

_DEFAULT_PREPROCESS = {
    "policy": "conservative-v1", "normalize": "auto", "blank_min_fold": "auto",
    "drift_correct": "auto", "max_qc_rsd": "auto", "impute": "half_min",
    "min_detection_rate": 0.0,
}

#: 個別comparisonの未指定時デフォルト(comparison_id/reference_group/test_groupは
#: 必須で既定を持たない)。既存dataset_analysis_tools.dataset_differentialの
#: 既定(q_threshold=0.05, log2fc_threshold=1.0, log_transform=True)と揃える。
_COMPARISON_DEFAULTS = {
    "q_threshold": 0.05, "log2fc_threshold": 1.0,
    "log_transform": True, "allow_confounded": False,
}

_DEFAULT_TIMEOUT_S = 21600  # 既存Consoleの既定（common-context.md）に合わせる

#: 明示nullで「無効化」を表せるトップレベルの省略可フィールドはこの1つだけ
#: （sample_manifest: 明示シートを解除し自動一覧生成へ切り替える）。
#: preprocess.blank_min_fold/max_qc_rsdは``_validate_preprocess``側で別途許可する。
#: （controller裁定R13、spec l.246「nullが無効化を意味する設定だけで許可する」/
#: l.421「他の不許可なnullはエラーにする」）。
_NULL_ALLOWED_TOP_LEVEL_KEYS = frozenset({"sample_manifest"})

#: 省略（未指定）は既定値扱いだが、明示nullを渡すとエラーになるトップレベル
#: フィールド。method_file等は「未確定」を``None``で表すが、それは省略時の
#: 既定値であって、明示nullを同じ意味で受け付けてよい理由にはならない
#: ——nullは「無効化」を意味する設定だけの特別な語彙（上のR13）。
_NULL_REJECTED_TOP_LEVEL_KEYS = frozenset({
    "method_file", "lbm_file", "polarity", "keep_extension", "output_root",
    "comparisons",
})


def _fail(message: str, **details) -> None:
    raise DomainError("PIPELINE_REQUEST_INVALID", message, details)


def _reject_disallowed_explicit_null(source: dict, keys) -> None:
    """``keys``のうち``source``に存在し、値が明示的に``None``のものを拒否する。

    キー自体が``source``に無い場合（＝未指定）は素通りする——未指定と明示null
    の区別（spec l.246/421、controller裁定R13）そのものがこの関数の核心。
    そのため呼び出し側は「キーが存在したかどうか」を保っている生の入力
    （``resolve_request``の``explicit``辞書、``merge_updates``の``updates``
    辞書）を渡す必要がある。既に全キーがそろっている検証済みdataを渡しても、
    そこでは省略も明示nullも同じ``None``に潰れていて区別できない。
    """
    for key in keys:
        if key in source and source[key] is None:
            _fail(f"{key}に明示nullは許可されません（無効化を意味しない設定のため、"
                  f"未指定にしてください）。", **{key: None})


def _is_finite_number(value: object) -> bool:
    """bool・非数値・NaN/Infinityを除いた「本物の有限数値」だけを真にする。

    Pythonの``bool``は``int``のサブクラスなので、素朴な``isinstance(v, (int, float))``
    だけでは``True``/``False``がしきい値として通ってしまう——timeout_sや
    各種閾値でboolを拒否するという契約の核心はここにある。
    """
    if isinstance(value, bool):
        return False
    if not isinstance(value, (int, float)):
        return False
    return math.isfinite(float(value))


# ---------- preprocess ----------

def _validate_preprocess(preprocess: object) -> dict:
    """conservative-v1の既定へ指定キーだけを重ね、全キーを検証して返す。

    ``preprocess``がNoneなら「何も指定されていない」ことを表す（resolve_requestの
    初回解決、またはmerge_updatesの対象外フィールド）。dict以外・未知キーは
    ここで拒否する。返り値は必ず7キー全部そろった新しいdict（呼び出し側の
    dictを書き換えない）。
    """
    if preprocess is None:
        preprocess = {}
    if not isinstance(preprocess, dict):
        _fail("preprocessはオブジェクトである必要があります。", value=preprocess)
    unknown = set(preprocess) - _PREPROCESS_KEYS
    if unknown:
        _fail(f"preprocessに未知のキーがあります: {sorted(unknown)}",
              unknown_keys=sorted(unknown))

    merged = dict(_DEFAULT_PREPROCESS)
    merged.update(preprocess)

    if merged["policy"] not in _POLICY_VALUES:
        _fail(f"preprocess.policyが不正です: {merged['policy']!r}",
              policy=merged["policy"])
    if merged["normalize"] not in _NORMALIZE_VALUES:
        _fail(f"preprocess.normalizeが不正です: {merged['normalize']!r}",
              normalize=merged["normalize"])
    if merged["impute"] not in _IMPUTE_VALUES:
        _fail(f"preprocess.imputeが不正です: {merged['impute']!r}",
              impute=merged["impute"])

    drift_correct = merged["drift_correct"]
    if not (drift_correct == "auto" or isinstance(drift_correct, bool)):
        _fail(f"preprocess.drift_correctが不正です: {drift_correct!r}",
              drift_correct=drift_correct)

    # blank_min_fold・max_qc_rsdだけは明示nullで無効化できる（自動閾値の停止）。
    # "auto"文字列・None・有限の正数、以外はすべて拒否する。
    for field in ("blank_min_fold", "max_qc_rsd"):
        value = merged[field]
        if not (value == "auto" or value is None
                or (_is_finite_number(value) and value > 0)):
            _fail(f"preprocess.{field}が不正です: {value!r}", **{field: value})

    min_detection_rate = merged["min_detection_rate"]
    if not (_is_finite_number(min_detection_rate)
            and 0.0 <= float(min_detection_rate) <= 1.0):
        _fail(f"preprocess.min_detection_rateが不正です: {min_detection_rate!r}",
              min_detection_rate=min_detection_rate)

    return merged


# ---------- comparisons ----------

def _validate_comparisons(comparisons: object) -> list:
    """配列全体を検証し、各要素へ未指定フィールドの既定値を埋めて返す。

    ``comparisons``は「置換」対象（部分更新はない）なので、ここでは常に
    配列全体を一から検証・正規化する。comparison_id重複・パスとして安全でない
    ID・reference_group==test_groupは1件でもあれば例外にする。
    """
    if comparisons is None:
        comparisons = []
    if not isinstance(comparisons, list):
        _fail("comparisonsは配列である必要があります。", value=comparisons)

    normalized: list = []
    seen_ids: set = set()
    for index, item in enumerate(comparisons):
        if not isinstance(item, dict):
            _fail(f"comparisons[{index}]はオブジェクトである必要があります。", index=index)
        unknown = set(item) - _COMPARISON_KEYS
        if unknown:
            _fail(f"comparisons[{index}]に未知のキーがあります: {sorted(unknown)}",
                  index=index, unknown_keys=sorted(unknown))
        missing = _COMPARISON_REQUIRED_KEYS - set(item)
        if missing:
            _fail(f"comparisons[{index}]に必須キーが不足しています: {sorted(missing)}",
                  index=index, missing_keys=sorted(missing))

        comparison_id = item["comparison_id"]
        if not isinstance(comparison_id, str) or not _SAFE_ID_RE.fullmatch(comparison_id):
            _fail(
                f"comparison_idが不正です（安全な文字だけのIDが必要・パス脱出不可）: "
                f"{comparison_id!r}",
                index=index, comparison_id=comparison_id,
            )
        if comparison_id in seen_ids:
            _fail(f"comparison_idが重複しています: {comparison_id!r}",
                  comparison_id=comparison_id)
        seen_ids.add(comparison_id)

        reference_group = item["reference_group"]
        test_group = item["test_group"]
        if not isinstance(reference_group, str) or not reference_group:
            _fail(f"reference_groupが不正です: {reference_group!r}", index=index)
        if not isinstance(test_group, str) or not test_group:
            _fail(f"test_groupが不正です: {test_group!r}", index=index)
        if reference_group == test_group:
            _fail(
                f"reference_groupとtest_groupが同一です（比較の向きを明示できません）: "
                f"{reference_group!r}",
                index=index, group=reference_group,
            )

        q_threshold = item.get("q_threshold", _COMPARISON_DEFAULTS["q_threshold"])
        if not (_is_finite_number(q_threshold) and 0 < float(q_threshold) <= 1):
            _fail(f"q_thresholdが不正です: {q_threshold!r}",
                  index=index, q_threshold=q_threshold)

        log2fc_threshold = item.get("log2fc_threshold", _COMPARISON_DEFAULTS["log2fc_threshold"])
        if not (_is_finite_number(log2fc_threshold) and float(log2fc_threshold) >= 0):
            _fail(f"log2fc_thresholdが不正です: {log2fc_threshold!r}",
                  index=index, log2fc_threshold=log2fc_threshold)

        log_transform = item.get("log_transform", _COMPARISON_DEFAULTS["log_transform"])
        if not isinstance(log_transform, bool):
            _fail(f"log_transformはboolのみ許可されます: {log_transform!r}", index=index)

        allow_confounded = item.get("allow_confounded", _COMPARISON_DEFAULTS["allow_confounded"])
        if not isinstance(allow_confounded, bool):
            _fail(f"allow_confoundedはboolのみ許可されます: {allow_confounded!r}", index=index)

        normalized.append({
            "comparison_id": comparison_id,
            "reference_group": reference_group,
            "test_group": test_group,
            "q_threshold": float(q_threshold),
            "log2fc_threshold": float(log2fc_threshold),
            "log_transform": log_transform,
            "allow_confounded": allow_confounded,
        })
    return normalized


# ---------- トップレベルフィールド ----------

def _validate_optional_str(data: dict, field: str) -> None:
    value = data.get(field)
    if value is not None and (not isinstance(value, str) or value == ""):
        _fail(f"{field}は空でない文字列またはnullである必要があります: {value!r}",
              **{field: value})


def _validate_fields(data: dict) -> None:
    """dataが13キー全部そろった形である前提で、型・enum・数値域を全検査する。

    ここで``data["preprocess"]``/``data["comparisons"]``を正規化済みの値へ
    書き換える（呼び出し側の意図した「返り値は常に完全な形」を満たすため）。
    """
    if data.get("schema") != SCHEMA:
        _fail(f"schemaは{SCHEMA!r}のみ許可されます。", schema=data.get("schema"))

    target = data.get("target")
    if target not in _TARGET_VALUES:
        _fail(f"targetが不正です: {target!r}", target=target)

    for field in ("method_file", "lbm_file", "output_root", "sample_manifest"):
        _validate_optional_str(data, field)

    polarity = data.get("polarity")
    if polarity is not None and polarity not in _POLARITY_VALUES:
        _fail(f"polarityが不正です: {polarity!r}", polarity=polarity)

    measure = data.get("measure")
    if measure not in _MEASURE_VALUES:
        _fail(f"measureが不正です: {measure!r}", measure=measure)

    _validate_optional_str(data, "keep_extension")

    timeout_s = data.get("timeout_s")
    if isinstance(timeout_s, bool) or not isinstance(timeout_s, int) or timeout_s <= 0:
        _fail(f"timeout_sは正の整数である必要があります（boolは不可）: {timeout_s!r}",
              timeout_s=timeout_s)

    save_project = data.get("save_project")
    if not isinstance(save_project, bool):
        _fail(f"save_projectはboolのみ許可されます: {save_project!r}",
              save_project=save_project)

    data["preprocess"] = _validate_preprocess(data.get("preprocess"))
    data["comparisons"] = _validate_comparisons(data.get("comparisons"))

    if target == "exploratory" and data["comparisons"]:
        # 明示的なexploratoryと比較指定は両立しない矛盾した要求。黙って
        # 比較を落とすと`engine.build_stages`が`differential:*`/`export:*`を
        # 計画から外し、`_mark_out_of_scope`はstage側のwarningsにしか
        # `STAGE_OUT_OF_SCOPE_FOR_TARGET`を書かない——`read_status`も品質
        # レポートの「Unverified Conditions」も`record["warnings"]`しか
        # 読まないため、利用者が頼んだ比較が1件も実行されないまま
        # `completed`が返り、痕跡がどこにも残らない。spec §8.2「明示的に
        # 要求した処理を実施できなければneeds_input」・§9.2「既知の不正入力は
        # Console起動前に拒否する」に従い、受付で止める。
        # `target="auto"`＋比較は矛盾ではない（`differential`へ昇格する）。
        _fail("target='exploratory'とcomparisonsは同時に指定できません"
              "（比較を実行するならtarget='differential'または'auto'にしてください）。",
              target=target,
              comparison_ids=[c["comparison_id"] for c in data["comparisons"]])


def _compute_effective_target(target: str, comparisons: list) -> str:
    """spec §9: target=autoは、有効な比較定義があればdifferential、なければ
    exploratoryへ固定する。target自体が明示されていればそれを採る。"""
    if target != "auto":
        return target
    return "differential" if comparisons else "exploratory"


def _promote_value_sources(value_sources: dict, updated_fields: dict) -> None:
    """updated_fieldsに挙げたフィールドの出所をexplicit_updateへ書き換える。

    preprocessだけは子キー単位（updates["preprocess"]に実際に含まれていた
    キーだけ）を更新し、挙げられなかった子キー・トップレベルフィールドの
    出所には一切触れない。
    """
    for key, value in updated_fields.items():
        if key == "preprocess":
            sub_sources = value_sources.setdefault("preprocess", {})
            for sub_key in value:
                sub_sources[sub_key] = "explicit_update"
        else:
            value_sources[key] = "explicit_update"


# ---------- 公開API ----------

def validate_request(data: dict, *, internal: bool = False,
                      updated_fields: dict | None = None) -> dict:
    """要求の形・値を検証し、``effective_target``を(再)計算して返す。

    ``internal=False``（既定）は外部入力の検査専用モード。``data``に
    ``effective_target``/``value_sources``という内部専用キーが含まれていたら
    （＝未知キーとして）拒否する。

    ``internal=True``は「検証済みの要求を更新したあとの再検証」専用。
    ``data``は既に``value_sources``/``effective_target``を持つ完全な要求で
    ある前提とし、``effective_target``を最新のtarget/comparisonsから
    再計算し、``updated_fields``に挙げたフィールドの``value_sources``を
    ``"explicit_update"``へ書き換える（挙げられなかったフィールドの出所は
    保持する）。
    """
    allowed = _TOP_LEVEL_KEYS | _INTERNAL_KEYS if internal else _TOP_LEVEL_KEYS
    unknown = set(data) - allowed
    if unknown:
        _fail(f"未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))

    if internal:
        missing = {"value_sources", "effective_target"} - set(data)
        if missing:
            _fail(f"internal=Trueの検証にはvalue_sources/effective_targetが必要です"
                  f"（不足: {sorted(missing)}）。", missing=sorted(missing))

    _validate_fields(data)

    data["effective_target"] = _compute_effective_target(data["target"], data["comparisons"])

    if internal and updated_fields:
        _promote_value_sources(data["value_sources"], updated_fields)

    return data


def read_request_file(source_root: Path) -> dict:
    """元フォルダ直下の``analysis-request.json``を読む（無ければ空dict、spec §7.1）。

    存在するのに読めない・JSONとして壊れている・オブジェクトでない・未知キーを
    持つ場合は``PIPELINE_REQUEST_INVALID``。「置いてあるのに黙って無視された」
    が一番危ない失敗——利用者は自分が書いた条件で走ったと思い込む。

    値そのものの検証（型・範囲・列挙値）はここではしない。明示値と重ねた
    あとに``validate_request``が1か所で行う。
    """
    path = Path(source_root) / REQUEST_FILE_NAME
    if not path.is_file():
        return {}
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        _fail(f"{REQUEST_FILE_NAME}を読めません: {exc}", path=str(path))
    try:
        data = json.loads(raw)
    except ValueError as exc:
        _fail(f"{REQUEST_FILE_NAME}がJSONとして壊れています: {exc}", path=str(path))
    if not isinstance(data, dict):
        _fail(f"{REQUEST_FILE_NAME}はオブジェクトである必要があります。",
              path=str(path), value=data)
    unknown = set(data) - _TOP_LEVEL_KEYS
    if unknown:
        _fail(f"{REQUEST_FILE_NAME}に未知のキーがあります: {sorted(unknown)}",
              path=str(path), unknown_keys=sorted(unknown))
    _reject_disallowed_explicit_null(data, _NULL_REJECTED_TOP_LEVEL_KEYS)
    return data


def _layer_sources(explicit: dict, from_file: dict) -> dict:
    """トップレベル各キーの出所（explicit/request_file/default）を返す。"""
    sources = {}
    for key in _TOP_LEVEL_KEYS:
        if key in explicit:
            sources[key] = "explicit"
        elif key in from_file:
            sources[key] = "request_file"
        else:
            sources[key] = "default"
    return sources


def _layer_preprocess(explicit: dict, from_file: dict) -> tuple:
    """preprocessを子キー単位で重ね、(値, 子キーごとの出所)を返す。

    どちらか一方でもdictでなければ重ねずに優先順位だけで選ぶ——不正な型は
    ``_validate_preprocess``が1か所で拒否する。
    """
    explicit_pp = explicit.get("preprocess")
    file_pp = from_file.get("preprocess")
    if isinstance(explicit_pp, dict) and isinstance(file_pp, dict):
        merged = {**file_pp, **explicit_pp}
    elif "preprocess" in explicit:
        merged = explicit_pp
    else:
        merged = file_pp

    explicit_keys = explicit_pp if isinstance(explicit_pp, dict) else {}
    file_keys = file_pp if isinstance(file_pp, dict) else {}
    sources = {}
    for key in _PREPROCESS_KEYS:
        if key in explicit_keys:
            sources[key] = "explicit"
        elif key in file_keys:
            sources[key] = "request_file"
        else:
            sources[key] = "default"
    return merged, sources


def _peek_schema_declared_by_request_file(source_root: Path) -> object:
    """``analysis-request.json``の``schema``フィールドだけを覗く。

    v1/v2どちらのkey setで本読込・検証すべきかを、内容を検証する前に決める
    ための最小限の先読み——本読み込み（v1の``read_request_file``、または
    v2の``request_v2.read_request_file``）が正しいkey setで担当する前に、
    このファイル自体をv1のkey setで一度でも検証してしまうと、v2形状の
    ファイルが「未知のキー」として誤って拒否される。

    ファイルが無い・JSONとして壊れている・オブジェクトでない場合は例外に
    せず``None``を返す——実際のエラー報告は、schemaが確定した後にどのみち
    同じファイルを読み直す本読み込み関数（v1なら``read_request_file``、v2
    なら``request_v2.read_request_file``）に一本化する。ここで二重にエラーを
    出さない。
    """
    path = Path(source_root) / REQUEST_FILE_NAME
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict):
        return None
    return data.get("schema")


def resolve_request(source_root: Path, explicit: dict | None = None, *,
                     profile: dict | None = None,
                     routine_overrides: dict | None = None) -> dict:
    """要求を解決し、``value_sources``/``effective_target``を付けて返す。

    優先順位は「MCPで明示した値 > 元フォルダ直下の``analysis-request.json``
    > 本specの既定」（spec §7.1）。``preprocess``だけは子キー単位で重ねる
    ——ファイルで``normalize``を決め、MCPでは``min_detection_rate``だけを
    上書きする、という現実的な使い方で、片方が丸ごと消える形にしない。

    ``source_root``は``analysis-request.json``の探索基準であり、Task 13/14以降が
    相対パス（method_file・sample_manifest等）を実解決する基準として持ち回る
    値でもある——このtask自体はパス脱出検証をsource_root基準では行わない
    （comparison_idを除く各パス系フィールドは、相対と絶対のどちらも許す。
    spec §10.1の例で``method_file``にラボ共有フォルダの絶対パスが使われている
    通り、source_rootの外を指す正当なケースがあるため）。

    ``explicit["schema"]``（省略時は``analysis-request.json``自身が宣言する
    ``schema``、それも無ければv1）が``"pipeline-request.v2"``のときは
    ``metabolomix.pipeline.request_v2.resolve``へ丸ごと委譲する（spec §6の
    schema dispatch）。**explicit/file双方がschemaを省略したときだけがv1**
    ——ファイル自身がv2を宣言している場合はそれ自体が明示的な指定であり、
    「省略」ではない（`docs/schema/pipeline-request-v2.md`「schema dispatch」
    節）。v2解決には検証済みprofile dict（``profile``引数、
    ``metabolomix.console.profile_schema.validate_profile``の戻り値相当）が要る
    ——読み込み自体はこの関数の責務ではない。``routine_overrides``（省略可）は
    ``metabolomix.console.profile_schema.validate_certificate``の戻り値の同名
    キー——``execution_purpose="routine"``のpreprocess上書きの明示許容集合。
    省略時（``None``）はfail-closed（routineでは何も上書きできない）。
    """
    source_root = Path(source_root)
    explicit = explicit if explicit is not None else {}
    if not isinstance(explicit, dict):
        _fail("requestはオブジェクトである必要があります。", value=explicit)

    schema = explicit["schema"] if "schema" in explicit \
        else _peek_schema_declared_by_request_file(source_root)
    if schema == request_v2.SCHEMA:
        from_file = request_v2.read_request_file(source_root)
        resolved = request_v2.resolve(
            explicit, profile, from_file=from_file, routine_overrides=routine_overrides)
        # workerのcwdに依存せず、受付で読んだ同じprofileを開けるようにする。
        profile_path = Path(resolved["profile_file"]).expanduser()
        if not profile_path.is_absolute():
            profile_path = source_root / profile_path
        resolved["profile_file"] = str(profile_path.resolve())
        return resolved

    unknown = set(explicit) - _TOP_LEVEL_KEYS
    if unknown:
        _fail(f"未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))

    # 未指定(キー無し)と明示null(値がNone)の区別はここでしか付かない
    # ——このあとdataを組み立てる時点で両者は同じNoneへ潰れる(R13)。
    _reject_disallowed_explicit_null(explicit, _NULL_REJECTED_TOP_LEVEL_KEYS)

    from_file = read_request_file(source_root)

    def _pick(key: str, default=None):
        """明示値 > ファイル値 > 既定値。キーの有無で選ぶ（値がNoneでも明示は明示）。"""
        if key in explicit:
            return explicit[key]
        if key in from_file:
            return from_file[key]
        return default

    preprocess, preprocess_sources = _layer_preprocess(explicit, from_file)
    data = {
        "schema": _pick("schema", SCHEMA),
        "target": _pick("target", "auto"),
        "method_file": _pick("method_file"),
        "lbm_file": _pick("lbm_file"),
        "polarity": _pick("polarity"),
        "measure": _pick("measure", "peak_height"),
        "keep_extension": _pick("keep_extension"),
        "timeout_s": _pick("timeout_s", _DEFAULT_TIMEOUT_S),
        "save_project": _pick("save_project", True),
        "output_root": _pick("output_root"),
        "sample_manifest": _pick("sample_manifest"),
        "preprocess": preprocess,
        "comparisons": _pick("comparisons"),
    }

    value_sources = _layer_sources(explicit, from_file)
    value_sources["preprocess"] = preprocess_sources

    validated = validate_request(data, internal=False)
    validated["value_sources"] = value_sources
    return validated


def merge_updates(request: dict, updates: dict, *, profile: dict | None = None,
                   routine_overrides: dict | None = None) -> dict:
    """resumeの入力訂正を反映し、再検証した要求を返す（spec §9/§10.1）。

    ``UPDATABLE``（target・sample_manifest・preprocess・comparisons）以外の
    キーを変えようとした場合は、上流の実行条件そのものの変更とみなし
    ``NEW_PIPELINE_REQUIRED``にする——effective_target/value_sourcesという
    内部キーも``UPDATABLE``に含まれないため、ここで同じ扱いになる
    （updates自身に内部キーを許可しないという契約を、UPDATABLEの外側として
    自然に満たす）。

    ``request["schema"] == "pipeline-request.v2"``のときは
    ``metabolomix.pipeline.request_v2.merge_updates``へ丸ごと委譲する（spec §6の
    schema dispatch）。v2解決同様、検証済みprofile dict・routine_overrides
    （証明書由来の明示許容集合）は呼び出し側が渡す。
    """
    if request.get("schema") == request_v2.SCHEMA:
        return request_v2.merge_updates(
            request, updates, profile, routine_overrides=routine_overrides)

    if set(updates) - UPDATABLE:
        raise DomainError("NEW_PIPELINE_REQUIRED", "上流条件の変更には新しい解析が必要です")
    # UPDATABLE個の中で明示nullが不許可なのはcomparisonsだけ(R13)——
    # sample_manifestは無効化として許可、preprocessは子キー単位で別途扱う。
    _reject_disallowed_explicit_null(updates, _NULL_REJECTED_TOP_LEVEL_KEYS & UPDATABLE)
    out = copy.deepcopy(request)
    for key, value in updates.items():
        if key == "preprocess":
            # out[key].update(value) より前に、更新値がdictであること・
            # そのキーが既知の子キーだけであることを確認する
            # （検証前にupdate()してしまうと不正な子キーが紛れ込む）。
            if not isinstance(value, dict):
                _fail("preprocessの更新値はオブジェクトである必要があります。", value=value)
            unknown = set(value) - _PREPROCESS_KEYS
            if unknown:
                _fail(f"preprocessに未知の更新キーがあります: {sorted(unknown)}",
                      unknown_keys=sorted(unknown))
            out[key].update(value)
        else:
            out[key] = copy.deepcopy(value)
    return validate_request(out, internal=True, updated_fields=updates)


def request_fingerprint(request: dict) -> str:
    """要求の内容hashを返す（`value_sources`/`effective_target`のような
    由来情報は対象から外す)。

    出所がexplicit/default/explicit_updateのどれであっても、最終的な
    フィールドの値が同じなら同一内容として扱う——Task 14がこのhashで
    run/request_idの再送・冪等性を判定する土台になるため、由来だけの違いで
    別内容と誤判定してはいけない。

    **対象キー集合はschemaごとに違う**（spec §6のschema dispatch）。v2要求を
    v1のキー集合で畳むと`statistics`・`feature_bindings`・`standard_assays`・
    `profile_file`・`omics`・`execution_purpose`が丸ごと落ち、**統計定義だけが
    違う2つの解析が同じ内容hash**になる——`store.find_or_create_run`はそれを
    「同じ要求の再送」として1本のrunへ畳み、`recovery.prepare_resume`は
    統計の訂正を「何も変わっていない」と読んでrevisionを上げない。spec §6.1が
    「下流結果IDは…変換、群・検定設定に依存する」と定めるのはこの逆である。

    v2のキー集合は`request_v2._TOP_LEVEL_KEYS`をそのまま使う——ここへ写しを
    作ると、v2にフィールドが増えたときに片方だけ腐って同じ衝突が戻ってくる
    （`request_v2`の内部定数を参照しているのは、この1箇所がv2キー集合の
    唯一の正準であることを崩さないため）。v1の対象キーと正規化は変えない。
    """
    key_set = (request_v2._TOP_LEVEL_KEYS
               if request.get("schema") == request_v2.SCHEMA else _TOP_LEVEL_KEYS)
    content = {key: request[key] for key in key_set if key in request}
    return canonical_hash(content)
