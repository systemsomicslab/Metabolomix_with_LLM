"""`pipeline-request.v2`（メタボロミクス）のstage計画と依存無効化表（spec §6.2）。

**このモジュールがv2 stage列の唯一の正準**。`store`（`record["stages"]`を作る側）・
`engine`（実行順を決める側）・`recovery`（再開時にpendingへ戻す側）は、いずれも
ここの`build_v2`/`invalidated_v2`を呼ぶ。順序の写しを2箇所に置くと、片方だけ直した
ときに「recordに無いstageを計画する」または「計画にあるのに実行されない」で静かに
ずれる——v1はまさにそれを避けるために`engine.build_stages`と`store._stage_ids_for`が
「同じ規則」であることを注記で守っていたが、v2では規則そのものを1箇所へ寄せる。

v1の`store._BASE_STAGE_IDS`は**変更しない**。v2はその隣に別の列として生える
（stage IDの綴りも異なる: v1 `prepare_input`/`upstream` に対し v2
`prepare_inputs`/`execute_console`）。

依存グラフ上は`metabolomix.pipeline.request_v2`（requestのschema定数）と
`metabolomix.core.atomic_io`（DomainError）だけを引く leaf 寄りのモジュール。
`store`/`engine`/`recovery`からimportされる側なので、それらと同じく
**グローバルsessionをimportしない**（`metabolomix.core.session_state` /
`metabolomix.core.mcp_core` / `metabolomix.tools.*`）。
"""
from __future__ import annotations

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline.request_v2 import SCHEMA as REQUEST_SCHEMA_V2

__all__ = [
    "BASE_V2_STAGE_IDS",
    "CHANGE_ENTRY_POINTS",
    "FINAL_V2_STAGE_ID",
    "PER_STATISTIC_HANDLERS",
    "REQUEST_SCHEMA_V2",
    "UPSTREAM_V2_STAGE_ID",
    "build_v2",
    "invalidated_v2",
    "is_v2_request",
    "stage_ids_v2",
    "statistic_id_from_stage_id",
]

#: spec §6.2のv2 stage列のうち、統計に依存しない前半部分。
#: `prepare_inputs → execute_console → validate_outputs → load_dataset →
#: resolve_metadata → load_assay_evidence → resolve_feature_bindings →
#: qc_raw → preprocess → qc_processed`。
BASE_V2_STAGE_IDS = (
    "prepare_inputs", "execute_console", "validate_outputs", "load_dataset",
    "resolve_metadata", "load_assay_evidence", "resolve_feature_bindings",
    "qc_raw", "preprocess", "qc_processed",
)

#: reportは必ず最後（spec §6.2の角括弧展開の後ろ）。
FINAL_V2_STAGE_ID = "report"

#: 統計1件につき1組ずつ生えるstageのhandlerキー。stage_idは`<handler>:<statistic_id>`。
PER_STATISTIC_HANDLERS = ("statistics", "export")

#: v2でConsoleを起動するstage（v1の`upstream`に相当）。再開時の
#: 「上流をやり直すか」の判定はこのstageの状態を見る。
UPSTREAM_V2_STAGE_ID = "execute_console"

#: 変更の種類 → 最初に無効化する基本stage（spec §6.2「metadata変更は
#: resolve_metadata以降、bindingまたはstandard_assays変更は
#: resolve_feature_bindings以降、recipe変更はpreprocess以降…元データQCの対象集合に
#: 影響するrecipe変更はqc_raw以降を無効化する」）。
#:
#: `sample_manifest`と`metadata`はどちらもresolve_metadata起点だが別の入口として
#: 残す——前者は要求の値そのものの変更、後者は「同じパスのシートを書き直した」
#: 内容変更（v1の`recovery._manifest_content_changed`と同じ信号）で、要求の
#: 差分だけを見ても検出できない。
#: `qc_raw_scope`は「元データQCの評価集合を動かすrecipe変更」を呼び出し側が
#: そう判定したときにだけ渡す token——preprocessの変更が常にqc_rawへ波及する
#: わけではないので、既定の`preprocess`とは分けてある。
CHANGE_ENTRY_POINTS = {
    "metadata": "resolve_metadata",
    "sample_manifest": "resolve_metadata",
    "standard_assays": "resolve_feature_bindings",
    "feature_bindings": "resolve_feature_bindings",
    "preprocess": "preprocess",
    "qc_raw_scope": "qc_raw",
}

#: 統計側だけを丸ごと無効化する変更。`target`（exploratory/differential）は
#: 統計集合の整合そのものを変えるため、全統計と同じ扱いにする。
_ALL_STATISTICS_CHANGES = frozenset({"statistics", "target"})

#: 個別統計の変更token（`statistics:<statistic_id>`）。追加・変更・削除のいずれも
#: この形で渡す——削除済み統計の成果物をcurrentのまま残さないため、削除でも
#: 該当stageとreportを無効化する（spec §6.2）。
_STATISTIC_CHANGE_PREFIX = "statistics:"


def _fail_request(message: str, **details) -> None:
    raise DomainError("PIPELINE_REQUEST_INVALID", message, details)


def is_v2_request(request: object) -> bool:
    """要求が`pipeline-request.v2`かを返す（schema省略は従来どおりv1）。"""
    return isinstance(request, dict) and request.get("schema") == REQUEST_SCHEMA_V2


def statistic_id_from_stage_id(stage_id: str) -> str | None:
    """`statistics:<id>` / `export:<id>` から`<id>`を取り出す（他はNone）。

    v1の`export:<comparison_id>`と綴りが重なるため、**v2の要求だと分かっている
    文脈でだけ**呼ぶこと（`engine.make_context`がschemaで振り分ける）。
    """
    for handler in PER_STATISTIC_HANDLERS:
        prefix = f"{handler}:"
        if stage_id.startswith(prefix):
            return stage_id[len(prefix):]
    return None


def _statistic_ids(request: dict) -> list[str]:
    """要求のstatisticsから`statistic_id`を順序どおりに取り出す。

    空・重複はここで拒否する。重複を通すと`record["stages"]`のdictキーとして
    潰れ、stageが1つ黙って消える（計画にはあるのに実行されない）。
    値そのものの検証（kind・閾値・feature_scope）は`request_v2`の責務で、
    ここでは再実装しない。
    """
    statistics = request.get("statistics") if isinstance(request, dict) else None
    if not isinstance(statistics, list) or not statistics:
        _fail_request("statisticsが空です（v2は非空の統計定義を必須とします）。",
                      statistics=statistics)
    ids: list[str] = []
    for item in statistics:
        if not isinstance(item, dict) or not item.get("statistic_id"):
            _fail_request("statisticsの各要素はstatistic_idを持つ必要があります。",
                          item=item)
        ids.append(item["statistic_id"])
    duplicated = sorted({sid for sid in ids if ids.count(sid) > 1})
    if duplicated:
        _fail_request(f"statistic_idが重複しています: {duplicated}",
                      duplicated=duplicated)
    return ids


def _statistic_stage_ids(statistic_id: str) -> list[str]:
    return [f"{handler}:{statistic_id}" for handler in PER_STATISTIC_HANDLERS]


def build_v2(request: dict) -> list[dict]:
    """v2のstage計画を組み立てる（spec §6.2の順序そのまま）。

    各要素は`{"stage_id", "handler", "statistic_id"}`。`statistics:<id>` /
    `export:<id>` は統計1件につき1組ずつ、要求に現れた順で展開し、`report`は
    必ず最後に置く。基本stageの`handler`はstage_idと同じ、`statistic_id`はNone。
    """
    statistic_ids = _statistic_ids(request)

    plan = [{"stage_id": stage_id, "handler": stage_id, "statistic_id": None}
            for stage_id in BASE_V2_STAGE_IDS]
    for statistic_id in statistic_ids:
        for handler in PER_STATISTIC_HANDLERS:
            plan.append({"stage_id": f"{handler}:{statistic_id}",
                         "handler": handler, "statistic_id": statistic_id})
    plan.append({"stage_id": FINAL_V2_STAGE_ID, "handler": FINAL_V2_STAGE_ID,
                 "statistic_id": None})
    return plan


def stage_ids_v2(request: dict) -> list[str]:
    """`build_v2`の計画順のstage ID列（`store`が`record["stages"]`を作るのに使う）。"""
    return [entry["stage_id"] for entry in build_v2(request)]


def _downstream_of(entry_stage_id: str, statistic_ids: list[str]) -> set[str]:
    """基本stage`entry_stage_id`から下流すべて（全統計とreportを含む）。"""
    index = BASE_V2_STAGE_IDS.index(entry_stage_id)
    dirty = set(BASE_V2_STAGE_IDS[index:])
    for statistic_id in statistic_ids:
        dirty.update(_statistic_stage_ids(statistic_id))
    dirty.add(FINAL_V2_STAGE_ID)
    return dirty


def invalidated_v2(changed: set[str], request: dict) -> set[str]:
    """変更内容から、やり直しが要るstage IDの集合を返す（spec §6.2）。

    `changed`は「何が変わったか」のtoken集合:

    * `CHANGE_ENTRY_POINTS`のキー（`metadata` / `sample_manifest` /
      `standard_assays` / `feature_bindings` / `preprocess` / `qc_raw_scope`）
      ——その入口stageから下流すべて。上流（`prepare_inputs` /
      `execute_console` / `validate_outputs`）はどの入口からも無効化されない
      ——binding訂正でConsoleを再実行しない（spec §6.2）。上流のやり直しは
      `rerun_upstream`の明示だけが起点で、それは呼び出し側が扱う。
    * `statistics` / `target` ——全統計stageとreport。
    * `statistics:<statistic_id>` ——その統計のstatistics/exportとreport。
      追加・変更・削除のいずれもこの形で渡す。削除された統計のstage IDは
      新しい計画には無いが、`record["stages"]`には残っているため戻り値には
      含める（呼び出し側がrecordのキーと突き合わせる）。

    未知のtokenは拒否する。黙って無視すると「無効化したつもりで何もしていない」
    ——古い結果をcurrentのまま残す、最も見つけにくい壊れ方になる。
    """
    statistic_ids = _statistic_ids(request)
    all_statistic_stages = {stage_id for sid in statistic_ids
                            for stage_id in _statistic_stage_ids(sid)}

    dirty: set[str] = set()
    for token in changed:
        if token in CHANGE_ENTRY_POINTS:
            dirty |= _downstream_of(CHANGE_ENTRY_POINTS[token], statistic_ids)
        elif token in _ALL_STATISTICS_CHANGES:
            dirty |= all_statistic_stages
            dirty.add(FINAL_V2_STAGE_ID)
        elif token.startswith(_STATISTIC_CHANGE_PREFIX):
            statistic_id = token[len(_STATISTIC_CHANGE_PREFIX):]
            if not statistic_id:
                raise DomainError(
                    "PIPELINE_STAGE_PLAN_INVALID",
                    f"statistic_idの無い変更tokenです: {token!r}", {"token": token})
            dirty.update(_statistic_stage_ids(statistic_id))
            dirty.add(FINAL_V2_STAGE_ID)
        else:
            raise DomainError(
                "PIPELINE_STAGE_PLAN_INVALID",
                f"未知の変更種別です: {token!r}"
                f"（既知: {sorted(CHANGE_ENTRY_POINTS)} / {sorted(_ALL_STATISTICS_CHANGES)}"
                f" / '{_STATISTIC_CHANGE_PREFIX}<statistic_id>'）",
                {"token": token})
    return dirty
