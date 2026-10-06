"""結果の同一性・来歴・派生結果の無効化（spec §6）。

数値層の純ロジックで、MCP にも session にも依存しない。ここが答えるのは 2 つの問い:

**この結果はどの入力から出たのか。** 図と TSV が別々の前処理から出た数字を並べても、
どちらも「直近の結果」を名乗る限り読み手には見分けがつかない。結果自身に
`provenance.result_id` と `input_fingerprint` を持たせ、あとから照合できるようにする。

**入力が変わったとき、どの結果が古くなるのか。** 前処理をやり直せば PCA も差次的解析も
無効になる。群の付け替えだけなら PCA は生きている。この依存関係を 1 か所
（`invalidate_results`）に書き、呼び出し側ごとに違う判断をさせない。

指紋は**バイト列**から作る。JSON 経由にすると NaN が null になって 0 と混ざり、
「欠測のまま」と「0 として埋めた」が同じ指紋になる——前処理の違いそのものが
消える。dtype と shape も混ぜる（転置は別の行列で、要素の並びだけでは区別できない）。
"""
from __future__ import annotations

import hashlib
import uuid

import numpy as np

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "PP_FIELDS",
    "STAT_FIELDS",
    "array_fingerprint",
    "assert_current",
    "dataset_fingerprint",
    "group_fingerprint",
    "invalidate_results",
    "is_current",
    "metadata_fingerprints",
    "new_provenance",
    "preprocess_fingerprint",
    "register_result",
]

#: 前処理の結果を左右するサンプルメタデータの列。ここが変われば前処理済み行列が
#: 変わるので、派生結果（PCA・差次的解析）まで無効になる。`group` は入っていない
#: ——群は比較の指定であって前処理の入力ではない。
PP_FIELDS = frozenset({"role", "batch", "injection_order", "qc_pool", "include",
                       "sample_id", "source_file"})

#: 比較・検定の設計を決めるサンプルメタデータの**列**。前処理の入力ではない
#: （群を付け替えただけで前処理をやり直す羽目にしない）が、ここが変われば
#: 検定結果は古い。`biological_sample_id` が入るのは、どの注入を独立した n と
#: 数えるかが検定そのものの前提だから（spec §7）。
STAT_FIELDS = frozenset({"group", "biological_sample_id"})

#: 差次的解析だけを無効化する変更。列（`STAT_FIELDS`）に加えて、列ではない
#: 比較指定そのもの（`comparison`）を含む。
_STAT_INVALIDATING = STAT_FIELDS | {"comparison"}

#: 前処理そのものをやり直す必要がある変更。`standard_assays`（内部標準の
#: 分母に使う標準品注入）・`bindings`（どのfeatureがその標準か）・
#: `matrix_recipe` を含めるのは、これらが変われば解析行列の値そのものが
#: 変わるため——行列を作り直さずに結果だけ残すと、次の解析が新しい設定の
#: 顔をした古い数字で走る。
_PP_INVALIDATING = PP_FIELDS | {"dataset", "detection", "recipe",
                                "standard_assays", "bindings", "matrix_recipe"}


def array_fingerprint(array) -> str:
    """NumPy 配列の指紋。dtype・shape・並びを含むバイト列から作る。

    `array` が None なら「無い」ことの指紋を返す（空配列と区別する）。
    """
    if array is None:
        return canonical_hash({"array": None})
    arr = np.ascontiguousarray(array)
    h = hashlib.sha256()
    h.update(str(arr.dtype).encode("ascii"))
    h.update(b"|")
    h.update(str(arr.shape).encode("ascii"))
    h.update(b"|")
    h.update(arr.tobytes())
    return h.hexdigest()


def dataset_fingerprint(ds) -> str:
    """解析の入力としてのデータセットの指紋。

    含めるのは「入力」だけ——定量行列・特徴 ID・検体の並び・検出マスク・出所。
    派生結果（pp_matrix / last_pca / last_differential）は含めない。含めると
    結果を 1 つ計算するたびに入力が変わったことになり、無限に再計算が要る。
    """
    return canonical_hash({
        "matrix": array_fingerprint(getattr(ds, "feature_matrix", None)),
        "features": list(getattr(ds, "feature_ids", []) or []),
        "samples": list(getattr(ds, "sample_names", []) or []),
        "assays": list(getattr(ds, "sample_assay_ids", []) or []),
        "measure": getattr(ds, "quantification_measure", None),
        "detected_mask": array_fingerprint(getattr(ds, "detected_mask", None)),
        "source_files": dict(getattr(ds, "source_files", {}) or {}),
    })


def metadata_fingerprints(rows: list[dict]) -> dict:
    """サンプルメタデータを**列ごとに**指紋化する。

    列ごとに分けるのは、何が変わったかで無効化の範囲が変わるから。全体を 1 つの
    指紋にすると、群ラベルを直しただけで前処理からやり直す羽目になる。
    行の並び順には依存させない（`sample_id` で整列してから見る）。
    """
    ordered = sorted(rows, key=lambda r: str(r.get("sample_id", "")))
    fields: set[str] = set()
    for row in ordered:
        fields.update(k for k in row if isinstance(k, str))
    return {field: canonical_hash([row.get(field) for row in ordered])
            for field in sorted(fields)}


def group_fingerprint(ds) -> str:
    """比較（差次的解析）の入力になる群割当の指紋。

    前処理の指紋（`preprocess_fingerprint`）には `group` を混ぜない——群を
    付け替えただけで前処理をやり直す羽目になるから。その代わり、群に依存する
    結果（差次的解析）が「今の群割当から出たものか」を言えるように、群だけを
    別の指紋にしてここへ分ける。

    明示メタデータ（実験情報シート）が無いデータセットでは「群という入力が
    無い」ことの指紋を返す。`sample_id` で整列してから見るので行の並び順には
    依存しない（`metadata_fingerprints` と同じ規則）。
    """
    rows = getattr(ds, "sample_metadata_rows", None)
    if not rows:
        return canonical_hash({"group": None, "biological_sample_id": None})
    fingerprints = metadata_fingerprints(rows)
    return canonical_hash({field: fingerprints.get(field)
                           for field in sorted(STAT_FIELDS)})


def preprocess_fingerprint(ds, recipe: dict, metadata_hash: str | None = None) -> str:
    """前処理の入力（データセット・レシピ・関与するメタデータ）の指紋。"""
    return canonical_hash({
        "dataset": dataset_fingerprint(ds),
        "recipe": recipe,
        "metadata": metadata_hash,
    })


def invalidate_results(ds, changed: set[str]) -> None:
    """変更内容に応じて、古くなった状態だけを捨てる。

    捨てる範囲は変更の性質で決まる:

    - 前処理の入力（`PP_FIELDS` ・`dataset` ・`detection` ・`recipe`）が変われば、
      前処理済み行列ごと捨てる。行列を残して結果だけ消すと、次の PCA が古い行列で
      走り、しかもそれが新しい設定の結果として記録される。
    - 群・比較の指定（`STAT_FIELDS`。どの注入を独立した n と数えるかを決める
      `biological_sample_id` を含む）が変われば差次的解析だけ。PCA は群を
      入力にしていない。
    - PCA の設定だけなら PCA だけ。

    「消す」を選ぶのは、古い数字を返し続けるより無いことにするほうが安全だから。
    """
    changed = set(changed)
    if changed & _PP_INVALIDATING:
        ds.pp_matrix = None
        ds.pp_sample_names = []
        ds.pp_feature_names = []
        ds.preprocessing_recipe = {}
        ds.preprocess_id = None
        ds.preprocess_metadata_hash = None
        ds.last_pca = None
        ds.last_differential = None
    elif changed & _STAT_INVALIDATING:
        ds.last_differential = None
    elif "pca_settings" in changed:
        ds.last_pca = None


def new_provenance(ds, *, kind: str, input_fingerprint: str,
                   effective_parameters: dict,
                   parent_ids: list[str] | None = None,
                   warnings: list[str] | None = None,
                   request_revision: int | None = None) -> dict:
    """結果に付ける来歴を作る（result_id はここでだけ発行する）。

    `request_revision` は pipeline 要求の版で、単体ツール経由では None。
    ID を発行するのは新しい計算をしたときだけ——古い結果の result_id を別の計算へ
    付け替えると、来歴が「同じ結果」と嘘をつく。
    """
    from metabolomix.core.version import server_version

    return {
        "result_id": f"res_{uuid.uuid4().hex}",
        "kind": kind,
        "dataset_id": getattr(ds, "dataset_id", None),
        "request_revision": request_revision,
        "input_fingerprint": input_fingerprint,
        "parent_ids": list(parent_ids or []),
        "code_version": server_version(),
        "effective_parameters": effective_parameters,
        "warnings": list(warnings or []),
    }


def register_result(ds, result: dict) -> dict:
    """結果を ID で引けるようにデータセットへ登録する。"""
    ds.results[result["provenance"]["result_id"]] = result
    return result


def is_current(ds, result: dict) -> bool:
    """その結果が、今のデータセットの状態から出たものとして通用するか。

    **来歴を持たない結果は True**（古いと判定しない）。来歴が無いのは「古い」の
    証拠ではなく「調べる材料が無い」ということで、そこを False にすると、来歴を
    付ける前から動いていた経路（ARF の旧結果・手で組んだ結果）の図が一律に
    描けなくなる。来歴があるなら厳密に照合する。

    照合は「その結果が何に依存しているか」で変わる。前処理・データセットの
    一致だけを見ると、**群依存**の結果を取りこぼす: 群だけを訂正した場合、
    前処理は正しく生き残る（`invalidate_results` が `last_differential` だけを
    捨てる）ので、古い差次的結果は親の preprocess_id も dataset_id も
    一致したままになる。`ds.results` に残っているそれを result_id で名指しすれば
    「現在の結果」として通ってしまい、旧群の数字が現在の前処理条件のラベル付きで
    書き出される——`dataset_export_differential` が拒否すると約束している当の
    ケースがこれ。そこで差次的結果は群の指紋（`group_fingerprint`）まで見る。

    群の指紋を持たない差次的結果（この検査より前に作られたもの・手で組んだ
    もの）は、来歴無しと同じく「調べる材料が無い」として通す。
    """
    if not (result or {}).get("provenance"):
        return True
    prov = result["provenance"]
    if prov.get("dataset_id") != getattr(ds, "dataset_id", None):
        return False
    parents = prov.get("parent_ids") or []
    if parents and getattr(ds, "preprocess_id", None) not in parents:
        return False
    if prov.get("kind") == "preprocess":
        return prov.get("result_id") == getattr(ds, "preprocess_id", None)
    if prov.get("kind") == "differential":
        recorded_groups = prov.get("group_fingerprint")
        if recorded_groups is not None and recorded_groups != group_fingerprint(ds):
            return False
    return True


def assert_current(ds, result: dict) -> None:
    """古い結果の利用を止める。理由は機械可読に返す。

    黙って使わせると、図と TSV が別々の前処理から出た数字を並べる（群だけを
    訂正した場合も同じで、旧群の数字が現在の条件の顔をして出て行く）。
    """
    if is_current(ds, result):
        return
    prov = (result or {}).get("provenance") or {}
    raise DomainError(
        "STALE_ANALYSIS_RESULT",
        "この結果は現在のデータセットの状態から出たものではありません"
        "（前処理のやり直し・別データセットの結果）。再計算してください。",
        {"result_id": prov.get("result_id"),
         "result_dataset_id": prov.get("dataset_id"),
         "dataset_id": getattr(ds, "dataset_id", None),
         "result_parent_ids": prov.get("parent_ids"),
         "preprocess_id": getattr(ds, "preprocess_id", None),
         "result_group_fingerprint": prov.get("group_fingerprint"),
         "group_fingerprint": group_fingerprint(ds)})
