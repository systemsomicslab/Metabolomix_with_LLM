"""結果の同一性・来歴・派生結果の無効化（spec §6）。

数値層の純ロジックで、MCP にも session にも依存しない。ここが答えるのは 2 つの問い:

**この結果はどの入力から出たのか。** 図と TSV が別々の前処理から出た数字を並べても、
どちらも「直近の結果」を名乗る限り読み手には見分けがつかない。結果自身に
`provenance.result_id` と `input_fingerprint` を持たせ、あとから照合できるようにする。

**入力が変わったとき、どの結果が古くなるのか。** 前処理の入力が変われば PCA も差次的解析も
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
    "STATS_KINDS",
    "STAT_FIELDS",
    "array_fingerprint",
    "assert_current",
    "dataset_fingerprint",
    "invalidate_results",
    "is_current",
    "is_stats_current",
    "metadata_fingerprints",
    "new_provenance",
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

#: 前処理の入力の変更（統計エンジンの行列 `stat_matrices` を捨てる）。
_PP_INVALIDATING = PP_FIELDS | {"dataset", "detection", "exclusion"}

#: 統計エンジンの結果の `provenance.kind`（`stats_results` が登録する）。古さは
#: `is_stats_current`（データの指紋と `metadata_revision`）で見る。
STATS_KINDS = ("compare", "stat_pca")


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
    派生結果（stat_matrices / last_pca / last_differential）は含めない。含めると
    結果を 1 つ計算するたびに入力が変わったことになり、無限に再計算が要る。
    """
    mask = getattr(ds, "detected_mask", None)
    return canonical_hash({
        "matrix": array_fingerprint(getattr(ds, "feature_matrix", None)),
        "features": list(getattr(ds, "feature_ids", []) or []),
        "samples": list(getattr(ds, "sample_names", []) or []),
        "assays": list(getattr(ds, "sample_assay_ids", []) or []),
        "measure": getattr(ds, "quantification_measure", None),
        "detected_mask": array_fingerprint(mask),
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


#: `dataset_plot` の kind のうち、PCA の結果・差次的結果から作るもの（`ds.figures` のキー）。
#: `species_pca`（選んだ分子種の PCA）は前処理済み行列を使わないが、除外・前処理の入力の変更で
#: 行列と同じく古くなるので `pca` と同じ組に入れる。
_PCA_FIGURES = ("pca", "loadings", "species_pca")
_DIFFERENTIAL_FIGURES = ("volcano",)


def _drop_figures(ds, kinds: tuple[str, ...]) -> None:
    figures = getattr(ds, "figures", None)
    if figures:
        for kind in kinds:
            figures.pop(kind, None)


def invalidate_results(ds, changed: set[str]) -> None:
    """変更内容に応じて、古くなった状態だけを捨てる。

    捨てる範囲は変更の性質で決まる:

    - 前処理の入力（`PP_FIELDS` ・`dataset` ・`detection`）と除外（`exclusion`。`dataset_exclude`）が
      変われば、統計エンジンの行列（`stat_matrices`）を空にし、直近の PCA・差次的結果とそれらの図を捨てる。
      `stat_matrices` は `stats_results.get_matrix` も古さを見るが、古い行列を抱え続けない。
    - 群・比較の指定（`STAT_FIELDS`。どの注入を独立した n と数えるかを決める
      `biological_sample_id` を含む）が変われば差次的解析だけ。PCA は群を
      入力にしていない。
    - PCA の設定だけなら PCA だけ。

    統計エンジンの結果（`ds.results` の `compare` / `stat_pca`）は消さず、`is_current` が
    古いと答える（`metadata_revision` とデータの指紋）。

    「消す」を選ぶのは、古い数字を返し続けるより無いことにするほうが安全だから。
    """
    changed = set(changed)
    # 結果から作った図（dataset_plot の置き場 ds.figures）も同じ依存で捨てる。残すと save_figure が
    # 古い結果の図を現在のものとして保存する。
    if changed & _PP_INVALIDATING:
        ds.stat_matrices = {}
        ds.last_pca = None
        ds.last_differential = None
        _drop_figures(ds, _PCA_FIGURES + _DIFFERENTIAL_FIGURES)
    elif changed & _STAT_INVALIDATING:
        ds.last_differential = None
        _drop_figures(ds, _DIFFERENTIAL_FIGURES)
    elif "pca_settings" in changed:
        ds.last_pca = None
        _drop_figures(ds, _PCA_FIGURES)


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


def is_stats_current(ds, record: dict, *, dataset_hash: str | None = None) -> bool:
    """統計エンジンの結果（`compare` / `stat_pca`）が、今のデータセットから出たものか。

    登録時のデータの指紋（`provenance.dataset_hash`）と `metadata_revision` が今の ds と同じなら
    True。除外（`dataset_exclude`）と試料の表の変更（群だけの訂正も）は `metadata_revision` を
    進めるので、どちらでも古くなる。結果は `dataset_id` に縛らない（同じデータなら同じ結果）。
    来歴・指紋が無ければ確かめられないので False。`dataset_hash` は呼び出し側が計算済みの
    `dataset_fingerprint(ds)`（任意）。
    """
    prov = (record or {}).get("provenance") or {}
    recorded = prov.get("dataset_hash")
    if recorded is None:
        return False
    if prov.get("metadata_revision") != getattr(ds, "metadata_revision", 0):
        return False
    current = dataset_hash if dataset_hash is not None else dataset_fingerprint(ds)
    return recorded == current


def is_current(ds, result: dict, *, dataset_hash: str | None = None) -> bool:
    """その結果が、今のデータセットの状態から出たものとして通用するか。

    **来歴を持たない結果は True**（古いと判定しない）。来歴が無いのは「古い」の
    証拠ではなく「調べる材料が無い」ということで、そこを False にすると、来歴を
    付ける前から動いていた経路（手で組んだ結果）の図が一律に描けなくなる。
    来歴があるなら厳密に照合する。

    統計エンジンの結果（`provenance.kind` が `STATS_KINDS`）は `is_stats_current` で見る
    （データの指紋と `metadata_revision`）。`dataset_hash` は呼び出し側が計算済みの
    `dataset_fingerprint(ds)`（任意）。指紋は行列全体の hash なので、多数の結果を続けて照合する
    呼び出し側（`feature_search.select_results`）は 1 回だけ計算して渡す。

    それ以外（分子種の PCA `species_pca`）は別のデータセットの結果なら古い。分子種の PCA は
    前処理済み行列を使わず群・除外を引数と試料の対応表から取るので、試料の情報
    （群・役割・除外。`metadata_revision`）が動いたら古い。キュレーションの判断の更新では古くしない。
    """
    if not (result or {}).get("provenance"):
        return True
    prov = result["provenance"]
    if prov.get("kind") in STATS_KINDS:
        return is_stats_current(ds, result, dataset_hash=dataset_hash)
    if prov.get("dataset_id") != getattr(ds, "dataset_id", None):
        return False
    if prov.get("kind") == "species_pca":
        recorded = prov.get("metadata_revision")
        if recorded is not None and recorded != getattr(ds, "metadata_revision", 0):
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
    if prov.get("kind") in STATS_KINDS:
        raise DomainError(
            "STALE_ANALYSIS_RESULT",
            "この結果は現在のデータセットの状態から出たものではありません"
            "（除外・試料の表・データが登録の後に変わった）。再計算してください。",
            {"result_id": prov.get("result_id"), "kind": prov.get("kind"),
             "result_dataset_hash": prov.get("dataset_hash"),
             "dataset_hash": dataset_fingerprint(ds),
             "result_metadata_revision": prov.get("metadata_revision"),
             "metadata_revision": getattr(ds, "metadata_revision", 0)})
    raise DomainError(
        "STALE_ANALYSIS_RESULT",
        "この結果は現在のデータセットの状態から出たものではありません"
        "（試料の情報の変更・別データセットの結果）。再計算してください。",
        {"result_id": prov.get("result_id"), "kind": prov.get("kind"),
         "result_dataset_id": prov.get("dataset_id"),
         "dataset_id": getattr(ds, "dataset_id", None),
         "result_metadata_revision": prov.get("metadata_revision"),
         "metadata_revision": getattr(ds, "metadata_revision", 0)})
