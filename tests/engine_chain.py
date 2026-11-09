"""統計エンジンの公開の入口（`dataset_preprocess` / `dataset_compare` / `dataset_pca`）で状態を作る足場。

書き出し・根拠バンドル・図・特徴検索など、差次的結果や PCA を読む利用者のテストが使う。
小さな fixture でも行が落ちないよう、既定では足切り・ドリフト補正・正規化を切った前処理
（`FLAT`）を当てる。前処理の設定を変えたいテストは `preprocess` に上書きを渡す。
"""
from __future__ import annotations

import json

from metabolomix.tools.stats_tools import dataset_compare, dataset_pca, dataset_preprocess

#: 足切り・ドリフト補正・正規化を切った前処理（小さな fixture の行を落とさない）。
FLAT = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
            max_qc_rsd=None, max_d_ratio=None, normalize="none")


def engine_preprocess(*, dataset_handle: str | None = None, **overrides) -> dict:
    """`dataset_preprocess` を `FLAT` に上書きを重ねて呼び、成功した戻り値（dict）を返す。"""
    out = json.loads(dataset_preprocess(**{**FLAT, **overrides}, dataset_handle=dataset_handle))
    assert out.get("status") == "success", out
    return out


def engine_compare(reference_group: str | None = None, test_group: str | None = None, *,
                   dataset_handle: str | None = None, preprocess: dict | None = None,
                   **compare_kw) -> dict:
    """前処理（`engine_preprocess`）から `dataset_compare` までを回し、比較の戻り値（dict）を返す。"""
    matrix = engine_preprocess(dataset_handle=dataset_handle, **(preprocess or {}))
    out = json.loads(dataset_compare(reference_group=reference_group, test_group=test_group,
                                     matrix_id=matrix["matrix_id"], dataset_handle=dataset_handle,
                                     **compare_kw))
    assert out.get("status") == "success", out
    return out


def engine_pca(*, dataset_handle: str | None = None, preprocess: dict | None = None,
               **pca_kw) -> dict:
    """前処理（`engine_preprocess`）から `dataset_pca` までを回し、PCA の戻り値（dict）を返す。"""
    matrix = engine_preprocess(dataset_handle=dataset_handle, **(preprocess or {}))
    out = json.loads(dataset_pca(matrix_id=matrix["matrix_id"], dataset_handle=dataset_handle,
                                 **pca_kw))
    assert out.get("status") == "success", out
    return out


def declare_groups(ds, group_of):
    """群と役割を MS-DIAL の試料の宣言（`ds.msdial_samples`）として与える。

    `group_of(name)` は試料名 → 群（None なら群なし）。名前が `QC` で始まる試料は役割 qc。
    試料名（`HC1` など）がトークンに割れず `group_levels` で群を解けない fixture に使う。
    """
    ds.msdial_samples = [{"name": name, "group": group_of(name),
                          "role": "qc" if name.startswith("QC") else "sample"}
                         for name in ds.sample_names]
    return ds


def hc_ra(name: str):
    """`tests.reanalysis_fixtures` の試料名（HC1… / RA1… / QC1…）の群。"""
    return name[:2] if name[:2] in ("HC", "RA") else None


#: spec v2 の前処理で `FLAT` と同じもの（ツールを通さずに登録するときに使う）。
FLAT_SPEC = {"drift": {"mode": "off"},
             "filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None},
             "normalize": {"method": "none"}}


def register_two_group(ds, reference: str, test: str, *, log2fc_threshold: float = 1.0) -> dict:
    """ツール（とセッション）を通さずに、統計エンジンの行列と 2 群比較の記録を `ds` に登録して記録を返す。

    `dataset_preprocess` → `dataset_compare` と同じ登録（`ds.stat_matrices` と `ds.results`）。群は
    `ds` の試料の表（試料の表の行・`msdial_samples`）から取る。前処理は `FLAT_SPEC`。
    `dataset_exclude` の除外（`ds.excluded_samples` / `excluded_features`）はツールと同じく反映する。
    """
    from metabolomix.analysis import stats_results
    from metabolomix.analysis.stats.spec_v2 import normalize_spec_v2, prepare_spec_v2, run_spec_v2

    spec = normalize_spec_v2({
        "schema": "analysis-spec.v2", "preprocess": FLAT_SPEC, "design": {"factors": ["group"]},
        "input": {"exclude_samples": sorted(getattr(ds, "excluded_samples", None) or {}),
                  "exclude_features": sorted(getattr(ds, "excluded_features", None) or {})},
        "comparisons": [{"id": f"{test}_vs_{reference}", "factor": "group",
                         "reference": reference, "test": test}],
        "analyses": [{"analysis_id": "compare", "kind": "compare"}]})
    prepared = prepare_spec_v2(spec, ds)
    stats_results.register_matrix(ds, prepared, spec["preprocess"], spec["input"]["group_overrides"])
    return stats_results.register_compare(ds, run_spec_v2(spec, ds, prepared=prepared),
                                          log2fc_threshold=log2fc_threshold, units=prepared.units)


def compare_halves(ds, reference: str = "A", test: str = "B", **kw) -> dict:
    """試料の前半を `reference`・後半を `test` の群に宣言し（`declare_groups`）、2 群比較を登録して記録を返す。"""
    half = len(ds.sample_names) // 2
    position = {name: i for i, name in enumerate(ds.sample_names)}
    declare_groups(ds, lambda name: reference if position[name] < half else test)
    return register_two_group(ds, reference, test, **kw)
