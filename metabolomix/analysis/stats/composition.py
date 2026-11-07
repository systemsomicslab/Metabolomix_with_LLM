"""クラスの組成（統計エンジン spec §11）。

個体ごとにクラス内の分子種を足した合計（`2**log2`、欠測は 0）から、
- CLR（クラス間の相対的な取り分: `log2(合計) − 除いたクラスを除く平均`）
- 絶対量（クラスの合計の log2）
の 2 つを作り、それぞれ既存の `compare` で比べる。総量だけが動く変化は絶対量にだけ出て、
取り分の変化は CLR に出る。合計が 0 以下の個体が 1 つでもあるクラスは CLR から除き、
値を補って検定することはしない。
"""
from __future__ import annotations

import dataclasses

import numpy as np

from metabolomix.analysis.stats.compare import compare
from metabolomix.analysis.stats.preprocess import PreparedMatrix

SCHEMA = "stats-composition.v1"
ZERO_RULE = "exclude_class_with_nonpositive_total"


def _class_prepared(prepared: PreparedMatrix, suffix: str, names: list[str],
                    values: np.ndarray) -> PreparedMatrix:
    """クラスを特徴に見立てた PreparedMatrix（検出マスク・注釈・QC は持たない）。"""
    return dataclasses.replace(
        prepared, matrix_id=f"{prepared.matrix_id}:{suffix}", log2=values,
        feature_ids=list(names), detected_units=None, qc_log2=np.empty((0, len(names))),
        annotated_features=frozenset(), feature_classes={})


def composition(prepared: PreparedMatrix, *, assay_kind: str, factors, comparisons,
                **compare_kwargs) -> dict:
    """クラスの合計の CLR と絶対量を、同じ比較条件で並べて返す。"""
    warnings: list[dict] = []
    if compare_kwargs.get("bh_population", "matrix") != "matrix":
        warnings.append({"code": "bh_population_overridden",
                         "message": f"組成の比較は特徴＝クラスなので bh_population を matrix に固定しました"
                                    f"（指定: {compare_kwargs['bh_population']}）。"})
    compare_kwargs["bh_population"] = "matrix"

    log2 = np.asarray(prepared.log2, dtype=float)
    classes: dict[str, list[str]] = {}
    members: dict[str, list[int]] = {}
    unclassified = 0
    for j, fid in enumerate(prepared.feature_ids):
        cls = prepared.feature_classes.get(fid)
        if cls is None:
            unclassified += 1
            continue
        classes.setdefault(cls, []).append(fid)
        members.setdefault(cls, []).append(j)
    names = sorted(classes)
    classes = {c: classes[c] for c in names}

    n = log2.shape[0]
    totals = np.zeros((n, len(names)))
    for k, c in enumerate(names):
        totals[:, k] = np.nansum(np.exp2(log2[:, members[c]]), axis=1)   # 欠測は 0
    positive = totals > 0
    with np.errstate(divide="ignore"):
        log_total = np.where(positive, np.log2(np.where(positive, totals, 1.0)), np.nan)

    excluded = {c: int((~positive[:, k]).sum()) for k, c in enumerate(names) if not positive[:, k].all()}
    keep = [k for k, c in enumerate(names) if c not in excluded]

    clr = None
    if len(keep) >= 2:
        sub = log_total[:, keep]
        clr_values = sub - sub.mean(axis=1, keepdims=True)
        clr = compare(_class_prepared(prepared, "clr", [names[k] for k in keep], clr_values),
                      factors=factors, comparisons=comparisons, **compare_kwargs)
    else:
        warnings.append({"code": "clr_not_computable",
                         "message": "CLR には、合計が全個体で正のクラスが 2 つ以上要ります"
                                    f"（残ったクラス: {len(keep)}）。",
                         "classes_remaining": [names[k] for k in keep]})

    absolute = None
    if names:
        absolute = compare(_class_prepared(prepared, "class_total", names, log_total),
                           factors=factors, comparisons=comparisons, **compare_kwargs)

    return {
        "schema": SCHEMA,
        "class_label_kind": "lipid_class" if assay_kind == "lipid" else "msdial_ontology",
        "classes": classes,
        "unclassified_features": unclassified,
        "zero_handling": {"rule": ZERO_RULE, "excluded_classes": excluded},
        "clr": clr,
        "absolute": absolute,
        "unit": prepared.quant.unit if prepared.quant is not None else None,
        "warnings": warnings,
    }
