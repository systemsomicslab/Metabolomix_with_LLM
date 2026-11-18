"""クラスの組成（統計エンジン spec §11）。

個体ごとにクラス内の分子種を足した合計（`2**log2`、欠測は 0）から、
- CLR（クラス間の相対的な取り分: `log2(合計) − 除いたクラスを除く平均`）
- 絶対量（クラスの合計の log2）
の 2 つを作り、それぞれ既存の `compare` で比べる。総量だけが動く変化は絶対量にだけ出て、
取り分の変化は CLR に出る。合計が 0 以下の個体が 1 つでもあるクラスは CLR から除き、
値を補って検定することはしない。

欠測のメンバーは合計に 0 として入る。spec §11 に従い、それを `zero_handling` に
（`member_missing_rule` と、クラス × 群ごとの欠測のメンバー × 個体のセル数 `missing_member_cells`）
書き、欠測のメンバーを含むクラスがあれば警告 `class_total_includes_missing_members` を出す。

脂質アッセイ（`assay_kind="lipid"`）では、脂質名の文法に乗らない注釈の特徴（研究室 MSP の天然物・医薬品の名前など）を
脂質クラスに数えない。乗るかは `lipid_identity.is_lipid_annotation`（クラスが MS-DIAL の脂質クラスか、名前を GOSLIN で
読めるか）。乗らない特徴は `non_lipid`（ラベル `non-lipid`）にまとめて件数と元のクラス（ontology）の内訳を返し、
クラスの合計・CLR・絶対量・欠測の数えには入れない。名前の無い特徴は判定しない（クラスがあればそのクラスに数える）。
ほかの種別では `non_lipid` は None。
"""
from __future__ import annotations

import dataclasses

import numpy as np

from metabolomix.analysis.stats.compare import compare
from metabolomix.analysis.stats.preprocess import PreparedMatrix
from metabolomix.msdial.lipid_identity import has_annotation_name, is_lipid_annotation

SCHEMA = "stats-composition.v1"
ZERO_RULE = "exclude_class_with_nonpositive_total"
MEMBER_MISSING_RULE = "missing_member_counts_as_zero"
_NO_GROUP = "(none)"
#: 脂質アッセイで脂質クラスに数えない特徴の入れ物のラベルと規則。
NON_LIPID_LABEL = "non-lipid"
NON_LIPID_RULE = "lipid_name_unparsable"
_NON_LIPID_EXAMPLES = 10


def _class_prepared(prepared: PreparedMatrix, suffix: str, names: list[str],
                    values: np.ndarray, raw_values: np.ndarray | None = None) -> PreparedMatrix:
    """クラスを特徴に見立てた PreparedMatrix（検出マスク・注釈・QC は持たない）。

    `raw_values` は注入ごとの読んだ値のクラスの合計（片群一定の印が読む）。渡さなければ持たない（compare は
    個体の値で調べる）。注入ごとの線形値（`sample_values`）は特徴の列なので持ち越さない。
    """
    return dataclasses.replace(
        prepared, matrix_id=f"{prepared.matrix_id}:{suffix}", log2=values,
        feature_ids=list(names), detected_units=None, qc_log2=np.empty((0, len(names))),
        annotated_features=frozenset(), feature_classes={}, feature_names={}, sample_values=None,
        sample_detected=None, sample_raw_values=raw_values)


def _raw_class_totals(prepared: PreparedMatrix, members: dict[str, list[int]],
                      names: list[str]) -> np.ndarray | None:
    """注入ごとの読んだ値（`sample_raw_values`）のクラスの合計（欠測・0 以下のメンバーは 0。合計の規則と同じ）。"""
    raw = getattr(prepared, "sample_raw_values", None)
    if raw is None:
        return None
    raw = np.asarray(raw, dtype=float)
    with np.errstate(invalid="ignore"):
        positive = np.where(np.isfinite(raw) & (raw > 0), raw, 0.0)
    return np.column_stack([positive[:, members[c]].sum(axis=1) for c in names]) if names else None


def _non_lipid_block(prepared: PreparedMatrix, ids: list[str]) -> dict:
    """脂質クラスから外した特徴の入れ物: 件数・元のクラス（ontology）ごとの件数・特徴 ID・例（クラスごとに 1 つ）。"""
    by_ontology: dict[str, int] = {}
    first: dict[str, str] = {}
    for fid in ids:
        cls = prepared.feature_classes[fid]
        by_ontology[cls] = by_ontology.get(cls, 0) + 1
        first.setdefault(cls, fid)
    order = sorted(by_ontology, key=lambda c: (-by_ontology[c], c))
    names = getattr(prepared, "feature_names", None) or {}
    return {"label": NON_LIPID_LABEL, "rule": NON_LIPID_RULE, "n_features": len(ids),
            "by_ontology": {c: by_ontology[c] for c in order},
            "feature_ids": list(ids),
            "examples": [{"feature_id": first[c], "name": names.get(first[c]), "ontology": c}
                         for c in order[:_NON_LIPID_EXAMPLES]]}


def composition(prepared: PreparedMatrix, *, assay_kind: str, factors, comparisons,
                **compare_kwargs) -> dict:
    """クラスの合計の CLR と絶対量を、同じ比較条件で並べて返す。

    脂質アッセイでは、脂質名の文法に乗らない注釈の特徴を `non_lipid` に分けてクラスから外す（モジュールの説明）。"""
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
    judge = assay_kind == "lipid"
    feature_names = getattr(prepared, "feature_names", None) or {}
    non_lipid_ids: list[str] = []
    for j, fid in enumerate(prepared.feature_ids):
        cls = prepared.feature_classes.get(fid)
        if cls is None:
            unclassified += 1
            continue
        if judge and has_annotation_name(name := feature_names.get(fid)):
            parses = is_lipid_annotation(name, cls)
            if parses is None:   # pygoslin が無い: 判定できないので外さない（警告を出す）
                judge = False
            elif not parses:
                non_lipid_ids.append(fid)
                continue
        classes.setdefault(cls, []).append(fid)
        members.setdefault(cls, []).append(j)
    non_lipid = None
    if assay_kind == "lipid":
        if judge:
            non_lipid = _non_lipid_block(prepared, non_lipid_ids)
        else:
            warnings.append({"code": "lipid_grammar_unavailable",
                             "message": "脂質名の文法（pygoslin）が使えないため、脂質名の文法に乗らない注釈を"
                                        "脂質クラスから外していません（non_lipid は判定していません）。"})
    names = sorted(classes)
    classes = {c: classes[c] for c in names}

    if not names:
        message = "クラスの付いた特徴が無いため、組成（CLR・絶対量）を計算できません。"
        if non_lipid is not None and non_lipid["n_features"]:
            message = (f"脂質クラスの付いた特徴が無いため、組成（CLR・絶対量）を計算できません"
                       f"（脂質名の文法に乗らない注釈の特徴 {non_lipid['n_features']} 個は non_lipid に外しました）。")
        warnings.append({"code": "no_feature_classes", "message": message})

    n = log2.shape[0]
    unit_groups = [u.group if u.group is not None else _NO_GROUP for u in prepared.units]
    group_levels = sorted(set(unit_groups))
    missing_cells: dict[str, dict[str, int]] = {}
    for c in names:
        miss = (~np.isfinite(log2[:, members[c]])).sum(axis=1)   # 個体ごとの欠測のメンバー数
        missing_cells[c] = {g: int(sum(m for m, ug in zip(miss, unit_groups) if ug == g))
                            for g in group_levels}
    with_missing = [c for c in names if any(missing_cells[c].values())]
    if with_missing:
        warnings.append({"code": "class_total_includes_missing_members",
                         "message": "欠測のメンバーを 0 として合計したクラスがあります"
                                    f"（{', '.join(with_missing)}）。合計は観測されたメンバーだけの和で、"
                                    "欠測の多い群では低めに出ます。",
                         "classes": with_missing})
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
        absolute = compare(_class_prepared(prepared, "class_total", names, log_total,
                                           _raw_class_totals(prepared, members, names)),
                           factors=factors, comparisons=comparisons, **compare_kwargs)

    return {
        "schema": SCHEMA,
        "class_label_kind": "lipid_class" if assay_kind == "lipid" else "msdial_ontology",
        "classes": classes,
        "unclassified_features": unclassified,
        "non_lipid": non_lipid,
        "zero_handling": {"rule": ZERO_RULE, "excluded_classes": excluded,
                          "member_missing_rule": MEMBER_MISSING_RULE,
                          "missing_member_cells": missing_cells},
        "clr": clr,
        "absolute": absolute,
        "unit": prepared.quant.unit if prepared.quant is not None else None,
        "warnings": warnings,
    }
