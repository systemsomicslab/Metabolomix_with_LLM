"""PLS-DA と VIP（統計エンジン。analysis-spec v2 の解析の種類 `plsda`）。

前処理済みの個体 × 特徴の log2（`PreparedMatrix.log2`）で、`design.factors` の水準の組（cell）を
クラスにした PLS 判別分析を当てはめ、次を**必ず**返す（過学習を数字で見えるようにするため）:

- 成分数（`n_components`。個体数と特徴数から決まる上限で頭打ちにしたら警告）
- 交差検証の Q2（層化グループ k 分割。`q2_cum` は成分 1〜A のそれぞれ。`Q2 = 1 − PRESS / TSS`、TSS は全個体の
  Y の平均の周りの平方和）。R2X・R2Y（当てはめた全個体で）も並べる
- 並べ替え検定の p（クラスの札を `n_permutations` 回入れ替え、同じ分割の交差検証で Q2 を出し直した
  帰無分布で `(1 + #{Q2_perm ≥ Q2}) / (1 + n_permutations)`）
- VIP（全個体の当てはめの変数重要度。`VIP_j = √(p·Σ_a SSY_a·(w_ja/‖w_a‖)² / Σ_a SSY_a)`、SSY_a は成分 a が
  説明する Y の平方和）

**独立な単位（グループ）。** 同じ個体の複数の単位（個体を属性で分けた単位・`collapse_replicates=False` の注入）は
互いに独立でない。単位ごとに分割・並べ替えると、同じ個体の単位が学習と予測の両方に入り、並べ替えは個体の共通の
札を壊すので、群差の無いデータでも Q2 と並べ替えの p が個体の違いで良く見える。分割と並べ替えは独立な単位で行う:
`block`（design.block）があればその値、無ければ `subject_id`、それも無い単位はそれぞれが 1 つ（`cv.group_by`。
`unit` は各単位が 1 グループ）。分割は同じグループの単位を同じ分割に入れ（グループのクラスの組で層化して配る）、
分割数の上限と成分数の上限はグループの数で決まる。並べ替えは、1 つのクラスだけのグループはグループの札ごと
入れ替え、複数のクラスにまたがるグループ（対応のある計画・交差した計画）はグループの中で入れ替える
（`permutation.scheme`: `between_groups` / `within_groups` / `mixed`）。どのクラスも 2 つ以上のグループに要る
（1 つなら学習と予測を分けられない。`STATISTIC_GROUP_TOO_SMALL`）。グループが少なく、入れ替えで作れる札の並びが
少ないと並べ替えの p は小さくなれない（`permutation.min_attainable_p`。0.05 を超えたら警告
`permutation_few_labelings`）。各単位が 1 グループなら、分割・並べ替えの乱数の使い方も結果も単位ごとの
層化 k 分割と同じ。

Y はクラスの指示変数（2 クラスなら後のクラスを 1 とした 1 列、3 クラス以上は 1 クラス 1 列）。X は
`scaling`（`autoscale`＝平均 0・SD 1 / `pareto`＝SD の平方根で割る）で尺度をそろえ、交差検証では
**学習の個体だけ**で中心と尺度を決める。欠測は PCA と同じ規則（`impute_for_pca`）で、検定前に 1 回だけ補う
（交差検証の外。`imputation` に件数）。全個体で値が変わらない特徴は外す（`n_features_dropped_constant`）。
QC は入れない（`PreparedMatrix.log2` は生体の個体だけ）。乱数（分割と並べ替え）は `seed` で決まる。

MCP・セッションには依存しない。数字は丸めない（非有限は None）。
"""
from __future__ import annotations

import math
import warnings as _warnings

import numpy as np

from metabolomix.analysis.stats.compare import _model_value as _value
from metabolomix.analysis.stats.preprocess import PreparedMatrix, impute_for_pca
from metabolomix.core.atomic_io import DomainError

__all__ = ["SCHEMA", "SCALINGS", "IMPUTES", "plsda"]

#: 並べ替えの p がおよそ下回れない値（`min_attainable_p`）がこれを超えたら警告する（独立な単位が少ない）。
_FEW_LABELINGS_P = 0.05
#: 帰無の Q2 を観測の Q2 以上と数える幅（同じ Q2 になる並び＝2 クラスの名前の入れ替えを丸めで取りこぼさない）。
_Q2_TIE_RTOL = 1e-12

SCHEMA = "stats-plsda.v1"
SCALINGS = ("autoscale", "pareto")
IMPUTES = ("half_min", "knn", "none")
_MIN_UNITS = 2
_SD_FLOOR = 1e-12


def _invalid(message: str, **details) -> DomainError:
    return DomainError("PLSDA_SPECIFICATION_INVALID", message, details)


def _num(x) -> float | None:
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def _label(unit, factors) -> str | None:
    values = [_value(unit, f) for f in factors]
    return None if any(v is None for v in values) else ":".join(values)


def _group_key(unit, block: str | None) -> str:
    """分割・並べ替えの独立な単位の鍵: block の値、無ければ subject_id、それも無ければ単位そのもの。"""
    if block is not None:
        return "b:" + str(_value(unit, block, block=True))
    subject = getattr(unit, "subject_id", None)
    if subject is not None and str(subject) != "":
        return "s:" + str(subject)
    return "u:" + str(unit.unit_id)


def _groups_of(groups: np.ndarray) -> dict[str, np.ndarray]:
    """グループの鍵 → その単位の位置（グループは初出の順）。"""
    return {g: np.flatnonzero(groups == g) for g in dict.fromkeys(groups.tolist())}


def _y(labels: np.ndarray, order: list[str]) -> np.ndarray:
    """クラスの札から Y（2 クラスなら後のクラスを 1 とした 1 列、3 クラス以上はクラスごとの指示変数）。"""
    if len(order) == 2:
        return np.column_stack([(labels == order[1]).astype(float)])
    return np.column_stack([(labels == c).astype(float) for c in order])


def _scale_params(x: np.ndarray, scaling: str) -> tuple[np.ndarray, np.ndarray]:
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=1) if x.shape[0] > 1 else np.ones(x.shape[1])
    sd = np.where(sd > _SD_FLOOR, sd, 1.0)
    return mean, (sd if scaling == "autoscale" else np.sqrt(sd))


def _fit(x: np.ndarray, y: np.ndarray, n_components: int):
    from sklearn.cross_decomposition import PLSRegression

    model = PLSRegression(n_components=n_components, scale=False)
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        model.fit(x, y)
    return model


def _predict_by_components(model, x_centered: np.ndarray, y_mean: np.ndarray, n_components: int):
    """成分 1〜A それぞれの予測（A 成分の当てはめの先頭 a 成分＝a 成分の当てはめ。NIPALS の P'W は上三角）。"""
    rot, q = model.x_rotations_, model.y_loadings_
    scores = x_centered @ rot
    return [y_mean + scores[:, :a] @ q[:, :a].T for a in range(1, n_components + 1)]


def _folds(labels: np.ndarray, k: int, rng: np.random.Generator,
           groups: np.ndarray | None = None) -> np.ndarray:
    """層化グループ k 分割の割り当て。同じグループの単位は同じ分割に入る。

    グループをそのクラスの組（1 クラスだけのグループはそのクラス）で層に分け、層ごと（組の順）に
    グループを並べ替えて分割へ順に配る（前の層の続きの番号から）。`groups` が無ければ単位ごと
    （各単位が 1 グループ。クラスごとに並べ替えて配る層化 k 分割と同じ）。
    """
    if groups is None:
        groups = np.arange(labels.size)
    members = _groups_of(np.asarray(groups))
    strata: dict[tuple, list] = {}
    for g, idx in members.items():
        strata.setdefault(tuple(sorted(set(labels[idx].tolist()))), []).append(g)
    fold = np.empty(labels.size, dtype=int)
    start = 0
    for key in sorted(strata):
        chosen = strata[key]
        for pos, pick in enumerate(rng.permutation(len(chosen))):
            fold[members[chosen[pick]]] = (pos + start) % k
        start += len(chosen)
    return fold


def _permutation_plan(labels: np.ndarray, groups: np.ndarray, order: list[str]) -> dict:
    """並べ替えの組み立て: 1 クラスだけのグループ（札ごと入れ替える）と、複数のクラスにまたがるグループ
    （中で入れ替える）。

    `n_labelings` は入れ替えで作れる区別できる札の並びの数、`min_attainable_p` は並べ替えの p がおよそ
    下回れない値（`観測と同じ Q2 になる並びの数 / n_labelings`）。2 クラスの Y は 1 列なので、2 つの
    クラスの名前を入れ替えた並びは同じ Q2 になる（その並びを入れ替えで作れるなら 2 つ）。3 クラス以上は
    名前の入れ替えでも Q2 が丸めの外で変わるので 1 つと数える（下限）。
    """
    members = _groups_of(groups)
    between = [g for g, idx in members.items() if len(set(labels[idx].tolist())) == 1]
    within = [g for g, idx in members.items() if len(set(labels[idx].tolist())) > 1]
    count = math.factorial(len(between))
    between_n = {c: sum(1 for g in between if labels[members[g][0]] == c) for c in order}
    for n_c in between_n.values():
        count //= math.factorial(n_c)
    swap = len(order) == 2 and between_n[order[0]] == between_n[order[1]]
    for g in within:
        idx = members[g]
        part = math.factorial(idx.size)
        inside = {c: int((labels[idx] == c).sum()) for c in order}
        for n_c in inside.values():
            part //= math.factorial(n_c)
        count *= part
        swap = swap and inside[order[0]] == inside[order[1]]
    scheme = ("within_groups" if not between else "between_groups" if not within else "mixed")
    return {"members": members, "between": between, "within": within, "scheme": scheme,
            "n_labelings": count, "min_attainable_p": (2 if swap else 1) / count}


def _permuted(labels: np.ndarray, plan: dict, rng: np.random.Generator) -> np.ndarray:
    """札を 1 回入れ替える。1 クラスのグループはグループの札を入れ替え、またがるグループは中で入れ替える
    （各単位が 1 グループなら、単位の札の並べ替え `labels[rng.permutation(n)]` と同じ）。"""
    members, between = plan["members"], plan["between"]
    out = labels.copy()
    if between:
        source = [labels[members[g][0]] for g in between]
        for g, pick in zip(between, rng.permutation(len(between))):
            out[members[g]] = source[pick]
    for g in plan["within"]:
        idx = members[g]
        out[idx] = labels[idx[rng.permutation(idx.size)]]
    return out


def _cv_q2(x: np.ndarray, y: np.ndarray, fold: np.ndarray, k: int, n_components: int,
           scaling: str) -> list[float]:
    """交差検証の Q2（成分 1〜A）。分割ごとに学習の個体だけで中心と尺度を決めて当てはめる。"""
    press = np.zeros(n_components)
    for f in range(k):
        test = fold == f
        if not test.any():
            continue
        train = ~test
        mean, scale = _scale_params(x[train], scaling)
        xs_train = (x[train] - mean) / scale
        xs_test = (x[test] - mean) / scale
        model = _fit(xs_train, y[train], n_components)
        x_center = xs_train.mean(axis=0)
        preds = _predict_by_components(model, xs_test - x_center, y[train].mean(axis=0), n_components)
        for a, pred in enumerate(preds):
            press[a] += float(((y[test] - pred) ** 2).sum())
    tss = float(((y - y.mean(axis=0)) ** 2).sum())
    return [1.0 - p / tss if tss > 0 else math.nan for p in press]


def _vip(model) -> np.ndarray:
    w = model.x_weights_                       # p × A（列は単位長）
    t = model.x_scores_                        # n × A
    q = model.y_loadings_                      # m × A
    ssy = (t ** 2).sum(axis=0) * (q ** 2).sum(axis=0)
    total = ssy.sum()
    if not total > 0:
        return np.full(w.shape[0], np.nan)
    wn = w / np.linalg.norm(w, axis=0, keepdims=True)
    return np.sqrt(w.shape[0] * ((wn ** 2) @ ssy) / total)


def _replicate_warning(prepared, group_by: str) -> dict | None:
    """単位に個体の鍵が無く、名前が反復注入らしい試料（`possible_technical_replicates`）があるときの注意。"""
    if group_by != "unit":
        return None
    stems = [w for w in (getattr(prepared, "warnings", None) or [])
             if w.get("code") == "possible_technical_replicates"]
    if not stems:
        return None
    return {"code": "plsda_units_may_not_be_independent",
            "message": "単位に subject_id が無く、名前が反復注入らしい試料があります"
                       f"（{len(stems)} 組。前処理の possible_technical_replicates）。交差検証と並べ替えは各単位を"
                       "独立として扱うので、同じ個体の反復なら Q2 と並べ替えの p は甘く出ます。試料の表の subject_id で"
                       "個体を指せば、反復は 1 個体にまとまり、分割と並べ替えも個体で行います。",
            "n_stems": len(stems)}


def plsda(prepared: PreparedMatrix, *, factors: list[str], classes: list[str] | None = None,
          n_components: int = 2, scaling: str = "autoscale", impute: str = "half_min",
          cv_folds: int = 7, n_permutations: int = 200, seed: int = 0,
          block: str | None = None) -> dict:
    """PLS-DA を当てはめ、成分数・交差検証の Q2・並べ替え検定の p・VIP を返す。

    `block`（design.block）があれば、その値が同じ単位を 1 つの独立な単位（グループ）として分割・並べ替えをする。
    無ければ subject_id で、それも無い単位はそれぞれが 1 グループ。
    """
    if scaling not in SCALINGS:
        raise _invalid(f"scaling は {SCALINGS} のどれか", field="scaling", value=repr(scaling))
    if impute not in IMPUTES:
        raise _invalid(f"impute は {IMPUTES} のどれか", field="impute", value=repr(impute))
    for name, value, low in (("n_components", n_components, 1), ("cv_folds", cv_folds, 2),
                             ("n_permutations", n_permutations, 1), ("seed", seed, 0)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < low:
            raise _invalid(f"{name} は {low} 以上の整数", field=name, value=repr(value))
    factors = list(factors)
    warnings: list[dict] = []

    labels_all = [_label(u, factors) for u in prepared.units]
    lacking, variables = [], []
    for u, lab in zip(prepared.units, labels_all):
        missing = [f for f in factors if _value(u, f) is None]
        if block is not None and _value(u, block, block=True) is None:
            missing.append(block)
        if missing:
            lacking.append(u.unit_id)
            variables.extend(v for v in missing if v not in variables)
    if lacking:
        warnings.append({"code": "units_missing_model_values",
                         "message": f"{', '.join(variables)} の値が欠けた {len(lacking)} 個体を外しました。",
                         "unit_ids": lacking, "variables": variables})
    lacking_set = set(lacking)
    present = sorted({lab for u, lab in zip(prepared.units, labels_all)
                      if lab is not None and u.unit_id not in lacking_set})
    if classes is not None:
        unknown = [c for c in classes if c not in present]
        if unknown:
            raise _invalid(f"classes に無い cell がある: {unknown}", field="classes",
                           unknown=unknown, available=present)
        order = list(dict.fromkeys(classes))
    else:
        order = present
    keep = [i for i, (u, lab) in enumerate(zip(prepared.units, labels_all))
            if lab is not None and lab in order and u.unit_id not in lacking_set]
    labels = np.array([labels_all[i] for i in keep])
    counts = {c: int((labels == c).sum()) for c in order}
    if len(order) < 2 or any(n < _MIN_UNITS for n in counts.values()):
        raise DomainError("STATISTIC_GROUP_TOO_SMALL",
                          f"PLS-DA には {_MIN_UNITS} 個体以上のクラスが 2 つ以上要る",
                          {"class_n": counts})
    groups = np.array([_group_key(prepared.units[i], block) for i in keep])
    with_subject = any(getattr(prepared.units[i], "subject_id", None) not in (None, "") for i in keep)
    group_by = block if block is not None else ("subject_id" if with_subject else "unit")
    members = _groups_of(groups)
    class_groups = {c: sum(1 for idx in members.values() if (labels[idx] == c).any()) for c in order}
    if any(n < _MIN_UNITS for n in class_groups.values()):
        raise DomainError("STATISTIC_GROUP_TOO_SMALL",
                          f"PLS-DA の交差検証には、どのクラスも {_MIN_UNITS} つ以上の独立な単位（{group_by}）に要る"
                          "（1 つだけのクラスは学習と予測を分けられない）",
                          {"class_n": counts, "class_groups": class_groups, "group_by": group_by})

    x = np.asarray(prepared.log2, dtype=float)[keep]
    n_missing = int((~np.isfinite(x)).sum())
    if n_missing:
        if impute == "none":
            raise DomainError("PLSDA_MISSING_VALUES", "欠測があるため PLS-DA を計算できません（補完なし）",
                              {"missing_cells": n_missing})
        x, rep = impute_for_pca(x, method=impute)
        imputation = {"method": rep["method"], "imputed_cells": rep["imputed_cells"], "used": True}
    else:
        imputation = {"method": impute, "imputed_cells": 0, "used": False}
    finite_cols = np.isfinite(x).all(axis=0)
    with np.errstate(invalid="ignore"):
        sd = np.where(finite_cols, np.nanstd(np.where(np.isfinite(x), x, 0.0), axis=0), 0.0)
        mean = np.where(finite_cols, np.nanmean(np.where(np.isfinite(x), x, 0.0), axis=0), 0.0)
    kept = finite_cols & (sd > _SD_FLOOR * np.maximum(1.0, np.abs(mean)))
    if int(kept.sum()) < 2:
        raise DomainError("PLSDA_TOO_FEW_FEATURES", "値が変わる特徴が 2 未満のため PLS-DA を計算できません",
                          {"kept": int(kept.sum())})
    x = x[:, kept]
    kept_ids = [fid for fid, k in zip(prepared.feature_ids, kept) if k]
    n, p = x.shape
    n_groups = len(members)

    y = _y(labels, order)
    k = min(int(cv_folds), n_groups)
    if k < cv_folds:
        unit_word = "個体数" if group_by == "unit" else f"独立な単位（{group_by}）の数"
        warnings.append({"code": "cv_folds_capped",
                         "message": f"分割数を{unit_word} {n_groups} に合わせて {k} にしました（1 つ抜き）。"})
    rng = np.random.default_rng(seed)
    fold = _folds(labels, k, rng, groups)
    smallest_train = n - int(np.bincount(fold, minlength=k).max())
    cap = max(1, min(smallest_train - 1, p))
    a = min(int(n_components), cap)
    if a < n_components:
        warnings.append({"code": "n_components_capped",
                         "message": f"成分数を {n_components} から {a} にしました（交差検証の学習の個体数と"
                                    "特徴数から決まる上限）。", "requested": int(n_components),
                         "used": a})
    q2_cum = _cv_q2(x, y, fold, k, a, scaling)

    mean_all, scale_all = _scale_params(x, scaling)
    xs = (x - mean_all) / scale_all
    model = _fit(xs, y, a)
    preds = _predict_by_components(model, xs - xs.mean(axis=0), y.mean(axis=0), a)
    yss = float(((y - y.mean(axis=0)) ** 2).sum())
    r2y = [1.0 - float(((y - pr) ** 2).sum()) / yss for pr in preds]
    xss = float((xs ** 2).sum())
    t, pl = model.x_scores_, model.x_loadings_
    r2x = [1.0 - float(((xs - xs.mean(axis=0) - t[:, :c] @ pl[:, :c].T) ** 2).sum()) / xss
           for c in range(1, a + 1)]
    vip = _vip(model)

    q2 = q2_cum[-1]
    plan = _permutation_plan(labels, groups, order)
    null = np.empty(int(n_permutations))
    for i in range(int(n_permutations)):
        null[i] = _cv_q2(x, _y(_permuted(labels, plan, rng), order), fold, k, a, scaling)[-1]
    finite_null = null[np.isfinite(null)]
    exceed = (int((finite_null >= q2 - _Q2_TIE_RTOL * max(1.0, abs(q2))).sum()) if math.isfinite(q2)
              else len(finite_null))
    perm_p = (1 + exceed) / (1 + len(finite_null)) if len(finite_null) else None
    if not (math.isfinite(q2) and q2 > 0):
        warnings.append({"code": "plsda_q2_not_positive",
                         "message": f"交差検証の Q2 が {q2:.3g} で正でありません（クラスを予測できていない）。"
                                    "VIP は当てはめの重要度で、予測できない模型の VIP は判別の根拠になりません。"})
    if plan["min_attainable_p"] > _FEW_LABELINGS_P:
        warnings.append({"code": "permutation_few_labelings",
                         "message": f"入れ替えで作れるクラスの札の並びは {plan['n_labelings']} 通りしか無く、並べ替えの p は"
                                    f"およそ {plan['min_attainable_p']:.3g} より小さくなれません（独立な単位 {group_by} が"
                                    f"少ない。{n_groups} グループ）。p が大きくても判別が無いとは言えません。",
                         "n_labelings": int(plan["n_labelings"]),
                         "min_attainable_p": float(plan["min_attainable_p"])})
    note = _replicate_warning(prepared, group_by)
    if note is not None:
        warnings.append(note)

    order_vip = np.argsort(-np.nan_to_num(vip, nan=-1.0), kind="stable")
    return {
        "schema": SCHEMA,
        "matrix_id": prepared.matrix_id,
        "factors": factors,
        "classes": counts,
        "y_coding": ("one_column_second_class" if len(order) == 2 else "indicator_per_class"),
        "n_units": int(n),
        "unit_ids": [prepared.units[i].unit_id for i in keep],
        "scaling": scaling,
        "imputation": imputation,
        "n_features": int(p),
        "n_features_dropped_constant": int(len(prepared.feature_ids) - p),
        "n_components_requested": int(n_components),
        "n_components": int(a),
        "r2x_cum": [_num(v) for v in r2x],
        "r2y_cum": [_num(v) for v in r2y],
        "q2_cum": [_num(v) for v in q2_cum],
        "q2": _num(q2),
        "cv": {"method": "stratified_group_kfold", "group_by": group_by, "n_groups": int(n_groups),
               "class_groups": class_groups, "folds": int(k), "seed": int(seed),
               "scaling_within_fold": True, "q2_definition": "1 - PRESS/TSS"},
        "permutation": {"n": int(n_permutations), "n_finite": int(len(finite_null)),
                        "p": _num(perm_p), "statistic": "q2", "scheme": plan["scheme"],
                        "min_attainable_p": float(plan["min_attainable_p"]),
                        "q2_null_mean": _num(finite_null.mean()) if len(finite_null) else None,
                        "q2_null_q95": (_num(np.quantile(finite_null, 0.95))
                                        if len(finite_null) else None)},
        "vip": [{"feature_id": kept_ids[j], "vip": _num(vip[j])} for j in order_vip],
        "scores": {"unit_ids": [prepared.units[i].unit_id for i in keep],
                   "values": [[_num(v) for v in row] for row in t]},
        "warnings": warnings,
    }
