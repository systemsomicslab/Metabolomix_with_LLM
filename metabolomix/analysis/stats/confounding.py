"""群とバッチ・注入順の交絡の検査（統計エンジン spec §9）。

単位（units）は属性 `group` `batch` `injection_order`（float | None）を持つ任意の
オブジェクト。生体試料だけを渡す（QC・ブランクは呼び出し側で除く）。
`groups` に入っていない群の単位は数えない。
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats as _sps

from metabolomix.core.atomic_io import DomainError

CRAMERS_V_THRESHOLD = 0.5
INJECTION_ORDER_ALPHA = 0.05


def _cramers_v(table: np.ndarray) -> float | None:
    """群 × バッチの分割表の Cramér の V（Yates 補正なし）。計算できなければ None。"""
    r, c = table.shape
    if r < 2 or c < 2:
        return None
    n = table.sum()
    if n <= 0:
        return None
    try:
        chi2 = float(_sps.chi2_contingency(table, correction=False)[0])
    except ValueError:
        return None
    v = math.sqrt(chi2 / (n * (min(r, c) - 1)))
    return float(v) if math.isfinite(v) else None


def _is_finite_number(x) -> bool:
    return x is not None and isinstance(x, (int, float, np.number)) and math.isfinite(float(x))


def _injection_order(counted: list, present_groups: list[str]) -> tuple[dict, list[dict]]:
    result: dict = {"status": "not_evaluable", "kruskal_p": None, "range_by_group": {}}
    warnings: list[dict] = []
    orders = {g: [float(u.injection_order) for u in counted
                  if u.group == g and _is_finite_number(u.injection_order)]
              for g in present_groups}
    for g, vals in orders.items():
        if vals:
            result["range_by_group"][g] = [min(vals), max(vals)]
    n_with = sum(len(v) for v in orders.values())
    if n_with == 0:
        return result, warnings  # 注入順が無いだけなら警告しない
    reason = None
    if n_with < len(counted):
        reason = "注入順が欠けた生体試料があるため、注入順の偏りは評価できません。"
    elif any(len(v) < 2 for v in orders.values()):
        reason = "注入順が 2 本未満の群があるため、注入順の偏りは評価できません。"
    elif len(present_groups) < 2:
        reason = "群が 1 つだけなので、注入順の偏りは評価できません。"
    elif len({x for v in orders.values() for x in v}) < 2:
        reason = "全試料の注入順が同一で退化しているため、注入順の偏りは評価できません。"
    if reason:
        warnings.append({"code": "injection_order_not_evaluable", "message": reason})
        return result, warnings
    try:
        p = float(_sps.kruskal(*[orders[g] for g in present_groups]).pvalue)
    except ValueError:
        warnings.append({"code": "injection_order_not_evaluable",
                         "message": "注入順の Kruskal–Wallis 検定を計算できませんでした。"})
        return result, warnings
    result["kruskal_p"] = p
    if p < INJECTION_ORDER_ALPHA:
        result["status"] = "biased"
        rng = "、".join(f"{g}: {result['range_by_group'][g][0]:g}–{result['range_by_group'][g][1]:g}"
                        for g in present_groups)
        warnings.append({
            "code": "injection_order_biased",
            "message": f"群ごとの注入順に偏りがあります（Kruskal–Wallis p = {p:.3g}）。"
                       f"群ごとの注入順の範囲: {rng}。ドリフトが群差に混ざる恐れがあります。",
        })
    else:
        result["status"] = "ok"
    return result, warnings


def confounding_report(units, groups: list[str]) -> dict:
    """群 × バッチの表と注入順の偏りを調べる。"""
    gset = set(groups)
    counted = [u for u in units if u.group in gset]
    present_groups = sorted({u.group for u in counted})
    warnings: list[dict] = []

    batches = {u.batch for u in counted}
    missing = any(b is None for b in batches)
    batch_names = sorted(b for b in batches if b is not None)

    raw: dict[str, dict[str, int]] = {g: {} for g in present_groups}
    for u in counted:
        if u.batch is None:
            continue
        raw[u.group][u.batch] = raw[u.group].get(u.batch, 0) + 1
    batch_table = {g: {b: raw[g][b] for b in sorted(raw[g])} for g in present_groups}

    cramers_v: float | None = None
    if not counted or missing:
        status = "not_evaluable"
        warnings.append({"code": "batch_not_evaluable",
                         "message": "バッチが不明な生体試料があるため、群とバッチの交絡は評価できません"
                                    "（交絡なしとは言えません）。"})
    elif len(batch_names) == 1:
        status = "single_batch"
    elif len(present_groups) < 2:
        status = "none"  # 群が 1 つだけなら交絡する相手がいない
    else:
        table = np.array([[batch_table[g].get(b, 0) for b in batch_names] for g in present_groups], float)
        cramers_v = _cramers_v(table)
        groups_per_batch = (table > 0).sum(axis=0)
        batches_per_group = (table > 0).sum(axis=1)
        if (groups_per_batch == 1).all() and (batches_per_group == 1).all():
            status = "complete"
            warnings.append({"code": "batch_group_complete_confounding",
                             "message": "どのバッチも 1 群だけで、どの群も 1 バッチだけです。"
                                        "群差とバッチ差を区別できません（完全交絡）。"})
        elif (groups_per_batch == 1).any() or (cramers_v is not None and cramers_v >= CRAMERS_V_THRESHOLD):
            status = "partial"
            vtxt = f"{cramers_v:.2f}" if cramers_v is not None else "計算不能"
            warnings.append({"code": "batch_group_partial_confounding",
                             "message": f"群とバッチが部分的に交絡しています（Cramér の V = {vtxt}、"
                                        f"1 群しか含まないバッチ: {int((groups_per_batch == 1).sum())} 個）。"
                                        "群差にバッチ差が混ざる恐れがあります。"})
        else:
            status = "none"

    inj, inj_warn = _injection_order(counted, present_groups)
    warnings.extend(inj_warn)
    return {"batch_table": batch_table, "batch_status": status, "cramers_v": cramers_v,
            "injection_order": inj, "warnings": warnings}


def assert_not_confounded(report: dict, *, allow_confounded: bool) -> None:
    """完全交絡なら、明示の許可がない限り STATISTIC_CONFOUNDED で止める。"""
    if report.get("batch_status") == "complete" and not allow_confounded:
        raise DomainError(
            "STATISTIC_CONFOUNDED",
            "群とバッチが完全に交絡しているため、群差を検定できません。"
            "それでも進めるなら allow_confounded を指定してください。",
            {"batch_table": report.get("batch_table")},
        )
