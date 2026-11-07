"""比較のエンジン（統計エンジン spec §5〜§9）。

入口は `compare` の 1 つ。moderated t / F（limma 型）・welch・student・mannwhitney をすべてここから回す。

- 量（§6・§7）: 個体 × 特徴の log2（`PreparedMatrix.log2`）を検定する。欠測は欠測のまま。
  モデルのどれかの cell で有限値の個体が 2 未満の特徴は全対比で `not_tested`（`insufficient_units`）とし、
  当てはめにも事前分布の推定にも入れない。
- 出方（§5）: `detected_units` があれば、単純な対比ごとに検出の例数の Fisher の正確検定。BH は量とは別に、
  対比ごとに Fisher の p を持つ特徴全体で掛ける。どちらかの群で実検出の割合が
  `min_detected_fraction` 未満なら `quantity_judgeable=false`（判定は `inconclusive`）。
- 判定（§7・§15）: ①検定していない・§5 の印 → `inconclusive`、②q ≤ 閾値 → `changed`、
  ③区間が ±幅の内 → `unchanged`、④それ以外 → `inconclusive`。区間が無ければ `unchanged` にしない。
- BH の母集団（§8）と交絡（§9）の表を結果に必ず書く。

MCP・セッションには依存しない。数字は丸めない（非有限は None）。
"""
from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
from scipy import stats

from metabolomix.analysis import spec_statistics
from metabolomix.analysis.stats.confounding import assert_not_confounded, confounding_report
from metabolomix.analysis.stats.design import Design, build_design
from metabolomix.analysis.stats.limma import contrasts_fit, ebayes, lm_fit
from metabolomix.analysis.stats.preprocess import PreparedMatrix
from metabolomix.core.atomic_io import DomainError

__all__ = ["TESTS", "EFFECT_SIZES", "BH_POPULATIONS", "SCHEMA", "CI_LEVEL", "compare"]

TESTS = ("moderated_t", "welch", "student", "mannwhitney")
EFFECT_SIZES = ("log2_mean_difference", "arithmetic_mean_ratio")
BH_POPULATIONS = ("matrix", "annotated", "reported")
SCHEMA = "stats-compare.v1"
CI_LEVEL = 0.95
_MIN_UNITS = 2

_EFFECT_TEXT = {
    ("log2_mean_difference", True): "線形モデルの対比の係数（log2 の上の平均差＝幾何平均比の log2）。"
                                    "区間はモデルの SE と moderated の自由度から。正なら test が高い。",
    ("log2_mean_difference", False): "群の log2 平均の差（test − reference。幾何平均比の log2）。正なら test が高い。",
    ("arithmetic_mean_ratio", True): "log2(mean(2**test) / mean(2**reference))（算術平均比の log2。"
                                     "交互作用では cell の算術平均の log2 に対比の係数を掛けて足す）。区間は出さない。",
}
_EFFECT_TEXT[("arithmetic_mean_ratio", False)] = _EFFECT_TEXT[("arithmetic_mean_ratio", True)]


def _spec_invalid(message: str, **details) -> DomainError:
    return DomainError("STATISTIC_SPECIFICATION_INVALID", message, details)


def _num(x) -> float | None:
    """有限なら float、それ以外（None・NaN・inf）は None。"""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def _model_value(unit, name: str, *, block: bool = False) -> str | None:
    """design._value と同じ読み方で値を読む。無ければ None。"""
    if block and name == "subject_id":
        raw = getattr(unit, "subject_id", None)
    elif name in ("group", "batch"):
        raw = getattr(unit, name, None)
    else:
        raw = (getattr(unit, "attributes", None) or {}).get(name)
    return None if raw is None or str(raw) == "" else str(raw)


def _validate(test, effect_size, bh_population, reported_feature_ids, q_threshold,
              equivalence_fc, min_detected_fraction, factors, comparisons, covariates, block):
    if test not in TESTS:
        raise _spec_invalid(f"test は {TESTS} のどれか", field="test", value=test)
    if effect_size not in EFFECT_SIZES:
        raise _spec_invalid(f"effect_size は {EFFECT_SIZES} のどれか", field="effect_size",
                            value=effect_size)
    if bh_population not in BH_POPULATIONS:
        raise _spec_invalid(f"bh_population は {BH_POPULATIONS} のどれか", field="bh_population",
                            value=bh_population)
    if bh_population == "reported" and reported_feature_ids is None:
        raise _spec_invalid("bh_population=reported には reported_feature_ids が要る",
                            field="reported_feature_ids")
    checks = (("q_threshold", q_threshold, lambda v: 0 < v <= 1),
              ("equivalence_fc", equivalence_fc, lambda v: v > 1),
              ("min_detected_fraction", min_detected_fraction, lambda v: 0 <= v <= 1))
    for name, value, ok in checks:
        if not (isinstance(value, (int, float)) and math.isfinite(value) and ok(value)):
            raise _spec_invalid(f"{name} の値が範囲の外", field=name, value=value)
    if not comparisons:
        raise _spec_invalid("比較（comparisons）が 1 つも無い", field="comparisons")
    if test != "moderated_t" and (len(factors) != 1 or covariates or block is not None):
        raise _spec_invalid(f"{test} は 1 因子・共変量なし・ブロックなしの比較だけ", field="test",
                            value=test)


def _not_tested(reason: str) -> dict:
    return {"status": "not_tested", "reason": reason, "log2fc": None, "ci_low": None,
            "ci_high": None, "ci_unavailable_reason": None, "statistic": None, "df": None,
            "p": None}


def _tested(log2fc, ci_low, ci_high, statistic, df, p, ci_reason=None) -> dict:
    return {"status": "tested", "reason": None, "log2fc": _num(log2fc), "ci_low": _num(ci_low),
            "ci_high": _num(ci_high), "ci_unavailable_reason": ci_reason,
            "statistic": _num(statistic), "df": _num(df), "p": _num(p)}


# --- 量の検定 -------------------------------------------------------------------

def _moderated(log2, design: Design, cids, eligible, want_omnibus):
    """moderated t（対比ごと）と、必要なら moderated F（omnibus）。"""
    g = log2.shape[1]
    out = {cid: [_not_tested("insufficient_units") for _ in range(g)] for cid in cids}
    omnibus = [None] * g if want_omnibus else None
    idx = np.flatnonzero(eligible)
    if idx.size == 0:
        return out, None, omnibus
    fit = lm_fit(log2[:, idx].T, design.matrix, design.columns)
    cmat = np.column_stack([design.contrasts[c] for c in cids])
    try:
        mod = ebayes(contrasts_fit(fit, cmat, list(cids)), ci_level=CI_LEVEL)
    except ValueError:  # どの特徴にも残差の自由度・有限の分散が無い
        for cid in cids:
            for j in idx:
                out[cid][j] = _not_tested("no_residual_degrees_of_freedom")
        return out, None, omnibus
    # 事後分散が 0 か浮動小数の丸めの大きさ（値の大きさ × 1e3·eps 以下の SD）なら退化とみなす。
    # 群の中で値が完全に同じ特徴は、残差が丸めの誤差だけになり t が 1e14 のような値になるため。
    scale = np.nanmax(np.abs(log2[:, idx]), axis=0)
    tol = (1e3 * np.finfo(float).eps * np.maximum(scale, 1.0)) ** 2
    with np.errstate(invalid="ignore"):
        ok = np.isfinite(mod.s2_post) & (mod.s2_post > tol)
    for k, cid in enumerate(cids):
        for pos, j in enumerate(idx):
            if not ok[pos]:
                out[cid][j] = _not_tested("degenerate_variance")
            elif not (math.isfinite(mod.coefficients[pos, k]) and math.isfinite(mod.t[pos, k])):
                out[cid][j] = _not_tested("not_estimable")
            else:
                out[cid][j] = _tested(mod.coefficients[pos, k], mod.ci_low[pos, k],
                                      mod.ci_high[pos, k], mod.t[pos, k], mod.df_total[pos],
                                      mod.p[pos, k])
    prior = {"s2_prior": float(mod.s2_prior), "df_prior": float(mod.df_prior)}
    if want_omnibus and design.omnibus.shape[1] > 0:
        names = [f"{c}-{design.cells[0]}" for c in design.cells[1:]]
        om = ebayes(contrasts_fit(fit, design.omnibus, names), ci_level=CI_LEVEL)
        for pos, j in enumerate(idx):
            if ok[pos] and _num(om.f[pos]) is not None and _num(om.f_p[pos]) is not None:
                omnibus[j] = (float(om.f[pos]), float(om.f_p[pos]))
    return out, prior, omnibus


def _two_group_ci(a, b, kind):
    """log2 の値の平均差・区間・自由度（welch は Welch–Satterthwaite、student はプールした SE）。"""
    d = b.mean() - a.mean()
    na, nb = a.size, b.size
    va, vb = a.var(ddof=1), b.var(ddof=1)
    if kind == "welch":
        sa, sb = va / na, vb / nb
        se = math.sqrt(sa + sb)
        df = (sa + sb) ** 2 / (sa ** 2 / (na - 1) + sb ** 2 / (nb - 1))
    else:
        df = na + nb - 2
        se = math.sqrt(((na - 1) * va + (nb - 1) * vb) / df * (1 / na + 1 / nb))
    half = stats.t.ppf((1 + CI_LEVEL) / 2, df) * se
    return d, d - half, d + half, df


def _two_group_tests(log2, cell_idx, info, cids, eligible, kind):
    g = log2.shape[1]
    out = {}
    for cid in cids:
        ref_rows = cell_idx[info[cid]["reference_cells"][0]]
        test_rows = cell_idx[info[cid]["test_cells"][0]]
        res = []
        for j in range(g):
            if not eligible[j]:
                res.append(_not_tested("insufficient_units"))
                continue
            a, b = log2[ref_rows, j], log2[test_rows, j]
            r = spec_statistics.two_group(a, b, kind=kind, transform="none")
            if r["status"] != "tested":
                res.append(_not_tested(r["reason"]))
                continue
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            if kind == "mannwhitney":
                res.append(_tested(b.mean() - a.mean(), None, None, r["statistic"], None, r["p"],
                                   ci_reason="rank_test"))
            else:
                d, lo, hi, df = _two_group_ci(a, b, kind)
                res.append(_tested(d, lo, hi, r["statistic"], df, r["p"]))
        out[cid] = res
    return out


def _log2_arith_mean(x: np.ndarray) -> float | None:
    """log2(mean(2**x))（有限値だけ。桁あふれしないよう最大値で揃える）。"""
    x = x[np.isfinite(x)]
    if not x.size:
        return None
    m = float(x.max())
    return m + math.log2(float(np.mean(np.exp2(x - m))))


def _apply_arithmetic_ratio(quant, log2, design: Design, cell_idx, cids):
    """検定した行の log2fc を算術平均比の log2 に置き換え、区間を外す。"""
    for cid in cids:
        vec = design.contrasts[cid]
        weights = [(c, vec[k]) for k, c in enumerate(design.cells) if vec[k] != 0]
        for j, row in enumerate(quant[cid]):
            if row["status"] != "tested":
                continue
            parts = [(_log2_arith_mean(log2[cell_idx[c], j]), w) for c, w in weights]
            row["log2fc"] = (None if any(v is None for v, _ in parts)
                             else _num(sum(v * w for v, w in parts)))
            row["ci_low"] = row["ci_high"] = None
            row["ci_unavailable_reason"] = "arithmetic_mean_ratio"


# --- 判定 -----------------------------------------------------------------------

def _verdict(row: dict, q_threshold: float, eq_log2: float) -> tuple[str, str]:
    if row["status"] != "tested":
        return "inconclusive", "not_tested"
    if row["quantity_judgeable"] is False:
        return "inconclusive", "quantity_not_judgeable"
    if row["q"] is not None and row["q"] <= q_threshold:
        return "changed", "q_below_threshold"
    if row["ci_low"] is None or row["ci_high"] is None:
        return "inconclusive", "no_ci"
    if -eq_log2 < row["ci_low"] and row["ci_high"] < eq_log2:
        return "unchanged", "ci_within_equivalence"
    return "inconclusive", "neither"


# --- 入口 -----------------------------------------------------------------------

def compare(prepared: PreparedMatrix, *, factors: list[str], comparisons: list[dict],
            covariates: list[dict] = (), block: str | None = None, test: str = "moderated_t",
            effect_size: str = "log2_mean_difference", bh_population: str = "matrix",
            reported_feature_ids: list[str] | None = None, q_threshold: float = 0.05,
            equivalence_fc: float = 1.2, min_detected_fraction: float = 0.5,
            allow_confounded: bool = False) -> dict:
    """前処理済みの行列で、指定の対比ごとに量と出方を検定し、判定まで付ける（spec §5〜§9）。"""
    factors, comparisons, covariates = list(factors), list(comparisons), [dict(c) for c in covariates]
    _validate(test, effect_size, bh_population, reported_feature_ids, q_threshold, equivalence_fc,
              min_detected_fraction, factors, comparisons, covariates, block)
    warnings = [dict(w) for w in prepared.warnings]

    # モデルの変数（因子・共変量・ブロック）の値が欠けた個体は解析から外す
    names = factors + [c.get("name") for c in covariates]
    keep = [i for i, u in enumerate(prepared.units)
            if all(_model_value(u, n) is not None for n in names)
            and (block is None or _model_value(u, block, block=True) is not None)]
    if len(keep) < len(prepared.units):
        warnings.append({"code": "units_missing_model_values",
                         "message": f"因子・共変量・ブロックの値が欠けた {len(prepared.units) - len(keep)} "
                                    "個体を比較から外しました。"})
    units = [prepared.units[i] for i in keep]
    log2 = np.asarray(prepared.log2, dtype=float)[keep]
    det = None if prepared.detected_units is None else np.asarray(prepared.detected_units, bool)[keep]
    feature_ids = list(prepared.feature_ids)
    n_feat = len(feature_ids)

    # cell と対比（共変量・ブロックなしで先に作り、群の個体数と交絡を当てはめの前に調べる）
    base = build_design(units, factors=factors, comparisons=comparisons)
    cids = [c["id"] for c in comparisons]
    info = base.contrast_info
    unit_cell = [":".join(_model_value(u, f) for f in factors) for u in units]
    cell_idx = {c: [i for i, uc in enumerate(unit_cell) if uc == c] for c in base.cells}
    group_n = {c: len(cell_idx[c]) for c in base.cells}
    used_cells = list(dict.fromkeys(c for cid in cids
                                    for c in info[cid]["reference_cells"] + info[cid]["test_cells"]))
    cell_n = {c: group_n[c] for c in used_cells}
    if any(n < _MIN_UNITS for n in cell_n.values()):
        raise DomainError("STATISTIC_GROUP_TOO_SMALL",
                          f"対比に出る群の個体数が {_MIN_UNITS} 未満（反復注入は個体でまとめた後で数える）",
                          {"cell_n": cell_n})

    proxies = [SimpleNamespace(group=unit_cell[i], batch=u.batch, injection_order=u.injection_order)
               for i, u in enumerate(units)]
    confounding = confounding_report(proxies, used_cells)
    assert_not_confounded(confounding, allow_confounded=allow_confounded)
    warnings.extend(dict(w) for w in confounding["warnings"])

    design = (build_design(units, factors=factors, comparisons=comparisons,
                           covariates=covariates, block=block)
              if covariates or block is not None else base)

    finite = np.isfinite(log2)
    eligible = np.ones(n_feat, dtype=bool)
    for c in base.cells:
        eligible &= finite[cell_idx[c]].sum(axis=0) >= _MIN_UNITS

    # 量
    want_omnibus = test == "moderated_t" and (len(base.cells) >= 3 or len(cids) >= 2)
    if test == "moderated_t":
        quant, prior, omnibus_raw = _moderated(log2, design, cids, eligible, want_omnibus)
    else:
        quant = _two_group_tests(log2, cell_idx, info, cids, eligible, test)
        prior, omnibus_raw = None, None
    if effect_size == "arithmetic_mean_ratio":
        _apply_arithmetic_ratio(quant, log2, design, cell_idx, cids)

    # BH の母集団
    if bh_population == "matrix":
        population = set(feature_ids)
    elif bh_population == "annotated":
        population = set(prepared.annotated_features)
    else:
        population = {str(f) for f in reported_feature_ids}
    eq_log2 = math.log2(equivalence_fc)

    rows: list[dict] = []
    bh_n: dict[str, int] = {}
    presence_n: dict[str, int] = {}
    for cid in cids:
        ref_cells, test_cells = info[cid]["reference_cells"], info[cid]["test_cells"]
        ref_rows = [i for c in ref_cells for i in cell_idx[c]]
        test_rows = [i for c in test_cells for i in cell_idx[c]]
        simple = info[cid]["kind"] == "simple"
        reference = info[cid]["reference"] if simple else " + ".join(ref_cells)
        test_label = info[cid]["test"] if simple else " + ".join(test_cells)

        ps = [quant[cid][j]["p"] if feature_ids[j] in population else None for j in range(n_feat)]
        qs = spec_statistics.bh(ps)
        bh_n[cid] = sum(p is not None for p in ps)

        pres_p: list = [None] * n_feat
        if det is not None and simple:
            for j in range(n_feat):
                pres_p[j] = spec_statistics.detection_fisher_p(
                    int(det[ref_rows, j].sum()), len(ref_rows),
                    int(det[test_rows, j].sum()), len(test_rows))
            presence_n[cid] = sum(p is not None for p in pres_p)
        pres_q = spec_statistics.bh(pres_p)

        for j, fid in enumerate(feature_ids):
            if det is None:
                det_ref = det_test = judgeable = None
            else:
                det_ref, det_test = int(det[ref_rows, j].sum()), int(det[test_rows, j].sum())
                judgeable = all(det[cell_idx[c], j].sum() / len(cell_idx[c]) >= min_detected_fraction
                                for c in ref_cells + test_cells)
            row = {"feature_id": fid, "contrast_id": cid, "reference": reference,
                   "test": test_label, "status": quant[cid][j]["status"],
                   "reason": quant[cid][j]["reason"],
                   "n_reference": int(finite[ref_rows, j].sum()),
                   "n_test": int(finite[test_rows, j].sum()),
                   "detected_reference": det_ref, "detected_test": det_test,
                   "quantity_judgeable": judgeable}
            for key in ("log2fc", "ci_low", "ci_high", "ci_unavailable_reason", "statistic", "df", "p"):
                row[key] = quant[cid][j][key]
            row["q"] = _num(qs[j])
            row["presence_p"], row["presence_q"] = _num(pres_p[j]), _num(pres_q[j])
            row["verdict"], row["verdict_reason"] = _verdict(row, q_threshold, eq_log2)
            rows.append(row)

    omnibus = None
    if want_omnibus:
        f_q = spec_statistics.bh([o[1] if o else None for o in omnibus_raw])
        omnibus = [{"feature_id": fid, "f": o[0] if o else None, "f_p": o[1] if o else None,
                    "f_q": _num(f_q[j])}
                   for j, (fid, o) in enumerate(zip(feature_ids, omnibus_raw))]

    presence = {"detection_available": det is not None, "bh_n_by_contrast": presence_n,
                "min_detected_fraction": min_detected_fraction}
    if det is None:
        presence["reason"] = "detected_mask_unavailable"

    return {
        "schema": SCHEMA,
        "matrix_id": prepared.matrix_id,
        "test": test,
        "effect_size_definition": {"name": effect_size,
                                   "description": _EFFECT_TEXT[(effect_size, test == "moderated_t")]},
        "ci_level": CI_LEVEL,
        "verdict_rule": {"q_threshold": q_threshold, "equivalence_fc": equivalence_fc,
                         "equivalence_log2": eq_log2},
        "bh_population": {"name": bh_population, "n_by_contrast": bh_n,
                          "q_interpretation": ("reproduction_only" if bh_population == "reported"
                                               else "discovery_fdr")},
        "presence": presence,
        "prior": prior,
        "design": {"columns": list(design.columns), "factors": factors, "covariates": covariates,
                   "block": block},
        "group_n": group_n,
        "subject_source": prepared.subject_source,
        "confounding": confounding,
        "warnings": warnings,
        "rows": rows,
        "omnibus": omnibus,
    }
