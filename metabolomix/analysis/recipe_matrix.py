"""行列の recipe v2（spec 2026-10-22 §4 の 4）: 論文の方法どおりの順序で前処理する。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。既存の前処理
（`preprocessing.preprocess` の固定順・`matrix_state` の analysis-matrix.v1）は
変えない——ここは再解析の spec から呼ぶ別の入口で、既定もここだけのもの。

**段の順序は recipe が決める**（`steps`）。論文によって「補完 → 正規化 → ドリフト補正」
だったり「正規化 → 補完」だったりし、順序が違えば数字が違う。挙げた段だけを、挙げた順に
走らせる。挙げない段の設定が recipe にあれば拒否する（書いたのに効かない設定を作らない）。

**標準の定義に合わせる。** KNN は R `impute::impute.knn`（`maxp` = 特徴数、すなわち
2 段の分割をしない）と同じ規則、LOESS は `statsmodels.nonparametric.smoothers_lowess.lowess`。
その場の実装の癖には合わせない。

行列は（試料 × 特徴）。入力の配列は書き換えない。
"""
from __future__ import annotations

import copy
import math

import numpy as np

from metabolomix.core.atomic_io import DomainError

__all__ = ["RECIPE_SCHEMA", "STEPS", "build_matrix", "normalize_recipe"]

RECIPE_SCHEMA = "matrix-recipe.v2"
STEPS = ("detection_filter", "impute", "normalize", "drift_correct", "qc_rsd_filter")
ROLES = ("qc", "sample", "blank")

_INVALID = "MATRIX_RECIPE_INVALID"

#: 段ごとの既定（正規化で埋めて返す）。
DEFAULTS: dict[str, dict] = {
    "detection_filter": {"basis": "detected", "min_rate": {}},
    "impute": {"method": "half_min", "missing": "zero", "k": 10, "scale": "log2",
               "rowmax": 0.5, "colmax": 0.8},
    "normalize": {"method": "pqn", "reference": "qc"},
    "drift_correct": {"method": "loess", "span": 0.75, "iterations": 0, "min_qc": 4},
    "qc_rsd_filter": {"max_rsd": 0.3},
}

_CHOICES = {
    ("detection_filter", "basis"): ("detected", "positive"),
    ("impute", "method"): ("half_min", "knn", "none"),
    ("impute", "missing"): ("zero", "nan"),
    ("impute", "scale"): ("log2", "none"),
    ("normalize", "method"): ("pqn", "tic", "median", "none"),
    ("normalize", "reference"): ("qc", "sample"),
    ("drift_correct", "method"): ("loess",),
}


def _fail(message: str, **details) -> DomainError:
    return DomainError(_INVALID, message, details)


def _number(value, label: str, *, low=None, high=None, integer=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise _fail(f"{label} は数で指定する: {value!r}", field=label)
    if integer and int(value) != value:
        raise _fail(f"{label} は整数で指定する: {value!r}", field=label)
    if (low is not None and value < low) or (high is not None and value > high):
        raise _fail(f"{label} は {low}〜{high} の範囲: {value!r}", field=label)
    return int(value) if integer else float(value)


def _normalize_step(step: str, raw) -> dict:
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise _fail(f"{step} は辞書で指定する", field=step)
    unknown = sorted(set(raw) - set(DEFAULTS[step]))
    if unknown:
        raise _fail(f"{step} の未知のキー: {unknown}", field=step, unknown=unknown)
    out = copy.deepcopy(DEFAULTS[step])
    out.update(copy.deepcopy(raw))
    for (name, key), choices in _CHOICES.items():
        if name == step and out[key] not in choices:
            raise _fail(f"{step}.{key} は {' / '.join(choices)}: {out[key]!r}",
                        field=f"{step}.{key}")
    if step == "detection_filter":
        rates = out["min_rate"]
        if not isinstance(rates, dict) or not set(rates) <= set(ROLES):
            raise _fail(f"detection_filter.min_rate は役割（{'/'.join(ROLES)}）→率の辞書",
                        field="detection_filter.min_rate")
        out["min_rate"] = {role: _number(rate, f"min_rate.{role}", low=0.0, high=1.0)
                           for role, rate in sorted(rates.items())}
    elif step == "impute":
        out["k"] = _number(out["k"], "impute.k", low=1, integer=True)
        out["rowmax"] = _number(out["rowmax"], "impute.rowmax", low=0.0, high=1.0)
        out["colmax"] = _number(out["colmax"], "impute.colmax", low=0.0, high=1.0)
    elif step == "drift_correct":
        out["span"] = _number(out["span"], "drift_correct.span", low=1e-9, high=1.0)
        out["iterations"] = _number(out["iterations"], "drift_correct.iterations", low=0,
                                    integer=True)
        out["min_qc"] = _number(out["min_qc"], "drift_correct.min_qc", low=2, integer=True)
    elif step == "qc_rsd_filter":
        out["max_rsd"] = _number(out["max_rsd"], "qc_rsd_filter.max_rsd", low=0.0)
    return out


def normalize_recipe(recipe) -> dict:
    """recipe を検証して既定を埋める。`steps` に挙げた段の設定だけを残す。"""
    if recipe is None:
        recipe = {}
    if not isinstance(recipe, dict):
        raise _fail("recipe は辞書で指定する")
    unknown = sorted(set(recipe) - {"schema", "steps", *STEPS})
    if unknown:
        raise _fail(f"recipe の未知のキー: {unknown}", unknown=unknown)
    if recipe.get("schema", RECIPE_SCHEMA) != RECIPE_SCHEMA:
        raise _fail(f"schema は {RECIPE_SCHEMA}: {recipe.get('schema')!r}")
    steps = recipe.get("steps", [])
    if not isinstance(steps, list) or not all(isinstance(s, str) for s in steps):
        raise _fail("steps は段の名前のリスト", field="steps")
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise _fail(f"未知の段: {bad}（使えるのは {list(STEPS)}）", field="steps")
    if len(set(steps)) != len(steps):
        raise _fail(f"steps に重複がある: {steps}", field="steps")
    orphan = sorted(set(recipe) & set(STEPS) - set(steps))
    if orphan:
        raise _fail(f"steps に無い段の設定がある（効かない設定は書かない）: {orphan}",
                    field="steps", orphan=orphan)
    out = {"schema": RECIPE_SCHEMA, "steps": list(steps)}
    for step in steps:
        out[step] = _normalize_step(step, recipe.get(step))
    return out


# ---------- 段 ----------

def _missing_mask(values: np.ndarray, definition: str) -> np.ndarray:
    if definition == "zero":
        return ~np.isfinite(values) | (values <= 0)
    return ~np.isfinite(values)


def _role_rows(sample_ids, roles, role) -> list[int]:
    return [i for i, s in enumerate(sample_ids) if roles.get(s, "sample") == role]


def _detection_filter(values, keep, sample_ids, feature_ids, roles, config, detected):
    """役割ごとの検出率の篩。(values, keep, report) を返す（values はそのまま）。"""
    if config["basis"] == "detected":
        if detected is None:
            raise DomainError(
                "DETECTION_UNAVAILABLE",
                "detection_filter.basis=detected には検出状態（gap-fill の区別）が要るが、"
                "このデータセットは取り込めていない。basis=positive（高さ>0）を使うか、"
                "検出状態を取り込めるデータセットを使う。", {})
        present = np.asarray(detected, dtype=bool)
    else:
        present = np.isfinite(values) & (values > 0)
    rates = {}
    passed = np.ones(values.shape[1], dtype=bool)
    for role, minimum in config["min_rate"].items():
        rows = _role_rows(sample_ids, roles, role)
        if not rows:
            raise DomainError("RECIPE_PRECONDITION",
                              f"detection_filter.min_rate に {role} があるが、その役割の試料が無い",
                              {"role": role})
        rate = present[rows].mean(axis=0)
        rates[role] = len(rows)
        passed &= rate >= minimum
    removed = [f for f, k, p in zip(feature_ids, keep, passed) if k and not p]
    return values, keep & passed, {"n_samples_by_role": rates, "removed": removed}


def _half_min(values, missing):
    out = values.copy()
    positive = np.where(np.isfinite(values) & (values > 0), values, np.nan)
    with np.errstate(all="ignore"):
        fill = np.nanmin(np.where(np.isnan(positive), np.inf, positive), axis=0) / 2.0
    fill = np.where(np.isfinite(fill), fill, np.nan)
    rows, cols = np.nonzero(missing)
    out[rows, cols] = fill[cols]
    return out


def knn_impute(data: np.ndarray, *, k: int, rowmax: float, colmax: float) -> np.ndarray:
    """R `impute::impute.knn`（maxp = 行数）と同じ規則の補完。`data` は（特徴 × 試料）、欠測は NaN。

    R の `knnimp` に倣い、行を 2 つに分ける:

    - 欠測の割合が `rowmax` を超える行（重い行）は、**重い行どうしの**列の平均で埋め
      （`meanimp`）、KNN の近傍にしない。重い行の中にその列の観測が無ければ、軽い行の列の平均。
    - それ以外の行（軽い行）は軽い行だけで KNN: 欠測のある行ごとに、その行が観測された列のうち
      相手も観測された列での平方差の平均を距離にし（共通の列が無い相手は除く）、近い `k` 個を
      選ぶ（同距離は行の順）。欠測の列を、近傍のうちその列を観測しているものの平均で埋める。
      全員が欠測なら**軽い行の**列の平均。
    - 欠測の割合が `colmax` を超える列（試料）があれば RECIPE_PRECONDITION。
    """
    data = np.asarray(data, dtype=float)
    observed = np.isfinite(data)
    n_rows, n_cols = data.shape
    col_missing = 1.0 - observed.mean(axis=0) if n_rows else np.zeros(n_cols)
    over = [int(j) for j in np.nonzero(col_missing > colmax)[0]]
    if over:
        raise DomainError("RECIPE_PRECONDITION",
                          f"欠測の割合が colmax={colmax} を超える試料がある（列 {over}）。"
                          "impute.knn と同じく補完しない。", {"columns": over})
    out = data.copy()
    row_missing = 1.0 - observed.mean(axis=1) if n_cols else np.zeros(n_rows)
    heavy_rows = np.nonzero(row_missing > rowmax)[0]
    light_rows = np.nonzero(row_missing <= rowmax)[0]

    def _column_means(rows):
        if not len(rows):
            return np.full(n_cols, np.nan)
        # 観測の無い列は NaN（numpy の「空の平均」の警告は出さない）
        values = np.where(observed[rows], data[rows], 0.0)
        counts = observed[rows].sum(axis=0)
        with np.errstate(all="ignore"):
            return np.where(counts > 0, values.sum(axis=0) / counts, np.nan)

    light_mean = _column_means(light_rows)
    if len(heavy_rows):
        heavy_mean = _column_means(heavy_rows)
        heavy_mean = np.where(np.isfinite(heavy_mean), heavy_mean, light_mean)
        for i in heavy_rows:
            out[i, ~observed[i]] = heavy_mean[~observed[i]]

    light = data[light_rows]
    light_observed = observed[light_rows]
    zeroed = np.where(light_observed, light, 0.0)
    mask = light_observed.astype(float)
    sq = zeroed ** 2
    targets = np.nonzero(~light_observed.all(axis=1))[0]
    block = 256
    for start in range(0, len(targets), block):
        rows = targets[start:start + block]
        o = mask[rows]
        shared = o @ mask.T
        with np.errstate(all="ignore"):
            dist = (sq[rows] @ mask.T - 2.0 * (zeroed[rows] @ zeroed.T) + o @ sq.T) / shared
        dist = np.where(shared > 0, np.maximum(dist, 0.0), np.inf)
        dist[np.arange(len(rows)), rows] = np.inf
        for r, i in enumerate(rows):
            neighbours = _nearest(dist[r], k)
            missing_cols = np.nonzero(~light_observed[i])[0]
            target_row = light_rows[i]
            if neighbours.size:
                block_values = np.where(light_observed[np.ix_(neighbours, missing_cols)],
                                        light[np.ix_(neighbours, missing_cols)], np.nan)
                counts = np.isfinite(block_values).sum(axis=0)
                with np.errstate(all="ignore"):
                    means = np.nansum(block_values, axis=0) / counts
                out[target_row, missing_cols] = np.where(counts > 0, means,
                                                         light_mean[missing_cols])
            else:
                out[target_row, missing_cols] = light_mean[missing_cols]
    return out


def _nearest(distances: np.ndarray, k: int) -> np.ndarray:
    """距離の小さい順に k 個（同距離は添字の小さい順）。有限の距離だけ。

    全体を並べ替えずに k 番目の値で切る。境界の同距離は添字の順に取るので、
    安定ソートの先頭 k 個と同じ集合になる（決定的）。
    """
    finite = np.isfinite(distances)
    k = min(k, int(finite.sum()))
    if k == 0:
        return np.array([], dtype=int)
    masked = np.where(finite, distances, np.inf)
    threshold = np.partition(masked, k - 1)[k - 1]
    below = np.nonzero(masked < threshold)[0]
    equal = np.nonzero(masked == threshold)[0][: k - below.size]
    return np.concatenate([below, equal])


def _impute(values, config):
    missing = _missing_mask(values, config["missing"])
    report = {"method": config["method"], "missing_cells": int(missing.sum())}
    if config["method"] == "none" or not missing.any():
        return values.copy(), report
    if config["method"] == "half_min":
        out = _half_min(values, missing)
    else:
        work = np.where(missing, np.nan, values)
        if config["scale"] == "log2":
            with np.errstate(all="ignore"):
                work = np.where(work > 0, np.log2(work), np.nan)
        filled = knn_impute(work.T, k=config["k"], rowmax=config["rowmax"],
                            colmax=config["colmax"]).T
        out = np.power(2.0, filled) if config["scale"] == "log2" else filled
        out = np.where(missing, out, values)
    report["filled_cells"] = int((missing & np.isfinite(out)).sum())
    return out, report


def _normalize(values, sample_ids, roles, config):
    method = config["method"]
    if method == "none":
        return values.copy(), {"method": "none"}
    positive = np.where(np.isfinite(values) & (values > 0), values, np.nan)
    if method == "pqn":
        rows = _role_rows(sample_ids, roles, config["reference"])
        if not rows:
            raise DomainError("RECIPE_PRECONDITION",
                              f"normalize.reference={config['reference']} の試料が無い",
                              {"reference": config["reference"]})
        with np.errstate(all="ignore"):
            reference = np.nanmedian(positive[rows], axis=0)
            factors = np.nanmedian(positive / reference, axis=1)
    else:
        with np.errstate(all="ignore"):
            raw = (np.nansum(np.where(np.isfinite(values), values, 0.0), axis=1) if method == "tic"
                   else np.nanmedian(positive, axis=1))
            factors = raw / np.nanmedian(raw[np.isfinite(raw) & (raw > 0)])
    bad = ~np.isfinite(factors) | (factors <= 0)
    factors = np.where(bad, 1.0, factors)
    report = {"method": method,
              "factors": {s: float(f) for s, f in zip(sample_ids, factors)},
              "unscaled_samples": [s for s, b in zip(sample_ids, bad) if b]}
    if method == "pqn":
        report["reference"] = config["reference"]
    return values / factors[:, None], report


def _drift_correct(values, sample_ids, roles, config, injection_order):
    from statsmodels.nonparametric.smoothers_lowess import lowess

    order = injection_order or {}
    lacking = [s for s in sample_ids if order.get(s) is None]
    if lacking:
        raise DomainError("INJECTION_ORDER_MISSING",
                          f"ドリフト補正には全試料の注入順が要る（無い試料: {lacking}）",
                          {"samples": lacking})
    x_all = np.array([float(order[s]) for s in sample_ids])
    qc = np.array(_role_rows(sample_ids, roles, "qc"), dtype=int)
    out = values.copy()
    uncorrected = 0
    for j in range(values.shape[1]):
        y = values[qc, j] if len(qc) else np.array([])
        ok = np.isfinite(y) & (y > 0)
        if ok.sum() < config["min_qc"]:
            uncorrected += 1
            continue
        fit = lowess(y[ok], x_all[qc][ok], frac=config["span"], it=config["iterations"],
                     delta=0.0, xvals=x_all)
        if not np.all(np.isfinite(fit)) or np.any(fit <= 0):
            uncorrected += 1
            continue
        out[:, j] = values[:, j] / fit * np.median(y[ok])
    return out, {"method": "loess", "n_qc": int(len(qc)), "uncorrected_features": uncorrected}


def _qc_rsd_filter(values, keep, sample_ids, feature_ids, roles, config):
    rows = _role_rows(sample_ids, roles, "qc")
    if len(rows) < 2:
        raise DomainError("RECIPE_PRECONDITION", "qc_rsd_filter には QC が 2 本以上要る",
                          {"n_qc": len(rows)})
    high, none = [], []
    passed = np.ones(values.shape[1], dtype=bool)
    for j, fid in enumerate(feature_ids):
        if not keep[j]:
            continue
        y = values[rows, j]
        y = y[np.isfinite(y)]
        rsd = (np.std(y, ddof=1) / np.mean(y)) if len(y) >= 2 and np.mean(y) > 0 else math.nan
        if not math.isfinite(rsd):
            none.append(fid)
            passed[j] = False
        elif rsd > config["max_rsd"]:
            high.append(fid)
            passed[j] = False
    return keep & passed, {"max_rsd": config["max_rsd"], "removed_high_rsd": high,
                           "removed_no_rsd": none}


def build_matrix(values, sample_ids: list[str], feature_ids: list[str], roles: dict[str, str],
                 recipe: dict, *, detected=None, injection_order: dict | None = None) -> dict:
    """正規化済みの recipe を（試料 × 特徴）の行列に当てる。

    篩の段（detection_filter・qc_rsd_filter）で落とした特徴は、**その段の直後に列から外す**。
    後の段（KNN の近傍・PQN の参照・TIC の和）に落とした特徴を混ぜると、論文の方法と
    違う数字になる。
    """
    current = np.array(values, dtype=float, copy=True)
    if current.shape != (len(sample_ids), len(feature_ids)):
        raise DomainError("RECIPE_PRECONDITION", "行列の形が試料・特徴の数と合わない",
                          {"shape": list(current.shape)})
    features = list(feature_ids)
    present = None if detected is None else np.asarray(detected, dtype=bool)
    reports = []
    for step in recipe["steps"]:
        config = recipe[step]
        keep = None
        if step == "detection_filter":
            _, keep, report = _detection_filter(current, np.ones(len(features), dtype=bool),
                                                sample_ids, features, roles, config, present)
        elif step == "impute":
            current, report = _impute(current, config)
        elif step == "normalize":
            current, report = _normalize(current, sample_ids, roles, config)
        elif step == "drift_correct":
            current, report = _drift_correct(current, sample_ids, roles, config,
                                             injection_order)
        else:
            keep, report = _qc_rsd_filter(current, np.ones(len(features), dtype=bool),
                                          sample_ids, features, roles, config)
        if keep is not None:
            columns = np.nonzero(keep)[0]
            current = current[:, columns]
            features = [features[j] for j in columns]
            if present is not None:
                present = present[:, columns]
        reports.append({"step": step, **report, "n_features_after": len(features)})
    return {
        "values": current,
        "sample_ids": list(sample_ids),
        "feature_ids": features,
        "steps": reports,
        "n_features_in": len(feature_ids),
        "n_features_out": len(features),
    }
