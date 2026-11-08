"""前処理（spec §4）: 選択 → ドリフト補正 → 品質の足切り → 正規化 → log2 → 反復のまとめ。

各段は `steps` に順に記録する（`{"step", "status": "applied"|"skipped", "reason"?, ...}`）。
0 以下・非有限は欠測で、pseudocount も clip も使わない。補完は `impute_for_pca` だけが行う。
MCP・セッションには依存しない。
"""
from __future__ import annotations

import copy
import hashlib
import math
import warnings
from dataclasses import dataclass, field

import numpy as np

from metabolomix.analysis.stats.data import (DETECTED, QuantState, StudyData, cell_states,
                                             replicate_name_warnings)
from metabolomix.core.atomic_io import DomainError, canonical_hash

PREPROCESS_DEFAULTS = {
    "drift": {"mode": "auto", "span": 0.75, "min_qc_per_batch": 4},
    "filters": {"blank_ratio": {"min": 3.0}, "detection": {"by": "group", "min_rate": 0.5},
                "qc_rsd": {"max": 0.30}, "d_ratio": {"max": 0.5}},
    "normalize": {"method": "auto"},
    "collapse_replicates": True,
    # 個体のまとめの鍵から外す試料の表の列（取得日・注入番号など注入ごとに違う属性）。
    "injection_columns": [],
}

_ENUMS = {("drift", "mode"): ("auto", "on", "off"),
          ("normalize", "method"): ("auto", "none", "tic", "median", "pqn"),
          ("filters", "detection", "by"): ("group", "all")}
_NULLABLE = {("filters", name) for name in PREPROCESS_DEFAULTS["filters"]}
_ROLES = ("sample", "qc", "blank")
_NO_BATCH = "(none)"


def _invalid(message: str, path: tuple, value=None) -> DomainError:
    return DomainError("PREPROCESS_CONFIG_INVALID", message,
                       {"field": ".".join(path), "value": repr(value)})


def _number(path: tuple, value, *, lo: float, hi: float = math.inf, lo_open: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise _invalid(f"{'.'.join(path)} は数値で指定する", path, value)
    if value < lo or value > hi or (lo_open and value == lo):
        raise _invalid(f"{'.'.join(path)} が範囲外", path, value)
    return float(value)


def _merge(default, raw, path: tuple):
    if raw is None and path in _NULLABLE:
        return None
    if isinstance(default, dict):
        if not isinstance(raw, dict):
            raise _invalid(f"{'.'.join(path) or '設定'} は object で指定する", path, raw)
        unknown = sorted(set(raw) - set(default))
        if unknown:
            raise _invalid(f"未知のキー: {unknown}", path + (unknown[0],), raw[unknown[0]])
        return {k: _merge(v, raw.get(k, copy.deepcopy(v)), path + (k,)) for k, v in default.items()}
    if path in _ENUMS:
        if raw not in _ENUMS[path]:
            raise _invalid(f"{'.'.join(path)} は {list(_ENUMS[path])} のどれか", path, raw)
        return raw
    if path == ("injection_columns",):
        if (not isinstance(raw, list)
                or any(not isinstance(c, str) or not c.strip() for c in raw)):
            raise _invalid("injection_columns は空でない文字列のリストで指定する", path, raw)
        return list(raw)
    if path == ("collapse_replicates",):
        if not isinstance(raw, bool):
            raise _invalid("collapse_replicates は真偽値で指定する", path, raw)
        return raw
    if path == ("drift", "min_qc_per_batch"):
        n = _number(path, raw, lo=2)
        if n != int(n):
            raise _invalid("min_qc_per_batch は整数で指定する", path, raw)
        return int(n)
    if path == ("drift", "span"):
        return _number(path, raw, lo=0.0, hi=1.0, lo_open=True)
    if path[-1] == "min_rate":
        return _number(path, raw, lo=0.0, hi=1.0)
    return _number(path, raw, lo=0.0)


def normalize_preprocess_config(raw: dict | None) -> dict:
    """既定を埋めた前処理設定を返す。未知のキー・不正な値は PREPROCESS_CONFIG_INVALID。"""
    return _merge(PREPROCESS_DEFAULTS, {} if raw is None else raw, ())


@dataclass
class UnitInfo:
    """検定の単位（個体）1 つ。design.py・confounding.py はこの属性名で読む。"""
    unit_id: str
    group: str | None
    batch: str | None
    subject_id: str | None
    injection_order: float | None
    n_injections: int
    sample_ids: list[str]
    attributes: dict[str, str] = field(default_factory=dict)


@dataclass
class PreparedMatrix:
    """前処理の済んだ個体 × 特徴の log2 行列（NaN = 欠測）と、その記録。"""
    matrix_id: str
    log2: np.ndarray
    units: list[UnitInfo]
    feature_ids: list[str]
    detected_units: np.ndarray | None
    qc_log2: np.ndarray
    qc_ids: list[str]
    steps: list[dict]
    config: dict
    quant: QuantState | None
    dataset_hash: str | None
    subject_source: str
    warnings: list[dict]
    feature_classes: dict
    annotated_features: frozenset
    # 生体（role=sample）の注入 × 残した特徴（feature_ids と同じ順）。正規化の後・log2 の前の線形値
    # （0 以下はそのまま、欠測は NaN）。注入ごとの値を見る配線層（箱ひげ・外れ値の根拠）が読む。
    sample_values: np.ndarray | None = None
    sample_ids: list[str] = field(default_factory=list)
    sample_detected: np.ndarray | None = None  # 検出の別が無い研究は None


@dataclass
class NeedsInput:
    """既定を決めずに利用者へ尋ねる封筒（spec §4 のドリフト補正の既定）。"""
    field: str
    reason: str
    options: list[dict]

    def to_envelope(self) -> dict:
        return {"status": "needs_input", "field": self.field, "reason": self.reason,
                "options": self.options}


def _positive(x: np.ndarray) -> np.ndarray:
    """有限の正値だけを残し、それ以外を NaN にした写し。"""
    x = np.asarray(x, dtype=float)
    return np.where(np.isfinite(x) & (x > 0), x, np.nan)


def _nanmedian(block: np.ndarray) -> np.ndarray:
    """列ごとの nanmedian。全欠測の列は NaN（numpy の警告は出さない）。"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(block, axis=0)


def _nanmin(block: np.ndarray) -> np.ndarray:
    """列ごとの有限値の最小。観測が無い列は NaN。"""
    finite = np.isfinite(block)
    low = np.where(finite, block, np.inf).min(axis=0) if block.shape[0] else np.full(block.shape[1], np.inf)
    return np.where(np.isfinite(low), low, np.nan)


def _matrix_id(study: StudyData, config: dict, normalize_method: str) -> str:
    table = [[s.sample_id, s.role, s.group, s.batch, s.injection_order, s.subject_id,
              bool(s.include), dict(sorted(s.attributes.items()))] for s in study.samples]
    source = study.dataset_hash
    # 特徴の一覧（除外・割れのまとめで変わる）は dataset_hash が同じでも ID に効かせる
    feature_hash = canonical_hash(list(study.feature_ids))
    if source is None:
        # 指紋が無いときは値そのものの bytes から（同じ入力なら同じ ID）
        values = np.ascontiguousarray(study.values, dtype=np.float64)
        det = b"" if study.detected is None else np.ascontiguousarray(study.detected, bool).tobytes()
        source = {"values_sha256": hashlib.sha256(values.tobytes()).hexdigest(),
                  "shape": list(values.shape), "detected_sha256": hashlib.sha256(det).hexdigest(),
                  "feature_ids": list(study.feature_ids)}
    # normalize=auto は quant_state で方法が決まるので、解決後の方法も指紋に入れる
    return "mx_" + canonical_hash({"dataset_hash": source, "sample_table_hash": canonical_hash(table),
                                   "config": config, "normalize_method": normalize_method,
                                   "feature_ids_hash": feature_hash})[:16]


# --- 1. 選択 -----------------------------------------------------------------

def _select(study: StudyData) -> tuple[list[int], dict]:
    keep, by_include, by_role = [], [], {}
    for i, s in enumerate(study.samples):
        if not s.include:
            by_include.append(s.sample_id)
        elif s.role not in _ROLES:
            by_role.setdefault(str(s.role), []).append(s.sample_id)
        else:
            keep.append(i)
    counts = {r: sum(study.samples[i].role == r for i in keep) for r in _ROLES}
    if counts["sample"] == 0:
        raise DomainError("PREPROCESS_PRECONDITION", "統計に使える生体試料（role=sample）が無い",
                          {"excluded_by_include": by_include, "excluded_by_role": by_role})
    step = {"step": "select", "status": "applied", "n_sample": counts["sample"],
            "n_qc": counts["qc"], "n_blank": counts["blank"],
            "excluded_by_include": by_include, "excluded_by_role": by_role}
    return keep, step


# --- 2. ドリフト補正 ---------------------------------------------------------

def _batch_key(s) -> str:
    return s.batch if s.batch is not None else _NO_BATCH


#: LOESS の近傍の最小数。statsmodels は近傍数を int(frac·n) で決め、最も遠い近傍の重みは 0。
#: 近傍が 3 だと等間隔の QC で等距離の 2 点がともに重み 0 になり、xvals 経路の当てはめが NaN になる。
_MIN_LOESS_NEIGHBOURS = 4


def _neighbours(span: float, n: int) -> int:
    return int(span * n + 1e-10)    # statsmodels の lowess と同じ数え方


def _drift_batches(samples) -> list[str]:
    return sorted({_batch_key(s) for s in samples if s.role in ("sample", "qc")})


def _drift_eligibility(samples, cfg: dict) -> list[str]:
    """補正できない理由（空なら全バッチで補正できる）。plan §15: 1 つでも欠ければ全体を飛ばす。"""
    if not any(s.role == "qc" for s in samples):
        return ["no_qc"]
    if any(s.injection_order is None for s in samples if s.role in ("sample", "qc")):
        return ["no_injection_order"]
    reasons = []
    for b in _drift_batches(samples):
        qc = [s.injection_order for s in samples if s.role == "qc" and _batch_key(s) == b]
        bio = [s.injection_order for s in samples if s.role == "sample" and _batch_key(s) == b]
        if len(qc) < cfg["min_qc_per_batch"]:
            reasons.append(f"too_few_qc:{b}")
        elif _neighbours(cfg["span"], len(qc)) < _MIN_LOESS_NEIGHBOURS:
            reasons.append(f"too_few_qc_for_span:{b}")
        elif bio and not (min(qc) <= min(bio) and max(qc) >= max(bio)):
            reasons.append(f"qc_not_interspersed:{b}")
    return reasons


def _drift(values: np.ndarray, samples, cfg: dict) -> dict:
    """バッチごとの LOESS（生の尺度）→ 特徴ごとのバッチ水準合わせ。

    全バッチで補正できた特徴だけを補正し（一部のバッチだけの補正はしない）、水準もその特徴だけ
    合わせる。補正できた特徴が 1 つも無いバッチがあれば全体を飛ばす。成功したときだけ values を書き換える。
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    span, min_qc = cfg["span"], cfg["min_qc_per_batch"]
    batches = _drift_batches(samples)
    n_feat = values.shape[1]
    work = values.copy()
    corrected = np.ones(n_feat, dtype=bool)
    too_few = {b: 0 for b in batches}
    fit_invalid = {b: 0 for b in batches}
    rows_of, qc_of, empty = {}, {}, []
    for b in batches:
        rows = np.array([i for i, s in enumerate(samples)
                         if _batch_key(s) == b and s.injection_order is not None], dtype=int)
        qc = np.array([i for i in rows if samples[i].role == "qc"], dtype=int)
        rows_of[b], qc_of[b] = rows, qc
        x_all = np.array([samples[i].injection_order for i in rows], dtype=float)
        x_qc = np.array([samples[i].injection_order for i in qc], dtype=float)
        in_batch = np.zeros(n_feat, dtype=bool)
        for j in range(n_feat):
            y = _positive(values[qc, j])
            ok = np.isfinite(y)
            if ok.sum() < min_qc or _neighbours(span, int(ok.sum())) < _MIN_LOESS_NEIGHBOURS:
                too_few[b] += 1
                continue
            fit = lowess(y[ok], x_qc[ok], frac=span, it=0, delta=0.0, xvals=x_all)
            if not np.all(np.isfinite(fit)) or np.any(fit <= 0):
                fit_invalid[b] += 1
                continue
            work[rows, j] = values[rows, j] / fit * float(np.median(y[ok]))
            in_batch[j] = True
        if not in_batch.any():
            empty.append(f"fit_invalid:{b}" if fit_invalid[b] else f"no_correctable_feature:{b}")
        corrected &= in_batch
    record = {"step": "drift", "method": "loess", "span": span, "batches": batches,
              "uncorrected_too_few_qc": too_few, "uncorrected_fit_invalid": fit_invalid}
    if not empty and not corrected.any():
        # バッチごとには補正できても、全バッチで補正できる特徴が 1 つも無い
        empty = ["no_feature_correctable_in_all_batches"]
    if empty:
        return {**record, "status": "skipped", "reason": empty[0], "reasons": empty}
    # 一部のバッチでしか補正できなかった特徴は、どのバッチでも補正しない
    cols = np.nonzero(corrected)[0]
    out = values.copy()
    out[:, cols] = work[:, cols]
    # 水準合わせ: 補正後の全 QC の中央値 / そのバッチの補正後の QC 中央値
    all_qc = np.concatenate([qc_of[b] for b in batches])
    global_med = _nanmedian(_positive(out[np.ix_(all_qc, cols)]))
    for b in batches:
        with np.errstate(all="ignore"):
            factor = global_med / _nanmedian(_positive(out[np.ix_(qc_of[b], cols)]))
        out[np.ix_(rows_of[b], cols)] *= factor
    values[:] = out
    done = {i for b in batches for i in rows_of[b]}
    not_corrected = [{"sample_id": s.sample_id,
                      "reason": "no_injection_order" if s.injection_order is None
                      else "batch_without_sample_or_qc"}
                     for i, s in enumerate(samples) if i not in done]
    return {**record, "status": "applied", "corrected_features": int(len(cols)),
            "uncorrected_features": int(n_feat - len(cols)), "rows_not_corrected": not_corrected}


_NEEDS_DRIFT_REASON = (
    "このデータは内部標準などで換算済み（quant_state.normalized）で、QC もある。"
    "内部標準の比はドリフトの一部をすでに打ち消しているため、QC の LOESS でドリフト補正を"
    "重ねると換算済みの値の尺度を動かしうる。補正するか（on）しないか（off）を選んでほしい。")


def _drift_step(values, samples, cfg, quant) -> dict | NeedsInput:
    mode = cfg["mode"]
    if mode == "off":
        return {"step": "drift", "status": "skipped", "reason": "mode_off"}
    has_qc = any(s.role == "qc" for s in samples)
    if mode == "auto" and has_qc and quant is not None and quant.normalized:
        return NeedsInput(field="preprocess.drift.mode", reason=_NEEDS_DRIFT_REASON, options=[
            {"value": "on", "label": "QC の LOESS でドリフト補正する（換算済みの値に重ねる）"},
            {"value": "off", "label": "ドリフト補正をしない（内部標準の換算だけに任せる）"}])
    reasons = _drift_eligibility(samples, cfg)
    step = ({"step": "drift", "status": "skipped", "reason": reasons[0], "reasons": reasons}
            if reasons else _drift(values, samples, cfg))
    if step["status"] == "skipped" and mode == "on":
        raise DomainError("PREPROCESS_PRECONDITION",
                          f"drift.mode=on だがドリフト補正の条件を満たさない: {step['reasons']}",
                          {"reasons": step["reasons"]})
    return step


# --- 3. 品質の足切り ---------------------------------------------------------

def _sd_mean(block: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """有限正値だけの (SD ddof=1, 平均, 有効数)。"""
    pos = _positive(block)
    n = np.isfinite(pos).sum(axis=0)
    with np.errstate(all="ignore"):
        mean = np.where(n > 0, np.nansum(pos, axis=0) / np.maximum(n, 1), np.nan)
        dev = np.where(np.isfinite(pos), pos - mean, 0.0)
        sd = np.where(n >= 2, np.sqrt((dev ** 2).sum(axis=0) / np.maximum(n - 1, 1)), np.nan)
    return sd, mean, n


def _filters(values, samples, detected, cfg, feature_ids) -> tuple[np.ndarray, dict]:
    """全マスクをドリフト後の値で同時に計算し、論理積で一度に当てる。"""
    roles = np.array([s.role for s in samples])
    bio, qc, blank = (np.nonzero(roles == r)[0] for r in _ROLES)
    n_feat = values.shape[1]
    out, keep = {}, np.ones(n_feat, dtype=bool)

    def record(name, mask, **extra):
        out[name] = {"status": "applied", "dropped": int((~mask).sum()),
                     "dropped_ids": [feature_ids[j] for j in np.nonzero(~mask)[0]], **extra}
        keep[:] &= mask

    def skip(name, reason):
        out[name] = {"status": "skipped", "reason": reason}

    c = cfg["blank_ratio"]
    if c is None:
        skip("blank_ratio", "disabled")
    elif not len(blank):
        skip("blank_ratio", "no_blank")
    else:
        # 平均（spec §4.3・preprocessing.blank_filter と同じ）。0 以下・非有限は欠測として除く。
        _, mean_bio, _ = _sd_mean(values[bio])
        _, mean_blank, _ = _sd_mean(values[blank])
        with np.errstate(all="ignore"):
            record("blank_ratio", ~np.isfinite(mean_blank) | (mean_bio >= c["min"] * mean_blank),
                   min=c["min"], statistic="mean")

    c = cfg["detection"]
    if c is None:
        skip("detection", "disabled")
    elif detected is None:
        skip("detection", "no_detected_mask")
    else:
        is_det = cell_states(values[bio], detected[bio]) == DETECTED
        if c["by"] == "all":
            parts = [np.arange(len(bio))]
        else:
            groups = [samples[i].group for i in bio]
            parts = [np.array([k for k, g in enumerate(groups) if g == lvl])
                     for lvl in dict.fromkeys(groups)]
        rates = np.vstack([is_det[p].mean(axis=0) for p in parts])
        record("detection", (rates >= c["min_rate"]).any(axis=0), by=c["by"], min_rate=c["min_rate"])

    qc_sd, qc_mean, _ = _sd_mean(values[qc]) if len(qc) else (None, None, None)
    c = cfg["qc_rsd"]
    if c is None:
        skip("qc_rsd", "disabled")
    elif len(qc) < 2:
        skip("qc_rsd", "too_few_qc")
    else:
        with np.errstate(all="ignore"):
            rsd = qc_sd / qc_mean
        ok = np.isfinite(rsd)
        record("qc_rsd", ~ok | (rsd <= c["max"]), max=c["max"], not_computable=int((~ok).sum()))

    c = cfg["d_ratio"]
    if c is None:
        skip("d_ratio", "disabled")
    elif len(qc) < 2:
        skip("d_ratio", "too_few_qc")
    elif len(bio) < 2:
        skip("d_ratio", "too_few_samples")
    else:
        bio_sd = _sd_mean(values[bio])[0]
        with np.errstate(all="ignore"):
            ratio = np.where(bio_sd > 0, qc_sd / bio_sd, np.nan)
        ok = np.isfinite(ratio)
        record("d_ratio", ~ok | (ratio <= c["max"]), max=c["max"], not_computable=int((~ok).sum()))

    applied = any(f["status"] == "applied" for f in out.values())
    step = {"step": "filter", "status": "applied" if applied else "skipped",
            "n_features_in": n_feat, "dropped": int((~keep).sum()), "filters": out}
    if not applied:
        step["reason"] = "all_filters_skipped"
    if not keep.any():
        raise DomainError("PREPROCESS_PRECONDITION", "品質の足切りで特徴が 1 つも残らなかった",
                          {"filters": {k: v.get("dropped") for k, v in out.items()}})
    return keep, step


# --- 4. 正規化 ---------------------------------------------------------------

def _normalize_step(values, samples, cfg, quant) -> tuple[np.ndarray, dict]:
    """係数は残った特徴（values の列）だけ・生体＋QC の行から。recipe_matrix._normalize を使う。"""
    from metabolomix.analysis.recipe_matrix import _normalize

    requested = cfg["method"]
    method = requested
    if requested == "auto":
        method = "none" if (quant is not None and quant.normalized) else "pqn"
    if method == "none":
        return values, {"step": "normalize", "status": "skipped", "requested": requested,
                        "method": "none",
                        "reason": "quant_normalized" if requested == "auto" else "method_none"}
    ids = [s.sample_id for s in samples]
    roles = {s.sample_id: s.role for s in samples}
    reference = "qc" if any(s.role == "qc" for s in samples) else "sample"
    out, report = _normalize(values, ids, roles, {"method": method, "reference": reference})
    step = {"step": "normalize", "status": "applied", "requested": requested, "method": method,
            "factors": [report["factors"][i] for i in ids], "factor_samples": ids,
            "factors_set_to_one": len(report["unscaled_samples"]),
            "unscaled_samples": report["unscaled_samples"]}
    if method == "pqn":
        step["reference"] = reference
    return out, step


# --- 5. log2 -----------------------------------------------------------------

def _log2(values: np.ndarray) -> tuple[np.ndarray, dict]:
    finite = np.isfinite(values)
    ok = finite & (values > 0)
    with np.errstate(all="ignore"):
        out = np.where(ok, np.log2(np.where(ok, values, 1.0)), np.nan)
    return out, {"non_positive": int((finite & ~ok).sum()), "non_finite": int((~finite).sum())}


# --- 6. 反復のまとめ ---------------------------------------------------------

def _unit_of(members, unit_id: str) -> UnitInfo:
    """まとめる注入（群・属性が一致する技術反復）から単位を作る。バッチが食い違えば矛盾。"""
    first = members[0]
    seen = sorted({str(s.batch) for s in members})
    if len(seen) > 1:
        raise DomainError("SUBJECT_INCONSISTENT",
                          f"同じ個体 {first.subject_id} の技術反復（群・属性が同じ注入）が別の batch を持つ: {seen}",
                          {"subject_id": first.subject_id, "field": "batch", "values": seen,
                           "samples": [s.sample_id for s in members]})
    orders = [s.injection_order for s in members]
    return UnitInfo(unit_id=unit_id, group=first.group, batch=first.batch,
                    subject_id=first.subject_id,
                    injection_order=None if any(o is None for o in orders) else float(np.mean(orders)),
                    n_injections=len(members), sample_ids=[s.sample_id for s in members],
                    attributes=dict(first.attributes))


def _unit_ids(keys: list[tuple], firsts: list) -> list[str]:
    """まとめた単位の ID（決定的）。

    個体の単位が 1 つなら subject_id。同じ個体から複数の単位ができたら `subject_id|group`、
    それでも重なれば、その個体・群の単位の間で値が違う属性を `|k=v` でキーの順に足す。
    subject_id の無い注入は試料 ID のまま。
    """
    ids: list[str | None] = [None] * len(keys)
    by_subject: dict[str, list[int]] = {}
    for k, key in enumerate(keys):
        if key[0] == "i":
            ids[k] = key[1]
        else:
            by_subject.setdefault(key[1], []).append(k)
    for subject, idx in by_subject.items():
        if len(idx) == 1:
            ids[idx[0]] = subject
            continue
        by_group: dict[str, list[int]] = {}
        for k in idx:
            by_group.setdefault(str(firsts[k].group), []).append(k)
        for group, sub in by_group.items():
            if len(sub) == 1:
                ids[sub[0]] = f"{subject}|{group}"
                continue
            names = sorted({n for k in sub for n in firsts[k].attributes})
            differing = [n for n in names
                         if len({firsts[k].attributes.get(n) for k in sub}) > 1]
            for k in sub:
                extra = "".join(f"|{n}={firsts[k].attributes.get(n, '')}" for n in differing)
                ids[k] = f"{subject}|{group}{extra}"
    return ids


def _split_warning(keys, firsts, injection_columns) -> dict | None:
    """同じ (個体, 群) が複数の単位に分かれたとき、分けた属性の列を名指しして注意にする。"""
    pairs: dict[tuple, list[int]] = {}
    for k, key in enumerate(keys):
        if key[0] == "s":
            pairs.setdefault((key[1], key[2]), []).append(k)
    split = {pair: idx for pair, idx in pairs.items() if len(idx) > 1}
    if not split:
        return None
    columns: set[str] = set()
    for idx in split.values():
        names = sorted({n for k in idx for n in firsts[k].attributes})
        columns.update(n for n in names if len({firsts[k].attributes.get(n) for k in idx}) > 1)
    subjects = [f"{subj}|{grp}" for subj, grp in split]
    shown = subjects[:20]
    more = f" ほか {len(subjects) - len(shown)} 件" if len(subjects) > len(shown) else ""
    cols = sorted(columns)
    return {"code": "units_split_by_attribute",
            "message": f"同じ個体・群なのに試料の表の列（{', '.join(cols)}）の値が違うため、{len(subjects)} 件が"
                       f"複数の単位に分かれました（{', '.join(shown)}{more}）。注入ごとに違う列（取得日・注入番号など）"
                       "なら、前処理の injection_columns にその列名を入れると 1 個体にまとまります。",
            "attributes": cols, "n_subjects": len(subjects), "subjects": subjects[:20],
            "suggestion": {"preprocess.injection_columns": cols}}


def _collapse(log2, det, samples, enabled: bool, injection_columns=()):
    """個体 × 群 × 属性が同じ注入を技術反復として log2 の上で平均する（spec §6 の対応のある計画のため、
    同じ個体でも群か属性が違えば別の単位。単位は subject_id を持ち続ける）。"""
    with_subject = sum(s.subject_id is not None for s in samples)
    source = ("subject_id" if with_subject == len(samples)
              else "injection" if with_subject == 0 else "mixed")
    groups: dict[tuple, list[int]] = {}
    for k, s in enumerate(samples):
        if enabled and s.subject_id is not None:
            key = ("s", s.subject_id, s.group,
                   tuple(sorted((k, v) for k, v in s.attributes.items()
                                if k not in injection_columns)))
        else:
            key = ("i", s.sample_id)
        groups.setdefault(key, []).append(k)
    keys = list(groups)
    firsts = [samples[groups[key][0]] for key in keys]
    unit_ids = _unit_ids(keys, firsts)
    units = [_unit_of([samples[k] for k in groups[key]], uid) for key, uid in zip(keys, unit_ids)]
    per_subject: dict[str, int] = {}
    for key in keys:
        if key[0] == "s":
            per_subject[key[1]] = per_subject.get(key[1], 0) + 1
    seen: dict[str, UnitInfo] = {}
    for u in units:
        if u.unit_id in seen:
            raise DomainError("SUBJECT_INCONSISTENT",
                              f"個体 ID {u.unit_id} が別の単位の ID と衝突する（subject_id と試料 ID が同じ）",
                              {"unit_id": u.unit_id,
                               "samples": seen[u.unit_id].sample_ids + u.sample_ids})
        seen[u.unit_id] = u
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out = np.vstack([np.nanmean(log2[idx], axis=0) if len(idx) > 1 else log2[idx[0]]
                         for idx in groups.values()])
    det_units = None if det is None else np.vstack([det[idx].any(axis=0) for idx in groups.values()])
    if not enabled:
        source = "injection"
    step = {"step": "collapse", "status": "applied" if enabled else "skipped",
            "subject_source": source, "n_injections": len(samples), "n_units": len(units),
            "collapsed_units": sum(u.n_injections > 1 for u in units),
            "subjects_with_multiple_units": sum(n > 1 for n in per_subject.values())}
    if not enabled:
        step["reason"] = "disabled"
    warn = _split_warning(keys, firsts, injection_columns) if enabled else None
    return out, det_units, units, source, step, warn


# --- 入口 --------------------------------------------------------------------

def preprocess(study: StudyData, config: dict) -> PreparedMatrix | NeedsInput:
    """spec §4 の 1〜6 を順に当てる。補完はしない。"""
    cfg = normalize_preprocess_config(config)
    steps = []
    keep_rows, step = _select(study)
    steps.append(step)
    samples = [study.samples[i] for i in keep_rows]
    values = np.array(np.asarray(study.values, dtype=float)[keep_rows], copy=True)
    detected = None if study.detected is None else np.asarray(study.detected, bool)[keep_rows]

    step = _drift_step(values, samples, cfg["drift"], study.quant)
    if isinstance(step, NeedsInput):
        return step
    steps.append(step)

    feature_ids = list(study.feature_ids)
    keep_cols, step = _filters(values, samples, detected, cfg["filters"], feature_ids)
    steps.append(step)
    cols = np.nonzero(keep_cols)[0]
    feature_ids = [feature_ids[j] for j in cols]
    # ブランクは足切りにだけ使い、行列には入れない
    rows = [k for k, s in enumerate(samples) if s.role != "blank"]
    samples = [samples[k] for k in rows]
    values = values[np.ix_(rows, cols)]
    detected = None if detected is None else detected[np.ix_(rows, cols)]

    values, step = _normalize_step(values, samples, cfg["normalize"], study.quant)
    steps.append(step)
    normalize_method = step["method"]

    log2, missing = _log2(values)
    bio = [k for k, s in enumerate(samples) if s.role == "sample"]
    qc = [k for k, s in enumerate(samples) if s.role == "qc"]
    _, qc_missing = _log2(values[qc])
    bio_missing = {k: missing[k] - qc_missing[k] for k in missing}
    steps.append({"step": "log2", "status": "applied", "missing_cells": bio_missing,
                  "qc_missing_cells": qc_missing})

    bio_samples = [samples[k] for k in bio]
    det_bio = None if detected is None else (cell_states(values[bio], detected[bio]) == DETECTED)
    unit_log2, det_units, units, source, step, split_warning = _collapse(
        log2[bio], det_bio, bio_samples, cfg["collapse_replicates"],
        frozenset(cfg["injection_columns"]))
    steps.append(step)

    notes = []
    if split_warning is not None:
        notes.append(split_warning)
    if not any(s.subject_id is not None for s in bio_samples):
        notes += replicate_name_warnings([s.sample_id for s in bio_samples])
    kept = set(feature_ids)
    return PreparedMatrix(
        matrix_id=_matrix_id(study, cfg, normalize_method), log2=unit_log2, units=units, feature_ids=feature_ids,
        detected_units=det_units, qc_log2=log2[qc] if qc else np.empty((0, len(cols))),
        qc_ids=[samples[k].sample_id for k in qc], steps=steps, config=cfg, quant=study.quant,
        dataset_hash=study.dataset_hash, subject_source=source, warnings=notes,
        feature_classes={f: c for f, c in study.feature_classes.items() if f in kept},
        annotated_features=frozenset(f for f in study.annotated_features if f in kept),
        sample_values=np.array(values[bio], copy=True),
        sample_ids=[samples[k].sample_id for k in bio],
        sample_detected=None if detected is None else np.array(detected[bio], copy=True))


def impute_for_pca(log2: np.ndarray, *, method: str = "half_min",
                   reference_rows: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    """PCA など完全な行列が要る計算のためだけの補完（検定には使わない。spec §4 の 7）。

    `half_min`: 特徴ごとに `reference_rows`（既定は全行）の最小の log2 − 1（＝半分の値）。
    参照の行に観測が無い特徴は全行の最小、それも無ければ NaN のまま残す。
    `knn`: `recipe_matrix.knn_impute(k=10, rowmax=0.5, colmax=0.8)` を log2 の上で。
    """
    data = np.array(np.asarray(log2, dtype=float), copy=True)
    missing = ~np.isfinite(data)
    if method == "half_min":
        ref = data if reference_rows is None else data[np.asarray(reference_rows)]
        fill = _nanmin(ref)
        fill = np.where(np.isfinite(fill), fill, _nanmin(data)) - 1.0
        out = np.where(missing, fill[None, :], data)
    elif method == "knn":
        from metabolomix.analysis.recipe_matrix import knn_impute
        out = knn_impute(data.T, k=10, rowmax=0.5, colmax=0.8).T
    else:
        raise DomainError("PREPROCESS_CONFIG_INVALID", f"未知の補完法: {method}",
                          {"field": "impute.method", "value": repr(method)})
    imputed = int((missing & np.isfinite(out)).sum())
    return out, {"method": method, "imputed_cells": imputed}
