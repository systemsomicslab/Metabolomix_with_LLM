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

#: 既定の dict に載せない任意のキー（path → 既定値）。既定値のときは設定に書かない（書けば matrix_id が
#: 変わり、既定どおりの前処理の ID がこのキーを足す前と食い違う）。
_OPTIONAL = {("drift", "loess_degree"): 1, ("normalize", "features"): None,
             ("normalize", "exclude_features"): None}

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
        optional = {key[-1]: value for key, value in _OPTIONAL.items() if key[:-1] == path}
        unknown = sorted(set(raw) - set(default) - set(optional))
        if unknown:
            raise _invalid(f"未知のキー: {unknown}", path + (unknown[0],), raw[unknown[0]])
        out = {k: _merge(v, raw.get(k, copy.deepcopy(v)), path + (k,)) for k, v in default.items()}
        for k, fallback in optional.items():
            if k in raw:
                value = _optional_value(path + (k,), raw[k])
                if value != fallback:
                    out[k] = value
        return out
    if path in _ENUMS:
        if raw not in _ENUMS[path]:
            raise _invalid(f"{'.'.join(path)} は {list(_ENUMS[path])} のどれか", path, raw)
        return raw
    if path == ("injection_columns",):
        if (not isinstance(raw, list)
                or any(not isinstance(c, str) or not c.strip() for c in raw)):
            raise _invalid("injection_columns は空でない文字列のリストで指定する", path, raw)
        return sorted({c.strip() for c in raw})
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


def _optional_value(path: tuple, raw):
    """任意のキーの値を検める（既定値なら呼び出し側が落とす）。"""
    if path == ("drift", "loess_degree"):
        if isinstance(raw, bool) or raw not in (1, 2):
            raise _invalid("drift.loess_degree は 1（局所線形。既定）か 2（局所 2 次）", path, raw)
        return int(raw)
    # normalize.features / exclude_features: 特徴 ID のリスト（None は指定なし）。空のリストは指定なしと同じ
    # （features の空は絞らない＝残った特徴すべて。normalize_preprocess_config が notes に書く）。
    if raw is None:
        return None
    if (not isinstance(raw, list)
            or any(not isinstance(f, str) or not f.strip() for f in raw)):
        raise _invalid(f"{'.'.join(path)} は空でない文字列（特徴 ID）のリストで指定する", path, raw)
    return sorted({f.strip() for f in raw}) or None


def _ignored(code: str, message: str, **extra) -> dict:
    return {"code": code, **extra, "message": message}


def _drop_inert(cfg: dict, raw, notes: list[dict]) -> None:
    """書いても効かないキーを設定から捨て、`notes` に書く（捨てたキーは matrix_id にも rule_hash にも入らない）。

    - `normalize.method="none"` の係数の特徴（`features` / `exclude_features`）: 係数を計算しない。
    - 空の `normalize.features`: 絞らない（指定なしと同じ。`_optional_value` が落としている）。
    - `drift.mode="off"` の `drift.loess_degree`: 補正しない。
    いずれも任意のキー（`_OPTIONAL`）だけで、既定の dict のキー（`drift.span` など）は従来どおり残す。
    """
    norm, drift = cfg["normalize"], cfg["drift"]
    raw_norm = raw.get("normalize") if isinstance(raw, dict) else None
    empty_features = isinstance(raw_norm, dict) and raw_norm.get("features") == []
    if norm["method"] == "none":
        given = [k for k in ("features", "exclude_features") if k in norm]
        for key in given:
            del norm[key]
        if given or empty_features:
            notes.append(_ignored(
                "normalization_features_ignored",
                "normalize=none は正規化の係数を計算しないので、normalization_features / "
                "normalization_exclude_features（normalize.features / exclude_features）は効かず、"
                "使っていない（行列と matrix_id は指定なしと同じ）。", reason="method_none"))
    elif empty_features:
        notes.append(_ignored(
            "normalization_features_empty",
            "normalization_features（normalize.features）が空のリストなので、指定なしとして扱った"
            "（係数は足切りの後に残ったすべての特徴から計算する）。内部標準などに絞るなら特徴 ID を渡す。"))
    if drift["mode"] == "off" and "loess_degree" in drift:
        del drift["loess_degree"]
        notes.append(_ignored(
            "loess_degree_ignored",
            "drift_mode=off はドリフトを補正しないので、loess_degree（drift.loess_degree）は効かず、"
            "使っていない（行列と matrix_id は指定なしと同じ）。", reason="drift_off"))


def normalize_preprocess_config(raw: dict | None, *, notes: list[dict] | None = None) -> dict:
    """既定を埋めた前処理設定を返す。未知のキー・不正な値は PREPROCESS_CONFIG_INVALID。

    書いても効かない任意のキー（`normalize.method="none"` の `features` / `exclude_features`、空の
    `features`、`drift.mode="off"` の `loess_degree`）は誤りにせず捨て、`notes`（渡せば）に
    `{"code", "reason"?, "message"}` を足す（`_drop_inert`）。`auto` が換算済みのデータで `none` に
    決まったときはデータで決まるので、前処理が `factor_features_ignored` を残す。
    """
    cfg = _merge(PREPROCESS_DEFAULTS, {} if raw is None else raw, ())
    _drop_inert(cfg, raw, [] if notes is None else notes)
    return cfg


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
    #: 注入順の出所（まとめた注入で 1 つに決まらなければ "mixed"。値が無ければ None）。
    injection_order_source: str | None = None
    #: 注入順が MS-DIAL の既定（ファイルの順）かもしれない注入を含む（`injection_order_may_be_file_order`）。
    injection_order_may_be_file_order: bool = False


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


def _matrix_id(study: StudyData, config: dict, normalize_method: str, extras: dict | None = None) -> str:
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
    body = {"dataset_hash": source, "sample_table_hash": canonical_hash(table),
            "config": config, "normalize_method": normalize_method,
            "feature_ids_hash": feature_hash,
            "feature_merges": [list(g) for g in study.feature_merges]}
    # 既定の判断を変えた事情（drift=auto を止めた理由など）。あるときだけ足す（無ければ ID は前と同じ）
    if extras:
        body["extras"] = extras
    return "mx_" + canonical_hash(body)[:16]


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


#: LOESS の近傍の最小数（局所線形）。statsmodels は近傍数を int(frac·n) で決め、最も遠い近傍の重みは 0。
#: 近傍が 3 だと等間隔の QC で等距離の 2 点がともに重み 0 になり、xvals 経路の当てはめが NaN になる。
_MIN_LOESS_NEIGHBOURS = 4


def _min_neighbours(degree: int) -> int:
    """局所 d 次の当てはめに要る近傍の数。重み 0 になりうる両端の 2 点を除いて d + 1 点が残る数。"""
    return _MIN_LOESS_NEIGHBOURS + (int(degree) - 1)


def _neighbours(span: float, n: int) -> int:
    return int(span * n + 1e-10)    # statsmodels の lowess と同じ数え方


def _local_poly(x, y, xvals, span: float, degree: int) -> np.ndarray:
    """局所多項式回帰（LOESS。頑健化の反復なし）を `xvals` で評価する。numpy だけで書く。

    近傍は statsmodels の lowess と同じ: 近い順に k = int(span·n) 点（2 ≤ k ≤ n。距離が同じなら左を先に）、
    半径はその最も遠い距離、重みは tricube (1 − (d/半径)³)³。近傍で重みが正の点が d + 1 未満、または
    計画行列の階数が足りなければ NaN。degree=1 は statsmodels の `lowess(it=0, delta=0, xvals=...)` と
    数値誤差の範囲で一致する（tests/test_stats_preprocess.py が縛る）。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xvals = np.asarray(xvals, dtype=float)
    n = len(x)
    k = min(max(_neighbours(span, n), 2), n)
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]
    out = np.full(len(xvals), np.nan)
    for t, xv in enumerate(xvals):
        dist = np.abs(xs - xv)
        near = np.sort(np.argsort(dist, kind="stable")[:k])
        radius = float(dist[near].max())
        if radius <= 0:
            continue
        weight = np.clip(1.0 - (dist[near] / radius) ** 3, 0.0, None) ** 3
        ok = weight > 1e-12
        if int(ok.sum()) < degree + 1:
            continue
        design = np.vander(xs[near][ok] - xv, degree + 1, increasing=True)
        root = np.sqrt(weight[ok])
        coef, _res, rank, _sv = np.linalg.lstsq(design * root[:, None], ys[near][ok] * root, rcond=None)
        if rank < degree + 1:
            continue
        out[t] = coef[0]
    return out


def _loess(y, x, xvals, span: float, degree: int) -> np.ndarray:
    """局所線形は statsmodels の lowess（従来どおり）、局所 2 次は `_local_poly`。"""
    if degree == 1:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        return lowess(y, x, frac=span, it=0, delta=0.0, xvals=xvals)
    return _local_poly(x, y, xvals, span, degree)


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
        elif _neighbours(cfg["span"], len(qc)) < _min_neighbours(cfg.get("loess_degree", 1)):
            reasons.append(f"too_few_qc_for_span:{b}")
        elif bio and not (min(qc) <= min(bio) and max(qc) >= max(bio)):
            reasons.append(f"qc_not_interspersed:{b}")
    return reasons


def _drift(values: np.ndarray, samples, cfg: dict) -> dict:
    """バッチごとの LOESS（生の尺度）→ 特徴ごとのバッチ水準合わせ。

    全バッチで補正できた特徴だけを補正し（一部のバッチだけの補正はしない）、水準もその特徴だけ
    合わせる。補正できた特徴が 1 つも無いバッチがあれば全体を飛ばす。成功したときだけ values を書き換える。
    """
    span, min_qc = cfg["span"], cfg["min_qc_per_batch"]
    degree = cfg.get("loess_degree", 1)
    min_near = _min_neighbours(degree)
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
            if ok.sum() < min_qc or _neighbours(span, int(ok.sum())) < min_near:
                too_few[b] += 1
                continue
            fit = _loess(y[ok], x_qc[ok], x_all, span, degree)
            if not np.all(np.isfinite(fit)) or np.any(fit <= 0):
                fit_invalid[b] += 1
                continue
            work[rows, j] = values[rows, j] / fit * float(np.median(y[ok]))
            in_batch[j] = True
        if not in_batch.any():
            empty.append(f"fit_invalid:{b}" if fit_invalid[b] else f"no_correctable_feature:{b}")
        corrected &= in_batch
    record = {"step": "drift", "method": "loess", "span": span, "loess_degree": degree, "batches": batches,
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


#: drift=auto が補正を止めた理由（補正の条件は満たしている。drift=on で呼び直せば補正する）。
_AUTO_HOLD_MESSAGES = {
    "quant_normalized": "内部標準などで換算済み（quant_state.normalized）。内部標準の比はドリフトの一部を"
                        "すでに打ち消しているので、QC の LOESS を重ねると換算済みの値の尺度を動かしうる",
    "injection_order_may_be_file_order": "注入順が MS-DIAL の既定（ファイルの順）かもしれない（またはシートの宣言で"
                                         "ファイル名の順）。実際の測定順でなければ、LOESS は誤った順でドリフトを当てはめる",
}


def _order_uncertain(samples) -> bool:
    """ドリフト補正に使う注入（生体・QC）の注入順が確かでない（MS-DIAL の既定かもしれない・ファイル名の順）。"""
    return any(s.injection_order_may_be_file_order or s.injection_order_source == "filename"
               for s in samples if s.role in ("sample", "qc"))


def _auto_hold(samples, quant) -> list[str]:
    """drift=auto が補正しない理由（補正の条件を満たすときだけ問う）。"""
    held = []
    if quant is not None and quant.normalized:
        held.append("quant_normalized")
    if _order_uncertain(samples):
        held.append("injection_order_may_be_file_order")
    return held


def _drift_step(values, samples, cfg, quant) -> dict:
    """ドリフト補正の段。mode=auto は条件を満たしても、換算済み・注入順が不確かなら補正しない
    （尋ねずに off で計算し、理由と drift=on の選択肢を `rerun_with` に残す）。mode=on は条件を満たせば
    必ず補正し、満たさなければ PREPROCESS_PRECONDITION。"""
    mode = cfg["mode"]
    if mode == "off":
        return {"step": "drift", "status": "skipped", "reason": "mode_off"}
    reasons = _drift_eligibility(samples, cfg)
    if not reasons and mode == "auto":
        held = _auto_hold(samples, quant)
        if held:
            return {"step": "drift", "status": "skipped", "reason": held[0], "reasons": held,
                    "rerun_with": {"preprocess.drift.mode": "on"}}
    step = ({"step": "drift", "status": "skipped", "reason": reasons[0], "reasons": reasons}
            if reasons else _drift(values, samples, cfg))
    if step["status"] == "skipped" and mode == "on":
        raise DomainError("PREPROCESS_PRECONDITION",
                          f"drift.mode=on だがドリフト補正の条件を満たさない: {step['reasons']}",
                          {"reasons": step["reasons"]})
    return step


def _drift_hold_warning(step: dict) -> dict | None:
    """drift=auto が補正を止めたときの注意（理由と、drift=on で呼び直す選択肢）。"""
    if step.get("rerun_with") is None:
        return None
    reasons = list(step.get("reasons") or [])
    return {"code": "drift_not_applied_by_default", "reasons": reasons,
            "message": "ドリフト補正は条件を満たすが、既定（drift=auto）では当てなかった: "
                       + "。".join(_AUTO_HOLD_MESSAGES[r] for r in reasons if r in _AUTO_HOLD_MESSAGES)
                       + "。補正するなら drift_mode=\"on\" で呼び直す。",
            "suggestion": dict(step["rerun_with"])}


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
        # 平均（spec §4.3）。0 以下・非有限は欠測として除く。
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

def _factor_columns(cfg, kept_ids: list[str], all_ids: list[str],
                    all_values: np.ndarray) -> tuple[np.ndarray | None, dict]:
    """正規化の係数を計算する列（試料 × 特徴）と、その記録。指定が無ければ (None, {})。

    既定は残った特徴（品質の足切りの後）すべて。`features` を指せばその特徴だけ（足切りで落ちた特徴も
    使える。内部標準は生体の分散が小さく D-ratio で落ちうるため。値はドリフト補正の後）。
    `exclude_features` はその中から外す。係数を計算する列から外した特徴も行列には残る。
    """
    chosen, excluded = cfg.get("features"), cfg.get("exclude_features")
    if chosen is None and excluded is None:
        return None, {}
    position = {fid: j for j, fid in enumerate(all_ids)}
    unknown = sorted(f for f in (chosen or []) + (excluded or []) if f not in position)
    if unknown:
        field_name = ("normalize.features" if chosen and set(unknown) & set(chosen)
                      else "normalize.exclude_features")
        raise DomainError("PREPROCESS_CONFIG_INVALID",
                          f"{field_name} に存在しない特徴 ID がある: {unknown[:20]}",
                          {"field": field_name, "unknown_features": unknown})
    base = list(chosen) if chosen is not None else list(kept_ids)
    drop = set(excluded or [])
    use = [f for f in base if f not in drop]
    if not use:
        raise DomainError("PREPROCESS_CONFIG_INVALID", "正規化の係数を計算する特徴が 1 つも残らない",
                          {"field": "normalize.exclude_features" if excluded else "normalize.features"})
    kept = set(kept_ids)
    record = {"factor_features": len(use),
              "factor_feature_basis": "features" if chosen is not None else "kept_features"}
    if excluded is not None:
        record["excluded_from_factors"] = len(drop & set(base))
    outside = [f for f in use if f not in kept]
    if outside:
        record["factor_features_filtered_out"] = outside
    return all_values[:, [position[f] for f in use]], record


def _normalize_step(values, samples, cfg, quant, factor_values=None,
                    factor_record=None) -> tuple[np.ndarray, dict]:
    """係数は残った特徴（values の列）だけ・生体＋QC の行から。recipe_matrix._normalize を使う。

    `factor_values` を渡せば係数はその列から計算し、`values` のすべての列に当てる（特徴の部分集合の正規化）。
    """
    from metabolomix.analysis.recipe_matrix import _normalize

    requested = cfg["method"]
    method = requested
    if requested == "auto":
        method = "none" if (quant is not None and quant.normalized) else "pqn"
    if method == "none":
        step = {"step": "normalize", "status": "skipped", "requested": requested, "method": "none",
                "reason": "quant_normalized" if requested == "auto" else "method_none"}
        if factor_record:
            # auto が換算済みのデータで none に決まった（明示の none は normalize_preprocess_config がキーごと
            # 捨てている）: 係数の特徴の指定は効かない（行列も ID も指定なしと同じ）
            step["factor_features_ignored"] = True
        return values, step
    ids = [s.sample_id for s in samples]
    roles = {s.sample_id: s.role for s in samples}
    reference = "qc" if any(s.role == "qc" for s in samples) else "sample"
    basis = values if factor_values is None else factor_values
    _scaled, report = _normalize(basis, ids, roles, {"method": method, "reference": reference})
    factors = np.array([report["factors"][i] for i in ids], dtype=float)
    out = values / factors[:, None]
    step = {"step": "normalize", "status": "applied", "requested": requested, "method": method,
            "factors": [report["factors"][i] for i in ids], "factor_samples": ids,
            "factors_set_to_one": len(report["unscaled_samples"]),
            "unscaled_samples": report["unscaled_samples"], **(factor_record or {})}
    if method == "pqn":
        step["reference"] = reference
        if reference == "sample":
            # QC が無いときは参照の取り方を明示する（生体試料・全群の中央値で、群の組成差を吸収しうる）。
            # QC があるときの段の形（reference="qc" だけ）は変えない
            step["reference_definition"] = "median_of_samples_all_groups"
    return out, step


def _pqn_reference_warning(step: dict) -> dict | None:
    """QC が無い PQN の参照（生体試料の特徴ごとの中央値）を明示する注意。"""
    if step.get("method") != "pqn" or step.get("reference") != "sample":
        return None
    return {"code": "pqn_reference_is_samples",
            "message": "QC が無いので、PQN の参照スペクトルは生体試料（全群）の特徴ごとの中央値。"
                       "群の間で全体の組成が大きく違うと、係数がその差の一部を吸収しうる"
                       "（係数を計算する特徴を絞るなら normalization_features）。"}



# --- 5. log2 -----------------------------------------------------------------

def _log2(values: np.ndarray) -> tuple[np.ndarray, dict]:
    finite = np.isfinite(values)
    ok = finite & (values > 0)
    with np.errstate(all="ignore"):
        out = np.where(ok, np.log2(np.where(ok, values, 1.0)), np.nan)
    return out, {"non_positive": int((finite & ~ok).sum()), "non_finite": int((~finite).sum())}


# --- 6. 反復のまとめ ---------------------------------------------------------

def _unit_of(members, unit_id: str, drop=frozenset()) -> UnitInfo:
    """まとめる注入（群・属性が一致する技術反復）から単位を作る。バッチが食い違えば矛盾。"""
    first = members[0]
    seen = sorted({str(s.batch) for s in members})
    sources = {s.injection_order_source for s in members if s.injection_order is not None}
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
                    attributes={k: v for k, v in first.attributes.items() if k not in drop},
                    injection_order_source=(None if not sources else sources.pop() if len(sources) == 1
                                            else "mixed"),
                    injection_order_may_be_file_order=any(s.injection_order_may_be_file_order
                                                          for s in members))


def _unit_ids(keys: list[tuple], firsts: list, ignore=frozenset()) -> list[str]:
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
            names = sorted({n for k in sub for n in firsts[k].attributes if n not in ignore})
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
        names = sorted({n for k in idx for n in firsts[k].attributes
                        if n not in injection_columns})
        columns.update(n for n in names if len({firsts[k].attributes.get(n) for k in idx}) > 1)
    subjects = [f"{subj}|{grp}" for subj, grp in split]
    shown = subjects[:20]
    more = f" ほか {len(subjects) - len(shown)} 件" if len(subjects) > len(shown) else ""
    cols = sorted(columns)
    return {"code": "units_split_by_attribute",
            "message": f"同じ個体・群なのに試料の表の列（{', '.join(cols)}）の値が違うため、{len(subjects)} 件が"
                       f"複数の単位に分かれました（{', '.join(shown)}{more}）。注入ごとに違う列（取得日・注入番号など）"
                       "なら、前処理の injection_columns にその列名を入れると 1 個体にまとまります。",
            "attributes": cols, "n_subjects": len(subjects), "subjects": shown,
            "subjects_n_more": len(subjects) - len(shown),
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
    drop = frozenset(injection_columns) if enabled else frozenset()
    unit_ids = _unit_ids(keys, firsts, drop)
    units = [_unit_of([samples[k] for k in groups[key]], uid, drop)
             for key, uid in zip(keys, unit_ids)]
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

def _batch_candidate_warning(study: StudyData, samples) -> dict | None:
    """バッチが不明な試料の名前に 8 桁が 2 種類以上あるとき、候補として返す（統計には当てていない）。"""
    found = {s.sample_id: study.batch_candidates[s.sample_id]["value"] for s in samples
             if s.batch is None and s.sample_id in study.batch_candidates}
    values = sorted(set(found.values()))
    if len(values) < 2:
        return None
    return {"code": "batch_candidates", "source": "filename_date", "n_samples": len(found),
            "values": values[:20], "values_n_more": max(0, len(values) - 20),
            "message": "試料名の 8 桁の数字（取得日らしい）が複数ある。バッチの候補として返すだけで、ドリフト補正・"
                       "反復のまとめ・交絡の判定には当てていない。バッチなら試料の表の batch 列に入れる"
                       "（dataset_set_sample_metadata）。"}


def preprocess(study: StudyData, config: dict) -> PreparedMatrix:
    """spec §4 の 1〜6 を順に当てる。補完はしない。"""
    from metabolomix.analysis.sample_manifest import file_order_warning
    config_notes: list[dict] = []
    cfg = normalize_preprocess_config(config, notes=config_notes)
    steps = []
    keep_rows, step = _select(study)
    steps.append(step)
    samples = [study.samples[i] for i in keep_rows]
    values = np.array(np.asarray(study.values, dtype=float)[keep_rows], copy=True)
    detected = None if study.detected is None else np.asarray(study.detected, bool)[keep_rows]

    notes = list(config_notes)     # 効かない引数を捨てた注記（normalize_preprocess_config）
    flagged = [s for s in samples if s.role in ("sample", "qc") and s.injection_order_may_be_file_order]
    if flagged:
        notes.append(file_order_warning(len(flagged)))
    step = _drift_step(values, samples, cfg["drift"], study.quant)
    steps.append(step)
    hold = _drift_hold_warning(step)
    if hold is not None:
        notes.append(hold)
    extras = {"drift_held": step["reasons"]} if hold is not None else None

    all_ids = list(study.feature_ids)
    feature_ids = list(all_ids)
    keep_cols, step = _filters(values, samples, detected, cfg["filters"], feature_ids)
    steps.append(step)
    cols = np.nonzero(keep_cols)[0]
    feature_ids = [feature_ids[j] for j in cols]
    # ブランクは足切りにだけ使い、行列には入れない
    rows = [k for k, s in enumerate(samples) if s.role != "blank"]
    samples = [samples[k] for k in rows]
    drifted = values[rows]
    values = values[np.ix_(rows, cols)]
    detected = None if detected is None else detected[np.ix_(rows, cols)]

    factor_values, factor_record = _factor_columns(cfg["normalize"], feature_ids, all_ids, drifted)
    values, step = _normalize_step(values, samples, cfg["normalize"], study.quant,
                                   factor_values=factor_values, factor_record=factor_record)
    steps.append(step)
    normalize_method = step["method"]
    pqn_note = _pqn_reference_warning(step)
    if pqn_note is not None:
        notes.append(pqn_note)
    id_cfg = cfg
    if step.get("factor_features_ignored"):
        notes.append({"code": "normalization_features_ignored", "reason": "quant_normalized",
                      "message": "normalize=auto は換算済みのデータ（quant_normalized）で正規化しないので、"
                                 "normalization_features / normalization_exclude_features は効いていない"
                                 "（行列は指定なしと同じ）。"})
        # 効かなかった指定で matrix_id を変えない（同じ行列は同じ ID）
        id_cfg = {**cfg, "normalize": {"method": cfg["normalize"]["method"]}}

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

    if split_warning is not None:
        notes.append(split_warning)
    if not any(s.subject_id is not None for s in bio_samples):
        notes += replicate_name_warnings([s.sample_id for s in bio_samples])
    batch_note = _batch_candidate_warning(study, samples)
    if batch_note is not None:
        notes.append(batch_note)
    kept = set(feature_ids)
    return PreparedMatrix(
        matrix_id=_matrix_id(study, id_cfg, normalize_method, extras), log2=unit_log2, units=units,
        feature_ids=feature_ids,
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
