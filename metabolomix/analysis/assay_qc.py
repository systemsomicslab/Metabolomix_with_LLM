"""固定母集団のQCと、解析filterの分離（spec §9）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

**QCと解析filterは別の軸。** 同じ mask を兼用すると、閾値で feature を落とすほど
バッチが合格に近づく——落とした分だけ fail が分母から消えるからで、これは
「品質が上がった」ではなく「測らなかった」。ここでは:

- 評価集合（母集団）を filter **前**に固定し、`population` として返す。処理後QCは
  同じ集合を受け取って再利用する（`evaluate_qc(..., population=前回のもの)`）。
- QCで落ちた feature は `eligibility_suggestion` として**提案するだけ**。適用は
  呼び出し側が解析行列の eligibility に対して行い、QC結果は再集計しない。
- 評価不能（U）を分母から消さない。`P/N≥t` なら pass、`(P+U)/N<t` なら fail、
  それ以外は not_evaluable（`aggregate_counts`）。「測れなかった」を
  「合格」にも「不合格」にも勝手に倒さない。

**前提が無い metric は not_evaluable。** 標準品しか無いバッチに pooled QC の合格を
付けない。検出 mask が無いのを「全部検出」と読まない。補完後の行列で QC-RSD を
合格にしない（補完は分散を縮める＝RSD を良く見せる）。

`+inf` は JSON に出さない。blank が 0 で試料に信号がある場合の blank_fold は
数学的には +inf だが、`Infinity` は JSON の数値ではない。`value=null,
special_value="positive_infinity"` として持ち、閾値比較では「どんな有限値より
大きい」として扱う——None に潰すと、最も明確な合格が評価不能に化ける。
"""
from __future__ import annotations

import numpy as np

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "QC_SCHEMA",
    "aggregate_counts",
    "evaluate_qc",
]

QC_SCHEMA = "assay-qc.v1"

_POLICY_INVALID = "QC_POLICY_INVALID"

#: spec §9 の表にある6 metric だけ。ここに無い名前は推測で計算しない。
_METRICS = ("pooled_qc_rsd", "blank_fold", "detection_rate",
            "standard_rt_error", "standard_mass_error_ppm",
            "internal_standard_valid_fraction")

#: pooled QC-RSD に必要な同pool・同batchの有効QC本数（spec §9 の表）。
_MIN_QC_FOR_RSD = 3

_POSITIVE_INFINITY = "positive_infinity"


def aggregate_counts(p: int, f: int, u: int, threshold: float) -> str:
    """集合判定。評価不能（u）を分母から除かない。

    `P/N≥t` なら pass。`(P+U)/N<t` なら——残り全部が pass でも届かないので
    ——fail。それ以外は「あと U 次第」なので not_evaluable。N=0 は not_evaluable
    （何も測っていない集合に合格を出さない）。
    """
    n = p + f + u
    if not n:
        return "not_evaluable"
    if p / n >= threshold:
        return "pass"
    return "fail" if (p + u) / n < threshold else "not_evaluable"


# ---------- 比較 ----------

def _compare(value, special_value, threshold: dict) -> str:
    """1要素の値を閾値と比べる。値が無ければ not_evaluable。"""
    operator = threshold.get("operator")
    limit = threshold.get("value")
    if special_value == _POSITIVE_INFINITY:
        # 有限のどんな閾値より大きい。`>=` `>` は pass、`<=` `<` は fail。
        return "pass" if operator in (">=", ">") else "fail"
    if value is None or not np.isfinite(value):
        return "not_evaluable"
    checks = {"<=": value <= limit, "<": value < limit,
              ">=": value >= limit, ">": value > limit,
              "==": value == limit}
    if operator not in checks:
        raise DomainError(
            _POLICY_INVALID, f"未知の比較演算子です: {operator!r}",
            {"operator": operator, "supported": sorted(checks)})
    return "pass" if checks[operator] else "fail"


def _element(element_id: str, value, *, special_value=None, reason=None,
             threshold: dict | None = None, absolute: bool = False) -> dict:
    """1要素の結果。判定に使う値と表示する値を分けられるようにする。"""
    judged = value
    if absolute and value is not None and np.isfinite(value):
        judged = abs(value)
    status = ("not_evaluable" if threshold is None
              else _compare(judged, special_value, threshold))
    return {"element_id": element_id,
            "value": None if value is None or not np.isfinite(value) else float(value),
            "special_value": special_value,
            "status": status,
            "reason": reason if reason else (None if status != "not_evaluable"
                                             else "value_not_available")}


# ---------- 母集団 ----------

def _fix_population(matrix: dict, entry: dict) -> dict:
    """metric 1件の評価集合を filter **前**の状態から作る。

    `scope="targets"` でも feature_id の解決（target → feature）は上流
    （binding）が済ませている前提で、ここは `population` に載って来たものを
    そのまま使う。解決できていない target を名前で拾いに行かない。
    """
    feature_ids = list(matrix["feature_ids"])
    assay_ids = list(matrix["assay_ids"])
    population = {"feature_ids": feature_ids, "assay_ids": assay_ids,
                  "targets": {}}
    population["hash"] = canonical_hash(
        {"feature_ids": feature_ids, "assay_ids": assay_ids})
    return population


# ---------- 集合の選択 ----------

def _included(rows: list[dict], *, role: str | None = None) -> list[int]:
    return [i for i, row in enumerate(rows)
            if row.get("include", True) and (role is None or row.get("role") == role)]


def _column(matrix, index: int, rows: list[int]):
    return np.asarray(matrix["values"], dtype=float)[rows, index]


# ---------- metric ごとの計算 ----------

def _pooled_qc_rsd(matrix, metadata, indices, entry, population):
    if np.asarray(matrix["imputed_mask"], dtype=bool).any():
        # 補完は分散を縮める。補完後の行列で RSD を合格にしない（spec §8）。
        return None, "imputed_matrix_not_allowed"

    groups: dict[tuple, list[int]] = {}
    for i in _included(metadata, role="qc"):
        pool = metadata[i].get("qc_pool")
        if not pool:
            # pool 不明はpooled QCの必須判定を満たさない（spec §7）。
            continue
        groups.setdefault((pool, metadata[i].get("batch")), []).append(i)
    if not groups:
        return None, "insufficient_qc_injections"
    # QC不足・QC不在のbatchも評価不能として母集団に残す。
    represented_batches = {batch for _, batch in groups}
    for i in _included(metadata, role="sample"):
        batch = metadata[i].get("batch")
        if batch not in represented_batches:
            groups.setdefault((None, batch), [])
    elements = []
    for index, feature_id in indices:
        for (pool, batch), rows in sorted(groups.items(), key=lambda kv: str(kv[0])):
            column = _column(matrix, index, rows)
            finite = column[np.isfinite(column)]
            if finite.size < _MIN_QC_FOR_RSD or float(np.mean(finite)) <= 0:
                elements.append(_element(f"{feature_id}@{pool}/{batch}", None,
                                         reason="insufficient_qc_values"))
                continue
            rsd = float(np.std(finite, ddof=1) / np.mean(finite) * 100.0)
            elements.append(_element(f"{feature_id}@{pool}/{batch}", rsd,
                                     threshold=entry["threshold"]))
    return elements, None


def _blank_fold(matrix, metadata, indices, entry, population):
    sample_rows = _included(metadata, role="sample")
    blank_rows = _included(metadata, role="blank")
    if not blank_rows or not sample_rows:
        return None, "blank_or_sample_missing"

    elements = []
    batches = sorted({metadata[i].get("batch") for i in sample_rows}, key=str)
    for index, feature_id in indices:
        for batch in batches:
            batch_samples = [i for i in sample_rows if metadata[i].get("batch") == batch]
            batch_blanks = [i for i in blank_rows if metadata[i].get("batch") == batch]
            element_id = f"{feature_id}@{batch}"
            samples = _column(matrix, index, batch_samples)
            blanks = _column(matrix, index, batch_blanks)
            sample_median = float(np.nanmedian(samples)) if samples.size else np.nan
            blank_median = float(np.nanmedian(blanks)) if blanks.size else np.nan
            if not np.isfinite(sample_median) or not np.isfinite(blank_median):
                elements.append(_element(element_id, None, reason="median_not_finite"))
            elif blank_median == 0 and sample_median > 0:
                elements.append(_element(element_id, None,
                                         special_value=_POSITIVE_INFINITY,
                                         threshold=entry["threshold"]))
            elif blank_median == 0:
                # 分子も 0。比は定義できない（0/0）。
                elements.append(_element(element_id, None, reason="blank_and_sample_zero"))
            else:
                elements.append(_element(element_id, sample_median / blank_median,
                                         threshold=entry["threshold"]))
    return elements, None


def _detection_rate(matrix, metadata, indices, entry, population):
    detected = matrix.get("detected_mask")
    if detected is None:
        return None, "detection_unknown"
    detected = np.asarray(detected, dtype=bool)
    rows = _included(metadata, role="sample")
    if not rows:
        return None, "no_included_samples"

    elements = []
    for index, feature_id in indices:
        column = detected[rows, index]
        # 分母は included の全注入。検出0の注入を分母から外すと率が常に1になる。
        elements.append(_element(feature_id, float(column.sum()) / len(rows),
                                 threshold=entry["threshold"]))
    return elements, None


def _standard_metric(matrix, metadata, indices, entry, population, evidence,
                     *, kind: str):
    if not (evidence or {}).get("availability"):
        return None, "assay_evidence_unavailable"
    targets = (population or {}).get("targets") or {}
    if not targets:
        return None, "target_expectations_missing"

    rows = {(r.get("feature_id"), r.get("assay_id")): r
            for r in (evidence.get("rows") or [])}
    allowed_assays = set((population or {}).get("assay_ids")
                         or matrix["assay_ids"])

    elements = []
    for target_id, spec in sorted(targets.items()):
        feature_id = spec.get("feature_id")
        expected_rt = spec.get("expected_rt_min")
        expected_mz = spec.get("expected_mz")
        for assay_id in matrix["assay_ids"]:
            if assay_id not in allowed_assays:
                continue
            row = rows.get((feature_id, assay_id))
            element_id = f"{target_id}@{assay_id}"
            if row is None or row.get("detection_status") != "detected":
                elements.append(_element(element_id, None,
                                         reason="standard_injection_not_detected"))
                continue
            if kind == "rt":
                observed = row.get("observed_rt_min")
                value = (None if observed is None or expected_rt is None
                         else abs(float(observed) - float(expected_rt)))
                elements.append(_element(element_id, value,
                                         threshold=entry["threshold"]))
            else:
                observed = row.get("observed_mz")
                value = (None if observed is None or not expected_mz
                         else (float(observed) - float(expected_mz))
                         / float(expected_mz) * 1e6)
                # 値の符号は残し、判定は絶対値で行う（spec §9 の表）。
                elements.append(_element(element_id, value,
                                         threshold=entry["threshold"],
                                         absolute=True))
    if not elements:
        return None, "no_standard_injections"
    return elements, None


def _internal_standard_valid_fraction(matrix, metadata, indices, entry, population):
    locked = np.asarray(matrix["locked_mask"], dtype=bool)
    rows = _included(metadata, role="sample")
    if not rows:
        return None, "no_included_samples"

    elements = []
    for index, feature_id in indices:
        column = locked[rows, index]
        elements.append(_element(feature_id,
                                 float((~column).sum()) / len(rows),
                                 threshold=entry["threshold"]))
    return elements, None


# ---------- 公開関数 ----------

def evaluate_qc(matrix: dict, evidence: dict, metadata: list[dict],
                policy: list[dict], population: dict | None) -> dict:
    """QC を評価し、metric結果・バッチ状態・eligibility提案を返す。

    `metadata` は解析行列の assay 列と**同じ順・同じ件数**の
    sample-manifest 行（`apply_metadata` の契約と同じ位置対応）。
    `population` が `None` なら、この呼び出しで母集団を固定して返す。以後の
    （処理後）QC は返ってきたものをそのまま渡して同じ集合で評価する。

    この関数は `matrix` を書き換えない。QCで落ちた feature は
    `eligibility_suggestion` として返すだけで、適用は呼び出し側の仕事。
    """
    assay_ids = list(matrix["assay_ids"])
    if len(metadata) != len(assay_ids):
        raise DomainError(
            _POLICY_INVALID,
            f"metadataの件数が解析行列のassay数と一致しません"
            f"（metadata={len(metadata)}, assays={len(assay_ids)}）。",
            {"metadata": len(metadata), "assays": len(assay_ids)})

    resolved_population = dict(population or {})
    metrics: dict[str, dict] = {}
    warnings: list[str] = []
    exclude: set[str] = set()

    for entry in policy:
        metric = entry.get("metric")
        if metric not in _METRICS:
            raise DomainError(
                _POLICY_INVALID,
                f"未知のQC metricです: {metric!r}",
                {"metric": metric, "supported": sorted(_METRICS)})

        if metric not in resolved_population:
            resolved_population[metric] = _fix_population(matrix, entry)
        metric_population = resolved_population[metric]

        index_of = {fid: i for i, fid in enumerate(matrix["feature_ids"])}
        indices = [(index_of[fid], fid)
                   for fid in metric_population["feature_ids"] if fid in index_of]

        if metric == "pooled_qc_rsd":
            elements, reason = _pooled_qc_rsd(matrix, metadata, indices, entry,
                                              metric_population)
        elif metric == "blank_fold":
            elements, reason = _blank_fold(matrix, metadata, indices, entry,
                                           metric_population)
        elif metric == "detection_rate":
            elements, reason = _detection_rate(matrix, metadata, indices, entry,
                                               metric_population)
        elif metric == "standard_rt_error":
            elements, reason = _standard_metric(matrix, metadata, indices, entry,
                                                metric_population, evidence,
                                                kind="rt")
        elif metric == "standard_mass_error_ppm":
            elements, reason = _standard_metric(matrix, metadata, indices, entry,
                                                metric_population, evidence,
                                                kind="mz")
        else:
            elements, reason = _internal_standard_valid_fraction(
                matrix, metadata, indices, entry, metric_population)

        if elements is None:
            metrics[metric] = {
                "status": "not_evaluable", "reason": reason, "elements": [],
                "counts": {"pass": 0, "fail": 0, "not_evaluable": 0, "total": 0},
                "required": bool(entry.get("required")),
                "population_hash": metric_population.get("hash"),
                "n_assays_used": 0,
            }
            continue

        counts = {
            "pass": sum(1 for e in elements if e["status"] == "pass"),
            "fail": sum(1 for e in elements if e["status"] == "fail"),
            "not_evaluable": sum(1 for e in elements
                                 if e["status"] == "not_evaluable"),
        }
        counts["total"] = sum(counts.values())
        fraction = entry.get("minimum_pass_fraction")
        # 集計指定が無い項目は「全要素が pass で初めて pass」として扱う
        # （t=1.0）。要素が1件なら spec の「その判定を直接使う」と同じになる。
        status = aggregate_counts(counts["pass"], counts["fail"],
                                  counts["not_evaluable"],
                                  1.0 if fraction is None else float(fraction))

        metrics[metric] = {
            "status": status,
            "reason": None,
            "elements": elements,
            "counts": counts,
            "required": bool(entry.get("required")),
            "minimum_pass_fraction": fraction,
            "population_hash": metric_population.get("hash"),
            "n_assays_used": len(metric_population.get("assay_ids") or []),
        }
        if status == "fail" and not entry.get("required"):
            warnings.append(
                f"任意QC項目 {metric} が fail です（バッチ合否には算入していません）。")
        for element in elements:
            if element["status"] == "fail":
                exclude.add(str(element["element_id"]).split("@", 1)[0])

    required = [m for m in metrics.values() if m["required"]]
    if not required:
        batch_status = "not_evaluable"
    elif any(m["status"] == "fail" for m in required):
        batch_status = "fail"
    elif any(m["status"] == "not_evaluable" for m in required):
        batch_status = "not_evaluable"
    else:
        batch_status = "pass"

    known_features = set(matrix["feature_ids"])
    return {
        "schema": QC_SCHEMA,
        "matrix_id": matrix.get("matrix_id"),
        "stage": matrix.get("stage"),
        "metrics": metrics,
        "batch_status": batch_status,
        "warnings": warnings,
        "population": resolved_population,
        # 提案であって適用ではない。適用しても QC は再集計しない
        # （除外で fail を分母から消すと、落とすほど合格に近づく）。
        "eligibility_suggestion": {
            "exclude_feature_ids": sorted(exclude & known_features),
            "applied": False,
        },
    }
