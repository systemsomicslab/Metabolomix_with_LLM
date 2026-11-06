"""2 つのアダプタ（解析フォルダ / mzTab-M）が同じアラインメントから同じ DatasetState を作るかの照合。

実データでの突き合わせ（B-1 Task 11）と、合成データのテストの両方で使う。
"""
from __future__ import annotations

import math

import numpy as np

from metabolomix.analysis.feature_search import _window_identity

__all__ = ["compare_states", "KNOWN_DIFFERENCES"]

_MAX_EXAMPLES = 5

KNOWN_DIFFERENCES = [
    "解析フォルダ側の SME 相当は、上流 ShouldWriteSmeLine の内部標準による正規化の条件を移植していない"
    "（パラメータを読まないため）。内部標準で正規化した解析の mzTab-M とは SME の有無が食い違いうる。",
]


def _identity(ds, fid: str) -> dict:
    # mzTab 側は全特徴に [] を持つ。キー欠落と [] は同じ「SME なし」として扱う。
    return _window_identity(ds.feature_candidates.get(fid) or [], ds.feature_annotations.get(fid))


def compare_states(folder_ds, mztab_ds, *, rel_tol: float = 1e-4) -> dict:
    ids_f, ids_m = list(folder_ds.feature_ids), list(mztab_ds.feature_ids)
    set_f, set_m = set(ids_f), set(ids_m)
    only_folder = [f for f in ids_f if f not in set_m]
    only_mztab = [f for f in ids_m if f not in set_f]
    common = [f for f in ids_f if f in set_m]

    names_f, names_m = list(folder_ds.sample_names), list(mztab_ds.sample_names)
    same_set = set(names_f) == set(names_m)
    same_order = names_f == names_m
    common_samples = [n for n in names_f if n in set(names_m)]

    f_row = {fid: i for i, fid in enumerate(ids_f)}
    m_row = {fid: i for i, fid in enumerate(ids_m)}
    f_col = {n: j for j, n in enumerate(names_f)}
    m_col = {n: j for j, n in enumerate(names_m)}

    n_cmp = n_over = n_det_cmp = n_det_diff = 0
    worst = 0.0
    f_mask, m_mask = folder_ds.detected_mask, mztab_ds.detected_mask
    compare_detected = f_mask is not None and m_mask is not None
    for fid in common:
        for name in common_samples:
            a = folder_ds.feature_matrix[f_row[fid], f_col[name]]
            b = mztab_ds.feature_matrix[m_row[fid], m_col[name]]
            a_nan, b_nan = bool(np.isnan(a)), bool(np.isnan(b))
            n_cmp += 1
            if a_nan != b_nan:
                n_over += 1
                worst = math.inf
            elif not a_nan:
                rel = abs(a - b) / max(abs(a), abs(b), 1.0)
                worst = max(worst, rel)
                if rel > rel_tol:
                    n_over += 1
            if compare_detected:
                n_det_cmp += 1
                if bool(f_mask[f_row[fid], f_col[name]]) != bool(m_mask[m_row[fid], m_col[name]]):
                    n_det_diff += 1

    name_examples, evidence_examples = [], []
    n_name_diff = n_evidence_diff = 0
    for fid in common:
        a, b = _identity(folder_ds, fid), _identity(mztab_ds, fid)
        if a["evidence"] != b["evidence"]:
            n_evidence_diff += 1
            if len(evidence_examples) < _MAX_EXAMPLES:
                evidence_examples.append(
                    {"feature_id": fid, "folder": a["evidence"], "mztab": b["evidence"]})
        if a["name"] != b["name"]:
            n_name_diff += 1
            if len(name_examples) < _MAX_EXAMPLES:
                name_examples.append({"feature_id": fid, "folder": a["name"], "mztab": b["name"]})

    worst_rel = worst if math.isinf(worst) else float(f"{worst:.6g}")
    checks = {
        "feature_ids": {"only_folder": only_folder[:_MAX_EXAMPLES],
                        "only_mztab": only_mztab[:_MAX_EXAMPLES], "n_common": len(common)},
        "samples": {"same_set": same_set, "same_order": same_order},
        "values": {"n_compared": n_cmp, "n_over": n_over, "worst_rel": worst_rel},
        "detected": {"n_compared": n_det_cmp, "n_diff": n_det_diff},
        "names": {"n_diff": n_name_diff, "examples": name_examples},
        "evidence": {"n_diff": n_evidence_diff, "examples": evidence_examples},
    }
    ok = (not only_folder and not only_mztab and same_set and n_over == 0 and n_det_diff == 0
          and n_name_diff == 0 and n_evidence_diff == 0)
    return {"ok": ok, "checks": checks, "known_differences": list(KNOWN_DIFFERENCES)}
