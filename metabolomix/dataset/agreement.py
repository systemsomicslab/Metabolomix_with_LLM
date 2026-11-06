"""2 つのアダプタ（解析フォルダ / mzTab-M）が同じアラインメントから同じ DatasetState を作るかの照合。

実データでの突き合わせ（B-1 Task 11）と、合成データのテストの両方で使う。
"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis.feature_search import _window_identity

__all__ = ["compare_states", "KNOWN_DIFFERENCES"]

_MAX_EXAMPLES = 5

KNOWN_DIFFERENCES = [
    "解析フォルダ側の SME 相当は、上流 ShouldWriteSmeLine の内部標準による正規化の条件を移植していない"
    "（パラメータを読まないため）。内部標準で正規化した解析の mzTab-M とは SME の有無が食い違いうる。",
    "Console の mzTab-M は Height を Math.Round(PeakHeightTop, 0)（C# 既定の偶数丸め）の整数で書く。"
    "解析フォルダのアダプタは .arf の丸める前の値を返すので、2 つの経路の値は最大 0.5 違う。"
    "比べるときは解析フォルダ側を偶数丸め（np.round）の整数にしてから見る（values.folder_rounded）。",
]


def _identity(ds, fid: str) -> dict:
    # mzTab 側は全特徴に [] を持つ。キー欠落と [] は同じ「SME なし」として扱う。
    return _window_identity(ds.feature_candidates.get(fid) or [], ds.feature_annotations.get(fid))


def _duplicates(names: list[str]) -> list[str]:
    seen, dups = set(), []
    for name in names:
        if name in seen and name not in dups:
            dups.append(name)
        seen.add(name)
    return dups


def _first_index(names: list[str]) -> dict[str, int]:
    index: dict[str, int] = {}
    for j, name in enumerate(names):
        index.setdefault(name, j)
    return index


def compare_states(folder_ds, mztab_ds, *, rel_tol: float = 1e-4) -> dict:
    ids_f, ids_m = list(folder_ds.feature_ids), list(mztab_ds.feature_ids)
    set_f, set_m = set(ids_f), set(ids_m)
    only_folder = [f for f in ids_f if f not in set_m]
    only_mztab = [f for f in ids_m if f not in set_f]
    common = [f for f in ids_f if f in set_m]

    names_f, names_m = list(folder_ds.sample_names), list(mztab_ds.sample_names)
    name_set_m = set(names_m)
    same_set = set(names_f) == name_set_m
    same_order = names_f == names_m
    same_length = len(names_f) == len(names_m)
    duplicates = _duplicates(names_f)
    duplicates += [n for n in _duplicates(names_m) if n not in duplicates]

    f_col, m_col = _first_index(names_f), _first_index(names_m)
    common_samples = [n for n in f_col if n in name_set_m]
    f_row, m_row = _first_index(ids_f), _first_index(ids_m)
    fr = np.array([f_row[f] for f in common], dtype=int)
    mr = np.array([m_row[f] for f in common], dtype=int)
    fc = np.array([f_col[n] for n in common_samples], dtype=int)
    mc = np.array([m_col[n] for n in common_samples], dtype=int)

    n_cmp = n_over = n_one_side = 0
    worst_rel = None
    if fr.size and fc.size:
        # 上流の書き出しは Height を偶数丸めの整数にする。アダプタの値は丸めないので、比べる前に
        # 解析フォルダ側だけを同じ丸め（np.round は偶数丸め）にそろえる。
        a = np.round(np.asarray(folder_ds.feature_matrix, dtype=float)[np.ix_(fr, fc)])
        b = np.asarray(mztab_ds.feature_matrix, dtype=float)[np.ix_(mr, mc)]
        n_cmp = int(a.size)
        a_nan, b_nan = np.isnan(a), np.isnan(b)
        one_side = a_nan != b_nan
        n_one_side = int(one_side.sum())
        finite = ~a_nan & ~b_nan
        if finite.any():
            fa, fb = a[finite], b[finite]
            rel = np.abs(fa - fb) / np.maximum(np.maximum(np.abs(fa), np.abs(fb)), 1.0)
            n_over = int((rel > rel_tol).sum())
            worst_rel = float(f"{float(rel.max()):.6g}")
        n_over += n_one_side

    f_mask, m_mask = folder_ds.detected_mask, mztab_ds.detected_mask
    skipped = f_mask is None or m_mask is None
    n_det_cmp = n_det_diff = 0
    if not skipped and fr.size and fc.size:
        da = np.asarray(f_mask, dtype=bool)[np.ix_(fr, fc)]
        db = np.asarray(m_mask, dtype=bool)[np.ix_(mr, mc)]
        n_det_cmp, n_det_diff = int(da.size), int((da != db).sum())

    name_examples, evidence_examples = [], []
    n_name_diff = n_evidence_diff = 0
    for fid in common:
        a_id, b_id = _identity(folder_ds, fid), _identity(mztab_ds, fid)
        if a_id["evidence"] != b_id["evidence"]:
            n_evidence_diff += 1
            if len(evidence_examples) < _MAX_EXAMPLES:
                evidence_examples.append(
                    {"feature_id": fid, "folder": a_id["evidence"], "mztab": b_id["evidence"]})
        if a_id["name"] != b_id["name"]:
            n_name_diff += 1
            if len(name_examples) < _MAX_EXAMPLES:
                name_examples.append(
                    {"feature_id": fid, "folder": a_id["name"], "mztab": b_id["name"]})

    checks = {
        "feature_ids": {"only_folder": only_folder[:_MAX_EXAMPLES],
                        "only_mztab": only_mztab[:_MAX_EXAMPLES], "n_common": len(common)},
        "samples": {"same_set": same_set, "same_order": same_order, "same_length": same_length,
                    "duplicates": duplicates[:_MAX_EXAMPLES]},
        "values": {"n_compared": n_cmp, "n_over": n_over, "n_nan_one_side": n_one_side,
                   "worst_rel": worst_rel, "folder_rounded": "half_to_even"},
        "detected": {"n_compared": n_det_cmp, "n_diff": n_det_diff, "skipped": skipped},
        "names": {"n_diff": n_name_diff, "examples": name_examples},
        "evidence": {"n_diff": n_evidence_diff, "examples": evidence_examples},
    }
    ok = bool(common and not only_folder and not only_mztab
              and same_set and same_length and not duplicates
              and n_over == 0 and n_det_diff == 0 and n_name_diff == 0 and n_evidence_diff == 0)
    return {"ok": ok, "checks": checks, "known_differences": list(KNOWN_DIFFERENCES)}
