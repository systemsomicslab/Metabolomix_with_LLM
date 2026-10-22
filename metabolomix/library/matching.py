"""参照ライブラリの照合の芯（`library_match_feature` と、再解析の MS/MS の合否・証拠の束が共有する）。

許容幅の解決（明示 → store の `search_params` → 既定）、候補の検索、採点、総合スコアでの
並べ替えを 1 か所に置く。セッションにも MCP にも依存しない（store と測定スペクトルを引数で
受ける）——並行のハンドルから呼んでも、直前の呼び出しの状態に左右されない。
"""
from __future__ import annotations

from metabolomix.analysis.spectral_match import cutoff_mask, match_spectrum, total_score
from metabolomix.library.defaults import DEFAULT_MS2_TOL, DEFAULT_MZ_TOL, DEFAULT_RT_TOL, pick_tol

__all__ = ["RANK_KEY", "score_candidates"]

#: 並べ替えの基準（MS-DIAL の総合スコア。library_match_feature の docstring を参照）。
RANK_KEY = "total_score"


def score_candidates(store, *, precursor_mz: float, rt: float | None, ion_mode: str | None,
                     measured: list, mz_tol: float | None = None, ms2_tol: float | None = None,
                     rt_tol: float | None = None, match_fn=None) -> dict:
    """候補を引いて採点し、総合スコアの降順（同点は候補順）に並べる。

    戻り値: `status`（`success` / `no_candidates`）、`query`（解決した許容幅と採点の前処理）、
    `scoring`（総合スコアの組み立て）、`scored`（`(record, result)` の list。全候補）、
    `unscored_mz`（足切りで採点から外れた測定ピーク）、`n_candidates`。

    `match_fn` は採点の関数（既定 `spectral_match.match_spectrum`）。`library_match_feature` は
    自分の module の名前を渡す（既存のテストがそこを差し替えて順位付けを確かめている）。
    """
    match_fn = match_fn or match_spectrum
    search_params = store.summary().get("search_params") or {}
    resolved_mz_tol = pick_tol(mz_tol, search_params, "ms1_tolerance", DEFAULT_MZ_TOL)
    resolved_ms2_tol = pick_tol(ms2_tol, search_params, "ms2_tolerance", DEFAULT_MS2_TOL)
    resolved_rt_tol = pick_tol(rt_tol, search_params, "rt_tolerance", DEFAULT_RT_TOL)
    mass_begin = pick_tol(None, search_params, "mass_range_begin", 0.0)
    mass_end = pick_tol(None, search_params, "mass_range_end", 2000.0)
    relative_amp_cutoff = pick_tol(None, search_params, "relative_amp_cutoff", 0.0)
    absolute_amp_cutoff = pick_tol(None, search_params, "absolute_amp_cutoff", 0.0)
    use_rt_scoring = bool(search_params.get("use_time_for_annotation_scoring", False))
    query = {
        "precursor_mz": precursor_mz, "rt": rt, "ion_mode": ion_mode,
        "mz_tol": resolved_mz_tol, "ms2_tol": resolved_ms2_tol, "rt_tol": resolved_rt_tol,
        "mass_begin": mass_begin, "mass_end": mass_end,
        "relative_amp_cutoff": relative_amp_cutoff, "absolute_amp_cutoff": absolute_amp_cutoff,
    }
    scoring = {
        "rule": "MS-DIAL GetTotalScore: rt + precursor_mz + (weighted+simple+reverse)/3 "
                "+ matched_peaks_percentage（正規化しない和。各項は > 0 のときだけ加算）",
        "use_rt": use_rt_scoring,
        "ms1_tol": resolved_mz_tol,
        "rt_tol": resolved_rt_tol if use_rt_scoring else None,
    }
    candidates = store.candidates(precursor_mz, mz_tol=resolved_mz_tol, ion_mode=ion_mode,
                                  rt=rt, rt_tol=resolved_rt_tol)
    out = {"query": query, "scoring": scoring, "scored": [], "unscored_mz": [],
           "n_candidates": len(candidates)}
    if not candidates:
        return {"status": "no_candidates", **out}

    kept = cutoff_mask(measured, relative_amp_cutoff=relative_amp_cutoff,
                       absolute_amp_cutoff=absolute_amp_cutoff)
    out["unscored_mz"] = [float(peak[0]) for peak, keep in zip(measured, kept) if not keep]
    scored = []
    for record in candidates:
        result = match_fn(
            measured, record["spectrum"], ms2_tol=resolved_ms2_tol,
            mass_begin=mass_begin, mass_end=mass_end,
            relative_amp_cutoff=relative_amp_cutoff, absolute_amp_cutoff=absolute_amp_cutoff,
        )
        result.update(total_score(
            result,
            precursor_mz=precursor_mz, reference_precursor_mz=record.get("precursor_mz"),
            ms1_tol=resolved_mz_tol,
            rt=rt, reference_rt=record.get("rt"),
            rt_tol=resolved_rt_tol, use_rt=use_rt_scoring,
        ))
        scored.append((record, result))
    scored.sort(key=lambda item: item[1][RANK_KEY], reverse=True)
    out["scored"] = scored
    return {"status": "success", **out}
