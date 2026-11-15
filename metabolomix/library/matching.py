"""参照ライブラリの照合の芯（`library_match_feature` と、再解析の MS/MS の合否・証拠の束が共有する）。

許容幅の解決（明示 → store の `search_params` → 既定）、候補の検索、採点、総合スコアでの
並べ替えを 1 か所に置く。セッションにも MCP にも依存しない（store と測定スペクトルを引数で
受ける）——並行のハンドルから呼んでも、直前の呼び出しの状態に左右されない。
"""
from __future__ import annotations

import re

from metabolomix.analysis.spectral_match import (
    cutoff_mask, fix_mass_tolerance, match_spectrum, total_score)
from metabolomix.library.defaults import DEFAULT_MS2_TOL, DEFAULT_MZ_TOL, DEFAULT_RT_TOL, pick_tol

__all__ = ["RANK_KEY", "dbs_cutoffs", "delta_ppm", "msdial_agreement", "normalize_ion_mode",
           "passes_dbs_cutoffs", "score_candidates"]

#: 並べ替えの基準（MS-DIAL の総合スコア。library_match_feature の docstring を参照）。
RANK_KEY = "total_score"


def score_candidates(store, *, precursor_mz: float, rt: float | None, ion_mode: str | None,
                     measured: list, mz_tol: float | None = None, ms2_tol: float | None = None,
                     rt_tol: float | None = None, mz_tol_ppm: float | None = None,
                     match_fn=None) -> dict:
    """候補を引いて採点し、総合スコアの降順（同点は候補順）に並べる。

    戻り値: `status`（`success` / `no_candidates`）、`query`（解決した許容幅と採点の前処理）、
    `scoring`（総合スコアの組み立て）、`scored`（`(record, result)` の list。全候補）、
    `unscored_mz`（足切りで採点から外れた測定ピーク）、`n_candidates`、`dbs_cutoffs`
    （`.dbs` の注釈の足切り。無ければ None。`passes_dbs_cutoffs` が使う）。

    前駆体の窓は MS-DIAL と同じ: `mz_tol`（m/z 500 以下での Da の幅）を、m/z 500 を超えたら同じ
    ppm に広げる（上流 `FixMassTolerance`。既定 0.01 Da なら `max(0.01 Da, m/z × 20 ppm)`）。
    実際に使った Da の幅は `query.mz_window`。
    `mz_tol_ppm`（任意）はその窓の**上に重ねる追加の絞り込み**で、窓を広げない:
    参照との差が ±`mz_tol_ppm` ppm（`(測定 - 参照) / 参照 * 1e6`）を超える候補を落とす。
    省略時は何も絞らない。

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
    mz_window = fix_mass_tolerance(resolved_mz_tol, precursor_mz)
    query = {
        "precursor_mz": precursor_mz, "rt": rt, "ion_mode": ion_mode,
        "mz_tol": resolved_mz_tol, "mz_window": mz_window, "mz_tol_ppm": mz_tol_ppm,
        "ms2_tol": resolved_ms2_tol, "rt_tol": resolved_rt_tol,
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
    candidates = store.candidates(precursor_mz, mz_tol=mz_window, ion_mode=ion_mode,
                                  rt=rt, rt_tol=resolved_rt_tol)
    if mz_tol_ppm is not None:
        # 参照の前駆体が 0 の候補は ppm を出せないので絞り込みの対象外（残す）。
        candidates = [r for r in candidates
                      if not r["precursor_mz"]
                      or abs(precursor_mz - r["precursor_mz"]) / r["precursor_mz"] * 1e6
                      <= mz_tol_ppm]
    out = {"query": query, "scoring": scoring, "scored": [], "unscored_mz": [],
           "n_candidates": len(candidates), "dbs_cutoffs": dbs_cutoffs(search_params)}
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


# ---------- 候補の行に添える材料 ----------

#: `.dbs` の注釈の足切り（`MsRefSearchParameterBase` の Key 9〜14）。dot の 3 つは二乗の値で持つ。
_CUTOFF_KEYS = ("squared_weighted_dot_cutoff", "squared_simple_dot_cutoff",
                "squared_reverse_dot_cutoff", "matched_peaks_percentage_cutoff",
                "total_score_cutoff", "minimum_spectrum_match")


def dbs_cutoffs(search_params: dict | None) -> dict | None:
    """store の `search_params` にある注釈の足切り。1 つも無ければ None（`.msp` 由来など）。"""
    found = {k: (search_params or {}).get(k) for k in _CUTOFF_KEYS}
    found = {k: v for k, v in found.items() if v is not None}
    return found or None


def passes_dbs_cutoffs(result: dict, cutoffs: dict | None) -> bool | None:
    """候補の採点（`match_spectrum` + `total_score`）が `.dbs` の足切りを全部通るか。足切りが無ければ None。

    MS-DIAL の `IsSpectrumMatch` の条件（weighted / simple / reverse dot・一致率・一致本数）と
    総合スコアの足切り。dot の足切りは二乗の値なので、こちらの dot（平方根）を二乗して比べる。
    比較不能の番兵 `-1` は通らない。**合否（msms-rules.v1）には当てない**材料の欄。
    """
    if not cutoffs:
        return None
    pairs = (("weighted_dot_product", "squared_weighted_dot_cutoff", True),
             ("simple_dot_product", "squared_simple_dot_cutoff", True),
             ("reverse_dot_product", "squared_reverse_dot_cutoff", True),
             ("matched_peaks_percentage", "matched_peaks_percentage_cutoff", False),
             ("matched_peaks_count", "minimum_spectrum_match", False),
             ("total_score", "total_score_cutoff", False))
    for score_key, cutoff_key, squared in pairs:
        cutoff = cutoffs.get(cutoff_key)
        if cutoff is None:
            continue
        value = result.get(score_key)
        if value is None or value < 0:
            return False
        if (value * value if squared else value) < cutoff:
            return False
    return True


def delta_ppm(query_mz: float, ref_mz) -> float | None:
    """`(測定 - 参照) / 参照 * 1e6`（測定が参照より高いと正）。参照が 0・無ければ None。"""
    if not ref_mz:
        return None
    return (float(query_mz) - float(ref_mz)) / float(ref_mz) * 1e6


def normalize_ion_mode(value, *, field: str = "ion_mode") -> str | None:
    """明示の極性を `positive` / `negative` にする（`pos` / `neg` / `+` / `-` も可）。None は None。

    読めない値は `ValueError`（黙って極性不明のレコードだけに絞らない。2026-11-13 の事実確認 A11）。
    """
    from metabolomix.core.polarity_names import normalize_polarity

    return None if value is None else normalize_polarity(value, field=field)


# ---------- MS-DIAL の注釈との食い違い ----------

_BLOCK_RE = re.compile(r"^[A-Z]{14}$")


def _key_block(value) -> str | None:
    block = value.strip()[:14].upper() if isinstance(value, str) else ""
    return block if _BLOCK_RE.match(block) else None


def _clean_name(name) -> str | None:
    """`|` の最後の区間から限定子を外し、大小を問わない比較の形にする（`Unknown` は None）。"""
    from metabolomix.msdial.lipid_identity import _clean_msdial_name

    if not isinstance(name, str) or not name.strip():
        return None
    clean, _ = _clean_msdial_name(name.split("|")[-1])
    clean = " ".join(clean.split()).casefold()
    return None if clean in ("", "unknown", "null") else clean


def _same_compound(record: dict, annotation: dict) -> bool:
    """ライブラリの記録が MS-DIAL の注釈と同じ物質か。InChIKey の先頭ブロック（両方にあれば）→
    脂質の種（両方が読めれば）→ 名前（限定子を外し大小を問わない）の順に比べる。"""
    from metabolomix.analysis.msms_rules import name_lipid_species, record_lipid_species

    a, b = _key_block(record.get("inchikey")), _key_block(annotation.get("inchikey"))
    if a and b:
        return a == b
    species_a = record_lipid_species(record)[0]
    species_b = name_lipid_species(annotation.get("name"))[0]
    if species_a and species_b:
        return species_a == species_b
    name_a, name_b = _clean_name(record.get("name")), _clean_name(annotation.get("name"))
    return bool(name_a and name_b and name_a == name_b)


def msdial_agreement(scored: list, annotation: dict | None) -> dict:
    """照合の上位と、そのスポットの MS-DIAL の注釈（`.arf2` の代表の名前・InChIKey）の食い違い。

    戻り値は `{msdial_name, msdial_rank, top_differs_from_msdial}`。`msdial_rank` は注釈と同じ
    物質の候補の最上位の順位（全候補の中。無ければ None）、`top_differs_from_msdial` は 1 位の
    候補が注釈と別の物質か。注釈が無い（未注釈・`Unknown`）・候補が無ければ空の辞書。
    """
    if not annotation or not scored or _clean_name(annotation.get("name")) is None:
        return {}
    rank = next((i + 1 for i, (record, _) in enumerate(scored)
                 if _same_compound(record, annotation)), None)
    return {"msdial_name": annotation["name"].split("|")[-1].strip(), "msdial_rank": rank,
            "top_differs_from_msdial": rank != 1}
