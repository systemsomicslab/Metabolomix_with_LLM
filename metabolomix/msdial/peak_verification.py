"""ピークアノテーションの生化学的検証（純ロジック層、MCP 非依存）。

分子式→質量、アダクト→m/z、精密質量誤差、アダクト-極性整合、脂質クラストークン
抽出を提供する。``server.py`` の ``verify_peak_annotation`` ツールがこれを使う。
``metabolomix.corpus.knowledge_store`` と同じく MCP に依存しない純関数群で、単体テスト可能。
"""

from __future__ import annotations

import re

from metabolomix.core import session_state
from metabolomix.library.defaults import DEFAULT_MS2_TOL as _LIBRARY_MS2_TOL
from metabolomix.library.defaults import DEFAULT_MZ_TOL as _LIBRARY_MZ_TOL
from metabolomix.library.defaults import DEFAULT_RT_TOL as _LIBRARY_RT_TOL
from metabolomix.library.defaults import pick_tol as _pick_tol

# session_state は leaf 側（arf/arf2/eic/pai2 の reader・msdial.classes/tags・
# analysis.preprocessing のみに依存）で、そのどれも peak_verification を import
# しないため循環しない。実際に `python -c "from metabolomix.core import
# session_state; ...sys.modules..."` で peak_verification が読み込まれない
# ことを確認済み（Task 11 レポート参照）。
#
# `metabolomix.library.defaults` は stdlib のみの leaf（Task 11 レビュー
# Important 1 で新設）。`metabolomix.library.tools`（MCP 面）と同じ既定許容幅を
# ここから読む——`library.tools` 自体は import しない（`@mcp.tool` の登録副作用と
# mcp_core / session_state 経由の上位層を引き込むため。peak_verification は
# 下位レイヤに留める）。

# --- 定数（モノアイソトピック質量） ---
PROTON_MASS = 1.00727646
ELECTRON_MASS = 0.00054858

ELEMENT_MASSES: dict[str, float] = {
    "H": 1.0078250319,
    "C": 12.0,
    "N": 14.0030740052,
    "O": 15.9949146221,
    "P": 30.97376151,
    "S": 31.97207069,
    "Na": 22.98976928,
    "Cl": 34.96885271,
    "K": 38.9637069,
    "D": 2.0141017779,
    "F": 18.9984031627,
    "Br": 78.9183376,
    "C13": 13.0033548378,
}

_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> dict[str, int]:
    """``"C42H82NO8P"`` のような分子式を 元素→個数 に分解する。

    数の省略は 1 とみなす。空/None は ``ValueError``。
    """
    if not formula or not formula.strip():
        raise ValueError("empty formula")
    counts: dict[str, int] = {}
    for element, digits in _FORMULA_TOKEN_RE.findall(formula.strip()):
        if not element:
            continue
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    if not counts:
        raise ValueError(f"unparseable formula: {formula!r}")
    return counts


def monoisotopic_mass(counts: dict[str, int]) -> float:
    """元素カウントから中性モノアイソトピック質量を返す。未知元素は ``KeyError``。"""
    return sum(ELEMENT_MASSES[element] * n for element, n in counts.items())


# --- アダクト表（(sign, shift, charge, n_mol); m/z = (n_mol*neutral + shift) / charge） ---
# sign はインデックス0のまま（adduct_consistency が entry[0] を参照）。多量体・多価に対応。
ADDUCT_SHIFTS: dict[str, tuple[str, float, int, int]] = {
    "[M+H]+": ("+", PROTON_MASS, 1, 1),
    "[M+NH4]+": ("+", 18.03382555, 1, 1),
    "[M+Na]+": ("+", 22.9892207, 1, 1),
    "[M-H2O+H]+": ("+", -17.00328823, 1, 1),
    "[M-H]-": ("-", -1.00727646, 1, 1),
    "[M+HCOO]-": ("-", 44.99820286, 1, 1),
    "[M+FA-H]-": ("-", 44.99820286, 1, 1),  # [M+HCOO]- の別名
    "[M+CH3COO]-": ("-", 59.01385292, 1, 1),
    "[M+Cl]-": ("-", 34.96940129, 1, 1),
    "[2M-H]-": ("-", -1.00727646, 1, 2),
    "[M-2H]2-": ("-", -2 * PROTON_MASS, 2, 1),
}


def adduct_mz(neutral_mass: float, adduct: str) -> float | None:
    """中性質量とアダクトから観測 m/z を返す。未知アダクトは ``None``。

    m/z = (n_mol * neutral + shift) / charge。多量体（[2M-H]-）・多価（[M-2H]2-）に対応。
    """
    entry = ADDUCT_SHIFTS.get(adduct)
    if entry is None:
        return None
    sign, shift, charge, n_mol = entry
    return (n_mol * neutral_mass + shift) / charge


def _band_for_ppm(ppm: float, pass_ppm: float, borderline_ppm: float) -> str:
    magnitude = abs(ppm)
    if magnitude <= pass_ppm:
        return "PASS"
    if magnitude <= borderline_ppm:
        return "BORDERLINE"
    return "FAIL"


def mass_error_ppm(
    observed_mz,
    formula,
    adduct,
    *,
    pass_ppm: float = 5.0,
    borderline_ppm: float = 10.0,
) -> dict:
    """実測 m/z と 分子式+アダクト の理論 m/z から ppm 誤差と帯を返す。

    分子式/アダクトが欠落・不明・解釈不能なら ``band="UNKNOWN"`` を返し例外は送出しない。
    """
    unknown = {"theoretical_mz": None, "ppm": None, "band": "UNKNOWN"}
    if observed_mz is None or not formula or formula == "Unknown":
        return unknown
    if not adduct or adduct == "Unknown":
        return unknown
    try:
        neutral = monoisotopic_mass(parse_formula(formula))
    except (ValueError, KeyError):
        return unknown
    theoretical = adduct_mz(neutral, adduct)
    if theoretical is None or theoretical == 0:
        return unknown
    ppm = (float(observed_mz) - theoretical) / theoretical * 1e6
    return {
        "theoretical_mz": round(theoretical, 4),
        "ppm": round(ppm, 2),
        "band": _band_for_ppm(ppm, pass_ppm, borderline_ppm),
    }


# --- クラス別の典型アダクト（助言のみ。判定はしない） ---
CLASS_TYPICAL_ADDUCTS: dict[str, list[str]] = {
    "pc": ["[M+H]+", "[M+HCOO]-", "[M+CH3COO]-"],
    "lpc": ["[M+H]+", "[M+HCOO]-"],
    "sm": ["[M+H]+", "[M+HCOO]-"],
    "pe": ["[M+H]+", "[M-H]-"],
    "pg": ["[M-H]-", "[M+NH4]+"],
    "pi": ["[M-H]-", "[M+NH4]+"],
    "ps": ["[M-H]-", "[M+H]+"],
    "tg": ["[M+NH4]+", "[M+Na]+"],
    "dg": ["[M+NH4]+", "[M+Na]+", "[M+H]+"],
    "cer": ["[M+H]+", "[M-H]-", "[M+HCOO]-"],
    "che": ["[M+NH4]+", "[M+Na]+"],
    "ce": ["[M+NH4]+", "[M+Na]+"],
    "fa": ["[M-H]-"],
}


def _class_key(ontology: str | None) -> str | None:
    if not ontology:
        return None
    token = ontology.strip().lower()
    if token in CLASS_TYPICAL_ADDUCTS:
        return token
    head = token.split()[0] if token.split() else token
    return head if head in CLASS_TYPICAL_ADDUCTS else None


def adduct_consistency(adduct, ion_mode, ontology=None) -> dict:
    """アダクトの電荷符号と実測極性の整合を判定し、クラス典型性を助言する。

    極性一致は決定的に PASS/FAIL。クラス典型性は助言のみで band には影響しない。
    """
    entry = ADDUCT_SHIFTS.get(adduct) if adduct else None
    if entry is None:
        return {
            "polarity_ok": None,
            "band": "UNKNOWN",
            "class_typical": None,
            "advisory": "アダクトが不明または未対応のため整合判定不可。",
        }

    sign = entry[0]
    mode = str(ion_mode).strip().lower()
    is_positive = mode.startswith("pos")
    is_negative = mode.startswith("neg")
    mode_known = is_positive or is_negative

    key = _class_key(ontology)
    if key is None:
        class_typical = None
        advisory = "クラス未知のため典型アダクト助言なし。"
    else:
        typical = CLASS_TYPICAL_ADDUCTS[key]
        class_typical = adduct in typical
        if class_typical:
            advisory = f"{key.upper()} で {adduct} は典型的。"
        else:
            advisory = f"{key.upper()} の典型は {', '.join(typical)}（{adduct} は非典型だが誤りとは限らない）。"

    if not mode_known:
        advisory = advisory + f" （ion_mode={ion_mode!r} が正/負極性として認識できないため極性整合を判定できませんでした。）"
        return {
            "polarity_ok": None,
            "band": "UNKNOWN",
            "class_typical": class_typical,
            "advisory": advisory,
        }

    polarity_ok = (sign == "+" and is_positive) or (sign == "-" and is_negative)
    band = "PASS" if polarity_ok else "FAIL"

    return {
        "polarity_ok": polarity_ok,
        "band": band,
        "class_typical": class_typical,
        "advisory": advisory,
    }


_ETHER_RE = re.compile(r"\b[POpo]-")


def extract_class_token(name: str | None, ontology: str | None) -> str | None:
    """脂質クラストークン（小文字）を返す。ontology 優先、無ければ name 先頭語。"""
    for source in (ontology, name):
        if source and source.strip():
            head = source.strip().lower().split()[0]
            if head:
                return head
    return None


def ether_caveats(name: str | None, ontology: str | None) -> list[str]:
    """エーテル脂質（P-/O- 表記）に該当する場合、注意ノートへのリンクを返す。"""
    text = f"{name or ''} {ontology or ''}"
    if _ETHER_RE.search(text):
        return [
            "エーテル脂質の P-/O- 表記混同に注意（[[pe-p-vs-pe-o-annotation]]）。"
            "P- は酸化ストレス仮説 [[plasmalogen-oxidation]] に関与するが O- は別物。"
        ]
    return []


def _spectral_match_for_feature(feature: dict, spectrum: list, store) -> dict:
    """band=PASS の実スペクトルを、読み込み済みの参照ライブラリと照合する。

    `library_match_feature` MCP ツール（`metabolomix/library/tools.py`）と同じ土俵
    （`store.candidates()` によるアラインメント検索 + `spectral_match.match_spectrum`
    採点 + `spectral_match.total_score` による順位付け）を使うが、ここは検証ドシエの
    一材料でしかないので、戻り値は最良候補のスコアだけに絞る。候補一覧が要る場合は
    `library_match_feature` を使うこと。**「最良」の基準は両者で必ず一致させること**
    ——食い違うと、同じ feature について 2 つのツールが別の候補を名指しする。

    どの段階で失敗しても（precursor m/z が無い・ライブラリ照会が失敗する等）
    例外は投げない——照合は msms_evidence の主目的である band 判定を止めてはいけない。
    """
    from metabolomix.analysis.spectral_match import match_spectrum, total_score

    precursor_mz = feature.get("m/z")
    if precursor_mz is None:
        return {"status": "unavailable", "reason": "precursor m/z が無いため照合できません。"}

    ion_mode = feature.get("ion_mode")
    # store.candidates() の ion_mode 比較は COLLATE NOCASE（Critical 1 修正後）なので
    # ここで小文字へ揃える必要は無い——PAI2 由来の Enum の `.name`（"Positive"）を
    # そのまま渡してよい。
    ion_mode_name = ion_mode.name if hasattr(ion_mode, "name") else (
        str(ion_mode) if ion_mode is not None else None)
    rt = (feature.get("time") or {}).get("rt")

    try:
        search_params = store.summary().get("search_params") or {}
    except Exception:  # noqa: BLE001 - ここも下の2箇所と同じ意図: 照合の失敗で
        search_params = {}  # msms_evidence 本来の band 判定を止めない（既定値へ）
    # `or` ではなく共有の `pick_tol`（`is None` 判定）を使う——`search_params` の
    # `0.0`（"足切りなし"のような正当な値）を偽値扱いで既定値に差し替えないため
    # （最終レビュー Important 5: ここが独自の `or` 判定を持っていて、
    # `metabolomix.library.tools._pick_tol` とロジックがずれていた）。
    mz_tol = _pick_tol(None, search_params, "ms1_tolerance", _LIBRARY_MZ_TOL)
    ms2_tol = _pick_tol(None, search_params, "ms2_tolerance", _LIBRARY_MS2_TOL)
    rt_tol = _pick_tol(None, search_params, "rt_tolerance", _LIBRARY_RT_TOL)
    # 採点前処理（`normalize_measured`）のパラメータも search_params から採る
    # （最終レビュー Important 2: 検証 CLI と同じ集合を使わないと、足切りを変えた
    # run の `.dbs` を読ませたときだけ `library_match_feature` と違う土俵で
    # 採点することになる）。
    mass_begin = _pick_tol(None, search_params, "mass_range_begin", 0.0)
    mass_end = _pick_tol(None, search_params, "mass_range_end", 2000.0)
    relative_amp_cutoff = _pick_tol(None, search_params, "relative_amp_cutoff", 0.0)
    absolute_amp_cutoff = _pick_tol(None, search_params, "absolute_amp_cutoff", 0.0)

    try:
        candidates = store.candidates(
            precursor_mz, mz_tol=mz_tol, ion_mode=ion_mode_name, rt=rt, rt_tol=rt_tol)
    except Exception as exc:  # noqa: BLE001 - ライブラリ照合の失敗は band 判定を止めない
        return {"status": "error", "reason": str(exc)}

    if not candidates:
        return {"status": "no_candidates", "n_candidates": 0}

    # 並び順は `library_match_feature` と**同じ**総合スコア（`GetTotalScore`）。
    # 揃えないと、同じ feature について検証ドシエと候補一覧が違う「最良候補」を
    # 名指しすることになる。RT 項の可否も同じく store の
    # `IsUseTimeForAnnotationScoring`（`.msp` 由来なら上流既定の False）に従う。
    use_rt_scoring = bool(search_params.get("use_time_for_annotation_scoring", False))

    scored = []
    for record in candidates:
        result = match_spectrum(
            spectrum, record["spectrum"], ms2_tol=ms2_tol,
            mass_begin=mass_begin, mass_end=mass_end,
            relative_amp_cutoff=relative_amp_cutoff, absolute_amp_cutoff=absolute_amp_cutoff,
        )
        # store のレコードは precursor m/z / RT を欠きうる（`.msp` 由来など）。
        # `total_score` は欠測を番兵で落とすので、ここで弾く必要はない。
        result.update(total_score(
            result,
            precursor_mz=precursor_mz, reference_precursor_mz=record.get("precursor_mz"),
            ms1_tol=mz_tol,
            rt=rt, reference_rt=record.get("rt"),
            rt_tol=rt_tol, use_rt=use_rt_scoring,
        ))
        scored.append((record, result))
    scored.sort(key=lambda item: item[1]["total_score"], reverse=True)
    best_record, best_result = scored[0]

    return {
        "status": "matched",
        "n_candidates": len(candidates),
        "best_match": {
            "name": best_record.get("name"),
            "ontology": best_record.get("ontology"),
            "adduct": best_record.get("adduct"),
            "total_score": round(best_result["total_score"], 6),
            "weighted_dot_product": round(best_result["weighted_dot_product"], 6),
            "simple_dot_product": round(best_result["simple_dot_product"], 6),
            "reverse_dot_product": round(best_result["reverse_dot_product"], 6),
            "matched_peaks_percentage": round(best_result["matched_peaks_percentage"], 6),
            "matched_peaks_count": best_result["matched_peaks_count"],
        },
    }


def msms_evidence(feature: dict, top_n: int = 5) -> dict:
    """MS/MS 証拠の band と主要フラグメントを返す（決定的判定）。

    MSI Level 2 は「MS/MS がある」ことを根拠にするが、PAI2 の ``has_msms`` は
    **取得参照があること**を示すだけで実スペクトルの有無を保証しない
    （output-format core §9-8）。実体を持つのは同名 ``.dcl``（MSDecResult）であり、
    ``dcl_reader.attach_msms_to_features`` が索引対応で ``msms_spectrum`` を書き込む。
    ここではその2状態を明示的に分ける:

    - ``PASS``      : 実スペクトルあり（source=``spectrum``）。確度主張の根拠に使える。
    - ``FLAG_ONLY`` : 取得フラグのみでスペクトル本体が無い（source=``flag``）。.dcl 未添付か
                      本当に空か区別できないため caveat を付す。
    - ``ABSENT``    : フラグもスペクトルも無い。

    この3状態の契約は変えない。``spectral_match`` はその内側に足すだけの追加情報:
    ``band == "PASS"`` かつ参照ライブラリが読み込み済み（``session.library.store``
    がある）ときだけ、実スペクトルをライブラリと照合した最良候補のスコアを載せる。
    ライブラリ未読み込みなら ``spectral_match`` キー自体を持たない
    （``dict.get("spectral_match")`` は ``None``）。

    ``top_fragments`` は強度降順の上位 ``top_n`` 本（``[[mz, intensity], ...]``）。
    ``n_peaks`` は間引き前の元本数を保つ。
    """
    spectrum = feature.get("msms_spectrum") or []
    declared = feature.get("n_msms_peaks")
    n_peaks = int(declared) if isinstance(declared, (int, float)) else len(spectrum)

    if spectrum:
        top = sorted(spectrum, key=lambda p: p[1], reverse=True)[:top_n]
        result = {
            "band": "PASS",
            "source": "spectrum",
            "n_peaks": n_peaks,
            "top_fragments": [[p[0], p[1]] for p in top],
            "caveat": None,
        }
        store = getattr(getattr(session_state.session, "library", None), "store", None)
        if store is not None:
            result["spectral_match"] = _spectral_match_for_feature(feature, spectrum, store)
        return result
    if feature.get("has_msms"):
        return {
            "band": "FLAG_ONLY",
            "source": "flag",
            "n_peaks": n_peaks,
            "top_fragments": [],
            "caveat": (
                "MS/MS 取得フラグ（has_msms）は立っているが、実スペクトルが添付されていない。"
                "PAI2 の has_msms は取得参照の有無を示すだけなので、これだけを根拠に同定確度を"
                "主張しないこと。実体は同名 .dcl にあり、dcl_parser で読み込むと付与される。"
            ),
        }
    return {"band": "ABSENT", "source": None, "n_peaks": n_peaks or 0,
            "top_fragments": [], "caveat": None}


def vocab_hits(class_token: str | None, vocab: dict[str, list[str]]) -> list[str]:
    """クラストークンに対応する語彙同義語を返す（無ければ空）。"""
    if not class_token:
        return []
    return list(vocab.get(class_token, []))
