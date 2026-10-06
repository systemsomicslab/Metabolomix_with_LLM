"""採点エンジン。上流 `MsScanMatching.cs` の定義の写し。

**ここで固定している値は上流の挙動そのもの**で、素直な実装に直すと落ちる。
瑕疵を含めて写しているのは mzTab の数値と比較可能にするため（spec §6.3）。
"""
import math

import pytest

from metabolomix.analysis import spectral_match as sm

# 参照窓のうち「正規化強度 > 0.1」が 5 個以上になるようにしてある。4 個以下だと
# weighted / reverse に peakCountPenalty が掛かり、同一スペクトルでも 1.0 にならない。
_A = [[100.0, 999.0], [200.0, 800.0], [300.0, 600.0], [400.0, 400.0], [500.0, 200.0]]


def test_an_identical_spectrum_scores_one():
    for fn in (sm.simple_dot_product, sm.weighted_dot_product, sm.reverse_dot_product):
        assert fn(_A, _A, bin_width=0.01) == pytest.approx(1.0, abs=1e-9)


def test_the_penalty_applies_even_to_an_identical_spectrum():
    """参照窓が 4 個以下なら同一でも 1.0 にならない。penalty は「似ていなさ」ではなく
    「参照の情報量の乏しさ」への減点なので、一致度とは独立に掛かる。"""
    three = [[100.0, 999.0], [200.0, 800.0], [300.0, 600.0]]
    assert sm.weighted_dot_product(three, three, bin_width=0.01) == pytest.approx(0.94, abs=1e-9)
    assert sm.simple_dot_product(three, three, bin_width=0.01) == pytest.approx(1.0, abs=1e-9)


def test_a_disjoint_spectrum_scores_zero():
    other = [[150.0, 999.0], [250.0, 500.0]]
    assert sm.simple_dot_product(_A, other, bin_width=0.01) == pytest.approx(0.0, abs=1e-9)


def test_an_empty_spectrum_returns_the_not_computed_sentinel():
    """上流は比較不能を 0 ではなく -1 で返す。0（＝合わない）と混同しない。"""
    assert sm.simple_dot_product([], _A, bin_width=0.01) == -1.0
    assert sm.weighted_dot_product(_A, [], bin_width=0.01) == -1.0
    assert sm.matched_peaks_scores(_A, [], bin_width=0.01) == (-1.0, -1.0)
    assert sm.spectral_entropy_similarity([], [], bin_width=0.01) == -1.0


def test_a_single_peak_reference_is_penalised():
    """正規化強度 > 0.1 の参照窓が 1 個なら penalty 0.75（weighted / reverse のみ）。"""
    one = [[100.0, 999.0]]
    assert sm.weighted_dot_product(one, one, bin_width=0.01) == pytest.approx(0.75, abs=1e-9)
    assert sm.reverse_dot_product(one, one, bin_width=0.01) == pytest.approx(0.75, abs=1e-9)
    # simple には penalty が無い。
    assert sm.simple_dot_product(one, one, bin_width=0.01) == pytest.approx(1.0, abs=1e-9)


def test_the_penalty_ladder_matches_upstream():
    for n, expected in ((1, 0.75), (2, 0.88), (3, 0.94), (4, 0.97), (5, 1.0)):
        spec = [[100.0 + 10 * i, 999.0] for i in range(n)]
        assert sm.weighted_dot_product(spec, spec, bin_width=0.01) == pytest.approx(expected, abs=1e-9)


def test_reverse_ignores_measured_peaks_absent_from_the_reference():
    """reverse は参照グリッドだけを歩く。測定側の余分なピークは効かない。"""
    reference = [[100.0, 999.0], [200.0, 999.0], [300.0, 999.0], [400.0, 999.0], [500.0, 999.0]]
    measured = reference + [[777.0, 999.0]]
    assert sm.reverse_dot_product(measured, reference, bin_width=0.01) == pytest.approx(
        sm.reverse_dot_product(reference, reference, bin_width=0.01), abs=1e-9)


def test_simple_is_reduced_by_an_extra_measured_peak():
    """simple は両方のグリッドを歩くので、測定側の余分なピークが効く。"""
    reference = [[100.0, 999.0], [200.0, 999.0]]
    measured = reference + [[777.0, 999.0]]
    assert sm.simple_dot_product(measured, reference, bin_width=0.01) < sm.simple_dot_product(
        reference, reference, bin_width=0.01)


def test_matched_peaks_counts_reference_windows_above_one_percent():
    reference = [[100.0, 999.0], [200.0, 999.0], [300.0, 1.0]]   # 3 本目は 1% 未満
    measured = [[100.0, 10.0]]
    percentage, count = sm.matched_peaks_scores(measured, reference, bin_width=0.01)
    assert count == 1.0
    assert percentage == pytest.approx(0.5)     # libCounter は 2


def test_entropy_similarity_of_an_identical_spectrum_is_one():
    assert sm.spectral_entropy_similarity(_A, _A, bin_width=0.01) == pytest.approx(1.0, abs=1e-9)


def test_match_spectrum_returns_the_square_rooted_values():
    """mzTab に出ているのは平方根側。比較の土俵を合わせる（spec §5.1）。"""
    result = sm.match_spectrum(_A, _A, ms2_tol=0.01)
    assert result["simple_dot_product"] == pytest.approx(
        math.sqrt(sm.simple_dot_product(_A, _A, bin_width=0.01)), abs=1e-9)
    assert result["matched_peaks_count"] == 5.0
    assert result["entropy_similarity"] == pytest.approx(1.0, abs=1e-9)


def test_the_not_computed_sentinel_survives_the_square_root():
    """-1 は「比較していない」。sqrt(-1) の NaN に化けさせない。"""
    result = sm.match_spectrum([], _A, ms2_tol=0.01)
    assert result["simple_dot_product"] == -1.0
    assert result["weighted_dot_product"] == -1.0
    assert result["reverse_dot_product"] == -1.0


def test_match_spectrum_reports_the_alignment():
    """数値には出ない「どの窓が合ったか」。対向プロットの注釈がこれを読む。"""
    reference = [[100.0, 999.0], [200.0, 500.0]]
    measured = [[100.0, 999.0]]
    alignment = sm.match_spectrum(measured, reference, ms2_tol=0.01)["alignment"]
    assert [(round(a["mz"], 3), a["matched"]) for a in alignment] == [(100.0, True), (200.0, False)]


def test_unsorted_input_is_sorted_before_scoring():
    """上流の走査は m/z 昇順を前提にしている。並べ替えずに渡すと静かに壊れる。"""
    shuffled = list(reversed(_A))
    assert sm.simple_dot_product(shuffled, _A, bin_width=0.01) == pytest.approx(1.0, abs=1e-9)


# --------------------------------------------------------------------------
# Important 4（最終レビュー）: 縮退した参照/測定スペクトルで例外を投げない。
# --------------------------------------------------------------------------
def test_entropy_similarity_with_a_reference_containing_a_zero_intensity_peak():
    """`spectral_entropy_similarity([[100,10]], [[100,0],[200,5]])` は以前
    `math.log2(0)` の定義域エラー（ValueError）を投げていた。総強度は非ゼロ
    （0+5=5）だが個々のピーク強度に 0 が混ざるケース。0*log2(0)=0 として
    寄与なし扱いにするだけで、比較不能（-1）にはしない。"""
    result = sm.spectral_entropy_similarity([[100.0, 10.0]], [[100.0, 0.0], [200.0, 5.0]], bin_width=0.01)
    assert result != -1.0
    assert not math.isnan(result)


def test_match_spectrum_with_an_all_zero_reference_does_not_raise():
    """`match_spectrum([[100,10]], [[100,0]])` は以前 entropy 計算内で
    ZeroDivisionError を投げていた（参照の総強度が 0）。空ではないので `-1`
    （比較不能）ではなく、他の 3 dot product が同じ状況（base_r==0）で返す
    `0.0`（＝合わなかった）に揃える。"""
    result = sm.match_spectrum([[100.0, 10.0]], [[100.0, 0.0]], ms2_tol=0.01)
    assert result["entropy_similarity"] == 0.0
    assert result["weighted_dot_product"] == 0.0
    assert result["simple_dot_product"] == 0.0
    assert result["reverse_dot_product"] == 0.0


def test_entropy_similarity_with_an_all_zero_measured_spectrum_does_not_raise():
    """測定側が縮退（総強度 0）の場合も同様に例外を投げず 0.0 を返す。"""
    result = sm.spectral_entropy_similarity([[100.0, 0.0]], [[100.0, 10.0], [200.0, 5.0]], bin_width=0.01)
    assert result == 0.0


def test_entropy_similarity_still_returns_the_sentinel_for_a_genuinely_empty_spectrum():
    """縮退（総強度 0 だが非空）と比較不能（空）を混同しない回帰テスト。"""
    assert sm.spectral_entropy_similarity([], [[100.0, 10.0]], bin_width=0.01) == -1.0
    assert sm.spectral_entropy_similarity([[100.0, 0.0]], [], bin_width=0.01) == -1.0


# --------------------------------------------------------------------------
# MS-DIAL の総合スコア（順位付けの正準）。
# 上流: `MsScanMatching.GetTotalScore` / `GetGaussianSimilarity`、
#       `MolecularFormulaUtility.FixMassTolerance`。
# --------------------------------------------------------------------------
_SCORES = {
    "simple_dot_product": 0.9, "weighted_dot_product": 0.96,
    "reverse_dot_product": 0.93, "matched_peaks_percentage": 1.0,
    "matched_peaks_count": 9,
}
_SPECTRUM_TERM = (0.96 + 0.9 + 0.93) / 3


def test_gaussian_similarity_peaks_at_one_when_the_values_agree():
    assert sm.gaussian_similarity(8.5, 8.5, 2.0) == 1.0


def test_gaussian_similarity_follows_the_upstream_formula():
    # exp(-0.5 * ((1 - 2) / 1)^2)
    assert sm.gaussian_similarity(1.0, 2.0, 1.0) == pytest.approx(0.6065306597126334)


def test_gaussian_similarity_returns_the_sentinel_for_missing_values():
    """`0` は「まったく似ていない」を意味する。欠測はそれと区別して `-1` を返す
    （上流の `out bool` 版と同じ規約。dot product 3 種の番兵とも揃う）。"""
    assert sm.gaussian_similarity(None, 1.0, 1.0) == -1.0
    assert sm.gaussian_similarity(1.0, None, 1.0) == -1.0
    assert sm.gaussian_similarity(0.0, 1.0, 1.0) == -1.0
    assert sm.gaussian_similarity(1.0, -1.0, 1.0) == -1.0


def test_mass_tolerance_is_left_alone_at_or_below_500():
    assert sm.fix_mass_tolerance(0.01, 400.0) == 0.01
    assert sm.fix_mass_tolerance(0.01, 500.0) == 0.01


def test_mass_tolerance_is_stretched_by_ppm_above_500():
    """500 での 0.01 は 20 ppm。その ppm を実測 m/z へ当て直す（上流 `FixMassTolerance`）。"""
    assert sm.fix_mass_tolerance(0.01, 1000.0) == pytest.approx(0.02)
    assert sm.fix_mass_tolerance(0.01, 885.54985) == pytest.approx(0.0177109970)


def test_total_score_is_an_unnormalised_sum_of_the_upstream_terms():
    """**1 に収まらない**。上流 `GetTotalScore` は平均でなく和なので、
    0〜1 に正規化すると別の量になり mzTab とも比較できなくなる。"""
    out = sm.total_score(
        _SCORES, precursor_mz=400.0, reference_precursor_mz=400.0, ms1_tol=0.01)

    assert out["mass_similarity"] == 1.0
    assert out["spectrum_score"] == pytest.approx(_SPECTRUM_TERM)
    assert out["rt_similarity"] == -1.0          # use_rt=False の既定
    assert out["total_score"] == pytest.approx(1.0 + _SPECTRUM_TERM + 1.0)
    assert out["total_score"] > 1.0


def test_total_score_adds_the_rt_term_only_when_the_library_says_to_use_it():
    """`IsUseTimeForAnnotationScoring` は既定 `False`。`.dbs` の実値で決まる。"""
    common = dict(precursor_mz=400.0, reference_precursor_mz=400.0, ms1_tol=0.01,
                  rt=8.0, reference_rt=8.0, rt_tol=2.0)
    off = sm.total_score(_SCORES, use_rt=False, **common)
    on = sm.total_score(_SCORES, use_rt=True, **common)

    assert off["rt_similarity"] == -1.0
    assert on["rt_similarity"] == 1.0
    assert on["total_score"] == pytest.approx(off["total_score"] + 1.0)


def test_total_score_skips_the_rt_term_when_the_reference_has_no_rt():
    """`.msp` 由来のレコードは RT を持たない。欠測を `0` として加算すると
    「RT が合わない」と読めてしまうので、項ごと落とす。"""
    out = sm.total_score(
        _SCORES, precursor_mz=400.0, reference_precursor_mz=400.0, ms1_tol=0.01,
        rt=8.0, reference_rt=None, rt_tol=2.0, use_rt=True)

    assert out["rt_similarity"] == -1.0
    assert out["total_score"] == pytest.approx(1.0 + _SPECTRUM_TERM + 1.0)


def test_total_score_does_not_add_the_incomparable_sentinel():
    """`-1`（比較不能）を素朴に足すと総合スコアが下がる。上流は「> 0 のときだけ加算」。"""
    incomparable = {
        "simple_dot_product": -1.0, "weighted_dot_product": -1.0,
        "reverse_dot_product": -1.0, "matched_peaks_percentage": -1.0,
        "matched_peaks_count": -1,
    }
    out = sm.total_score(
        incomparable, precursor_mz=400.0, reference_precursor_mz=400.0, ms1_tol=0.01)

    assert out["spectrum_score"] == -1.0
    assert out["total_score"] == pytest.approx(1.0)   # m/z 項だけが残る


def test_total_score_survives_a_reference_without_a_precursor_mz():
    """store のレコードが precursor m/z を欠いていても例外にしない（番兵で落とす）。"""
    out = sm.total_score(
        _SCORES, precursor_mz=400.0, reference_precursor_mz=None, ms1_tol=0.01)

    assert out["mass_similarity"] == -1.0
    assert out["total_score"] == pytest.approx(_SPECTRUM_TERM + 1.0)


def test_total_score_refuses_to_use_rt_without_a_tolerance():
    """`rt_tol=None` のまま RT を使うと 0 除算になる。黙って縮退させない
    （`store.candidates` の `rt_tol` 欠落と同じ扱い）。"""
    with pytest.raises(ValueError):
        sm.total_score(
            _SCORES, precursor_mz=400.0, reference_precursor_mz=400.0, ms1_tol=0.01,
            rt=8.0, reference_rt=8.0, use_rt=True)


# --------------------------------------------------------------------------
# 足切りの判定（`cutoff_mask`）。採点に入るピークと、図にだけ出るピークを
# 分ける唯一の根拠。`normalize_measured` はこの mask の上に載る。
# --------------------------------------------------------------------------
def test_the_relative_threshold_is_measured_against_the_pre_cutoff_max():
    """分母は**足切り前**の最大強度。上流 `maxIntensity` は足切りループの外で
    一度だけ計算される（生き残ったピークの最大値ではない）。"""
    spectrum = [[100.0, 1000.0], [200.0, 200.0], [300.0, 50.0]]

    mask = sm.cutoff_mask(spectrum, relative_amp_cutoff=0.1, absolute_amp_cutoff=0.0)

    assert mask == [True, True, False]


def test_a_peak_sitting_exactly_on_the_threshold_is_excluded():
    """上流の条件は `>`。境界上は「外」。"""
    spectrum = [[100.0, 1000.0], [200.0, 100.0]]

    mask = sm.cutoff_mask(spectrum, relative_amp_cutoff=0.1, absolute_amp_cutoff=0.0)

    assert mask == [True, False]


def test_normalize_measured_keeps_exactly_the_peaks_the_mask_marks():
    """判定規則を 2 箇所に書かない。`normalize_measured` は mask の上に載る。"""
    spectrum = [[100.0, 1000.0], [200.0, 120.0], [300.0, 20.0]]
    cutoffs = {"relative_amp_cutoff": 0.05, "absolute_amp_cutoff": 30.0}

    kept = sm.normalize_measured(spectrum, **cutoffs)
    mask = sm.cutoff_mask(spectrum, **cutoffs)

    assert [mz for mz, _ in kept] == [mz for (mz, _), keep in zip(spectrum, mask) if keep]


def test_an_empty_spectrum_has_an_empty_mask():
    assert sm.cutoff_mask([], relative_amp_cutoff=0.1, absolute_amp_cutoff=0.0) == []
