"""MS/MS スペクトルの照合スコア（純関数）。

**上流 `MsScanMatching.cs` の写しである。** 定義を揃えているのは mzTab-M の
`id_confidence_measure[4..8]` と同じ土俵で数値を比較するためで、次の「明らかに変な点」は
**意図的にそのまま残している**。直すと比較の土俵が消える（spec §6.3）。

- 上流には `wM` / `wR` の計算があるがどこにも使われていない。移植していない。
- 上流の weighted には中身が同一の if/else 分岐がある。条件ごと落としている。
- simple の `× 999` は比になる時点で打ち消える。意味は無い。
- 窓の走査は固定幅でなく、**同じピークが隣り合う 2 つの窓に二重計上されうる**。
  素直な 1 対 1 アラインメントに書き直してはいけない。
- entropy に Li et al. 2021 の低エントロピー重み変換は入っていない。足さない。

比較不能（どちらかのスペクトルが空）は **0 ではなく -1** を返す。0 は「合わない」を
意味するので、混同すると「照合していない」が「合わなかった」に化ける。

非空だが総強度が 0 の縮退スペクトル（例: `[[100.0, 0.0]]`）は「比較不能」ではない
（空ではないため）ので `-1` にはしない。dot product 3 種は窓合算後の信号が無い
（`base_m == 0` / `base_r == 0`）ときに既に `0.0` を返しており、entropy も同じ
規約に揃えて `0.0` を返す（最終レビュー Important 4。以前は `_entropy` が
ゼロ除算・`log2(0)` の定義域エラーで例外を投げていた）。

前提: 入力は m/z 昇順。`_prepare()` が並べ替えるので呼び側は気にしなくてよい。

上流ソース: `MsdialWorkbench` (master, HEAD afd5f9522)
`src/Common/CommonStandard/Algorithm/Scoring/MsScanMatching.cs`
  - `GetWeightedDotProduct`  (~4244行目)
  - `GetSimpleDotProduct`    (~4365行目)
  - `GetReverseDotProduct`   (~3958行目)
  - `GetMatchedPeaksScores`  (~731行目)
  - `GetSpectralEntropySimilarity` / `GetSpectralEntropy` (~795行目)
  - `IsComparedAvailable`    (~26行目)
前処理: `src/MSDIAL5/MsdialCore/Utility/DataAccess.cs` の `GetNormalizedMs2Spectra`。
"""
import math

__all__ = [
    "weighted_dot_product",
    "simple_dot_product",
    "reverse_dot_product",
    "matched_peaks_scores",
    "spectral_entropy_similarity",
    "cutoff_mask",
    "normalize_measured",
    "match_spectrum",
    "gaussian_similarity",
    "fix_mass_tolerance",
    "total_score",
]

_PEAK_COUNT_PENALTY = {1: 0.75, 2: 0.88, 3: 0.94, 4: 0.97}


def _prepare(spectrum):
    """None を空に潰し、m/z 昇順に並べ替えた [[mz, intensity], ...] を返す。"""
    if not spectrum:
        return []
    return sorted(([float(p[0]), float(p[1])] for p in spectrum), key=lambda p: p[0])


def _is_compared_available(peaks1, peaks2):
    """上流 `IsComparedAvailable`: どちらも None でなく、どちらも要素数 0 でないこと。
    ここでは `_prepare()` 済みの入力を受け取るので None は既に [] に潰れている。"""
    return bool(peaks1) and bool(peaks2)


def _scan_sticky(peaks, start_index, focused_mz, bin_width):
    """weighted / reverse / matched-peaks が使う走査。

    上流の `for (int i = remaindIndexM; i < peaks1.Count; i++) { ... else { remaindIndexM = i; break; } }`
    の写し。カーソル（返す index）は「focusedMz + bin 以上の質量」に出会って break したときだけ
    進む。ループが break せずに配列末尾まで達した場合、カーソルは呼び出し時の値のまま変わらない
    （＝次の外側ループでも同じ位置から再スキャンする。これが「隣接窓への二重計上」の元）。
    """
    total = 0.0
    n = len(peaks)
    i = start_index
    while i < n:
        mz = peaks[i][0]
        if mz < focused_mz - bin_width:
            i += 1
            continue
        if mz < focused_mz + bin_width:
            total += peaks[i][1]
            i += 1
            continue
        return total, i
    return total, start_index


def _scan_advancing(peaks, start_index, focused_mz, bin_width):
    """simple だけが使う走査。

    上流の `for (int i = remaindIndexM; i < peaks1.Count; remaindIndexM = ++i) { ... else { break; } }`
    の写し。C# の for 文は `continue` でも増分式 `remaindIndexM = ++i` を必ず実行するため、
    ここではカーソルは要素を処理するたび（continue でも）進み、配列末尾まで達すると
    `len(peaks)` になり得る（sticky 版とはここが違う）。simple 側で末尾到達の分岐が
    明示的に要るのはこのため。
    """
    total = 0.0
    n = len(peaks)
    i = start_index
    while i < n:
        mz = peaks[i][0]
        if mz < focused_mz - bin_width:
            i += 1
            continue
        if mz < focused_mz + bin_width:
            total += peaks[i][1]
            i += 1
            continue
        break
    return total, i


def weighted_dot_product(measured, reference, *, bin_width, mass_begin=0.0, mass_end=2000.0):
    """`GetWeightedDotProduct` の写し（**二乗値**）。カーソルは両スペクトルの和集合を進む。"""
    peaks1 = _prepare(measured)
    peaks2 = _prepare(reference)
    if not _is_compared_available(peaks1, peaks2):
        return -1.0

    last1 = peaks1[-1][0]
    last2 = peaks2[-1][0]

    min_mz = min(peaks1[0][0], peaks2[0][0])
    max_mz = max(last1, last2)
    if mass_begin > min_mz:
        min_mz = mass_begin
    if max_mz > mass_end:
        max_mz = mass_end

    focused_mz = min_mz
    idx_m = idx_l = 0
    measured_buf = []
    reference_buf = []
    base_m = base_r = -math.inf

    while focused_mz <= max_mz:
        sum_m, idx_m = _scan_sticky(peaks1, idx_m, focused_mz, bin_width)
        sum_r, idx_l = _scan_sticky(peaks2, idx_l, focused_mz, bin_width)

        # 上流に `if (sumM <= 0 && sumR > 0) {...} else {...}` の分岐があるが両枝の中身は
        # 同一（brief 記載の瑕疵）。分岐ごと落として共通処理だけ残す。
        measured_buf.append([focused_mz, sum_m])
        if sum_m > base_m:
            base_m = sum_m
        reference_buf.append([focused_mz, sum_r])
        if sum_r > base_r:
            base_r = sum_r

        if focused_mz + bin_width > max(last1, last2):
            break
        next_m = peaks1[idx_m][0]
        next_l = peaks2[idx_l][0]
        if focused_mz + bin_width > next_l and focused_mz + bin_width <= next_m:
            focused_mz = next_m
        elif focused_mz + bin_width <= next_l and focused_mz + bin_width > next_m:
            focused_mz = next_l
        else:
            focused_mz = min(next_m, next_l)

    if base_m == 0 or base_r == 0:
        return 0.0

    sum_measure = 0.0
    sum_reference = 0.0
    l_counter = 0
    for entry_m, entry_r in zip(measured_buf, reference_buf):
        entry_m[1] = entry_m[1] / base_m
        entry_r[1] = entry_r[1] / base_r
        sum_measure += entry_m[1]
        sum_reference += entry_r[1]
        if entry_r[1] > 0.1:
            l_counter += 1

    penalty = _PEAK_COUNT_PENALTY.get(l_counter, 1.0)
    # 上流はここで wM = 1/(sumMeasure-0.5), wR = 1/(sumReference-0.5) を計算するが、
    # どちらも以降どこにも使われていない（brief 記載の瑕疵）。移植しない。

    cutoff = 0.01
    scalar_m = scalar_r = covariance = 0.0
    for entry_m, entry_r in zip(measured_buf, reference_buf):
        if entry_m[1] < cutoff:
            continue
        scalar_m += entry_m[1] * entry_m[0]
        scalar_r += entry_r[1] * entry_r[0]
        covariance += math.sqrt(entry_m[1] * entry_r[1]) * entry_m[0]

    if scalar_m == 0 or scalar_r == 0:
        return 0.0
    return (covariance ** 2) / scalar_m / scalar_r * penalty


def reverse_dot_product(measured, reference, *, bin_width, mass_begin=0.0, mass_end=2000.0):
    """`GetReverseDotProduct` の写し（**二乗値**）。カーソルは参照グリッドだけを進む。

    weighted との違いは 2 つ: (1) カーソルが `peaks2[idx_l]`（参照）だけで進む、
    (2) 最終合算のカットオフ判定を測定側でなく**参照側**の強度で行う。
    """
    peaks1 = _prepare(measured)
    peaks2 = _prepare(reference)
    if not _is_compared_available(peaks1, peaks2):
        return -1.0

    last2 = peaks2[-1][0]

    min_mz = peaks2[0][0]
    max_mz = last2
    if mass_begin > min_mz:
        min_mz = mass_begin
    if max_mz > mass_end:
        max_mz = mass_end

    focused_mz = min_mz
    idx_m = idx_l = 0
    measured_buf = []
    reference_buf = []
    base_m = base_r = -math.inf

    while focused_mz <= max_mz:
        sum_l, idx_l = _scan_sticky(peaks2, idx_l, focused_mz, bin_width)
        sum_m, idx_m = _scan_sticky(peaks1, idx_m, focused_mz, bin_width)

        # 上流の `if (sumM <= 0) {...} else { counter++; }` は counter を使うが、
        # counter 自体は戻り値に影響しない（バッファへの記録は両枝で同一）ので、
        # 分岐ごと落として共通処理だけ残す。
        measured_buf.append([focused_mz, sum_m])
        if sum_m > base_m:
            base_m = sum_m
        reference_buf.append([focused_mz, sum_l])
        if sum_l > base_r:
            base_r = sum_l

        if focused_mz + bin_width > last2:
            break
        focused_mz = peaks2[idx_l][0]

    if base_m == 0 or base_r == 0:
        return 0.0

    sum_measure = 0.0
    sum_reference = 0.0
    l_counter = 0
    for entry_m, entry_r in zip(measured_buf, reference_buf):
        entry_m[1] = entry_m[1] / base_m
        entry_r[1] = entry_r[1] / base_r
        sum_measure += entry_m[1]
        sum_reference += entry_r[1]
        if entry_r[1] > 0.1:
            l_counter += 1

    penalty = _PEAK_COUNT_PENALTY.get(l_counter, 1.0)
    # wM / wR は weighted と同様、上流でも未使用。移植しない。

    cutoff = 0.01
    scalar_m = scalar_r = covariance = 0.0
    for entry_m, entry_r in zip(measured_buf, reference_buf):
        if entry_r[1] < cutoff:      # weighted と違い、参照側でカットオフを掛ける。
            continue
        scalar_m += entry_m[1] * entry_m[0]
        scalar_r += entry_r[1] * entry_r[0]
        covariance += math.sqrt(entry_m[1] * entry_r[1]) * entry_m[0]

    if scalar_m == 0 or scalar_r == 0:
        return 0.0
    return (covariance ** 2) / scalar_m / scalar_r * penalty


def simple_dot_product(measured, reference, *, bin_width, mass_begin=0.0, mass_end=2000.0):
    """`GetSimpleDotProduct` の写し（**二乗値**）。penalty も m/z 重みも無い。

    カーソルは `_scan_advancing` を使うので、両方の配列末尾に達すると
    `len(peaks)` に等しくなり得る。上流同様、その越境を明示的に処理する。
    """
    peaks1 = _prepare(measured)
    peaks2 = _prepare(reference)
    if not _is_compared_available(peaks1, peaks2):
        return -1.0

    last1 = peaks1[-1][0]
    last2 = peaks2[-1][0]

    max_mz = max(last1, last2)
    if max_mz > mass_end:
        max_mz = mass_end

    idx_m = 0
    while idx_m < len(peaks1) and peaks1[idx_m][0] < mass_begin - bin_width:
        idx_m += 1
    idx_l = 0
    while idx_l < len(peaks2) and peaks2[idx_l][0] < mass_begin - bin_width:
        idx_l += 1

    focused_mz = min(
        peaks1[idx_m][0] if idx_m < len(peaks1) else math.inf,
        peaks2[idx_l][0] if idx_l < len(peaks2) else math.inf,
    )

    measured_buf = []
    reference_buf = []
    base_m = base_r = -math.inf

    while focused_mz <= max_mz:
        sum_m, idx_m = _scan_advancing(peaks1, idx_m, focused_mz, bin_width)
        sum_r, idx_l = _scan_advancing(peaks2, idx_l, focused_mz, bin_width)

        measured_buf.append([focused_mz, sum_m])
        if sum_m > base_m:
            base_m = sum_m
        reference_buf.append([focused_mz, sum_r])
        if sum_r > base_r:
            base_r = sum_r

        if focused_mz + bin_width > max(last1, last2):
            break
        if idx_m >= len(peaks1) or idx_l >= len(peaks2):
            focused_mz = peaks1[idx_m][0] if idx_l >= len(peaks2) else peaks2[idx_l][0]
            continue
        next_m = peaks1[idx_m][0]
        next_l = peaks2[idx_l][0]
        if focused_mz + bin_width > next_l and focused_mz + bin_width <= next_m:
            focused_mz = next_m
        elif focused_mz + bin_width <= next_l and focused_mz + bin_width > next_m:
            focused_mz = next_l
        else:
            focused_mz = min(next_m, next_l)

    if base_m == 0 or base_r == 0:
        return 0.0

    for entry in measured_buf:
        # 上流の `* 999` は比になる時点で打ち消えるが、忠実に写す（意味は無い）。
        entry[1] = entry[1] / base_m * 999
    for entry in reference_buf:
        entry[1] = entry[1] / base_r * 999

    scalar_m = sum(entry[1] for entry in measured_buf)
    scalar_r = sum(entry[1] for entry in reference_buf)
    covariance = sum(
        math.sqrt(entry_m[1] * entry_r[1]) for entry_m, entry_r in zip(measured_buf, reference_buf)
    )

    if scalar_m == 0 or scalar_r == 0:
        return 0.0
    return (covariance ** 2) / scalar_m / scalar_r


def _matched_peaks_walk(measured, reference, bin_width, mass_begin, mass_end):
    """`GetMatchedPeaksScores` の走査本体。matched_peaks_scores と match_spectrum の
    alignment が共有する。比較不能なら (None, None, []) を返す。"""
    peaks1 = _prepare(measured)
    peaks2 = _prepare(reference)
    if not _is_compared_available(peaks1, peaks2):
        return None, None, []

    last2 = peaks2[-1][0]

    min_mz = peaks2[0][0]
    max_mz = last2
    if mass_begin > min_mz:
        min_mz = mass_begin
    if max_mz > mass_end:
        max_mz = mass_end

    focused_mz = min_mz
    idx_m = idx_l = 0
    counter = 0
    lib_counter = 0
    max_lib_intensity = max(p[1] for p in peaks2)
    records = []

    while focused_mz <= max_mz:
        sum_l, idx_l = _scan_sticky(peaks2, idx_l, focused_mz, bin_width)
        is_lib_hit = sum_l >= 0.01 * max_lib_intensity
        if is_lib_hit:
            lib_counter += 1

        sum_m, idx_m = _scan_sticky(peaks1, idx_m, focused_mz, bin_width)
        matched = sum_m > 0 and is_lib_hit
        if matched:
            counter += 1

        records.append({"mz": focused_mz, "measured": sum_m, "reference": sum_l, "matched": matched})

        if focused_mz + bin_width > last2:
            break
        focused_mz = peaks2[idx_l][0]

    return counter, lib_counter, records


def matched_peaks_scores(measured, reference, *, bin_width, mass_begin=0.0, mass_end=2000.0):
    """`GetMatchedPeaksScores` の写し。戻り値は `(percentage, count)`。"""
    counter, lib_counter, _records = _matched_peaks_walk(measured, reference, bin_width, mass_begin, mass_end)
    if counter is None:
        return (-1.0, -1.0)
    if lib_counter == 0:
        return (0.0, 0.0)
    return (counter / lib_counter, float(counter))


def _normalized_by_total(peaks):
    """`SpectrumHandler.GetNormalizedByTotalIntensityPeaks` の写し。"""
    total = sum(intensity for _, intensity in peaks)
    return [(mz, intensity / total) for mz, intensity in peaks]


def _binned_spectrum(peaks, bin_width, halve=False):
    """`SpectrumHandler.GetBinnedSpectrum` / `GetCombinedSpectrum` が使う質量フレーム
    への集約の写し。`massframe = int(mz / bin_width)`（C# の `(int)` キャストは 0 方向への
    切り捨てだが、質量は正なので floor と同じ）ごとにグループ化し、代表 m/z はそのグループ内で
    強度最大のピークの m/z、強度は合算（`halve=True` のときだけ 0.5 倍。`GetCombinedSpectrum`
    が両スペクトルを混ぜたときの平均化に使う 0.5 倍で、単独スペクトルの binning には掛からない）。
    """
    buckets = {}
    for mz, intensity in peaks:
        frame = int(mz / bin_width)
        buckets.setdefault(frame, []).append((mz, intensity))

    result = []
    for group in buckets.values():
        max_mz = max(group, key=lambda p: p[1])[0]
        total = sum(intensity for _, intensity in group)
        if halve:
            total *= 0.5
        result.append((max_mz, total))
    return result


def _entropy(peaks):
    """`GetSpectralEntropy` の写し。

    強度 0 のピークは `0 * log2(0) = 0`（情報理論の慣例的な極限値）として寄与なし
    に扱う——素の `math.log2(0)` は定義域エラーになる（最終レビュー Important 4:
    `spectral_entropy_similarity([[100,10]], [[100,0],[200,5]])` が
    `ValueError` を投げていた）。合計強度が 0（全ピーク強度 0）の縮退入力は
    `spectral_entropy_similarity` 側で事前にガードするので通常ここには来ないが、
    単体で呼ばれても例外を投げないよう `0.0` を返す。
    """
    total = sum(intensity for _, intensity in peaks)
    if total <= 0:
        return 0.0
    entropy = 0.0
    for _, intensity in peaks:
        if intensity <= 0:
            continue
        p = intensity / total
        entropy -= p * math.log2(p)
    return entropy


def spectral_entropy_similarity(measured, reference, *, bin_width):
    """`GetSpectralEntropySimilarity` の写し。Li et al. 2021 の低エントロピー重み変換は
    入っていない（上流にも無い）。足さない。

    **縮退入力の扱い（最終レビュー Important 4）**: 測定・参照どちらかが
    「非空だが総強度が 0」（例: `[[100.0, 0.0]]`）だと `_normalized_by_total` が
    ゼロ除算で例外を投げていた（上流は同じ状況で NaN になる）。この関数の既存の
    番兵規約では `-1.0` は「比較不能＝どちらかが空」専用（`_is_compared_available`）
    であり、この入力は空ではないので `-1.0` を流用すると意味が変わる。
    weighted/reverse/simple の 3 dot product が同じ「窓を合算したら信号が無かった」
    状況（`base_m == 0` / `base_r == 0`）で `0.0`（＝合わなかった）を返しているのに
    倣い、ここも `0.0` を返す（controller裁定 2026-09-20: 上流の NaN をそのまま
    模倣する先例が本モジュールに無い一方、「総強度 0 は無情報」を「合わなかった」
    として扱う `0.0` は同モジュール内の既存規約と一貫する）。
    """
    peaks1 = _prepare(measured)
    peaks2 = _prepare(reference)
    if not _is_compared_available(peaks1, peaks2):
        return -1.0

    total1 = sum(intensity for _, intensity in peaks1)
    total2 = sum(intensity for _, intensity in peaks2)
    if total1 <= 0 or total2 <= 0:
        return 0.0

    norm1 = _normalized_by_total(peaks1)
    norm2 = _normalized_by_total(peaks2)
    combined = _binned_spectrum(norm1 + norm2, bin_width, halve=True)

    entropy12 = _entropy(combined)
    entropy1 = _entropy(_binned_spectrum(peaks1, bin_width))
    entropy2 = _entropy(_binned_spectrum(peaks2, bin_width))

    return 1 - (2 * entropy12 - entropy1 - entropy2) * 0.5


def cutoff_mask(spectrum, *, relative_amp_cutoff=0.0, absolute_amp_cutoff=0.0):
    """足切りを通るピークを入力と同じ並びの真偽リストで返す。

    判定規則そのもの（`>` の厳密不等号・分母は**足切り前**の最大強度）はここにしか
    無い。`normalize_measured` はこの mask の上に載り、描画層は「図には出るが採点に
    入らなかったピーク」を同じ mask から割り出す——規則を 2 箇所に書くと、片方だけが
    上流の変更に追従して静かにずれる（`defaults.py` を作ったのと同じ理由）。
    """
    if not spectrum:
        return []
    intensities = [float(p[1]) for p in spectrum]
    max_intensity = max(intensities)
    return [
        intensity > max_intensity * relative_amp_cutoff and intensity > absolute_amp_cutoff
        for intensity in intensities
    ]


def normalize_measured(spectrum, *, relative_amp_cutoff=0.0, absolute_amp_cutoff=0.0):
    """`DataAccess.GetNormalizedMs2Spectra` の写し。**測定側にだけ**掛ける前処理。

    1. 足切り: `intensity > max(original) * relative_amp_cutoff and intensity > absolute_amp_cutoff`
       を満たすピークだけを残す。**並べ替えはしない**（元の並び順を保つ。上流も同じ）。
       閾値の分母は足切り前の元の最大強度であって、生き残ったピークの最大強度ではない
       （上流 `maxIntensity` は足切りループの外で一度だけ計算される）。
    2. 再スケール: 残ったピークを `intensity / max(original) * 100` にする。

    **再スケールは本モジュールの 5 種のスコア関数すべてに対して数学的に no-op である。**
    dot product 3 種は各窓の合算値を自分自身の最大値で割ってから使うので測定強度の定数倍は
    打ち消える。matched peaks は `sumM > 0` の真偽しか見ない。entropy は総和で正規化する。
    どの経路も入力を定数倍しても結果は変わらない。**忠実性のため実装はするが、数値が
    合わないときにここ（再スケール）を疑って時間を溶かさないこと** — 実際に効くのは
    足切りのほうだけである。

    既定 `relative_amp_cutoff=0.0` / `absolute_amp_cutoff=0.0`
    （`MsRefSearchParameterBase` の Key 7 / 8 の既定値）。空スペクトルは空のまま返す。
    """
    if not spectrum:
        return []
    peaks = [(float(p[0]), float(p[1])) for p in spectrum]
    max_intensity = max(intensity for _, intensity in peaks)
    mask = cutoff_mask(spectrum, relative_amp_cutoff=relative_amp_cutoff,
                       absolute_amp_cutoff=absolute_amp_cutoff)

    return [
        [mz, intensity / max_intensity * 100.0]
        for (mz, intensity), keep in zip(peaks, mask) if keep
    ]


def _sqrt_or_sentinel(value):
    """-1（比較不能の番兵）はそのまま通す。sqrt(-1) の NaN に化けさせない。"""
    if value < 0:
        return -1.0
    return math.sqrt(value)


def match_spectrum(measured, reference, *, ms2_tol, mass_begin=0.0, mass_end=2000.0,
                    relative_amp_cutoff=0.0, absolute_amp_cutoff=0.0):
    """5 種のスコアを一括計算する。前処理（`normalize_measured`）は測定側にのみ掛ける。

    dot product 3 種は mzTab の `id_confidence_measure` と揃えるため平方根を返す
    （このモジュールの単体関数は二乗値を返す。spec §5.1）。`-1` の番兵は sqrt を経ても
    そのまま `-1` で通す。`alignment` は参照グリッドの窓ごとの記録
    （`{"mz", "measured", "reference", "matched"}`）で、対向プロットの注釈に使う。

    **`-1` は mzTab の生値とは一致しない（比較時に注意）。**
    この関数は比較不能（測定・参照どちらかが空）を `-1` で返す。`0` は「合わなかった」を
    意味するので、意図的に区別している。ところが上流が mzTab に書き出す値はこの区別を
    保っていない:

    - `id_confidence_measure[4..6]`（simple / weighted / reverse dot product）は
      `MsScanMatchResult` の非二乗 getter（`WeightedDotProduct` など）が
      `Math.Sqrt(Math.Max(Squared*, 0f))` として `-1` を `0` にクランプしてから
      `sqrt` を取る（`MsScanMatchResult.cs:36-38,44-46,51-53`）。`MztabFormatExport.cs`
      はこのクランプ済み getter をそのまま書き出す（`:1412-1414` 付近）。
      **したがって mzTab 上のこれら 3 列では、比較不能は `-1` ではなく `0` として出る。**
    - `id_confidence_measure[7..8]`（`matched_peaks_count` / `matched_peaks_percentage`）
      は素のフィールドで、クランプを経ないので `-1` のまま出る。

    Task 7 で mzTab の値と突き合わせる際は、**dot product 3 種についてのみ `-1` と `0`
    を同一視**すること。同一視せずに素朴に比較すると、比較不能ケース（測定または参照が
    空）で「存在しない乖離」を検出してしまう。matched peaks の 2 列はそのまま `-1` 同士で
    比較してよい。
    """
    normalized = normalize_measured(
        measured, relative_amp_cutoff=relative_amp_cutoff, absolute_amp_cutoff=absolute_amp_cutoff)

    simple = simple_dot_product(normalized, reference, bin_width=ms2_tol,
                                 mass_begin=mass_begin, mass_end=mass_end)
    weighted = weighted_dot_product(normalized, reference, bin_width=ms2_tol,
                                     mass_begin=mass_begin, mass_end=mass_end)
    reverse = reverse_dot_product(normalized, reference, bin_width=ms2_tol,
                                   mass_begin=mass_begin, mass_end=mass_end)
    percentage, count = matched_peaks_scores(normalized, reference, bin_width=ms2_tol,
                                              mass_begin=mass_begin, mass_end=mass_end)
    entropy = spectral_entropy_similarity(normalized, reference, bin_width=ms2_tol)
    _, _, alignment = _matched_peaks_walk(normalized, reference, ms2_tol, mass_begin, mass_end)

    return {
        "simple_dot_product": _sqrt_or_sentinel(simple),
        "weighted_dot_product": _sqrt_or_sentinel(weighted),
        "reverse_dot_product": _sqrt_or_sentinel(reverse),
        "matched_peaks_percentage": percentage,
        "matched_peaks_count": count,
        "entropy_similarity": entropy,
        "alignment": alignment,
    }


def gaussian_similarity(actual, reference, tolerance):
    """RT / precursor m/z の一致度。上流 `MsScanMatching.GetGaussianSimilarity`。

    `exp(-0.5 * ((actual - reference) / tolerance)^2)`。出典は上流の
    docstring が挙げる Tsugawa, H. et al. Anal. Chem. 85, 5191-5199 (2013)。

    **欠測は `0` ではなく `-1` を返す**（上流の `out bool` 版と同じ規約）。
    `0` は「まったく似ていない」を意味するので、混同すると「値が無い」が
    「一致しない」に化ける——このモジュールの dot product 3 種の番兵と同じ考え方。
    上流は非正の値（`<= 0`）も欠測として扱う（RT も m/z も正のはずのため）ので、
    ここもそれに倣う。
    """
    if actual is None or reference is None:
        return -1.0
    if actual <= 0 or reference <= 0:
        return -1.0
    return math.exp(-0.5 * ((actual - reference) / tolerance) ** 2)


def fix_mass_tolerance(tolerance, mass):
    """precursor m/z の許容幅を高質量側へ伸ばす。上流
    `MolecularFormulaUtility.FixMassTolerance`。

    500 以下はそのまま。500 超は「500 における `tolerance` が何 ppm か」を求め直し、
    その ppm を実測 m/z へ当てる（絶対幅ではなく相対幅で効かせる）。
    上流の `PpmCalculator` は小数 4 桁で丸めるので、その丸めも写す——
    丸めないと高質量側で許容幅がわずかにずれ、`gaussian_similarity` の値が
    mzTab と合わなくなる。
    """
    if mass <= 500:
        return tolerance
    ppm = abs(round(((500.0 + tolerance) - 500.0) / 500.0 * 1000000, 4))
    return ppm * mass / 1000000.0


def total_score(scores, *, precursor_mz, reference_precursor_mz, ms1_tol,
                rt=None, reference_rt=None, rt_tol=None, use_rt=False):
    """MS-DIAL の総合スコア。上流 `MsScanMatching.GetTotalScore`。

    **正規化されていない和である。**

        RtSimilarity + AcurateMassSimilarity + (Weighted + Simple + Reverse)/3
        + MatchedPeaksPercentage

    各項は **`> 0` のときだけ加算する**（`-1` の番兵＝比較不能・欠測を足して
    総合スコアを下げないため。上流も同じ条件で加算する）。0〜1 に正規化しては
    いけない——平均にすると別の量になり、MS-DIAL の順位と比較できなくなる。

    上流の `GetTotalScore` はこのほかに `CcsSimilarity` / `IsotopeSimilarity` /
    `AndromedaScore` も足すが、**この経路には存在しない**ので扱わない
    （CCS は IM-MS、isotope は MS1 の同位体パターン、Andromeda はプロテオミクス）。

    `use_rt` は `.dbs` の `IsUseTimeForAnnotationScoring`（`Key(16)`）に対応する。
    **上流の既定は `False`** で、`.msp` のようにフラグを持たないライブラリでは
    RT 項が入らない。`rt_tol` は `use_rt=True` のとき必須——`None` のまま
    使うと 0 除算になるので、黙って縮退させず例外にする
    （`store.candidates` の `rt_tol` 欠落と同じ扱い）。

    上流が順位付けに使うのはこの値**単独ではない**（`MsScanMatchResultContainer.
    ResultOrder` は `IsReferenceMatched` 等を先に見る辞書順タプル）。その
    ブール判定は脂質クラス固有の判定（`Lipidomics/` 66,932 行）に依存するため
    移植していない。詳細は `docs/output_format/library.md` §14.9。
    """
    if use_rt and rt_tol is None:
        raise ValueError(
            "use_rt=True のときは rt_tol を指定してください "
            "（None のままだと gaussian_similarity が 0 除算になります）。"
        )

    rt_similarity = gaussian_similarity(rt, reference_rt, rt_tol) if use_rt else -1.0

    if precursor_mz is None or precursor_mz <= 0:
        mass_similarity = -1.0
    else:
        mass_similarity = gaussian_similarity(
            precursor_mz, reference_precursor_mz, fix_mass_tolerance(ms1_tol, precursor_mz))

    weighted = scores.get("weighted_dot_product", -1.0)
    simple = scores.get("simple_dot_product", -1.0)
    reverse = scores.get("reverse_dot_product", -1.0)
    spectrum_score = (weighted + simple + reverse) / 3.0
    percentage = scores.get("matched_peaks_percentage", -1.0)

    total = 0.0
    if rt_similarity > 0:
        total += rt_similarity
    if mass_similarity > 0:
        total += mass_similarity
    if weighted > 0:          # 上流も weighted だけを見て 3 種の平均を足す
        total += spectrum_score
    if percentage > 0:
        total += percentage

    return {
        "total_score": total,
        "rt_similarity": rt_similarity,
        "mass_similarity": mass_similarity,
        "spectrum_score": spectrum_score,
    }
