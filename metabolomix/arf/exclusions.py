"""ARF スポットリストに対するサンプル/スポット手動除外（純ロジック層・MCP 非依存）。

metabolomix.analysis.preprocessing / metabolomix.analysis.differential と同格の leaf モジュール。行列を組む直前に
filtered_features から除外集合を適用した派生リストを作るために使う。
サンプルキーは各 AlignedPeakProperties エントリの file_name、スポットキーは
MasterAlignmentID。build_pca_matrix と同じ file_name 導出を用いる。
"""
from __future__ import annotations


def _entry_file_name(entry) -> str | None:
    """AlignedPeakProperties の1エントリ（生 list）から file_name を導出する。

    `reader.file_name_of` は `_convert_to_alignment_feature` と同じ導出を共有する
    ので答えは変わらない。feature dict を組まないぶんだけ速い——ここは spot ×
    注入ぶん呼ばれる（`roster` が spot ごとに全エントリを舐める）。
    """
    from metabolomix.arf.reader import file_name_of
    try:
        return file_name_of(entry)
    except Exception:
        return None


def prune_spots(spots, excluded_samples, excluded_spots):
    """除外集合を適用したスポットリストの非破壊コピーを返す。

    - MasterAlignmentID in excluded_spots のスポットを丸ごと除外。
    - 残スポットの AlignedPeakProperties から file_name in excluded_samples の
      エントリを除去する（スポット dict は浅いコピーし、AlignedPeakProperties を
      フィルタ済みリストへ差し替え。元の spots / エントリは変更しない）。
    - 除外集合が両方空なら入力をそのまま返す（コピー不要・恒等）。
    """
    if not excluded_samples and not excluded_spots:
        return spots
    excluded_samples = set(excluded_samples or ())
    excluded_spots = set(excluded_spots or ())
    out = []
    for spot in spots:
        if spot.get("MasterAlignmentID") in excluded_spots:
            continue
        if not excluded_samples:
            out.append(spot)
            continue
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list):
            out.append(spot)
            continue
        kept = [e for e in aligned if _entry_file_name(e) not in excluded_samples]
        new_spot = dict(spot)
        new_spot["AlignedPeakProperties"] = kept
        out.append(new_spot)
    return out


def roster(spots):
    """現データに存在する (サンプル file_name 集合, MasterAlignmentID 集合) を返す。

    除外指定の未一致検出と list 表示に使う。
    """
    names: set[str] = set()
    ids: set[int] = set()
    for spot in spots or []:
        mid = spot.get("MasterAlignmentID")
        if mid is not None:
            ids.add(mid)
        aligned = spot.get("AlignedPeakProperties")
        if isinstance(aligned, list):
            for e in aligned:
                fn = _entry_file_name(e)
                if fn:
                    names.add(fn)
    return names, ids
