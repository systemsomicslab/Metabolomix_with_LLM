"""アラインメント 1 バッチのファイル群と、スポット ID（MasterAlignmentID）からの引き当て。

MS-DIAL のアラインメント出力は `AlignmentResult_<timestamp>` を共通の stem に持つ
（`.arf2` / `.EIC.aef` / `.dcl`）。スポット ID はこの 3 つで共通で、
アラインメント `.dcl` は `dcl[MasterAlignmentID]` が代表試料の MS/MS
（`curation/evidence.py` と同じ対応）。索引の崩れに備え、precursor m/z が
`PRECURSOR_TOLERANCE` を、RT が `RT_TOLERANCE` を超えて ARF2 と食い違う MS/MS は使わない。

`.dcl` の precursor と RT は**代表試料の実測値**で、ARF2 の `MassCenter`（最大強度の
試料の値）・`RT`（平均）とは一致しない。実データ 6 バッチ・約 2.1 万スポット
（liver / feces / kidney、2026-10-15 測定）で差は最大 0.027 Da（44 ppm）・0.106 分
だったので、許容はその約 2 倍に置く。±0.01 Da では正しい対応の約 1 割を落とした。
RT も見るのは、隣り合う AlignmentID は m/z がほぼ同じことがあり、m/z だけでは
1 つずれた対応を見分けられないため。

deps: arf2.reader（load_catalog）と dcl.reader だけ。tools 層は import しない。
"""

from __future__ import annotations

from pathlib import Path

from metabolomix.arf2.reader import load_catalog
from metabolomix.dcl.reader import deserialize_dcl, read_dcl_entry

PRECURSOR_TOLERANCE = 0.05  # Da
RT_TOLERANCE = 0.2  # min（dcl.reader の既定 RT 窓と同じ）
_ALIGNMENT_SUFFIXES = (".EIC.aef", ".arf2", ".dcl")

_counts_cache: dict[tuple[str, int, int], list[dict]] = {}


def alignment_stem(path: str | Path) -> Path | None:
    """`<dir>/<stem>`（拡張子なし）。アラインメント出力でなければ None。"""
    path = Path(path)
    lowered = path.name.lower()
    for suffix in _ALIGNMENT_SUFFIXES:
        if lowered.endswith(suffix.lower()):
            return path.with_name(path.name[: -len(suffix)])
    return None


def alignment_sibling(path: str | Path, suffix: str) -> Path | None:
    """同じバッチの `<stem><suffix>`。無ければ None（別バッチへは落ちない）。"""
    stem = alignment_stem(path)
    if stem is None:
        return None
    candidate = stem.with_name(stem.name + suffix)
    return candidate if candidate.is_file() else None


def spot_record(arf2_path: str | Path, spot_id: int) -> dict | None:
    """ARF2 で `AlignmentID == spot_id` のレコード。無ければ None。"""
    for record in load_catalog(Path(arf2_path)):
        if record.get("AlignmentID") == spot_id:
            return record
    return None


def msms_peak_counts(dcl_path: str | Path) -> list[dict]:
    """アラインメント `.dcl` 全件の `{precursor_mz, rt, n_peaks}`（スペクトルは読まない）。

    一覧表示用。ファイルの mtime とサイズでキャッシュする。
    """
    path = Path(dcl_path)
    stat = path.stat()
    key = (str(path.resolve()), stat.st_mtime_ns, stat.st_size)
    cached = _counts_cache.get(key)
    if cached is None:
        cached = [
            {"precursor_mz": entry["precursor_mz"], "rt": entry["rt"],
             "n_peaks": entry["n_msms_peaks"]}
            for entry in deserialize_dcl(str(path), include_spectrum=False)
        ]
        _counts_cache.clear()
        _counts_cache[key] = cached
    return cached


def precursor_agrees(precursor_mz: float, expected_mz: float | None,
                     rt: float | None = None, expected_rt: float | None = None) -> bool:
    """`.dcl` の (precursor, RT) が ARF2 のスポットと同じとみなせるか。

    期待値が無い（None、または m/z が欠損の -1）軸は判定しない。
    """
    if expected_mz is not None and expected_mz >= 0:
        if abs(float(precursor_mz) - float(expected_mz)) > PRECURSOR_TOLERANCE:
            return False
    if rt is not None and expected_rt is not None:
        if abs(float(rt) - float(expected_rt)) > RT_TOLERANCE:
            return False
    return True


def alignment_msms(dcl_path: str | Path, spot_id: int, *, expected_mz: float | None,
                   expected_rt: float | None = None) -> dict:
    """スポットの MS/MS を引く。`status` で結果を区別する。

    - `ok`: `spectrum`（全ピーク）・`n_peaks`・`precursor_mz`・`rt` を持つ
    - `no_msms`: エントリはあるがピークが無い（MS/MS 未取得）
    - `precursor_mismatch`: precursor か RT が期待値と合わない（索引対応の崩れ）
    - `out_of_range`: `.dcl` にその番号のエントリが無い
    """
    try:
        entry = read_dcl_entry(str(dcl_path), spot_id)
    except IndexError:
        return {"status": "out_of_range"}
    result = {
        "precursor_mz": entry["precursor_mz"],
        "rt": entry["rt"],
        "n_peaks": entry["n_msms_peaks"],
    }
    if not precursor_agrees(entry["precursor_mz"], expected_mz, entry["rt"], expected_rt):
        return {"status": "precursor_mismatch", **result}
    if not entry["msms_spectrum"]:
        return {"status": "no_msms", **result}
    return {"status": "ok", **result, "spectrum": entry["msms_spectrum"]}
