"""DatasetState の極性（ジョブの宣言 → 測定側の多数）。

測定側の多数は、mzTab-M のアダプタでは SML 付加体の多数決、解析フォルダのアダプタでは
`.arf2` の ion_mode の多数（`ion_mode_majority`）。どちらのアダプタでもジョブの宣言を優先し、
測定側の多数は裏取りに使う。

`dataset_feature_search` の目標照合と、再解析のハンドル（`load_dataset`）が
同じ規則で極性を決めるための唯一の実装（spec 2026-10-19 §7）。
"""
from __future__ import annotations

from pathlib import Path

__all__ = ["dataset_polarity"]

_POLARITIES = ("positive", "negative")

#: 宣言と測定側の多数が食い違ったときの注記（どちらのアダプタでも同じ文面。`evidence` だけが違う）。
_MISMATCH_NOTE = "ジョブの宣言極性（{declared}）と {evidence}（{majority}）が食い違う。宣言を使った。"


def dataset_polarity(ds) -> tuple[str | None, str, list[str]]:
    """データセットの極性・その出所（`job_declared` / `adduct_majority` / `arf2_ion_mode` /
    `unknown`）・注記。

    ジョブの宣言極性を優先し、測定側の多数（mzTab-M の SML 付加体の多数決、解析フォルダなら
    `.arf2` の ion_mode の多数）は裏取りに使う（spec §7）。両方あって食い違えば宣言を使い、
    食い違いを注記する。
    """
    notes: list[str] = []
    declared = None
    if getattr(ds, "job_path", None):
        from metabolomix.console import job_manager
        try:
            declared = job_manager.load_job(Path(ds.job_path)).polarity
        except (OSError, ValueError, KeyError, TypeError) as exc:
            notes.append(f"ジョブ記録を読めず、宣言極性を使っていない（{exc}）。")
        if declared is not None and declared not in _POLARITIES:
            notes.append(f"ジョブの宣言極性 {declared!r} を解釈できず、使っていない。")
            declared = None

    if getattr(ds, "source_format", None) == "msdial_folder":
        # 解析フォルダには mzTab-M の付加体の多数決が無い。.arf2 の ion_mode の多数を使う。
        majority = getattr(ds, "ion_mode_majority", None)
        if majority not in _POLARITIES:
            majority = None
        evidence, origin = ".arf2 の ion_mode の多数", "arf2_ion_mode"
        unknown = ("データセットの極性を決められなかった（.arf2 の ion_mode が陽陰で決まらない）。"
                   "目標の polarity による適用判定（not_applicable）はしていない。")
    else:
        from metabolomix.console.output_collector import read_adduct_polarity
        majorities = {read_adduct_polarity(path)["adduct_majority"]
                      for path in (getattr(ds, "source_files", None) or {})}
        majorities.discard(None)
        majority = next(iter(majorities)) if len(majorities) == 1 else None
        if len(majorities) > 1:
            notes.append("読み込んだ mzTab-M の付加体の多数決が陽陰に割れたので、多数決は使っていない。")
        evidence, origin = "mzTab-M の付加体の多数決", "adduct_majority"
        unknown = ("データセットの極性を決められなかった（ジョブの宣言も mzTab-M の付加体の多数決も無い）。"
                   "目標の polarity による適用判定（not_applicable）はしていない。")

    if declared is not None:
        if majority is not None and majority != declared:
            notes.append(_MISMATCH_NOTE.format(declared=declared, evidence=evidence,
                                               majority=majority))
        return declared, "job_declared", notes
    if majority is not None:
        return majority, origin, notes
    notes.append(unknown)
    return None, "unknown", notes
