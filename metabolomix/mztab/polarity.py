"""DatasetState の極性（ジョブの宣言 → mzTab-M の SML 付加体の多数決）。

`dataset_feature_search` の目標照合と、再解析のハンドル（`dataset_handle_open`）が
同じ規則で極性を決めるための唯一の実装（spec 2026-10-19 §7）。
"""
from __future__ import annotations

from pathlib import Path

__all__ = ["dataset_polarity"]

_POLARITIES = ("positive", "negative")


def dataset_polarity(ds) -> tuple[str | None, str, list[str]]:
    """データセットの極性・その出所（`job_declared` / `adduct_majority` / `unknown`）・注記。

    ジョブの宣言極性を優先し、mzTab-M の SML 付加体の多数決は裏取りに使う（spec §7）。
    両方あって食い違えば宣言を使い、食い違いを注記する。
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

    if declared is None and getattr(ds, "source_format", None) == "msdial_folder":
        # 解析フォルダには mzTab-M の付加体の多数決が無い。.arf2 の ion_mode の多数を使う。
        folder_majority = getattr(ds, "ion_mode_majority", None)
        if folder_majority in _POLARITIES:
            return folder_majority, "arf2_ion_mode", notes
        notes.append("データセットの極性を決められなかった（.arf2 の ion_mode が陽陰で決まらない）。"
                     "目標の polarity による適用判定（not_applicable）はしていない。")
        return None, "unknown", notes

    from metabolomix.console.output_collector import read_adduct_polarity
    majorities = {read_adduct_polarity(path)["adduct_majority"]
                  for path in (getattr(ds, "source_files", None) or {})}
    majorities.discard(None)
    majority = next(iter(majorities)) if len(majorities) == 1 else None
    if len(majorities) > 1:
        notes.append("読み込んだ mzTab-M の付加体の多数決が陽陰に割れたので、多数決は使っていない。")

    if declared is not None:
        if majority is not None and majority != declared:
            notes.append(
                f"ジョブの宣言極性（{declared}）と mzTab-M の付加体の多数決（{majority}）が"
                "食い違う。宣言を使った。")
        return declared, "job_declared", notes
    if majority is not None:
        return majority, "adduct_majority", notes
    notes.append(
        "データセットの極性を決められなかった（ジョブの宣言も mzTab-M の付加体の多数決も無い）。"
        "目標の polarity による適用判定（not_applicable）はしていない。")
    return None, "unknown", notes
