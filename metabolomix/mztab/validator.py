"""mzTab-M 2.0 構造・定量種別バリデータ（純関数）。spec §8, §9 参照。"""
from __future__ import annotations
import re

_HEIGHT_RE = re.compile(r"\bheight\b", re.IGNORECASE)
_AREA_RE = re.compile(r"\barea\b", re.IGNORECASE)
_HEIGHT_PREFIX_RE = re.compile(r"^Height_", re.IGNORECASE)
_AREA_PREFIX_RE = re.compile(r"^Area_", re.IGNORECASE)
_SUPPORTED_VERSION = "2.0.0-M"


def validate_mztab(parse_result: dict) -> dict:
    """mzTab-M 2.0 の構造を検証する。

    戻り値: {"ok": bool, "errors": list[str], "warnings": list[str]}
    """
    errors: list[str] = []
    warnings: list[str] = []
    meta = parse_result.get("metadata", {})
    sections = parse_result.get("sections", {})

    # バージョン検査
    version = meta.get("mzTab-version", "")
    if version != _SUPPORTED_VERSION:
        errors.append(
            f"mzTab-version が '{version}' です。'{_SUPPORTED_VERSION}' のみ対応しています。"
        )

    # 必須セクション
    if "SMF" not in sections or not sections["SMF"].get("rows"):
        errors.append("SMF セクションが存在しないか行が 0 件です（特徴量行列が空）。")

    # assay 定義（値も確認）
    has_assay = any(re.match(r"^assay\[\d+\]-ms_run_ref$", k) and v for k, v in meta.items())
    if not has_assay:
        errors.append("assay[] の定義が MTD にありません。")

    # SME 末尾空欄 warning を伝播
    for section_name, sec in sections.items():
        for w in sec.get("warnings", []):
            warnings.append(f"[{section_name}] {w}")

    # ファイル全体 warning
    warnings.extend(parse_result.get("warnings", []))

    return {"ok": len(errors) == 0, "errors": errors, "warnings": warnings}


def detect_quantification_measure(
    parse_result: dict, filename: str
) -> tuple[str | None, str]:
    """ファイル名と MTD から定量種別を推定する。

    戻り値: (measure, confidence)
      measure: "peak_height" / "peak_area_above_zero" / None
      confidence: "verified" / "inferred" / "unknown" / "conflict"
    """
    meta = parse_result.get("metadata", {})
    mtd_unit = " ".join(
        v for k, v in meta.items()
        if "quantification_unit" in k or "quantification_method" in k
    )
    fname = filename

    height_in_fname = bool(_HEIGHT_PREFIX_RE.match(fname))
    area_in_fname = bool(_AREA_PREFIX_RE.match(fname))
    height_in_mtd = bool(_HEIGHT_RE.search(mtd_unit))
    area_in_mtd = bool(_AREA_RE.search(mtd_unit))

    if height_in_fname and area_in_fname:
        return None, "conflict"
    if height_in_fname:
        confidence = "verified" if height_in_mtd and not area_in_mtd else "inferred"
        return "peak_height", confidence
    if area_in_fname:
        confidence = "verified" if area_in_mtd and not height_in_mtd else "inferred"
        return "peak_area_above_zero", confidence
    if height_in_mtd and not area_in_mtd:
        return "peak_height", "inferred"
    if area_in_mtd and not height_in_mtd:
        return "peak_area_above_zero", "inferred"
    return None, "unknown"
