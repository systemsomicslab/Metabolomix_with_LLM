"""全feature定量表と統計結果の書き出し（spec §11）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

**注釈のない行を捨てない。** 同定できた feature だけを書き出すと、実際には大半を
占める未同定 feature が「存在しないもの」になり、検出率も多重比較の母数も
その表からは再現できなくなる。定量表のキーは feature_id で、annotation は
feature_id で繋ぐ**別表**にする（1 feature に複数候補が付くので、定量表へ
畳み込むと行が増えるか候補が消えるかのどちらかになる）。

**best hit を確定同定と呼ばない。** SME の候補は rank 付きで全部残し、状態は
`candidate` と書く。「一番上の候補」を確定として書き出すと、下流はそれを
同定済みとして扱う。

**欠損と ±inf を数値として書かない。** 空欄にして、理由は別列に持つ。TSV に
`inf` や `nan` を書くと、読み手のパーサが 0 や巨大値として通す。

既存の15列 volcano/pathway TSV（`analysis/export_contract.py`）はここでは触らない
——あれは別リポジトリとの契約で、ANOVA の列を足す場所ではない。
"""
from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path

import numpy as np

__all__ = ["export_features", "export_statistic"]

#: 全feature定量表の列（spec §11 の long TSV）。
_FEATURE_COLUMNS = ("dataset_id", "feature_id", "assay_id", "matrix_result_id",
                    "value", "unit", "detected", "gap_filled", "imputed",
                    "exclusion_reason")

#: annotation 側表の列。定量表とは feature_id で繋ぐ。
_ANNOTATION_COLUMNS = ("feature_id", "candidate_rank", "identification_status",
                       "name", "database_identifier", "inchikey", "adduct",
                       "charge", "confidence_measure", "confidence_value",
                       "mz", "rt_min")


def _atomic_write_text(path: Path, text: str) -> Path:
    """同じ親の一時ファイルへ書いてから置換する（`dataset_export.py` と同じ流儀）。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="\n", dir=path.parent,
                prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            tmp_path = Path(handle.name)
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
        tmp_path = None
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
    return path


def _cell(value) -> str:
    """TSV の1セル。欠損・非有限は**空欄**にする（0 や inf と読ませない）。"""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return "" if not math.isfinite(float(value)) else repr(float(value))
    text = str(value)
    return text.replace("\t", " ").replace("\n", " ")


def _table(columns, rows) -> str:
    lines = ["\t".join(columns)]
    lines.extend("\t".join(_cell(row.get(c)) for c in columns) for row in rows)
    return "\n".join(lines) + "\n"


def export_features(ds, matrices: list[dict], path: Path) -> dict:
    """全 feature × assay × matrix の long TSV と、annotation 側表を書く。

    `matrices` に複数の解析行列（recipe 違い）を渡してよい。行は
    `matrix_result_id` で区別され、単位は行列ごとの `units` をそのまま載せる
    ——未補正 height と内部標準比が同じ列に並んでも、単位を見れば区別できる。
    """
    path = Path(path)
    dataset_id = getattr(ds, "dataset_id", None)

    rows: list[dict] = []
    for matrix in matrices:
        values = np.asarray(matrix["values"], dtype=float)
        eligibility = np.asarray(matrix["eligibility_mask"], dtype=bool)
        imputed = np.asarray(matrix["imputed_mask"], dtype=bool)
        locked = np.asarray(matrix["locked_mask"], dtype=bool)
        detected = matrix.get("detected_mask")
        detected = None if detected is None else np.asarray(detected, dtype=bool)
        units = list(matrix["units"])
        missing_reasons = matrix.get("missing_reasons") or {}

        for column, feature_id in enumerate(matrix["feature_ids"]):
            for row, assay_id in enumerate(matrix["assay_ids"]):
                # 検出状態が分からないことを "false"（未検出）に倒さない。
                is_detected = (None if detected is None
                               else bool(detected[row, column]))
                reasons = []
                if not eligibility[column]:
                    reasons.append("not_eligible")
                if locked[row, column]:
                    reasons.append((missing_reasons.get(feature_id, {})
                                    .get(assay_id, {})
                                    .get("code", "denominator_invalid")))
                rows.append({
                    "dataset_id": matrix.get("dataset_id") or dataset_id,
                    "feature_id": feature_id,
                    "assay_id": assay_id,
                    "matrix_result_id": matrix["matrix_id"],
                    "value": float(values[row, column]),
                    "unit": units[column],
                    "detected": "unknown" if is_detected is None else is_detected,
                    "gap_filled": ("unknown" if is_detected is None
                                   else not is_detected),
                    "imputed": bool(imputed[row, column]),
                    "exclusion_reason": "+".join(reasons),
                })

    _atomic_write_text(path, _table(_FEATURE_COLUMNS, rows))

    annotation_path = path.with_name(f"{path.stem}.annotations.tsv")
    metadata = getattr(ds, "feature_metadata", {}) or {}
    candidates_by_feature = getattr(ds, "feature_candidates", {}) or {}
    # SML 由来のラベル。**証拠ではない**ので、読むのは ms1_annotation 行だけ。
    annotations_by_feature = getattr(ds, "feature_annotations", {}) or {}
    feature_ids = list(getattr(ds, "feature_ids", []) or [])

    annotation_rows: list[dict] = []
    for feature_id in feature_ids:
        meta = metadata.get(feature_id) or {}
        candidates = candidates_by_feature.get(feature_id) or []
        annotation = annotations_by_feature.get(feature_id) or {}
        if not candidates and annotation and not annotation.get("ambiguous"):
            # MS1 の照合だけで付いた名前。candidate と混ぜると MS/MS 裏付けが
            # あるように読まれ、unidentified と混ぜると同定できなかったと読まれる。
            annotation_rows.append({
                "feature_id": feature_id, "candidate_rank": None,
                "identification_status": "ms1_annotation",
                "name": annotation.get("name"),
                "database_identifier": annotation.get("database_identifier"),
                "inchikey": annotation.get("inchikey"),
                "adduct": annotation.get("adduct"), "charge": None,
                "confidence_measure": annotation.get("confidence_measure"),
                "confidence_value": annotation.get("confidence_value"),
                "mz": meta.get("mz"), "rt_min": meta.get("rt")})
            continue
        if not candidates:
            # 未同定でも1行残す。表から消すと「同定できたものだけの世界」になる。
            annotation_rows.append({
                "feature_id": feature_id, "candidate_rank": None,
                "identification_status": "unidentified",
                "name": meta.get("name"), "database_identifier": None,
                "inchikey": meta.get("inchikey"), "adduct": None, "charge": None,
                "confidence_measure": None, "confidence_value": None,
                "mz": meta.get("mz"), "rt_min": meta.get("rt")})
            continue
        for candidate in candidates:
            annotation_rows.append({
                "feature_id": feature_id,
                "candidate_rank": candidate.get("rank"),
                # rank 1 でも "candidate"。確定同定はここでは名乗らない。
                "identification_status": "candidate",
                "name": candidate.get("chemical_name"),
                "database_identifier": candidate.get("database_identifier"),
                "inchikey": meta.get("inchikey"),
                "adduct": candidate.get("adduct"),
                "charge": candidate.get("charge"),
                "confidence_measure": candidate.get("confidence_measure"),
                "confidence_value": candidate.get("confidence_value"),
                "mz": meta.get("mz"), "rt_min": meta.get("rt")})

    _atomic_write_text(annotation_path, _table(_ANNOTATION_COLUMNS, annotation_rows))

    return {
        "path": str(path),
        "annotation_file": annotation_path.name,
        "n_rows": len(rows),
        "n_features": len(feature_ids),
        "n_matrices": len(matrices),
        "matrix_result_ids": [m["matrix_id"] for m in matrices],
    }


# ---------- 統計結果 ----------

_WELCH_COLUMNS = ("feature_id", "status", "reason", "n_reference", "n_test",
                  "t_statistic", "p_value", "q_value", "log2fc")

_ANOVA_COLUMNS = ("feature_id", "status", "reason", "f_statistic", "df_between",
                  "df_within", "p_value", "q_value", "group_n")

_TUKEY_COLUMNS = ("feature_id", "reference_group", "test_group",
                  "mean_difference", "ci_low", "ci_high", "p_adjusted", "alpha")


def _meta_lines(result: dict, extra: dict | None = None) -> str:
    """表の先頭に付ける `# key = value` 行。どの行列・どの定義かを表に残す。"""
    meta = {
        "schema": result.get("schema"),
        "statistic_id": result.get("statistic_id"),
        "kind": result.get("kind"),
        "status": result.get("status"),
        "reason": result.get("reason"),
        "matrix_result_id": result.get("matrix_id"),
        "matrix_recipe_id": result.get("matrix_recipe_id"),
        "transform": result.get("transform"),
        "bh_population": result.get("bh_population"),
        "groups": result.get("groups"),
        **(extra or {}),
    }
    return "".join(f"# {key} = {_cell(value)}\n" for key, value in meta.items()
                   if value is not None)


def export_statistic(result: dict, path: Path) -> dict:
    """統計結果1件を TSV に書く（kind ごとに列が違う）。

    Welch と ANOVA/Tukey を同じ表にしない。Welch の15列 volcano 契約へ ANOVA の
    F や df を足すと、同じ列名が別の意味を持つ表が2種類できる。
    """
    path = Path(path)
    kind = result.get("kind")
    features = result.get("features") or []

    if kind == "anova_tukey":
        rows = [{**f, "group_n": ",".join(f"{k}={v}" for k, v in
                                          sorted((f.get("group_n") or {}).items()))}
                for f in features]
        text = _meta_lines(result, {"alpha": result.get("alpha")}) \
            + _table(_ANOVA_COLUMNS, rows)
        _atomic_write_text(path, text)

        tukey_path = path.with_name(f"{path.stem}.tukey.tsv")
        pairs = [{"feature_id": f["feature_id"], **pair, "alpha": result.get("alpha")}
                 for f in features for pair in (f.get("tukey") or [])]
        tukey_text = _meta_lines(result, {
            # feature 横断の FDR ではないことを表自身に書く（spec §10）。
            "correction_scope": result.get("tukey_correction_scope",
                                           "within_feature_group_pairs")}) \
            + _table(_TUKEY_COLUMNS, pairs)
        _atomic_write_text(tukey_path, tukey_text)
        extra = {"tukey_file": tukey_path.name, "n_tukey_pairs": len(pairs)}
    else:
        text = _meta_lines(result, {
            "effect_size_definition": result.get("effect_size_definition"),
            "reference_group": result.get("reference_group"),
            "test_group": result.get("test_group")}) \
            + _table(_WELCH_COLUMNS, features)
        _atomic_write_text(path, text)
        extra = {"effect_size_definition": result.get("effect_size_definition")}

    return {
        "path": str(path),
        "statistic_id": result.get("statistic_id"),
        "kind": kind,
        "status": result.get("status"),
        "reason": result.get("reason"),
        "n_features": len(features),
        "n_tested": sum(1 for f in features if f.get("status") == "tested"),
        **extra,
    }
