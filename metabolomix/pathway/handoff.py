# metabolomix/pathway/handoff.py
"""ARF / DatasetState → 上流 write_mztab_m と同形の最小 mzTab-M（純ロジック）。

群・除外・前処理の決定権はこちらに置く（spec §2）。上流は mzTab-M の
`study_variable` を群に、`opt_global_Ontology` をサブクラスに読む。
session_state / mcp_core を import しない。
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_UNANNOTATED = {"", "unknown", "others", "null"}
_SPOT_RE = re.compile(r"^Spot_(\d+)_(.+)$")
_MIN_GROUP = 2


class HandoffError(Exception):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details or {}


@dataclass
class HandoffTable:
    sample_names: list[str]
    groups: list[str]
    names: list[str]
    ontologies: list[str]
    matrix: np.ndarray            # (n_features, n_samples)。NaN = 欠損
    report: dict = field(default_factory=dict)


def _check_group_sizes(groups: list[str], label_a: str, label_b: str) -> None:
    n_a, n_b = groups.count(label_a), groups.count(label_b)
    if n_a < _MIN_GROUP or n_b < _MIN_GROUP:
        raise HandoffError("PATHWAY_BAD_REQUEST",
                           f"各群 2 例以上が要ります（{label_a}={n_a}, {label_b}={n_b}）。",
                           {"n_a": n_a, "n_b": n_b})


def _annotated(name, ontology) -> bool:
    return (str(name or "").strip().lower() not in _UNANNOTATED
            and str(ontology or "").strip().lower() not in _UNANNOTATED)


def _dedup_by_name(names, ontologies, rows):
    """同名の行は平均強度が最大の 1 行だけ残す（sum の二重計上を避ける。spec §4）。"""
    best: dict[str, int] = {}
    for i, name in enumerate(names):
        mean = float(np.nanmean(rows[i])) if np.isfinite(rows[i]).any() else -math.inf
        j = best.get(name)
        if j is None:
            best[name] = i
            continue
        mean_j = float(np.nanmean(rows[j])) if np.isfinite(rows[j]).any() else -math.inf
        if mean > mean_j:
            best[name] = i
    keep = sorted(best.values())
    dup_names = sorted({n for n in names if names.count(n) > 1})
    return keep, len(names) - len(keep), dup_names[:5]


def _finish(sample_names, groups, names, ontologies, rows, report) -> HandoffTable:
    keep, n_dup, examples = _dedup_by_name(names, ontologies, rows)
    report["dropped_duplicate_names"] = n_dup
    report["duplicate_name_examples"] = examples
    matrix = (np.vstack([rows[i] for i in keep]) if keep
              else np.empty((0, len(sample_names))))
    if not keep:
        raise HandoffError("PATHWAY_BAD_REQUEST",
                           "サブクラス注釈の付いた特徴が 1 件もありません。", dict(report))
    return HandoffTable(sample_names=sample_names, groups=groups,
                        names=[names[i] for i in keep],
                        ontologies=[ontologies[i] for i in keep],
                        matrix=matrix, report=report)


def build_from_arf(*, matrix, sample_names, feature_names, sample_meta,
                   excluded_samples, excluded_spots, catalog, props,
                   group_a: str, group_b: str, compared: dict) -> HandoffTable:
    matrix = np.asarray(matrix, dtype=float)          # (n_samples, n_features)
    prop_list = [props] if isinstance(props, str) else list(props or [])
    prop = prop_list[0] if prop_list else None
    a_set, b_set = set(compared["group_a"]), set(compared["group_b"])

    col_idx, sample_out, groups = [], [], []
    for i, name in enumerate(sample_names):
        entry = sample_meta.get(name) or {}
        if entry.get("role", "sample") != "sample" or name in excluded_samples:
            continue
        if name in a_set:
            label = group_a
        elif name in b_set:
            label = group_b
        else:
            label = entry.get("group")
            # 比較外の試料の群名が比較の指定名と同じなら書かない（比較群に紛れ込む）
            if not label or str(label) in (group_a, group_b):
                continue
        col_idx.append(i)
        sample_out.append(name)
        groups.append(str(label))
    _check_group_sizes(groups, group_a, group_b)

    names, ontologies, rows = [], [], []
    report = {"source": "arf", "prop": prop, "dropped_unannotated": 0,
              "dropped_excluded_spots": 0}
    for j, feat in enumerate(feature_names):
        m = _SPOT_RE.match(str(feat))
        if not m or (prop is not None and m.group(2) != prop):
            continue
        spot_id = int(m.group(1))
        if spot_id in excluded_spots:
            report["dropped_excluded_spots"] += 1
            continue
        spot = catalog.get(spot_id) or {}
        name, onto = spot.get("Name"), spot.get("Ontology")
        if not _annotated(name, onto):
            report["dropped_unannotated"] += 1
            continue
        names.append(str(name).strip())
        ontologies.append(str(onto).strip())
        rows.append(matrix[col_idx, j])
    return _finish(sample_out, groups, names, ontologies, rows, report)


def build_from_dataset(*, pp_matrix, pp_sample_names, pp_feature_names, roles,
                       sample_names, sample_ids, metadata_rows, annotations,
                       group_a, group_b, group_a_label: str,
                       group_b_label: str) -> HandoffTable:
    matrix = np.asarray(pp_matrix, dtype=float)       # (n_samples, n_features)
    a_set, b_set = set(group_a), set(group_b)
    meta_group: dict[str, str] = {}
    if metadata_rows:
        id_of = (dict(zip(sample_names, sample_ids)) if sample_ids
                 else {n: n for n in sample_names})
        by_id = {str(r.get("sample_id")): r for r in metadata_rows}
        for name in sample_names:
            row = by_id.get(str(id_of.get(name, name)))
            if (row and row.get("include", True) and row.get("role", "sample") == "sample"
                    and row.get("group")):
                meta_group[name] = str(row["group"])

    col_idx, sample_out, groups = [], [], []
    for i, name in enumerate(pp_sample_names):
        if roles.get(name, "sample") != "sample":
            continue
        if name in a_set:
            label = group_a_label
        elif name in b_set:
            label = group_b_label
        elif name in meta_group and meta_group[name] not in (group_a_label, group_b_label):
            label = meta_group[name]
        else:
            continue
        col_idx.append(i)
        sample_out.append(name)
        groups.append(label)
    _check_group_sizes(groups, group_a_label, group_b_label)

    names, ontologies, rows = [], [], []
    report = {"source": "dataset", "dropped_unannotated": 0}
    for j, fid in enumerate(pp_feature_names):
        ann = annotations.get(str(fid)) or {}
        if ann.get("ambiguous") or not _annotated(ann.get("name"), ann.get("ontology")):
            report["dropped_unannotated"] += 1
            continue
        names.append(str(ann["name"]).strip())
        ontologies.append(str(ann["ontology"]).strip())
        rows.append(matrix[col_idx, j])
    return _finish(sample_out, groups, names, ontologies, rows, report)


def _cell(v) -> str:
    if v is None or not math.isfinite(float(v)):
        return "null"
    return repr(float(v))


def _clean(s: str) -> str:
    return re.sub(r"[\t\r\n]+", " ", str(s))


def write_mztab(table: HandoffTable, path: Path, *, normalization: str) -> None:
    """上流 R `write_mztab_m` と同じ列構成で書く（R/mztab_m_import.R）。"""
    levels: list[str] = []
    for g in table.groups:
        if g not in levels:
            levels.append(g)
    lines = [["MTD", "mzTab-version", "2.0.0-M"],
             ["MTD", "mzTab-ID", "metabolomix-pathway"],
             ["MTD", "software[1]", "[,,metabolomix,]"],
             ["MTD", "small_molecule-quantification_unit", "[,,abundance,]"],
             ["MTD", "comment[1]", f"Normalization: {_clean(normalization)}"]]
    for i, name in enumerate(table.sample_names, start=1):
        lines.append(["MTD", f"assay[{i}]", _clean(name)])
    for k, level in enumerate(levels, start=1):
        refs = "|".join(f"assay[{i}]" for i, g in enumerate(table.groups, start=1) if g == level)
        lines += [["MTD", f"study_variable[{k}]", _clean(level)],
                  ["MTD", f"study_variable[{k}]-assay_refs", refs],
                  ["MTD", f"study_variable[{k}]-description", _clean(level)]]
    n = len(table.sample_names)
    lines.append(["SMH", "SML_ID", "SMF_ID_REFS", "database_identifier",
                  "chemical_formula", "smiles", "inchi", "chemical_name", "uri",
                  "theoretical_neutral_mass", "adduct_ions", "reliability",
                  "best_id_confidence_measure", "best_id_confidence_value",
                  *[f"abundance_assay[{i}]" for i in range(1, n + 1)],
                  "opt_global_Ontology"])
    for r, (name, onto) in enumerate(zip(table.names, table.ontologies), start=1):
        lines.append(["SML", str(r), str(r), "null", "null", "null", "null",
                      _clean(name), "null", "null", "null", "null", "null", "null",
                      *[_cell(v) for v in table.matrix[r - 1]], _clean(onto)])
    path.write_text("\n".join("\t".join(l) for l in lines) + "\n", encoding="utf-8")
