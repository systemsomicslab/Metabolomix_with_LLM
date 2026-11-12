"""sample-manifest v2 の後ろに足した因子の列（Genotype・Treatment など）を通す（裁定 A1）。"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from metabolomix.analysis.sample_manifest import (
    FIELDS_V2,
    OPTIONAL_FIELDS,
    apply_metadata,
    parse_manifest,
)
from metabolomix.analysis.stats.spec_v2 import normalize_spec_v2, run_spec_v2
from metabolomix.core.atomic_io import DomainError
from tests.pipeline_fixtures import make_dataset

_V2 = "\t".join(FIELDS_V2)


def _write(tmp_path: Path, header: str, lines: list[str]) -> Path:
    sheet = tmp_path / "sample-manifest.tsv"
    sheet.write_text("# schema = sample-manifest.v2\n" + header + "\n" + "\n".join(lines) + "\n",
                     encoding="utf-8")
    return sheet


def _line(i: int, *extra: str, bio: str | None = None) -> str:
    cells = [f"s{i}", f"S{i}.wiff", "sample", "A", "B1", str(i + 1), "", "true",
             bio if bio is not None else f"bio{i}", *extra]
    return "\t".join(cells)


def _parse(tmp_path: Path, header: str, lines: list[str], n: int = 2) -> list[dict]:
    root = tmp_path / "source"
    root.mkdir(exist_ok=True)
    names = [f"S{i}.wiff" for i in range(n)]
    for name in names:
        (root / name).write_bytes(b"fake")
    sheet = _write(tmp_path, header, lines)
    return parse_manifest(sheet, source_root=root, expected_sources=names)


def test_manifest_accepts_trailing_factor_columns(tmp_path):
    rows = _parse(tmp_path, _V2 + "\tGenotype\tTreatment",
                  [_line(0, "WT", "Veh"), _line(1, "KO", "")])
    assert [r["Genotype"] for r in rows] == ["WT", "KO"]
    assert [r["Treatment"] for r in rows] == ["Veh", None]
    assert rows[0]["provenance"]["Genotype"]["value"] == "WT"
    assert rows[1]["provenance"]["Treatment"]["value"] is None


def test_manifest_accepts_subject_id_then_factor_columns(tmp_path):
    rows = _parse(tmp_path, _V2 + "\t" + OPTIONAL_FIELDS[0] + "\tGenotype",
                  [_line(0, "m1", "WT"), _line(1, "m2", "KO")])
    assert rows[0]["subject_id"] == "m1" and rows[1]["Genotype"] == "KO"


@pytest.mark.parametrize("extra", ["group", "biological_sample_id", "provenance", "conflicts", "", "   ", "Geno\tGeno"])
def test_manifest_rejects_factor_column_clashing_with_known_column(tmp_path, extra):
    # 既知の列と重なる・空の列名・足した列名どうしの重複は拒む。
    header = _V2 + "\t" + extra
    with pytest.raises(DomainError) as exc:
        _parse(tmp_path, header, [_line(0, "x"), _line(1, "y")])
    assert exc.value.code == "SAMPLE_MANIFEST_INVALID"


def test_v1_sheet_with_trailing_extra_column_is_still_rejected(tmp_path):
    from metabolomix.analysis.sample_manifest import FIELDS
    root = tmp_path / "source"
    root.mkdir()
    (root / "S0.wiff").write_bytes(b"fake")
    sheet = tmp_path / "v1.tsv"
    sheet.write_text("# schema = sample-manifest.v1\n" + "\t".join(FIELDS) + "\tGenotype\n"
                     "s0\tS0.wiff\tsample\tA\tB1\t1\t\ttrue\tWT\n", encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert exc.value.code == "SAMPLE_MANIFEST_INVALID"


def test_manifest_rejects_known_column_among_extra_columns(tmp_path):
    with pytest.raises(DomainError) as exc:
        _parse(tmp_path, _V2 + "\tGenotype\tsubject_id", [_line(0, "WT", "m"), _line(1, "KO", "n")])
    assert exc.value.code == "SAMPLE_MANIFEST_INVALID"


def test_manifest_rejects_short_row_with_factor_columns(tmp_path):
    with pytest.raises(DomainError) as exc:
        _parse(tmp_path, _V2 + "\tGenotype", [_line(0, "WT"), _line(1)])
    assert exc.value.code == "SAMPLE_MANIFEST_INVALID"


def test_old_manifest_without_factor_columns_still_parses(tmp_path):
    rows = _parse(tmp_path, _V2, [_line(0), _line(1)])
    assert [r["group"] for r in rows] == ["A", "A"]
    assert set(rows[0]) == set(FIELDS_V2) | {"provenance", "conflicts"}


def _factor_ds(tmp_path: Path, n_per=4):
    """Genotype × Treatment の 2×2 を持つ sample-manifest を当てたデータセット。"""
    cells = [(g, t) for g in ("WT", "KO") for t in ("Veh", "Drug")]
    n = len(cells) * n_per
    root = tmp_path / "source"
    root.mkdir()
    lines = []
    for i in range(n):
        g, t = cells[i % len(cells)]
        (root / f"S{i}.wiff").write_bytes(b"fake")
        lines.append("\t".join([f"s{i}", f"S{i}.wiff", "sample", f"{g}_{t}", "B1", str(i + 1), "",
                                "true", f"bio{i}", g, t]))
    sheet = _write(tmp_path, _V2 + "\tGenotype\tTreatment", lines)
    rows = parse_manifest(sheet, source_root=root, expected_sources=[f"S{i}.wiff" for i in range(n)])
    rng = np.random.default_rng(5)
    fm = np.exp2(rng.normal(20, 1, (8, n)))
    for j, r in enumerate(rows):
        if r["Genotype"] == "KO" and r["Treatment"] == "Drug":
            fm[:3, j] *= 4
    names = [r["sample_id"] for r in rows]
    ds = SimpleNamespace(
        feature_matrix=fm, detected_mask=np.ones_like(fm, dtype=bool),
        sample_names=names, sample_ids=names, sample_metadata_rows=rows,
        feature_ids=[f"f{k}" for k in range(8)],
        feature_metadata={f"f{k}": {"ontology": "PC"} for k in range(8)},
        feature_candidates={})
    return ds


def test_design_factors_resolve_from_manifest_factor_columns(tmp_path):
    ds = _factor_ds(tmp_path)
    spec = {"schema": "analysis-spec.v2",
            "design": {"factors": ["Genotype", "Treatment"]},
            "preprocess": {"drift": {"mode": "off"},
                           "filters": {"blank_ratio": None, "d_ratio": None, "qc_rsd": None}},
            "comparisons": [
                {"id": "drug_in_ko", "factor": "Treatment", "reference": "Veh", "test": "Drug",
                 "at": {"Genotype": "KO"}},
                {"id": "gxt", "kind": "interaction", "factors": ["Genotype", "Treatment"],
                 "levels": {"Genotype": ["WT", "KO"], "Treatment": ["Veh", "Drug"]}}],
            "analyses": [{"analysis_id": "a1", "kind": "compare"}]}
    result = run_spec_v2(normalize_spec_v2(spec), ds, assay_kind="lipid")
    res = result["analyses"][0]["result"]
    for contrast in ("drug_in_ko", "gxt"):
        mine = [r for r in res["rows"] if r["contrast_id"] == contrast]
        assert len(mine) == 8
        assert all(r["status"] == "tested" for r in mine)
    assert res["design"]["factors"] == ["Genotype", "Treatment"]
    assert res["omnibus"]


def test_factor_column_change_invalidates_results():
    ds = make_dataset()

    def rows(treat):
        out = []
        for i, name in enumerate(ds.sample_names):
            row = {"sample_id": name, "source_file": f"{name}.wiff", "role": "sample",
                   "group": "A", "batch": "B1", "injection_order": i + 1, "qc_pool": None,
                   "include": True, "biological_sample_id": f"bio{i}",
                   "Treatment": treat(i), "conflicts": []}
            row["provenance"] = {f: {"value": row[f], "source": "user_manifest",
                                     "confidence": "confirmed"}
                                 for f in (*FIELDS_V2, "Treatment")}
            out.append(row)
        return out

    first = apply_metadata(ds, rows(lambda i: "Veh" if i % 2 else "Drug"))
    assert "Treatment" in first["changed_fields"]
    ds.stat_matrices = {"mx": {}}
    ds.last_differential = {"x": 1}
    same = apply_metadata(ds, rows(lambda i: "Veh" if i % 2 else "Drug"))
    assert same["changed_fields"] == []
    second = apply_metadata(ds, rows(lambda i: "Drug" if i % 2 else "Veh"))
    assert second["changed_fields"] == ["Treatment"]
    assert second["metadata_revision"] == first["metadata_revision"] + 1
    assert ds.last_differential is None
