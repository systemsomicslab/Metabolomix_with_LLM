"""sample-manifest.v2 — standard役割と生物試料ID（spec §7）。

v1が縛るのは「注入の台帳として矛盾がないか」だった。v2が追加で縛るのは
**統計の母集団が何か**の2点:

1. `role=standard`（標準品注入）は測定の道具であって生物試料ではない。統計群・
   pooled QC・正規化参照へ自動投入しない。roleを増やす以上、v1のシートへ
   `standard`が書かれたら受け付けない——v1の読み手（既存経路）はこの値を
   知らないため、素通りさせるより読込時に止めるほうが安全。
2. `biological_sample_id` は反復注入を束ねる生物単位。注入数をnと読むと、同じ
   個体の2注入が「n=2」に化けて有意性が水増しされる。初期版は自動平均も
   混合モデルもせず、`REPEATED_MEASURES_UNSUPPORTED` で止める。

列選択（include/role/group）はここの共通関数へ集約する。同じ規則を統計・QC・
標準注入解決が各々書き直すと、片方だけ `include=false` を見落とす。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from metabolomix.analysis.result_state import invalidate_results
from metabolomix.analysis.sample_manifest import (
    FIELDS,
    FIELDS_V2,
    apply_metadata,
    parse_manifest,
    select_statistical_samples,
    validate_independent_samples,
    validate_standard_assays,
)
from metabolomix.core.atomic_io import DomainError
from tests.pipeline_fixtures import make_dataset, metadata_rows

_HEADER_V1 = "\t".join(FIELDS)
_HEADER_V2 = "\t".join(FIELDS_V2)


def _sheet(tmp_path: Path, body: str, *, schema: str = "sample-manifest.v2",
           header: str | None = None) -> Path:
    sheet = tmp_path / "sample-manifest.tsv"
    head = _HEADER_V2 if header is None else header
    sheet.write_text(f"# schema = {schema}\n{head}\n{body}", encoding="utf-8")
    return sheet


def _line(sample_id="s0", source_file="S0.wiff", role="sample", group="A",
          batch="B1", injection_order="1", qc_pool="", include="true",
          biological_sample_id="bio0") -> str:
    return "\t".join([sample_id, source_file, role, group, batch,
                      injection_order, qc_pool, include,
                      biological_sample_id]) + "\n"


def _touch(root: Path, *names: str) -> None:
    for name in names:
        (root / name).write_bytes(b"fake")


def _rows(*specs: dict) -> list[dict]:
    """選択関数・apply_metadataへ渡す合成v2行（provenanceは最小形）。"""
    rows = []
    for index, spec in enumerate(specs):
        row = {"sample_id": f"s{index}", "source_file": f"S{index}.wiff",
               "role": "sample", "group": "A", "batch": "B1",
               "injection_order": index + 1, "qc_pool": None, "include": True,
               "biological_sample_id": f"bio{index}"}
        row.update(spec)
        row["provenance"] = {
            field: {"value": row[field], "source": "user_manifest",
                    "confidence": "confirmed"}
            for field in FIELDS_V2
        }
        row["conflicts"] = []
        rows.append(row)
    return rows


# ---------- 9列とschema dispatch ----------

def test_v2_sheet_parses_nine_columns(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, _line())
    rows = parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert rows[0]["biological_sample_id"] == "bio0"
    assert rows[0]["provenance"]["biological_sample_id"]["source"] == "user_manifest"


def test_v2_header_with_v1_columns_is_rejected(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, "\t".join([
        "s0", "S0.wiff", "sample", "A", "B1", "1", "", "true"]) + "\n",
        header=_HEADER_V1)
    with pytest.raises(DomainError) as caught:
        parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert caught.value.code == "SAMPLE_MANIFEST_INVALID"


def test_v1_header_with_v2_columns_is_rejected(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, _line(), schema="sample-manifest.v1")
    with pytest.raises(DomainError) as caught:
        parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert caught.value.code == "SAMPLE_MANIFEST_INVALID"


def test_unknown_schema_version_is_rejected(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, _line(), schema="sample-manifest.v3")
    with pytest.raises(DomainError) as caught:
        parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert caught.value.code == "SAMPLE_MANIFEST_INVALID"


def test_v1_sheet_still_parses_eight_columns(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, "\t".join([
        "s0", "S0.wiff", "sample", "A", "B1", "1", "", "true"]) + "\n",
        schema="sample-manifest.v1", header=_HEADER_V1)
    rows = parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert "biological_sample_id" not in rows[0]


# ---------- standard role ----------

def test_standard_role_is_accepted_in_v2(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, _line(role="standard", group="",
                               biological_sample_id=""))
    rows = parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert rows[0]["role"] == "standard"
    assert rows[0]["biological_sample_id"] is None


def test_standard_role_is_rejected_in_v1(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff")
    sheet = _sheet(root, "\t".join([
        "s0", "S0.wiff", "standard", "", "B1", "1", "", "true"]) + "\n",
        schema="sample-manifest.v1", header=_HEADER_V1)
    with pytest.raises(DomainError) as caught:
        parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert caught.value.code == "SAMPLE_MANIFEST_INVALID"


def test_standard_is_not_a_statistical_sample():
    rows = _rows({}, {"role": "standard", "biological_sample_id": None},
                 {"role": "qc", "biological_sample_id": None})
    selected = select_statistical_samples(rows, ["A"])
    assert [r["sample_id"] for r in selected] == ["s0"]


# ---------- include / group の選択 ----------

def test_excluded_rows_leave_the_statistical_population():
    rows = _rows({}, {"include": False}, {"group": "B"})
    assert [r["sample_id"] for r in select_statistical_samples(rows, ["A"])] == ["s0"]
    assert [r["sample_id"] for r in select_statistical_samples(rows, ["A", "B"])] \
        == ["s0", "s2"]


def test_rows_without_group_are_not_silently_grouped():
    rows = _rows({"group": None})
    assert select_statistical_samples(rows, ["A"]) == []


# ---------- 生物試料ID ----------

def test_duplicate_biological_id_is_not_independent_n():
    rows = [{"sample_id": x, "biological_sample_id": "bio1", "group": "A",
             "role": "sample", "include": True} for x in ("inj1", "inj2")]
    with pytest.raises(DomainError) as caught:
        validate_independent_samples(rows, ["A"])
    assert caught.value.code == "REPEATED_MEASURES_UNSUPPORTED"


def test_duplicate_biological_id_across_groups_is_also_repeated_measures():
    rows = _rows({"biological_sample_id": "bio1"},
                 {"group": "B", "biological_sample_id": "bio1"})
    with pytest.raises(DomainError) as caught:
        validate_independent_samples(rows, ["A", "B"])
    assert caught.value.code == "REPEATED_MEASURES_UNSUPPORTED"


def test_duplicate_biological_id_outside_the_population_is_allowed():
    """除外行・standard・対象外群の重複は統計のnに入らないので止めない。"""
    rows = _rows({"biological_sample_id": "bio1"},
                 {"include": False, "biological_sample_id": "bio1"},
                 {"role": "standard", "biological_sample_id": "bio1"},
                 {"group": "Z", "biological_sample_id": "bio1"})
    validate_independent_samples(rows, ["A"])


def test_missing_biological_id_blocks_the_test():
    rows = _rows({"biological_sample_id": None}, {})
    with pytest.raises(DomainError) as caught:
        validate_independent_samples(rows, ["A"])
    assert caught.value.code == "BIOLOGICAL_SAMPLE_ID_REQUIRED"
    assert caught.value.details["sample_ids"] == ["s0"]


def test_v1_rows_without_the_column_block_the_test():
    """v1シートは列自体が無い。無い＝空欄と同じく「独立性未確認」で止める。"""
    ds = make_dataset()
    rows = metadata_rows(ds, n_qc=0)
    for row in rows:
        row["group"] = "A"
    with pytest.raises(DomainError) as caught:
        validate_independent_samples(rows, ["A"])
    assert caught.value.code == "BIOLOGICAL_SAMPLE_ID_REQUIRED"


# ---------- standard_assays ----------

def test_standard_assays_must_point_at_standard_rows():
    rows = _rows({}, {"role": "standard", "sample_id": "std-1",
                      "biological_sample_id": None})
    validate_standard_assays(rows, {"IS1": ["std-1"]})
    with pytest.raises(DomainError) as caught:
        validate_standard_assays(rows, {"IS1": ["s0"]})
    assert caught.value.code == "STANDARD_ASSAY_INVALID"


def test_standard_assays_reject_excluded_injections():
    rows = _rows({"role": "standard", "sample_id": "std-1", "include": False,
                  "biological_sample_id": None})
    with pytest.raises(DomainError) as caught:
        validate_standard_assays(rows, {"IS1": ["std-1"]})
    assert caught.value.code == "STANDARD_ASSAY_INVALID"


def test_standard_assays_reject_unknown_sample_id():
    rows = _rows({"role": "standard", "sample_id": "std-1",
                  "biological_sample_id": None})
    with pytest.raises(DomainError) as caught:
        validate_standard_assays(rows, {"IS1": ["std-9"]})
    assert caught.value.code == "STANDARD_ASSAY_INVALID"


# ---------- rawとの対応（v2でもv1と同じ厳しさ） ----------

def test_missing_source_row_is_rejected_in_v2(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff", "S1.wiff")
    sheet = _sheet(root, _line())
    with pytest.raises(DomainError) as caught:
        parse_manifest(sheet, source_root=root,
                       expected_sources=["S0.wiff", "S1.wiff"])
    assert caught.value.code == "SAMPLE_MANIFEST_INVALID"


def test_extra_source_row_is_rejected_in_v2(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff", "S1.wiff")
    sheet = _sheet(root, _line() + _line(sample_id="s1", source_file="S1.wiff",
                                         injection_order="2",
                                         biological_sample_id="bio1"))
    with pytest.raises(DomainError) as caught:
        parse_manifest(sheet, source_root=root, expected_sources=["S0.wiff"])
    assert caught.value.code == "SAMPLE_MANIFEST_INVALID"


def test_duplicate_biological_id_is_allowed_in_the_sheet(tmp_path):
    """シート自体は反復注入を書けてよい（止めるのは検定の直前）。"""
    root = tmp_path / "source"
    root.mkdir()
    _touch(root, "S0.wiff", "S1.wiff")
    sheet = _sheet(root, _line() + _line(sample_id="s1", source_file="S1.wiff",
                                         injection_order="2",
                                         biological_sample_id="bio0"))
    rows = parse_manifest(sheet, source_root=root,
                          expected_sources=["S0.wiff", "S1.wiff"])
    assert [r["biological_sample_id"] for r in rows] == ["bio0", "bio0"]


# ---------- 指紋と無効化 ----------

def test_biological_id_change_is_detected_and_keeps_preprocessing():
    ds = make_dataset()
    rows = _rows(*[{"sample_id": n, "source_file": f"{n}.wiff",
                    "injection_order": i + 1, "biological_sample_id": f"bio{i}"}
                   for i, n in enumerate(ds.sample_names)])
    apply_metadata(ds, rows)
    ds.pp_matrix = object()
    ds.last_differential = {"x": 1}

    changed = _rows(*[{"sample_id": n, "source_file": f"{n}.wiff",
                       "injection_order": i + 1,
                       "biological_sample_id": f"donor{i}"}
                      for i, n in enumerate(ds.sample_names)])
    report = apply_metadata(ds, changed)

    assert "biological_sample_id" in report["changed_fields"]
    assert ds.last_differential is None
    assert ds.pp_matrix is not None


def test_standard_assay_change_invalidates_the_preprocessed_matrix():
    ds = make_dataset()
    ds.pp_matrix = object()
    ds.last_pca = {"x": 1}
    invalidate_results(ds, {"standard_assays"})
    assert ds.pp_matrix is None
    assert ds.last_pca is None
