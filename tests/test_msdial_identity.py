import pytest

from metabolomix.arf2.spot_fields import spot_fields
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.dataset.msdial_identity import (
    IdentityCache, annotation_source, identity_columns, is_sme_like, sme_candidate,
    sml_annotation)
from metabolomix.mztab.loading import load_dataset_state
from tests.curation_fixtures import arf2_spot_raw, match_result
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.reanalysis_fixtures import POS_FEATURES, write_mztab

TAG_XML = """<?xml version="1.0" encoding="utf-8"?>
<PeakSpotTags>
  <Definitions><Tag><Id>1</Id><Label>Confirmed</Label></Tag></Definitions>
  <Peaks><Peak Id="0"><Tag>1</Tag></Peak></Peaks>
</PeakSpotTags>
"""


def _samples(n=3):
    return [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(n)]


@pytest.mark.parametrize("name,overrides,msms,expected", [
    ("PC 34:1", {}, True, True),
    ("PC 34:1", {}, False, False),                 # (a) MS/MS 未取得
    ("PC 34:1", {26: 4 | 64}, True, False),        # (b) 手で直した代表
    ("PC 34:1", {26: 16}, True, False),            # (c) Text DB
    ("", {}, True, False),                         # (d) 名前が空
    ("no MS2: PC 34:1", {}, True, False),          # (f)
    # `w/o MS2:` は上流の ShouldWriteSmeLine が除外していない（`no MS2` だけ）。MS/MS が付いて
    # いれば SME 相当（mzTab-M の書き出しと同じ）。付いていなければ (a) で落ちる。
    ("w/o MS2: PC 34:1", {}, True, True),
    ("w/o MS2: PC 34:1", {}, False, False),
])
def test_sme_rule_ported_from_mztab_export(name, overrides, msms, expected):
    raw = arf2_spot_raw(spot_id=0, name=name, matches=[match_result(overrides)])
    assert is_sme_like(spot_fields(raw), msms_assigned=msms) is expected


def test_blank_filtered_spot_is_not_sme():
    raw = arf2_spot_raw(spot_id=0, matches=[match_result()])
    raw[58] = True
    assert is_sme_like(spot_fields(raw), msms_assigned=True) is False


def test_candidate_and_annotation_keys_match_mztab_path(tmp_path):
    features = [dict(f, matched=True) for f in POS_FEATURES]
    ds_mz = load_dataset_state(mztab_path=write_mztab(tmp_path / "x.mzTab", features=features))
    fid = next(iter(ds_mz.feature_candidates))
    raw = arf2_spot_raw(spot_id=0, matches=[match_result()])
    assert set(sme_candidate(spot_fields(raw), "0")) == set(ds_mz.feature_candidates[fid][0])
    fid_a = next(iter(ds_mz.feature_annotations))
    assert set(sml_annotation(spot_fields(raw))) == set(ds_mz.feature_annotations[fid_a])


def test_unknown_has_no_annotation():
    assert sml_annotation(spot_fields(arf2_spot_raw(spot_id=0, name="Unknown"))) is None


def test_sme_candidate_name_is_representative_last_segment():
    raw = arf2_spot_raw(spot_id=0, name="low score: PC 34:1",
                        matches=[match_result({0: "PC 34:1|PC 16:0_18:1"})])
    fields = spot_fields(raw)
    assert sme_candidate(fields, "0")["chemical_name"] == "PC 16:0_18:1"
    assert sml_annotation(fields)["name"] == "low score: PC 34:1"
    piped = spot_fields(arf2_spot_raw(spot_id=0, name="PC 34:1|PC 16:0_18:1",
                                      matches=[match_result()]))
    assert sml_annotation(piped)["name"] == "PC 16:0_18:1"


def test_sme_candidate_without_representative_has_no_name():
    fields = spot_fields(arf2_spot_raw(spot_id=0, name="PC 34:1"))
    assert fields["representative"] is None
    assert is_sme_like(fields, msms_assigned=True) is True
    candidate = sme_candidate(fields, "0")
    assert candidate["chemical_name"] is None and candidate["confidence_value"] is None


def test_unknown_spot_with_msms_and_representative_is_not_sme():
    fields = spot_fields(arf2_spot_raw(spot_id=0, name="Unknown", matches=[match_result()]))
    assert is_sme_like(fields, msms_assigned=True) is False


@pytest.mark.parametrize("name,expected", [
    ("null", False), ("Unknown", False), ("", False), (None, False),
    ("unknown", True)])   # 上流は大文字小文字を区別する（spot.Name != "Unknown"）
def test_sml_name_rule_matches_upstream_case_sensitivity(name, expected):
    fields = spot_fields(arf2_spot_raw(spot_id=0, name="x"))
    fields["name"] = name
    assert (sml_annotation(fields) is not None) is expected


def test_spectrum_matched_is_none_without_representative_in_both_slots():
    fields = spot_fields(arf2_spot_raw(spot_id=0, name="PC 34:1"))
    assert sml_annotation(fields)["spectrum_matched"] is None


def test_unknown_source_bits_are_other_not_none():
    fields = spot_fields(arf2_spot_raw(spot_id=0, matches=[match_result({26: 2})]))
    assert annotation_source(fields) == "other"


def test_generated_lipid_beats_msp():
    fields = spot_fields(arf2_spot_raw(spot_id=0, matches=[match_result({26: 4 | 32})]))
    assert annotation_source(fields) == "generated_lipid"


@pytest.mark.parametrize("source,expected", [
    (4, "msp"), (16, "text_db"), (32, "generated_lipid"), (64, "manual"),
    (4 | 64, "manual"), (4 | 16, "text_db"), (0, "none")])
def test_annotation_source_bits_and_precedence(source, expected):
    raw = arf2_spot_raw(spot_id=0, matches=[match_result({26: source})])
    fields = spot_fields(raw)
    if source == 0:
        fields["representative"] = None
    assert annotation_source(fields) == expected


def test_identity_is_computed_once_per_name_ontology_msms():
    cache = IdentityCache()
    fields = spot_fields(arf2_spot_raw(spot_id=0))
    other = spot_fields(arf2_spot_raw(spot_id=1))
    cache.columns(fields, msms_assigned=True)
    cache.columns(other, msms_assigned=True)
    assert cache.computed == 1
    cache.columns(other, msms_assigned=False)
    assert cache.computed == 2
    assert set(identity_columns(fields, msms_assigned=True, tables=cache.tables)) == {
        "normalized", "refmet", "lipid_maps_category", "msi_level"}


def test_adapter_fills_identity_tags_and_flags(tmp_path):
    write_msdial_folder(tmp_path, spots=[
        {"name": "PC 34:1", "ontology": "PC", "msms": [True, False, False], "matches": [match_result()]},
        {"name": "Unknown"}], samples=_samples())
    selection = select_alignment(tmp_path)
    (tmp_path / f"{selection.arf2.stem}_tags.xml").write_text(TAG_XML, encoding="utf-8")
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_candidates["0"][0]["chemical_name"] == "PC 34:1"
    assert ds.feature_candidates["0"][0]["confidence_value"] == pytest.approx(3.2)  # match_result の Key2（total_score）
    assert "1" not in ds.feature_candidates and "1" not in ds.feature_annotations
    assert ds.feature_annotations["0"]["name"] == "PC 34:1"
    meta = ds.feature_metadata["0"]
    assert meta["identity"]["msi_level"] in (2, 3)
    assert meta["annotation_source"] == "msp" and meta["spectrum_matched"] is True
    assert meta["inchikey_source"] == "arf2" and meta["inchi"] is None
    assert ds.feature_tags == {"0": ["Confirmed"]}
    assert ds.curation_flags["n"] == 0
    assert meta["annotator_id"] == "Msp20260101000000_lib_1"
    assert meta["total_score"] == 3.2  # 丸めない


def test_adapter_without_tags_has_empty_tags(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "PC 34:1"}], samples=_samples())
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_tags == {}
    assert ds.feature_metadata["0"]["annotation_source"] == "none"
    assert ds.feature_metadata["0"]["spectrum_matched"] is None


def test_broken_flag_file_becomes_a_load_note(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "PC 34:1"}], samples=_samples())
    (tmp_path / "curation").mkdir()
    (tmp_path / "curation" / "flags.jsonl").write_text("{not json\n", encoding="utf-8")
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert any("1 行目" in n for n in ds.load_notes)
    assert ds.feature_annotations  # 読み込みは止まらない
    assert ds.curation_flags is None


def test_corrupt_tag_file_becomes_a_load_note(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "PC 34:1"}], samples=_samples())
    selection = select_alignment(tmp_path)
    (tmp_path / f"{selection.arf2.stem}_tags.xml").write_text("<PeakSpotTags>", encoding="utf-8")
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_tags == {}
    assert any("タグファイル" in n for n in ds.load_notes)
    assert ds.feature_annotations


@pytest.mark.parametrize("name", ["   ", "\t", " Unknown ", " null "])
def test_whitespace_name_is_not_sml(name):
    fields = spot_fields(arf2_spot_raw(spot_id=0, name="x"))
    fields["name"] = name
    assert sml_annotation(fields) is None


@pytest.mark.parametrize("name", ["unknown", "UNKNOWN", " Unknown "])
def test_unknown_name_case_insensitive(name):
    fields = spot_fields(arf2_spot_raw(spot_id=0, name="x"))
    fields["name"] = name
    cols = identity_columns(fields, msms_assigned=False, tables=IdentityCache().tables)
    ref = identity_columns({**fields, "name": "Unknown"}, msms_assigned=False,
                           tables=IdentityCache().tables)
    assert cols == ref
