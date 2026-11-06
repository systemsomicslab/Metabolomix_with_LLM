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
