from lipidmix.arf2.ion_features import LINK_KINDS, decode_ion_features, load_ion_features
from tests.curation_fixtures import arf2_spot_raw, write_arf2


def test_decode_links_and_isotope_fields():
    raw = arf2_spot_raw(spot_id=3, peak_links=[(0, 4), (1, 5), (2, 9)])
    decoded = decode_ion_features(raw)
    assert decoded["links"] == [{"spot_id": 0, "kind": "found_in_upper_msms"},
                                {"spot_id": 1, "kind": "correl_similar"},
                                {"spot_id": 2, "kind": "unknown_9"}]
    assert decoded["isotope_weight"] == 0
    assert decoded["peak_group"] == -1


def test_decode_tolerates_missing_key_10():
    raw = arf2_spot_raw(spot_id=0)
    raw[10] = None
    assert decode_ion_features(raw) == {"links": [], "isotope_weight": None,
                                        "isotope_parent": None, "peak_group": None}


def test_load_ion_features_is_keyed_by_master_alignment_id(tmp_path):
    path = write_arf2(tmp_path / "a.arf2", [arf2_spot_raw(spot_id=0),
                                            arf2_spot_raw(spot_id=1, peak_links=[(0, 3)])])
    features = load_ion_features(path)
    assert features[1]["links"] == [{"spot_id": 0, "kind": "chrom_similar"}]
    assert features[0]["links"] == []


def test_link_kinds_match_the_upstream_enum():
    assert LINK_KINDS == {0: "same_feature", 1: "isotope", 2: "adduct", 3: "chrom_similar",
                          4: "found_in_upper_msms", 5: "correl_similar"}
