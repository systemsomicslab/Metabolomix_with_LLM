"""MatchResults（.arf2 Key 56）の復号と、GUI が表示する代表の選択。"""
import io
import math

from metabolomix.arf2 import match_results as mr
from metabolomix.arf2 import reader as arf2_reader
from tests.curation_fixtures import arf2_spot_raw, match_result, write_arf2


def test_decode_names_every_key_from_the_schema():
    decoded = mr.decode_match_result(match_result())
    assert decoded["name"] == "PC 34:1"
    assert decoded["library_id"] == 5
    assert decoded["annotator_id"] == "Msp20260101000000_lib_1"
    assert decoded["is_reference_matched"] is True
    assert decoded["squared_weighted_dot_product"] == 0.81
    assert decoded["priority"] == 1
    assert decoded["is_manually_modified"] is False
    assert decoded["has_msms"] is True


def test_decode_marks_msms_absent_by_the_minus_one_sentinel():
    decoded = mr.decode_match_result(match_result({3: -1.0, 4: -1.0, 5: -1.0, 6: -1.0}))
    assert decoded["has_msms"] is False


def test_decode_rejects_rows_that_are_not_match_results():
    assert mr.decode_match_result(None) is None
    assert mr.decode_match_result([1, 2]) is None


def test_representative_follows_upstream_result_order():
    low = match_result({0: "A", 2: 9.0, 33: False, 34: True})
    matched = match_result({0: "B", 2: 1.0, 33: True})
    manual = match_result({0: "C", 2: 0.1, 33: False, 26: 4 | 64})
    assert mr.representative([[low, matched], {}, []])["name"] == "B"
    assert mr.representative([[low, matched, manual], {}, []])["name"] == "C"


def test_representative_ignores_decoys_and_unknown_source():
    decoy = match_result({0: "D", 30: True, 2: 99.0})
    unknown = match_result({0: "U", 26: 1, 2: 99.0})
    real = match_result({0: "R", 2: 0.5})
    assert mr.representative([[decoy, unknown, real], {}, []])["name"] == "R"
    assert mr.representative([[decoy, unknown], {}, []]) is None
    assert mr.representative([[], {}, []]) is None
    assert mr.representative(None) is None


def test_representative_priority_breaks_ties_before_total_score():
    a = match_result({0: "A", 31: 2, 2: 0.1})
    b = match_result({0: "B", 31: 1, 2: 9.0})
    assert mr.representative([[a, b], {}, []])["name"] == "A"


def test_name_prefix_reads_msdial_qualifiers():
    assert mr.name_prefix("low score: PC 34:1") == "low score"
    assert mr.name_prefix("no MS2: PC 34:1") == "no MS2"
    assert mr.name_prefix("Unsettled: PC 34:1") == "unsettled"
    assert mr.name_prefix("w/o MS2: PC 34:1") == "w/o MS2"
    assert mr.name_prefix("PC 34:1") is None
    assert mr.name_prefix(None) is None


def test_name_prefix_handles_extra_whitespace():
    # Test with extra whitespace between qualifier tokens
    assert mr.name_prefix("no  MS2: PC 34:1") == "no MS2"
    assert mr.name_prefix("low   score: PC 34:1") == "low score"
    assert mr.name_prefix("w/o  MS2: PC 34:1") == "w/o MS2"


def test_load_spot_annotations_keys_by_master_alignment_id(tmp_path):
    path = write_arf2(tmp_path / "AlignmentResult_x.arf2", [
        arf2_spot_raw(spot_id=0, matches=[match_result()], representative_file_id=2),
        arf2_spot_raw(spot_id=1, name="Unknown", matches=[]),
    ])
    annotations = mr.load_spot_annotations(path)
    assert set(annotations) == {0, 1}
    assert annotations[0]["representative"]["name"] == "PC 34:1"
    assert annotations[0]["representative_file_id"] == 2
    assert annotations[0]["n_candidates"] == 1
    assert annotations[1]["representative"] is None


def test_iter_raw_spots_keeps_deserialize_behaviour(tmp_path):
    path = write_arf2(tmp_path / "a.arf2", [arf2_spot_raw(spot_id=0), arf2_spot_raw(spot_id=1)])
    spots = arf2_reader.deserialize(io.BytesIO(path.read_bytes()))
    assert [s["MasterAlignmentID"] for s in spots] == [0, 1]
    assert len(arf2_reader.load_raw_spots(path)) == 2


def test_nan_scores_are_passed_through_as_none():
    decoded = mr.decode_match_result(match_result({37: float("nan")}))
    assert decoded["enhanced_dot_product"] is None
    assert not any(isinstance(v, float) and math.isnan(v) for v in decoded.values())


def test_usable_candidates_drop_decoys_and_unknown_source():
    container = [[match_result({0: "A", 31: 2}), match_result({0: "B", 30: True}),
                  match_result({0: "C", 26: 1})], {}, []]
    assert [c["name"] for c in mr.usable_candidates(container)] == ["A"]
