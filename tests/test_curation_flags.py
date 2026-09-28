import json

import pytest

from lipidmix.curation import flags

ALIGN = {"alignment_file": "AlignmentResult_x.arf2", "alignment_sha256": "aa" * 32}


def test_latest_flag_wins_and_clear_removes(tmp_path):
    store = flags.FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "suspect"}, {"spot_id": 2, "flag": "wrong"}],
                 alignment=ALIGN, review_id="cr-1", source="user")
    store.append([{"spot_id": 1, "flag": "wrong", "note": "EIC が二峰"},
                  {"spot_id": 2, "flag": "clear"}],
                 alignment=ALIGN, review_id="cr-2", source="user")
    effective = store.effective(ALIGN["alignment_sha256"])
    assert set(effective) == {1}
    assert effective[1]["flag"] == "wrong"
    assert effective[1]["note"] == "EIC が二峰"


def test_file_is_append_only(tmp_path):
    store = flags.FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=ALIGN, review_id="r", source="user")
    first = (tmp_path / "flags.jsonl").read_text(encoding="utf-8")
    store.append([{"spot_id": 1, "flag": "clear"}], alignment=ALIGN, review_id="r", source="user")
    assert (tmp_path / "flags.jsonl").read_text(encoding="utf-8").startswith(first)


def test_flags_of_another_alignment_do_not_apply(tmp_path):
    store = flags.FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=ALIGN, review_id="r", source="user")
    assert store.effective("bb" * 32) == {}


def test_invalid_entry_rejects_the_whole_submission():
    with pytest.raises(ValueError):
        flags.validate_entries([{"spot_id": 1, "flag": "wrong"}, {"spot_id": 2, "flag": "maybe"}],
                               allowed_spot_ids={1, 2})
    with pytest.raises(ValueError):
        flags.validate_entries([{"spot_id": 9, "flag": "wrong"}], allowed_spot_ids={1, 2})
    with pytest.raises(ValueError):
        flags.validate_entries([{"spot_id": "x", "flag": "wrong"}], allowed_spot_ids=None)


def test_submission_text_round_trips_through_surrounding_chat_text():
    text = flags.build_submission_text("cr-1", [{"spot_id": 3, "flag": "wrong", "note": "a\tb"}])
    pasted = "これを送ります\n\n" + text + "\nよろしく"
    parsed = flags.parse_submission_text(pasted)
    assert parsed["review_id"] == "cr-1"
    assert parsed["flags"] == [{"spot_id": 3, "flag": "wrong", "note": "a\tb"}]


def test_submission_text_without_marker_is_rejected():
    with pytest.raises(ValueError):
        flags.parse_submission_text(json.dumps({"review_id": "x", "flags": []}))


def test_digest_changes_with_effective_flags(tmp_path):
    store = flags.FlagStore(tmp_path)
    empty = store.digest(ALIGN["alignment_sha256"])
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=ALIGN, review_id="r", source="user")
    assert store.digest(ALIGN["alignment_sha256"]) != empty


def test_alignment_key_hashes_the_file(tmp_path):
    path = tmp_path / "AlignmentResult_x.arf2"
    path.write_bytes(b"abc")
    key = flags.alignment_key(path)
    assert key["alignment_file"] == "AlignmentResult_x.arf2"
    assert len(key["alignment_sha256"]) == 64
    assert flags.curation_dir(path) == tmp_path / "curation"
