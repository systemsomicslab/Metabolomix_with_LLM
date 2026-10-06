import json

import pytest

from metabolomix.curation import flags

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


# ---------- I4: 送信用テキストは .arf2 の絶対パスを運ぶ ----------

def test_submission_text_carries_the_arf2_path(tmp_path):
    arf2 = str((tmp_path / "AlignmentResult_x.arf2").resolve())
    text = flags.build_submission_text("cr-1", [{"spot_id": 3, "flag": "wrong"}], arf2)
    assert json.loads(text[len(flags.SUBMISSION_PREFIX):])["arf2_path"] == arf2
    assert flags.parse_submission_text(text)["arf2_path"] == arf2


def test_submission_text_without_arf2_path_still_parses():
    text = flags.SUBMISSION_PREFIX + json.dumps({"review_id": "cr-1", "flags": []})
    assert flags.parse_submission_text(text)["arf2_path"] is None


# ---------- I6: 壊れた行は名前付きの例外 ----------

def test_a_truncated_line_raises_flag_file_error_naming_the_line(tmp_path):
    store = flags.FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=ALIGN, review_id="r", source="user")
    with open(store.path, "a", encoding="utf-8") as handle:
        handle.write('{"spot_id": 2, "flag": "wro\n')
    with pytest.raises(flags.FlagFileError) as info:
        store.effective(ALIGN["alignment_sha256"])
    assert info.value.path == store.path
    assert info.value.line_no == 2
    assert "flags.jsonl" in str(info.value) and "2" in str(info.value)


def test_a_line_without_an_integer_spot_id_is_a_flag_file_error(tmp_path):
    store = flags.FlagStore(tmp_path)
    store.path.write_text('{"flag": "wrong", "alignment_sha256": "x"}\n', encoding="utf-8")
    with pytest.raises(flags.FlagFileError) as info:
        store.effective("x")
    assert info.value.line_no == 1


# ---------- M2: メモの改行・タブは空白へ ----------

def test_notes_have_control_whitespace_normalised():
    cleaned = flags.validate_entries([{"spot_id": 1, "flag": "wrong", "note": "a\r\nb\tc"}],
                                     allowed_spot_ids=None)
    assert cleaned[0]["note"] == "a  b c"


def test_cleared_spots_are_those_whose_latest_row_is_clear(tmp_path):
    # likely_wrong の「間違い」プリセットを、人が明示的に取り消したスポットへは掛け直さないため。
    store = flags.FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "wrong"}, {"spot_id": 2, "flag": "wrong"}],
                 alignment=ALIGN, review_id="cr-1", source="user")
    store.append([{"spot_id": 1, "flag": "clear"}], alignment=ALIGN, review_id="cr-2", source="user")
    other = {"alignment_file": "AlignmentResult_x.arf2", "alignment_sha256": "bb" * 32}
    store.append([{"spot_id": 3, "flag": "clear"}], alignment=other, review_id="cr-3", source="user")
    assert flags.cleared_spots(store.rows(), ALIGN["alignment_sha256"]) == {1}


from metabolomix.curation.flags import (
    FlagStore, effective_flags, flags_digest, split_decisions, validate_entries)

_ALIGN = {"alignment_file": "a.arf2", "alignment_sha256": "s1"}


def test_assign_after_wrong_wins_and_clear_removes_it(tmp_path):
    store = FlagStore(tmp_path)
    store.append([{"spot_id": 4, "flag": "wrong", "note": ""}], alignment=_ALIGN, review_id="cr-1", source="user")
    store.append([{"spot_id": 4, "flag": "assign", "name": "PC 34:1", "level": "sum", "note": ""}],
                 alignment=_ALIGN, review_id="cs-1", source="user")
    effective = store.effective("s1")
    assert effective[4]["flag"] == "assign"
    assert split_decisions(effective)["assign"][4]["name"] == "PC 34:1"
    assert split_decisions(effective)["wrong"] == set()
    store.append([{"spot_id": 4, "flag": "clear", "note": ""}], alignment=_ALIGN, review_id="cs-1", source="user")
    assert store.effective("s1") == {}


def test_rows_accept_the_new_flag_values(tmp_path):
    store = FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "redundant", "of": 0, "relation": "isotope_M+2", "note": ""}],
                 alignment=_ALIGN, review_id="cs-1", source="user")
    assert store.rows()[0]["flag"] == "redundant"


def test_digest_for_wrong_and_suspect_is_unchanged():
    # 既存ユーザーの curation_flags_sha256 を変えない（メタ行が現行と同じであること）
    import hashlib, json
    effective = {1: {"flag": "wrong"}, 2: {"flag": "suspect"}}
    legacy = hashlib.sha256(json.dumps(sorted([(1, "wrong"), (2, "suspect")])).encode()).hexdigest()
    assert flags_digest(effective) == legacy


def test_digest_changes_when_the_assigned_name_changes():
    a = {1: {"flag": "assign", "name": "PC 34:1", "level": "sum"}}
    b = {1: {"flag": "assign", "name": "PC 34:2", "level": "sum"}}
    assert flags_digest(a) != flags_digest(b)


def test_digest_changes_when_only_the_assigned_record_changes():
    a = {1: {"flag": "assign", "name": "PC 34:1", "level": "sum", "inchikey": "KEY-A"}}
    b = {1: {"flag": "assign", "name": "PC 34:1", "level": "sum", "inchikey": "KEY-B"}}
    assert flags_digest(a) != flags_digest(b)


def test_review_submissions_still_reject_assign():
    with pytest.raises(ValueError, match="flag"):
        validate_entries([{"spot_id": 1, "flag": "assign"}], allowed_spot_ids={1})
