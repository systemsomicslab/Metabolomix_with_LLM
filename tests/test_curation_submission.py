import pytest

from metabolomix.curation import flags, submission
from metabolomix.library import store as library_store
from metabolomix.msdial import tags as msdial_tags
from tests.curation_fixtures import build_saved_review


@pytest.fixture()
def saved(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    return build_saved_review(tmp_path / "neg")


def _load(saved):
    r = saved["review"]
    return submission.load_saved(flags.curation_dir(r["arf2_path"]), r["review_id"])


def test_submit_flags_records_and_updates_the_tags_file(saved):
    r = saved["review"]
    body = submission.submit_flags(_load(saved), [{"spot_id": 1, "flag": "wrong", "note": "x"}],
                                   review_id=r["review_id"], source="user")
    assert body["status"] == "ok" and body["recorded"] == 1 and body["n_wrong"] == 1
    assert body["tags_xml"]["added"] == [1]
    rows = flags.FlagStore(flags.curation_dir(r["arf2_path"])).rows()
    assert rows[-1]["source"] == "user" and rows[-1]["review_id"] == r["review_id"]
    tag_path = msdial_tags.alignment_tag_path(saved["paths"]["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {1: frozenset({3})}


def test_invalid_entries_raise_invalid_without_writing(saved):
    r = saved["review"]
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 1, "flag": "nope"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "invalid"
    assert not flags.FlagStore(flags.curation_dir(r["arf2_path"])).exists()


def test_duplicate_spot_ids_are_invalid(saved):
    r = saved["review"]
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 0, "flag": "wrong"}, {"spot_id": 0, "flag": "suspect"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "invalid" and "[0]" in str(info.value)


def test_alignment_change_is_its_own_kind(saved):
    r = saved["review"]
    with open(saved["paths"]["arf2"], "ab") as handle:
        handle.write(b"\x00")
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 0, "flag": "wrong"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "alignment_changed"


def test_a_corrupt_flags_file_is_flag_file_with_details(saved):
    r = saved["review"]
    path = flags.curation_dir(r["arf2_path"]) / flags.FLAGS_FILENAME
    path.write_text("not json\n", encoding="utf-8")
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 0, "flag": "wrong"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "flag_file"
    assert info.value.details == {"flags_file": str(path), "line": 1}
    assert path.read_text(encoding="utf-8") == "not json\n"


def test_load_saved_raises_file_not_found_for_a_missing_review(saved):
    r = saved["review"]
    with pytest.raises(FileNotFoundError):
        submission.load_saved(flags.curation_dir(r["arf2_path"]), "cr-20000101-000000-abcd")
