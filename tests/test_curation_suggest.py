import json

import pytest

from lipidmix.arf2.reader import load_catalog
from lipidmix.curation import evidence, flags, judge, review, suggest
from lipidmix.library import store as library_store
from tests.curation_fixtures import write_suggest_set

TH = judge.resolve_thresholds(None)
OPTIONS = {"wrong": "flagged_or_likely", "unannotated": True, "include_decided": False, "top_n": 5,
           "rt_window": None, "relation_mz_tol": None, "relation_min_r": None}


def _setup(tmp_path, monkeypatch, **kw):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_suggest_set(tmp_path / "neg", **kw)
    s = library_store.open_store(paths["msp"])
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    base = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025, th=TH, file_ids=None,
                             max_traces=12, selection={})
    review.save_review(base)
    flags.FlagStore(flags.curation_dir(paths["arf2"])).append(
        [{"spot_id": 1, "flag": "wrong", "note": ""}], alignment=base["alignment"],
        review_id=base["review_id"], source="user")
    return paths, s, base


@pytest.fixture()
def built(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    result = suggest.run_suggestion(paths["arf2"], base_review=base, store=s, th=TH, options=OPTIONS)
    yield paths, s, base, result
    s.close()


def test_targets_are_the_flagged_and_the_unannotated(built):
    _, _, _, result = built
    kinds = {sp["spot_id"]: sp["target_kind"] for sp in result["spots"]}
    assert kinds == {1: "flagged", 2: "unannotated", 3: "unannotated", 4: "unannotated"}
    assert result["counts"]["targets"] == {"flagged": 1, "likely_wrong": 0, "unannotated": 3}


def test_isotope_explanation_is_preset_for_the_wrong_spot(built):
    _, _, _, result = built
    spot1 = next(sp for sp in result["spots"] if sp["spot_id"] == 1)
    assert spot1["strong"] is True
    preset = next(r for r in spot1["relations"] if r["candidate_id"] == spot1["preset"])
    assert preset["relation"] == "isotope_M+2" and preset["of"] == 0
    # 代表（PC 16:0_18:1）は除かれ、store に無い 2 位（PG 34:1）は reference_unresolved で残る
    names = [c["name"] for c in spot1["candidates"]]
    assert "PC 16:0_18:1" not in names and "PG 34:1" in names


def test_relation_carries_the_partner_eic(built):
    _, _, _, result = built
    spot1 = next(sp for sp in result["spots"] if sp["spot_id"] == 1)
    trace = spot1["relations"][0]["partner_eic"]
    assert trace["file_id"] == 0 and trace["left"] < trace["top"] < trace["right"]
    assert 0 < len(trace["points"]) <= evidence.EIC_MAX_POINTS


def test_research_candidate_for_an_unannotated_spot(built):
    _, _, _, result = built
    spot2 = next(sp for sp in result["spots"] if sp["spot_id"] == 2)
    assert spot2["candidates"][0]["name"] == "PE 18:0_18:2"
    assert spot2["candidates"][0]["source"] == "research"
    assert spot2["candidates"][0]["sum_name"] == "PE 36:2"
    assert spot2["preset"] is None


def test_insource_link_is_strong(built):
    _, _, _, result = built
    spot3 = next(sp for sp in result["spots"] if sp["spot_id"] == 3)
    assert spot3["relations"][0]["relation"] == "insource:-HCOOCH3" and spot3["relations"][0]["strong"]


def test_spot_without_msms_and_relation_has_nothing(built):
    _, _, _, result = built
    spot4 = next(sp for sp in result["spots"] if sp["spot_id"] == 4)
    assert spot4["candidates"] == [] and spot4["relations"] == []


def test_flagged_only_mode(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    result = suggest.run_suggestion(paths["arf2"], base_review=base, store=s, th=TH,
                                    options={**OPTIONS, "wrong": "flagged", "unannotated": False})
    assert [sp["spot_id"] for sp in result["spots"]] == [1]
    s.close()


def test_run_suggestion_rejects_a_base_review_of_another_alignment(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    stale = {**base, "alignment": {**base["alignment"], "alignment_sha256": "other"}}
    with pytest.raises(ValueError, match="アラインメント"):
        suggest.run_suggestion(paths["arf2"], base_review=stale, store=s, th=TH, options=OPTIONS)
    s.close()


def test_run_suggestion_without_param_file_and_arf_uses_defaults(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch, with_param=False, with_arf=False)
    result = suggest.run_suggestion(paths["arf2"], base_review=base, store=s, th=TH, options=OPTIONS)
    assert result["analysis_params"]["source"] == "default"
    assert any("param" in w for w in result["warnings"])
    spot3 = next(sp for sp in result["spots"] if sp["spot_id"] == 3)
    assert spot3["relations"][0]["r"] is None and spot3["relations"][0]["strong"]   # リンクだけで裏付け
    s.close()


def test_latest_review_finds_the_base(built):
    paths, _, base, _ = built
    found = suggest.latest_review(paths["arf2"], base["alignment"]["alignment_sha256"])
    assert found["review_id"] == base["review_id"]


def test_save_load_and_summary(built):
    paths, _, _, result = built
    saved = suggest.save_suggestion(result)
    assert saved["json"].name == f"suggest-{result['suggestion_id']}.json"
    loaded = suggest.load_suggestion(paths["arf2"], result["suggestion_id"])
    assert loaded["suggestion_id"] == result["suggestion_id"]
    lines = suggest.summary_tsv(loaded, max_rows=10).splitlines()
    assert lines[0].split("\t") == suggest.TSV_COLUMNS and len(lines) == 5
    with pytest.raises(ValueError):
        suggest.load_suggestion(paths["arf2"], "../evil")


def test_expand_entries_builds_assign_and_redundant(built):
    _, _, _, result = built
    rows = suggest.expand_entries([
        {"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "sum", "note": "ok"},
        {"spot_id": 1, "flag": "redundant", "candidate": "R1"},
        {"spot_id": 4, "flag": "clear"}], result)
    assign, redundant, clear = rows
    assert assign["name"] == "PE 36:2" and assign["species_name"] == "PE 18:0_18:2"
    assert assign["inchikey"] == "KEY-PE362" and assign["candidate_source"] == "research"
    assert assign["suggestion_id"] == result["suggestion_id"]
    assert redundant["of"] == 0 and redundant["relation"] == "isotope_M+2"
    assert clear == {"spot_id": 4, "flag": "clear", "note": ""}


def test_expand_entries_species_level_keeps_the_library_name(built):
    _, _, _, result = built
    row = suggest.expand_entries([{"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "species"}],
                                 result)[0]
    assert row["name"] == "PE 18:0_18:2" and row["level"] == "species"


def test_expand_entries_uses_msdial_identity_when_reference_unresolved(built):
    _, _, _, result = built
    spot1 = next(sp for sp in result["spots"] if sp["spot_id"] == 1)
    pg = next(c for c in spot1["candidates"] if c["name"] == "PG 34:1")
    row = suggest.expand_entries([{"spot_id": 1, "flag": "assign", "candidate": pg["candidate_id"],
                                   "level": "sum"}], result)[0]
    assert row["name"] == "PG 34:1" and row["inchikey"] == "KEY-PG341"


@pytest.mark.parametrize("entry, message", [
    ({"spot_id": 2, "flag": "assign", "candidate": "L9", "level": "sum"}, "candidate"),
    ({"spot_id": 0, "flag": "assign", "candidate": "L1", "level": "sum"}, "対象外"),
    ({"spot_id": 2, "flag": "wrong"}, "flag"),
    ({"spot_id": 2, "flag": "redundant", "candidate": "L1"}, "candidate"),
    ({"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "chain"}, "level"),
])
def test_expand_entries_rejects_unknown_candidate_and_spot(built, entry, message):
    _, _, _, result = built
    with pytest.raises(ValueError, match=message):
        suggest.expand_entries([entry], result)


def test_decided_spots_are_skipped_by_default(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    flags.FlagStore(flags.curation_dir(paths["arf2"])).append(
        [{"spot_id": 2, "flag": "assign", "name": "PE 36:2", "level": "sum", "note": ""}],
        alignment=base["alignment"], review_id="cs-x", source="user")
    ids = [sp["spot_id"] for sp in suggest.run_suggestion(
        paths["arf2"], base_review=base, store=s, th=TH, options=OPTIONS)["spots"]]
    assert 2 not in ids
    ids = [sp["spot_id"] for sp in suggest.run_suggestion(
        paths["arf2"], base_review=base, store=s, th=TH, options={**OPTIONS, "include_decided": True})["spots"]]
    assert 2 in ids
    s.close()
