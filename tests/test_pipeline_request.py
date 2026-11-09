"""pipeline-request.v3 の厳密な解決・更新契約（plan 2026-11-09「要求」）を検証する。

resolve_request / merge_updates / validate_request / request_fingerprint / unit_request は
すべて metabolomix.pipeline.request にある。実 raw や既存成果物は使わない
（tmp_path 上に空の source_root を作るだけ）。
"""
import copy
import json

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline.request import (
    UPDATABLE,
    merge_updates,
    pin_polarity,
    request_fingerprint,
    resolve_request,
    unit_request,
    validate_request,
)


def _root(tmp_path, name="source"):
    root = tmp_path / name
    root.mkdir()
    return root


def _source_with_request_file(tmp_path, content) -> object:
    root = _root(tmp_path)
    text = content if isinstance(content, str) else json.dumps(content)
    (root / "analysis-request.json").write_text(text, encoding="utf-8")
    return root


def _code(fn) -> DomainError:
    with pytest.raises(DomainError) as info:
        fn()
    return info.value


# ---------- 2 つの形 → 1 つの正規形 ----------

def test_single_polarity_form_normalizes_to_one_unit(tmp_path):
    req = resolve_request(_root(tmp_path), {"polarity": "negative",
                                            "method_file": "lab_param.txt",
                                            "sample_manifest": "sheet.tsv"})
    assert req["schema"] == "pipeline-request.v3"
    assert req["units"] == [{
        "polarity": "negative", "raw_dir": ".", "method_file": "lab_param.txt",
        "lbm_file": None, "msp_file": None, "sample_manifest": "sheet.tsv",
        "keep_extension": None,
    }]
    # 単位のキーはトップに残らない（読む側は units だけを見る）。
    for key in ("polarity", "method_file", "lbm_file", "msp_file", "sample_manifest"):
        assert key not in req
    assert req["value_sources"]["units"][0]["polarity"] == "explicit"
    assert req["value_sources"]["units"][0]["lbm_file"] == "default"


def test_two_units_are_accepted_in_request_order(tmp_path):
    req = resolve_request(_root(tmp_path), {"units": [
        {"polarity": "negative", "raw_dir": "raw/neg"},
        {"polarity": "positive", "raw_dir": "raw\\pos", "method_file": "m_pos.txt"},
    ]})
    assert [u["polarity"] for u in req["units"]] == ["negative", "positive"]
    assert [u["raw_dir"] for u in req["units"]] == ["raw/neg", "raw/pos"]
    assert req["units"][1]["method_file"] == "m_pos.txt"
    assert req["value_sources"]["units"][1]["method_file"] == "explicit"


def test_two_units_require_distinct_polarity_and_raw_dir(tmp_path):
    root = _root(tmp_path)
    same_polarity = _code(lambda: resolve_request(root, {"units": [
        {"polarity": "negative", "raw_dir": "neg"},
        {"polarity": "negative", "raw_dir": "neg2"}]}))
    assert same_polarity.code == "PIPELINE_REQUEST_INVALID"
    same_dir = _code(lambda: resolve_request(root, {"units": [
        {"polarity": "negative", "raw_dir": "raw"},
        {"polarity": "positive", "raw_dir": "raw/"}]}))
    assert same_dir.code == "PIPELINE_REQUEST_INVALID"
    for missing in ({"raw_dir": "neg"}, {"polarity": "negative"}):
        err = _code(lambda: resolve_request(root, {"units": [missing]}))
        assert err.code == "PIPELINE_REQUEST_INVALID"


@pytest.mark.parametrize("units", [
    [],
    [{"polarity": "negative", "raw_dir": "a"}, {"polarity": "positive", "raw_dir": "b"},
     {"polarity": "positive", "raw_dir": "c"}],
    [{"polarity": "both", "raw_dir": "a"}],
    [{"polarity": "negative", "raw_dir": "../outside"}],
    [{"polarity": "negative", "raw_dir": "C:/abs"}],
    [{"polarity": "negative", "raw_dir": "a", "unknown": 1}],
    [{"polarity": "negative", "raw_dir": None}],
    "negative",
])
def test_malformed_units_are_rejected(tmp_path, units):
    err = _code(lambda: resolve_request(_root(tmp_path), {"units": units}))
    assert err.code == "PIPELINE_REQUEST_INVALID"


def test_units_and_top_level_unit_keys_cannot_be_mixed(tmp_path):
    err = _code(lambda: resolve_request(_root(tmp_path), {
        "polarity": "negative", "units": [{"polarity": "negative", "raw_dir": "neg"}]}))
    assert err.code == "PIPELINE_REQUEST_INVALID"


# ---------- 下流のキー ----------

@pytest.mark.parametrize("key, value", [
    ("target", "differential"),
    ("preprocess", {"normalize": "tic"}),
    ("comparisons", [{"comparison_id": "a", "reference_group": "x", "test_group": "y"}]),
])
def test_downstream_keys_are_rejected_with_replacements(tmp_path, key, value):
    err = _code(lambda: resolve_request(_root(tmp_path), {key: value}))
    assert err.code == "PIPELINE_DOWNSTREAM_REMOVED"
    assert err.details["removed_keys"] == [key]
    assert "dataset_compare" in err.details["replacements"]
    assert "load_dataset" in err.details["replacements"]


def test_request_file_downstream_keys_are_ignored_with_warning(tmp_path):
    """reanalysis-study の古い要求（`target`・`comparisons` を持つ analysis-request.json）。"""
    root = _source_with_request_file(tmp_path, {
        "target": "differential", "polarity": "negative", "timeout_s": 60,
        "comparisons": [{"comparison_id": "a", "reference_group": "x", "test_group": "y"}],
        "preprocess": {"normalize": "tic"}})

    req = resolve_request(root)

    assert req["units"][0]["polarity"] == "negative"  # 上流の値は生きる
    assert req["timeout_s"] == 60
    for key in ("target", "comparisons", "preprocess"):
        assert key not in req
    [warning] = req["intake_warnings"]
    assert warning["code"] == "REQUEST_FILE_DOWNSTREAM_IGNORED"
    assert warning["details"]["ignored_keys"] == ["comparisons", "preprocess", "target"]


def test_no_request_file_has_no_intake_warning(tmp_path):
    assert resolve_request(_root(tmp_path))["intake_warnings"] == []


# ---------- schema ----------

@pytest.mark.parametrize("explicit", [{}, {"schema": "pipeline-request.v1"},
                                      {"schema": "pipeline-request.v3"}])
def test_v1_v3_and_omitted_schema_are_saved_as_v3(tmp_path, explicit):
    assert resolve_request(_root(tmp_path), explicit)["schema"] == "pipeline-request.v3"


def test_v2_schema_is_rejected(tmp_path):
    err = _code(lambda: resolve_request(_root(tmp_path), {"schema": "pipeline-request.v2"}))
    assert err.code == "PIPELINE_REQUEST_INVALID"
    assert err.details["removed_schema"] == "pipeline-request.v2"


def test_v2_schema_in_the_request_file_is_rejected(tmp_path):
    root = _source_with_request_file(tmp_path, {"schema": "pipeline-request.v2"})
    err = _code(lambda: resolve_request(root))
    assert err.details["removed_schema"] == "pipeline-request.v2"


def test_unknown_schema_is_rejected(tmp_path):
    err = _code(lambda: resolve_request(_root(tmp_path), {"schema": "pipeline-request.v9"}))
    assert err.code == "PIPELINE_REQUEST_INVALID"


# ---------- 既定値・出所 ----------

def test_defaults_are_fully_populated(tmp_path):
    req = resolve_request(_root(tmp_path))
    assert req["library_mode"] == "auto"
    assert req["measure"] == "peak_height"
    assert req["keep_extension"] is None
    assert req["timeout_s"] == 21600
    assert req["save_project"] is True
    assert req["output_root"] is None
    assert req["alignment_reference"] == "median_injection"
    assert req["units"] == [{"polarity": None, "raw_dir": ".", "method_file": None,
                             "lbm_file": None, "msp_file": None, "sample_manifest": None,
                             "keep_extension": None}]
    assert req["value_sources"]["timeout_s"] == "default"
    assert all(v == "default" for v in req["value_sources"]["units"][0].values())


def test_an_explicit_value_beats_the_request_file(tmp_path):
    root = _source_with_request_file(tmp_path, {"polarity": "positive", "timeout_s": 60,
                                                "method_file": "from_file.txt"})
    req = resolve_request(root, {"polarity": "negative"})
    unit = req["units"][0]
    assert unit["polarity"] == "negative"
    assert unit["method_file"] == "from_file.txt"    # キーごとに重なる（1 極性の形）
    assert req["timeout_s"] == 60
    assert req["value_sources"]["units"][0]["polarity"] == "explicit"
    assert req["value_sources"]["units"][0]["method_file"] == "request_file"
    assert req["value_sources"]["timeout_s"] == "request_file"


def test_units_from_the_request_file_are_taken_whole(tmp_path):
    root = _source_with_request_file(tmp_path, {"units": [
        {"polarity": "negative", "raw_dir": "neg"}, {"polarity": "positive", "raw_dir": "pos"}]})
    req = resolve_request(root, {"timeout_s": 30})
    assert [u["raw_dir"] for u in req["units"]] == ["neg", "pos"]
    assert req["value_sources"]["units"][0]["raw_dir"] == "request_file"
    # ファイルが units を持つのに単位のキーを明示すると、どの極性に効くか決められない。
    err = _code(lambda: resolve_request(root, {"method_file": "m.txt"}))
    assert err.code == "PIPELINE_REQUEST_INVALID"
    # 明示の units はファイルの units を丸ごと置き換える。
    req = resolve_request(root, {"units": [{"polarity": "positive", "raw_dir": "pos"}]})
    assert [u["polarity"] for u in req["units"]] == ["positive"]


@pytest.mark.parametrize("content", ["{ not json", ["x"], {"targett": 1},
                                     {"method_file": None}])
def test_a_broken_request_file_is_rejected_not_ignored(tmp_path, content):
    root = _source_with_request_file(tmp_path, content)
    assert _code(lambda: resolve_request(root)).code == "PIPELINE_REQUEST_INVALID"


def test_a_request_file_value_participates_in_the_fingerprint(tmp_path):
    with_file = resolve_request(_source_with_request_file(tmp_path, {"timeout_s": 60}))
    plain = resolve_request(_root(tmp_path, "plain"))
    assert request_fingerprint(with_file) != request_fingerprint(plain)


# ---------- 値の検査 ----------

@pytest.mark.parametrize("explicit", [
    {"timeout_s": True}, {"timeout_s": 0}, {"timeout_s": 1.5},
    {"polarity": "both"}, {"measure": "peak_area"}, {"keep_extension": ""},
    {"save_project": "true"}, {"not_a_real_key": 1},
    {"value_sources": {}}, {"intake_warnings": []},
])
def test_invalid_values_are_rejected(tmp_path, explicit):
    assert _code(lambda: resolve_request(_root(tmp_path), explicit)).code \
        == "PIPELINE_REQUEST_INVALID"


@pytest.mark.parametrize("field", ["method_file", "lbm_file", "polarity", "keep_extension",
                                   "output_root", "library_mode", "msp_file",
                                   "alignment_reference", "units"])
def test_explicit_null_is_rejected_for_non_disable_fields(tmp_path, field):
    assert _code(lambda: resolve_request(_root(tmp_path), {field: None})).code \
        == "PIPELINE_REQUEST_INVALID"


def test_sample_manifest_explicit_null_is_permitted(tmp_path):
    req = resolve_request(_root(tmp_path), {"sample_manifest": None})
    assert req["units"][0]["sample_manifest"] is None
    assert req["value_sources"]["units"][0]["sample_manifest"] == "explicit"


def test_absolute_paths_outside_the_source_root_are_accepted(tmp_path):
    outside = str(tmp_path / "elsewhere")
    req = resolve_request(_root(tmp_path), {"output_root": outside,
                                            "method_file": "C:/lab/methods/lipid_NEG.txt"})
    assert req["output_root"] == outside
    assert req["units"][0]["method_file"] == "C:/lab/methods/lipid_NEG.txt"


# ---------- library_mode / alignment_reference（v8） ----------

@pytest.mark.parametrize("explicit", [
    {"library_mode": "lbm_only"},
    {"library_mode": "msp_only", "msp_file": ""},
    {"library_mode": "msp_only", "lbm_file": "x.lbm2"},
    {"library_mode": "msp_then_lbm", "lbm_file": "x.lbm2"},
    {"msp_file": "lab.msp"},  # auto では使われないので黙って受けない
    {"library_mode": "auto",
     "units": [{"polarity": "negative", "raw_dir": "n", "msp_file": "lab.msp"}]},
    {"alignment_reference": "first"}, {"alignment_reference": 3},
])
def test_invalid_library_and_alignment_settings_are_rejected(tmp_path, explicit):
    assert _code(lambda: resolve_request(_root(tmp_path), explicit)).code \
        == "PIPELINE_REQUEST_INVALID"


def test_library_mode_and_alignment_reference_fingerprints(tmp_path):
    root = _root(tmp_path)
    fingerprints = {request_fingerprint(resolve_request(root, explicit)) for explicit in
                    ({}, {"library_mode": "msp_only"}, {"library_mode": "msp_then_lbm"},
                     {"alignment_reference": "method"})}
    assert len(fingerprints) == 4
    # 既定値のキーは数えない（キーを足しただけで hash を動かさない）。
    req = resolve_request(root)
    stripped = {k: v for k, v in req.items() if k not in ("library_mode", "alignment_reference")}
    assert request_fingerprint(stripped) == request_fingerprint(req)


def test_request_fingerprint_ignores_provenance_only_differences(tmp_path):
    root = _root(tmp_path)
    explicit = resolve_request(root, {"timeout_s": 21600})
    default = resolve_request(root)
    assert explicit["value_sources"] != default["value_sources"]
    assert request_fingerprint(explicit) == request_fingerprint(default)


# ---------- resume の訂正（試料の対応表だけ） ----------

def test_sample_manifest_can_be_updated_on_a_single_unit_run(tmp_path):
    req = resolve_request(_root(tmp_path), {"sample_manifest": "a.tsv"})
    updated = merge_updates(req, {"sample_manifest": "b.tsv"})
    assert updated["units"][0]["sample_manifest"] == "b.tsv"
    assert updated["value_sources"]["units"][0]["sample_manifest"] == "explicit_update"
    released = merge_updates(updated, {"sample_manifest": None})
    assert released["units"][0]["sample_manifest"] is None
    assert UPDATABLE == {"sample_manifest", "units"}


def test_unit_sample_manifest_is_updated_by_polarity(tmp_path):
    req = resolve_request(_root(tmp_path), {"units": [
        {"polarity": "negative", "raw_dir": "neg"}, {"polarity": "positive", "raw_dir": "pos"}]})
    updated = merge_updates(req, {"units": {"positive": {"sample_manifest": "pos/s.tsv"}}})
    assert updated["units"][1]["sample_manifest"] == "pos/s.tsv"
    assert updated["units"][0]["sample_manifest"] is None
    # 2 極性ではトップの sample_manifest はどちらの極性か決められない。
    assert _code(lambda: merge_updates(req, {"sample_manifest": "x.tsv"})).code \
        == "PIPELINE_REQUEST_INVALID"
    assert _code(lambda: merge_updates(
        req, {"units": {"neutral": {"sample_manifest": "x"}}})).code == "PIPELINE_REQUEST_INVALID"


@pytest.mark.parametrize("updates", [
    {"method_file": "other.txt"}, {"library_mode": "msp_only"},
    {"alignment_reference": "method"}, {"polarity": "positive"},
    {"value_sources": {}}, {"units": {"negative": {"method_file": "x.txt"}}},
])
def test_resume_cannot_change_upstream(tmp_path, updates):
    req = resolve_request(_root(tmp_path), {"polarity": "negative"})
    before = copy.deepcopy(req)
    assert _code(lambda: merge_updates(req, updates)).code == "NEW_PIPELINE_REQUIRED"
    assert req == before


@pytest.mark.parametrize("updates", [
    {"comparisons": []}, {"target": "differential"}, {"preprocess": {"impute": "knn"}},
    {"units": {"negative": {"comparisons": []}}},
])
def test_resume_downstream_updates_are_refused_with_replacements(tmp_path, updates):
    req = resolve_request(_root(tmp_path), {"polarity": "negative"})
    err = _code(lambda: merge_updates(req, updates))
    assert err.code == "PIPELINE_DOWNSTREAM_REMOVED"
    assert "dataset_compare" in err.details["replacements"]


# ---------- 単位の平らな要求・極性の確定 ----------

def test_unit_request_flattens_shared_and_unit_keys(tmp_path):
    root = _root(tmp_path)
    req = resolve_request(root, {"keep_extension": "wiff", "units": [
        {"polarity": "negative", "raw_dir": "neg", "method_file": "methods/m.txt",
         "sample_manifest": "neg/sheet.tsv"},
        {"polarity": "positive", "raw_dir": "pos", "keep_extension": "raw",
         "method_file": "C:/abs/m.txt"}]})
    neg = unit_request(req, 0, root)
    assert neg["polarity"] == "negative" and neg["keep_extension"] == "wiff"
    # raw_dir 起点で読まれるので、相対パスは dataset_root 基準の絶対パスにする。
    assert neg["method_file"] == str(root / "methods/m.txt")
    assert neg["sample_manifest"] == str(root / "neg/sheet.tsv")
    pos = unit_request(req, 1, root)
    assert pos["keep_extension"] == "raw" and pos["method_file"] == "C:/abs/m.txt"
    # 1 極性の形（raw_dir "."）は今までどおりそのまま。
    single = resolve_request(root, {"method_file": "m.txt"})
    assert unit_request(single, 0, root)["method_file"] == "m.txt"


def test_pin_polarity_fills_only_an_omitted_polarity(tmp_path):
    req = resolve_request(_root(tmp_path))
    pinned = pin_polarity(req, 0, "negative")
    assert pinned["units"][0]["polarity"] == "negative"
    assert pinned["value_sources"]["units"][0]["polarity"] == "method_declaration"
    explicit = resolve_request(_root(tmp_path, "b"), {"polarity": "positive"})
    assert pin_polarity(explicit, 0, "negative")["units"][0]["polarity"] == "positive"


def test_validate_request_rejects_unknown_key_directly(tmp_path):
    req = resolve_request(_root(tmp_path))
    req.pop("value_sources")
    req.pop("intake_warnings")
    validate_request(req)
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_request({**req, "target": "auto"})
