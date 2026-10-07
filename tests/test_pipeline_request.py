"""pipeline-request.v1 の厳密な解決・更新契約（spec §10.1）を検証する。

resolve_request/merge_updates/validate_request/request_fingerprint はすべて
metabolomix.pipeline.request にある。実rawや既存成果物は使わない
（tmp_path上に空のsource_rootを作るだけ）。
"""
import copy
import math

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline.request import (
    UPDATABLE,
    merge_updates,
    request_fingerprint,
    resolve_request,
    validate_request,
)


# ---------- brief記載のRED ----------

def test_null_disable_is_distinct_from_omitted(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    original = resolve_request(root)
    assert original["preprocess"]["blank_min_fold"] == "auto"
    changed = merge_updates(original, {"preprocess": {"blank_min_fold": None}})
    assert changed["preprocess"]["blank_min_fold"] is None
    assert changed["preprocess"]["normalize"] == "auto"


def test_resume_cannot_change_upstream(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(resolve_request(root), {"method_file": "other.txt"})


# ---------- 既定値・出所 ----------

def test_defaults_are_fully_populated(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)

    assert req["schema"] == "pipeline-request.v1"
    assert req["target"] == "auto"
    assert req["effective_target"] == "exploratory"  # comparisons空 => exploratory
    assert req["method_file"] is None
    assert req["lbm_file"] is None
    assert req["library_mode"] == "auto"
    assert req["msp_file"] is None
    assert req["polarity"] is None
    assert req["measure"] == "peak_height"
    assert req["keep_extension"] is None
    assert req["timeout_s"] == 21600
    assert req["save_project"] is True
    assert req["output_root"] is None
    assert req["sample_manifest"] is None
    assert req["comparisons"] == []
    assert req["preprocess"] == {
        "policy": "conservative-v1", "normalize": "auto", "blank_min_fold": "auto",
        "drift_correct": "auto", "max_qc_rsd": "auto", "impute": "half_min",
        "min_detection_rate": 0.0,
    }
    # 未指定は全フィールドdefault
    assert req["value_sources"]["target"] == "default"
    assert req["value_sources"]["sample_manifest"] == "default"
    assert all(v == "default" for v in req["value_sources"]["preprocess"].values())
    assert req["value_sources"]["comparisons"] == "default"


def test_explicit_overrides_default_and_source_is_explicit(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"polarity": "negative", "preprocess": {"impute": "knn"}})
    assert req["polarity"] == "negative"
    assert req["value_sources"]["polarity"] == "explicit"
    # 未指定の他preprocessキーはdefaultのまま
    assert req["preprocess"]["impute"] == "knn"
    assert req["preprocess"]["normalize"] == "auto"
    assert req["value_sources"]["preprocess"]["impute"] == "explicit"
    assert req["value_sources"]["preprocess"]["normalize"] == "default"


def test_explicit_gt_file_gt_default_resolution_order(tmp_path):
    """明示値 > (前revisionの)ファイル的な既存値 > 既定、の3段の優先順位。"""
    root = tmp_path / "source"
    root.mkdir()
    base = resolve_request(root, {"preprocess": {"normalize": "tic"}})
    assert base["preprocess"]["normalize"] == "tic"  # 既定autoではなく明示値

    # merge_updatesの`request`引数は「保存済みファイルから読み戻した値」を表す。
    # ここでの明示updatesがそれをさらに上書きすることを確認する。
    updated = merge_updates(base, {"preprocess": {"normalize": "median"}})
    assert updated["preprocess"]["normalize"] == "median"
    assert updated["value_sources"]["preprocess"]["normalize"] == "explicit_update"
    # 更新していない他のpreprocessキーの出所は維持される
    assert updated["value_sources"]["preprocess"]["impute"] == "default"


# ---------- bool/NaN/Infinityの拒否 ----------

def test_timeout_s_rejects_bool(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"timeout_s": True})


def test_preprocess_threshold_rejects_bool(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"preprocess": {"blank_min_fold": True}})


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_preprocess_rejects_nan_and_infinity(tmp_path, bad):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"preprocess": {"max_qc_rsd": bad}})


def test_comparison_threshold_rejects_nan(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"comparisons": [{
            "comparison_id": "a_vs_b", "reference_group": "a", "test_group": "b",
            "q_threshold": math.nan,
        }]})


# ---------- target 固定/更新 ----------

def test_target_auto_resolves_to_differential_when_comparisons_present(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"comparisons": [{
        "comparison_id": "t_vs_c", "reference_group": "control", "test_group": "treated",
    }]})
    assert req["target"] == "auto"
    assert req["effective_target"] == "differential"


def test_target_can_be_updated_via_resume(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    assert req["effective_target"] == "exploratory"
    updated = merge_updates(req, {"target": "differential"})
    assert updated["target"] == "differential"
    assert updated["effective_target"] == "differential"
    assert updated["value_sources"]["target"] == "explicit_update"


# ---------- sample_manifest: null(明示解除) と 未指定(既定探索) の区別 ----------

def test_sample_manifest_omitted_means_default_search(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    assert req["sample_manifest"] is None
    assert req["value_sources"]["sample_manifest"] == "default"


def test_sample_manifest_explicit_path_then_explicit_null_release(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"sample_manifest": "sample-manifest.tsv"})
    assert req["sample_manifest"] == "sample-manifest.tsv"
    assert req["value_sources"]["sample_manifest"] == "explicit"

    released = merge_updates(req, {"sample_manifest": None})
    assert released["sample_manifest"] is None
    assert released["value_sources"]["sample_manifest"] == "explicit_update"


# ---------- preprocessの深い更新 ----------

def test_merge_updates_preprocess_updates_only_named_child_keys(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"preprocess": {
        "normalize": "tic", "impute": "knn", "min_detection_rate": 0.1,
    }})
    updated = merge_updates(req, {"preprocess": {"impute": "column_mean"}})
    assert updated["preprocess"]["impute"] == "column_mean"
    # 指定していないキーは前の値のまま(既定ではなく、resolve_request時の明示値)
    assert updated["preprocess"]["normalize"] == "tic"
    assert updated["preprocess"]["min_detection_rate"] == 0.1
    assert updated["value_sources"]["preprocess"]["impute"] == "explicit_update"
    assert updated["value_sources"]["preprocess"]["normalize"] == "explicit"


def test_merge_updates_preprocess_value_must_be_dict(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        merge_updates(req, {"preprocess": "tic"})


def test_merge_updates_preprocess_rejects_unknown_child_key(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        merge_updates(req, {"preprocess": {"unknown_field": 1}})


def test_merge_updates_preprocess_invalid_update_checked_before_mutation(tmp_path):
    """`.update()`呼出しより先に検証することを固定する（brief: 「上記update呼出し
    より先に検証する」）。

    ``merge_updates``自身の事前検査だけが使う文言「未知の更新キー」でmatchする
    ことが肝——これは``_validate_preprocess``（``.update()``適用後の値を検査
    するため、事前検査が後回しになっても最終的に同じDomainError/同じcodeで
    捕まえてしまう）が使う文言「未知のキー」とは異なる。事前検査を``.update()``
    より後ろへ移動すると、この一致は壊れて別の文言でDomainErrorが起きるため、
    「例外型・codeが同じだから合格」という偽陽性を防ぐ。
    あわせて、検証に失敗しても呼び出し側が渡した``request``自体（``out``へ
    deepcopyする前の入力）は一切書き換わらないことも確認する。
    """
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    before = copy.deepcopy(req)

    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        merge_updates(req, {"preprocess": "tic"})
    assert req == before

    with pytest.raises(DomainError, match="未知の更新キー"):
        merge_updates(req, {"preprocess": {"unknown_field": 1}})
    assert req == before


# ---------- comparisons: 方向・重複ID・path脱出 ----------

def test_comparison_reference_and_test_must_differ(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"comparisons": [{
            "comparison_id": "x", "reference_group": "control", "test_group": "control",
        }]})


def test_comparison_duplicate_id_rejected(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"comparisons": [
            {"comparison_id": "a", "reference_group": "control", "test_group": "treated"},
            {"comparison_id": "a", "reference_group": "control", "test_group": "treated2"},
        ]})


@pytest.mark.parametrize("bad_id", ["../escape", "a/b", "a\\b", "..", ""])
def test_comparison_id_path_escape_rejected(tmp_path, bad_id):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"comparisons": [
            {"comparison_id": bad_id, "reference_group": "control", "test_group": "treated"},
        ]})


def test_comparison_defaults_are_filled(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"comparisons": [{
        "comparison_id": "t_vs_c", "reference_group": "control", "test_group": "treated",
    }]})
    comparison = req["comparisons"][0]
    assert comparison["q_threshold"] == 0.05
    assert comparison["log2fc_threshold"] == 1.0
    assert comparison["log_transform"] is True
    assert comparison["allow_confounded"] is False


def test_comparisons_replace_whole_array_on_resume(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"comparisons": [
        {"comparison_id": "a", "reference_group": "control", "test_group": "treated"},
    ]})
    updated = merge_updates(req, {"comparisons": [
        {"comparison_id": "b", "reference_group": "control", "test_group": "treated2"},
    ]})
    ids = [c["comparison_id"] for c in updated["comparisons"]]
    assert ids == ["b"]


# ---------- 未知キー・内部キーの拒否 ----------

def test_resolve_request_rejects_unknown_top_level_key(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"not_a_real_key": 1})


@pytest.mark.parametrize("internal_key", ["effective_target", "value_sources"])
def test_resolve_request_rejects_internal_keys_from_external_input(tmp_path, internal_key):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {internal_key: "x"})


@pytest.mark.parametrize("internal_key", ["effective_target", "value_sources"])
def test_merge_updates_rejects_internal_keys_in_updates(tmp_path, internal_key):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(req, {internal_key: "x"})


def test_merge_updates_rejects_keys_outside_updatable(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root)
    assert UPDATABLE == {"target", "sample_manifest", "preprocess", "comparisons"}
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(req, {"polarity": "positive"})


# ---------- request_fingerprint ----------

def test_request_fingerprint_stable_for_same_content(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req1 = resolve_request(root, {"polarity": "negative"})
    req2 = resolve_request(root, {"polarity": "negative"})
    assert request_fingerprint(req1) == request_fingerprint(req2)


def test_request_fingerprint_differs_for_different_content(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req1 = resolve_request(root, {"polarity": "negative"})
    req2 = resolve_request(root, {"polarity": "positive"})
    assert request_fingerprint(req1) != request_fingerprint(req2)


def test_request_fingerprint_ignores_provenance_only_differences(tmp_path):
    """value_sources/effective_targetのような由来情報だけが違っても内容が同じなら同一hash。"""
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"polarity": "negative"})
    updated = merge_updates(req, {"target": "auto"})  # 値は変わらないが出所はexplicit_updateへ
    assert request_fingerprint(req) == request_fingerprint(updated)


# ---------- validate_request 単体 ----------

def test_validate_request_rejects_unknown_key_directly():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_request({"not_a_real_key": 1})


# ---------- 明示無効化(auto/既定と同値でも出所で区別) ----------

def test_normalize_none_and_drift_correct_false_are_explicit_not_default(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"preprocess": {"normalize": "none", "drift_correct": False}})
    assert req["preprocess"]["normalize"] == "none"
    assert req["preprocess"]["drift_correct"] is False
    assert req["value_sources"]["preprocess"]["normalize"] == "explicit"
    assert req["value_sources"]["preprocess"]["drift_correct"] == "explicit"


def test_min_detection_rate_explicit_zero_distinct_from_default_zero(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    default_req = resolve_request(root)
    explicit_req = resolve_request(root, {"preprocess": {"min_detection_rate": 0.0}})
    assert default_req["preprocess"]["min_detection_rate"] == 0.0
    assert explicit_req["preprocess"]["min_detection_rate"] == 0.0
    assert default_req["value_sources"]["preprocess"]["min_detection_rate"] == "default"
    assert explicit_req["value_sources"]["preprocess"]["min_detection_rate"] == "explicit"


def test_blank_min_fold_explicit_null_is_disable_not_auto(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"preprocess": {"blank_min_fold": None}})
    assert req["preprocess"]["blank_min_fold"] is None
    assert req["value_sources"]["preprocess"]["blank_min_fold"] == "explicit"


def test_max_qc_rsd_explicit_null_is_disable_not_auto(tmp_path):
    """blank_min_foldと対になるもう一方の明示null無効化フィールド（spec l.419）。

    request.pyの``_validate_preprocess``は``for field in ("blank_min_fold",
    "max_qc_rsd")``という共有ループでこの2フィールドを対称に扱っている
    ——blank_min_fold側にしかテストが無いと、max_qc_rsd側だけ壊れる回帰
    （例えばループをblank_min_foldだけの専用分岐へ書き換えて片方だけ
    handlingが漏れる）を拾えない。
    """
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"preprocess": {"max_qc_rsd": None}})
    assert req["preprocess"]["max_qc_rsd"] is None
    assert req["value_sources"]["preprocess"]["max_qc_rsd"] == "explicit"


# ---------- 個別フィールドのenum/型 ----------

def test_polarity_rejects_unknown_value(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"polarity": "both"})


def test_measure_rejects_value_other_than_peak_height(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"measure": "peak_area"})


def test_keep_extension_rejects_empty_string(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"keep_extension": ""})


def test_output_root_accepts_absolute_path_outside_source_root(tmp_path):
    """output_rootは書込先の変更用途があるためsource_root外の絶対パスも許可する。"""
    root = tmp_path / "source"
    root.mkdir()
    outside = str(tmp_path / "elsewhere")
    req = resolve_request(root, {"output_root": outside})
    assert req["output_root"] == outside


def test_method_file_accepts_absolute_path_outside_source_root(tmp_path):
    """spec §10.1の例に合わせ、method_fileはラボ共有フォルダの絶対パスを許可する。"""
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"method_file": "C:/lab/methods/lipid_NEG.txt"})
    assert req["method_file"] == "C:/lab/methods/lipid_NEG.txt"


def test_save_project_rejects_non_bool(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"save_project": "true"})


def test_timeout_s_rejects_non_positive(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"timeout_s": 0})


# ---------- 明示nullの許可範囲（controller裁定R13、spec l.246/421） ----------
#
# 明示nullが許されるのは sample_manifest・preprocess.blank_min_fold・
# preprocess.max_qc_rsd の3つだけ（無効化を意味する設定）。それ以外のキーへ
# 明示nullを渡すのは「未指定」の代用ではなく、常にエラーにする
# （spec l.421「他の不許可なnullはエラーにする」/ l.246「nullが無効化を
# 意味する設定だけで許可する」）。特に comparisons: null を空配列へ黙って
# 変換するのは、spec l.417が戒める「差次的解析の要求を捏造」そのものになる
# ため最も鋭い事例。

@pytest.mark.parametrize("field", [
    "method_file", "lbm_file", "polarity", "keep_extension", "output_root",
])
def test_explicit_null_is_rejected_for_non_disable_top_level_fields(tmp_path, field):
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {field: None})


def test_explicit_null_comparisons_is_rejected_in_resolve_request(tmp_path):
    """comparisons: null を既定の空配列へ黙って読み替えてはいけない。"""
    root = tmp_path / "source"
    root.mkdir()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"comparisons": None})


def test_explicit_null_comparisons_is_rejected_in_merge_updates(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {"comparisons": [
        {"comparison_id": "a", "reference_group": "control", "test_group": "treated"},
    ]})
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        merge_updates(req, {"comparisons": None})


def test_explicit_exploratory_target_with_comparisons_is_rejected(tmp_path):
    """最終レビュー指摘5: `target="exploratory"`＋比較指定は矛盾した要求。

    `engine.build_stages`は`effective_target == "exploratory"`のとき
    `differential:*`/`export:*`を計画から外し、`_mark_out_of_scope`が
    `STAGE_OUT_OF_SCOPE_FOR_TARGET`を**stage側のwarningsにだけ**書く
    ——`read_status`もレポートの「Unverified Conditions」も
    `record["warnings"]`しか読まないので、利用者が明示的に頼んだ比較が
    1件も実行されないまま`completed`が返り、痕跡がどこにも出ない。
    spec §8.2「明示的に要求した処理を実施できなければneeds_input」に反する。

    `target="auto"`＋比較は矛盾ではない（`differential`へ昇格する）ので、
    拒否されるのは**明示的なexploratory**との組合せだけ。
    """
    root = tmp_path / "source"
    root.mkdir()
    comparisons = [{"comparison_id": "treated_vs_control",
                    "reference_group": "control", "test_group": "treated"}]

    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root, {"target": "exploratory", "comparisons": comparisons})

    # autoは昇格するので通る（拒否の的が広すぎないことの確認）。
    promoted = resolve_request(root, {"target": "auto", "comparisons": comparisons})
    assert promoted["effective_target"] == "differential"

    # resumeで後から比較だけを足す経路も同じ理由で拒否する。
    exploratory = resolve_request(root, {"target": "exploratory"})
    assert exploratory["effective_target"] == "exploratory"
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        merge_updates(exploratory, {"comparisons": comparisons})

    # 同じresumeでtargetも一緒に直せば通る（利用者に出口がある）。
    fixed = merge_updates(exploratory, {"target": "differential", "comparisons": comparisons})
    assert fixed["effective_target"] == "differential"


def test_permitted_explicit_nulls_still_work(tmp_path):
    """R13の許可リスト3つ（sample_manifest・blank_min_fold・max_qc_rsd）は
    引き続き明示nullで無効化できる。"""
    root = tmp_path / "source"
    root.mkdir()
    req = resolve_request(root, {
        "sample_manifest": None,
        "preprocess": {"blank_min_fold": None, "max_qc_rsd": None},
    })
    assert req["sample_manifest"] is None
    assert req["value_sources"]["sample_manifest"] == "explicit"
    assert req["preprocess"]["blank_min_fold"] is None
    assert req["preprocess"]["max_qc_rsd"] is None
    assert req["value_sources"]["preprocess"]["blank_min_fold"] == "explicit"
    assert req["value_sources"]["preprocess"]["max_qc_rsd"] == "explicit"

    released = merge_updates(req, {"sample_manifest": "sheet.tsv"})
    released_again = merge_updates(released, {"sample_manifest": None})
    assert released_again["sample_manifest"] is None
    assert released_again["value_sources"]["sample_manifest"] == "explicit_update"


# ---------- analysis-request.json の探索と3段の優先順位（spec §7.1） ----------

def _source_with_request_file(tmp_path, payload) -> "Path":
    import json
    from pathlib import Path

    from metabolomix.pipeline.request import REQUEST_FILE_NAME

    root = tmp_path / "source"
    root.mkdir(exist_ok=True)
    text = payload if isinstance(payload, str) else json.dumps(payload)
    (root / REQUEST_FILE_NAME).write_text(text, encoding="utf-8")
    return Path(root)


def test_the_request_file_is_used_when_nothing_is_specified(tmp_path):
    """元フォルダ直下に置いた要求ファイルが、既定値より優先される。"""
    root = _source_with_request_file(tmp_path, {
        "polarity": "positive", "timeout_s": 60,
        "preprocess": {"normalize": "pqn", "impute": "knn"},
        "comparisons": [{"comparison_id": "cmp1", "reference_group": "control",
                         "test_group": "treated"}]})

    request = resolve_request(root)

    assert request["polarity"] == "positive"
    assert request["timeout_s"] == 60
    assert request["preprocess"]["normalize"] == "pqn"
    assert request["effective_target"] == "differential"   # 比較があるので自動判定も効く
    assert request["value_sources"]["polarity"] == "request_file"
    assert request["value_sources"]["preprocess"]["normalize"] == "request_file"
    assert request["value_sources"]["measure"] == "default"


def test_an_explicit_value_beats_the_request_file(tmp_path):
    """MCPで明示した値 > analysis-request.json（spec §7.1の優先順位）。"""
    root = _source_with_request_file(tmp_path, {"polarity": "positive", "timeout_s": 60})

    request = resolve_request(root, {"polarity": "negative"})

    assert request["polarity"] == "negative"
    assert request["timeout_s"] == 60           # 明示していない側はファイルのまま
    assert request["value_sources"]["polarity"] == "explicit"
    assert request["value_sources"]["timeout_s"] == "request_file"


def test_preprocess_is_layered_per_sub_key(tmp_path):
    """preprocessは子キー単位で重ねる（片方が丸ごと消えない）。"""
    root = _source_with_request_file(
        tmp_path, {"preprocess": {"normalize": "pqn", "impute": "knn"}})

    request = resolve_request(root, {"preprocess": {"impute": "none"}})

    assert request["preprocess"]["normalize"] == "pqn"     # ファイル由来が生き残る
    assert request["preprocess"]["impute"] == "none"       # 明示が勝つ
    assert request["value_sources"]["preprocess"]["normalize"] == "request_file"
    assert request["value_sources"]["preprocess"]["impute"] == "explicit"
    assert request["value_sources"]["preprocess"]["max_qc_rsd"] == "default"


def test_a_request_file_value_participates_in_the_fingerprint(tmp_path):
    """ファイル由来でも最終的な値が違えば別の要求（受付の冪等性判定の土台）。"""
    with_file = resolve_request(_source_with_request_file(tmp_path, {"timeout_s": 60}))
    plain_root = tmp_path / "plain"
    plain_root.mkdir()
    assert request_fingerprint(with_file) != request_fingerprint(resolve_request(plain_root))


def test_a_broken_request_file_is_rejected_not_ignored(tmp_path):
    """壊れたJSONを黙って無視すると、利用者は自分の指定で走ったと思い込む。"""
    root = _source_with_request_file(tmp_path, "{ not json")
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root)


def test_an_unknown_key_in_the_request_file_is_rejected(tmp_path):
    root = _source_with_request_file(tmp_path, {"targett": "exploratory"})
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root)


def test_a_non_object_request_file_is_rejected(tmp_path):
    root = _source_with_request_file(tmp_path, ["exploratory"])
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root)


def test_a_disallowed_null_in_the_request_file_is_rejected(tmp_path):
    """nullは「無効化」を意味する設定だけの語彙（明示値と同じ規則を課す）。"""
    root = _source_with_request_file(tmp_path, {"method_file": None})
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(root)


def test_no_request_file_keeps_the_previous_defaults(tmp_path):
    """ファイルが無ければ従来どおり既定値だけで解決する。"""
    root = tmp_path / "source"
    root.mkdir()
    request = resolve_request(root)
    assert request["timeout_s"] == 21600
    assert request["value_sources"]["timeout_s"] == "default"


# ---------- library_mode / msp_file（spec 2026-10-07） ----------

def _root(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    return root


def test_library_mode_msp_only_is_accepted_with_an_msp_file(tmp_path):
    req = resolve_request(_root(tmp_path), {"library_mode": "msp_only", "msp_file": "lab.msp"})
    assert req["library_mode"] == "msp_only"
    assert req["msp_file"] == "lab.msp"
    assert req["value_sources"]["library_mode"] == "explicit"


@pytest.mark.parametrize("explicit", [
    {"library_mode": "lbm_only"},
    {"library_mode": None},
    {"msp_file": None},
    {"library_mode": "msp_only", "msp_file": ""},
    {"library_mode": "msp_only", "lbm_file": "x.lbm2"},
    {"msp_file": "lab.msp"},  # auto では使われないので黙って受けない
])
def test_invalid_library_mode_combinations_are_rejected(tmp_path, explicit):
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(_root(tmp_path), explicit)


def test_library_mode_cannot_be_changed_on_resume(tmp_path):
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(resolve_request(_root(tmp_path)), {"library_mode": "msp_only"})


def test_default_library_mode_does_not_change_the_request_fingerprint(tmp_path):
    from metabolomix.core.atomic_io import canonical_hash
    from metabolomix.pipeline import request as request_mod
    req = resolve_request(_root(tmp_path))
    legacy_keys = request_mod._TOP_LEVEL_KEYS - {"library_mode", "msp_file"}
    assert request_fingerprint(req) == canonical_hash({k: req[k] for k in legacy_keys if k in req})


def test_msp_only_changes_the_request_fingerprint(tmp_path):
    root = _root(tmp_path)
    assert (request_fingerprint(resolve_request(root, {"library_mode": "msp_only"}))
            != request_fingerprint(resolve_request(root)))


def test_v2_rejects_library_mode_with_a_migration_hint(tmp_path):
    from metabolomix.pipeline import request_v2
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID") as info:
        request_v2._check_known_keys({"library_mode": "msp_only"})
    assert "library_mode" in str(info.value)


def test_resume_of_a_request_saved_before_library_mode_still_validates(tmp_path):
    """最終レビュー #1: 変更前に保存された要求（library_mode / msp_file キー無し）でも
    resume の updates が通り、内容 hash も変わらない。"""
    req = resolve_request(_root(tmp_path))
    legacy = {k: v for k, v in req.items() if k not in ("library_mode", "msp_file")}
    legacy["value_sources"] = {k: v for k, v in req["value_sources"].items()
                               if k not in ("library_mode", "msp_file")}
    assert request_fingerprint(legacy) == request_fingerprint(req)
    updated = merge_updates(legacy, {"target": "exploratory"})
    assert updated["target"] == "exploratory"
