"""conservative-v1: 自動前処理レシピの解決と適用状況の監査（spec §8）。

自動前処理は「証拠が無いのに補正した」と「補正したと言いながら実は補正していない」の
両方を防がなければならない。ここで確認するのは3つ:

1. `resolve_policy` — include=true のメタデータ（confidence を含む）とQC健全性から、
   何を適用してよいかを判断する。未確認の batch/injection_order/qc_pool は
   auto のドリフト補正・正規化・QC-RSDフィルタを一切有効化しない。
2. `check_applied_policy` — resolve_policy の判断が、実際の数値層の実行結果
   （report）と一致しているかを再照合する。数値層の caveat だけに委ねない。
3. `preprocess_auto` — 上記2つとTask 6の `preprocess_dataset` を結び、
   `applied_steps` は計算が完了するまで空のままにする。
"""
from __future__ import annotations

import pytest

from metabolomix.analysis.dataset_service import preprocess_auto
from metabolomix.analysis.preprocess_policy import (
    POLICY_VERSION,
    check_applied_policy,
    resolve_policy,
)
from metabolomix.core.atomic_io import DomainError
from tests.pipeline_fixtures import make_dataset, metadata_rows


def _with_overrides(rows, overrides_by_index):
    """metadata_rows() の合成行を、テストケースごとに部分的に上書きする。

    `overrides_by_index` は `{index: {field: value}}`。値・provenance の両方を
    書き換える必要がある場合は `{index: {field: (value, source, confidence)}}`
    の3要素タプルで渡す。
    """
    rows = [dict(r, provenance={k: dict(v) for k, v in r["provenance"].items()})
            for r in rows]
    for index, fields in overrides_by_index.items():
        row = rows[index]
        for field, value in fields.items():
            if isinstance(value, tuple):
                actual, source, confidence = value
                row[field] = actual
                row["provenance"][field] = {
                    "value": actual, "source": source, "confidence": confidence,
                }
            else:
                row[field] = value
                row["provenance"][field]["value"] = value
    return rows


# ---------- 代表ケース（brief step 1） ----------

def test_unverified_order_never_enables_auto_drift():
    ds = make_dataset()
    rows = metadata_rows(ds, confirmed=False)
    plan = resolve_policy(ds, {"policy": "conservative-v1", "drift_correct": "auto"}, rows)
    assert plan["resolved_recipe"]["drift_correct"] is False
    assert "drift_correct" in plan["skipped_steps"]


def test_confirmed_pool_with_enclosing_qc_enables_drift():
    ds = make_dataset()
    plan = resolve_policy(ds, {"policy": "conservative-v1"}, metadata_rows(ds))
    assert plan["resolved_recipe"]["normalize"] == "pqn"
    assert plan["resolved_recipe"]["drift_correct"] is True
    assert plan["resolved_recipe"]["max_qc_rsd"] == 0.30


# ---------- resolve_policy: 追加ケース（brief step 4） ----------

def test_plan_has_all_required_keys():
    ds = make_dataset()
    plan = resolve_policy(ds, {"policy": "conservative-v1"}, metadata_rows(ds))
    for key in ("requested_recipe", "resolved_recipe", "applied_steps",
               "skipped_steps", "reasons", "policy_version", "assumptions"):
        assert key in plan
    assert plan["applied_steps"] == []  # 計算完了まで空のまま
    assert plan["policy_version"] == POLICY_VERSION


@pytest.mark.parametrize("n_qc", [0, 2, 3, 4])
def test_qc_count_thresholds_gate_normalize_and_rsd(n_qc):
    ds = make_dataset()
    plan = resolve_policy(ds, {}, metadata_rows(ds, n_qc=n_qc))
    resolved = plan["resolved_recipe"]
    if n_qc >= 3:
        assert resolved["normalize"] == "pqn"
        assert resolved["max_qc_rsd"] == 0.30
    else:
        assert resolved["normalize"] == "none"
        assert resolved["max_qc_rsd"] is None
        assert "normalize" in plan["skipped_steps"]
        assert "qc_rsd_filter" in plan["skipped_steps"]


@pytest.mark.parametrize("n_qc", [3, 4])
def test_drift_correct_needs_four_healthy_qc(n_qc):
    ds = make_dataset()
    plan = resolve_policy(ds, {}, metadata_rows(ds, n_qc=n_qc))
    resolved = plan["resolved_recipe"]
    if n_qc >= 4:
        assert resolved["drift_correct"] is True
    else:
        assert resolved["drift_correct"] is False
        assert "drift_correct" in plan["skipped_steps"]


def test_unknown_batch_disables_auto_qc_steps():
    ds = make_dataset()
    rows = metadata_rows(ds)
    # 全件の batch を「未確認」に落とす（値そのものはあるが provenance が unverified）。
    rows = _with_overrides(rows, {i: {"batch": ("B1", "mztab", "unverified")}
                                  for i in range(8)})
    plan = resolve_policy(ds, {}, rows)
    resolved = plan["resolved_recipe"]
    assert resolved["normalize"] == "none"
    assert resolved["drift_correct"] is False
    assert resolved["max_qc_rsd"] is None
    assert plan["assumptions"]["batch_ok"] is False


def test_multi_batch_disables_auto_qc_steps():
    ds = make_dataset()
    rows = metadata_rows(ds)
    # QCの半分をB2へ動かす（単一confirmed batchでなくなる）。
    rows = _with_overrides(rows, {0: {"batch": ("B2", "user_manifest", "confirmed")}})
    plan = resolve_policy(ds, {}, rows)
    resolved = plan["resolved_recipe"]
    assert resolved["normalize"] == "none"
    assert resolved["drift_correct"] is False
    assert plan["assumptions"]["batch_ok"] is False


def test_two_qc_pools_disables_auto_and_rejects_explicit_pqn():
    ds = make_dataset()
    rows = metadata_rows(ds)
    # QC 4件のうち2件をpool2へ（単一poolでなくなる）。
    rows = _with_overrides(rows, {
        0: {"qc_pool": ("pool2", "user_manifest", "confirmed")},
        1: {"qc_pool": ("pool2", "user_manifest", "confirmed")},
    })
    plan = resolve_policy(ds, {}, rows)
    assert plan["resolved_recipe"]["normalize"] == "none"
    assert plan["assumptions"]["same_pool"] is False

    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"normalize": "pqn"}, rows)
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"


def test_explicit_pqn_with_zero_qc_is_prerequisite_missing():
    """参照に使える健全なQCが1件も無い明示PQNは実施不能として止める(spec §8.2)。

    pool_ambiguous には触れない条件（QCが1件も無ければpoolも定義しようがない）
    で、n_qc==0 だけを理由にPREPREQUISITE_MISSINGになることを確認する。
    """
    ds = make_dataset()
    rows = metadata_rows(ds, n_qc=0)
    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"normalize": "pqn"}, rows)
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"


@pytest.mark.parametrize("n_qc", [1, 2])
def test_explicit_pqn_with_insufficient_qc_proceeds_and_records_shortfall(n_qc):
    """auto閾値(3件)未満のQCでも明示PQNは実施するが、不足を記録する(spec §8.1/8.2)。

    「実施しない」ではなく「実施した上でQC補正済みと誤認させない」ことが論点
    なので、raiseしないこと・resolved_recipeがpqnのままであること・assumptions
    に不足件数と対象sample_idが載ることの3点を確認する。
    """
    ds = make_dataset()
    rows = metadata_rows(ds, n_qc=n_qc)
    plan = resolve_policy(ds, {"normalize": "pqn"}, rows)
    assert plan["resolved_recipe"]["normalize"] == "pqn"
    weak = plan["assumptions"]["normalize_weak_qc_reference"]
    assert weak["n_qc"] == n_qc
    assert weak["sample_ids"] == [f"S{i}" for i in range(n_qc)]
    assert "normalize" in plan["reasons"]


def test_explicit_pqn_with_failed_qc_proceeds_and_names_affected_samples():
    """健全なQCが3件以上あっても、失敗疑いのQCを含んだままの明示PQNは記録する。

    n_qc=4(auto閾値以上)なので「件数不足」ではなく「健全性」だけが理由になる
    ケースを、既存test_failed_qc_disables_auto_qc_stepsと同じ壊し方(S0の総強度を
    1/10以下に落とす)で作る。対象sample_id(S0)がassumptionsへ載ることを確認する
    —— これが無いと品質レポートは「QC補正済み」と誤認させ得る(spec §8.2)。
    """
    ds = make_dataset()
    ds.feature_matrix = ds.feature_matrix.copy()
    ds.feature_matrix[:, 0] *= 0.01
    rows = metadata_rows(ds)  # n_qc既定4 → 件数条件は満たす
    plan = resolve_policy(ds, {"normalize": "pqn"}, rows)
    assert plan["resolved_recipe"]["normalize"] == "pqn"
    weak = plan["assumptions"]["normalize_weak_qc_reference"]
    assert weak["n_qc"] == 4  # 件数自体はauto閾値を満たしている
    assert "S0" in weak["failed_qc"]
    assert "normalize" in plan["reasons"]


def test_failed_qc_disables_auto_qc_steps():
    ds = make_dataset()
    # QC(S0)の総強度をQC中央値の1/10以下へ落とす(detect_failed_qcが失敗判定する)。
    ds.feature_matrix = ds.feature_matrix.copy()
    ds.feature_matrix[:, 0] *= 0.01
    plan = resolve_policy(ds, {}, metadata_rows(ds))
    assert plan["resolved_recipe"]["normalize"] == "none"
    assert plan["resolved_recipe"]["drift_correct"] is False
    assert "S0" in plan["assumptions"]["failed_qc"]


def test_sample_outside_qc_span_disables_drift_but_not_normalize():
    ds = make_dataset()
    rows = metadata_rows(ds)
    # QCの注入順区間は[1,8]。1件のsampleを区間外(order=99)へ動かす。
    outside_index = 4  # role="sample" (n_qc既定4なのでindex4以降がsample)
    assert rows[outside_index]["role"] == "sample"
    rows = _with_overrides(rows, {outside_index: {"injection_order": 99}})
    rows[outside_index]["provenance"]["injection_order"] = {
        "value": 99, "source": "user_manifest", "confidence": "confirmed",
    }
    plan = resolve_policy(ds, {}, rows)
    assert plan["resolved_recipe"]["drift_correct"] is False
    assert plan["assumptions"]["all_samples_enclosed"] is False
    # 正規化・RSDフィルタはQC区間外サンプルの影響を受けない。
    assert plan["resolved_recipe"]["normalize"] == "pqn"


def test_no_blank_leaves_blank_filter_disabled():
    ds = make_dataset()
    plan = resolve_policy(ds, {}, metadata_rows(ds))  # フィクスチャにblankは無い
    assert plan["resolved_recipe"]["blank_min_fold"] is None
    assert "blank_filter" in plan["skipped_steps"]


def test_blank_and_sample_present_enables_blank_filter():
    ds = make_dataset()
    rows = metadata_rows(ds)
    # role="sample"の1件をblankへ差し替える(blank/sampleが両方揃う)。
    rows = _with_overrides(rows, {7: {"role": "blank"}})
    plan = resolve_policy(ds, {}, rows)
    assert plan["resolved_recipe"]["blank_min_fold"] == 3.0
    assert "blank_filter" not in plan["skipped_steps"]


def test_explicit_false_and_null_disable_without_error():
    ds = make_dataset()
    requested = {
        "normalize": "none", "blank_min_fold": None,
        "drift_correct": False, "max_qc_rsd": None,
    }
    plan = resolve_policy(ds, requested, metadata_rows(ds))
    resolved = plan["resolved_recipe"]
    assert resolved["normalize"] == "none"
    assert resolved["blank_min_fold"] is None
    assert resolved["drift_correct"] is False
    assert resolved["max_qc_rsd"] is None
    # 明示的な無効化は「条件不足のskip」ではない。
    assert plan["skipped_steps"] == []
    assert plan["reasons"] == {}


def test_explicit_drift_without_prerequisite_is_prerequisite_missing():
    ds = make_dataset()
    rows = metadata_rows(ds, confirmed=False)
    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"drift_correct": True}, rows)
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"


def test_explicit_max_qc_rsd_without_enough_qc_is_prerequisite_missing():
    ds = make_dataset()
    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"max_qc_rsd": 0.30}, metadata_rows(ds, n_qc=2))
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"


def test_explicit_blank_min_fold_without_blank_is_prerequisite_missing():
    ds = make_dataset()
    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"blank_min_fold": 3.0}, metadata_rows(ds))
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"


def test_no_detection_mask_with_explicit_min_detection_rate_is_prerequisite_missing():
    ds = make_dataset()
    assert getattr(ds, "detected_mask", None) is None
    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"min_detection_rate": 0.2}, metadata_rows(ds))
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"


def test_unknown_policy_version_is_rejected():
    ds = make_dataset()
    with pytest.raises(DomainError) as exc:
        resolve_policy(ds, {"policy": "conservative-v2"}, metadata_rows(ds))
    assert exc.value.code == "PREPROCESS_POLICY_UNKNOWN"


# ---------- check_applied_policy ----------

def test_check_applied_policy_passes_when_report_matches_plan():
    plan = {"resolved_recipe": {"normalize": "pqn", "blank_min_fold": None,
                                "drift_correct": False, "max_qc_rsd": 0.30}}
    report = {
        "recipe_applied": ["normalize", "qc_rsd_filter", "impute"],
        "recipe_skipped": [],
        "steps": {"normalize": {"method": "pqn", "n_samples": 8}},
        "features_after": 5, "features_before": 6,
    }
    check_applied_policy(plan, report)  # 例外が出なければ良い


def test_check_applied_policy_rejects_step_claimed_but_not_applied():
    """resolved_recipeがdrift_correctを要求したのに、reportでは未適用（数値層都合）。"""
    plan = {"resolved_recipe": {"normalize": "none", "blank_min_fold": None,
                                "drift_correct": True, "max_qc_rsd": None}}
    report = {
        "recipe_applied": ["impute"],
        "recipe_skipped": ["drift_correct"],
        "steps": {"drift_correct": {"status": "skipped"}},
        "features_after": 6, "features_before": 6,
    }
    with pytest.raises(DomainError) as exc:
        check_applied_policy(plan, report)
    assert exc.value.code == "PREPROCESS_STEP_NOT_APPLIED"


def test_check_applied_policy_rejects_degenerate_normalization():
    plan = {"resolved_recipe": {"normalize": "pqn", "blank_min_fold": None,
                                "drift_correct": False, "max_qc_rsd": None}}
    report = {
        "recipe_applied": ["normalize", "impute"],
        "recipe_skipped": [],
        "steps": {"normalize": {"method": "pqn", "unscaled_samples": 2, "n_samples": 8}},
        "features_after": 6, "features_before": 6,
    }
    with pytest.raises(DomainError) as exc:
        check_applied_policy(plan, report)
    assert exc.value.code == "NORMALIZATION_DEGENERATE"


def test_check_applied_policy_rejects_zero_features_remaining():
    plan = {"resolved_recipe": {"normalize": "none", "blank_min_fold": None,
                                "drift_correct": False, "max_qc_rsd": 0.0}}
    report = {
        "recipe_applied": ["qc_rsd_filter", "impute"],
        "recipe_skipped": [],
        "steps": {},
        "features_after": 0, "features_before": 6,
    }
    with pytest.raises(DomainError) as exc:
        check_applied_policy(plan, report)
    assert exc.value.code == "PREPROCESS_FEATURES_EXHAUSTED"


# ---------- preprocess_auto ----------

def test_preprocess_auto_commits_and_reports_applied_steps():
    ds = make_dataset()
    plan = preprocess_auto(ds, {"policy": "conservative-v1"}, metadata_rows(ds))
    assert plan["resolved_recipe"]["normalize"] == "pqn"
    assert "normalize" in plan["applied_steps"]
    assert "impute" in plan["applied_steps"]
    assert ds.pp_matrix is not None
    assert ds.preprocessing_recipe["normalize"] == "pqn"


def test_preprocess_auto_raises_and_leaves_dataset_untouched_on_prerequisite_missing():
    ds = make_dataset()
    with pytest.raises(DomainError) as exc:
        preprocess_auto(ds, {"max_qc_rsd": 0.30}, metadata_rows(ds, n_qc=2))
    assert exc.value.code == "PREPROCESS_PREREQUISITE_MISSING"
    assert ds.pp_matrix is None
    assert ds.preprocess_id is None


def test_preprocess_auto_raises_on_degenerate_normalization_and_stays_uncommitted():
    """正規化係数が0/非有限の試料が混じる場合、実データ経由でも自動進行しない。

    check_applied_policyの単体テストは合成reportを使うが、こちらは
    resolve_policy→run_dataset_preprocess→check_applied_policyの一連が
    実際の前処理経路で本当に噛み合っているかを確認する（`ds.feature_matrix`の
    1試料を全特徴ゼロにし、PQNの正規化係数を退化させる）。
    """
    ds = make_dataset()
    ds.feature_matrix = ds.feature_matrix.copy()
    ds.feature_matrix[:, 4] = 0.0  # S4(role=sample) を全特徴ゼロにする
    with pytest.raises(DomainError) as exc:
        preprocess_auto(ds, {}, metadata_rows(ds))
    assert exc.value.code == "NORMALIZATION_DEGENERATE"
    # 検査で落ちた場合、Task 6のcommitへは到達していない。
    assert ds.pp_matrix is None
    assert ds.preprocess_id is None



def test_preprocess_auto_computes_the_matrix_only_once(monkeypatch):
    """検査用とcommit用で同じ行列を2回計算しない（レビュー指摘P3）。

    `run_dataset_preprocess`の呼出しの間に``ds``の入力は何も変わらないので、
    2度目は結果まで同じ——検体数×特徴量数の全量に対する重複計算で、
    得られるものが無い。
    """
    from metabolomix.analysis import dataset_service
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess

    calls = {"count": 0}

    def _counting(ds_arg, recipe):
        calls["count"] += 1
        return run_dataset_preprocess(ds_arg, recipe)

    monkeypatch.setattr(dataset_service, "run_dataset_preprocess", _counting)

    ds = make_dataset()
    plan = preprocess_auto(ds, {"policy": "conservative-v1"}, metadata_rows(ds))

    assert calls["count"] == 1
    # 1回に畳んでも、commitされた結果と検査は従来どおり成立している。
    assert ds.pp_matrix is not None
    assert plan["result"]["provenance"]["result_id"] == ds.preprocess_id
    assert plan["applied_steps"] == list(plan["result"].get("recipe_applied", []))


def test_preprocess_auto_reuses_an_identical_preprocessing_without_recomputing(monkeypatch):
    """同じ入力・同じレシピの再実行では、計算そのものを1回も行わない。"""
    from metabolomix.analysis import dataset_service
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess

    ds = make_dataset()
    rows = metadata_rows(ds)
    first = preprocess_auto(ds, {"policy": "conservative-v1"}, rows)

    calls = {"count": 0}

    def _counting(ds_arg, recipe):
        calls["count"] += 1
        return run_dataset_preprocess(ds_arg, recipe)

    monkeypatch.setattr(dataset_service, "run_dataset_preprocess", _counting)
    second = preprocess_auto(ds, {"policy": "conservative-v1"}, rows)

    assert calls["count"] == 0
    assert (second["result"]["provenance"]["result_id"]
            == first["result"]["provenance"]["result_id"])


def test_dataset_preprocess_default_arguments_are_unchanged():
    """通常のdataset_preprocess(単体ツール経路)の既定引数はpipeline由来のautoで変えない。"""
    import inspect

    from metabolomix.analysis.dataset_service import preprocess_dataset

    sig = inspect.signature(preprocess_dataset)
    assert list(sig.parameters) == ["ds", "recipe", "request_revision"]
    assert sig.parameters["request_revision"].default is None
