import numpy as np
import pytest

from metabolomix.mztab.dataset_state import DatasetState


def _make_ds(n_features=20, n_samples=8, with_blank=False):
    """最小限の DatasetState。ロールは dataset_analysis がサンプル名から推定する。"""
    rng = np.random.default_rng(42)
    ds = DatasetState()
    names = [f"ctrl_{i}" for i in range(n_samples // 2)]
    names += [f"treat_{i}" for i in range(n_samples - len(names))]
    if with_blank:
        names[-1] = "blank_1"
    ds.feature_matrix = rng.random((n_features, n_samples)) * 1000.0
    ds.sample_names = names
    ds.feature_ids = [f"f{i}" for i in range(n_features)]
    return ds


def _ds_with_assay_meta(names, assay_meta):
    """sample_names と assay_metadata を持つ最小の DatasetState を作る。"""
    ds = DatasetState()
    ds.feature_matrix = np.arange(3 * len(names), dtype=float).reshape(3, len(names)) + 1.0
    ds.sample_names = list(names)
    ds.feature_ids = ["1", "2", "3"]
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(len(names))]
    ds.assay_metadata = assay_meta
    return ds


def _preprocessed(ds, recipe=None):
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess
    (ds.pp_matrix, ds.pp_sample_names, ds.pp_feature_names,
     ds.roles, ds.sample_meta, report) = run_dataset_preprocess(ds, recipe or {})
    return report


# ---------- build_dataset_pp_inputs ----------

def test_build_dataset_pp_inputs_transposes():
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _make_ds(n_features=20, n_samples=8)
    matrix, sample_names, feature_names, roles, sample_meta = build_dataset_pp_inputs(ds)
    assert matrix.shape == (8, 20)          # (n_samples, n_features)
    assert sample_names == ds.sample_names
    assert feature_names == ds.feature_ids
    assert set(roles) == set(ds.sample_names)
    assert set(sample_meta) == set(ds.sample_names)


def test_build_dataset_pp_inputs_detects_blank_role():
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _make_ds(with_blank=True)
    _, _, _, roles, _ = build_dataset_pp_inputs(ds)
    assert roles["blank_1"] == "blank"


def test_build_dataset_pp_inputs_uses_mztab_injection_order():
    """assay 順序を全件 None にすると QC ドリフト補正が常に skipped になる。"""
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _ds_with_assay_meta(["s1", "s2"], {
        "assay[1]": {"run_order": 7, "batch": "1"},
        "assay[2]": {"run_order": 3, "batch": "1"},
    })
    _, _, _, _, sample_meta = build_dataset_pp_inputs(ds)
    assert sample_meta["s1"]["run_order"] == 7
    assert sample_meta["s2"]["run_order"] == 3
    assert sample_meta["s1"]["run_order_source"] == "mztab_injection_sequence"


def test_build_dataset_pp_inputs_prefers_filename_date_when_batch_is_constant():
    """既定の定数バッチを採ると日付由来の交絡を見落とす。"""
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _ds_with_assay_meta(["20220901_a", "20220902_b"], {
        "assay[1]": {"run_order": 1, "batch": "1"},
        "assay[2]": {"run_order": 2, "batch": "1"},
    })
    _, _, _, _, sample_meta = build_dataset_pp_inputs(ds)
    assert sample_meta["20220901_a"]["batch"] == "20220901"
    assert sample_meta["20220901_a"]["batch_source"] == "filename_date"


def test_build_dataset_pp_inputs_uses_mztab_batch_when_it_varies():
    """変化する assay バッチを捨てると真の実験バッチを使えない。"""
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _ds_with_assay_meta(["20220901_a", "20220901_b"], {
        "assay[1]": {"run_order": 1, "batch": "B1"},
        "assay[2]": {"run_order": 2, "batch": "B2"},
    })
    _, _, _, _, sample_meta = build_dataset_pp_inputs(ds)
    assert sample_meta["20220901_a"]["batch"] == "B1"
    assert sample_meta["20220901_b"]["batch"] == "B2"
    assert sample_meta["20220901_a"]["batch_source"] == "mztab_batch_label"


def test_build_dataset_pp_inputs_rejects_empty_matrix():
    from metabolomix.analysis.dataset_analysis import (
        PreconditionError, build_dataset_pp_inputs,
    )
    ds = DatasetState()
    with pytest.raises(PreconditionError):
        build_dataset_pp_inputs(ds)


# ---------- run_dataset_preprocess ----------

def test_run_dataset_preprocess_keeps_feature_count():
    ds = _make_ds(n_features=20, n_samples=8)
    report = _preprocessed(ds)
    assert ds.pp_matrix.shape == (8, 20)
    assert len(ds.pp_feature_names) == 20
    assert report["features_before"] == 20


def test_run_dataset_preprocess_drops_blank_samples():
    """ブランクは背景除去の参照として使った後、解析行列から外す（ARF と同じ）。"""
    ds = _make_ds(n_features=20, n_samples=8, with_blank=True)
    report = _preprocessed(ds)
    assert "blank_1" not in ds.pp_sample_names
    assert report["excluded_from_matrix"]["blank"] == ["blank_1"]


def test_run_dataset_preprocess_warns_no_run_order():
    """注入順を持たない mzTab では drift_correct は必ず skipped になる。"""
    ds = _ds_with_assay_meta(["s1", "s2"], {"assay[1]": {}, "assay[2]": {}})
    report = _preprocessed(ds, {"drift_correct": True})
    assert report["steps"]["drift_correct"]["status"] == "skipped"
    assert any("注入順" in c for c in report["caveats"])


def test_run_dataset_preprocess_caveats_injection_order_provenance():
    """注入順を読めたときに「未読」と誤報すると補正の妥当性確認を誤らせる。"""
    ds = _ds_with_assay_meta(["s1", "s2"], {
        "assay[1]": {"run_order": 1, "batch": "1"},
        "assay[2]": {"run_order": 2, "batch": "1"},
    })
    report = _preprocessed(ds, {
        "normalize": "none", "drift_correct": True, "impute": "half_min",
    })
    caveats = " ".join(report["caveats"])
    assert "読み取っていない" not in caveats
    assert "ファイル読み込み順" in caveats


def _drift_fixture_dataset():
    """移動中央値を手計算できる、QC 4件と試料2件の単一特徴量。"""
    names = ["QC_1", "QC_2", "QC_3", "QC_4", "sample_A", "sample_B"]
    values = [1.0, 3.0, 5.0, 7.0, 2.0, 4.0]
    orders = [1, 3, 5, 7, 2, 4]
    ds = _ds_with_assay_meta(names, {
        f"assay[{i + 1}]": {"run_order": order, "batch": "1"}
        for i, order in enumerate(orders)
    })
    ds.feature_matrix = np.asarray([values], dtype=float)
    ds.feature_ids = ["f1"]
    return ds, names, values, orders


def test_run_dataset_preprocess_applies_mztab_qc_drift_correction():
    """run_order 配線を外すと applied にならず、手計算済み補正値にも到達しない。"""
    ds, names, _, _ = _drift_fixture_dataset()
    matrix, pp_names, _, _, _, report = _run_dataset_preprocess_for_test(ds)

    # QC [1, 3, 5, 7] の window=3 移動中央値は [2, 3, 5, 6]、中央値は 4。
    # 実装は window=5 を QC数4へ縮小後に奇数化するため、実際の trend は [2, 3, 5, 6]。
    # したがって各注入順の補正値はこのリテラルになる。
    assert pp_names == names
    np.testing.assert_allclose(matrix[:, 0], [2.0, 4.0, 4.0, 14.0 / 3.0, 3.2, 4.0])
    assert report["steps"]["drift_correct"]["status"] == "applied"
    assert not np.allclose(matrix[:, 0], [1.0, 3.0, 5.0, 7.0, 2.0, 4.0])


def _run_dataset_preprocess_for_test(ds):
    """テスト中に DatasetState へ副作用を書かず、純アダプタの戻り値を得る。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess
    return run_dataset_preprocess(ds, {
        "normalize": "none", "drift_correct": True, "impute": "half_min",
    })


def test_run_dataset_preprocess_matches_arf_preprocess_for_qc_drift(monkeypatch):
    """DatasetState 側だけで補正式を複製すると、ARF 経路との数値乖離を見逃す。"""
    import json

    from metabolomix.arf.tools import arf_preprocess
    from metabolomix.core import session_state
    from metabolomix.core.session_state import AnalysisSession

    ds, names, values, orders = _drift_fixture_dataset()
    dataset_matrix, dataset_names, _, _, _, dataset_report = _run_dataset_preprocess_for_test(ds)

    monkeypatch.setattr(session_state, "session", AnalysisSession())
    rows = []
    for i, (name, value) in enumerate(zip(names, values)):
        row = [0] * 40
        row[1] = name
        row[2] = 101
        row[18] = value
        rows.append(row)
    session_state.session.arf.filtered_features = [{
        "MasterAlignmentID": 0,
        "Name": "Lipid_0",
        "MassCenter": 700.0,
        "RT": 5.0,
        "AlignedPeakProperties": rows,
    }]
    session_state.session.arf.class_index = {"records": [
        {"file_name": name, "class_id": "QC" if name.startswith("QC_") else "sample",
         "analytical_order": order}
        for name, order in zip(names, orders)
    ]}

    arf_report = json.loads(arf_preprocess(
        normalize="none", drift_correct=True, impute="half_min",
    ))
    arf_by_name = dict(zip(session_state.session.arf.pp_sample_names,
                           session_state.session.arf.feature_matrix[:, 0]))
    dataset_by_name = dict(zip(dataset_names, dataset_matrix[:, 0]))

    assert dataset_report["steps"]["drift_correct"]["status"] == "applied"
    assert arf_report["steps"]["drift_correct"]["status"] == "applied"
    np.testing.assert_allclose(
        [dataset_by_name[name] for name in names],
        [arf_by_name[name] for name in names],
    )


def test_run_dataset_preprocess_rejects_unknown_normalize():
    from metabolomix.analysis.dataset_analysis import PreconditionError
    ds = _make_ds()
    with pytest.raises(PreconditionError):
        _preprocessed(ds, {"normalize": "not_a_method"})


def test_run_dataset_preprocess_warns_multi_batch_qc():
    """QC が複数バッチ（ファイル名日付）に分かれる場合、arf_preprocess と同じ caveat を出す。"""
    rng = np.random.default_rng(7)
    ds = DatasetState()
    names = [f"ctrl_{i}" for i in range(3)] + [f"treat_{i}" for i in range(3)]
    names += ["20260901_QC_1", "20260902_QC_2"]
    ds.feature_matrix = rng.random((20, len(names))) * 1000.0
    ds.sample_names = names
    ds.feature_ids = [f"f{i}" for i in range(20)]
    report = _preprocessed(ds)
    assert any("複数バッチ" in c for c in report["caveats"])


def test_run_dataset_preprocess_no_multi_batch_caveat_without_dates():
    """既定フィクスチャ（サンプル名に日付なし）では複数バッチ caveat は出ない。"""
    ds = _make_ds()
    report = _preprocessed(ds)
    assert not any("複数バッチ" in c for c in report["caveats"])


# ---------- run_dataset_pca ----------

def test_run_dataset_pca_requires_preprocess():
    from metabolomix.analysis.dataset_analysis import PreconditionError, run_dataset_pca
    ds = _make_ds()
    with pytest.raises(PreconditionError):
        run_dataset_pca(ds, n_components=2)


def test_run_dataset_pca_returns_scores_per_sample():
    from metabolomix.analysis.dataset_analysis import run_dataset_pca
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    result = run_dataset_pca(ds, n_components=2)
    assert len(result["explained_variance_ratio"]) == 2
    assert [s["name"] for s in result["scores"]] == ds.pp_sample_names
    assert "PC1" in result["scores"][0]


def test_run_dataset_pca_reports_run_order_correlation(monkeypatch):
    """QC が無くても、注入順があればドリフトの有無は言える。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_pca
    names = [f"s{i}" for i in range(8)]
    ds = _ds_with_assay_meta(names, {
        f"assay[{i + 1}]": {"run_order": i + 1, "batch": None} for i in range(8)})
    _preprocessed(ds)
    result = run_dataset_pca(ds, n_components=2)
    assert len(result["run_order_correlation"]) == 2
    assert all(v is None or -1.0 <= v <= 1.0 for v in result["run_order_correlation"])


def test_run_dataset_pca_run_order_correlation_is_null_without_injection_order():
    """注入順が無いのに相関を数字で返してはいけない。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_pca
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    result = run_dataset_pca(ds, n_components=2)
    assert result["run_order_correlation"] == [None, None]


def test_run_dataset_pca_caveats_strong_drift_component():
    """主成分が注入順と強く相関したら、生物学と読む前に警告する。"""
    import numpy as np
    from metabolomix.analysis.dataset_analysis import run_dataset_pca
    names = [f"s{i}" for i in range(8)]
    ds = _ds_with_assay_meta(names, {
        f"assay[{i + 1}]": {"run_order": i + 1, "batch": None} for i in range(8)})
    # 注入順に沿って単調に増える強度＝典型的な分析ドリフト。
    ds.feature_matrix = np.array([[100.0 * (i + 1) + f for i in range(8)] for f in range(20)])
    ds.feature_ids = [f"f{i}" for i in range(20)]
    _preprocessed(ds)
    result = run_dataset_pca(ds, n_components=2)
    assert any("注入順" in c for c in result["caveats"])


# ---------- run_dataset_differential ----------

def test_run_dataset_differential_summarizes():
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if n.startswith("ctrl")],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )
    assert result["summary"]["n_tested"] > 0
    assert len(result["results"]) == 20
    assert "q" in result["results"][0]
    assert result["contract_version"] == 1
    assert result["log2fc_sign"] == "positive means group_b is higher"


def test_run_dataset_differential_reports_unknown_samples():
    """存在しないサンプル名を黙って捨てず、caveat で名指しする。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    result = run_dataset_differential(
        ds,
        group_a_samples=["ctrl_0", "ctrl_1", "ctrl_2", "nope_1"],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )
    assert any("nope_1" in c for c in result["caveats"])
    assert result["n_a"] == 3


def test_run_dataset_differential_rejects_small_groups():
    from metabolomix.analysis.dataset_analysis import (
        PreconditionError, run_dataset_differential,
    )
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    with pytest.raises(PreconditionError):
        run_dataset_differential(ds, group_a_samples=["ctrl_0"],
                                 group_b_samples=["treat_0"])


def test_run_dataset_differential_deduplicates_requested_samples():
    """group_a=["S1","S1"] は同じ行を2回数えるだけで、分散0の縮退比較を
    len(idx_a)<2 の拒否をすり抜けて通してしまう。重複は1回に丸め、
    その旨を caveat で名指しする。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    result = run_dataset_differential(
        ds,
        group_a_samples=["ctrl_0", "ctrl_0", "ctrl_1"],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )
    assert result["n_a"] == 2
    assert result["samples_a"] == ["ctrl_0", "ctrl_1"]
    assert any("ctrl_0" in c and "重複" in c for c in result["caveats"])


def test_run_dataset_differential_rejects_group_that_is_all_duplicates():
    """重複除去後に n<2 まで縮む場合、水増しされた見かけの n ではなく
    実際のサイズ不足として拒否されなければならない。"""
    from metabolomix.analysis.dataset_analysis import (
        PreconditionError, run_dataset_differential,
    )
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    with pytest.raises(PreconditionError):
        run_dataset_differential(
            ds,
            group_a_samples=["ctrl_0", "ctrl_0"],
            group_b_samples=["treat_0", "treat_1"],
        )


# ---------- run_dataset_differential: ARF parity caveats (fix round 1) ----------

def _make_batch_ds(names, n_features=20):
    """バッチ判定テスト用。サンプル名の8桁日付から sample_meta.batch が決まる。"""
    rng = np.random.default_rng(11)
    ds = DatasetState()
    ds.feature_matrix = rng.random((n_features, len(names))) * 1000.0
    ds.sample_names = names
    ds.feature_ids = [f"f{i}" for i in range(n_features)]
    return ds


def test_run_dataset_differential_warns_when_not_normalized():
    """arf_differential:846-848 と同じ「正規化未適用」警告。既定 recipe は normalize=none。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)  # ds.preprocessing_recipe は既定の {} のまま
    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if n.startswith("ctrl")],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )
    assert any("正規化が未適用" in c for c in result["caveats"])


def test_run_dataset_differential_no_normalize_caveat_when_normalized():
    """normalize が none 以外なら「正規化未適用」警告は出ない。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    ds.preprocessing_recipe = {"normalize": "median"}
    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if n.startswith("ctrl")],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )
    assert not any("正規化が未適用" in c for c in result["caveats"])


def test_run_dataset_differential_unassessable_confounding_without_dates():
    """既定フィクスチャ（サンプル名に日付なし）は全 batch=None → 交絡評価不可。

    check_confounding のドキュメント（metabolomix/analysis/differential.py:256）どおり、
    「バッチが1つ（この場合は判明ゼロ）＝交絡なし」と誤読させないための assessable=False 経路。
    """
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    ds = _make_ds(n_features=20, n_samples=8)
    _preprocessed(ds)
    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if n.startswith("ctrl")],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )
    assert any("交絡評価不可" in c for c in result["caveats"])
    # バッチ情報自体が無い（filename_date 由来ではない）ので、出所注記は付かない。
    assert not any("バッチはファイル名の日付から推定" in c for c in result["caveats"])


def test_run_dataset_differential_confounded_caveat_when_groups_are_single_batch():
    """group_a が丸ごとバッチA、group_b が丸ごとバッチBなら交絡と判定する。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    names = ([f"20260901_ctrl_{i}" for i in range(3)]
             + [f"20260902_treat_{i}" for i in range(3)])
    ds = _make_batch_ds(names)
    _preprocessed(ds)
    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if "ctrl" in n],
        group_b_samples=[n for n in ds.pp_sample_names if "treat" in n],
    )
    assert any("交絡:" in c for c in result["caveats"])
    assert any("バッチはファイル名の日付から推定" in c for c in result["caveats"])


def test_run_dataset_differential_no_confounding_caveat_when_batches_mixed():
    """両群にバッチA・バッチBが混在していれば、交絡系の caveat はどちらも出ない。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_differential
    names = ["20260901_ctrl_0", "20260902_ctrl_1", "20260901_ctrl_2",
             "20260901_treat_0", "20260902_treat_1", "20260902_treat_2"]
    ds = _make_batch_ds(names)
    _preprocessed(ds)
    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if "ctrl" in n],
        group_b_samples=[n for n in ds.pp_sample_names if "treat" in n],
    )
    assert not any(c.startswith("交絡") for c in result["caveats"])


# ---------- 検出率フィルタ（evidence sidecar 由来） ----------

def test_min_detection_rate_drops_features_below_threshold():
    """gap-fill を除いた実検出率で足切りする（ARF 経路の同名引数と同義）。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess
    ds = _make_ds(n_features=4, n_samples=4)
    # 特徴 0: 4/4 検出、特徴 1: 2/4、特徴 2: 1/4、特徴 3: 0/4
    ds.detected_mask = np.array([
        [True, True, True, True],
        [True, True, False, False],
        [True, False, False, False],
        [False, False, False, False],
    ])
    ds.feature_qc = {"source": "arf"}

    _, _, feature_names, _, _, report = run_dataset_preprocess(
        ds, {"impute": "none", "min_detection_rate": 0.5})

    assert feature_names == ["f0", "f1"]
    assert report["detection_filter"]["min_detection_rate"] == 0.5
    assert report["detection_filter"]["features_removed"] == 2


def test_min_detection_rate_zero_keeps_everything():
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess
    ds = _make_ds(n_features=3, n_samples=4)
    ds.detected_mask = np.zeros((3, 4), dtype=bool)
    ds.feature_qc = {"source": "arf"}
    _, _, feature_names, _, _, report = run_dataset_preprocess(
        ds, {"impute": "none", "min_detection_rate": 0.0})
    assert len(feature_names) == 3
    assert "detection_filter" not in report


def test_min_detection_rate_without_detection_state_is_a_bad_request():
    """検出状態が無いのに足切りを求められたら、黙って 0 扱いにせず拒否する。"""
    from metabolomix.analysis.dataset_analysis import (
        PreconditionError, run_dataset_preprocess)
    ds = _make_ds(n_features=3, n_samples=4)
    assert ds.detected_mask is None
    with pytest.raises(PreconditionError) as excinfo:
        run_dataset_preprocess(ds, {"min_detection_rate": 0.5})
    assert excinfo.value.kind == "bad_request"
    assert "検出" in str(excinfo.value)


def test_detection_rate_is_reported_even_without_filtering():
    """フィルタを掛けなくても、検出状態があれば実検出率を報告する。"""
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess
    ds = _make_ds(n_features=2, n_samples=4)
    ds.detected_mask = np.array([[True, True, True, True],
                                 [True, False, False, False]])
    ds.feature_qc = {"source": "arf"}
    _, _, _, _, _, report = run_dataset_preprocess(ds, {"impute": "none"})
    assert report["detection"]["gap_filled_rate"] == 0.375
