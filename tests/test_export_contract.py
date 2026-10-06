"""差次的エクスポート契約（metabolomix.analysis.export_contract）の単体テスト。

ARF 経路（arf_export_differential）と DatasetState 経路
（dataset_export_differential）の両方がこのモジュールを参照する。契約の列・版・
メタ行順序はここで固定し、二重定義を作らない（arf/tools.py 側は同一性を確認する）。
"""


def test_contract_constants_are_shared_with_arf():
    """ARF 側が契約モジュールの定数を参照していること（二重定義を作らない）。"""
    from metabolomix.analysis import export_contract
    from metabolomix.arf import tools as arf_tools
    assert arf_tools._EXPORT_COLUMNS is export_contract.EXPORT_COLUMNS
    assert arf_tools._DIFFERENTIAL_CONTRACT_VERSION == export_contract.CONTRACT_VERSION
    assert arf_tools._LOG2FC_SIGN == export_contract.LOG2FC_SIGN


def test_export_columns_are_frozen():
    """列と順序は下流との契約。変更は contract_version の上げ方とセット。"""
    from metabolomix.analysis.export_contract import EXPORT_COLUMNS
    assert EXPORT_COLUMNS == [
        "spot_id", "name", "name_source", "ontology", "inchikey",
        "inchikey_source", "msi_level", "mz", "rt", "log2fc",
        "p_value", "q_value", "mean_a", "mean_b", "significant",
    ]


def test_format_number_blanks_non_finite():
    from metabolomix.analysis.export_contract import format_number
    assert format_number(None, ".4f") == ""
    assert format_number(float("nan"), ".4f") == ""
    assert format_number(float("inf"), ".4f") == ""
    assert format_number(1.23456, ".4f") == "1.2346"


def _meta(**over):
    from metabolomix.analysis.export_contract import build_meta
    kwargs = dict(group_a="ctrl", n_a=4, group_b="treat", n_b=4,
                  q_threshold=0.05, log2fc_threshold=1.0, log_transform=True,
                  n_features_total=100, n_with_inchikey=80, n_unannotated=20,
                  msi_note="# msi_level は注釈確度。MS/MS の有無ではない")
    kwargs.update(over)
    return build_meta(**kwargs)


def test_build_meta_has_mandatory_lines():
    lines = _meta()
    joined = "\n".join(lines)
    assert "# contract_version = 1" in joined
    assert "# log2fc_sign = positive means group_b is higher" in joined
    assert "n_a = 4" in joined and "n_b = 4" in joined
    assert "n_unannotated = 20" in joined
    assert all(line.startswith("#") for line in lines)


def test_build_meta_places_optional_slots_in_contract_order():
    """行の順序も契約。source は exported_at の直後、preprocess は閾値の直後。"""
    lines = _meta(source_lines=["# source_arf = a.arf", "# source_arf2 = a.arf2"],
                  preprocess_line="# preprocess = {'normalize': 'none'}")
    keys = [l.split(" = ")[0].split("\t")[0] for l in lines]
    assert keys == [
        "# contract_version", "# exported_at",
        "# source_arf", "# source_arf2",
        "# group_a", "# group_b", "# log2fc_sign", "# q_threshold",
        "# preprocess", "# n_features_total",
        "# msi_level は注釈確度。MS/MS の有無ではない",
    ]


def test_build_meta_omits_preprocess_when_absent():
    assert not any("preprocess" in l for l in _meta())


def test_run_dataset_differential_tracks_contract_module(monkeypatch):
    """dataset_analysis.py がバージョン/符号ラベルをローカル複製せず、
    export_contract を都度参照していることを行動で確認する。

    値の一致（==）だけを見るテストは、dataset_analysis.py 側に値の等しい
    ローカル定数（旧 _CONTRACT_VERSION/_LOG2FC_SIGN）が残っていても偶然パスして
    しまう。ここでは export_contract 側だけを書き換え、その変更が
    run_dataset_differential の戻り値へ伝播することを見て、参照がコピーでなく
    リンクであることを確認する（コントローラ裁定: これを直さないと、
    CONTRACT_VERSION を上げた瞬間 dataset_export_differential が永久に
    「契約非互換」で拒否し続ける無限ループになる）。
    """
    import numpy as np

    from metabolomix.analysis import export_contract
    from metabolomix.analysis.dataset_analysis import (
        run_dataset_differential, run_dataset_preprocess,
    )
    from metabolomix.mztab.dataset_state import DatasetState

    monkeypatch.setattr(export_contract, "CONTRACT_VERSION", 999)
    monkeypatch.setattr(export_contract, "LOG2FC_SIGN", "patched sign for linkage test")

    rng = np.random.default_rng(0)
    ds = DatasetState()
    ds.feature_matrix = rng.random((20, 8)) * 1000.0
    ds.sample_names = ([f"ctrl_{i}" for i in range(4)]
                       + [f"treat_{i}" for i in range(4)])
    ds.feature_ids = [f"f{i}" for i in range(20)]
    (ds.pp_matrix, ds.pp_sample_names, ds.pp_feature_names,
     ds.roles, ds.sample_meta, _) = run_dataset_preprocess(ds, {})

    result = run_dataset_differential(
        ds,
        group_a_samples=[n for n in ds.pp_sample_names if n.startswith("ctrl")],
        group_b_samples=[n for n in ds.pp_sample_names if n.startswith("treat")],
    )

    assert result["contract_version"] == 999
    assert result["log2fc_sign"] == "patched sign for linkage test"


# --- significant 判定（契約 15 列目）---

def test_is_significant_at_thresholds_is_inclusive():
    """境界ちょうどは有意。q <= 閾値、|log2fc| >= 閾値。"""
    from metabolomix.analysis.export_contract import is_significant
    assert is_significant(q=0.05, log2fc=1.0, q_threshold=0.05, log2fc_threshold=1.0)
    assert is_significant(q=0.05, log2fc=-1.0, q_threshold=0.05, log2fc_threshold=1.0)
    assert not is_significant(q=0.051, log2fc=1.0, q_threshold=0.05, log2fc_threshold=1.0)
    assert not is_significant(q=0.05, log2fc=0.999, q_threshold=0.05, log2fc_threshold=1.0)


def test_is_significant_rejects_missing_and_non_finite():
    """欠測・NaN・inf は有意にしない（下流が真偽を鵜呑みにするため）。"""
    from metabolomix.analysis.export_contract import is_significant
    nan, inf = float("nan"), float("inf")
    assert not is_significant(q=None, log2fc=1.0, q_threshold=0.05, log2fc_threshold=1.0)
    assert not is_significant(q=0.01, log2fc=None, q_threshold=0.05, log2fc_threshold=1.0)
    assert not is_significant(q=nan, log2fc=1.0, q_threshold=0.05, log2fc_threshold=1.0)
    assert not is_significant(q=0.01, log2fc=nan, q_threshold=0.05, log2fc_threshold=1.0)
    assert not is_significant(q=0.01, log2fc=inf, q_threshold=0.05, log2fc_threshold=1.0)


def test_both_export_paths_use_the_contract_significance():
    """ARF 経路と DatasetState 経路が別実装を持たないこと。

    キー名が `q_value` / `q` で違うだけの同一ロジックが 2 箇所にあると、
    片方だけ直した状態で下流は同一契約として読む。
    """
    from metabolomix.arf import tools as arf_tools
    from metabolomix.tools import dataset_analysis_tools
    assert not hasattr(arf_tools, "_is_significant")
    assert not hasattr(dataset_analysis_tools, "_is_significant")
