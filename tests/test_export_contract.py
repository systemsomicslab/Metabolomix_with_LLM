"""差次的エクスポート契約（metabolomix.analysis.export_contract）の単体テスト。

書き出し（dataset_export_differential と pipeline の書き出し）がこのモジュールを参照する。
契約の列・版・メタ行順序はここで固定し、二重定義を作らない。
"""


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
    """書き出し側が有意判定を別実装として持たないこと（契約モジュールだけが持つ）。

    同一ロジックが 2 箇所にあると、片方だけ直した状態で下流は同一契約として読む
    （旧 ARF 経路は B-2b で消えた。ここでは残った書き出し側を縛る）。
    """
    from metabolomix.tools import dataset_analysis_tools
    assert not hasattr(dataset_analysis_tools, "_is_significant")


def test_build_meta_appends_stats_lines_after_msi_note():
    """統計エンジンの行は msi_note の後ろだけ。足さなければ今までの出力と同じ。"""
    base = _meta()
    lines = _meta(stats_lines=["# test = moderated_t\tci_level = 0.95", "# unit = subject"])
    assert lines[-3].startswith("# msi_level") and lines[-2:] == [
        "# test = moderated_t\tci_level = 0.95", "# unit = subject"]
    assert [l for l in lines if "exported_at" not in l][:-2] == [l for l in base if "exported_at" not in l]


def test_stats_meta_lines_empty_for_v1_result_and_keys_do_not_collide():
    from metabolomix.analysis.export_contract import stats_meta_lines
    assert stats_meta_lines({"kind": "two_group", "a": "x"}) == []
    result = {"effect_size_definition": {"name": "log2_mean_difference"}, "matrix_id": "m", "rule_hash": "h", "contrast_id": "c", "provenance": {"kind": "compare"},
              "statistics": {"test": "welch", "ci_level": 0.95, "subject_source": "injection",
                             "verdict_rule": {"q_threshold": 0.05, "equivalence_fc": 1.2},
                             "bh_population": {"name": "matrix", "n": 7, "q_interpretation": "discovery_fdr"},
                             "n_units_a": 3, "n_units_b": 4}}
    lines = stats_meta_lines(result)
    keys = _meta_keys(lines)
    assert len(keys) == len(set(keys))
    # 既存のメタ行は**タブで割ったすべての** `key = value` を見る（1 行に複数のキーがある。下流の
    # parse_differential はそれを全部拾い、同じキーは後勝ちで消える）
    existing = _meta_keys(_meta(source_lines=_DATASET_SOURCE_LINES,
                                preprocess_line="# preprocess = {}"))
    assert {"n_a", "n_b", "log2fc_threshold", "log_transform", "n_with_inchikey", "n_unannotated",
            "curation_flags", "unadjusted_confounded_note"} <= set(existing)   # 先頭以外のキーも拾えている
    assert not set(keys) & set(existing)
    assert "# unit = injection\tn_units_a = 3\tn_units_b = 4" in lines
    assert "significant_rule" in keys


#: dataset_export が source_lines に差す行（キーの形だけ。値は問わない）
_DATASET_SOURCE_LINES = [
    "# source_mztab = a.mzTab", "# source_job = j", "# id_space = mztab_smf_id", "# result_id = r",
    "# preprocess_id = p", "# source_verification = v", "# unadjusted_confounded = true",
    "# unadjusted_confounded_note = x",
    "# curation = applied\tcuration_flags = 1\tcuration_wrong_excluded = 0\tcuration_suspect = 0"
    "\tcuration_assigned = 0\tcuration_redundant_excluded = 0\tcuration_flags_sha256 = d",
]


def _meta_keys(lines) -> list[str]:
    """massbank-context の parse_differential と同じ読み方で、各行のタブ区切りの全キーを並べる。"""
    return [chunk.partition("=")[0].strip() for line in lines
            for chunk in line.lstrip("#").split("\t") if "=" in chunk]


def test_stats_meta_lines_write_zero_for_an_empty_bh_population():
    """BH の母集団が空（対比の n が無い）なら `bh_n = 0`。"None" という文字列は書かない。"""
    from metabolomix.analysis.export_contract import stats_meta_lines
    result = {"effect_size_definition": {"name": "log2_mean_difference"}, "matrix_id": "m",
              "rule_hash": "h", "contrast_id": "c", "provenance": {"kind": "compare"},
              "statistics": {"test": "moderated_t", "ci_level": 0.95, "subject_source": "subject_id",
                             "verdict_rule": {"q_threshold": 0.05, "equivalence_fc": 1.2},
                             "bh_population": {"name": "annotated", "n": None,
                                               "q_interpretation": "discovery_fdr"},
                             "n_units_a": 3, "n_units_b": 4}}
    text = "\n".join(stats_meta_lines(result))
    assert "\tbh_n = 0\t" in text and "None" not in text


def test_row_significant_honours_verdict_and_row_flag():
    from metabolomix.analysis.export_contract import row_significant
    thr = dict(q_threshold=0.05, log2fc_threshold=1.0)
    assert row_significant({"q": 0.01, "log2fc": 2.0}, **thr)                       # v1 の行
    assert not row_significant({"q": 0.01, "log2fc": 0.5}, **thr)
    assert row_significant({"q": 0.01, "log2fc": 2.0, "verdict": "changed"}, **thr)
    for verdict in ("inconclusive", "unchanged"):
        assert not row_significant({"q": 0.01, "log2fc": 2.0, "verdict": verdict}, **thr)
    assert not row_significant({"q": 0.01, "log2fc": 2.0, "significant": False}, **thr)
    assert row_significant({"q": 0.9, "log2fc": 0.0, "significant": True}, **thr)   # 行の印が勝つ
