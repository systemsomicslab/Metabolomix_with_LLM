"""比較前提（reference/test群の解決・交絡判定）とメタデータ訂正MCPの契約（spec §7.4）。

resolve_comparison / run_comparison は数値・ドメイン層の純関数で、群名から向きを
推測しない・QC/blank/unknown/include=falseを比較対象から除外する・完全交絡を
allow_confounded明示なしには通さない、という前提検証だけを持つ。
dataset_set_sample_metadata はTask 9のparse_manifest/resolve_metadata/apply_metadata
をMCPへ薄く繋ぐだけで、失敗時にsessionを一切変えない契約を壊さない。
"""
from __future__ import annotations

import pytest

from tests.pipeline_fixtures import make_dataset, metadata_rows

from metabolomix.core.atomic_io import DomainError


def test_group_names_do_not_determine_reference_implicitly():
    """comparisonにreference_group/test_groupが無ければ、群名がどれだけそれらしくても
    比較は実行しない（spec §7.4「群が2種類あるだけでは対照と処置の向きを決めない」）。"""
    from metabolomix.analysis.dataset_service import resolve_comparison

    ds = make_dataset()
    rows = metadata_rows(ds, n_qc=0)
    for i, row in enumerate(rows):
        row["group"] = "control" if i < 4 else "treated"
    with pytest.raises(DomainError, match="COMPARISON_REQUIRED"):
        resolve_comparison(ds, {}, rows)


def _labelled_rows(ds, *, n_control=4):
    """全8検体にcontrol/treatedの群ラベルを振っただけの行を返す（role=sample・include=true）。"""
    rows = metadata_rows(ds, n_qc=0)
    for i, row in enumerate(rows):
        row["group"] = "control" if i < n_control else "treated"
    return rows


def _confirm(row: dict, field: str, value) -> None:
    """メタデータ1行の1フィールドを`user_manifest`確定値として書き換える。"""
    row[field] = value
    row["provenance"][field] = {
        "value": value, "source": "user_manifest", "confidence": "confirmed"}


def test_resolve_comparison_excludes_qc_blank_unknown_and_include_false():
    """QC/blank/role=unknown/include=falseは、比較対象の群ラベルに合致していても
    比較には混ざらない（brief step3・spec §7.4）。"""
    from metabolomix.analysis.dataset_service import resolve_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    _confirm(rows[0], "role", "qc")        # control: 除外
    _confirm(rows[1], "role", "unknown")   # control: 除外
    _confirm(rows[4], "include", False)    # treated: 除外

    resolved = resolve_comparison(
        ds, {"reference_group": "control", "test_group": "treated"}, rows)

    assert set(resolved["reference_samples"]) == {ds.sample_names[2], ds.sample_names[3]}
    assert set(resolved["test_samples"]) == {
        ds.sample_names[5], ds.sample_names[6], ds.sample_names[7]}


def test_resolve_comparison_rejects_duplicate_sample_id():
    from metabolomix.analysis.dataset_service import resolve_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    rows[1]["sample_id"] = rows[0]["sample_id"]

    with pytest.raises(DomainError, match="SAMPLE_MANIFEST_INVALID"):
        resolve_comparison(
            ds, {"reference_group": "control", "test_group": "treated"}, rows)


def test_resolve_comparison_rejects_group_with_fewer_than_two_samples():
    from metabolomix.analysis.dataset_service import resolve_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    # controlに該当するのはidx0の1件だけにする(idx1-3をtreatedへ回す)。
    for row in rows[1:4]:
        row["group"] = "treated"

    with pytest.raises(DomainError, match="COMPARISON_GROUP_TOO_SMALL"):
        resolve_comparison(
            ds, {"reference_group": "control", "test_group": "treated"}, rows)


def test_run_comparison_direction_positive_log2fc_means_test_group_higher():
    """reference_group→group_a、test_group→group_bへ渡す（正のlog2FCはtest_groupが高い）。"""
    from metabolomix.analysis.dataset_service import preprocess_dataset, run_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})

    result = run_comparison(
        ds, {"reference_group": "control", "test_group": "treated"}, rows)

    assert result["a"] == "control"
    assert result["b"] == "treated"
    assert result["log2fc_sign"] == "positive means group_b is higher"
    assert result["provenance"]["comparison"]["reference_group"] == "control"
    assert result["provenance"]["comparison"]["test_group"] == "treated"


def test_run_comparison_stops_on_complete_confounding_without_override():
    """各群が単一バッチに完全に一致する（群⟂バッチが分離不能）場合は既定で止める。"""
    from metabolomix.analysis.dataset_service import preprocess_dataset, run_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    for row in rows[:4]:
        _confirm(row, "batch", "B1")
    for row in rows[4:]:
        _confirm(row, "batch", "B2")
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})

    with pytest.raises(DomainError, match="CONFOUNDED_COMPARISON"):
        run_comparison(ds, {"reference_group": "control", "test_group": "treated"}, rows)

    # 探索結果（last_differential）は残さない——確認しないまま「直近の結果」を語らせない。
    assert ds.last_differential is None


def test_run_comparison_continues_with_explicit_allow_confounded():
    """allow_confounded=trueを明示すれば継続し、未調整である旨を来歴とcaveatsへ残す。

    spec §7.4「未調整であることを図・TSV付随メタ・レポートへ記録して継続できる」の
    記録先は永続化されたprovenance（`result["provenance"]["warnings"]`）でなければ
    ならない——`compare_dataset`は`result["caveats"]`を先にスナップショットして
    provenance.warningsを作るため、run_comparisonがそのあとに追記したこの注記が
    provenance側に反映されるかは別に確認する必要がある（レビュー Finding 2）。
    """
    from metabolomix.analysis.dataset_service import preprocess_dataset, run_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    for row in rows[:4]:
        _confirm(row, "batch", "B1")
    for row in rows[4:]:
        _confirm(row, "batch", "B2")
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})

    result = run_comparison(
        ds,
        {"reference_group": "control", "test_group": "treated", "allow_confounded": True},
        rows,
    )

    assert result["provenance"]["comparison"]["unadjusted_confounded"] is True
    assert result["provenance"]["comparison"]["confounding"]["confounded"] is True
    assert any("未調整" in c for c in result["caveats"])
    # 永続化されたprovenance側にも同じ注記が残っていること（caveatsだけでは
    # 図/TSV/レポートの来歴側から読めない——Task 17がこのprovenanceを読む）。
    assert any("未調整" in w for w in result["provenance"]["warnings"])


def test_run_comparison_insufficient_batch_info_is_not_confounded():
    """batch情報が無い（欠落）場合は「評価不可」であって「交絡あり」ではない
    ——完全交絡と決め付けて止めない（spec §7.4・common-context 参照）。"""
    from metabolomix.analysis.dataset_service import preprocess_dataset, run_comparison

    ds = make_dataset()
    rows = _labelled_rows(ds)
    for row in rows:
        row["batch"] = None
        row["provenance"]["batch"] = {"value": None, "source": None, "confidence": None}
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})

    result = run_comparison(ds, {"reference_group": "control", "test_group": "treated"}, rows)

    assert result["provenance"]["comparison"]["confounding"]["confounded"] is False
    assert result["provenance"]["comparison"]["confounding"]["assessable"] is False


# ---------- dataset_set_sample_metadata（MCP、Task 9のparse/resolve/applyを接続） ----------

def _job_dataset_with_raws(tmp_path):
    """assay_sources（実行時rawの絶対パス）を持つDatasetStateを合成する。

    parse_manifestのsource_root/expected_sources逆算(_raw_manifest_layout)は
    ds.assay_sourcesから行うため、実ファイルを伴うrawパスが要る。
    """
    ds = make_dataset()
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    ds.assay_sources = {}
    for i, name in enumerate(ds.sample_names):
        raw_path = raw_dir / f"{name}.wiff"
        raw_path.write_text("x", encoding="utf-8")
        ds.assay_sources[f"abundance_assay[{i + 1}]"] = str(raw_path)
    return ds


def _manifest_text(ds, groups: list[str]) -> str:
    header = "# schema = sample-manifest.v1\n"
    header += "sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude\n"
    lines = [header.rstrip("\n")]
    for i, name in enumerate(ds.sample_names):
        lines.append(f"{name}\t{name}.wiff\tsample\t{groups[i]}\tB1\t{i + 1}\t\ttrue")
    return "\n".join(lines) + "\n"


def test_dataset_set_sample_metadata_reuses_pca_and_invalidates_differential(tmp_path):
    """groupだけを訂正した2回目の適用は、PCAを生かしたまま差次的解析だけを無効化する
    （group以外の列は前回と同一——Task 6 invalidate_resultsの「group変更だけならPCAは
    無効化しない」規則、common-context 参照）。1回目の適用は全列確定として扱われる
    （Task 9 apply_metadataの契約）ので、その後にPCA・差次的解析を作ってから2回目を送る。
    """
    from metabolomix.analysis.dataset_service import (
        apply_sample_manifest, compare_dataset, pca_dataset, preprocess_dataset,
    )

    ds = _job_dataset_with_raws(tmp_path)
    groups_v1 = ["control"] * 4 + ["treated"] * 4
    manifest_v1 = tmp_path / "sample-manifest-v1.tsv"
    manifest_v1.write_text(_manifest_text(ds, groups_v1), encoding="utf-8")
    apply_sample_manifest(ds, str(manifest_v1))

    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})
    pca_dataset(ds, n_components=2)
    old_pca_id = ds.last_pca["provenance"]["result_id"]
    compare_dataset(ds, ds.sample_names[:4], ds.sample_names[4:])
    assert ds.last_differential is not None

    # group以外は全列同一のまま、idx3だけをtreatedへ付け替えた2回目の適用。
    groups_v2 = list(groups_v1)
    groups_v2[3] = "treated"
    manifest_v2 = tmp_path / "sample-manifest-v2.tsv"
    manifest_v2.write_text(_manifest_text(ds, groups_v2), encoding="utf-8")
    apply_sample_manifest(ds, str(manifest_v2))

    assert ds.last_pca is not None
    assert ds.last_pca["provenance"]["result_id"] == old_pca_id
    assert ds.last_differential is None
    assert ds.sample_metadata_rows[3]["group"] == "treated"


def test_dataset_set_sample_metadata_does_not_mutate_session_on_failure(tmp_path):
    """シートが不正なら、apply_metadataの契約どおりdsを一切変更しない。"""
    from metabolomix.analysis.dataset_service import apply_sample_manifest, preprocess_dataset

    ds = _job_dataset_with_raws(tmp_path)
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})
    before_revision = ds.metadata_revision
    before_rows = ds.sample_metadata_rows

    manifest = tmp_path / "bad-manifest.tsv"
    # role列が不正な値("invalid_role")を持つ壊れたシート。
    groups = ["control"] * 4 + ["treated"] * 4
    text = _manifest_text(ds, groups).replace("sample\tcontrol", "invalid_role\tcontrol", 1)
    manifest.write_text(text, encoding="utf-8")

    with pytest.raises(DomainError, match="SAMPLE_MANIFEST_INVALID"):
        apply_sample_manifest(ds, str(manifest))

    assert ds.metadata_revision == before_revision
    assert ds.sample_metadata_rows is before_rows
