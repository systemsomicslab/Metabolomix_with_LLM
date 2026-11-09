"""メタデータ訂正の契約（spec §7.4）。

`apply_sample_manifest` は parse_manifest / resolve_metadata / apply_metadata をこの順に繋ぐだけで、
失敗時にデータセットを一切変えない契約を壊さない。群だけの訂正は PCA を生かし、差次的解析だけを捨てる。
比較の前提（群の向き・試料の選択・交絡）は統計エンジンの `dataset_compare` が持つ（`tests/test_stats_*`）。
"""
from __future__ import annotations

import pytest

from tests.pipeline_fixtures import make_dataset

from metabolomix.core.atomic_io import DomainError


# ---------- dataset_set_sample_metadata（MCP、Task 9のparse/resolve/applyを接続） ----------

def _job_dataset_with_raws(tmp_path):
    """assay_sources（実行時rawの絶対パス）を持つDatasetStateを合成する。

    parse_manifestのsource_root/expected_sources逆算(sample_manifest._raw_manifest_layout)は
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
    from metabolomix.analysis.sample_manifest import apply_sample_manifest
    from tests.engine_chain import register_two_group

    ds = _job_dataset_with_raws(tmp_path)
    groups_v1 = ["control"] * 4 + ["treated"] * 4
    manifest_v1 = tmp_path / "sample-manifest-v1.tsv"
    manifest_v1.write_text(_manifest_text(ds, groups_v1), encoding="utf-8")
    apply_sample_manifest(ds, str(manifest_v1))

    register_two_group(ds, "control", "treated")
    pca = ds.last_pca = {"kind": "stat_pca", "provenance": {"result_id": "res_pca"}}
    matrices = dict(ds.stat_matrices)
    assert ds.last_differential is not None and matrices

    # group以外は全列同一のまま、idx3だけをtreatedへ付け替えた2回目の適用。
    groups_v2 = list(groups_v1)
    groups_v2[3] = "treated"
    manifest_v2 = tmp_path / "sample-manifest-v2.tsv"
    manifest_v2.write_text(_manifest_text(ds, groups_v2), encoding="utf-8")
    apply_sample_manifest(ds, str(manifest_v2))

    assert ds.last_pca is pca
    assert ds.stat_matrices == matrices
    assert ds.last_differential is None
    assert ds.sample_metadata_rows[3]["group"] == "treated"


def test_dataset_set_sample_metadata_does_not_mutate_session_on_failure(tmp_path):
    """シートが不正なら、apply_metadataの契約どおりdsを一切変更しない。"""
    from metabolomix.analysis.sample_manifest import apply_sample_manifest

    ds = _job_dataset_with_raws(tmp_path)
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
