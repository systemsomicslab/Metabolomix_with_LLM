"""注入単位のRT/m/z証拠（spec §9.1 `assay-feature-evidence.v1`）。

QCと同定の根拠は**注入ごと**の実測値でなければならない。アライメント代表値
（SMFのRT/m/z、`.arf`スポットの代表行）を各注入へ複製すると、注入間のばらつきが
消えたまま「全注入でRTが一致した」という嘘になる。ここで縛るのは3点:

1. 値は `AlignedPeakProperties` の**行ごと**に読む。2注入の値は別物として残る。
2. 単位変換は明示された元単位からだけ行う。不明単位は推測せず
   `EVIDENCE_UNIT_UNSUPPORTED` で止める。
3. feature軸（`.arf`スポット ↔ mzTab SMF_ID）とassay軸（注入 ↔ job の
   source→assay対応）のどちらかでも照合できなければ、部分的な表を作らず
   `availability=false` と理由を返す。名前一致や行順の一致で代用しない。

fixtureは全て合成。実データ・実Consoleの検証はTask 15。
"""
from __future__ import annotations

import lz4.block
import msgpack
import pytest

from metabolomix.analysis.assay_evidence import (
    EVIDENCE_SCHEMA,
    build_assay_evidence,
    normalize_cell,
)
from metabolomix.console.profile_adapter import adapter_capabilities
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState


# ---------- 合成 `.arf`（tests/test_arf_multiblock.py と同じ組み立て） ----------

def _peak_row(file_id: int, sample_name: str, *, mz: float, rt: float,
              height: float = 1000.0, gap_filled: bool = False) -> list:
    row = [None] * 38
    row[0] = file_id
    row[1] = sample_name
    row[2] = -2 if gap_filled else file_id + 100
    row[3] = file_id + 1000
    row[10] = {}
    row[15] = [[1, [rt]], [3, [mz]]]
    row[18] = height
    row[20] = height * 2
    row[21] = height * 1.5
    row[22] = mz
    row[23] = 0
    row[24] = ""
    row[37] = [1.0, 20.0]
    return row


def _arf_bytes(groups: list[list]) -> bytes:
    inner = b"".join(msgpack.packb(item, use_bin_type=True)
                     for item in [[0, 0, *groups]])
    payload = (msgpack.packb(len(inner), use_bin_type=True)
               + lz4.block.compress(inner, store_size=False))
    return msgpack.packb(msgpack.ExtType(99, payload), use_bin_type=True)


def _write_arf(tmp_path, groups: list[list]):
    path = tmp_path / "AlignmentResult_PeakProperties.arf"
    path.write_bytes(_arf_bytes(groups))
    return path


def _dataset(tmp_path, *, n_features=2, raw_names=("S1", "S2", "S3"), mzs=None):
    """`.arf` と同じ特徴数・注入数を持つ最小のDatasetState。

    `assay_sources` は job の source→assay 対応（abundance列 → 予定raw の絶対パス）。
    これが assay 軸の唯一の照合材料で、assay 表示名では代用しない。

    注入は3件。`.arf` の既存リーダは「グループ=スポット」を行数3以上で判定する
    （`arf/reader.py` `extract_arf_data`）ので、2注入のスポットは読めない。
    """
    ds = DatasetState()
    ds.feature_ids = [str(i) for i in range(n_features)]
    ds.sample_names = list(raw_names)
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(len(raw_names))]
    mz_list = mzs if mzs is not None else [700.5 + i for i in range(n_features)]
    ds.feature_metadata = {str(i): {"mz": mz_list[i], "rt": 1.0}
                           for i in range(n_features)}
    ds.assay_sources = {f"abundance_assay[{i + 1}]": str(tmp_path / f"{name}.wiff")
                        for i, name in enumerate(raw_names)}
    return ds


def _adapter():
    return adapter_capabilities("msdial5")


# ---------- normalize_cell ----------

def test_per_injection_rt_is_not_representative_rt():
    a = normalize_cell({"rt": 60., "m_z": 104.}, "second")
    b = normalize_cell({"rt": 66., "m_z": 104.001}, "second")
    assert (a["observed_rt_min"], b["observed_rt_min"]) == (1., 1.1)
    assert a["observed_mz"] != b["observed_mz"]


def test_minute_source_is_not_rescaled():
    cell = normalize_cell({"rt": 3.5, "m_z": 700.5}, "minute")
    assert cell["observed_rt_min"] == 3.5


def test_missing_rt_stays_missing_with_a_reason():
    cell = normalize_cell({"rt": None, "m_z": 700.5}, "minute")
    assert cell["observed_rt_min"] is None
    assert cell["missing_reason"] == "rt_absent_in_artifact"


def test_unknown_unit_raises_value_error_in_the_pure_layer():
    with pytest.raises(ValueError):
        normalize_cell({"rt": 1.0, "m_z": 1.0}, "hour")


# ---------- build_assay_evidence: 正常系 ----------

def test_evidence_rows_are_per_feature_and_per_assay(tmp_path):
    ds = _dataset(tmp_path)
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0),
         _peak_row(1, "S2", mz=700.52, rt=3.1),
         _peak_row(2, "S3", mz=700.49, rt=3.05)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0),
         _peak_row(1, "S2", mz=701.5, rt=4.0, gap_filled=True),
         _peak_row(2, "S3", mz=701.51, rt=4.02)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())

    assert result["schema"] == EVIDENCE_SCHEMA
    assert result["availability"] is True
    assert len(result["rows"]) == 6          # 特徴2 × 注入3

    by_key = {(r["feature_id"], r["assay_id"]): r for r in result["rows"]}
    assert by_key[("0", "assay[1]")]["observed_rt_min"] == 3.0
    assert by_key[("0", "assay[2]")]["observed_rt_min"] == 3.1
    # 代表値の複製なら両注入が同値になる。別物として残っていることを見る。
    assert (by_key[("0", "assay[1]")]["observed_mz"]
            != by_key[("0", "assay[2]")]["observed_mz"])
    assert by_key[("1", "assay[2]")]["detection_status"] == "gap_filled"
    assert by_key[("1", "assay[1]")]["detection_status"] == "detected"


def test_rows_carry_artifact_hash_and_reader_version(tmp_path):
    ds = _dataset(tmp_path)
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0), _peak_row(1, "S2", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0), _peak_row(1, "S2", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    row = result["rows"][0]
    assert row["source_artifact_hash"] == result["source_artifact_hash"]
    assert len(row["source_artifact_hash"]) == 64
    assert row["reader_version"] == result["reader_version"]
    assert row["dataset_id"] == ds.dataset_id
    assert result["rt_unit"] == {"source": "minute", "normalized": "minute"}


def test_assay_axis_follows_the_job_mapping_not_the_row_order(tmp_path):
    """`.arf`の行順（FileID順）とmzTabのassay順が逆でも、値は注入に従う。"""
    ds = _dataset(tmp_path, raw_names=("S1", "S2", "S3"))
    artifact = _write_arf(tmp_path, [
        [_peak_row(7, "S2", mz=700.5, rt=9.0),     # 先に S2 が来る
         _peak_row(9, "S3", mz=700.5, rt=5.0),
         _peak_row(3, "S1", mz=700.5, rt=3.0)],
        [_peak_row(7, "S2", mz=701.5, rt=9.5),
         _peak_row(9, "S3", mz=701.5, rt=5.5),
         _peak_row(3, "S1", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    by_key = {(r["feature_id"], r["assay_id"]): r for r in result["rows"]}
    assert by_key[("0", "assay[1]")]["observed_rt_min"] == 3.0
    assert by_key[("0", "assay[2]")]["observed_rt_min"] == 9.0


# ---------- build_assay_evidence: 照合できない場合 ----------

def _reasons(result) -> list[str]:
    return [r["code"] for r in result["reasons"]]


def test_unresolved_assay_source_yields_no_rows(tmp_path):
    ds = _dataset(tmp_path, raw_names=("S1", "S2", "S3"))
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0),
         _peak_row(1, "OTHER", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0),
         _peak_row(1, "OTHER", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert "assay_source_unresolved" in _reasons(result)
    assert result["rows"] == []


def test_missing_job_source_mapping_is_not_replaced_by_names(tmp_path):
    """assay_sourcesが無いとき、表示名一致へ落ちない（v2の厳密照合の代用禁止）。"""
    ds = _dataset(tmp_path)
    ds.assay_sources = {}
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0), _peak_row(1, "S2", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0), _peak_row(1, "S2", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert "assay_source_unavailable" in _reasons(result)


def test_duplicate_assay_key_in_one_spot_is_rejected(tmp_path):
    ds = _dataset(tmp_path)
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0),
         _peak_row(1, "S1", mz=700.5, rt=3.4),      # 同じrawが2行
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0),
         _peak_row(1, "S1", mz=701.5, rt=4.4),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert "duplicate_assay_key" in _reasons(result)


def test_feature_count_mismatch_is_rejected(tmp_path):
    ds = _dataset(tmp_path, n_features=3)
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0), _peak_row(1, "S2", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0), _peak_row(1, "S2", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert "feature_count_mismatch" in _reasons(result)


def test_feature_axis_mz_mismatch_is_rejected(tmp_path):
    """1つずらした誤接合はm/zがDaオーダーで食い違う。名前では気付けない。"""
    ds = _dataset(tmp_path, n_features=2, mzs=[701.5, 700.5])
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0), _peak_row(1, "S2", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0), _peak_row(1, "S2", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert "feature_axis_mz_mismatch" in _reasons(result)


def test_representative_only_artifact_has_no_per_injection_evidence(tmp_path):
    """代表行だけで注入行を持たない成果物から、代表値を複製して埋めない。"""
    ds = _dataset(tmp_path)
    artifact = tmp_path / "representative.arf"
    inner = b"".join(msgpack.packb(item, use_bin_type=True)
                     for item in [[0, 0, ["header-only"], ["also-not-a-group"]]])
    payload = (msgpack.packb(len(inner), use_bin_type=True)
               + lz4.block.compress(inner, store_size=False))
    artifact.write_bytes(msgpack.packb(msgpack.ExtType(99, payload),
                                       use_bin_type=True))
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert result["rows"] == []


def test_unreadable_artifact_is_reported_not_raised(tmp_path):
    ds = _dataset(tmp_path)
    artifact = tmp_path / "broken.arf"
    artifact.write_bytes(b"not msgpack at all")
    result = build_assay_evidence(ds, artifact, _adapter())
    assert result["availability"] is False
    assert "artifact_unreadable" in _reasons(result)


def test_unsupported_rt_unit_becomes_a_domain_error(tmp_path):
    ds = _dataset(tmp_path)
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0), _peak_row(1, "S2", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0), _peak_row(1, "S2", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    adapter = _adapter()
    adapter["evidence_rt_unit"] = "hour"
    with pytest.raises(DomainError) as caught:
        build_assay_evidence(ds, artifact, adapter)
    assert caught.value.code == "EVIDENCE_UNIT_UNSUPPORTED"


def test_unsupported_evidence_reader_is_a_domain_error(tmp_path):
    ds = _dataset(tmp_path)
    artifact = _write_arf(tmp_path, [
        [_peak_row(0, "S1", mz=700.5, rt=3.0), _peak_row(1, "S2", mz=700.5, rt=3.0),
         _peak_row(2, "S3", mz=700.5, rt=3.0)],
        [_peak_row(0, "S1", mz=701.5, rt=4.0), _peak_row(1, "S2", mz=701.5, rt=4.0),
         _peak_row(2, "S3", mz=701.5, rt=4.0)],
    ])
    adapter = _adapter()
    adapter["evidence_reader"] = "pai2"
    with pytest.raises(DomainError) as caught:
        build_assay_evidence(ds, artifact, adapter)
    assert caught.value.code == "EVIDENCE_READER_UNSUPPORTED"


# ---------- 全SME候補（Task 7 のbindingが読む） ----------

def test_all_sme_candidates_are_kept_with_library_identity_and_score():
    """rank最上位だけを残すと、bindingがadduct/charge/scoreで選び直せない。"""
    from metabolomix.mztab.dataset_state import build_dataset_state

    parse_result = {
        "metadata": {"assay[1]": "S1"},
        "sections": {
            "SMF": {"rows": [{"SMF_ID": "1", "SME_ID_REFS": "1|2",
                              "exp_mass_to_charge": "104.0706",
                              "retention_time_in_seconds": "60.0",
                              "abundance_assay[1]": "1000"}]},
            "SME": {"rows": [
                {"SME_ID": "1", "rank": "1", "chemical_name": "GABA",
                 "database_identifier": "HMDB:HMDB0000112",
                 "adduct_ion": "[M+H]1+", "charge": "1",
                 "identification_method": "matched to library",
                 "best_id_confidence_measure": "MS-DIAL total score",
                 "best_id_confidence_value": "0.92"},
                {"SME_ID": "2", "rank": "2", "chemical_name": "other",
                 "database_identifier": "HMDB:HMDB0000999",
                 "adduct_ion": "[M+Na]1+", "charge": "1",
                 "best_id_confidence_value": "0.41"},
            ]},
            "SML": {"rows": []},
        },
    }
    ds = build_dataset_state(parse_result, "Height_x.mzTab", "Height_x.mzTab")

    candidates = ds.feature_candidates["1"]
    assert [c["sme_id"] for c in candidates] == ["1", "2"]
    assert candidates[0]["rank"] == 1
    assert candidates[0]["database_identifier"] == "HMDB:HMDB0000112"
    assert candidates[0]["adduct"] == "[M+H]1+"
    assert candidates[0]["charge"] == 1
    assert candidates[0]["confidence_value"] == pytest.approx(0.92)
    assert candidates[1]["adduct"] == "[M+Na]1+"


def test_features_without_identification_keep_an_empty_candidate_list():
    from metabolomix.mztab.dataset_state import build_dataset_state

    parse_result = {
        "metadata": {"assay[1]": "S1"},
        "sections": {
            "SMF": {"rows": [{"SMF_ID": "1", "SME_ID_REFS": "",
                              "exp_mass_to_charge": "104.0706",
                              "abundance_assay[1]": "1000"}]},
            "SME": {"rows": []},
            "SML": {"rows": []},
        },
    }
    ds = build_dataset_state(parse_result, "Height_x.mzTab", "Height_x.mzTab")
    assert ds.feature_candidates["1"] == []
