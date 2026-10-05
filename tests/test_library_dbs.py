"""`.dbs` / `.lbm2` リーダ。fixture はテスト自身が組む（実データに依存しない）。"""
import struct
import zipfile

import lz4.block
import msgpack
import pytest

from lipidmix.library import dbs
from tests.dbs_fixture import pack_chunk as _pack_chunk
from tests.dbs_fixture import record as _record
from tests.dbs_fixture import write_dbs


def test_a_short_array_header_is_read_correctly(tmp_path):
    """要素数が小さいと array ヘッダが短形式になる。常に 5 バイトと仮定すると壊れる。"""
    chunk = _pack_chunk([_record()] * 3)
    dec = next(dbs.iter_decompressed_chunks(chunk))
    assert dbs.read_array_count(dec) == 3


def test_records_come_back_from_a_multi_chunk_lbm2(tmp_path):
    path = tmp_path / "lib.lbm2"
    path.write_bytes(_pack_chunk([_record(name="A")] * 2) + _pack_chunk([_record(name="B")]))

    records = list(dbs.iter_records(path))

    assert [r["name"] for r in records] == ["A", "A", "B"]
    assert records[0]["precursor_mz"] == pytest.approx(104.0706)
    assert records[0]["ion_mode"] == "negative"      # IonMode=1
    assert records[0]["adduct"] == "[M+H]+"
    assert records[0]["inchikey"] == "BTCSSZJGUNDROE-UHFFFAOYSA-N"
    assert records[0]["compound_class"] == "AminoAcid"
    assert records[0]["spectrum"] == [[87.04, 999.0], [69.03, 500.0]]


def test_a_dbs_zip_is_read_from_its_database_entry(tmp_path):
    path = tmp_path / "Project_Loaded.msp2.dbs"
    storage = {"MetabolomicsDataBases": [{"DataBase": ["mylib", 4, 2, "C:/x/mylib.lbm2"],
                                          "Pairs": [[0, {"SerializableAnnotatorKey":
                                                         [3, {"Parameter": [0.0, 2000.0, 2.0, 100.0, 20.0,
                                                                            0.01, 0.025, 0.0, 0.0, 0.0225,
                                                                            0.0225, 0.09, 0.0, 0.8, 1.0,
                                                                            True, True, False, False, 0.1],
                                                              "SourceType": 4, "Key": "mylib_1",
                                                              "Priority": 1}]}]]}]}
    packed = msgpack.packb(storage, use_bin_type=True)
    comp = lz4.block.compress(packed, store_size=False)
    wrapped = (b"\xc9" + struct.pack(">I", len(comp) + 5) + b"\x63"
               + b"\xd2" + struct.pack(">i", len(packed)) + comp)
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("MetabolomicsDB/mylib/DataBase", _pack_chunk([_record(name="Z")]))
        z.writestr("MetabolomicsDB/MS-FINDER/DataBase", b"")
        z.writestr("Storage", wrapped)

    assert [r["name"] for r in dbs.iter_records(path)] == ["Z"]

    meta = dbs.read_storage_meta(path)
    assert meta["library_name"] == "mylib"
    assert meta["source_path"] == "C:/x/mylib.lbm2"
    assert meta["search_params"]["ms2_tolerance"] == pytest.approx(0.025)
    assert meta["search_params"]["mass_range_end"] == pytest.approx(2000.0)
    assert meta["search_params"]["rt_tolerance"] == pytest.approx(2.0)
    # Key 7 / 8。採点前の足切りで、MS-DIAL の正規化の唯一の可変部分（Task 1）。
    assert meta["search_params"]["relative_amp_cutoff"] == pytest.approx(0.0)
    assert meta["search_params"]["absolute_amp_cutoff"] == pytest.approx(0.0)


def test_an_empty_msfinder_entry_is_not_mistaken_for_the_database(tmp_path):
    """MS-FINDER エントリは空で出る。これを DataBase と取り違えない。"""
    path = tmp_path / "P_Loaded.msp2.dbs"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("MetabolomicsDB/MS-FINDER/DataBase", b"")
        z.writestr("MetabolomicsDB/real/DataBase", _pack_chunk([_record(name="R")]))
    assert [r["name"] for r in dbs.iter_records(path)] == ["R"]


def test_a_proteomics_entry_is_not_read_as_a_metabolomics_record(tmp_path):
    """`ProteomicsDB/` / `EadLipidomicsDB/` は対象外。中身があっても拾わない
    （拾うと MoleculeMsReference の Key 配置とは限らずフィールドがずれる）。"""
    path = tmp_path / "P_Loaded.msp2.dbs"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("MetabolomicsDB/real/DataBase", _pack_chunk([_record(name="R")]))
        z.writestr("ProteomicsDB/other/DataBase", _pack_chunk([_record(name="P")]))
        z.writestr("EadLipidomicsDB/other/DataBase", _pack_chunk([_record(name="E")]))
    assert [r["name"] for r in dbs.iter_records(path)] == ["R"]


def test_rt_is_matched_by_chromxs_type_code_not_by_position(tmp_path):
    """ChromXs は `[[種別, [値, ...]], ...]`。種別 1 が RT、3 が m/z など。
    先頭要素が RT とは限らないので、位置ではなく種別コードで一致を取る。
    実データでは末尾に入れ子でないスカラーが付くこともあり、それは無視する。"""
    r = _record()
    r[2] = [[3, [104.07, 0, 0]], [1, [1.23, 0, 0]], 0, 2.4764, 0, 0]
    path = tmp_path / "lib.lbm2"
    path.write_bytes(_pack_chunk([r]))
    records = list(dbs.iter_records(path))
    assert records[0]["rt"] == pytest.approx(1.23)


def test_rt_is_none_for_the_unset_sentinel(tmp_path):
    """`-1` は MS-DIAL の「未設定」番兵。RT として -1 をそのまま返さない。"""
    r = _record()
    r[2] = [[1, [-1.0, 0, 0]]]
    path = tmp_path / "lib.lbm2"
    path.write_bytes(_pack_chunk([r]))
    records = list(dbs.iter_records(path))
    assert records[0]["rt"] is None


def _write_dbs(path, parameter):
    """`Storage` だけを差し替えられる最小の `.dbs`（fixture はテスト自身が作る規約）。"""
    storage = {"MetabolomicsDataBases": [{"DataBase": ["mylib", 4, 2, "C:/x/mylib.lbm2"],
                                          "Pairs": [[0, {"SerializableAnnotatorKey":
                                                         [3, {"Parameter": parameter,
                                                              "SourceType": 4, "Key": "mylib_1",
                                                              "Priority": 1}]}]]}]}
    packed = msgpack.packb(storage, use_bin_type=True)
    comp = lz4.block.compress(packed, store_size=False)
    wrapped = (b"\xc9" + struct.pack(">I", len(comp) + 5) + b"\x63"
               + b"\xd2" + struct.pack(">i", len(packed)) + comp)
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("MetabolomicsDB/mylib/DataBase", _pack_chunk([_record(name="Z")]))
        z.writestr("Storage", wrapped)
    return path


# Key 0..19。15〜18 が注釈スコアリングの可否（既定は全部 False）。
_PARAMETER = [0.0, 2000.0, 2.0, 100.0, 20.0, 0.01, 0.025, 0.0, 0.0, 0.0225,
              0.0225, 0.09, 0.0, 0.8, 1.0, True, True, False, False, 0.1]


def test_annotation_scoring_flags_are_read_from_the_storage(tmp_path):
    """Key 15〜18。総合スコアに RT / CCS 項を入れるかはここで決まる。
    **既定は全部 False だが実 run では True のことがある**ので、推測せず読む。"""
    meta = dbs.read_storage_meta(_write_dbs(tmp_path / "on.msp2.dbs", _PARAMETER))

    assert meta["search_params"]["use_time_for_annotation_filtering"] is True
    assert meta["search_params"]["use_time_for_annotation_scoring"] is True
    assert meta["search_params"]["use_ccs_for_annotation_filtering"] is False
    assert meta["search_params"]["use_ccs_for_annotation_scoring"] is False


def test_the_scoring_flags_are_read_and_not_hardcoded(tmp_path):
    off = list(_PARAMETER)
    off[15] = off[16] = False
    off[17] = off[18] = True
    meta = dbs.read_storage_meta(_write_dbs(tmp_path / "off.msp2.dbs", off))

    assert meta["search_params"]["use_time_for_annotation_scoring"] is False
    assert meta["search_params"]["use_ccs_for_annotation_scoring"] is True


def test_the_storage_maps_every_annotator_key_to_its_library(tmp_path):
    """照合結果の AnnotatorID はライブラリ名と同じとは限らない。MS-DIAL Console は
    LBM ならファイルのパス（後の版は `LbmDB: <stem>`）、MSP なら `.msp` のパスや
    設定ファイルで付けた任意の名前を使う。`Storage` が保存している Key → ライブラリ名の
    対応を、全ライブラリ・全注釈器について読む。"""
    path = write_dbs(tmp_path / "P_Loaded.msp2.dbs", [
        ("LbmDB", [_record(name="L")], ["C:/lib/NCDK_conventional.lbm2"]),
        ("MspDB_lab_1", [_record(name="M")], ["lab-strict", "lab-loose"]),
        ("TextDB", [_record(name="T")], ["C:/lib/targets.txt"]),
    ])

    meta = dbs.read_storage_meta(path)

    assert meta["annotator_libraries"] == {
        "C:/lib/NCDK_conventional.lbm2": "LbmDB",
        "lab-strict": "MspDB_lab_1", "lab-loose": "MspDB_lab_1",
        "C:/lib/targets.txt": "TextDB",
    }
    assert meta["library_name"] == "LbmDB"   # 先頭ライブラリの扱いは変えない


def test_the_storage_keeps_search_params_per_annotator(tmp_path):
    """RT を使うかは注釈器ごとに違いうる（先頭の注釈器だけで代表させない）。"""
    off = list(_PARAMETER)
    off[15] = off[16] = False
    path = tmp_path / "P_Loaded.msp2.dbs"
    write_dbs(path, [("LbmDB", [_record(name="L")], ["C:/lib/a.lbm2"])], parameter=off)
    meta = dbs.read_storage_meta(path)
    params = meta["annotator_search_params"]["C:/lib/a.lbm2"]
    assert params["use_time_for_annotation_filtering"] is False
    assert params["use_time_for_annotation_scoring"] is False
    assert params["rt_tolerance"] == pytest.approx(2.0)
