"""MS-DIAL バイナリの復号経路を、テスト自身が合成したバイト列だけで縛る。

実データ（`.arf` `.arf2` `.pai2` `.dcl` `.EIC.aef`）は追跡外で巨大なため
リポジトリに置けない。その結果、復号本体はこれまでほぼ未実行のまま残っていた
（実測で `_convert_to_peakfeature` 3%、`parse_eic_aef_css1` 2%、
`deserialize_lz4_packed_msgpack` 9〜17%）。**MessagePack の Key 番号を 1 つずらしても
テストは全部緑のまま**という状態であり、CLAUDE.md が「推測で直すな」と最も強く
警告している箇所に自動検知が無かった。

ここでは msgpack + LZ4 のコンテナと CSS1 バイナリをテスト側で組み立て、
`docs/schema/*.md`（MS-DIAL の `[Key(N)]` から抽出した Key 番号表）から読んだ Key 番号に
番兵値を置いて、リーダが**同じ番号から**拾うことを検証する。Key 番号の正解は
実装ではなくスキーマ文書から取るため、実装をなぞるだけのテストにはならない。

「fixture はテスト自身が作る」規約に沿い、実ファイルには一切依存しない。
"""
import io
import re
import struct
import unittest
from pathlib import Path

import lz4.block
import msgpack

from metabolomix.arf import reader as arf_reader
from metabolomix.arf2 import reader as arf2_reader
from metabolomix.dcl import reader as dcl_reader
from metabolomix.eic import reader as eic_reader
from metabolomix.pai2 import reader as pai2_reader

REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_DIR = REPO_ROOT / "docs" / "schema"


# --------------------------------------------------------------------------
# docs/schema/*.md から Key 番号 → プロパティ名を読む
# --------------------------------------------------------------------------
def schema_keys(doc: str, cls: str) -> dict:
    """`docs/schema/*.md` の Key 番号表から Key 番号 → メンバ名の対応を取る。

    表は `| Key | 型 | メンバ | 備考 |` の行が並ぶ形。1 ファイルに複数クラスの表が
    入っており（`AlignmentSpotVariableCorrelation` `LinkedPeakFeature`）、**後続クラスも
    Key(0) から振り直す**ため、`## 同ファイル内の別クラス: <名前>` 見出しでクラスの範囲を
    切らないと Key(0) が別クラスの値で上書きされる。主クラスの表はファイル冒頭にあり、
    クラス名はファイル名と同じ。

    「欠番」表（`| Key | 状態 |`）と入れ子型の表（1 列目が型名）は列の形が違うため、
    「1 列目が数字・3 列目がバッククォート付きの識別子」という条件で自然に外れる。
    """
    text = (SCHEMA_DIR / f"{doc}.md").read_text(encoding="utf-8")

    keys: dict = {}
    current = doc  # 主クラスの表はファイル名と同じクラス名で始まる
    for line in text.splitlines():
        heading = re.match(r"##\s*同ファイル内の別クラス:\s*(\w+)", line)
        if heading:
            current = heading.group(1)
            continue
        if current != cls:
            continue
        row = re.match(r"\|\s*(\d+)\s*\|\s*`[^`]+`\s*\|\s*`(\w+)`\s*\|", line)
        if row:
            keys.setdefault(int(row.group(1)), row.group(2))

    if not keys:
        raise AssertionError(
            f"{doc}.md から class {cls} の Key 番号表を読めなかった。"
            "表の形（| Key | 型 | メンバ | 備考 |）が変わっていないか確認する。"
        )
    return keys


SPOT_KEYS = schema_keys("AlignmentSpotProperty", "AlignmentSpotProperty")
PEAK_KEYS = schema_keys("ChromatogramPeakFeature", "ChromatogramPeakFeature")
ROW_KEYS = schema_keys("AlignmentChromPeakFeature", "AlignmentChromPeakFeature")


# --------------------------------------------------------------------------
# バイナリの組み立て
# --------------------------------------------------------------------------
def at_keys(size: int, values: dict) -> list:
    """`values` の Key 番号の位置にだけ値を置いた、長さ `size` の行を作る。"""
    row = [None] * size
    for key, value in values.items():
        row[key] = value
    return row


def chromxs(rt=None, ri=None, mz=None, dt=None) -> list:
    """ChromXs（Key4 TimesCenter / Key15 ChromXsTop）の入れ子表現。

    `[[種別, [値, ...]], ...]` で種別は 1=RT 2=RI 3=m/z 4=ドリフト時間。
    """
    out = []
    for kind, value in ((1, rt), (2, ri), (3, mz), (4, dt)):
        if value is not None:
            out.append([kind, [value]])
    return out


def lz4_container(inner: bytes, header="hdr") -> bytes:
    """MS-DIAL の外側コンテナ: msgpack([header, msgpack(非圧縮長) + LZ4 ブロック])。"""
    payload = msgpack.packb(len(inner), use_bin_type=True) + lz4.block.compress(
        inner, store_size=False
    )
    return msgpack.packb([header, payload], use_bin_type=True)


def msgpack_stream(objects) -> bytes:
    """msgpack オブジェクトを連結したストリーム（`.arf` / `.arf2` の内側）。"""
    return b"".join(msgpack.packb(obj, use_bin_type=True) for obj in objects)


def css1_bytes(spots: list) -> bytes:
    """`.EIC.aef`（CSS1）を組み立てる。

    レイアウト: 'CSS1' + 6 バイト予約 + スポット数(i) + シークポインタ(q)×N、
    各スポットは RT/RI/Mass/Drift(f×4) + MainType(b) + サンプル数(i)、
    各サンプルは FileID(i)/ピーク数(i)/Top/Left/Right(f×3) + (横軸(f), 強度(f))×ピーク数。
    """
    header = b"CSS1" + b"\x00" * 6 + struct.pack("<i", len(spots))
    table_size = 8 * len(spots)
    body = b""
    pointers = []
    for spot in spots:
        pointers.append(len(header) + table_size + len(body))
        body += struct.pack(
            "<ffffbi",
            spot["rt"], spot["ri"], spot["mz"], spot["drift"],
            spot["main_type"], len(spot["samples"]),
        )
        for sample in spot["samples"]:
            points = sample["points"]
            body += struct.pack(
                "<iifff",
                sample["file_id"], len(points),
                sample["top"], sample["left"], sample["right"],
            )
            for horiz, intensity in points:
                body += struct.pack("<ff", horiz, intensity)
    return header + b"".join(struct.pack("<q", p) for p in pointers) + body


# --------------------------------------------------------------------------
# スキーマ文書の読み取り自体を縛る
# --------------------------------------------------------------------------
class TestSchemaKeyExtraction(unittest.TestCase):
    """Key 番号の正解表を実装ではなくスキーマ文書から取る、という前提を固定する。"""

    def test_alignment_spot_property_keys(self):
        expected = {
            0: "MasterAlignmentID", 1: "AlignmentID", 4: "TimesCenter",
            5: "MassCenter", 11: "IonMode", 12: "Name", 13: "Formula",
            14: "Ontology", 15: "SMILES", 16: "InChIKey", 31: "HeightAverage",
            32: "HeightMin", 33: "HeightMax", 43: "MassMin", 44: "MassMax",
            # MS-DIAL 側の綴り（Percentage ではない）。直さないこと。
            49: "FillParcentage", 51: "MonoIsotopicPercentage", 54: "AdductType",
        }
        for key, name in expected.items():
            with self.subTest(key=key):
                self.assertEqual(SPOT_KEYS.get(key), name)

    def test_chromatogram_peak_feature_keys(self):
        expected = {
            3: "ChromXsLeft", 4: "ChromXsTop", 5: "ChromXsRight",
            6: "PeakHeightLeft", 7: "PeakHeightTop", 8: "PeakHeightRight",
            9: "PeakAreaAboveZero", 10: "PeakAreaAboveBaseline",
            11: "MasterPeakID", 18: "MS2RawSpectrumID", 19: "MS2RawSpectrumID2CE",
            22: "IonMode", 24: "Spectrum", 25: "Name", 26: "Formula",
            27: "Ontology", 28: "SMILES", 29: "InChIKey", 30: "AdductType",
            31: "CollisionCrossSection", 38: "Comment", 40: "PeakShape", 43: "Mass",
        }
        for key, name in expected.items():
            with self.subTest(key=key):
                self.assertEqual(PEAK_KEYS.get(key), name)

    def test_alignment_chrom_peak_feature_keys(self):
        expected = {
            0: "FileID", 2: "MasterPeakID", 3: "PeakID", 9: "MS2RawSpectrumID",
            10: "MS2RawSpectrumID2CE", 15: "ChromXsTop", 18: "PeakHeightTop",
            20: "PeakAreaAboveZero", 21: "PeakAreaAboveBaseline", 22: "Mass",
            23: "IonMode", 24: "Name", 37: "PeakShape",
        }
        for key, name in expected.items():
            with self.subTest(key=key):
                self.assertEqual(ROW_KEYS.get(key), name)

    def test_second_class_in_the_same_file_does_not_leak(self):
        # AlignmentSpotProperty.md には AlignmentSpotVariableCorrelation が続く。
        # 範囲を切らないと Key(0) が CorrelateAlignmentID で上書きされる。
        self.assertNotEqual(SPOT_KEYS.get(0), "CorrelateAlignmentID")


# --------------------------------------------------------------------------
# .arf2 = AlignmentSpotProperty（スポット代表）
# --------------------------------------------------------------------------
def arf2_spot(overrides: dict | None = None) -> list:
    values = {
        0: 7, 1: 3,
        4: chromxs(rt=12.5, mz=700.5),
        5: 700.5001, 11: 1, 12: "PC 34:1",
        13: ["C42H82NO8P", 759.5778],
        14: "PC", 15: "CCCCC", 16: "SUBSTITUTE-KEY-A",
        31: 12345.0, 32: 100.0, 33: 20000.0,
        34: 0.25, 35: 30.0, 36: 55.0, 37: 12.0,
        43: 700.4, 44: 700.6, 49: 0.75, 51: 0.9,
        54: [1.00794, 1, "[M+H]+"],
    }
    values.update(overrides or {})
    return at_keys(60, values)


class TestArf2Decode(unittest.TestCase):
    def one(self, overrides: dict | None = None) -> dict:
        blob = lz4_container(msgpack_stream([arf2_spot(overrides)]))
        spots = arf2_reader.deserialize(io.BytesIO(blob))
        self.assertEqual(len(spots), 1)
        return spots[0]

    def test_container_roundtrips_multiple_spots(self):
        blob = lz4_container(msgpack_stream([arf2_spot(), arf2_spot({0: 8})]))
        spots = arf2_reader.deserialize(io.BytesIO(blob))
        self.assertEqual([s["MasterAlignmentID"] for s in spots], [7, 8])

    def test_values_land_in_the_fields_named_by_the_schema(self):
        spot = self.one()
        self.assertEqual(spot["MasterAlignmentID"], 7)          # Key0
        self.assertEqual(spot["AlignmentID"], 3)                # Key1
        self.assertAlmostEqual(spot["RT"], 12.5)                # Key4 TimesCenter
        self.assertAlmostEqual(spot["MassCenter"], 700.5001)    # Key5
        self.assertEqual(spot["IonMode"], "Negative")           # Key11 (1=Negative)
        self.assertEqual(spot["Name"], "PC 34:1")               # Key12
        self.assertEqual(spot["Ontology"], "PC")                # Key14
        self.assertEqual(spot["SMILES"], "CCCCC")               # Key15
        self.assertEqual(spot["InChIKey"], "SUBSTITUTE-KEY-A")  # Key16
        self.assertAlmostEqual(spot["HeightAverage"], 12345.0)  # Key31
        self.assertAlmostEqual(spot["HeightMin"], 100.0)        # Key32
        self.assertAlmostEqual(spot["HeightMax"], 20000.0)      # Key33
        self.assertAlmostEqual(spot["MassMin"], 700.4)          # Key43
        self.assertAlmostEqual(spot["MassMax"], 700.6)          # Key44
        self.assertAlmostEqual(spot["FillPercentage"], 0.75)         # Key49
        self.assertAlmostEqual(spot["MonoIsotopicPercentage"], 0.9)  # Key51

    def test_formula_and_adduct_come_from_inside_the_nested_object(self):
        # Key13 Formula = [組成式, 質量, ...] の先頭
        # Key54 AdductType = [質量差, charge, 表記, ...] の index 2
        spot = self.one()
        self.assertEqual(spot["Formula"], "C42H82NO8P")
        self.assertEqual(spot["AdductType"], "[M+H]+")

    def test_missing_nested_objects_fall_back_instead_of_raising(self):
        spot = self.one({13: None, 54: None})
        self.assertEqual(spot["Formula"], "")
        self.assertEqual(spot["AdductType"], "")

    def test_non_numeric_value_falls_back_to_zero(self):
        # _to_float は TypeError/ValueError を既定値に倒す。壊れた行で全体を落とさない。
        self.assertEqual(self.one({31: "N/A"})["HeightAverage"], 0.0)

    def test_ion_mode_maps_by_number(self):
        for value, name in ((0, "Positive"), (1, "Negative"), (2, "Both"), (99, "Unknown")):
            with self.subTest(value=value):
                self.assertEqual(self.one({11: value})["IonMode"], name)

    def test_empty_name_becomes_unknown(self):
        self.assertEqual(self.one({12: ""})["Name"], "Unknown")

    def test_times_center_reads_only_its_own_kind(self):
        # 種別 3 は m/z であって RT ではない。取り違えると RT が m/z になる。
        spot = self.one({4: chromxs(mz=700.5)})
        self.assertEqual(spot["RT"], 0.0)

    def test_short_lists_are_treated_as_containers_of_spots(self):
        # deserialize は len<14 の項目を「スポット群を内包する箱」として展開する
        blob = lz4_container(msgpack_stream([[[arf2_spot(), arf2_spot({0: 9})]]]))
        spots = arf2_reader.deserialize(io.BytesIO(blob))
        self.assertEqual([s["MasterAlignmentID"] for s in spots], [7, 9])

    def test_too_short_rows_are_dropped(self):
        blob = lz4_container(msgpack_stream([[1, 2, 3]]))
        self.assertEqual(arf2_reader.deserialize(io.BytesIO(blob)), [])


class TestArf2Summary(unittest.TestCase):
    def test_summary_counts_annotated_spots(self):
        spots = [
            arf2_reader.extract_arf2_data(arf2_spot()),
            arf2_reader.extract_arf2_data(arf2_spot({12: "Unknown"})),
        ]
        summary = arf2_reader.summarize_arf2_data(spots)
        self.assertEqual(summary["total_spots"], 2)
        self.assertEqual(summary["annotated_count"], 1)
        self.assertAlmostEqual(summary["annotation_rate"], 50.0)
        self.assertEqual(summary["ontology_top"], {"PC": 2})

    def test_summary_of_nothing_reports_an_error(self):
        self.assertIn("error", arf2_reader.summarize_arf2_data([]))

    def test_text_summary_mentions_the_counts(self):
        spots = [arf2_reader.extract_arf2_data(arf2_spot())]
        text = arf2_reader.generate_text_summary(spots)
        self.assertIn("ARF2 カタログ要約", text)
        self.assertIn("PC", text)


# --------------------------------------------------------------------------
# .pai2 = ChromatogramPeakFeature（個別測定のピーク）
# --------------------------------------------------------------------------
def pai2_feature(overrides: dict | None = None) -> list:
    values = {
        3: chromxs(rt=11.9), 4: chromxs(rt=12.0), 5: chromxs(rt=12.1),
        6: 800.0, 7: 9000.0, 8: 700.0,
        9: 45000.0, 10: 43000.0,
        11: 42, 18: 17, 19: {0: 20.0, 1: 35.0},
        22: 1, 24: [[184.0733, 999.0], [104.1075, 120.0]],
        25: "PC 34:1", 26: ["C42H82NO8P", 759.5778],
        27: "PC", 28: "CCCCC", 29: "SUBSTITUTE-KEY-A",
        30: [1.00794, 1, "[M+H]+"], 31: 250.5,
        38: "note", 40: [12.0, 55.5], 43: 760.5851,
    }
    values.update(overrides or {})
    return at_keys(50, values)


class TestPai2Decode(unittest.TestCase):
    def one(self, overrides: dict | None = None) -> dict:
        blob = lz4_container(msgpack.packb([pai2_feature(overrides)], use_bin_type=True))
        features = pai2_reader.deserialize(io.BytesIO(blob))
        self.assertEqual(len(features), 1)
        return features[0]

    def test_values_land_in_the_fields_named_by_the_schema(self):
        feature = self.one()
        self.assertAlmostEqual(feature["time"]["rt"], 12.0)          # Key4 ChromXsTop
        self.assertAlmostEqual(feature["time_left"]["rt"], 11.9)     # Key3 ChromXsLeft
        self.assertAlmostEqual(feature["time_right"]["rt"], 12.1)    # Key5 ChromXsRight
        self.assertAlmostEqual(feature["peak_height"], 9000.0)       # Key7
        self.assertAlmostEqual(feature["peak_height_left"], 800.0)   # Key6
        self.assertAlmostEqual(feature["peak_height_right"], 700.0)  # Key8
        self.assertAlmostEqual(feature["peak_area"], 45000.0)                 # Key9
        self.assertAlmostEqual(feature["peak_area_above_baseline"], 43000.0)  # Key10
        self.assertEqual(feature["id"], 42)                          # Key11 MasterPeakID
        self.assertEqual(feature["ion_mode"], pai2_reader.IonMode.Negative)   # Key22
        self.assertEqual(feature["name"], "PC 34:1")                 # Key25
        self.assertEqual(feature["ontology"], "PC")                  # Key27
        self.assertEqual(feature["smiles"], "CCCCC")                 # Key28
        self.assertEqual(feature["inchikey"], "SUBSTITUTE-KEY-A")    # Key29
        self.assertAlmostEqual(feature["collision_cross_section"], 250.5)  # Key31
        self.assertEqual(feature["comment"], "note")                 # Key38
        self.assertAlmostEqual(feature["m/z"], 760.5851)             # Key43 Mass

    def test_signal_to_noise_is_the_second_element_of_peak_shape(self):
        # Key40 PeakShape = [EstimatedNoise, SignalToNoise, ...]。先頭を取ると
        # ノイズ推定値を S/N として報告してしまう。
        self.assertAlmostEqual(self.one()["S/N"], 55.5)

    def test_formula_and_adduct_come_from_inside_the_nested_object(self):
        feature = self.one()
        self.assertEqual(feature["formula"], "C42H82NO8P")   # Key26 の先頭
        self.assertEqual(feature["adduct"], "[M+H]+")        # Key30 の index 2

    def test_collision_energies_are_deduplicated_and_sorted(self):
        feature = self.one({19: {0: 35.0, 1: 20.0, 2: 35.0}})
        self.assertEqual(feature["collision_energies"], [20.0, 35.0])

    def test_has_msms_is_true_when_a_collision_energy_map_exists(self):
        self.assertTrue(self.one()["has_msms"])

    def test_has_msms_is_true_on_a_raw_id_alone(self):
        self.assertTrue(self.one({19: {}, 18: 5})["has_msms"])

    def test_has_msms_is_false_without_either(self):
        feature = self.one({19: {}, 18: -1})
        self.assertFalse(feature["has_msms"])
        self.assertEqual(feature["collision_energies"], [])

    def test_spectrum_peaks_are_trimmed_to_mass_and_intensity(self):
        feature = self.one({24: [[184.0733, 999.0, "extra"], [104.1075, 120.0]]})
        self.assertEqual(feature["msms_peak_count"], 2)
        self.assertEqual(feature["msms_spectrum"], [[184.0733, 999.0], [104.1075, 120.0]])

    def test_negative_times_are_ignored(self):
        # MS-DIAL は未設定の軸に負値を入れる。そのまま採ると RT が負になる。
        self.assertEqual(self.one({4: [[1, [-1.0]]]})["time"], {})


# --------------------------------------------------------------------------
# .arf = AlignmentChromPeakFeature の行がスポットごとに束ねられたもの
# --------------------------------------------------------------------------
def arf_row(overrides: dict | None = None) -> list:
    values = {
        0: 3, 1: "Sample_A.abf", 2: 11, 3: 5,
        9: 17, 10: {0: 20.0},
        15: chromxs(rt=12.0), 18: 9000.0, 20: 45000.0, 21: 43000.0,
        22: 760.5851, 23: 1, 24: "PC 34:1", 37: [12.0, 55.5],
    }
    values.update(overrides or {})
    return at_keys(50, values)


def arf_group(n_rows: int = 3, first: dict | None = None) -> list:
    rows = [arf_row(first)]
    rows += [arf_row({0: i, 1: "Sample_" + str(i) + ".abf"}) for i in range(1, n_rows)]
    return rows


class TestArfBlockDecode(unittest.TestCase):
    def test_container_object_yields_its_groups(self):
        container = [0, 1, arf_group(), arf_group()]
        blob = lz4_container(msgpack_stream([container]))
        spots = arf_reader.deserialize(io.BytesIO(blob))
        self.assertEqual([s["MasterAlignmentID"] for s in spots], [0, 1])
        self.assertEqual([s["SourceLocalIndex"] for s in spots], [0, 1])

    def test_bare_group_outside_a_container_is_still_read(self):
        blob = lz4_container(msgpack_stream([arf_group()]))
        self.assertEqual(len(arf_reader.deserialize(io.BytesIO(blob))), 1)

    def test_ext_type_payload_is_accepted(self):
        inner = msgpack_stream([arf_group()])
        payload = msgpack.packb(len(inner), use_bin_type=True) + lz4.block.compress(
            inner, store_size=False
        )
        blob = msgpack.packb(msgpack.ExtType(1, payload), use_bin_type=True)
        self.assertEqual(len(arf_reader.deserialize(io.BytesIO(blob))), 1)

    def test_unsupported_top_level_object_is_skipped_not_fatal(self):
        good = lz4_container(msgpack_stream([arf_group()]))
        blob = msgpack.packb("メタデータ", use_bin_type=True) + good
        self.assertEqual(len(arf_reader.deserialize(io.BytesIO(blob))), 1)

    def test_groups_shorter_than_three_rows_are_dropped(self):
        blob = lz4_container(msgpack_stream([[arf_row(), arf_row()]]))
        self.assertEqual(arf_reader.deserialize(io.BytesIO(blob)), [])

    def test_representative_row_supplies_the_spot_fields(self):
        blob = lz4_container(msgpack_stream([arf_group()]))
        spot = arf_reader.deserialize(io.BytesIO(blob))[0]
        self.assertAlmostEqual(spot["RT"], 12.0)               # Key15 ChromXsTop
        self.assertAlmostEqual(spot["HeightAverage"], 9000.0)  # Key18
        self.assertAlmostEqual(spot["MassCenter"], 760.5851)   # Key22
        self.assertEqual(spot["IonMode"], "Negative")          # Key23
        self.assertEqual(spot["Name"], "PC 34:1")              # Key24
        self.assertEqual(len(spot["AlignedPeakProperties"]), 3)


class TestArfRowConversion(unittest.TestCase):
    def row(self, overrides: dict | None = None) -> dict:
        return arf_reader._convert_to_alignment_feature(arf_row(overrides))

    def test_values_land_in_the_fields_named_by_the_schema(self):
        row = self.row()
        self.assertEqual(row["file_id"], 3)                          # Key0 FileID
        self.assertEqual(row["master_peak_id"], 11)                  # Key2
        self.assertEqual(row["peak_id"], 5)                          # Key3
        self.assertEqual(row["ms2_raw_id"], 17)                      # Key9
        self.assertAlmostEqual(row["height"], 9000.0)                # Key18
        self.assertAlmostEqual(row["area"], 45000.0)                 # Key20
        self.assertAlmostEqual(row["area_above_baseline"], 43000.0)  # Key21
        self.assertAlmostEqual(row["m_z"], 760.5851)                 # Key22
        self.assertAlmostEqual(row["rt"], 12.0)                      # Key15
        self.assertEqual(row["file_name"], "Sample_A.abf")

    def test_peak_shape_supplies_noise_and_signal_to_noise(self):
        # Key37 PeakShape = [EstimatedNoise, SignalToNoise]。順序を取り違えない。
        row = self.row()
        self.assertAlmostEqual(row["estimated_noise"], 12.0)
        self.assertAlmostEqual(row["signal_to_noise"], 55.5)

    def test_negative_master_peak_id_marks_a_gap_filled_value(self):
        # ギャップフィル＝未検出を補間した値。検出ピークと同一視すると偽陽性になる。
        self.assertTrue(self.row({2: -2})["is_gap_filled"])
        self.assertFalse(self.row({2: 11})["is_gap_filled"])

    def test_msms_flag_follows_the_collision_energy_map(self):
        self.assertTrue(self.row()["is_msms"])
        self.assertFalse(self.row({10: {}})["is_msms"])


# --------------------------------------------------------------------------
# .EIC.aef = CSS1（msgpack ではない独自バイナリ）
# --------------------------------------------------------------------------
def eic_spot(spot_id=0, rt=12.0, mz=760.5, samples=None) -> dict:
    return {
        "rt": rt, "ri": 0.0, "mz": mz, "drift": 0.0, "main_type": 1,
        "samples": samples if samples is not None else [
            {"file_id": 0, "top": rt, "left": rt - 0.1, "right": rt + 0.1,
             "points": [(rt - 0.1, 10.0), (rt, 900.0), (rt + 0.1, 20.0)]},
        ],
    }


class TestEicCss1Decode(unittest.TestCase):
    def parse(self, spots, include_chromatogram=False):
        path = Path(self.tmp.name) / "synthetic.EIC.aef"
        path.write_bytes(css1_bytes(spots))
        return eic_reader.parse_eic_aef_css1(str(path), include_chromatogram)

    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_header_and_pointer_table_roundtrip(self):
        parsed = self.parse([eic_spot(0, rt=1.0, mz=100.0), eic_spot(1, rt=9.0, mz=900.0)])
        self.assertEqual([s["spot_id"] for s in parsed], [0, 1])
        self.assertAlmostEqual(parsed[1]["rt"], 9.0, places=4)
        self.assertAlmostEqual(parsed[1]["mz"], 900.0, places=3)

    def test_a_file_without_the_css1_magic_is_rejected(self):
        path = Path(self.tmp.name) / "bad.EIC.aef"
        path.write_bytes(b"NOPE" + b"\x00" * 20)
        with self.assertRaises(ValueError):
            eic_reader.parse_eic_aef_css1(str(path))

    def test_sample_block_reports_max_and_mean_intensity(self):
        parsed = self.parse([eic_spot()])
        sample = parsed[0]["samples"][0]
        self.assertEqual(sample["num_peaks"], 3)
        self.assertAlmostEqual(sample["max_intensity"], 900.0, places=3)
        self.assertAlmostEqual(sample["mean_intensity"], 310.0, places=3)

    def test_peak_top_is_a_horizontal_coordinate_not_an_intensity(self):
        # peak_top は RT 座標。強度と取り違えると「遅い RT の一覧」を強度上位と誤認する。
        parsed = self.parse([eic_spot(rt=12.0)])
        self.assertAlmostEqual(parsed[0]["samples"][0]["peak_top"], 12.0, places=4)

    def test_chromatogram_points_are_omitted_unless_requested(self):
        self.assertNotIn("chromatogram", self.parse([eic_spot()])[0]["samples"][0])
        with_points = self.parse([eic_spot()], include_chromatogram=True)
        self.assertEqual(len(with_points[0]["samples"][0]["chromatogram"]), 3)

    def test_summary_aggregates_across_spots_and_samples(self):
        parsed = self.parse([eic_spot(0, rt=1.0, mz=100.0), eic_spot(1, rt=9.0, mz=900.0)])
        summary = eic_reader.summarize_eic_data(parsed)
        self.assertEqual(summary["total_spots"], 2)
        self.assertEqual(summary["total_samples"], 2)
        self.assertEqual(summary["total_peaks"], 6)
        self.assertEqual(summary["unique_file_ids"], [0])
        self.assertAlmostEqual(summary["rt_range"][0], 1.0, places=4)
        self.assertAlmostEqual(summary["rt_range"][1], 9.0, places=4)

    def test_summary_of_nothing_is_zeroed_not_an_error(self):
        summary = eic_reader.summarize_eic_data([])
        self.assertEqual(summary["total_spots"], 0)
        self.assertIsNone(summary["rt_range"])

    def test_range_searches_filter_inclusively(self):
        parsed = self.parse([eic_spot(0, rt=1.0, mz=100.0), eic_spot(1, rt=9.0, mz=900.0)])
        self.assertEqual(
            [s["spot_id"] for s in eic_reader.search_eic_by_mz_range(parsed, 50.0, 500.0)], [0]
        )
        self.assertEqual(
            [s["spot_id"] for s in eic_reader.search_eic_by_rt_range(parsed, 5.0, 10.0)], [1]
        )

    def test_ranking_uses_intensity_not_retention_time(self):
        """旧実装は peak_top（RT 座標）で並べており、実質「遅い RT 順」だった。

        遅く溶出するが弱いスポットと、早く溶出するが強いスポットを並べ、
        強い方が先頭に来ることで回帰を防ぐ。
        """
        weak_late = eic_spot(0, rt=18.0, mz=100.0, samples=[
            {"file_id": 0, "top": 18.0, "left": 17.9, "right": 18.1,
             "points": [(18.0, 5.0)]},
        ])
        strong_early = eic_spot(1, rt=2.0, mz=200.0, samples=[
            {"file_id": 0, "top": 2.0, "left": 1.9, "right": 2.1,
             "points": [(2.0, 5000.0)]},
        ])
        ranked = eic_reader.top_eic_spots_by_max_intensity(self.parse([weak_late, strong_early]))
        self.assertEqual([r["spot_id"] for r in ranked], [1, 0])
        self.assertAlmostEqual(ranked[0]["max_intensity"], 5000.0, places=2)

    def test_truncated_stream_is_reported_with_context(self):
        with self.assertRaises(ValueError) as caught:
            eic_reader._read_exact(io.BytesIO(b"ab"), 8, "spot header")
        self.assertIn("spot header", str(caught.exception))


# --------------------------------------------------------------------------
# .dcl = MSDecResult（msgpack ではなく BinaryWriter 由来の独自バイナリ）
# --------------------------------------------------------------------------
def dcl_bytes(results: list) -> bytes:
    """`.dcl` を組み立てる。レイアウトは metabolomix/dcl/reader.py の冒頭に準拠。"""
    header = b"DC" + struct.pack("<i", 1) + b"\x00" + struct.pack("<i", len(results))
    table_size = 8 * len(results)
    body = b""
    pointers = []
    for result in results:
        pointers.append(len(header) + table_size + len(body))
        spectrum = result["spectrum"]
        body += struct.pack(
            "<qiididdd d",
            0, result["scan_id"], result["raw_spec_id"], result["precursor_mz"],
            result["ion_mode"], result["rt"], 0.0, 0.0, result["precursor_mz"],
        )
        body += struct.pack("<5d", 0.0, result["model_height"], 0.0, 0.0, 0.0)
        body += struct.pack("<5f", 0.0, 0.0, 0.0, result["sn"], 0.0)
        body += struct.pack("<3i", len(spectrum), 0, 0)
        for mass, intensity in spectrum:
            body += struct.pack("<2di", mass, intensity, 0)
    return header + b"".join(struct.pack("<q", p) for p in pointers) + body


def dcl_result(precursor_mz=760.5851, rt=12.0, spectrum=None, scan_id=1,
               raw_spec_id=17, ion_mode=1, model_height=9000.0, sn=55.5) -> dict:
    return {
        "scan_id": scan_id, "raw_spec_id": raw_spec_id, "precursor_mz": precursor_mz,
        "ion_mode": ion_mode, "rt": rt, "model_height": model_height, "sn": sn,
        "spectrum": [[184.0733, 999.0], [104.1075, 120.0]] if spectrum is None else spectrum,
    }


class TestDclDecode(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, results) -> str:
        path = Path(self.tmp.name) / "synthetic.dcl"
        path.write_bytes(dcl_bytes(results))
        return str(path)

    def test_header_and_seek_table_roundtrip(self):
        parsed = dcl_reader.deserialize_dcl(self.write([dcl_result(), dcl_result(precursor_mz=500.0)]))
        self.assertEqual([r["dcl_index"] for r in parsed], [0, 1])
        self.assertAlmostEqual(parsed[1]["precursor_mz"], 500.0)

    def test_a_file_without_the_dc_magic_is_rejected(self):
        path = Path(self.tmp.name) / "bad.dcl"
        path.write_bytes(b"XX" + b"\x00" * 30)
        with self.assertRaises(ValueError):
            dcl_reader.deserialize_dcl(str(path))

    def test_scan_block_fields_land_in_the_documented_offsets(self):
        parsed = dcl_reader.deserialize_dcl(self.write([dcl_result()]))[0]
        self.assertEqual(parsed["scan_id"], 1)
        self.assertEqual(parsed["raw_spec_id"], 17)
        self.assertEqual(parsed["ion_mode"], 1)
        self.assertAlmostEqual(parsed["rt"], 12.0)
        self.assertAlmostEqual(parsed["model_peak_height"], 9000.0)
        self.assertAlmostEqual(parsed["signal_to_noise"], 55.5, places=4)
        self.assertEqual(parsed["n_msms_peaks"], 2)
        self.assertEqual(parsed["msms_spectrum"], [[184.0733, 999.0], [104.1075, 120.0]])

    def test_spectrum_can_be_skipped(self):
        parsed = dcl_reader.deserialize_dcl(self.write([dcl_result()]), include_spectrum=False)[0]
        self.assertEqual(parsed["msms_spectrum"], [])
        self.assertEqual(parsed["n_msms_peaks"], 2)

    def test_top_n_keeps_the_strongest_peaks_but_restores_mz_order(self):
        spectrum = [[100.0, 10.0], [200.0, 900.0], [300.0, 500.0]]
        parsed = dcl_reader.deserialize_dcl(self.write([dcl_result(spectrum=spectrum)]), top_n_peaks=2)[0]
        self.assertEqual(parsed["msms_spectrum"], [[200.0, 900.0], [300.0, 500.0]])

    def test_lookup_matches_within_the_mass_tolerance(self):
        results = dcl_reader.deserialize_dcl(self.write([dcl_result(precursor_mz=760.5851)]))
        self.assertEqual(len(dcl_reader.get_msms_by_precursor(results, 760.5851)), 1)
        self.assertEqual(len(dcl_reader.get_msms_by_precursor(results, 760.60)), 0)
        self.assertEqual(len(dcl_reader.get_msms_by_precursor(results, 760.60, tol=0.05)), 1)

    def test_lookup_can_separate_co_eluting_peaks_by_retention_time(self):
        results = dcl_reader.deserialize_dcl(self.write([
            dcl_result(precursor_mz=760.5851, rt=3.0),
            dcl_result(precursor_mz=760.5851, rt=12.0),
        ]))
        hits = dcl_reader.get_msms_by_precursor(results, 760.5851, rt=12.0)
        self.assertEqual([h["dcl_index"] for h in hits], [1])

    def test_entries_without_peaks_are_not_returned_as_evidence(self):
        # 「MS/MS が取れていない」を「フラグメントが無い」と混同させないための線。
        results = dcl_reader.deserialize_dcl(self.write([dcl_result(spectrum=[])]))
        self.assertEqual(dcl_reader.get_msms_by_precursor(results, 760.5851), [])


# --------------------------------------------------------------------------
# 一度も実行されていなかった MCP ツールのスモーク
# --------------------------------------------------------------------------
class TestNeverExecutedToolsSmoke(unittest.TestCase):
    """登録・注釈・文書の検査からは名前として参照されるだけで、実際には
    一度も呼ばれていなかったツール群（実測カバレッジ 4〜10%）。

    これらは例外を握り潰してエラー文字列を返すため、「例外が出ない」だけでは
    何も確かめたことにならない。合成ファイルに入れた値が戻り値に現れることまで
    見て、成功経路を通ったことを確かめる。
    """

    def setUp(self):
        import tempfile
        from metabolomix.core import session_state as state

        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        # セッションはパーサ別スロットを持つ可変グローバル。試験間で持ち越さない。
        state.session = state.AnalysisSession()

    def write(self, name: str, payload: bytes) -> str:
        path = Path(self.tmp.name) / name
        path.write_bytes(payload)
        return str(path)

    def eic_file(self) -> str:
        spots = [
            eic_spot(0, rt=1.5, mz=100.25, samples=[
                {"file_id": 0, "top": 1.5, "left": 1.4, "right": 1.6,
                 "points": [(1.5, 50.0)]},
            ]),
            eic_spot(1, rt=12.5, mz=760.5, samples=[
                {"file_id": 1, "top": 12.5, "left": 12.4, "right": 12.6,
                 "points": [(12.5, 8000.0)]},
            ]),
        ]
        return self.write("synthetic.EIC.aef", css1_bytes(spots))

    def test_eic_parser_reports_the_summary(self):
        from metabolomix.eic.tools import eic_parser

        result = eic_parser(file_path=self.eic_file())
        self.assertIn("EIC解析完了", result)
        self.assertIn("total_spots", result)

    def test_eic_rank_by_max_intensity_puts_the_strong_spot_first(self):
        from metabolomix.eic.tools import eic_rank_by_max_intensity

        result = eic_rank_by_max_intensity(file_path=self.eic_file(), top_n=2)
        self.assertIn("EIC強度上位スポット", result)
        lines = [line for line in result.splitlines() if "spot" in line.lower()]
        self.assertTrue(lines, f"スポット行が出ていない: {result}")
        self.assertIn("760", lines[0])

    def test_eic_search_by_mz_range_filters(self):
        from metabolomix.eic.tools import eic_search_by_mz_range

        result = eic_search_by_mz_range(file_path=self.eic_file(), min_mz=700.0, max_mz=800.0)
        self.assertIn("一致件数: 1", result)

    def test_eic_search_by_rt_range_filters(self):
        from metabolomix.eic.tools import eic_search_by_rt_range

        result = eic_search_by_rt_range(file_path=self.eic_file(), min_rt=0.0, max_rt=5.0)
        self.assertIn("一致件数: 1", result)

    def test_arf2_annotate_identities_returns_a_table(self):
        from metabolomix.arf2.tools import arf2_annotate_identities

        path = self.write("synthetic.arf2", lz4_container(msgpack_stream([arf2_spot()])))
        result = arf2_annotate_identities(file_path=path, max_rows=10)
        self.assertIn("PC 34:1", result)

    def test_arf2_annotate_identities_shows_the_curation_flag(self):
        """Task 9: wrong フラグが付いたスポットは curation_flag 列に出る。"""
        from metabolomix.arf2.tools import arf2_annotate_identities
        from metabolomix.curation import flags as curation_flags

        path = Path(self.write(
            "synthetic.arf2", lz4_container(msgpack_stream([arf2_spot()]))))
        curation_flags.FlagStore(curation_flags.curation_dir(path)).append(
            [{"spot_id": 7, "flag": "wrong"}],
            alignment=curation_flags.alignment_key(path),
            review_id="r", source="user")

        result = arf2_annotate_identities(file_path=str(path), max_rows=10)

        header, *body = result.splitlines()[-2:]
        self.assertIn("curation_flag", header)
        self.assertIn("wrong", body[0])

    def test_arf2_annotate_identities_shows_assign(self):
        from metabolomix.arf2.tools import arf2_annotate_identities
        from metabolomix.curation import flags as curation_flags

        path = Path(self.write(
            "synthetic.arf2", lz4_container(msgpack_stream([arf2_spot()]))))
        curation_flags.FlagStore(curation_flags.curation_dir(path)).append(
            [{"spot_id": 7, "flag": "assign", "name": "PC 34:2", "level": "sum"}],
            alignment=curation_flags.alignment_key(path), review_id="cs-x", source="user")

        result = arf2_annotate_identities(file_path=str(path), max_rows=10)

        header, *body = result.splitlines()[-2:]
        column = header.split("\t").index("curation_flag")
        self.assertEqual(body[0].split("\t")[column], "assign:PC 34:2")

    def test_arf2_annotate_identities_survives_a_malformed_flags_file(self):
        """I6: 壊れた flags.jsonl でも一覧は返し、curation_flag を空にして 1 行で知らせる。"""
        from metabolomix.arf2.tools import arf2_annotate_identities
        from metabolomix.curation import flags as curation_flags

        path = Path(self.write(
            "synthetic.arf2", lz4_container(msgpack_stream([arf2_spot()]))))
        flags_path = curation_flags.FlagStore(curation_flags.curation_dir(path)).path
        flags_path.parent.mkdir(parents=True, exist_ok=True)
        flags_path.write_text('{"spot_id": 7, "fl\n', encoding="utf-8")

        result = arf2_annotate_identities(file_path=str(path), max_rows=10)

        self.assertIn("PC 34:1", result)
        notes = [l for l in result.splitlines() if l.startswith("#") and "flags.jsonl" in l]
        self.assertEqual(len(notes), 1)
        header, *body = result.splitlines()[-2:]
        self.assertEqual(body[0].split("\t")[header.split("\t").index("curation_flag")], "")

    def test_dcl_find_msms_returns_the_matching_spectrum(self):
        from metabolomix.dcl.tools import dcl_find_msms

        path = self.write("synthetic.dcl", dcl_bytes([dcl_result(precursor_mz=760.5851)]))
        result = dcl_find_msms(precursor_mz=760.5851, file_path=path)
        self.assertIn("184.0733", result)

    def test_dcl_find_msms_reports_not_found_rather_than_an_empty_hit(self):
        # not_found は「MS/MS が取得されていない」であって「期待フラグメントが無い」
        # ではない。この区別が MSI レベルの主張を左右する。
        from metabolomix.dcl.tools import dcl_find_msms

        path = self.write("synthetic.dcl", dcl_bytes([dcl_result(precursor_mz=760.5851)]))
        result = dcl_find_msms(precursor_mz=200.0, file_path=path)
        self.assertIn("not_found", result)


if __name__ == "__main__":
    unittest.main()
