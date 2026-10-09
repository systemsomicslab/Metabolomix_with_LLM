"""キュレーション用の合成入力（.arf2 / .dcl / .EIC.aef / .arf / .msp）。

実データはリポジトリに含めないので、各形式のレイアウトをここで組み立てる。
Key 番号の正準は docs/schema/*.md（MsScanMatchResult.md / AlignmentSpotProperty.md）。
"""
from __future__ import annotations

from pathlib import Path

import lz4.block
import msgpack

from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula


def at_keys(size: int, values: dict) -> list:
    row = [None] * size
    for key, value in values.items():
        row[key] = value
    return row


def chromxs(rt: float, mz: float) -> list:
    """ChromXs の入れ子表現（1=RT, 3=m/z）。"""
    return [[1, [rt]], [3, [mz]]]


def match_result(overrides: dict | None = None) -> list:
    """MsScanMatchResult の生配列（39 要素）。既定は「参照一致した MS/MS あり」。

    overrides: dict with numeric keys for Value at that Key position.
    """
    if overrides is None:
        overrides = {}
    values = {
        0: "PC 34:1", 1: "KEY-PC341", 2: 3.2,
        3: 0.81, 4: 0.64, 5: 0.9, 6: 12.0, 7: 0.7, 8: 0.0,
        9: 0.99, 10: 0.0, 11: 0.0, 12: -1.0, 13: 0.98,
        14: 5, 15: True, 16: True, 17: True, 18: False,
        19: True, 20: True, 21: False, 22: False, 23: False,
        24: -1, 26: 4, 27: "Msp20260101000000_lib_1", 28: -1,
        29: 0.0, 30: False, 31: 1, 32: 0.0, 33: True, 34: False,
        35: False, 36: 0.0, 37: float("nan"), 38: -1.0,
    }
    values.update(overrides)
    return at_keys(39, values)


def arf2_spot_raw(*, spot_id: int, name: str = "PC 34:1", mz: float = 760.5851,
                  rt: float = 12.0, ontology: str = "PC", adduct: str = "[M+H]+",
                  formula: str = "C42H82NO8P", ion_mode: int = 0,
                  representative_file_id: int = 0,
                  matches: list | None = None,
                  peak_links: list | None = None,
                  isotopes: list | None = None) -> list:
    """AlignmentSpotProperty の生配列（Key 0..59）。matches は match_result() の並び。
    peak_links は Key 10（IonFeatureCharacter）の PeakLinks で、(リンク先スポット ID, 種類の整数) の並び。
    isotopes は Key 53 IsotopicPeaks（要素は [相対強度 %, m/z, 単同位体からの差, コメント, 絶対強度]）。
    既定は M=100 / M+1=46 / M+2=12。"""
    if isotopes is None:
        isotopes = [[100.0, mz, 0.0, "", 10000.0], [46.0, mz + 1.00467, 1.00467, "", 4600.0],
                    [12.0, mz + 2.00934, 2.00934, "", 1200.0]]
    values = {
        0: spot_id, 1: spot_id, 3: representative_file_id,
        4: chromxs(rt, mz), 5: mz,
        10: [[-1.00782503207, 1, adduct, 1, 1 if ion_mode == 1 else 0, True, 0.0, 0.0, False, True],
             [0.0, 0, "", 0, 0, False, 0.0, 0.0, False, False], 1,
             [[int(i), int(c)] for i, c in (peak_links or [])], 0, -1, -1, bool(peak_links), -1],
        11: ion_mode, 12: name,
        13: [formula, 0.0], 14: ontology, 15: "", 16: "",
        31: 10000.0, 32: 100.0, 33: 20000.0, 34: 0.2,
        35: 30.0, 36: 50.0, 37: 10.0, 43: mz - 0.001, 44: mz + 0.001,
        49: 1.0, 51: 1.0, 53: isotopes, 54: [0.0, 1, adduct],
        56: [list(matches or []), {}, []],
        59: -1,
    }
    return at_keys(60, values)


def write_arf2(path: Path, spots: list[list]) -> Path:
    """MS-DIAL の外側コンテナ msgpack([header, msgpack(非圧縮長) + LZ4]) で書く。"""
    inner = b"".join(msgpack.packb(spot, use_bin_type=True) for spot in spots)
    payload = msgpack.packb(len(inner), use_bin_type=True) + lz4.block.compress(
        inner, store_size=False)
    path = Path(path)
    path.write_bytes(msgpack.packb(["hdr", payload], use_bin_type=True))
    return path


import struct
import textwrap

from tests.dcl_fixture import build_dcl_bytes


def css1_bytes(spots: list[dict]) -> bytes:
    """`.EIC.aef`（CSS1）。spots[i] = {"rt","mz","samples":[{"file_id","top","left","right","points"}]}"""
    chunks = []
    for spot in spots:
        chunk = bytearray(struct.pack("<ffffbi", spot["rt"], 0.0, spot["mz"], 0.0, 0,
                                      len(spot["samples"])))
        for sample in spot["samples"]:
            chunk.extend(struct.pack("<iifff", sample["file_id"], len(sample["points"]),
                                     sample["top"], sample["left"], sample["right"]))
            for x, y in sample["points"]:
                chunk.extend(struct.pack("<ff", x, y))
        chunks.append(bytes(chunk))
    header = 14 + 8 * len(chunks)
    offsets, cursor = [], header
    for chunk in chunks:
        offsets.append(cursor)
        cursor += len(chunk)
    body = bytearray(b"CSS1" + b"\x00" * 6 + struct.pack("<i", len(chunks)))
    for offset in offsets:
        body.extend(struct.pack("<q", offset))
    for chunk in chunks:
        body.extend(chunk)
    return bytes(body)


def arf_row(*, file_id: int, mz: float, rt: float, height: float, gap_filled: bool) -> list:
    """AlignmentChromPeakFeature の行（docs/schema/AlignmentChromPeakFeature.md）。"""
    row = [None] * 50
    row[0] = file_id
    row[1] = f"sample_{file_id}"
    row[2] = -2 if gap_filled else file_id + 100
    row[3] = file_id + 1000
    row[10] = {}
    row[15] = chromxs(rt, mz)
    row[18] = height
    row[20] = height * 2
    row[21] = height * 1.5
    row[22] = mz
    row[23] = 0
    row[24] = ""
    row[37] = [1.0, 20.0]
    return row


def write_arf(path: Path, groups: list[list[list]]) -> Path:
    inner = b"".join(msgpack.packb(item, use_bin_type=True) for item in [[0, 0, *groups]])
    payload = msgpack.packb(len(inner), use_bin_type=True) + lz4.block.compress(
        inner, store_size=False)
    path = Path(path)
    path.write_bytes(msgpack.packb(msgpack.ExtType(99, payload), use_bin_type=True))
    return path


def gaussian_points(center: float, *, height: float = 1000.0, sigma: float = 0.03,
                    n: int = 31, step: float = 0.01) -> list[list[float]]:
    import math
    start = center - (n // 2) * step
    return [[start + i * step, height * math.exp(-((start + i * step - center) ** 2) / (2 * sigma ** 2))]
            for i in range(n)]


LIBRARY_MSP = textwrap.dedent("""\
    NAME: PC 34:1
    PRECURSORMZ: 760.5851
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    RETENTIONTIME: 12.1
    INCHIKEY: KEY-PC341
    Num Peaks: 2
    184.07 999
    760.58 200

    NAME: PC 36:2
    PRECURSORMZ: 786.6007
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    RETENTIONTIME: 12.4
    INCHIKEY: KEY-PC362
    Num Peaks: 1
    184.07 999
""")


#: `write_alignment_set` の各スポットの match_result 上書き（`match_overrides` で足すときの土台）。
_KEYS_BY_SPOT = {0: {0: "PC 34:1", 14: 0, 27: "lib_1"},
                 1: {0: "PC 36:2", 1: "KEY-PC362", 14: 1, 27: "lib_1", 33: False, 34: True}}


def write_alignment_set(folder: Path, *, n_files: int = 3, eic_points: int = 31,
                        match_overrides: dict[int, dict] | None = None) -> dict:
    """スポット 3 件の一式を `AlignmentResult_x.*` として書く。

    spot 0: PC 34:1、参照一致・きれいなピーク（ok になるべき）
    spot 1: low score: PC 36:2、参照 RT から 1.5 分ずれ・ピークがギザギザ
    spot 2: Unknown（選択されない）

    spot 1 の match_result は Key1（InChIKey）を参照レコード（PC 36:2 = KEY-PC362）に
    合わせて上書きする。`match_result()` の既定値はどの候補も "KEY-PC341"（PC 34:1 の
    InChIKey）なので、上書きしないと `evidence._reference` の InChIKey 不一致ガード
    （別ライブラリのレコードを ScanID が指してしまった場合の安全弁）に spot 1 自身が
    引っかかり、参照が解決できなくなる（実装時に判明。ここは fixture 側の記述漏れの
    修正であり、ガードのロジックは変えていない）。

    `match_overrides` は {spot_id: {Key: 値}}。spot 0/1 の照合結果の Key を追加で上書きする
    （脂質規則フラグ Key 19/20/22 を試すため）。

    `eic_points`（既定 31）は各サンプルの EIC 点数。既定の `step=0.01` のまま増やすと、
    `evidence._trim` の窓（`(right-left) * EIC_WINDOW_FACTOR`、既定で ±0.3 分）を
    超えた点が捨てられるだけなので、payload の間引き（`EIC_MAX_POINTS`）を試すテストは
    61（オフセットがちょうど窓幅 0.3 分に収まる点数）のように窓に収まる値を渡すこと。
    """
    folder.mkdir(parents=True, exist_ok=True)
    specs = [
        {"spot_id": 0, "name": "PC 34:1", "mz": 760.5851, "rt": 12.1, "scan": 0,
         "matches": [match_result({0: "PC 34:1", 14: 0, 27: "lib_1"})]},
        {"spot_id": 1, "name": "low score: PC 36:2", "mz": 786.6007, "rt": 13.9, "scan": 1,
         "matches": [match_result({0: "PC 36:2", 1: "KEY-PC362", 14: 1, 27: "lib_1",
                                    33: False, 34: True})]},
        {"spot_id": 2, "name": "Unknown", "mz": 500.0, "rt": 5.0, "scan": None, "matches": []},
    ]
    for spot_id, keys in (match_overrides or {}).items():
        specs[spot_id]["matches"][0][:] = match_result({**_KEYS_BY_SPOT[spot_id], **keys})
    write_arf2(folder / "AlignmentResult_x.arf2", [
        arf2_spot_raw(spot_id=s["spot_id"], name=s["name"], mz=s["mz"], rt=s["rt"],
                      matches=s["matches"], representative_file_id=0)
        for s in specs])
    (folder / "AlignmentResult_x.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": s["mz"], "rt": s["rt"],
         "spectrum": [(184.07, 1000.0), (s["mz"], 150.0)] if s["scan"] is not None else []}
        for s in specs]))
    eic_spots, groups = [], []
    for s in specs:
        samples, rows = [], []
        for file_id in range(n_files):
            points = gaussian_points(s["rt"], n=eic_points)
            if s["spot_id"] == 1:
                for i in range(1, len(points) - 1, 2):
                    points[i][1] *= 0.3
            samples.append({"file_id": file_id, "top": s["rt"], "left": s["rt"] - 0.1,
                            "right": s["rt"] + 0.1, "points": points})
            rows.append(arf_row(file_id=file_id, mz=s["mz"], rt=s["rt"],
                                height=1000.0 * (file_id + 1), gap_filled=(file_id == n_files - 1)))
        eic_spots.append({"rt": s["rt"], "mz": s["mz"], "samples": samples})
        groups.append(rows)
    (folder / "AlignmentResult_x.EIC.aef").write_bytes(css1_bytes(eic_spots))
    write_arf(folder / "AlignmentResult_x_PeakProperties.arf", groups)
    (folder / "lib.msp").write_text(LIBRARY_MSP, encoding="utf-8")
    return {"arf2": folder / "AlignmentResult_x.arf2", "msp": folder / "lib.msp"}


SUGGEST_MSP = textwrap.dedent("""\
    NAME: PC 16:0_18:2
    PRECURSORMZ: {pc342:.4f}
    PRECURSORTYPE: [M+HCOO]-
    IONMODE: Negative
    INCHIKEY: KEY-PC342
    FORMULA: C42H80NO8P
    Num Peaks: 2
    255.23 999
    279.23 800

    NAME: PC 16:0_18:1
    PRECURSORMZ: {pc341:.4f}
    PRECURSORTYPE: [M+HCOO]-
    IONMODE: Negative
    INCHIKEY: KEY-PC341
    FORMULA: C42H82NO8P
    Num Peaks: 2
    255.23 999
    281.25 800

    NAME: PE 18:0_18:2
    PRECURSORMZ: {pe362:.4f}
    PRECURSORTYPE: [M-H]-
    IONMODE: Negative
    INCHIKEY: KEY-PE362
    FORMULA: C41H78NO8P
    Num Peaks: 2
    283.26 999
    279.23 700
""")


def _mz(formula: str, adduct: str) -> float:
    return mz_from_neutral(monoisotopic_mass(parse_formula(formula)), parse_adduct(adduct))


def write_suggest_set(folder: Path, *, with_param: bool = True, with_arf: bool = True) -> dict:
    """候補付け用の一式（負イオン、6 試料）。

    spot 0: PC 34:2 [M+HCOO]-、注釈付き（相手 Y）
    spot 1: PC 34:1 と注釈されているが、実は spot 0 の M+2（強度 0.1 倍・同時溶出）。wrong フラグの対象。
            MatchResults は代表（PC 16:0_18:1、ライブラリ 1 番）と、store に無い 2 位（PG 34:1、LibraryID 99）
    spot 2: 未注釈、MS/MS が PE 18:0_18:2 [M-H]- と合う（② で候補が出る）
    spot 3: 未注釈、spot 0 の [M+HCOO]- → [M-CH3]-（-60.0211）、spot 0 への FoundInUpperMsMs リンク
    spot 4: 未注釈、MS/MS なし、関係なし（候補は出ない）
    """
    folder.mkdir(parents=True, exist_ok=True)
    pc342 = _mz("C42H80NO8P", "[M+HCOO]-")
    pc341 = _mz("C42H82NO8P", "[M+HCOO]-")
    pe362 = _mz("C41H78NO8P", "[M-H]-")
    unknown = match_result({0: None, 1: "", 2: 0.0, 3: 0.0, 14: -1, 26: 1, 33: False})
    specs = [
        {"spot_id": 0, "name": "PC 34:2", "mz": pc342, "rt": 10.00, "adduct": "[M+HCOO]-",
         "formula": "C42H80NO8P", "ontology": "PC", "scale": 1.0, "spectrum": [(255.23, 999.0), (279.23, 800.0)],
         "matches": [match_result({0: "PC 16:0_18:2", 1: "KEY-PC342", 14: 0, 27: "lib_1"})], "links": []},
        {"spot_id": 1, "name": "PC 34:1", "mz": pc342 + 2.006710, "rt": 10.01, "adduct": "[M+HCOO]-",
         "formula": "C42H82NO8P", "ontology": "PC", "scale": 0.1, "spectrum": [(255.23, 999.0), (279.23, 700.0)],
         "matches": [match_result({0: "PC 16:0_18:1", 1: "KEY-PC341", 14: 1, 27: "lib_1", 31: 2}),
                     match_result({0: "PG 34:1", 1: "KEY-PG341", 14: 99, 27: "lib_1", 31: 1, 2: 2.0})],
         "links": [(0, 5)]},
        {"spot_id": 2, "name": "Unknown", "mz": pe362, "rt": 11.50, "adduct": "[M-H]-", "formula": "",
         "ontology": "", "scale": 0.5, "spectrum": [(283.26, 999.0), (279.23, 650.0)], "matches": [unknown],
         "links": []},
        {"spot_id": 3, "name": "Unknown", "mz": pc342 - 60.021129, "rt": 10.00, "adduct": "[M-H]-",
         "formula": "", "ontology": "", "scale": 0.3, "spectrum": [(255.23, 999.0)], "matches": [unknown],
         "links": [(0, 4)]},
        {"spot_id": 4, "name": "Unknown", "mz": 432.1234, "rt": 3.00, "adduct": "[M-H]-", "formula": "",
         "ontology": "", "scale": 0.7, "spectrum": [], "matches": [unknown], "links": []},
    ]
    write_arf2(folder / "AlignmentResult_x.arf2", [
        arf2_spot_raw(spot_id=s["spot_id"], name=s["name"], mz=s["mz"], rt=s["rt"], ontology=s["ontology"],
                      adduct=s["adduct"], formula=s["formula"], ion_mode=1, representative_file_id=0,
                      matches=s["matches"], peak_links=s["links"])
        for s in specs])
    (folder / "AlignmentResult_x.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": s["mz"], "rt": s["rt"], "spectrum": s["spectrum"]} for s in specs]))
    eic_spots, groups = [], []
    for s in specs:
        samples, rows = [], []
        for file_id in range(6):
            height = 1000.0 * (file_id + 1) * s["scale"]
            samples.append({"file_id": file_id, "top": s["rt"], "left": s["rt"] - 0.1, "right": s["rt"] + 0.1,
                            "points": gaussian_points(s["rt"], height=height)})
            rows.append(arf_row(file_id=file_id, mz=s["mz"], rt=s["rt"], height=height, gap_filled=False))
        eic_spots.append({"rt": s["rt"], "mz": s["mz"], "samples": samples})
        groups.append(rows)
    (folder / "AlignmentResult_x.EIC.aef").write_bytes(css1_bytes(eic_spots))
    if with_arf:
        write_arf(folder / "AlignmentResult_x_PeakProperties.arf", groups)
    if with_param:
        (folder / "Dataset_x_param_1.txt").write_text(
            "Ion mode: Negative\nSearched adduct ions: [M-H]-,[M+HCOO]-,[M+CH3COO]-\n"
            "MS1 tolerance for centroid: 0.01\nRetention time tolerance for alignment: 0.1\n", encoding="utf-8")
    (folder / "lib.msp").write_text(SUGGEST_MSP.format(pc342=pc342, pc341=pc341, pe362=pe362), encoding="utf-8")
    return {"arf2": folder / "AlignmentResult_x.arf2", "msp": folder / "lib.msp"}


def build_saved_review(folder: Path) -> dict:
    """`write_alignment_set` → `review.run_review` → `review.save_review` まで済ませる（送信系のテスト用）。
    `LIBRARY_CACHE_ENV` は呼び出し側が tmp に向けておく。"""
    from metabolomix.arf2.reader import load_catalog
    from metabolomix.curation import evidence, judge, review
    from metabolomix.library import store as library_store

    paths = write_alignment_set(folder)
    store = library_store.open_store(paths["msp"])
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        result = review.run_review(paths["arf2"], spots, store=store, ms2_tol=0.025,
                                   th=judge.resolve_thresholds(None), file_ids=None, max_traces=12,
                                   selection={"kind": "annotated"})
    finally:
        store.close()
    review.save_review(result)
    return {"paths": paths, "review": result}
