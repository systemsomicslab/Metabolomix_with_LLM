"""MS-DIAL の解析フォルダを tmp に作るテスト用 fixture（後続の Task が足していく）。"""
from __future__ import annotations

from pathlib import Path

import lz4.block
import msgpack


def pack_msdial(value) -> bytes:
    unpacked = msgpack.packb(value, use_bin_type=True)
    size_header = msgpack.packb(len(unpacked))
    compressed = lz4.block.compress(unpacked, store_size=False)
    extension = msgpack.ExtType(99, size_header + compressed)
    return msgpack.packb(extension, use_bin_type=True)


def analysis_file_row(
    file_id: int,
    name: str,
    class_id: str,
    *,
    file_type: int = 0,
    order: int | None = None,
    batch: int = 1,
) -> list:
    """AnalysisFileBean 1 行（Key0〜11。docs/schema/AnalysisFileBean.md）。"""
    return [
        f"C:/data/{name}.wiff",  # Key0
        name,  # Key1
        file_type,  # Key2
        class_id,  # Key3
        file_id + 1 if order is None else order,  # Key4
        file_id,  # Key5
        True,  # Key6
        "", [], "", "",  # Key7〜10（未使用）
        batch,  # Key11
    ]


def write_mddata(path, rows: list[list]) -> Path:
    path = Path(path)
    path.write_bytes(pack_msdial([rows, [], [], [], [], None, None, None, None]))
    return path


def write_msdial_folder(
    folder,
    *,
    spots: list[dict],
    samples: list[dict],
    with_mddata: bool = True,
    stem: str = "AlignmentResult_2026_02_01_00_00_00",
) -> dict:
    """`.arf` / `.arf2`（と `.mddata`）を 1 組、フォルダへ書く。

    spots[i]: name, mz, rt, ontology, adduct, heights, gap_filled, msms, matches, ion_mode,
    blank_filtered（`.arf2` Key58。すべて任意）。
    samples[j]: file_id, name, class_id（必須）と file_type / order / batch、
    `.arf` に書く名前 arf_name（既定は name。改名後の `.mddata` との突き合わせ用）。
    """
    from tests.curation_fixtures import arf2_spot_raw, arf_row, write_arf, write_arf2

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    groups, raws = [], []
    for spot_id, spot in enumerate(spots):
        mz = spot.get("mz", 500.0)
        rt = spot.get("rt", 12.0)
        heights = spot.get("heights") or [100.0] * len(samples)
        gap = spot.get("gap_filled") or [False] * len(samples)
        msms = spot.get("msms") or [False] * len(samples)
        rows = []
        for sample, height, gap_filled, has_msms in zip(samples, heights, gap, msms):
            row = arf_row(file_id=sample["file_id"], mz=mz, rt=rt, height=height,
                          gap_filled=gap_filled)
            row[1] = sample.get("arf_name", sample["name"])
            if has_msms:
                row[10] = {1: 30.0}
            rows.append(row)
        groups.append(rows)
        raw = arf2_spot_raw(
            spot_id=spot_id, name=spot.get("name", "Unknown"), mz=mz, rt=rt,
            ontology=spot.get("ontology", ""), adduct=spot.get("adduct", "[M+H]+"),
            ion_mode=spot.get("ion_mode", 0), matches=spot.get("matches"))
        if spot.get("blank_filtered"):
            raw[58] = True
        raws.append(raw)
    paths = {
        "arf": write_arf(folder / f"{stem}_PeakProperties.arf", groups),
        "arf2": write_arf2(folder / f"{stem}.arf2", raws),
    }
    if with_mddata:
        paths["mddata"] = write_mddata(folder / f"{stem}.mddata", [
            analysis_file_row(s["file_id"], s["name"], s["class_id"],
                              file_type=s.get("file_type", 0), order=s.get("order"),
                              batch=s.get("batch", 1))
            for s in samples])
    return paths


def _last_segment(name: str) -> str:
    return name.split("|")[-1]


def _has_name(name) -> bool:
    return name is not None and name not in ("Unknown", "null", "") and bool(str(name).strip())


#: MsScanMatchResult Key26 のビット（上流の SourceType）。
_SOURCE_TEXT_DB, _SOURCE_MANUAL = 16, 64


def _upstream_writes_sme(spot: dict) -> bool:
    """MS-DIAL の exporter が SME を書く条件（MztabFormatExport.cs の ShouldWriteSmeLine ∧
    SME_ID_REFS の条件）を、アダプタの `is_sme_like` を呼ばずに fixture 側で独立に書き下したもの。
    内部標準の正規化の条件（IsNormalizeIS / IsNormalizeSplash）は移植しない（既定は無効）。"""
    name = spot.get("name", "Unknown")
    matches = spot.get("matches") or []
    if not matches:
        return False
    rep_name, source = matches[0][0], int(matches[0][26] or 0)
    return bool(
        any(spot.get("msms") or [])                    # (a) MS/MS assigned
        and not source & _SOURCE_MANUAL                # (b) 手動修正でない
        and not source & _SOURCE_TEXT_DB               # (c) Text DB 由来でない
        and _has_name(name) and _has_name(rep_name)    # (d) 名前あり
        and not spot.get("blank_filtered")             # (e) ブランク除外でない
        and "no MS2" not in name)                      # (f)


def write_matching_mztab(path, folder_spec: dict) -> Path:
    """`write_msdial_folder` に渡したのと同じ spots / samples から、同じアラインメントの mzTab-M を書く。

    MS-DIAL 自身の書き出しに合わせる: SMF_ID = スポットの MasterAlignmentID（fixture では連番）、
    assay は FileID 昇順で同じ名前、値はピーク高、SME は MS/MS があり代表（matches の先頭）を
    持つスポットのうち上流の条件（`_upstream_writes_sme`）を満たすものだけ（名前は代表名の最後の `|` 区間）、SML は名前付きスポット全部
    （名前はスポット名の最後の区間）。
    """
    spots, samples = folder_spec["spots"], folder_spec["samples"]
    ordered = sorted(range(len(samples)), key=lambda j: samples[j]["file_id"])
    lines = [
        "MTD\tmzTab-version\t2.0.0-M",
        "MTD\tmzTab-mode\tComplete",
        "MTD\tmzTab-type\tQuantification",
        "MTD\tquantification_method\t[MS, MS:1001829, label-free raw feature quantitation, ]",
        "MTD\tsmall_molecule-quantification_unit\t[PRIDE, PRIDE:0000429, Abundance, ]",
    ]
    for k, j in enumerate(ordered, start=1):
        lines += [f"MTD\tms_run[{k}]-location\tfile:///C:/data/{samples[j]['name']}.wiff",
                  f"MTD\tassay[{k}]-ms_run_ref\tms_run[{k}]",
                  f"MTD\tassay[{k}]\t{samples[j]['name']}"]
    abundance = "\t".join(f"abundance_assay[{k}]" for k in range(1, len(samples) + 1))
    sml, smf, sme = [], [], []
    for spot_id, spot in enumerate(spots):
        name = spot.get("name", "Unknown")
        is_sme = _upstream_writes_sme(spot)
        refs = ""
        mz = spot.get("mz", 500.0)
        adduct = spot.get("adduct", "[M+H]+")
        if is_sme:
            refs = f"{spot_id}_1"
            sme.append(f"SME\t{refs}\tSMF:{spot_id}\tnull\tnull\t{_last_segment(spot["matches"][0][0])}\t"
                       f"null\tnull\t{adduct}\t{mz}\t1\t{mz}\t1\tMS-DIAL total score\t0.9")
        if _has_name(name):
            sml.append(f"SML\t{spot_id}\t{spot_id}\tnull\t{_last_segment(name)}\t{adduct}")
        heights = spot.get("heights") or [100.0] * len(samples)
        # 上流の Height は Math.Round(PeakHeightTop, 0)（C# 既定の偶数丸め）の整数。Python の round も偶数丸め。
        values = "\t".join(str(int(round(float(heights[j])))) for j in ordered)
        smf.append(f"SMF\t{spot_id}\t{spot_id if _has_name(name) else 'null'}\t{refs or 'null'}\t"
                   f"{mz}\t{spot.get('rt', 12.0) * 60.0}\t{values}")
    lines.append("SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_name\tadduct_ions")
    lines += sml
    lines.append("SFH\tSMF_ID\tSML_ID_REFS\tSME_ID_REFS\texp_mass_to_charge\t"
                 f"retention_time_in_seconds\t{abundance}")
    lines += smf
    lines.append("SEH\tSME_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_formula\tchemical_name\t"
                 "smiles\tinchi\tadduct_ion\texp_mass_to_charge\tcharge\ttheoretical_mass_to_charge\t"
                 "rank\tbest_id_confidence_measure\tbest_id_confidence_value")
    lines += sme
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
