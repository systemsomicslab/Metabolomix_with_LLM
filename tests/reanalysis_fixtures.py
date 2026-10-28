"""再解析の道具（ハンドル・統計の spec・証拠の束）のテスト用の合成 mzTab-M。

実物の `parse_mztab` / `build_dataset_state` を通る小さなファイルを書く。研究室の
ライブラリや実データは使わない。
"""
from __future__ import annotations

from pathlib import Path

#: 既定の試料（群 A=HC 3 人・群 B=RA 3 人・QC 2 本）。
SAMPLES = ("HC1", "HC2", "HC3", "RA1", "RA2", "RA3", "QC1", "QC2")

#: 既定の特徴。values は SAMPLES の順。sme は SME 候補（rank 順）。
POS_FEATURES = (
    {"id": "1", "mz": 104.1070, "rt": 1.20, "values": [10, 12, 11, 20, 22, 21, 15, 16],
     "sme": [{"name": "Choline", "inchikey": "OEYIOHPDSNJKLS-UHFFFAOYSA-N", "formula": "C5H14NO",
              "adduct": "[M]1+", "theo": 104.1070}]},
    {"id": "2", "mz": 184.0733, "rt": 1.25, "values": [5, 6, 5, 5, 6, 5, 5, 5],
     "sme": []},
    {"id": "3", "mz": 500.2000, "rt": 5.00, "values": [100, 110, 90, 100, 105, 95, 100, 101],
     "sme": []},
)
NEG_FEATURES = (
    {"id": "1", "mz": 179.0561, "rt": 2.00, "values": [10, 11, 12, 13, 14, 15, 12, 12],
     "sme": [{"name": "Glucose", "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N", "formula": "C6H12O6",
              "adduct": "[M-H]1-", "theo": 179.0561}]},
    {"id": "2", "mz": 300.0000, "rt": 3.00, "values": [1, 2, 3, 4, 5, 6, 3, 3], "sme": []},
)


def write_mztab(path: Path, *, features=POS_FEATURES, samples=SAMPLES,
                sml_adduct: str = "[M+H]1+") -> Path:
    """mzTab-M 2.0.0-M を書く。SML の adduct_ions は極性の多数決に使われる。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "MTD\tmzTab-version\t2.0.0-M",
        "MTD\tmzTab-mode\tComplete",
        "MTD\tmzTab-type\tQuantification",
        "MTD\tquantification_method\t[MS, MS:1001829, label-free raw feature quantitation, ]",
        "MTD\tsmall_molecule-quantification_unit\t[PRIDE, PRIDE:0000429, Abundance, ]",
    ]
    for i, name in enumerate(samples, start=1):
        lines.append(f"MTD\tms_run[{i}]-location\tfile:///{name}.raw")
        lines.append(f"MTD\tassay[{i}]-ms_run_ref\tms_run[{i}]")
        lines.append(f"MTD\tassay[{i}]\t{name}")
    abundance = "\t".join(f"abundance_assay[{i}]" for i in range(1, len(samples) + 1))
    # 特徴に `matched`（True / False / None）があれば、SML に MS-DIAL の拡張列
    # `opt_global_ms2_presence` / `opt_global_spectrum_matched` を足し、SML に同定（rank 1 の SME の
    # 名前と InChIKey）を書く（同定の無い SML 行は注釈にならないため）。
    with_match = any("matched" in feature for feature in features)
    lines.append("SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_name\tadduct_ions"
                 + ("\topt_global_ms2_presence\topt_global_spectrum_matched" if with_match else ""))
    for feature in features:
        first = (feature.get("sme") or [None])[0] if with_match else None
        identity = f"{first['inchikey']}\t{first['name']}" if first else "null\tnull"
        match = ""
        if with_match:
            flag = feature.get("matched")
            match = ("\tnull\tnull" if flag is None
                     else f"\tTrue\t{'True' if flag else 'False'}")
        lines.append(f"SML\t{feature['id']}\t{feature['id']}\t{identity}\t{sml_adduct}{match}")
    lines.append("SFH\tSMF_ID\tSML_ID_REFS\tSME_ID_REFS\texp_mass_to_charge\t"
                 f"retention_time_in_seconds\t{abundance}")
    sme_rows = []
    for feature in features:
        refs = []
        for rank, sme in enumerate(feature.get("sme") or [], start=1):
            sme_id = f"{feature['id']}_{rank}"
            refs.append(sme_id)
            sme_rows.append((sme_id, feature, sme, rank))
        values = "\t".join(str(float(v)) if v is not None else "null" for v in feature["values"])
        lines.append(f"SMF\t{feature['id']}\t{feature['id']}\t{'|'.join(refs) or 'null'}\t"
                     f"{feature['mz']}\t{feature['rt'] * 60.0}\t{values}")
    lines.append("SEH\tSME_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_formula\tchemical_name\t"
                 "smiles\tinchi\tadduct_ion\texp_mass_to_charge\tcharge\ttheoretical_mass_to_charge\t"
                 "rank\tbest_id_confidence_measure\tbest_id_confidence_value")
    for sme_id, feature, sme, rank in sme_rows:
        charge = -1 if sme["adduct"].endswith("-") else 1
        lines.append(
            f"SME\t{sme_id}\tSMF:{feature['id']}\t{sme['inchikey']}\t{sme['formula']}\t{sme['name']}\t"
            f"null\tnull\t{sme['adduct']}\t{sme.get('exp', feature['mz'])}\t{charge}\t{sme['theo']}\t"
            f"{rank}\tMS-DIAL total score\t0.9")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
