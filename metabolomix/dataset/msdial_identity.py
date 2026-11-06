"""`.arf2` のスポットの同定を、mzTab-M の SME / SML と同じ 2 つの欄へ写す。

SME 相当（feature_candidates）の規則は上流 `MztabFormatExport.cs` の
`ShouldWriteSmeLine` の移植。上流の内部標準による正規化の条件は、パラメータを
読まないので移植しない。SML 相当（feature_annotations）は名前が Unknown / null / 空でない
全スポット（SME 相当のスポットを含む。mzTab-M と同じ）。
"""
from __future__ import annotations

from metabolomix.core.tool_helpers import _identity_tables
from metabolomix.msdial import lipid_identity

__all__ = ["is_sme_like", "sme_candidate", "sml_annotation", "identity_columns",
           "annotation_source", "IdentityCache"]

_SOURCE_MSP, _SOURCE_TEXT_DB, _SOURCE_GENERATED, _SOURCE_MANUAL = 4, 16, 32, 64


def _has_name(name) -> bool:
    return bool(name) and str(name).strip().lower() not in ("unknown", "null", "")


def is_sme_like(fields: dict, *, msms_assigned: bool) -> bool:
    """上流 ShouldWriteSmeLine の (a)〜(f)。"""
    rep = fields.get("representative")
    name = fields.get("name")
    return bool(
        msms_assigned                                              # (a)
        and rep is not None
        and not int(rep.get("source") or 0) & _SOURCE_MANUAL       # (b)
        and not fields.get("is_text_db_representative")            # (c)
        and name and str(name).strip()                             # (d)
        and not fields.get("blank_filtered")                       # (e)
        and "no MS2" not in name)                                  # (f)


def annotation_source(fields: dict) -> str:
    rep = fields.get("representative")
    if rep is None:
        return "none"
    source = int(rep.get("source") or 0)
    # 複数のビットが立つときの優先順: manual > text_db > generated_lipid > msp。
    for bit, label in ((_SOURCE_MANUAL, "manual"), (_SOURCE_TEXT_DB, "text_db"),
                       (_SOURCE_GENERATED, "generated_lipid"), (_SOURCE_MSP, "msp")):
        if source & bit:
            return label
    return "none"


def _inchikey(fields: dict) -> str | None:
    rep = fields.get("representative") or {}
    return fields.get("inchikey") or rep.get("inchikey") or None


def sme_candidate(fields: dict, fid: str) -> dict:
    """mzTab の `_all_candidates` の 1 件と同じキーの集合。mzTab にしか無い欄は None。"""
    rep = fields.get("representative") or {}
    inchikey = _inchikey(fields)
    return {
        "sme_id": f"arf2:{fid}", "rank": 1, "charge": None,
        "exp_mass_to_charge": fields.get("mz"), "theoretical_mass_to_charge": None,
        "confidence_value": rep.get("total_score"),
        "chemical_name": fields.get("name"), "database_identifier": None, "inchi": None,
        "smiles": fields.get("smiles"), "adduct": fields.get("adduct"),
        "identification_method": None, "confidence_measure": None, "spectra_ref": None,
        "inchikey": inchikey, "inchikey_source": "arf2" if inchikey else None,
    }


def sml_annotation(fields: dict, *, msms_assigned: bool | None = None) -> dict | None:
    """`_build_feature_annotations` の 1 件と同じキーの集合。曖昧さは持たない。"""
    if not _has_name(fields.get("name")):
        return None
    rep = fields.get("representative")
    inchikey = _inchikey(fields)
    return {
        "sml_id": f"arf2:{fields['master_id']}", "ambiguous": False,
        "name": fields["name"], "database_identifier": None,
        "chemical_formula": fields.get("formula"), "smiles": fields.get("smiles"),
        "adduct": fields.get("adduct"), "reliability": None, "confidence_measure": None,
        "ontology": fields.get("ontology"),
        "confidence_value": (rep or {}).get("total_score"),
        "ms2_presence": msms_assigned,
        "spectrum_matched": None if rep is None else bool(rep.get("is_spectrum_match")),
        "inchikey": inchikey, "inchikey_source": "arf2" if inchikey else None,
    }


def identity_columns(fields: dict, *, msms_assigned: bool, tables) -> dict:
    name = fields.get("name") or ""
    feature = {"name": "" if name.strip().lower() == "unknown" else name,
               "ontology": fields.get("ontology") or "", "has_msms": bool(msms_assigned)}
    block = lipid_identity.build_identity_block(
        feature, tables, mass_error_band="UNKNOWN", adduct_band="UNKNOWN")
    return {"normalized": block["goslin"]["normalized"],
            "refmet": block["reference"]["refmet_name"],
            "lipid_maps_category": block["reference"]["lipid_maps_category"],
            "msi_level": block["msi"]["level"]}


class IdentityCache:
    """同じ (name, ontology, msms_assigned) の標準化を 1 回だけ計算する。"""

    def __init__(self, tables=None):
        self.tables = tables if tables is not None else _identity_tables()
        self.computed = 0
        self._cache: dict[tuple, dict] = {}

    def columns(self, fields: dict, *, msms_assigned: bool) -> dict:
        key = (fields.get("name") or "", fields.get("ontology") or "", bool(msms_assigned))
        if key not in self._cache:
            self.computed += 1
            self._cache[key] = identity_columns(fields, msms_assigned=msms_assigned,
                                                tables=self.tables)
        return dict(self._cache[key])
