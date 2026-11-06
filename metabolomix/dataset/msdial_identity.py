"""`.arf2` のスポットの同定を、mzTab-M の SME / SML と同じ 2 つの欄へ写す。

SME 相当（feature_candidates）の規則は上流 `MztabFormatExport.cs` の
`ShouldWriteSmeLine` の移植。上流の内部標準による正規化の条件は、パラメータを
読まないので移植しない。SML 相当（feature_annotations）は名前が Unknown / null / 空でない
全スポット（SME 相当のスポットを含む。mzTab-M と同じ）。
"""
from __future__ import annotations

from functools import lru_cache

from metabolomix.arf2.match_results import SOURCE_GENERATED_LIPID, SOURCE_MANUAL, SOURCE_MSP
from metabolomix.arf2.spot_fields import SOURCE_TEXT_DB
from metabolomix.msdial import lipid_identity

__all__ = ["is_sme_like", "sme_candidate", "sml_annotation", "identity_columns",
           "annotation_source", "resolve_inchikey", "spectrum_matched", "IdentityCache"]

@lru_cache(maxsize=1)
def _identity_tables():
    """同梱の RefMet / LIPID MAPS 対応表（core.tool_helpers を引かないよう、この層で持つ）。"""
    return lipid_identity.load_reference_tables()


def _has_name(name) -> bool:
    """上流の SML / SME の名前条件。大文字小文字を区別する（`spot.Name != "Unknown"` と同じ）。"""
    return name is not None and name not in ("Unknown", "null", "")


def _last_segment(name) -> str | None:
    """MS-DIAL の名前は `|` 区切りで、mzTab-M の chemical_name は最後の区間（上流 :350-351, :559）。"""
    return None if name is None else name.split("|")[-1]


def is_sme_like(fields: dict, *, msms_assigned: bool) -> bool:
    """上流 ShouldWriteSmeLine の (a)〜(f)。"""
    name = fields.get("name")
    return bool(
        msms_assigned                                              # (a)
        and not fields.get("is_manually_modified")                 # (b)
        and not fields.get("is_text_db_representative")            # (c)
        and _has_name(name) and str(name).strip()                  # (d) ＋ SME_ID_REFS の名前条件
        and not fields.get("blank_filtered")                       # (e)
        and "no MS2" not in name)                                  # (f)


def annotation_source(fields: dict) -> str:
    rep = fields.get("representative")
    if rep is None:
        return "none"
    source = int(rep.get("source") or 0)
    # 複数のビットが立つときの優先順: manual > text_db > generated_lipid > msp。
    for bit, label in ((SOURCE_MANUAL, "manual"), (SOURCE_TEXT_DB, "text_db"),
                       (SOURCE_GENERATED_LIPID, "generated_lipid"), (SOURCE_MSP, "msp")):
        if source & bit:
            return label
    return "other"  # 既知のビットが無い。"none" は代表が無いときだけ


def resolve_inchikey(fields: dict) -> str | None:
    rep = fields.get("representative") or {}
    return fields.get("inchikey") or rep.get("inchikey") or None


def spectrum_matched(fields: dict) -> bool | None:
    """代表の is_spectrum_match。代表が無ければ None（偽ではなく不明。mzTab の列欠落と同じ）。"""
    rep = fields.get("representative")
    return None if rep is None else bool(rep.get("is_spectrum_match"))


def sme_candidate(fields: dict, fid: str) -> dict:
    """mzTab の `_all_candidates` の 1 件と同じキーの集合。mzTab にしか無い欄は None。"""
    rep = fields.get("representative") or {}
    inchikey = resolve_inchikey(fields)
    return {
        "sme_id": f"arf2:{fid}", "rank": 1, "charge": None,
        "exp_mass_to_charge": fields.get("mz"), "theoretical_mass_to_charge": None,
        "confidence_value": rep.get("total_score"),
        "chemical_name": _last_segment(rep.get("name")) if rep else None,
        "database_identifier": None, "inchi": None,
        "smiles": fields.get("smiles"), "adduct": fields.get("adduct"),
        "identification_method": None, "confidence_measure": None, "spectra_ref": None,
        "inchikey": inchikey, "inchikey_source": "arf2" if inchikey else None,
    }


def sml_annotation(fields: dict, *, msms_assigned: bool | None = None) -> dict | None:
    """`_build_feature_annotations` の 1 件と同じキーの集合。曖昧さは持たない。"""
    if not _has_name(fields.get("name")):
        return None
    rep = fields.get("representative")
    inchikey = resolve_inchikey(fields)
    return {
        "sml_id": f"arf2:{fields['master_id']}", "ambiguous": False,
        "name": _last_segment(fields["name"]), "database_identifier": None,
        "chemical_formula": fields.get("formula"), "smiles": fields.get("smiles"),
        "adduct": fields.get("adduct"), "reliability": None, "confidence_measure": None,
        "ontology": fields.get("ontology"),
        "confidence_value": (rep or {}).get("total_score"),
        "ms2_presence": msms_assigned,
        "spectrum_matched": spectrum_matched(fields),
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
