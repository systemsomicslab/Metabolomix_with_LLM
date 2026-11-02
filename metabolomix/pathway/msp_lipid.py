# metabolomix/pathway/msp_lipid.py
"""研究室 MSP の脂質行の ClassyFire Ontology と組成式から、MS-DIAL の脂質クラスと和の組成の名前を作る。

名前で読めない行（慣用名・IUPAC 名・ID）の経路（設計書 4.3）。推測はしない: Ontology が表に無い・
MS-DIAL で複数のクラスに分かれる、組成式が無い・そのクラスの形に合わない行は、作らずに理由を返す。

クラスごとに基準の種（`ref_n`:`ref_db` の組成式 `ref_formula`）を持つ。組成式の C・H から
    鎖の炭素 n  = C − (C₀ − n₀)
    二重結合 db = db₀ + ((H₀ + 2(n − n₀)) − H) / 2
を出す。C・H 以外の元素の数が基準と違う、db が整数でない・負、n < 1 の行は形に合わない。
鎖ごとの内訳は出さない（和の組成のみ）。
"""
from __future__ import annotations

import re
from dataclasses import dataclass

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str | None) -> dict[str, int] | None:
    """`C57H104O6` → `{"C": 57, "H": 104, "O": 6}`。元素記号＋数だけの式でなければ None（`[2H]`・電荷・空）。"""
    if not formula or not formula.strip():
        return None
    text = formula.strip()
    counts: dict[str, int] = {}
    pos = 0
    while pos < len(text):
        m = _FORMULA_TOKEN.match(text, pos)
        if m is None:
            return None
        counts[m.group(1)] = counts.get(m.group(1), 0) + (int(m.group(2)) if m.group(2) else 1)
        pos = m.end()
    return counts


@dataclass(frozen=True)
class ClassForm:
    """名前は `f"{prefix}{n}:{db}{suffix}"`。基準の種 `ref_n`:`ref_db` の組成式が `ref_formula`。"""
    lbm_class: str
    prefix: str
    suffix: str
    ref_n: int
    ref_db: int
    ref_formula: str


CLASS_FORMS: dict[str, ClassForm] = {
    "PC": ClassForm("PC", "PC ", "", 34, 1, "C42H82NO8P"),
    "PE": ClassForm("PE", "PE ", "", 34, 1, "C39H76NO8P"),
    "PS": ClassForm("PS", "PS ", "", 34, 1, "C40H76NO10P"),
    "PG": ClassForm("PG", "PG ", "", 34, 1, "C40H75O10P"),
    "PI": ClassForm("PI", "PI ", "", 34, 1, "C43H81O13P"),
    "PA": ClassForm("PA", "PA ", "", 34, 1, "C37H71O8P"),
    "LPC": ClassForm("LPC", "LPC ", "", 16, 0, "C24H50NO7P"),
    "LPE": ClassForm("LPE", "LPE ", "", 16, 0, "C21H44NO7P"),
    "PE P-": ClassForm("EtherPE", "PE P-", "", 32, 1, "C37H72NO7P"),
    "PE O-": ClassForm("EtherPE", "PE O-", "", 32, 1, "C37H74NO7P"),
    "PC P-": ClassForm("EtherPC", "PC P-", "", 34, 1, "C42H82NO7P"),
    "PC O-": ClassForm("EtherPC", "PC O-", "", 34, 1, "C42H84NO7P"),
    "TG": ClassForm("TG", "TG ", "", 54, 3, "C57H104O6"),
    "DG": ClassForm("DG", "DG ", "", 34, 1, "C37H70O5"),
    "MG": ClassForm("MG", "MG ", "", 18, 1, "C21H40O4"),
    "FA": ClassForm("FA", "FA ", "", 18, 1, "C18H34O2"),
    "CAR": ClassForm("CAR", "CAR ", "", 16, 0, "C23H45NO4"),
    "NAE": ClassForm("NAE", "NAE ", "", 16, 0, "C18H37NO2"),
    "SM": ClassForm("SM", "SM ", ";O2", 34, 1, "C39H79N2O6P"),
}

# ClassyFire の名前 → CLASS_FORMS の鍵。MS-DIAL で 1 つのクラスに決まるものだけ。
CLASSYFIRE_TO_CLASS: dict[str, str] = {
    "Phosphatidylcholines": "PC",
    "Phosphatidylethanolamines": "PE",
    "Phosphatidylserines": "PS",
    "Phosphatidylglycerols": "PG",
    "Phosphatidylinositols": "PI",
    "1-phosphatidyl-1D-myo-inositols": "PI",
    "1,2-diacylglycerol-3-phosphates": "PA",
    "1-acyl-sn-glycero-3-phosphocholines": "LPC",
    "2-acyl-sn-glycero-3-phosphocholines": "LPC",
    "1-acyl-sn-glycero-3-phosphoethanolamines": "LPE",
    "1-(1Z-alkenyl),2-acylglycerophosphoethanolamines": "PE P-",
    "1-alkyl,2-acylglycerophosphoethanolamines": "PE O-",
    "1-(1Z-alkenyl),2-acyl-glycerophosphocholines": "PC P-",
    "1-alkyl,2-acylglycero-3-phosphocholines": "PC O-",
    "Triacylglycerols": "TG",
    "1,2-diacylglycerols": "DG",
    "1,3-diacylglycerols": "DG",
    "1-monoacylglycerols": "MG",
    "2-monoacylglycerols": "MG",
    "Long-chain fatty acids": "FA",
    "Very long-chain fatty acids": "FA",
    "Medium-chain fatty acids": "FA",
    "Lineolic acids and derivatives": "FA",
    "Acyl carnitines": "CAR",
    "N-acylethanolamines": "NAE",
    "Sphingomyelins": "SM",
}

# MS-DIAL で複数のクラスに分かれる（Cer_NS / Cer_NDS / HexCer_NS …）ので、この経路では決めない。
AMBIGUOUS_CLASSYFIRE: frozenset[str] = frozenset({
    "Ceramides", "Long-chain ceramides", "Very long-chain ceramides",
    "Phosphosphingolipids", "Glycosyl-N-acylsphingosines", "Cholesterols and derivatives",
})


def from_class_formula_key(key: str, formula: str | None) -> tuple[str, str] | str:
    """クラスの鍵と組成式から `(名前（限定子なし）, LbmClass)`。作れなければ理由 `no_formula` / `formula_mismatch`。"""
    if formula is None or not formula.strip():
        return "no_formula"
    counts = parse_formula(formula)
    form = CLASS_FORMS[key]
    ref = parse_formula(form.ref_formula)
    if counts is None or ref is None:
        return "formula_mismatch"
    if any(counts.get(el, 0) != ref.get(el, 0) for el in (counts.keys() | ref.keys()) - {"C", "H"}):
        return "formula_mismatch"
    n = counts.get("C", 0) - (ref.get("C", 0) - form.ref_n)
    if n < 1:
        return "formula_mismatch"
    twice_db_shift = (ref.get("H", 0) + 2 * (n - form.ref_n)) - counts.get("H", 0)
    if twice_db_shift % 2:
        return "formula_mismatch"
    db = form.ref_db + twice_db_shift // 2
    if db < 0:
        return "formula_mismatch"
    return f"{form.prefix}{n}:{db}{form.suffix}", form.lbm_class


def from_class_formula(ontology: str | None, formula: str | None) -> tuple[str, str] | str:
    """ClassyFire の Ontology と組成式から `(名前（限定子なし）, LbmClass)`。

    作れなければ理由（この順に見る）: `unmapped_ontology` / `ambiguous_ontology` / `no_formula` / `formula_mismatch`。
    """
    name = (ontology or "").strip()
    if name in AMBIGUOUS_CLASSYFIRE:
        return "ambiguous_ontology"
    key = CLASSYFIRE_TO_CLASS.get(name)
    if key is None:
        return "unmapped_ontology"
    return from_class_formula_key(key, formula)
