# metabolomix/pathway/msp_lipid.py
"""研究室 MSP の脂質行の ClassyFire Ontology と組成式から、MS-DIAL の脂質クラスと和の組成の名前を作る。

名前で読めない行（慣用名・IUPAC 名・ID）の経路（設計書 4.3）。推測はしない: Ontology が表に無い・
MS-DIAL で複数のクラスに分かれる、組成式が無い・そのクラスの形に合わない行は、作らずに理由を返す。

クラスごとに基準の種（`ref_n`:`ref_db` の組成式 `ref_formula`）を持つ。組成式の C・H から
    鎖の炭素 n  = C − (C₀ − n₀)
    二重結合 db = db₀ + ((H₀ + 2(n − n₀)) − H) / 2
を出す。C・H 以外の元素の数が基準と違う、db が整数でない・負、n < 1 の行は形に合わない。
鎖ごとの内訳は出さない（和の組成のみ）。PC と PE は C3H6 しか違わず組成式では見分けられないので、
名前に頭部の語（`cholin` / `ethanolamin` / `serin` / `inositol`）があればクラスと突き合わせる。

名前で読む経路（設計書 4.2）は `from_name`: 古い MS-DIAL の書式・英語のクラス名を MS-DIAL 5 の略記に直し、
GOSLIN で読めたものだけ、先頭の頭部から LbmClass を決める。
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from metabolomix.msdial.lipid_identity import _MSDIAL_QUALIFIER_RE, normalize_lipid_name

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
    """名前は `f"{prefix}{n}:{db}{suffix}"`。基準の種 `ref_n`:`ref_db` の組成式が `ref_formula`。

    `head_word` は名前に出うる頭部の語（無ければ None）。
    """
    lbm_class: str
    prefix: str
    suffix: str
    ref_n: int
    ref_db: int
    ref_formula: str
    head_word: str | None = None


CLASS_FORMS: dict[str, ClassForm] = {
    "PC": ClassForm("PC", "PC ", "", 34, 1, "C42H82NO8P", "cholin"),
    "PE": ClassForm("PE", "PE ", "", 34, 1, "C39H76NO8P", "ethanolamin"),
    "PS": ClassForm("PS", "PS ", "", 34, 1, "C40H76NO10P", "serin"),
    "PG": ClassForm("PG", "PG ", "", 34, 1, "C40H77O10P"),
    "PI": ClassForm("PI", "PI ", "", 34, 1, "C43H81O13P", "inositol"),
    "PA": ClassForm("PA", "PA ", "", 34, 1, "C37H71O8P"),
    "LPC": ClassForm("LPC", "LPC ", "", 16, 0, "C24H50NO7P", "cholin"),
    "LPE": ClassForm("LPE", "LPE ", "", 16, 0, "C21H44NO7P", "ethanolamin"),
    "PE P-": ClassForm("EtherPE", "PE P-", "", 32, 1, "C37H72NO7P", "ethanolamin"),
    "PE O-": ClassForm("EtherPE", "PE O-", "", 32, 1, "C37H74NO7P", "ethanolamin"),
    "PC P-": ClassForm("EtherPC", "PC P-", "", 34, 1, "C42H82NO7P", "cholin"),
    "PC O-": ClassForm("EtherPC", "PC O-", "", 34, 1, "C42H84NO7P", "cholin"),
    "TG": ClassForm("TG", "TG ", "", 54, 3, "C57H104O6"),
    "DG": ClassForm("DG", "DG ", "", 34, 1, "C37H70O5"),
    "MG": ClassForm("MG", "MG ", "", 18, 1, "C21H40O4"),
    "FA": ClassForm("FA", "FA ", "", 18, 1, "C18H34O2"),
    "CAR": ClassForm("CAR", "CAR ", "", 16, 0, "C23H45NO4"),
    "NAE": ClassForm("NAE", "NAE ", "", 16, 0, "C18H37NO2", "ethanolamin"),
    "SM": ClassForm("SM", "SM ", ";O2", 34, 1, "C39H79N2O6P", "cholin"),
}

_HEAD_WORDS = ("cholin", "ethanolamin", "serin", "inositol")

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


def class_name_conflict(key: str, name: str | None) -> bool:
    """名前（限定子なし）の頭部の語がクラスの鍵と食い違うか。名前に頭部の語が 1 つも無ければ False。"""
    text = (name or "").lower()
    words = {w for w in _HEAD_WORDS if w in text}
    return bool(words) and CLASS_FORMS[key].head_word not in words


def from_class_formula(ontology: str | None, formula: str | None,
                       name: str | None = None) -> tuple[str, str] | str:
    """ClassyFire の Ontology と組成式から `(名前（限定子なし）, LbmClass)`。`name` は元の名前（限定子なし）。

    作れなければ理由（この順に見る）: `ambiguous_ontology` / `unmapped_ontology` / `class_name_conflict`
    （名前の頭部の語がクラスと食い違う）/ `no_formula` / `formula_mismatch`。
    """
    onto = (ontology or "").strip()
    if onto in AMBIGUOUS_CLASSYFIRE:
        return "ambiguous_ontology"
    key = CLASSYFIRE_TO_CLASS.get(onto)
    if key is None:
        return "unmapped_ontology"
    if class_name_conflict(key, name):
        return "class_name_conflict"
    return from_class_formula_key(key, formula)


# 同位体の標識: 重水素 `d5` / `(d7)` / `-d31`（長鎖塩基の `d18:1` は `:` が続くので当たらない）と `13C`。
_ISOTOPE_LABEL = re.compile(r"(?<![A-Za-z0-9])d\d+(?![\d:])|(?<!\d)13C(?![a-z])")


def has_isotope_label(name: str | None) -> bool:
    """名前に同位体の標識（内部標準）があるか。"""
    return bool(_ISOTOPE_LABEL.search(name or ""))


# --- 名前から読む（経路 1） ---

AMBIGUOUS_CLASS = "ambiguous_class"

# 古い MS-DIAL の書式: `Cer[NS] d30:1` / `SM d34:2` / `SM d18:1/16:0`（d は ;O2、t は ;O3）。括弧は Cer / HexCer だけ。
_OLD_FORMAT = re.compile(r"^(Cer|HexCer|SM)(?:\[(NS|NDS)\])?\s+([dt])(\d+:\d+)(?:/(.+))?$")
_OLD_FORMAT_BASE = {"d": ";O2", "t": ";O3"}

_ENGLISH_HEADS = {
    "phosphatidylcholine": "PC",
    "phosphatidylethanolamine": "PE",
    "phosphatidylserine": "PS",
    "phosphatidylglyceride": "PG",
    "phosphatidylglycerol": "PG",
    "phosphatidylinositol": "PI",
    "phosphatidic acid": "PA",
}
_ENGLISH_NAME = re.compile(
    rf"^({'|'.join(_ENGLISH_HEADS)})(?:\s+(lyso|alkyl|alkenyl))?\s+(\d.*)$", re.IGNORECASE)
_ENGLISH_MODIFIER = {"alkyl": "O-", "alkenyl": "P-"}

# `;` の後ろが空白を含む注記（`LPC 16:0; PlaSMA ID-1486`）のときだけ切る。`;O2` の修飾は残す。
_ANNOTATION = re.compile(r";(?=[^;]*\s)")
# 酸化の修飾 `;O` / `;O2` …。Cer / HexCer / SM（長鎖塩基の ;O2 / ;O3）以外では経路 1 で読まない。
_OXIDIZED = re.compile(r";O(?![a-z])")
_SPHINGO_HEADS = ("Cer", "HexCer", "SM")

# 頭部 → LbmClass（同名のもの）。PC / PE の `O-` / `P-` は EtherPC / EtherPE。
_SAME_NAME_CLASSES = frozenset({
    "PC", "PE", "PS", "PG", "PI", "PA", "LPC", "LPE", "LPG", "LPI", "LPS", "LPA",
    "TG", "DG", "MG", "FA", "CAR", "CE", "NAE", "SM",
})
_ETHER_CLASSES = {"PC": "EtherPC", "PE": "EtherPE"}

_HEAD = re.compile(r"^([A-Za-z0-9]+)\s+([OP]-)?")
# 長鎖塩基 ;O2 ＋ 水酸基の注記の無い脂肪酸だけが Cer_NS / Cer_NDS に決まる（塩基の二重結合 ≥ 1 なら NS）。
_CERAMIDE_NS_NDS = re.compile(r"^(?:Cer|HexCer) \d+:(\d+);O2/\d+:\d+$")


def split_qualifier(name: str) -> tuple[str, str]:
    """MS-DIAL の限定子接頭辞を分ける。`("low score: ", "PC 34:1")`。無ければ `("", name.strip())`。"""
    m = _MSDIAL_QUALIFIER_RE.match(name)
    if m is None:
        return "", name.strip()
    return m.group(0), name[m.end():].strip()


def _to_msdial5(name: str) -> tuple[str, str | None]:
    """古い書式・英語のクラス名を MS-DIAL 5 の略記に直す。`(直した名前, 括弧のクラス NS/NDS または None)`。"""
    m = _OLD_FORMAT.match(name)
    if m:
        head, bracket, base, chain, fatty = m.groups()
        if bracket and head == "SM":
            return name, None
        return f"{head} {chain}{_OLD_FORMAT_BASE[base]}" + (f"/{fatty}" if fatty else ""), bracket
    m = _ENGLISH_NAME.match(name)
    if m:
        head, modifier, chains = m.groups()
        abbrev = _ENGLISH_HEADS[head.lower()]
        modifier = (modifier or "").lower()
        if modifier == "lyso":
            return f"L{abbrev} {chains.replace('-', '_')}", None
        return f"{abbrev} {_ENGLISH_MODIFIER.get(modifier, '')}{chains.replace('-', '_')}", None
    return name, None


def _ceramide_class(head: str, name: str, bracket: str | None) -> str:
    """Cer / HexCer の LbmClass。決められない形（和の組成だけ・;O3・水酸基付き）は `ambiguous_class`。"""
    if bracket:
        return f"{head}_{bracket}"
    m = _CERAMIDE_NS_NDS.match(name)
    if m is None:
        return AMBIGUOUS_CLASS
    return f"{head}_NS" if int(m.group(1)) >= 1 else f"{head}_NDS"


def from_name(name: str | None) -> tuple[str, str] | str | None:
    """限定子を外した名前から `(MS-DIAL 5 の名前, LbmClass)`。

    クラスが一意でなければ `"ambiguous_class"`、読めなければ None（経路 2 に回す）。
    返す名前は直した文字列であり、GOSLIN の正規化形ではない。Cer / HexCer のクラスの曖昧さは
    GOSLIN の可否より先に決める（水酸基付きの脂肪酸など、GOSLIN が読めない形でも曖昧と分かる）。
    Cer / HexCer / SM 以外の頭部で `;O…` の修飾がある名前（`FA 18:1;O (9-HODE)` も。注記を切る前に見る）も None。
    """
    if not name or not name.strip():
        return None
    text = name.split("|", 1)[0].strip()
    oxidized = _OXIDIZED.search(text) is not None
    m = _ANNOTATION.search(text)
    if m:
        text = text[:m.start()].strip()
    text, bracket = _to_msdial5(text)
    head_match = _HEAD.match(text)
    if head_match is None:
        return None
    head, ether = head_match.groups()
    if oxidized and head not in _SPHINGO_HEADS:
        return None
    if head in ("Cer", "HexCer"):
        if ether:
            return None
        lbm_class = _ceramide_class(head, text, bracket)
        if lbm_class == AMBIGUOUS_CLASS:
            return AMBIGUOUS_CLASS
    elif head in _SAME_NAME_CLASSES:
        if ether:
            lbm_class = _ETHER_CLASSES.get(head)
            if lbm_class is None:
                return None
        else:
            lbm_class = head
    else:
        return None
    if not normalize_lipid_name(text)["parse_ok"]:
        return None
    return text, lbm_class
