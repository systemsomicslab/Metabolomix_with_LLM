"""脂質同定の標準化・信頼度（純ロジック層、MCP 非依存、オフライン）。

GOSLIN による脂質名正規化、同梱表による RefMet/LIPID MAPS ID 付与、MSI レベル推定。
ネットワークは一切使わない（pygoslin と同梱 TSV のみ）。
"""

from __future__ import annotations

import re
import threading
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple

from metabolomix.msdial import peak_verification as pv


# MS-DIAL は Name に信頼度の限定子接頭辞を付ける（実データ NEG で 317/856=37%）。
# 例: "no MS2: FA 5:0" / "low score: PC 35:2" / "w/o MS2:PC 34:1"。GOSLIN 解釈前に除く。
_MSDIAL_QUALIFIER_RE = re.compile(
    r'^\s*(?:no\s*ms2|w/?o\s*ms[12]|low\s*score|unsettled)\s*:\s*',
    re.IGNORECASE,
)


def _clean_msdial_name(name: str):
    """MS-DIAL の限定子接頭辞と候補区切り '|' を除き、GOSLIN が解釈できる単一の
    脂質ショートハンドに整える。戻り値 (clean_name, stripped: bool)。"""
    s = str(name).strip()
    stripped = False
    m = _MSDIAL_QUALIFIER_RE.match(s)
    if m:
        s = s[m.end():].strip()
        stripped = True
    if "|" in s:  # 複数候補（"Cer 24:1;O2|Cer 12:0;O2/12:1"）は先頭（species 表記）を採用
        s = s.split("|", 1)[0].strip()
        stripped = True
    return s, stripped


#: pygoslin のパーサは構築が重い（1 回 0.1 秒前後。文法 6 本）ので、モジュールで 1 つ遅延生成して
#: 使い回す。`LipidParser` は parse のたびに自分の状態を書き換えるので、ロックで排他する
#: （ハンドルの読み込みはワーカースレッドで走る。`curation/trend.py` と同じ流儀）。
_SHARED_PARSER = None
_PARSER_LOCK = threading.Lock()
#: 文法ごとの「決して読めない文字列」を見分ける表（`_grammar_screen`。`_SHARED_PARSER.parser_list` と同じ順）。
#: pygoslin の内部の形が想定と違えば None（表を使わず `LipidParser.parse` をそのまま呼ぶ）。
_SCREENS = None


def _grammar_screen(parser):
    """文法 1 本が受け付けうる文字列の必要条件の表: `(小文字にするか, 文字の集合, 隣り合う 2 文字の集合,
    先頭の文字の集合, 末尾の文字の集合)`。全文一致でない（EOF を使わない）文法は None（表で除かない）。

    pygoslin は文法をチョムスキー標準形（`TtoNT`: 文字 → それを導く規則、`NTtoNT`: 規則の対 → それを導く規則）
    にして CYK で読み、全文一致の文法は末尾に `EOF_SIGN` を足して読む。導ける文字列の隣り合う 2 文字は、どれも
    ある規則の対（B, C）について「B が導く文字列の末尾の文字 × C が導く文字列の先頭の文字」に入る（2 文字を
    両方含む最も低い節で分かれる）。各規則の先頭・末尾の文字の集合を不動点で求め、この集合を作る（実際に
    現れる 2 文字の上位集合なので、ここに無い 2 文字を含む文字列はその文法では決して読めない）。
    """
    from pygoslin.parser.Parser import FattyAcidParser

    if not getattr(parser, "used_eof", False):
        return None
    mask = (1 << 32) - 1
    first: dict = {}
    last: dict = {}
    for ch, rules in parser.TtoNT.items():
        for rule in rules:
            first.setdefault(rule, set()).add(ch)
            last.setdefault(rule, set()).add(ch)
    pairs = [(key >> 32, key & mask, rules) for key, rules in parser.NTtoNT.items()]
    changed = True
    while changed:
        changed = False
        for left, right, rules in pairs:
            head, tail = first.get(left), last.get(right)
            for rule in rules:
                if head and not head <= first.setdefault(rule, set()):
                    first[rule] |= head
                    changed = True
                if tail and not tail <= last.setdefault(rule, set()):
                    last[rule] |= tail
                    changed = True
    bigrams = {x + y for left, right, _ in pairs for x in last.get(left, ()) for y in first.get(right, ())}
    start = parser.START_RULE
    return (isinstance(parser, FattyAcidParser), frozenset(parser.TtoNT), frozenset(bigrams),
            frozenset(first.get(start, ())), frozenset(last.get(start, ())))


def _may_parse(screen, parser, text: str) -> bool:
    """その文法が `text` を読みうるか（False なら決して読めない。`parser.parse` と同じ前処理の後で表と照らす）。"""
    if screen is None:
        return True
    lower, alphabet, bigrams, heads, tails = screen
    t = (text.lower() if lower else text).strip(" ") + parser.EOF_SIGN
    if t[0] not in heads or t[-1] not in tails or not alphabet.issuperset(t):
        return False
    return all(t[k:k + 2] in bigrams for k in range(len(t) - 1))


def _parse_lipid(clean: str):
    """共有パーサで 1 つ解析する（pygoslin が無ければ ImportError、解析不能なら例外）。

    `LipidParser.parse` と同じく文法を順に試し、最初に読めた文法の結果を返す。ただし、その文法では決して読めない
    文字列（`_may_parse` が False）ではその文法を走らせない（結果は同じ。読めない名前は 6 本の文法をすべて CYK で
    走らせるので、天然物・医薬品の名前が多いと遅い）。
    """
    global _SHARED_PARSER, _SCREENS
    from pygoslin.parser.Parser import LipidParser

    with _PARSER_LOCK:
        if _SHARED_PARSER is None:
            _SHARED_PARSER = LipidParser()
            try:
                _SCREENS = [_grammar_screen(p) for p in _SHARED_PARSER.parser_list]
            except Exception:   # pygoslin の内部の形が違う: 表を使わない
                _SCREENS = None
        if _SCREENS is None:
            return _SHARED_PARSER.parse(clean)
        from pygoslin.domain.LipidExceptions import LipidException

        _SHARED_PARSER.parser_event_handler = None
        for screen, parser in zip(_SCREENS, _SHARED_PARSER.parser_list):
            if not _may_parse(screen, parser, clean):
                continue
            lipid = parser.parse(clean, raise_error=False)
            if lipid is not None:
                _SHARED_PARSER.parser_event_handler = parser.parser_event_handler
                return lipid
        raise LipidException("Lipid not found")


def normalize_lipid_name(name: str) -> dict:
    """脂質ショートハンド名を GOSLIN で正規化する（オフライン）。

    先に MS-DIAL の限定子接頭辞（'no MS2:'/'low score:' 等）と候補区切り '|' を除く。
    pygoslin が無い/解析不能でも例外を投げず parse_ok=False を返す。
    返り値: parse_ok / normalized / level / lipid_maps_category / stripped / error。
    """
    result = {"parse_ok": False, "normalized": None, "level": None,
              "lipid_maps_category": None, "stripped": False, "error": None}
    if not name or not str(name).strip():
        result["error"] = "empty name"
        return result
    clean, stripped = _clean_msdial_name(name)
    result["stripped"] = stripped
    if not clean:
        result["error"] = "empty after stripping qualifier"
        return result
    try:
        import pygoslin.parser.Parser  # noqa: F401
    except Exception as exc:  # pygoslin 未導入
        result["error"] = f"pygoslin unavailable: {exc}"
        return result
    try:
        lipid = _parse_lipid(clean)
        result["parse_ok"] = True
        result["normalized"] = lipid.get_lipid_string()
        # 構造レベル（species/molecular species/sn-position 等）。版差に強い経路を優先。
        level = None
        try:
            level = lipid.lipid.info.level.name
        except Exception:
            lvl = getattr(lipid, "level", None)
            level = lvl.name if hasattr(lvl, "name") else (str(lvl) if lvl else None)
        result["level"] = level
        try:
            result["lipid_maps_category"] = lipid.lipid.headgroup.lipid_category.name
        except Exception:
            result["lipid_maps_category"] = None
    except Exception as exc:
        result["error"] = f"parse error: {exc}"
    return result


@lru_cache(maxsize=4096)
def lipid_species_of(name) -> str | None:
    """脂質名を GOSLIN の種のレベル（``LipidLevel.SPECIES``）の文字列に直す（照合の鍵）。

    ``PC 18:3_16:1`` → ``PC 34:4``、``PE 16:0_18:1;O`` → ``PE 34:1;O``（修飾は残る）、
    ``PC O-34:1`` はエーテルのまま（``PC 34:1`` とは別）。``PC P-34:1`` は GOSLIN の規約で
    ``PC O-34:2`` になる。付加体の表記は落とす。限定子接頭辞と ``|`` は
    `normalize_lipid_name` と同じく除く。読めない・種より粗い（クラスだけ）名前は None。
    """
    if not isinstance(name, str) or not name.strip():
        return None
    clean, _ = _clean_msdial_name(name)
    if not clean:
        return None
    try:
        from pygoslin.domain.LipidLevel import LipidLevel
        import pygoslin.parser.Parser  # noqa: F401
    except Exception:  # pygoslin 未導入
        return None
    try:
        species = _parse_lipid(clean).lipid
        if species.info.level.value < LipidLevel.SPECIES.value:
            return None
        return species.get_lipid_string(LipidLevel.SPECIES)
    except Exception:
        return None


#: 同位体の標識（内部標準）: `(d7)` `[d9]` `(13C16)` `-d31`。脂質名かどうかを読むときだけ外す
#: （長鎖塩基の `d18:1` は `:` が続くので当たらない）。
_ISOTOPE_TAG_RE = re.compile(r"\s*[(\[](?:d\d+|13C\d*)[)\]]|-d\d+(?![\d:])")

#: 注釈の名前が無いことを表す MS-DIAL の値（`dataset.msdial_identity._has_name` と同じ。大小は問わない）。
_NO_NAME = ("", "unknown", "null")


def has_annotation_name(name) -> bool:
    """注釈の名前があるか（空・`Unknown`・`null` でない。限定子接頭辞だけの名前も無いとみなす）。"""
    if not isinstance(name, str):
        return False
    clean, _ = _clean_msdial_name(name)
    return clean.strip().casefold() not in _NO_NAME


@lru_cache(maxsize=65536)
def _parses(text: str) -> bool:
    """共有パーサで読めるか（pygoslin が無ければ ImportError を通す。例外は覚えない）。名前の形は名前どうしで
    重なる（限定子だけが違う名前など）ので、形ごとに覚える。"""
    try:
        _parse_lipid(text)
    except ImportError:
        raise
    except Exception:
        return False
    return True


#: 内部標準（SPLASH など）の名前の鎖の区切りの `-`（`SM d18:1-18:1` `PC 15:0-18:1` `TG 15:0-18:1-15:0`）。
#: 鎖と鎖の間（`数:数` の後で、次も `数:数`）の `-` だけを `/` に直す（`PC O-16:0` の `O-` は当たらない）。
_CHAIN_DASH_RE = re.compile(r"(\d+:\d+)-(?=[dt]?\d+:\d+)")
#: 脂質の略記の鎖（`34:1`）。名前から MS-DIAL の脂質クラスを決める（`msp_lipid.from_name`）のは、これを含む名前だけ。
_CHAIN_TOKEN_RE = re.compile(r"\d+:\d+")


@lru_cache(maxsize=65536)
def _name_forms(name: str) -> tuple[str, ...]:
    """名前を脂質名の文法（GOSLIN）で読むときに試す形（重ねない。この順）。

    限定子接頭辞と `|` を外した名前 → 注記・付加体の表記を切り、古い MS-DIAL の書式・英語のクラス名を略記に
    直した名前（`pathway.msp_lipid.shorthand`。`msms_rules.name_lipid_species` の `name_cleaned` と同じ直し）
    → それぞれから同位体の標識（`(d7)` など。内部標準）を外した名前 → その鎖の区切りの `-` を `/` に直した名前と、
    その略記（`SM d18:1-18:1(d9)` → `SM d18:1/18:1` → `SM 18:1;O2/18:1`）。
    """
    from metabolomix.pathway import msp_lipid  # msp_lipid がこのモジュールを import するので遅延

    forms: list[str] = []

    def add(text):
        if text and text not in forms:
            forms.append(text)

    clean, _ = _clean_msdial_name(name)
    add(clean)
    add(msp_lipid.shorthand(clean))
    for text in list(forms):
        bare = _ISOTOPE_TAG_RE.sub("", text).strip()
        add(bare)
        joined = _CHAIN_DASH_RE.sub(r"\1/", bare)
        if joined != bare:
            add(joined)
            add(msp_lipid.shorthand(joined))
    return tuple(forms)


@lru_cache(maxsize=65536)
def lipid_name_parses(name) -> bool | None:
    """注釈の名前が脂質名の文法（GOSLIN）で読めるか。脂質アッセイでクラスに数えてよいかの判定に使う。

    `_name_forms` の形のどれかが読めれば True、どれも読めなければ False（天然物・医薬品などの名前。
    `Ginsenoside Rb1` など）。GOSLIN は慣用名（`Oleic acid` `Cholesterol`）も読む。名前が無ければ False
    （呼び出し側が `has_annotation_name` で先に分ける）。pygoslin が無ければ None（判定できない）。
    """
    if not has_annotation_name(name):
        return False
    try:
        return any(_parses(text) for text in _name_forms(name))
    except ImportError:  # pygoslin 未導入
        return None


class LipidCall(NamedTuple):
    """脂質アッセイでの注釈 1 つの扱い（`lipid_annotation_class`）。

    `is_lipid`: 脂質クラスに数えるか（pygoslin が無ければ None＝判定できない）。
    `lipid_class`: 数えるクラス（脂質でなければ None）。元のクラスと違えば、そのクラスへ移して数える。
    `route`: `lbm_class` / `name` / `class_formula` / `name_grammar` / `unmatched` / `grammar_unavailable`。
    `rescued`: 規則 3（クラスと組成式）だけで脂質になったか（名前は脂質名の文法で読めない。研究室 MSP の IUPAC 名・ID）。
    """
    is_lipid: bool | None
    lipid_class: str | None
    route: str
    rescued: bool


#: `LipidCall.route` の値。
ROUTE_LBM_CLASS = "lbm_class"
ROUTE_NAME = "name"
ROUTE_CLASS_FORMULA = "class_formula"
ROUTE_NAME_GRAMMAR = "name_grammar"
ROUTE_UNMATCHED = "unmatched"
ROUTE_GRAMMAR_UNAVAILABLE = "grammar_unavailable"


def _text_or_none(value) -> str | None:
    if not isinstance(value, str):
        return None
    value = value.strip()
    return value if value and value.casefold() != "null" else None


def lipid_annotation_class(name, ontology, formula=None) -> LipidCall:
    """脂質アッセイで、この注釈（名前・クラス・組成式）を脂質クラスに数えてよいかと、数えるクラス。

    次のどれかなら脂質（上から順に見る）:
    1. クラス（ontology）が MS-DIAL の脂質クラス（`pathway.msdial_lipid_classes.MSDIAL_LIPID_CLASSES`）。クラスは
       そのまま（`lbm_class`）。MS-DIAL の略記はクラストークン＝Ontology の文法（output-format core §2.1）で、GOSLIN が
       読めない形（`Cer 12:0;O3/8:0(2OH)`・`AHexCer (O-14:1)12:1;O2/12:1`・`NAGly 17:0`）も含む。
    2. 名前を GOSLIN で読める（`lipid_name_parses`）。名前の頭部から MS-DIAL の脂質クラスが決まれば
       （`pathway.msp_lipid.from_name`）そのクラスへ移す（`name`）。
    3. クラスが研究室 MSP の ClassyFire の脂質クラス（`msp_lipid.CLASSYFIRE_TO_CLASS`）で、クラスと組成式から作った
       和の組成の略記（`msp_lipid.from_class_formula`）を GOSLIN で読める。その略記のクラス（`PC` `TG` `EtherPC` …）
       へ移す（`class_formula`）。名前が IUPAC 名・ID で読めない研究室 MSP の脂質を救う。名前で読めたが 2. で
       クラスが決まらない注釈も、ここでクラスが決まれば移す。
    2. で読めて 3. でクラスが決まらなければ、元のクラスのまま（`name_grammar`）。どれでもなければ脂質でない
    （`unmatched`。`Ginsenoside Rb1`・`Enniatin B`・ClassyFire のクラスと組成式が合わない行など）。
    注釈の名前・クラスそのものは変えない（数えるクラスを返すだけ）。pygoslin が無ければ判定しない（`is_lipid` が None）。
    """
    onto = ontology.strip() if isinstance(ontology, str) else ""
    return _lipid_call(name if isinstance(name, str) else None, onto, _text_or_none(formula))


@lru_cache(maxsize=65536)
def _lipid_call(name: str | None, onto: str, formula: str | None) -> LipidCall:
    from metabolomix.pathway import msp_lipid  # msp_lipid がこのモジュールを import するので遅延
    from metabolomix.pathway.msdial_lipid_classes import MSDIAL_LIPID_CLASSES

    if onto in MSDIAL_LIPID_CLASSES:
        return LipidCall(True, onto, ROUTE_LBM_CLASS, False)
    unavailable = LipidCall(None, None, ROUTE_GRAMMAR_UNAVAILABLE, False)
    parses = lipid_name_parses(name)
    if parses is None:
        return unavailable
    named = has_annotation_name(name)
    if parses and _CHAIN_TOKEN_RE.search(name):
        for text in _name_forms(name):
            by_name = msp_lipid.from_name(text)
            if isinstance(by_name, tuple):
                return LipidCall(True, by_name[1], ROUTE_NAME, False)
    rest = msp_lipid.split_qualifier(name)[1] if named else None
    built = msp_lipid.from_class_formula(onto, formula, rest)
    if isinstance(built, tuple):
        try:
            if _parses(built[0]):
                return LipidCall(True, built[1], ROUTE_CLASS_FORMULA, not parses)
        except ImportError:  # pygoslin 未導入
            return unavailable
    if parses:
        return LipidCall(True, onto or None, ROUTE_NAME_GRAMMAR, False)
    return LipidCall(False, None, ROUTE_UNMATCHED, False)


def may_move_lipid_class(name, ontology) -> bool:
    """`lipid_annotation_class` がこの注釈を元のクラスと違うクラスへ移しうるか（GOSLIN を使わない安い前段）。

    移すのは、クラスが MS-DIAL の脂質クラスでなく、研究室 MSP の ClassyFire の対応表（`msp_lipid.CLASSYFIRE_TO_CLASS`）に
    あるか、名前に鎖（`34:1`）があるときだけ。False なら、脂質と読めてもクラスは元のまま（図のクラスの展開が、移った
    特徴を探す範囲を絞るのに使う）。
    """
    from metabolomix.pathway import msp_lipid  # msp_lipid がこのモジュールを import するので遅延
    from metabolomix.pathway.msdial_lipid_classes import MSDIAL_LIPID_CLASSES

    onto = ontology.strip() if isinstance(ontology, str) else ""
    if onto in MSDIAL_LIPID_CLASSES:
        return False
    return onto in msp_lipid.CLASSYFIRE_TO_CLASS or bool(_CHAIN_TOKEN_RE.search(name if isinstance(name, str) else ""))


def is_lipid_annotation(name, ontology, formula=None) -> bool | None:
    """脂質アッセイで、この注釈を脂質クラスに数えてよいか（`lipid_annotation_class` の `is_lipid`）。

    名前の無い注釈は、クラスと組成式（規則 1・3）だけで決まる（呼び出し側が `has_annotation_name` で先に分ける）。
    pygoslin が無ければ None（判定できない）。
    """
    return lipid_annotation_class(name, ontology, formula).is_lipid


def _read_tsv(path: Path) -> list[dict]:
    """タブ区切り表を dict のリストで読む。存在しなければ空リスト。"""
    rows = []
    if not path.is_file():
        return rows
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            if len(values) == len(header):
                rows.append(dict(zip(header, values)))
    return rows


#: 同梱表の置き場所（リポジトリルート直下の reference/）。起動元の作業
#: ディレクトリに依存させない——相対パスで探すと表が見つからず `_read_tsv` が
#: 空を返し、全クラスが「同梱表に無い」扱いになって ID 付与が黙って全滅する。
REFERENCE_DIR = Path(__file__).resolve().parents[2] / "reference"


def load_reference_tables(reference_dir=None) -> dict:
    """同梱の RefMet / LIPID MAPS クラス対応表を読み込む（オフライン）。

    `reference_dir` 省略時は `REFERENCE_DIR`（作業ディレクトリに依存しない）。
    """
    base = Path(reference_dir) if reference_dir is not None else REFERENCE_DIR
    lm = {r["class_token"]: r for r in _read_tsv(base / "lipidmaps_classes.tsv")}
    rm = {r["class_token"]: r for r in _read_tsv(base / "refmet_map.tsv")}
    return {"lipidmaps": lm, "refmet": rm}


def map_to_reference(class_token, tables) -> dict:
    """クラストークンを同梱表で RefMet 名・LIPID MAPS カテゴリ/メインクラスに写像する。

    表に無ければ matched=False＋caveat（推測はしない）。
    """
    token = (class_token or "").strip().lower()
    lm = tables.get("lipidmaps", {}).get(token)
    rm = tables.get("refmet", {}).get(token)
    if not lm and not rm:
        return {"matched": False, "lipid_maps_category": None,
                "lipid_maps_main_class": None, "refmet_name": None,
                "caveat": f"クラス '{class_token}' は同梱マッピング表に無いため ID 未付与。"}
    return {
        "matched": True,
        "lipid_maps_category": lm["lipid_maps_category"] if lm else None,
        "lipid_maps_main_class": lm["lipid_maps_main_class"] if lm else None,
        "refmet_name": rm["refmet_name"] if rm else None,
        "caveat": None,
    }


def msi_level(*, name, ontology, has_msms, mass_error_band, adduct_band,
              msms_band=None) -> dict:
    """決定論的シグナルから MSI 同定信頼度レベルを推定する（ヒューリスティック）。

    Level 1（標準品照合）は主張しない。
    - 2: 名称あり かつ MS/MS取得 かつ 精密質量整合（putative annotated compound）
    - 3: クラス（ontology）は分かるが上記を満たさない（putative class）
    - 4: 名称もクラスも無い（unknown）

    ``msms_band``（MS/MS の証拠の band）を渡すと、Level 2 の
    根拠が**実スペクトル**（``PASS``）なのか**取得フラグのみ**（``FLAG_ONLY``）なのかを
    rationale で開示する。PAI2 の has_msms は取得参照の有無を示すだけなので、
    フラグ由来の Level 2 は弱い根拠であることを下流に伝える必要がある。
    """
    has_name = bool(name and str(name).strip())
    has_class = bool(ontology and str(ontology).strip() and str(ontology).strip() != "Unknown")
    mass_ok = mass_error_band == "PASS"
    adduct_ok = adduct_band in ("PASS", "UNKNOWN")  # FAIL は同定を疑う

    if has_name and has_msms and mass_ok and adduct_ok:
        if msms_band == "FLAG_ONLY":
            rationale = (
                "名称あり＋精密質量整合。ただし MS/MS は取得フラグのみで実スペクトルを"
                "確認できていない（.dcl 未添付）ため、この Level 2 の根拠は弱い。"
                "dcl_find_msms で実フラグメントを確認すること。"
            )
        elif msms_band == "PASS":
            rationale = (
                "名称あり＋実 MS/MS スペクトル確認＋精密質量整合。"
                "標準品照合ではないため Level 2。"
            )
        else:
            rationale = "名称あり＋MS/MS取得＋精密質量整合。標準品照合ではないため Level 2。"
        return {"level": 2, "label": "putative annotated compound",
                "rationale": rationale, "heuristic": True}
    if has_class:
        return {"level": 3, "label": "putative class-level",
                "rationale": "クラス（ontology）は判別できるが MS/MS または質量整合が不十分。",
                "heuristic": True}
    return {"level": 4, "label": "unknown",
            "rationale": "名称・クラスとも無し。", "heuristic": True}


def build_identity_block(feature, tables, *, mass_error_band, adduct_band,
                         msms_band=None):
    """1 feature の同定標準化ブロック（GOSLIN + reference + MSI）を組み立てる。

    feature は pai2/arf2 のキー name / ontology / has_msms を用いる。``msms_band`` を
    渡せば、実スペクトル（PASS）と取得フラグのみ（FLAG_ONLY）を区別して MSI の
    rationale に反映する（ARF2 のように .dcl を引けない経路では None のまま）。
    """
    name = feature.get("name") or ""
    ontology = feature.get("ontology") or ""
    has_msms = bool(feature.get("has_msms")) or msms_band == "PASS"
    class_token = pv.extract_class_token(name, ontology)
    goslin = normalize_lipid_name(name) if name.strip() else {
        "parse_ok": False, "normalized": None, "level": None,
        "lipid_maps_category": None, "error": "no name"}
    reference = map_to_reference(class_token, tables)
    msi = msi_level(name=name, ontology=ontology, has_msms=has_msms,
                    mass_error_band=mass_error_band, adduct_band=adduct_band,
                    msms_band=msms_band)
    return {"class_token": class_token, "goslin": goslin,
            "reference": reference, "msi": msi}
