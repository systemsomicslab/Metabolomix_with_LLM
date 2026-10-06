"""脂質同定の標準化・信頼度（純ロジック層、MCP 非依存、オフライン）。

GOSLIN による脂質名正規化、同梱表による RefMet/LIPID MAPS ID 付与、MSI レベル推定。
ネットワークは一切使わない（pygoslin と同梱 TSV のみ）。
"""

from __future__ import annotations

import re
from pathlib import Path

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
        from pygoslin.parser.Parser import LipidParser
    except Exception as exc:  # pygoslin 未導入
        result["error"] = f"pygoslin unavailable: {exc}"
        return result
    try:
        lipid = LipidParser().parse(clean)
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

    ``msms_band``（``peak_verification.msms_evidence`` の band）を渡すと、Level 2 の
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
                "dcl_parser / dcl_find_msms で実フラグメントを確認すること。"
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
