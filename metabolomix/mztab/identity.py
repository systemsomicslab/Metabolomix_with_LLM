"""SME同定情報からのInChIKey導出（純関数）。

優先順位: database_identifier 正規表現一致 → inchi → smiles (RDKit) → None。
RDKit 非インストール時は smiles/inchi 経路が silent fallback する。
spec §23 参照。

RDKit の解決結果はプロセス内で 1 度だけ確かめて memo する。失敗した import は
`sys.modules` に載らないため、memo が無いと**壊れた RDKit の機械では特徴量ごとに
DLL ロードを再試行する**（実測 60 ms/回。35,803 特徴量で理論 2,147 秒。実測でも
`dataset_load` が検体数と無関係に 15〜22 分かかっていた）。健全な機械では import が
キャッシュされるのでこの病理は見えない——任意依存が「壊れている」状態でだけ出る形。
memo を差し替えるテストは `_rdkit_cache` を `_UNSET` へ monkeypatch する。
"""
from __future__ import annotations
import re
from typing import Any

_INCHIKEY_RE = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')

# RDKit 解決結果の memo。_UNSET = 未判定、None = 使えない。
_UNSET: Any = object()
_rdkit_cache: Any = _UNSET


def _rdkit() -> tuple[Any, Any, Any] | None:
    """(Chem, MolToInchi, InchiToInchiKey) を返す。使えなければ None。

    ImportError だけでなく全例外を使用不可として扱う。Windows のアプリケーション
    制御がネイティブ DLL を塞ぐ環境では `import rdkit` は通るのに `rdkit.Chem` が
    落ちるなど、壊れ方が ImportError に収まらない。
    """
    global _rdkit_cache
    if _rdkit_cache is _UNSET:
        try:
            from rdkit import Chem
            from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
        except Exception:
            _rdkit_cache = None
        else:
            _rdkit_cache = (Chem, MolToInchi, InchiToInchiKey)
    return _rdkit_cache


def derive_inchikey(
    database_identifier: str | None,
    inchi: str | None,
    smiles: str | None,
) -> tuple[str | None, str]:
    """SME フィールドから InChIKey を導出する。

    戻り値: (inchikey, source)。source は
    "database_identifier" / "inchi_derived" / "smiles_derived" / "none"。
    """
    if database_identifier:
        candidate = database_identifier.strip()
        if _INCHIKEY_RE.match(candidate):
            return candidate, "database_identifier"

    if not inchi and not smiles:
        return None, "none"

    rdkit = _rdkit()
    if rdkit is None:
        return None, "none"
    Chem, MolToInchi, InchiToInchiKey = rdkit

    try:
        if inchi:
            ik = InchiToInchiKey(inchi.strip())
            if ik:
                return ik, "inchi_derived"
        if smiles:
            mol = Chem.MolFromSmiles(smiles.strip())
            if mol:
                inchi_str = MolToInchi(mol)
                if inchi_str:
                    ik = InchiToInchiKey(inchi_str)
                    if ik:
                        return ik, "smiles_derived"
    except Exception:
        # 文書化された silent fallback（spec §23）。診断を潰す点は
        # docs/task.md 2026-09-19(2) の HOLD として別途判断する。
        pass

    return None, "none"


def rdkit_available() -> bool:
    """RDKit で SMILES / InChI から InChIKey を導出できるか返す。"""
    return _rdkit() is not None
