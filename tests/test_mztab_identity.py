# tests/test_mztab_identity.py
import re
import sys
from unittest.mock import patch, MagicMock

import pytest

from metabolomix.mztab import identity
from metabolomix.mztab.identity import derive_inchikey, _INCHIKEY_RE


@pytest.fixture(autouse=True)
def _fresh_rdkit_memo(monkeypatch):
    """RDKit 解決の memo をテストごとに捨てる。

    `derive_inchikey` はプロセス内で 1 度しか RDKit を解決しないので、
    `sys.modules` を差し替えるテストは memo を無効化しないと前のテストの
    結論を引き継いでしまう。
    """
    monkeypatch.setattr(identity, "_rdkit_cache", identity._UNSET)

# PC 36:2 の本物の InChIKey
_VALID_IK = "IPCSVZSSVZVIGE-UHFFFAOYSA-N"
# PC 36:2 の SMILES（テスト用に簡略化してもよいが、実際の値を使う）
_PC362_SMILES = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"

def test_inchikey_regex_accepts_valid():
    assert _INCHIKEY_RE.match(_VALID_IK)

def test_inchikey_regex_rejects_short():
    assert not _INCHIKEY_RE.match("IPCSVZSSVZVIGE")

def test_derive_from_database_identifier():
    ik, src = derive_inchikey(_VALID_IK, None, None)
    assert ik == _VALID_IK
    assert src == "database_identifier"

def test_derive_ignores_non_inchikey_database_identifier():
    ik, src = derive_inchikey("HMDB:HMDB0000001", None, None)
    assert ik is None
    assert src == "none"

def test_derive_returns_none_when_all_absent():
    ik, src = derive_inchikey(None, None, None)
    assert ik is None
    assert src == "none"

def test_derive_from_smiles_with_rdkit(monkeypatch):
    # RDKit をモックして SMILES 経路を検証する
    mock_rdkit = MagicMock()
    mock_mol = MagicMock()
    mock_rdkit.Chem.MolFromSmiles.return_value = mock_mol
    mock_rdkit.Chem.inchi.MolToInchi.return_value = "InChI=1S/test"
    mock_rdkit.Chem.inchi.InchiToInchiKey.return_value = _VALID_IK
    monkeypatch.setitem(sys.modules, "rdkit", mock_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", mock_rdkit.Chem)
    monkeypatch.setitem(sys.modules, "rdkit.Chem.inchi", mock_rdkit.Chem.inchi)

    ik, src = derive_inchikey(None, None, _PC362_SMILES)
    assert ik == _VALID_IK
    assert src == "smiles_derived"

def test_derive_returns_none_when_rdkit_absent():
    with patch.dict(sys.modules, {"rdkit": None, "rdkit.Chem": None, "rdkit.Chem.inchi": None}):
        ik, src = derive_inchikey(None, None, _PC362_SMILES)
    assert ik is None
    assert src == "none"

def test_derive_returns_none_on_invalid_smiles(monkeypatch):
    mock_rdkit = MagicMock()
    mock_rdkit.Chem.MolFromSmiles.return_value = None  # 無効 SMILES → None mol
    monkeypatch.setitem(sys.modules, "rdkit", mock_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", mock_rdkit.Chem)
    monkeypatch.setitem(sys.modules, "rdkit.Chem.inchi", mock_rdkit.Chem.inchi)

    ik, src = derive_inchikey(None, None, "INVALID_SMILES_XYZ")
    assert ik is None
    assert src == "none"

def test_derive_from_inchi_with_rdkit(monkeypatch):
    mock_rdkit = MagicMock()
    mock_rdkit.Chem.inchi.InchiToInchiKey.return_value = _VALID_IK
    monkeypatch.setitem(sys.modules, "rdkit", mock_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", mock_rdkit.Chem)
    monkeypatch.setitem(sys.modules, "rdkit.Chem.inchi", mock_rdkit.Chem.inchi)

    ik, src = derive_inchikey(None, "InChI=1S/test", None)
    assert ik == _VALID_IK
    assert src == "inchi_derived"

def test_database_identifier_takes_priority_over_smiles():
    ik, src = derive_inchikey(_VALID_IK, None, _PC362_SMILES)
    assert src == "database_identifier"


def test_rdkit_is_resolved_only_once_when_it_is_unusable(monkeypatch):
    """壊れた RDKit の失敗 import は `sys.modules` に載らないので、memo が無いと
    特徴量ごとに DLL ロードを再試行する（実測 60 ms/回。35,803 特徴量で理論 2,147 秒）。
    """
    import builtins
    real_import = builtins.__import__
    attempts = []

    def broken_rdkit_import(name, *args, **kwargs):
        if name == "rdkit":
            attempts.append(name)
            raise ImportError("DLL load failed while importing rdchem")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", broken_rdkit_import)
    for _ in range(5):
        assert derive_inchikey(None, None, _PC362_SMILES) == (None, "none")

    assert len(attempts) == 1, f"rdkit の解決が {len(attempts)} 回走った"
