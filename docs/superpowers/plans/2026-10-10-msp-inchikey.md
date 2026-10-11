# MSP の INCHIKEY で mzTab-M の同定に InChIKey を付ける Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Console に渡した研究室 MSP の `INCHIKEY` を `MspDB:<名前>` で引き、RDKit が使えない機械でも mzTab-M の同定に InChIKey を付ける。

**Architecture:** 新モジュール `metabolomix/mztab/msp_inchikey.py` が MSP の `NAME`/`INCHIKEY` 行だけを読んで名前→InChIKey 表を作る（stdlib のみ、プロセス内 memo）。`build_dataset_state` は表を任意引数で受け取り、「InChIKey 形の database_identifier → MSP 表 → RDKit」の順で解決する。`load_dataset_state(job_path=...)` はジョブの実効メソッドの `Msp file path` から表を作って渡す。

**Tech Stack:** Python 3.14（`C:/Python314/python.exe`）、pytest。

**Spec:** `docs/superpowers/specs/2026-10-10-msp-inchikey-design.md`（承認 2026-10-10）

## Global Constraints

- `metabolomix/mztab/loading.py` は `metabolomix.core.session_state` / `metabolomix.core.mcp_core` / `metabolomix.tools` を import しない（`tests/test_dataset_loading_service.py::test_loading_does_not_import_session_or_mcp_core`）。新モジュール `msp_inchikey.py` も stdlib と `metabolomix.mztab.identity` だけを import する（`library.store` は使わない）。
- 差次的エクスポートの列は変えない。`metabolomix/analysis/export_contract.py` の `CONTRACT_VERSION` は上げない。
- `inchikey_source` の新しい値は `msp_library`（この綴りで固定）。
- `inchikey_coverage["msp_lookup"]` のキーは `status` / `file` / `resolved` / `conflicting_names` / `unmatched`。`status` は `used` / `no_msp` / `unavailable` / `not_applicable`。
- MSP の置き場所（パス）を戻り値・警告・文書・テストに書かない。出すのはファイル名だけ（研究室ライブラリの非開示）。テストの MSP は `tmp_path` に作る。
- `derive_inchikey`（`metabolomix/mztab/identity.py`）のシグネチャと挙動は変えない。
- テストは `C:/Python314/python.exe -m pytest tests -q` をリポジトリ（この worktree）ルートから。
- コミットで全テストが走る（`.githooks/pre-commit`、数分）。自動化した `git commit` はバックグラウンドで実行し、出力をファイルへ落とす。実行中は作業ツリーを編集しない。
- 文書・コメントは日本語、既存の文体（だ・である調）に合わせる。

## Review Focus

1. MSP が CRLF 改行・先頭 BOM 付きでも名前と InChIKey が正しく読める（Windows で作られた MSP が普通）。→ Task 1 のテスト `test_crlf_and_bom_are_tolerated`。
2. 名前にコロンや空白を含む（例 `MspDB:LPC 18:0`）`database_identifier` を、先頭の `MspDB:` だけ外して引く。→ Task 1 `test_lookup_keeps_colons_in_the_name`。
3. INCHIKEY を持たないエントリの直後のエントリへ、前の名前が漏れない（空行で区切りを確定する）。→ Task 1 `test_entry_without_inchikey_does_not_borrow_the_next_key`。
4. 実効メソッドの `Msp file path` が相対パスでもメソッドファイル基準で解決する。→ Task 3 `test_relative_msp_path_is_resolved_against_the_method_file`。
5. 陽・陰の MSP を交互に読んでも memo が互いを追い出さない（1 GB 級の再読込を避ける）。→ Task 1 `test_tables_for_two_files_are_memoized_independently`。

---

### Task 1: MSP の名前→InChIKey 表（`msp_inchikey.py`）

**Files:**
- Create: `metabolomix/mztab/msp_inchikey.py`
- Test: `tests/test_msp_inchikey.py`

**Interfaces:**
- Consumes: `metabolomix.mztab.identity._INCHIKEY_RE`（`re.Pattern`、`^[A-Z]{14}-[A-Z]{10}-[A-Z]$`）
- Produces:
  - `MSP_DB_PREFIX: str = "MspDB:"`
  - `@dataclass(frozen=True) class MspInchikeyTable: file_name: str; by_name: dict[str, str]; conflicting: frozenset[str]`
  - `load_msp_inchikey_table(path: str | Path) -> MspInchikeyTable`（`OSError` はそのまま投げる）
  - `lookup(table: MspInchikeyTable, database_identifier: str | None) -> tuple[str | None, str]`。第 2 要素は `"resolved"` / `"conflicting"` / `"unmatched"` / `"not_msp"`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_msp_inchikey.py`:

```python
"""研究室 MSP の NAME → INCHIKEY 表（spec 2026-10-10 msp-inchikey-design §3.1）。

研究室ライブラリには触れない。MSP は tmp_path に数エントリだけ書く。
"""
from __future__ import annotations

import os

import pytest

from metabolomix.mztab import msp_inchikey
from metabolomix.mztab.msp_inchikey import (
    MSP_DB_PREFIX, MspInchikeyTable, load_msp_inchikey_table, lookup)

CIT = "RHGKLRLOHDJJDR-BYPYZUCNSA-N"   # L-citrulline
ORN = "AHLPHDHHMVZTML-BYPYZUCNSA-N"   # L-ornithine
SUC = "KDYFGRWQOYBRFD-UHFFFAOYSA-N"   # succinic acid
OTHER = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"


def _entry(name, inchikey=None):
    lines = [f"NAME: {name}", "PRECURSORMZ: 176.1030", "PRECURSORTYPE: [M+H]+"]
    if inchikey is not None:
        lines.append(f"INCHIKEY: {inchikey}")
    lines += ["Num Peaks: 1", "70.065\t1000", ""]
    return "\n".join(lines) + "\n"


def _write(tmp_path, text, name="lib.msp", newline="\n"):
    p = tmp_path / name
    p.write_bytes(text.replace("\n", newline).encode("utf-8"))
    return p


@pytest.fixture(autouse=True)
def _fresh_memo(monkeypatch):
    monkeypatch.setattr(msp_inchikey, "_memo", {})


def test_names_map_to_their_inchikey(tmp_path):
    p = _write(tmp_path, _entry("Citrulline", CIT) + _entry("Ornithine", ORN))
    table = load_msp_inchikey_table(p)
    assert table.by_name == {"Citrulline": CIT, "Ornithine": ORN}
    assert table.conflicting == frozenset()
    assert table.file_name == "lib.msp"


def test_empty_and_malformed_inchikeys_are_skipped(tmp_path):
    p = _write(tmp_path, _entry("Empty", "") + _entry("Short", "ABC-DEF")
               + _entry("Lower", CIT.lower()) + _entry("Citrulline", CIT))
    assert load_msp_inchikey_table(p).by_name == {"Citrulline": CIT}


def test_same_name_same_key_is_kept_once(tmp_path):
    p = _write(tmp_path, _entry("Citrulline", CIT) + _entry("Citrulline", CIT))
    table = load_msp_inchikey_table(p)
    assert table.by_name == {"Citrulline": CIT}
    assert table.conflicting == frozenset()


def test_same_name_different_keys_is_conflicting(tmp_path):
    p = _write(tmp_path, _entry("Leucine", CIT) + _entry("Leucine", OTHER)
               + _entry("Leucine", CIT))
    table = load_msp_inchikey_table(p)
    assert "Leucine" not in table.by_name
    assert table.conflicting == frozenset({"Leucine"})


def test_field_names_are_case_insensitive(tmp_path):
    p = _write(tmp_path, f"Name: Citrulline\nInChIKey: {CIT}\nNum Peaks: 0\n\n")
    assert load_msp_inchikey_table(p).by_name == {"Citrulline": CIT}


def test_crlf_and_bom_are_tolerated(tmp_path):
    text = "\ufeff" + _entry("Citrulline", CIT) + _entry("Ornithine", ORN)
    p = _write(tmp_path, text, newline="\r\n")
    assert load_msp_inchikey_table(p).by_name == {"Citrulline": CIT, "Ornithine": ORN}


def test_entry_without_inchikey_does_not_borrow_the_next_key(tmp_path):
    # 2 件目は NAME を持たない壊れたエントリ。1 件目の名前を引き継がせない。
    text = _entry("NoKey") + f"PRECURSORMZ: 1\nINCHIKEY: {CIT}\n\n"
    assert load_msp_inchikey_table(_write(tmp_path, text)).by_name == {}


def test_lookup_resolves_only_mspdb_identifiers(tmp_path):
    table = MspInchikeyTable("lib.msp", {"Citrulline": CIT}, frozenset({"Leucine"}))
    assert lookup(table, "MspDB:Citrulline") == (CIT, "resolved")
    assert lookup(table, "MspDB:Leucine") == (None, "conflicting")
    assert lookup(table, "MspDB:Unknown") == (None, "unmatched")
    assert lookup(table, "TextDB:Citrulline") == (None, "not_msp")
    assert lookup(table, "Citrulline") == (None, "not_msp")
    assert lookup(table, None) == (None, "not_msp")
    assert lookup(table, "") == (None, "not_msp")
    assert MSP_DB_PREFIX == "MspDB:"


def test_lookup_keeps_colons_in_the_name(tmp_path):
    table = MspInchikeyTable("lib.msp", {"LPC 18:0": CIT, "Citrulline": ORN}, frozenset())
    assert lookup(table, "MspDB:LPC 18:0") == (CIT, "resolved")
    assert lookup(table, "  MspDB:Citrulline  ") == (ORN, "resolved")


def test_table_is_memoized_until_the_file_changes(tmp_path, monkeypatch):
    p = _write(tmp_path, _entry("Citrulline", CIT))
    calls = []
    real_scan = msp_inchikey._scan
    monkeypatch.setattr(msp_inchikey, "_scan",
                        lambda path: calls.append(path) or real_scan(path))
    first = load_msp_inchikey_table(p)
    assert load_msp_inchikey_table(p) is first
    assert len(calls) == 1
    p.write_text(_entry("Ornithine", ORN), encoding="utf-8")
    st = p.stat()
    os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    assert load_msp_inchikey_table(p).by_name == {"Ornithine": ORN}
    assert len(calls) == 2


def test_tables_for_two_files_are_memoized_independently(tmp_path, monkeypatch):
    pos = _write(tmp_path, _entry("Citrulline", CIT), name="pos.msp")
    neg = _write(tmp_path, _entry("Succinate", SUC), name="neg.msp")
    calls = []
    real_scan = msp_inchikey._scan
    monkeypatch.setattr(msp_inchikey, "_scan",
                        lambda path: calls.append(path) or real_scan(path))
    for _ in range(2):
        assert load_msp_inchikey_table(pos).by_name == {"Citrulline": CIT}
        assert load_msp_inchikey_table(neg).by_name == {"Succinate": SUC}
    assert len(calls) == 2


def test_missing_file_raises_oserror(tmp_path):
    with pytest.raises(OSError):
        load_msp_inchikey_table(tmp_path / "absent.msp")


def test_module_imports_no_mcp_or_library_store():
    import ast
    from pathlib import Path
    tree = ast.parse(Path(msp_inchikey.__file__).read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert not any(m.startswith(("metabolomix.core.mcp_core", "metabolomix.core.session_state",
                                 "metabolomix.library", "metabolomix.tools"))
                   for m in imported), imported
```

- [ ] **Step 2: テストが落ちることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_msp_inchikey.py -q`
Expected: FAIL（`ModuleNotFoundError: No module named 'metabolomix.mztab.msp_inchikey'`）

- [ ] **Step 3: 実装を書く**

`metabolomix/mztab/msp_inchikey.py`:

```python
"""研究室 MSP の NAME → INCHIKEY 表（spec 2026-10-10 msp-inchikey-design §3.1）。

`library_mode="msp_only"` の Console 実行では、mzTab-M の同定は
`database_identifier = MspDB:<MSP の NAME>` で出て、InChIKey の列は無い。RDKit が
使えない機械（Windows のアプリケーション制御がネイティブ DLL を拒否する）でも
InChIKey に届くよう、Console に渡した MSP から名前で引く。

読むのは各エントリの `NAME:` と `INCHIKEY:` 行だけで、ピーク行は解さない（研究室の
MSP は 1 GB 級・数百万エントリで、全フィールドを解すと遅い）。同じ名前に別の
InChIKey があるものは、名前では決められないので表から外して `conflicting` に入れる。

stdlib と `metabolomix.mztab.identity` だけを import する。`mztab/loading.py`
（pipeline の worker も通る）から呼ぶので、`mcp_core` を引き込む `library.store` の
SQLite キャッシュは使わない。代わりに同じプロセスでは (パス, サイズ, 更新時刻) を
鍵に表を使い回す。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from metabolomix.mztab.identity import _INCHIKEY_RE

__all__ = ["MSP_DB_PREFIX", "MspInchikeyTable", "load_msp_inchikey_table", "lookup"]

#: MS-DIAL が mzTab-M の `database_identifier` に書く MSP 由来の接頭辞。
MSP_DB_PREFIX = "MspDB:"

_BOM = b"\xef\xbb\xbf"


@dataclass(frozen=True)
class MspInchikeyTable:
    """名前 → InChIKey。`file_name` は来歴用のファイル名だけ（置き場所は持たない）。"""
    file_name: str
    by_name: dict[str, str]
    conflicting: frozenset[str]


#: (絶対パス, サイズ, 更新時刻 ns) → 表。陽・陰の 2 本を交互に読んでも追い出さない。
_memo: dict[tuple[str, int, int], MspInchikeyTable] = {}


def load_msp_inchikey_table(path: str | Path) -> MspInchikeyTable:
    """MSP を読み、名前 → InChIKey 表を返す。読めなければ OSError をそのまま投げる。"""
    p = Path(path).resolve()
    st = p.stat()
    key = (str(p), st.st_size, st.st_mtime_ns)
    hit = _memo.get(key)
    if hit is None:
        hit = _scan(p)
        _memo[key] = hit
    return hit


def _scan(path: Path) -> MspInchikeyTable:
    by_name: dict[str, str] = {}
    conflicting: set[str] = set()
    name: str | None = None
    with path.open("rb") as f:
        for i, raw in enumerate(f):
            if i == 0 and raw.startswith(_BOM):
                raw = raw[len(_BOM):]
            line = raw.strip()
            if not line:
                # 空行がエントリの区切り。名前を次のエントリへ持ち越さない。
                name = None
                continue
            head = line[:9].upper()
            if head.startswith(b"NAME:"):
                name = line[5:].decode("utf-8", "replace").strip() or None
            elif head == b"INCHIKEY:" and name is not None:
                inchikey = line[9:].decode("ascii", "replace").strip()
                if not _INCHIKEY_RE.match(inchikey) or name in conflicting:
                    continue
                known = by_name.get(name)
                if known is None:
                    by_name[name] = inchikey
                elif known != inchikey:
                    del by_name[name]
                    conflicting.add(name)
    return MspInchikeyTable(path.name, by_name, frozenset(conflicting))


def lookup(table: MspInchikeyTable, database_identifier: str | None) -> tuple[str | None, str]:
    """`MspDB:<名前>` を表で引く。

    戻り値の第 2 要素: `resolved`（付いた）/ `conflicting`（同名に別キーがあり決めない）/
    `unmatched`（MSP 由来なのに名前が表に無い）/ `not_msp`（MSP 由来の同定ではない）。
    名前は先頭の接頭辞だけを外して引く（`LPC 18:0` のようにコロンを含みうる）。
    """
    text = (database_identifier or "").strip()
    if not text.startswith(MSP_DB_PREFIX):
        return None, "not_msp"
    name = text[len(MSP_DB_PREFIX):].strip()
    inchikey = table.by_name.get(name)
    if inchikey is not None:
        return inchikey, "resolved"
    if name in table.conflicting:
        return None, "conflicting"
    return None, "unmatched"
```

- [ ] **Step 4: テストが通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_msp_inchikey.py -q`
Expected: PASS（全件）

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/mztab/msp_inchikey.py tests/test_msp_inchikey.py
git commit -m "feat(mztab): 研究室 MSP の NAME → INCHIKEY 表を足す" > "$TEMP/msp_task1_commit.log" 2>&1; echo "exit=$?" >> "$TEMP/msp_task1_commit.log"
```

コミットメッセージの末尾に `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` を付ける。

---

### Task 2: `build_dataset_state` で MSP 表を使う

**Files:**
- Modify: `metabolomix/mztab/dataset_state.py`（`build_dataset_state` と `_build_feature_annotations`、import 行）
- Test: `tests/test_dataset_state.py`（末尾に追記）

**Interfaces:**
- Consumes: Task 1 の `MspInchikeyTable`、`lookup`
- Produces:
  - `build_dataset_state(parse_result: dict, filename: str, source_path: str | Path, *, msp_table: MspInchikeyTable | None = None) -> DatasetState`
  - `ds.inchikey_coverage["by_source"]` に `msp_library` キー（初期値 0）
  - `ds.inchikey_coverage["msp_lookup"] = {"status": "used" | "not_applicable", "file": str | None, "resolved": int, "conflicting_names": int, "unmatched": int}`。数は**照合の試行回数**（特徴ごとの SME 最良証拠 1 件 + SML 行 × 参照特徴）。`status` は表があれば `used`、無ければ `not_applicable`（Task 3 の loader が `no_msp` / `unavailable` に上書きする）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_dataset_state.py` の末尾に追記:

```python
# ---------- MSP の INCHIKEY（spec 2026-10-10 msp-inchikey-design）----------

from metabolomix.mztab.msp_inchikey import MspInchikeyTable

_CIT = "RHGKLRLOHDJJDR-BYPYZUCNSA-N"
_ORN = "AHLPHDHHMVZTML-BYPYZUCNSA-N"
_DBIK = "IPCSVZSSVZVIGE-UHFFFAOYSA-N"

# msp_only の Console が書く実形状: database_identifier は MspDB:<名前>、smiles あり、
# inchi は null。特徴 1 は SME（MS/MS 証拠）、特徴 2 は SML だけ（MS1 注釈）、
# 特徴 3 は名前が衝突、特徴 4 は MSP に無い名前、特徴 5 は InChIKey 形の識別子。
_MSPDB_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tmzTab-mode\tComplete
    MTD\tmzTab-type\tQuantification
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_formula\tsmiles\tinchi\tchemical_name\tadduct_ions\treliability\tbest_id_confidence_measure\tbest_id_confidence_value
    SML\t1\t1\tMspDB:Citrulline\tnull\tNC(=O)NCCCC(N)C(O)=O\tnull\tCitrulline\t[M+H]1+\t530\tnull\tnull
    SML\t2\t2\tMspDB:Ornithine\tnull\tNCCCC(N)C(O)=O\tnull\tlow score: Ornithine\t[M+H]1+\t530\tnull\tnull
    SML\t3\t3\tMspDB:Leucine\tnull\tCC(C)CC(N)C(O)=O\tnull\tLeucine\t[M+H]1+\t530\tnull\tnull
    SML\t4\t4\tMspDB:Unlisted\tnull\tCCO\tnull\tUnlisted\t[M+H]1+\t530\tnull\tnull
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\t1\t176.1030\t300.0\t1000.0
    SMF\t2\tnull\t133.0972\t310.0\t900.0
    SMF\t3\tnull\t132.1019\t320.0\t800.0
    SMF\t4\tnull\t47.0491\t330.0\t700.0
    SMF\t5\t2\t786.6\t340.0\t600.0
    SEH\tSME_ID\tdatabase_identifier\tsmiles\tinchi\tchemical_name\trank
    SME\t1\tMspDB:Citrulline\tNC(=O)NCCCC(N)C(O)=O\tnull\tCitrulline\t1
    SME\t2\tIPCSVZSSVZVIGE-UHFFFAOYSA-N\tnull\tnull\tPC 36:2\t1
""")

_TABLE = MspInchikeyTable("lab.msp", {"Citrulline": _CIT, "Ornithine": _ORN},
                          frozenset({"Leucine"}))


def _build_msp(tmp_path, msp_table=_TABLE):
    p = tmp_path / "Height_mspdb.mzTab"
    p.write_text(_MSPDB_MZTAB, encoding="utf-8")
    return build_dataset_state(parse_mztab(p), p.name, p, msp_table=msp_table)


def test_mspdb_evidence_gets_the_msp_inchikey(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    assert ds.feature_metadata["1"]["inchikey"] == _CIT
    assert ds.feature_metadata["1"]["inchikey_source"] == "msp_library"


def test_mspdb_sml_annotation_gets_the_msp_inchikey(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    assert ds.feature_annotations["2"]["inchikey"] == _ORN
    assert ds.feature_annotations["2"]["inchikey_source"] == "msp_library"


def test_conflicting_and_unmatched_names_get_no_inchikey(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    assert ds.feature_annotations["3"]["inchikey"] is None
    assert ds.feature_annotations["4"]["inchikey"] is None


def test_inchikey_shaped_identifier_wins_over_the_msp(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    assert ds.feature_metadata["5"]["inchikey"] == _DBIK
    assert ds.feature_metadata["5"]["inchikey_source"] == "database_identifier"


def test_msp_wins_over_rdkit(tmp_path, monkeypatch):
    """RDKit が使えても MSP の値を採る（機械によって結果を変えない）。"""
    def rdkit_derive(database_identifier, inchi, smiles):
        if smiles:
            return "ZZZZZZZZZZZZZZ-ZZZZZZZZZZ-Z", "smiles_derived"
        return None, "none"
    monkeypatch.setattr("metabolomix.mztab.dataset_state.derive_inchikey", rdkit_derive)
    ds = _build_msp(tmp_path)
    assert ds.feature_metadata["1"]["inchikey"] == _CIT
    assert ds.feature_metadata["1"]["inchikey_source"] == "msp_library"
    # MSP で決まらない名前は従来どおり RDKit へ落ちる。
    assert ds.feature_annotations["4"]["inchikey_source"] == "smiles_derived"


def test_by_source_counts_every_feature_once(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    cov = ds.inchikey_coverage
    assert sum(cov["by_source"].values()) == cov["total_features"] == 5
    assert cov["by_source"]["msp_library"] == 2      # 特徴 1（SME）と 2（SML）
    assert cov["by_source"]["database_identifier"] == 1
    assert cov["by_source"]["none"] == 2
    assert cov["with_inchikey"] == 3


def test_msp_lookup_reports_attempts(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    look = ds.inchikey_coverage["msp_lookup"]
    assert look["status"] == "used"
    assert look["file"] == "lab.msp"
    # 試行: SME 最良証拠（特徴 1: Citrulline）+ SML 4 行（Citrulline/Ornithine/Leucine/Unlisted）
    assert look["resolved"] == 3
    assert look["conflicting_names"] == 1
    assert look["unmatched"] == 1


def test_without_a_table_the_lookup_is_not_applicable(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path, msp_table=None)
    assert ds.feature_metadata["1"]["inchikey"] is None
    assert ds.inchikey_coverage["by_source"]["msp_library"] == 0
    assert ds.inchikey_coverage["msp_lookup"] == {
        "status": "not_applicable", "file": None,
        "resolved": 0, "conflicting_names": 0, "unmatched": 0}


def test_msp_resolved_features_are_not_reported_as_rdkit_misses(tmp_path, monkeypatch):
    _without_rdkit(monkeypatch)
    ds = _build_msp(tmp_path)
    # SME 由来で SMILES を持ち InChIKey が無いのは無い（特徴 1 は MSP で付いた）。
    assert not any("RDKit" in w for w in ds.validation_result.get("warnings", []))
```

- [ ] **Step 2: テストが落ちることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_dataset_state.py -q -k "msp"`
Expected: FAIL（`TypeError: build_dataset_state() got an unexpected keyword argument 'msp_table'`）

- [ ] **Step 3: 実装を書く**

`metabolomix/mztab/dataset_state.py`:

(a) import に足す（既存の `from metabolomix.mztab.identity import derive_inchikey` の下）:

```python
from metabolomix.mztab import msp_inchikey
from metabolomix.mztab.identity import _INCHIKEY_RE
```

(b) `build_dataset_state` の直前に関数を足す:

```python
#: `msp_inchikey.lookup` の結果 → `inchikey_coverage["msp_lookup"]` の数のキー。
_MSP_LOOKUP_COUNTERS = {"resolved": "resolved", "conflicting": "conflicting_names",
                        "unmatched": "unmatched"}


def _resolve_inchikey(database_identifier, inchi, smiles, msp_table, msp_counts: dict):
    """InChIKey 形の識別子 → 研究室 MSP の表 → RDKit の順で解決する。

    MSP を RDKit より前に置くのは、研究室が整備した値をそのまま使い、RDKit の有無
    （機械のアプリケーション制御）で結果を変えないため（spec §3.2）。`derive_inchikey`
    はモジュール属性として呼ぶ（既存テストがここを差し替えて RDKit 不在を再現する）。
    """
    text = str(database_identifier or "").strip()
    if text and _INCHIKEY_RE.match(text):
        return text, "database_identifier"
    if msp_table is not None:
        inchikey, outcome = msp_inchikey.lookup(msp_table, database_identifier)
        counter = _MSP_LOOKUP_COUNTERS.get(outcome)
        if counter is not None:
            msp_counts[counter] += 1
        if inchikey:
            return inchikey, "msp_library"
    return derive_inchikey(database_identifier, inchi, smiles)
```

(c) `build_dataset_state` のシグネチャを変える:

```python
def build_dataset_state(
    parse_result: dict,
    filename: str,
    source_path: str | Path,
    *,
    msp_table: "msp_inchikey.MspInchikeyTable | None" = None,
) -> DatasetState:
    """parse_mztab() の戻り値から DatasetState を構築する。

    `msp_table` は Console に渡した研究室 MSP の名前 → InChIKey 表（`mztab/loading.py`
    がジョブの実効メソッドから作る）。あれば `MspDB:<名前>` の同定に MSP の InChIKey を付ける。
    """
```

(d) `by_source` の初期化と、SME ループの `derive_inchikey` 呼び出しを置き換える:

```python
    by_source: dict[str, int] = {"database_identifier": 0, "msp_library": 0,
                                 "inchi_derived": 0, "smiles_derived": 0, "none": 0}
    msp_counts = {"resolved": 0, "conflicting_names": 0, "unmatched": 0}
```

SME ループ内:

```python
        ik, src = _resolve_inchikey(
            evidence.get("database_identifier"),
            evidence.get("inchi"),
            evidence.get("smiles"),
            msp_table, msp_counts,
        )
```

(e) `ds.inchikey_coverage = {...}` の辞書に足すのではなく、SML 由来の注釈を作ったあと（`identified_by` を入れる直前）で入れる。まず `_build_feature_annotations` 呼び出しを変える:

```python
    ds.feature_annotations, annotation_warnings = _build_feature_annotations(
        ds.sml_rows, set(ds.feature_ids), msp_table=msp_table, msp_counts=msp_counts)
```

`identified_by` を入れる行の直前に:

```python
    ds.inchikey_coverage["msp_lookup"] = {
        "status": "used" if msp_table is not None else "not_applicable",
        "file": msp_table.file_name if msp_table is not None else None,
        **msp_counts,
    }
```

(f) `_build_feature_annotations` のシグネチャと、その中の `derive_inchikey` 呼び出しを置き換える:

```python
def _build_feature_annotations(sml_rows, known_feature_ids: set, *,
                               msp_table=None, msp_counts: dict | None = None
                               ) -> tuple[dict, list[str]]:
```

関数の先頭（`by_feature` の前）に:

```python
    if msp_counts is None:
        msp_counts = {"resolved": 0, "conflicting_names": 0, "unmatched": 0}
```

ループ内の呼び出し:

```python
            ik, src = _resolve_inchikey(row.get("database_identifier"),
                                        row.get("inchi"), row.get("smiles"),
                                        msp_table, msp_counts)
```

- [ ] **Step 4: テストが通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_dataset_state.py tests/test_mztab_identity.py -q`
Expected: PASS（既存テストを含め全件）

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/mztab/dataset_state.py tests/test_dataset_state.py
git commit -m "feat(mztab): MspDB の同定に研究室 MSP の InChIKey を付ける（RDKit より優先）" > "$TEMP/msp_task2_commit.log" 2>&1; echo "exit=$?" >> "$TEMP/msp_task2_commit.log"
```

---

### Task 3: ジョブの実効メソッドから MSP を見つける（`loading.py`）

**Files:**
- Modify: `metabolomix/mztab/loading.py`（import、`_load_from_job`、新関数 `_msp_table_for_job`）
- Test: `tests/test_dataset_loading_service.py`（末尾に追記）

**Interfaces:**
- Consumes: Task 1 の `load_msp_inchikey_table`、Task 2 の `build_dataset_state(..., msp_table=...)`、`metabolomix.console.method_file.read_method_keys(path: Path) -> dict[str, str]`（キーは小文字）、`AnalysisJob.method_file: str`
- Produces: `_msp_table_for_job(job) -> tuple[MspInchikeyTable | None, str, str | None, str | None]`（表、`status`、ファイル名、警告文）。`load_dataset_state(job_path=...)` の戻り値の `ds.inchikey_coverage["msp_lookup"]["status"]` が `used` / `no_msp` / `unavailable`。直接読み（`mztab_path`）は `not_applicable`（Task 2 の既定のまま）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_dataset_loading_service.py` の末尾に追記:

```python
# ---------- MSP の INCHIKEY（spec 2026-10-10 msp-inchikey-design §3.3）----------

from metabolomix.handoff.schema import AnalysisJob
from metabolomix.mztab import msp_inchikey

_CIT = "RHGKLRLOHDJJDR-BYPYZUCNSA-N"


@pytest.fixture(autouse=False)
def _fresh_msp_memo(monkeypatch):
    monkeypatch.setattr(msp_inchikey, "_memo", {})


def _write_msp(path: Path) -> Path:
    path.write_text(f"NAME: Citrulline\nINCHIKEY: {_CIT}\nNum Peaks: 0\n\n",
                    encoding="utf-8")
    return path


def _write_mspdb_mztab(path: Path, raws) -> None:
    """msp_only の Console が書く形の同定（MspDB:<名前>、inchi は null）。"""
    lines = ["MTD\tmzTab-version\t2.0.0-M", "MTD\tmzTab-mode\tComplete",
             "MTD\tmzTab-type\tQuantification"]
    for i, raw in enumerate(raws, start=1):
        lines.append(f"MTD\tms_run[{i}]-location\t{Path(raw).resolve().as_uri()}")
        lines.append(f"MTD\tassay[{i}]-ms_run_ref\tms_run[{i}]")
        lines.append(f"MTD\tassay[{i}]\tS{i}")
    cols = "\t".join(f"abundance_assay[{i}]" for i in range(1, len(raws) + 1))
    vals = "\t".join("1000.0" for _ in raws)
    lines.append(f"SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\t{cols}")
    lines.append(f"SMF\t1\t1\t176.1030\t300.0\t{vals}")
    lines.append("SEH\tSME_ID\tdatabase_identifier\tsmiles\tinchi\tchemical_name\trank")
    lines.append("SME\t1\tMspDB:Citrulline\tNC(=O)NCCCC(N)C(O)=O\tnull\tCitrulline\t1")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _job_declaring_msp(tmp_path, msp_line: str):
    job_path, raws = _job_with_outputs(tmp_path)
    job = AnalysisJob.load(job_path)
    Path(job.method_file).write_text(f"Ion mode: Positive\n{msp_line}\n", encoding="ascii")
    _write_mspdb_mztab(Path(job.run_dir) / "msdial" / "Height_AlignmentResult_1.mzTab", raws)
    return job_path


def _no_rdkit(monkeypatch):
    monkeypatch.setattr("metabolomix.mztab.identity._rdkit_cache", None)


def test_job_with_a_declared_msp_gets_msp_inchikeys(tmp_path, monkeypatch, _fresh_msp_memo):
    _no_rdkit(monkeypatch)
    msp = _write_msp(tmp_path / "lab.msp")
    ds = load_dataset_state(job_path=_job_declaring_msp(tmp_path, f"Msp file path: {msp}"))
    assert ds.feature_metadata["1"]["inchikey"] == _CIT
    assert ds.feature_metadata["1"]["inchikey_source"] == "msp_library"
    look = ds.inchikey_coverage["msp_lookup"]
    assert look["status"] == "used"
    assert look["file"] == "lab.msp"
    assert look["resolved"] == 1


def test_relative_msp_path_is_resolved_against_the_method_file(tmp_path, monkeypatch,
                                                                _fresh_msp_memo):
    _no_rdkit(monkeypatch)
    job_path = _job_declaring_msp(tmp_path, "Msp file path: lab.msp")
    method_dir = Path(AnalysisJob.load(job_path).method_file).parent
    _write_msp(method_dir / "lab.msp")
    ds = load_dataset_state(job_path=job_path)
    assert ds.inchikey_coverage["msp_lookup"]["status"] == "used"
    assert ds.feature_metadata["1"]["inchikey"] == _CIT


def test_job_without_an_msp_declaration_is_no_msp(tmp_path, monkeypatch, _fresh_msp_memo):
    _no_rdkit(monkeypatch)
    ds = load_dataset_state(job_path=_job_declaring_msp(tmp_path, "Msp file path: "))
    assert ds.inchikey_coverage["msp_lookup"]["status"] == "no_msp"
    assert ds.feature_metadata["1"]["inchikey"] is None


def test_declared_but_missing_msp_is_unavailable_with_a_warning(tmp_path, monkeypatch,
                                                                 _fresh_msp_memo):
    _no_rdkit(monkeypatch)
    missing = tmp_path / "nowhere" / "lab.msp"
    ds = load_dataset_state(job_path=_job_declaring_msp(tmp_path, f"Msp file path: {missing}"))
    look = ds.inchikey_coverage["msp_lookup"]
    assert look["status"] == "unavailable"
    assert look["file"] == "lab.msp"
    hits = [w for w in ds.validation_result["warnings"] if "MSP" in w]
    assert len(hits) == 1
    assert str(missing.parent) not in hits[0]     # 置き場所は出さない


def test_direct_read_is_not_applicable(tmp_path, _fresh_msp_memo):
    job_path, _ = _job_with_outputs(tmp_path)
    mztab = Path(AnalysisJob.load(job_path).run_dir) / "msdial" / "Height_AlignmentResult_1.mzTab"
    ds = load_dataset_state(mztab_path=mztab)
    assert ds.inchikey_coverage["msp_lookup"]["status"] == "not_applicable"
```

- [ ] **Step 2: テストが落ちることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_dataset_loading_service.py -q -k "msp"`
Expected: FAIL（`msp_lookup` の `status` が `not_applicable`、または `inchikey` が None）

- [ ] **Step 3: 実装を書く**

`metabolomix/mztab/loading.py`:

(a) import に足す:

```python
from metabolomix.console.method_file import read_method_keys
from metabolomix.mztab.msp_inchikey import MspInchikeyTable, load_msp_inchikey_table
```

(b) `_load_from_job` の `ds = build_dataset_state(parse_result, mztab_abs.name, mztab_abs)` を置き換える:

```python
    msp_table, msp_status, msp_file, msp_warning = _msp_table_for_job(job)
    ds = build_dataset_state(parse_result, mztab_abs.name, mztab_abs, msp_table=msp_table)
    ds.inchikey_coverage["msp_lookup"]["status"] = msp_status
    ds.inchikey_coverage["msp_lookup"]["file"] = msp_file
    if msp_warning:
        ds.validation_result.setdefault("warnings", []).insert(0, msp_warning)
```

(c) `_parse_and_validate` の直前に関数を足す:

```python
def _msp_table_for_job(job) -> tuple[MspInchikeyTable | None, str, str | None, str | None]:
    """ジョブの実効メソッドが宣言する研究室 MSP から名前 → InChIKey 表を作る。

    戻り値は (表, status, ファイル名, 警告文)。status は `used`（表を作れた）/
    `no_msp`（実効メソッドが MSP を宣言していない、またはメソッド自体を読めない）/
    `unavailable`（宣言はあるがファイルを読めない）。読めなくても読み込みは止めない
    ——MSP 表は InChIKey の付け足しで、無くても従来どおり RDKit 経路で読める。
    警告と戻り値に出すのはファイル名だけ（研究室ライブラリの置き場所は出さない）。
    """
    method_text = (getattr(job, "method_file", "") or "").strip()
    if not method_text:
        return None, "no_msp", None, None
    method = Path(method_text)
    try:
        keys = read_method_keys(method)
    except (OSError, DomainError):
        return None, "no_msp", None, None
    declared = (keys.get("msp file path") or "").strip()
    if not declared:
        return None, "no_msp", None, None
    msp = Path(declared)
    if not msp.is_absolute():
        msp = method.parent / msp
    try:
        table = load_msp_inchikey_table(msp)
    except OSError:
        return None, "unavailable", msp.name, (
            f"実効メソッドが宣言する MSP（{msp.name}）を読めないため、MSP の INCHIKEY で"
            "InChIKey を付けられませんでした（RDKit 経路だけで読みました）。")
    return table, "used", table.file_name, None
```

- [ ] **Step 4: テストが通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_dataset_loading_service.py tests/test_dataset_state.py -q`
Expected: PASS（`test_loading_does_not_import_session_or_mcp_core` を含め全件）

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/mztab/loading.py tests/test_dataset_loading_service.py
git commit -m "feat(mztab): ジョブの実効メソッドの MSP から InChIKey 表を作って読み込みに渡す" > "$TEMP/msp_task3_commit.log" 2>&1; echo "exit=$?" >> "$TEMP/msp_task3_commit.log"
```

---

### Task 4: エクスポートの通り抜けと文書

**Files:**
- Test: `tests/test_msp_inchikey.py`（末尾に追記）
- Modify: `docs/output_format/mztab.md`（§13.5 の表、§13.6 の表）
- Modify: `metabolomix/mztab/identity.py`（モジュール docstring の優先順位の文）
- Modify: `docs/workflow/mztab.md`（`load_dataset_state` の呼び出し順を書いた箇所）

**Interfaces:**
- Consumes: Task 2 の `inchikey_source == "msp_library"`、`metabolomix.analysis.dataset_export.export_dataset_result(ds, result, path, curation=None) -> dict`、`tests.pipeline_fixtures.make_dataset()`、`metabolomix.analysis.dataset_service.preprocess_dataset` / `compare_dataset`、`metabolomix.analysis.export_contract.EXPORT_COLUMNS`（列名のリスト）
- Produces: 文書の更新のみ（コードの公開面は増えない）

- [ ] **Step 1: 失敗するかもしれないテストを書く（通り抜けの固定）**

`tests/test_msp_inchikey.py` の末尾に追記:

```python
def test_export_writes_msp_library_rows_with_the_contract_columns(tmp_path):
    """エクスポートは inchikey_source を通すだけ。列は契約のまま（CONTRACT_VERSION 不変）。"""
    from metabolomix.analysis import export_contract
    from metabolomix.analysis.dataset_export import export_dataset_result
    from metabolomix.analysis.dataset_service import compare_dataset, preprocess_dataset
    from tests.pipeline_fixtures import make_dataset

    ds = make_dataset()
    fid = ds.feature_ids[0]
    ds.feature_metadata.setdefault(fid, {}).update(
        {"name": "Citrulline", "inchikey": CIT, "inchikey_source": "msp_library"})
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})
    result = compare_dataset(ds, ds.sample_names[:4], ds.sample_names[4:])

    export_dataset_result(ds, result, tmp_path / "out.tsv", curation=None)

    lines = [l for l in (tmp_path / "out.tsv").read_text(encoding="utf-8").splitlines()
             if not l.startswith("#")]
    header = lines[0].split("\t")
    assert header == export_contract.EXPORT_COLUMNS
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    hit = [r for r in rows if r["spot_id"] == str(fid)]
    assert hit and hit[0]["inchikey_source"] == "msp_library"
    assert hit[0]["inchikey"] == CIT
```

Run: `C:/Python314/python.exe -m pytest tests/test_msp_inchikey.py -q -k export`
Expected: PASS（エクスポートは無改修で通る想定。通らなければ原因を調べ、仕様（列不変）を保つ最小の修正を入れる）

- [ ] **Step 2: `docs/output_format/mztab.md` を直す**

§13.5 の表の `inchikey_source` 行を次に置き換える:

```markdown
| `inchikey_source` | `database_identifier` / `msp_library` / `inchi_derived` / `smiles_derived` | どこから導出したか。`msp_library` は Console に渡した研究室 MSP の `INCHIKEY` を `MspDB:<名前>` で引いた値（RDKit 不要。spec 2026-10-10 `msp-inchikey-design`） |
```

§13.6 の表の `by_source` 行を次に置き換え、表の末尾に `msp_lookup` 行を足す:

```markdown
| `by_source` | `database_identifier` / `msp_library` / `inchi_derived` / `smiles_derived` / `none` の内訳。全特徴を 1 回ずつ数える |
| `msp_lookup` | MSP 表の使用状況。`status`（`used` / `no_msp` / `unavailable` / `not_applicable`）、`file`（MSP のファイル名だけ）、`resolved` / `conflicting_names` / `unmatched`（**照合の試行回数**: 特徴ごとの SME 最良証拠 1 件と、SML 行 × 参照特徴）。`conflicting_names` は同じ名前に別の InChIKey が MSP にあって付けなかった数 |
```

§13.6 の表の直後に段落を足す:

```markdown
**解決の順序**: `database_identifier` が InChIKey 形 → **研究室 MSP の表**（`MspDB:<名前>`）→ RDKit（InChI / SMILES）。
MSP を RDKit より前に置くのは、RDKit の有無（Windows のアプリケーション制御がネイティブ DLL を
拒否する機械がある）で結果を変えないため。MSP 表はジョブ経由の読み込み（`job_path`）でだけ作られ、
実効メソッドの `Msp file path` を読む。直接読み（`mztab_path`）では `not_applicable`。
```

- [ ] **Step 3: `metabolomix/mztab/identity.py` の docstring を直す**

先頭 docstring の 3 行目 `優先順位: database_identifier 正規表現一致 → inchi → smiles (RDKit) → None。` を次に置き換える:

```python
優先順位: database_identifier 正規表現一致 → inchi → smiles (RDKit) → None。
（`build_dataset_state` は、正規表現一致と RDKit の間に研究室 MSP の表
 `mztab/msp_inchikey.py` を挟む。この関数自体は MSP を知らない純関数のまま。）
```

- [ ] **Step 4: `docs/workflow/mztab.md` を直す**

`grep -n "load_dataset_state\|build_dataset_state" docs/workflow/mztab.md` で `_load_from_job` → `build_dataset_state` の呼び出し順を書いた箇所を探し、その間に次の 1 行（同じ字下げ・番号付けの書式に合わせる）を足す:

```markdown
metabolomix/mztab/loading.py  _msp_table_for_job()  → metabolomix/mztab/msp_inchikey.py  load_msp_inchikey_table()
```

`docs/workflow/pipeline.md` にも `load_dataset_state` の呼び出しが書かれていれば、同じ 1 行を足す。

Run: `C:/Python314/python.exe -m pytest tests/test_workflow_docs.py tests/test_readme_links.py -q`
Expected: PASS（書いた関数名・パスが実在すること）

- [ ] **Step 5: 全テストを回してコミット（バックグラウンド）**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: PASS（全件）

```bash
git add tests/test_msp_inchikey.py docs/output_format/mztab.md metabolomix/mztab/identity.py docs/workflow/mztab.md docs/workflow/pipeline.md
git commit -m "docs(mztab): inchikey_source=msp_library と msp_lookup を出力形式と workflow に書く" > "$TEMP/msp_task4_commit.log" 2>&1; echo "exit=$?" >> "$TEMP/msp_task4_commit.log"
```

---

## 完了後（計画外の運用。コントローラが行う）

1. `main` へ `--no-ff` マージ: `git merge --no-ff feat/msp-inchikey -m "Merge feat/msp-inchikey: 研究室 MSP の INCHIKEY で mzTab-M の同定に InChIKey を付ける（RDKit 不要）"`。worktree を片付け、ローカルブランチを削除。
2. `docs/HISTRY.md` に経緯を追記、`docs/task.md` を更新（どちらも追跡外。メイン作業ツリーで）。フロー（入口・stage 列・前提状態）は変わらないので vault の流れ図は直さない。
3. MCP サーバの再起動はユーザーに依頼する。
4. reanalysis-study の 4 本の run は、`comparison_id` を変えた同じ比較で `pipeline_resume`（spec §3.5）。
