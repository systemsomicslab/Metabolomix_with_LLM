# MSP の INCHIKEY で mzTab-M の同定に InChIKey を付ける

- 日付: 2026-10-10
- 状態: 設計合意済み（ユーザー承認 2026-10-10。書面レビュー待ち）
- 由来: reanalysis-study（`C:\Users\yuu18\reanalysis-study`）の慣らし運転 ST003841。`issues.md` 2026-10-10
- 前提: spec 2026-10-07 `library-mode-msp-only-design.md`（マージ `f41f415`）

## 1. 目的

`library_mode="msp_only"` の Console 実行（研究室 MSP だけで同定）では、mzTab-M の同定行が
`database_identifier = MspDB:<MSP の NAME>`、`smiles` あり、`inchi` は常に空、で出る。InChIKey は
`derive_inchikey` が RDKit で SMILES から作るしかない。ところが Windows のアプリケーション制御
（Smart App Control）が RDKit のネイティブ DLL を拒否する機械では RDKit が使えず、**全特徴の InChIKey が
空**になる。結果として:

- `dataset_export_differential`（pipeline の `export:<id>` 工程）が `EXPORT_BACKGROUND_EMPTY` で何も書かない。
- 下流の massbank-context（InChIKey を結合キーにするパスウェイ解析）に何も渡らない。

実測（ST003841、陽イオン、37,156 特徴）: SML 29,374 行・SME 1,165 行が SMILES を持つが、InChIKey は 0 件。

研究室 MSP は各エントリに `INCHIKEY:` を持つ（pos で 261 万エントリ中ほぼ全件）。Console に渡した MSP から
名前で InChIKey を引けば、RDKit なしで同じ情報に届く。

成功の条件:

1. Console の実効メソッドが `Msp file path` を持つジョブを `load_dataset_state(job_path=...)` で読むと、
   `MspDB:<名前>` の同定に MSP の InChIKey が付き、`inchikey_source = "msp_library"` になる。
2. RDKit の有無で 1 の結果が変わらない（機械に依存しない）。
3. MSP が見つからない・読めないときは、今まで通りの動き（RDKit 経路）に戻り、理由を警告に残す。
4. 差次的エクスポートの列は変わらない（`CONTRACT_VERSION` は上げない）。

## 2. 確認した事実

- `derive_inchikey`（`metabolomix/mztab/identity.py`）の順序は「`database_identifier` が InChIKey 形 → InChI →
  SMILES（RDKit）→ None」。`build_dataset_state`（`metabolomix/mztab/dataset_state.py`）が SME の最良証拠と SML
  行の両方でこれを呼び、`inchikey_coverage.by_source` を数える。
- MS-DIAL の mzTab-M には InChIKey の列が無い（SML: `database_identifier` `chemical_formula` `smiles` `inchi`
  `chemical_name` …、SME も同様）。`database_identifier` は `MspDB:<名前>` 形式で、名前に確度接頭辞
  （`no MS2:` `low score:`）は付かない（`chemical_name` には付く）。
- ジョブ（`analysis-job.json`）の `software.method_file` は実効メソッド（pipeline では
  `<pipeline_root>/inputs/effective-method.txt`）を指し、その `Msp file path:` が Console に渡した MSP の
  絶対パス。
- `metabolomix/mztab/loading.py` は `session_state` も `mcp_core` も import しない契約
  （`tests/test_dataset_loading_service.py::test_loading_does_not_import_session_or_mcp_core`）。
  `metabolomix/library/store.py` は `mcp_core` を import するので、ここから SQLite キャッシュは使えない。
- 研究室 MSP（pos）の `NAME` と `INCHIKEY` 行だけを読む走査は 43 秒（261 万エントリ、29.4 万の名前）。
  同じ名前に別の InChIKey があるものが 4,930 名（1.7%）。
- massbank-context は `inchikey_source` の値を解釈しない（列の存在だけを契約にしている）。既存でも
  `curation` `arf2` という文書外の値が使われている。

## 3. 設計

### 3.1 名前 → InChIKey の表（新規 `metabolomix/mztab/msp_inchikey.py`）

```python
@dataclass(frozen=True)
class MspInchikeyTable:
    by_name: dict[str, str]        # 一意に決まる名前だけ
    conflicting: frozenset[str]    # 同じ名前に別の InChIKey があった名前
    sha256: str | None             # 読んだ MSP の内容 hash（来歴用。計算しない場合は None）

def load_msp_inchikey_table(path: Path) -> MspInchikeyTable: ...
def lookup(table: MspInchikeyTable, database_identifier: str | None) -> str | None: ...
```

- 読むのは各エントリの `NAME:` と `INCHIKEY:` 行だけ（キーは大文字小文字を区別しない、値は前後の空白を除く）。
  ピーク行は解さない。文字コードは `library/msp.py` と同じく UTF-8 を試し、壊れたバイトは置換する。
- `INCHIKEY` が空、または `_INCHIKEY_RE`（`identity.py`）に合わないものは表に入れない。
- 同じ名前に別の InChIKey が現れたら `by_name` から外して `conflicting` に入れる（名前で決めない）。
- `lookup` は `database_identifier` が `MspDB:` で始まるときだけ、その後ろを名前として引く。それ以外は None。
- 同じプロセスでは `(絶対パス, サイズ, 更新時刻)` を鍵に表を使い回す（pipeline の worker は 1 回、MCP サーバは
  同じ MSP を何度も読まない）。
- stdlib だけで書き、`mcp_core` / `session_state` / `library.store` を import しない。
- `sha256` は計算しない（1 GB 級を 2 回読まないため。来歴の MSP hash は pipeline が既に固定している）。

### 3.2 解決の順序（`build_dataset_state`）

`build_dataset_state(parse_result, filename, source_path, *, msp_table: MspInchikeyTable | None = None)`。
SME 証拠と SML 行のそれぞれで:

1. `database_identifier` が InChIKey 形 → `database_identifier`（今まで通り）
2. **`msp_table` があり、`lookup` が InChIKey を返す → `msp_library`**
3. InChI / SMILES を RDKit で → `inchi_derived` / `smiles_derived`（今まで通り）
4. どれも無い → `none`

2 を RDKit より前に置くのは、研究室が整備した値をそのまま使い、RDKit の有無で結果を変えないため
（成功条件 2）。`derive_inchikey` 自体は変えない（純関数のまま）。順序の合成は `dataset_state.py` に
小さな関数を 1 つ置いて行い、SME と SML の 2 か所から呼ぶ。

`inchikey_coverage.by_source` の初期キーに `msp_library` を足す（全特徴を 1 回ずつ数える不変条件は保つ）。
さらに `inchikey_coverage["msp_lookup"]` を足す:

| キー | 意味 |
|---|---|
| `status` | `used`（表を使った）/ `no_msp`（ジョブにMSPの宣言が無い）/ `unavailable`（宣言はあるが読めない）/ `not_applicable`（直接読み） |
| `file` | MSP のファイル名だけ（置き場所は載せない。研究室ライブラリの非開示） |
| `resolved` | `MspDB:` の同定のうち、表で InChIKey が付いた数 |
| `conflicting_names` | 名前が衝突して付けなかった数 |
| `unmatched` | `MspDB:` なのに表に名前が無かった数 |

RDKit 不在の警告（`derivable_but_missing`）は、MSP 経路で付いたものを差し引いて数える（付いたものを
「取りこぼし」と言わない）。

### 3.3 MSP の見つけ方（`metabolomix/mztab/loading.py`）

`_load_from_job` が、ジョブの `software.method_file` を読み、`Msp file path:` の値（空でない絶対パス）を取る。

- 宣言が無い・空 → `msp_lookup.status = "no_msp"`、表なし（今まで通り）。
- 宣言はあるがファイルが無い・読めない → `status = "unavailable"`、表なし、警告を 1 件。
- 読めた → `load_msp_inchikey_table` → `build_dataset_state(..., msp_table=table)`。

`_load_from_mztab`（直接読み）は MSP を知らないので `not_applicable`。引数で MSP を渡す口は今回足さない（YAGNI）。

メソッドファイルの読み取りは `metabolomix.console.method_file.read_method_keys` を使う（キーは小文字化される）。
`console.method_file` の依存は `core.atomic_io` と `core.user_config`（どちらも stdlib だけの leaf）で、
import しても `mcp_core` / `session_state` は読み込まれない（2026-10-10 に確認）。既存の
`test_loading_does_not_import_session_or_mcp_core` がこれを縛り続ける。

### 3.4 エクスポートと文書

- エクスポートの列は不変。`inchikey_source` に新しい値 `msp_library` が出るだけ。`name_source` は今まで通り
  （`mztab_sme` / `mztab_sml`）。
- `docs/output_format/mztab.md`（MCP リソース `lipidmix://docs/output-format/mztab`）の §13.5 と §13.6 に
  `msp_library` と `msp_lookup` を足す。`identity.py` の docstring の優先順位の記述も直す。
- `docs/workflow/` で `load_dataset_state` の呼び出し順を書いている箇所に、MSP 表の構築を足す。
- `docs/HISTRY.md` に経緯、`docs/task.md` は該当があれば更新。

### 3.5 既存 run への適用（運用。コードは変えない）

`pipeline_resume` は `load_dataset`〜`pca` を毎回作り直すが、成功済みの `export:<id>` は作り直さない。
改修前に `EXPORT_BACKGROUND_EMPTY` で終わった run は、`comparisons` の `comparison_id` だけを変えた
同じ比較で `pipeline_resume` すると、新しい ID の差次的解析と書き出しが作られる（Console は回し直さない）。

## 4. 対象外

- RDKit を使えるようにすること（機械のセキュリティ設定。ユーザーの判断）。
- Text DB・LBM 由来の同定の InChIKey（`MspDB:` 以外は今まで通り）。
- 名前の衝突を構造（SMILES・式）で解くこと。衝突は数えるだけ。
- 直接読み（`mztab_path`）に MSP を渡す引数。
- resume が成功済みの export を作り直す仕組み。

## 5. テスト（TDD）

研究室ライブラリには触れない。MSP は `tmp_path` に数エントリの小さなファイルを作る。

1. `load_msp_inchikey_table`: 名前→InChIKey、空・不正形の除外、同名同キーの重複は 1 件、同名別キーは
   `conflicting`、キーの大文字小文字。
2. `lookup`: `MspDB:<名前>` だけを引く。接頭辞なし・別 DB・None は None。
3. `build_dataset_state(msp_table=...)`: SME・SML の両方で `msp_library` が付く。InChIKey 形の
   `database_identifier` が MSP より優先、MSP が RDKit より優先（RDKit の memo を偽物に差し替えて確かめる）、
   `by_source` の総数が特徴数と一致、`msp_lookup` の各数。
4. `load_dataset_state(job_path=...)`: 合成ジョブの実効メソッドの `Msp file path` から表を作る。宣言なし →
   `no_msp`、ファイルなし → `unavailable` と警告、直接読み → `not_applicable`。
5. `loading.py` が `mcp_core` / `session_state` を import しない既存テストが緑のまま。
6. エクスポート: `msp_library` の行が書き出され、列は `export_contract` のまま。
