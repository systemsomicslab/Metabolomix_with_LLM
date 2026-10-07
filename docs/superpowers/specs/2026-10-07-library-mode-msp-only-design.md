# library_mode="msp_only": LBM を使わず研究室 MSP だけで Console を回す

- 日付: 2026-10-07
- 状態: 設計合意済み（ユーザー承認 2026-10-07）
- 由来: reanalysis-study（`C:\Users\yuu18\reanalysis-study`）plan `docs/plans/2026-10-05-phase0-setup.md` の Task 9（要求 R1〜R5）
- 前提: task_a4b18144 のマージ `27884ee`（相対パス宣言の絶対化）と、設定ファイル（spec 2026-10-06 `user-config-file-design.md`、マージ `9b6935b`）

## 1. 目的

reanalysis-study のライブラリ契約は「同定はリピドミクス・メタボロミクスとも研究室 MSP（VS21、測定の極性で
pos / neg を選ぶ）だけ。LBM は使わない」。いまの ms-data-parser は、リピドミクスでは `resolve_lbm` が LBM を
必ず解決して実効メソッドに書き込み、止める手段が無い。MSP は、メソッドファイルが宣言しているときだけ
絶対パスにして Console に渡しており、設定ファイルの `[library] msp_*` は `library_load` しか使っていない。

この spec は、契約どおりの Console 実行を 1 つの切り替えで選べるようにする。成功の条件:

1. 実効メソッドの `Lbm file path` が空。
2. `Msp file path` が、測定の極性に合う研究室 MSP の絶対パス。
3. 1・2 を満たせないときは Console を起動せず、機械可読な封筒で止まる（同定 0 件のまま黙って走らせない）。

既定の動き（`library_mode="auto"`）は一切変えない。

## 2. 確認した事実（この設計の根拠）

- Console は `Msp file path` を `Target omics` に関係なく読み、LBM の照合は `Target omics` に関係なく脂質モードで
  走る（MsdialWorkbench afd5f95 の `ConfigParser.cs` / `CommonProcess.cs` `ParseLibraries` / `LcmsProcess.cs`）。
  したがってリピドミクスでも、`Lbm file path` を空にして `Msp file path` だけを渡せば MSP だけで同定される。
- 同定用ライブラリのキーは「実在するときだけ読む」。空・不在なら警告なしで同定なし（`IsFileExist`）。
- LC-MS の Console は相対パスをメソッドファイル基準で解決しない。パスは絶対かつ ASCII で書く。
- `write_effective_method_file` は、値が空の上書きを `key: ` の行として書き、原本に行が無いキーは末尾に足す。
  同じキーの行が複数あれば全部差し替える。
- pipeline v2（`pipeline-request.v2`）はプロファイルで依存を明示する方式で、`method_file` / `lbm_file` を受け付けない。

## 3. インタフェース

### 3.1 新しい引数

`console_plan`・`console_method_template`・pipeline v1 の要求（`pipeline_plan` / `pipeline_run` の `request`、
`analysis-request.json`）に次の 2 つを足す。

| 引数 | 値 | 既定 |
|---|---|---|
| `library_mode` | `"auto"` / `"msp_only"` | `"auto"` |
| `msp_file` | MSP の明示パス（`msp_only` のときだけ意味を持つ） | 無し |

- `"auto"`: 今の動きのまま（`resolve_lbm` による LBM の解決と自動補完、メソッドの宣言の絶対化）。
- `"msp_only"`: §3.2 の動き。
- pipeline v1 の要求では、`library_mode` と `msp_file` の明示の null を拒否する（`PIPELINE_REQUEST_INVALID`）。
  不正な `library_mode` も同じく拒否する。`library_mode` は resume で変更できない（`UPDATABLE` に入れない）。
- `library_mode="auto"` で `msp_file` を渡したら止める（Console 系は `LIBRARY_MODE_CONFLICT`、pipeline v1 は
  `PIPELINE_REQUEST_INVALID`）。黙って無視すると、研究室 MSP を使ったつもりで LBM の解析が走る。
- 要求の内容 hash（`request_fingerprint`）は、`library_mode="auto"` と `msp_file=None`（どちらも既定値）を
  hash の入力から外す。既定値のキーが増えただけで既存の要求の hash が変わると、過去の run への再送が
  `IDEMPOTENCY_CONFLICT` になり、完了済み run の再利用も外れるため（`tests/test_metabolomics_stages.py`
  `test_v1_request_fingerprint_is_unchanged` が実測値で縛っている）。
- pipeline-request.v2 に `library_mode` / `msp_file` が来たら、`method_file` / `lbm_file` と同じ移行案内付きで拒否する。

### 3.2 `library_mode="msp_only"` の動き

1. **LBM を切る。** `resolve_lbm` を呼ばない（ビルドツリー・`[msdial] lbm` / `MSDIAL_LBM`・実行体のフォルダ
   からの自動補完もしない）。実効メソッドの `Lbm file path` を空にする（原本が宣言していても消す）。
   `lbm_file` と同時に指定されたら `LIBRARY_MODE_CONFLICT` で止める。
2. **MSP を決める。** 解決順は「`msp_file` 引数 → 極性の設定」。極性の設定は
   `user_config.get_setting("library.msp_positive")`（polarity=positive）/ `"library.msp_negative"`（negative）
   で引く（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → `lipidmix.local.toml`）。
   メソッドファイルの `Msp file path` の宣言は**使わない**（研究室 MSP 以外の MSP を宣言した GUI パラメータを
   流用しても、契約から外れないようにする。ユーザー決定 2026-10-07）。
   解決したパスは絶対パスで実効メソッドの `Msp file path` に書く。`msp_file` の扱いは既存の `lbm_file` と
   同じにする: Console 系は `Path(msp_file).expanduser()` で実在を確かめてから絶対化し、pipeline v1 は相対なら
   `source_root` 基準で解く。
3. **ほかの同定用の宣言を消す。** 原本が宣言している次のキーを、実効メソッドで空にする:
   `Text DB file path`、`SETTINGS_PATH_KEYS` の 7 キー（MSP / Text の注釈器設定表）。消したキーは戻り値の
   `removed_declarations`（原本での綴りのリスト）に出す。同定に使わないキー（`Isotope text DB file path`、
   ターゲット検出用・RT 補正用のキー）は消さない（`auto` と同じく絶対化だけする）。
   `removed_declarations` に出すのはキーの正準の綴り（`TEXT_DB_KEY` / `SETTINGS_PATH_KEYS` の値。原本の行は
   大文字小文字を問わず差し替わる）。原本が `Lbm file path` を空でなく宣言していたら、それも含める。
   pipeline v1 では、メソッドの宣言参照の厳格検査（`resolve_method_references`。対象は LBM）を `msp_only` では
   行わない（使わない LBM の宣言が実在しないだけで止めない）。
4. **止める条件**（どれも Console を起動しない・ジョブを作らない）:

| コード | 条件 |
|---|---|
| `MSP_NOT_CONFIGURED` | `msp_file` が無く、その極性の設定（環境変数・設定ファイル）も無い |
| `MSP_NOT_FOUND` | `msp_file` か設定が指すファイルが無い |
| `LIBRARY_MODE_CONFLICT` | `msp_only` と `lbm_file` の同時指定 |
| `CONFIG_INVALID`（既存） | 設定ファイルが読めない |
| `METHOD_ENCODING_UNSUPPORTED`（既存） | MSP のパスが ASCII でない |

   メッセージには、設定ファイルの節とキー（`[library] msp_positive` など）と環境変数名を示す
   （`user_config.missing_hint` / `describe_missing` を使う）。

### 3.3 戻り値

**研究室 MSP の置き場所（ディレクトリ）は戻り値・エラー文に載せない。** ms-data-parser の既存の規約
（`core/path_resolvers.py` `LibraryPathError`: 戻り値は LLM の文脈に入るので、ファイル名・設定キー・環境変数名・
設定ファイルのパスだけを載せる）に合わせる。絶対パスは実効メソッドと、データ側（`pipeline_root`）の入力記録にだけ置く。

`console_plan` と `console_method_template`:

- `"library_mode"` を足す。
- `msp_only` のとき `"lbm": {"path": null, "source": "disabled"}`。
- `"msp": {"file": <ファイル名>, "source": "argument" | "env" | "config_file"}` を足す（`auto` では
  `{"file": null, "source": "not_used"}`。`auto` で原本が `Msp file path` を宣言していても、ここは
  `not_used` のままにする——`auto` の MSP は今までどおり宣言の絶対化だけで、解決はしていないため）。
- `"removed_declarations": [...]` を足す（`auto` では空リスト）。

pipeline v1:

- 入力計画（`inspect_inputs` の戻り値。データ側に保存される）に `"library_mode"`、
  `"msp": {"path", "file", "sha256", "source"}`、`"removed_declarations"` を足す。
  `msp_only` でないときは `{"path": null, "file": null, "sha256": null, "source": "not_used"}`。
- `pipeline_plan` の resolved には `library_mode` と `msp: {"file", "sha256", "source"}` を載せる（R4。`path` は載せない）。
- 計画の `lbm` は、`msp_only` では `{"path": null, "sha256": null}`。
- 入力の固定（snapshot）にも `msp` を写し、`verify_inputs` の内容検査（method / lbm / exe）に `msp` を足す。
  計画後に同じパスの MSP の中身が変われば `INPUT_CHANGED`（`which: "msp"`）。
- `_plan_fingerprint` は、`msp.sha256` が**あるときだけ** `msp_sha256` をハッシュの入力に足す。
  これで既存の `auto` の計画の指紋は変わらず、受け付けの冪等性と再開が壊れない。
- 実行レポート（`report.py` の Method / LBM / Version 節）に `library_mode`、`msp_file`、`msp_sha256`、
  `msp_source` の行を足す（レポートは LLM が読むので、ここもファイル名だけ）。

### 3.4 コスト

研究室 MSP は NAS 上で約 2 GB（pos）/ 約 0.5 GB（neg）。pipeline v1 は計画時と実行時の検査で sha256 を計算する
ので、それぞれ数十秒かかりうる。キャッシュは作らない（ユーザー了承 2026-10-07）。`console_plan` /
`console_method_template` は sha256 を計算しない。

## 4. 実装の単位

| 単位 | 変更 |
|---|---|
| `metabolomix/core/user_config.py` | `MSP_SETTING_KEYS = {"positive": "library.msp_positive", "negative": "library.msp_negative"}`（極性 → 設定キーの正準。`core/path_resolvers.py` の `LIBRARY_SETTING_KEYS` はこれを参照する形に変える。pipeline の worker が `mcp_core` を import しないよう、leaf の user_config に置く） |
| `metabolomix/console/method_file.py` | `MspResolution`（`path / source / error_code / message`）、`resolve_msp(override, polarity, msp_setting)`（純粋関数。`resolve_lbm` と同じ形）、`IDENTIFICATION_CLEAR_KEYS`（`TEXT_DB_KEY` と `SETTINGS_PATH_KEYS`）、`msp_only_overrides(method_keys, msp_path) -> (overrides, removed)` |
| `metabolomix/tools/console_tools.py` | `console_plan` / `console_method_template` に引数を足し、`msp_only` のとき `resolve_lbm` の代わりに上の 2 関数を使う。MSP の設定を引く小関数（`_lbm_setting` と同じ形で `CONFIG_INVALID` を封筒にする） |
| `metabolomix/pipeline/request.py` | `_TOP_LEVEL_KEYS`・`_NULL_REJECTED_TOP_LEVEL_KEYS`・検証・`_pick` の既定（`library_mode="auto"`） |
| `metabolomix/pipeline/inputs.py` | `inspect_inputs` の分岐、`_resolve_msp_pinned`（sha256 付き）、snapshot への写し、`verify_inputs` |
| `metabolomix/pipeline/service.py` | `_plan_fingerprint`、resolved への `msp` |
| `metabolomix/pipeline/report.py` | Method / LBM / Version 節 |
| `metabolomix/tools/pipeline_tools.py` | docstring（要求キーの説明） |

`overrides` の組み立ては `msp_only` でも `relative_path_overrides` を先に作り、その上に
`msp_only_overrides` の結果を重ねる（消すキー・LBM・MSP が勝つ）。ASCII 検査（`_non_ascii_override_error`）は
重ねた後に行う。

## 5. テスト（TDD）

- 研究室 MSP の実パス・実ファイル名はテストに書かない（R5）。tmp に偽の `.msp` を作り、設定は tmp の toml を
  `LIPIDMIX_CONFIG` で指す（`tests/conftest.py` の `_isolate_external_asset_settings` の流儀）。
- `method_file`（純粋関数）:
  - `resolve_msp`: 引数優先 / 極性で pos・neg を選ぶ / 設定が env と config_file のどちらから来たかを `source` に出す /
    `MSP_NOT_CONFIGURED` / `MSP_NOT_FOUND`。
  - `msp_only_overrides`: `Lbm file path` が空・`Msp file path` が絶対パス・宣言された Text DB と注釈器設定表が空・
    宣言されていないキーは足さない・`removed` が原本の綴り。
  - 空の値の上書きが `write_effective_method_file` で `key: ` になり、重複行も全部空になる。
- `console_plan` / `console_method_template`（fake Console の既存流儀）:
  - リピドミクスで、ビルドツリーや `[msdial] lbm` に LBM があっても、実効メソッドの `Lbm file path` が空。
  - メソッドが別の MSP を宣言していても、実効メソッドの `Msp file path` は設定の MSP。
  - 止める条件 4 つ（`MSP_NOT_CONFIGURED` / `MSP_NOT_FOUND` / `LIBRARY_MODE_CONFLICT` / 非 ASCII）で
    ジョブが作られない。
  - `auto` の既存テストは変更なしで全部通る。
- pipeline v1:
  - 要求の検証（不正値・null・resume での変更不可）。
  - 計画に `msp`（sha256 付き）が載り、`lbm` が空。
  - 計画後に MSP の中身を変えると `INPUT_CHANGED`（`which: "msp"`）。
  - `auto` の計画の指紋が変更前と同じ（既存の指紋テスト、無ければ追加）。
- `tests/test_gitignore_library.py` を含め全テストが緑。

## 6. 文書

- `USAGE.md`: `console_plan` / `console_method_template` の引数、pipeline v1 の要求キー。
- `docs/workflow/` の Console・pipeline の節（`tests/test_workflow_docs.py` が名前の実在を検証する）。
- vault の流れ図は直さない（CLAUDE.md: ツールの引数だけの変更は対象外。既定の動きも変えない）。
- 記録は main ツリーの `docs/HISTRY.md` と `docs/task.md`。task.md から Task 9 の項目を消す。
- マージ後、reanalysis-study の手順書の工程 5 に `library_mode="msp_only"` を書き、版を v3 に上げる
  （reanalysis-study 側の作業）。

## 7. 対象外

- pipeline v2（プロファイル方式）。必要になれば、プロファイルの依存に `msp` を足す別 spec にする。
- `library_mode` を設定ファイルで既定にすること（この機械の他の解析まで巻き込むため、ユーザーが選ばなかった）。
- MSP の sha256 のキャッシュ。
- `library_load` の変更（既に設定ファイルを読んでいる）。
