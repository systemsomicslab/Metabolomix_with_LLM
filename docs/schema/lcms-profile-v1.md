# lcms-profile.v1 — 検証済みLC–MSプロファイル契約

検証済みLC–MSメタボロミクス実行経路（spec
[2026-09-15-validated-lcms-metabolomics-design.md](../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md)
§5, §5.1）が使う測定法・処理・統計・QCの固定プロファイルと、その検証証明書
`lcms-profile-validation.v1` の契約。実装は
[`metabolomix/console/profile_schema.py`](../../metabolomix/console/profile_schema.py)。

> このモジュールはファイルを読まない。profileファイル・証明書ファイルの読込、
> 依存ファイル・method・証明書出力の実体hash計算は呼び出し側（後続タスク）の
> 責務で、ここでは渡されたdictを構造的に検証するだけ。

## 全体像

- `validate_profile(data: dict) -> dict`: `lcms-profile.v1` を検証し、正規化した
  新しいdictを返す。未知キー・型違反・NaN/Infinity・参照整合性違反（存在しない
  target_id/recipe_id参照、内部標準対応の重複・循環、既定matrix recipeの欠落等）
  はすべて `DomainError("PROFILE_INVALID", ...)`。
- `profile_content_hash(data: dict) -> str`: `validation` キーを除いたprofile本体の
  canonical hash（`metabolomix.core.atomic_io.canonical_hash`）。証明書の付与や
  `status` の変更（draft→validated）はprofile本体の内容identityを変えない。
- `validate_certificate(profile: dict, certificate: dict, observed_hashes: dict) -> None`:
  証明書が「このprofileを正当にvalidatedと名乗らせる」ものかを検証する。不一致・
  不足はすべて `DomainError("PROFILE_VALIDATION_INVALID", ...)`。署名認証基盤は
  作らない（spec §5「署名認証基盤は作らない」）。

**draft と validated の違い**は `validation.status` だけで決まらない。
`validate_profile` は `status="validated"` のprofileに `certificate_path` /
`certificate_sha256` が構造的に埋まっていることは要求するが、それらが実在の
正当な証明書を指しているかどうかは検証しない。実在性・整合性は
`validate_certificate` が別途、証明書の中身と `observed_hashes`（呼び出し側が
実測したhash）を突き合わせて検証する。**profileを手書きで `validated` にしても、
対応する証明書が検証を通らなければ信頼されない。**

## 「未記載の測定条件」の表現

spec は測定条件の未記載値に「値nullと理由」を持たせると定めるが、その
value/reasonの対を実際に保持する場所は `acquisition`/`processing`/`software`
自身ではなく **`evidence`**（field pathキーのdict）である。条件フィールド自身は
素の値（nullable な項目は素の `null`）を持ち、「なぜnullなのか」「どの根拠tier
（`paper_explicit` 等）に基づくか」は `evidence["<field path>"]` 側に一元化する。
条件値と出典証跡を分離することで、未記載理由を書き忘れたまま値だけ埋める事故を
「evidenceの`reason`必須（valueがnullのとき）」で機械的に防げる。

## トップレベル

未知キー・欠落キーはどちらも拒否（`PROFILE_INVALID`）。13キーちょうど。

| フィールド | 型 | null可否 | 内容 |
|---|---|---|---|
| `schema` | str | 不可 | 固定値 `"lcms-profile.v1"` |
| `profile_id` | str | 不可 | 安全なslug（`^[A-Za-z0-9_-]+$`） |
| `revision` | int | 不可 | 正の整数（`bool`不可、`1.5`や`"1"`は拒否）。更新時は新revision |
| `omics` | str | 不可 | 初期v1は `"metabolomics"` のみ（lipidomicsは互換adapter側で解決し、このschemaには乗らない） |
| `acquisition` | object | 不可 | 下記参照 |
| `software` | object | 不可 | 下記参照 |
| `processing` | object | 不可 | 下記参照 |
| `analysis_recipe` | object | 不可 | 下記参照 |
| `feature_targets` | object | 不可 | `target_id`キーのdict。下記参照 |
| `matrix_recipes` | object | 不可 | `recipe_id`キーのdict。`"default"`キー必須。下記参照 |
| `qc_policy` | object | 不可 | metricキーのdict。下記参照 |
| `evidence` | object | 不可 | field pathキーのdict。下記参照 |
| `validation` | object | 不可 | 下記参照 |

## `acquisition`

8キーちょうど必須（値は個別にnull可）。

| フィールド | 型 | null可否 | 内容 |
|---|---|---|---|
| `separation` | str | 不可 | 固定値 `"lc"` |
| `acquisition_type` | str | 不可 | 固定値 `"dda"` |
| `polarity` | str | 不可 | `"positive"` / `"negative"` |
| `instrument` | str | 可 | 装置名。未記載は`null`（`evidence`側に理由） |
| `lc` | object | 不可 | 下記。個々のキーはnull可・未指定キーは`null`埋め |
| `ms_range` | object | 不可 | 下記 |
| `sample_matrix` | str | 可 | 試料マトリクス |
| `scope` | str | 不可 | 適用範囲の記述（例:「単一LC条件・positive・DDA・peak height限定」）。**`validate_certificate`が突合するのは`validation.scope`であり、`acquisition.scope`はどことも機械比較されない**（詳細は「`validation`」節・「証明書 `lcms-profile-validation.v1`」節を参照） |

`acquisition.lc`（キーは未指定なら`null`埋め。未知キーのみ拒否）:

| フィールド | 型 | 内容 |
|---|---|---|
| `column` | str/null | カラム仕様 |
| `mobile_phase_a` / `mobile_phase_b` | str/null | 移動相 |
| `flow_rate_ul_min` | float/null | 流量（>0） |
| `column_temperature_c` | float/null | カラム温度（>0） |
| `gradient_profile` | str/null | 溶媒グラジエントの記述 |

`acquisition.ms_range`（4キーちょうど必須）:

| フィールド | 型 | 内容 |
|---|---|---|
| `ms1_low_mz` / `ms1_high_mz` | float | MS1範囲（>0、high>low） |
| `ms2_low_mz` / `ms2_high_mz` | float/null | MS2範囲。両方null、または両方指定（high>low）のいずれか |

## `software`

4キーちょうど必須。すべて空でない文字列。

| フィールド | 内容 |
|---|---|
| `msdial_version` | MS-DIALのバージョン文字列 |
| `executable_path` | Console実行ファイルのパス |
| `executable_sha256` | 実行ファイルのSHA-256（64桁小文字16進） |
| `adapter_version` | このプロファイルを解釈するadapterのバージョン |

## `processing`

5キーちょうど必須。

| フィールド | 型 | 内容 |
|---|---|---|
| `method_path` | str | メソッドファイルの原本パス |
| `method_sha256` | str | メソッドファイルのSHA-256 |
| `measure` | str | 固定値 `"peak_height"` |
| `dependencies` | list | 下記。`dependency_id`は一意 |
| `effective_settings` | object | 実効設定snapshot（**開いた辞書**。MS-DIALの実パラメータキー集合を固定しない唯一の例外。NaN/Infinity・非JSON型は依然として拒否される） |

`processing.dependencies[]`（6キーちょうど必須。spec §5.1）:

| フィールド | 型 | 内容 |
|---|---|---|
| `dependency_id` | str | 安全なslug。プロファイル内で一意 |
| `kind` | str | `msp` / `text_identification` / `rt_reference` / `lbm` のいずれか |
| `method_key` | str | メソッドファイル側のキー名（例: `"Lbm file path"`） |
| `path` | str | 依存ファイルの原本パス |
| `sha256` | str | 依存ファイルのSHA-256 |
| `required` | bool | 起動前必須か |

metabolomicsでLBM（`kind="lbm"`）を必須にしない（`dependencies`に含めなくてよい）。

## `analysis_recipe`

**`statistics` と `internal_standards` の2キーだけ**を持つ（数値前処理そのものは
`matrix_recipes` が唯一の保存先で、ここには含めない）。

`analysis_recipe.statistics[]`（spec §6.2の既定統計定義。共通5キー
`statistic_id` / `kind` / `matrix_recipe_id` / `transform` / `feature_scope` に
`kind`別の追加キーがちょうど加わる。未知キーは`kind`ごとに拒否）:

| フィールド | 型 | 内容 |
|---|---|---|
| `statistic_id` | str | 安全なslug。プロファイル内で一意（重複は`PROFILE_INVALID`） |
| `kind` | str | `pca` / `welch` / `anova_tukey` |
| `matrix_recipe_id` | str | `matrix_recipes`に実在するrecipe_id（不存在参照は`PROFILE_INVALID`） |
| `transform` | str | `none` / `log2` |
| `feature_scope` | object | `{"mode": "all_eligible"}` または `{"mode": "targets", "target_ids": [...]}`（`target_ids`は`feature_targets`のキーを参照） |

`kind`別追加キー:

| `kind` | 追加キー |
|---|---|
| `pca` | `scaling`（`none`/`autoscale`）, `n_components`（正整数） |
| `welch` | `reference_group`, `test_group`（異なる2群）, `q_threshold`（0<x<1）, `log2fc_threshold`（非負有限） |
| `anova_tukey` | `groups`（重複のない3群以上）, `alpha`（0<x<1） |

`analysis_recipe.internal_standards[]`（`target_id` / `standard_target_id`の
2キーちょうど）:

- 両方とも `feature_targets` に実在する `target_id` を参照する（不存在参照は
  `PROFILE_INVALID`）。
- `target_id == standard_target_id`（自己参照）は拒否。
- 同じ `target_id` を2度対応付ける（重複マッピング）ことは拒否。
- 対応付けが循環する（例: A→B→A）ことは拒否。

## `feature_targets`

`target_id`（安全なslug）をキーとするdict。各要素は8キーちょうど（spec §8.1）:

| フィールド | 型 | 内容 |
|---|---|---|
| `compound_identifiers` | object | 名前以外の識別子を最低1つ含む必要がある（`name`だけは不可）。許可キー: `name` / `inchikey` / `formula` / `cas` / `hmdb_id` / `kegg_id` / `pubchem_cid` |
| `adduct` | str | 例: `"[M+H]+"` |
| `charge` | int | 0でない整数 |
| `expected_mz` | float | >0 |
| `mz_tolerance_ppm` | float | >0（負値・0は`PROFILE_INVALID`） |
| `expected_rt_min` | float | >0 |
| `rt_tolerance_min` | float | >0 |
| `required_evidence` | list | 非空。各要素は下記。`kind`の重複は不可 |

`required_evidence[]`（`kind`ごとに許可キーが決まる。未知キーは拒否）:

| `kind` | 追加キー |
|---|---|
| `mass_rt` | なし（`kind`のみ） |
| `library_match` | `library_id`, `library_sha256`, `score_field`, `score_threshold` |
| `authentic_standard_match` | `mz_tolerance_ppm`（>0）, `rt_tolerance_min`（>0） |

## `matrix_recipes`

`recipe_id`（安全なslug）をキーとするdict。**`"default"`キーが必須**
（欠落は`PROFILE_INVALID`）。各要素は5キーちょうど（spec §8.2）:

| フィールド | 型 | 内容 |
|---|---|---|
| `base` | str | `peak_height` / `internal_standard_ratio` |
| `normalize` | str | `none` / `tic` / `median` / `pqn`。**`base="internal_standard_ratio"`では`none`固定**（二重正規化の禁止。他の値は`PROFILE_INVALID`） |
| `drift_correct` | bool | |
| `filter` | object/null | `null`または`{"min_detection_rate": 0..1}` |
| `impute` | str | `none` / `half_min` / `knn` / `column_mean` |

## `qc_policy`

metricをキーとするdict。metricは次の6つのいずれかのみ（spec §9表）:
`pooled_qc_rsd` / `blank_fold` / `detection_rate` / `standard_rt_error` /
`standard_mass_error_ppm` / `internal_standard_valid_fraction`。

各要素は6キーちょうど:

| フィールド | 型 | 内容 |
|---|---|---|
| `scope` | str | `targets` / `all_features` |
| `target_ids` | list/null | `scope="targets"`なら`feature_targets`参照の非空一意リスト必須、`all_features`なら`null`必須 |
| `required` | bool | バッチQC合否に必須算入するか |
| `threshold` | object | `{"operator": "<="/"<"/">="/">"/"==",  "value": <number>}` |
| `evidence_requirement` | bool | 注入ごとの証拠（§9.1 `assay-feature-evidence.v1`）を要求するか |
| `minimum_pass_fraction` | float/null | 集合判定の閾値（0..1）。単一項目判定なら`null` |

## `evidence`

field path（ドット区切りの識別子列、例 `"acquisition.lc.column"`）をキーと
するdict。各要素は6キーちょうど:

| フィールド | 型 | 内容 |
|---|---|---|
| `value` | any/null | 記録された値。未記載条件は`null` |
| `reason` | str/null | `value`が`null`のとき**必須**（空文字不可）。`value`が非nullなら任意 |
| `tier` | str | `paper_explicit` / `paper_cross_reference` / `supplied` / `locally_validated` / `proposed` |
| `source_uri` | str/null | 出典URI |
| `source_hash` | str/null | ローカル資料のSHA-256 |
| `location` | str/null | 出典中の所在（節・表・ページ等） |

## `validation`

4キーちょうど:

| フィールド | 型 | 内容 |
|---|---|---|
| `status` | str | `draft` / `validated` |
| `scope` | str | 適用範囲の記述。`validated`ではこれと一致する`scope`を持つ証明書が必要 |
| `certificate_path` | str/null | `status="draft"`では`null`必須、`validated`では非空文字列必須 |
| `certificate_sha256` | str/null | `status="draft"`では`null`必須、`validated`では証明書のcanonical hashに一致する64桁小文字16進文字列必須 |

初期の論文由来プロファイルは`draft`。ソフトウェアが動くことを理由に自動昇格しない。

## 証明書 `lcms-profile-validation.v1`

`validate_certificate(profile, certificate, observed_hashes) -> dict`が検証し、
全検証に通った場合のみ正規化済み証明書を返す（失敗時は`None`を返さず
`DomainError`を送出する——「副作用なし」は変わらないが、戻り値は
`dict | 例外`の二択になった）。10キー必須＋`routine_overrides`が任意の
**11キーまで**:

| フィールド | 型 | 内容 |
|---|---|---|
| `schema` | str | 固定値 `"lcms-profile-validation.v1"` |
| `profile_content_sha256` | str | `profile_content_hash(profile)`と一致必須 |
| `dependency_hashes` | object | `{dependency_id: sha256}`。`observed_hashes["dependencies"]`と一致必須 |
| `fixed_input_hashes` | object | `{id: sha256}`。`observed_hashes["fixed_inputs"]`と一致必須 |
| `reference_file_hashes` | object | `{id: sha256}`。`observed_hashes["reference_files"]`と一致必須 |
| `validation_output_hashes` | object | `{id: sha256}`。`observed_hashes["validation_outputs"]`と一致必須 |
| `criteria` | list | 非空。各要素`{criterion_id, status, required, detail}`。`status`は`pass`/`fail`/`not_evaluable`。`required=true`の要素は全て`status="pass"`必須 |
| `performed_by` | str | 実施者 |
| `performed_at` | str | 実施日時 |
| `scope` | str | `profile.validation.scope`と一致必須 |
| `routine_overrides` | object | **任意**（省略可）。`execution_purpose="routine"`が上書きできる範囲の明示宣言。下記「`routine_overrides`」節 |

`observed_hashes`は呼び出し側が実測した参照値（`{"dependencies": {...},
"fixed_inputs": {...}, "reference_files": {...}, "validation_outputs": {...}}`）。
このモジュール自身はファイルにアクセスしないため、実測は必ず呼び出し側が行う。

**検証手順**（`validate_certificate`の内部順序）:

1. `profile.validation.status == "validated"`（draftには適用不可）。
2. `canonical_hash(certificate)` が `profile.validation.certificate_sha256` と
   一致すること——**証明書のどの一要素を改変してもこのhashが変わる**ため、
   実質的な改竄検知はここで完結する。
3. `certificate.profile_content_sha256` が `profile_content_hash(profile)` と一致。
4. `certificate.scope` が `profile.validation.scope` と一致（適用範囲の検証）。
5. 4種のhash辞書が対応する`observed_hashes`のカテゴリと完全一致。
6. `required=true`の`criteria`が全て`status="pass"`。

**`routine_overrides`は上記いずれの手順よりも先に「構造だけ」検証される**
（証明書全体の構造検証`_validate_certificate_shape`の一部）が、その中身が
呼び出し側にとって**意味を持つ**のは、上記6手順すべてを通過し
`validate_certificate`が正常に戻り値を返した後だけである。つまり、
`validation.status='validated'`というラベルだけでは何も証明しない——
証明書自体がpinされたhashと一致しない・profile内容と一致しない等の理由で
手順2〜6のどこかで拒否されれば、`routine_overrides`がどれだけ許容的な
内容であっても、その内容は呼び出し側に一切渡らない
（`test_hand_written_validated_with_mismatched_certificate_rejected_before_allowance_consulted`
がこの順序を固定する）。

署名認証基盤（鍵管理・電子署名）は実装しない。

### `routine_overrides`（spec §6「証明書の明示許容集合」）

`execution_purpose="routine"`のrequestが`preprocess`をどこまで上書きできるかを
証明書自身が明示的に宣言する——**fail-closed**: このキー自体が無ければ
`{}`として扱われ、routineでは何も上書きできない（Task 1が書いた既存の
証明書fixtureはこのキーを持たないが、そのまま有効であり続ける。マイグレー
ション不要）。

```json
"routine_overrides": {
  "preprocess": {
    "<recipe_id>": {
      "<field>": {"mode": "any"}
      // または
      "<field>": {"mode": "values", "values": [<v1>, <v2>, ...]}
    }
  }
}
```

- トップレベルのカテゴリキーは現状`preprocess`だけ（未知のカテゴリは
  `PROFILE_VALIDATION_INVALID`）。将来別カテゴリ（例:
  検定の比較群変更）を追加する余地を残すための開いた構造。
- `<recipe_id>`は`matrix_recipes`のキーである必要はない（証明書は
  profile本体を跨いで再利用されうるため、参照整合性はここでは検査しない
  ——`request_v2`側が実際の`profile.matrix_recipes`と突き合わせて判定する）。
- `<field>`は`normalize` / `drift_correct` / `filter` / `impute`のいずれか
  （未知のfieldは拒否）。
- 許容は`{"mode": "any"}`（どんな値への上書きも許可）または
  `{"mode": "values", "values": [...]}`（列挙した値だけ許可。`values`は
  非空必須）の**どちらか一方**——2つのキー以外は許可しない
  （`{"mode": "any", "values": [...]}`のような混在は拒否）。
- **明示のみ**: この宣言に無いrecipe_id×fieldへの上書きは、profile本体の
  他の場所に同じ値が現れていても許可されない——「profile自身が既に使って
  いる値だから」という推測による許容は行わない（controller裁定。以前の
  `request_v2`独自ルールがこの理由で差し戻された）。
- `execution_purpose="validation"`ではこの宣言を一切参照しない
  （validation目的はまだ証明書が検証していない値を試すためのものだから）。

**`acquisition.scope`は証明書と機械比較されない**。上の手順4が突合するのは
`validation.scope`だけである。`acquisition.scope`（測定法が原理的に対応する
適用範囲）と`validation.scope`（実際にこの証明書が検証した、それより狭いか
同じ範囲）が食い違うこと自体は正当な状態でありうる——例えば測定法が広い
質量範囲をカバーしていても、検証はその一部の化合物クラスだけに限定して
行われた、という場合。spec はこの2つのscopeの一致を要求しておらず、
実行時に値を推測して補うような一致チェックも追加しない。

## 合成profile例（`schema`のみ抜粋、完全体はテスト参照）

```json
{
  "schema": "lcms-profile.v1",
  "profile_id": "kanzo-lcms-metabolomics",
  "revision": 1,
  "omics": "metabolomics",
  "acquisition": {
    "separation": "lc",
    "acquisition_type": "dda",
    "polarity": "positive",
    "instrument": "ExionLC AD / ZenoTOF 7600",
    "lc": {
      "column": "ACQUITY UPLC BEH 1.7um 100x2.1mm",
      "mobile_phase_a": "0.1% formic acid in water",
      "mobile_phase_b": "methanol",
      "flow_rate_ul_min": 150.0,
      "column_temperature_c": 30.0,
      "gradient_profile": null
    },
    "ms_range": {
      "ms1_low_mz": 50.0, "ms1_high_mz": 1500.0,
      "ms2_low_mz": null, "ms2_high_mz": null
    },
    "sample_matrix": "licorice extract",
    "scope": "single LC method, positive polarity, DDA acquisition, peak height only"
  },
  "software": {
    "msdial_version": "5.x",
    "executable_path": "C:/tools/MSDIAL5/MSDIALCUI.exe",
    "executable_sha256": "aaaa...(64桁)",
    "adapter_version": "1.0.0"
  },
  "processing": {
    "method_path": "method/kanzo_pos_param.txt",
    "method_sha256": "bbbb...(64桁)",
    "measure": "peak_height",
    "dependencies": [
      {"dependency_id": "msp-lib-1", "kind": "msp", "method_key": "MSP file",
       "path": "library/lib1.msp", "sha256": "cccc...(64桁)", "required": true}
    ],
    "effective_settings": {"Minimum peak height": "1000"}
  },
  "analysis_recipe": {
    "statistics": [
      {"statistic_id": "pca", "kind": "pca", "matrix_recipe_id": "default",
       "transform": "none", "feature_scope": {"mode": "all_eligible"},
       "scaling": "autoscale", "n_components": 2}
    ],
    "internal_standards": [
      {"target_id": "gaba", "standard_target_id": "gaba_d6"}
    ]
  },
  "feature_targets": {
    "gaba": {
      "compound_identifiers": {"name": "GABA", "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"},
      "adduct": "[M+H]+", "charge": 1,
      "expected_mz": 104.0706, "mz_tolerance_ppm": 10.0,
      "expected_rt_min": 1.2, "rt_tolerance_min": 0.1,
      "required_evidence": [{"kind": "mass_rt"}]
    }
  },
  "matrix_recipes": {
    "default": {"base": "peak_height", "normalize": "none",
                "drift_correct": false, "filter": null, "impute": "none"}
  },
  "qc_policy": {
    "pooled_qc_rsd": {
      "scope": "all_features", "target_ids": null, "required": true,
      "threshold": {"operator": "<=", "value": 30.0},
      "evidence_requirement": false, "minimum_pass_fraction": 0.8
    }
  },
  "evidence": {
    "acquisition.instrument": {
      "value": "ExionLC AD / ZenoTOF 7600", "reason": null,
      "tier": "paper_explicit",
      "source_uri": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12738721/",
      "source_hash": null, "location": "Methods section"
    }
  },
  "validation": {
    "status": "draft",
    "scope": "single LC method, positive polarity, DDA acquisition, peak height only",
    "certificate_path": null, "certificate_sha256": null
  }
}
```

## routine実行が上書きできる範囲

`matrix_recipes` / `feature_targets` / `analysis_recipe` / `qc_policy` はprofileの
固定契約であり、値を変えるには新しいprofile revisionと再検証が要る。
`execution_purpose=routine`のrequestが実行時に指定できるのは、あくまで
「どのrecipe/statisticを選ぶか」（`matrix_recipe_id`・`statistic_id`の選択）や
`preprocess`の明示上書きに限られ、証明書が明示的に許容した範囲
（上記「証明書 `lcms-profile-validation.v1`」節の`routine_overrides`）の外に
出る変更は`PROFILE_SCOPE_MISMATCH`で拒否される（spec §6「routineで許される
比較群の変更等は証明書の適用範囲内に限定する」）。この判定は
`metabolomix.pipeline.request_v2`（`docs/schema/pipeline-request-v2.md`）が
`validate_certificate`の戻り値の`routine_overrides`を消費して行う——
このモジュール自体はMCP requestを扱わない。
