# 定量の規則プラグイン

定量の工程（内部標準による正規化・極性の統合・付加体の選別）の**規則**は、ms-data-parser の外の
パッケージ `msdial_rules` が持つ。ms-data-parser はこれを依存として宣言せず、設定ファイルで
指されたときだけ読み込む。このファイルは、その設定の書式と、プラグインが満たす面（API）を定める。
規則の中身（参照表の値・研究室の標準名）はここに書かない。

## 設定

`lipidmix.local.toml`（雛形 `lipidmix.example.toml`）の `[quant]` 節。**環境変数は無い**。

```toml
[quant]
rules = 'D:\tools\msdial-rules\msdial_rules'   # パッケージ msdial_rules のフォルダ（__init__.py のある所）
adducts = 'D:\tables\adduct-selection.csv'      # 付加体の採用表（任意）
```

- どちらもパスとして解く（相対パスは設定ファイルのあるフォルダ基準）。書き換えは次の呼び出しから効く。
- 内部標準の参照表は設定に既定を置かず、定量ツールの呼び出しごとに指定する。
- `rules` が未設定でも、定量以外のツールは今まで通り動く。

## 読み込み

`rules` のフォルダを固有の名前（`_lipidmix_quant_rules_<指紋の先頭 8 桁>`）で読み込み、
`<パッケージ>.api` を使う。`sys.modules` に `msdial_rules` は残らない。パッケージ内の import は
相対 import にすること（どのフォルダから読まれても動くように）。

指紋（digest）は、フォルダ配下の `*.py`（`__pycache__` を除く）を相対パスの名前順に、
`相対パス NUL 中身 NUL` で連結した sha256。中身が変わると指紋が変わり、次の呼び出しで読み直す。

## API（`API_VERSION = 1`）

引数と戻り値は辞書・list・str・float・None だけ（ms-data-parser は規則側の型を知らない）。

```python
API_VERSION = 1
RULES_VERSION: str          # 規則の版

def normalize(features: list[dict], n_samples: int, reference_text: str,
              method: str = "cell_coefficient", bindings: dict[str, int] | None = None) -> dict
    # features[i]: {feature_id:int, name, lipid_class, mz, rt, adduct, heights:[float|None], detected:[bool]}
    # -> {unit, normalization_method, values:[[float|None]], statuses:[[str]], roles:[str],
    #     assignments:{target_class: {standard_name, feature_id, coefficient, concentration,
    #                  standard_adduct, unit, status, candidates:[{feature_id,name,adduct,mz,rt}]}},
    #     any_others_classes:[str], fallback_by_feature:[bool],
    #     target_class_by_feature:[str|None]}  # 各特徴を割った `assignments` のキー（None = クラス未対応）

def coefficient_reference_text(reference_text: str, concentrations: dict[str, float],
                               solvent_volume: float, cells_millions: float,
                               unit: str | None = None) -> str

def match_samples(names_by_mode: dict[str, list[str]], mapping_rows: list[dict] | None = None) -> dict
    # -> {rule, ids_by_mode:{mode:[共通 ID（元の順）]}, mismatched:{mode:[ID]},
    #     rows_by_mode:{mode:{元の名前: 対応表の行}}|None}

def select(features: list[dict], adduct_text: str | None) -> list[dict]
    # features[i]: {feature_id:str, name, comment, msms_assigned:"True"|"False", ontology, adduct, mode, role}
    # -> [{feature_id, action, reason, original_class, canonical_class, adduct, mode, class_rule_version}]
```

規則が入力を拒むときは `ValueError` の派生を投げる。`.code`（機械可読な種別）と `.details`
（辞書）を持たせると、ms-data-parser が `QUANT_RULES_REJECTED` の封筒の
`details.rule_code` と詳細に写す。

## 記録に残るもの

各工程の記録の `settings.rules` に `{"rules_version", "digest", "api_version"}` を書く。
記録から計算し直すとき、指紋が違えばその工程から先を「古い」として戻さない
（規則の版を固定する代わりに、記録が使った版を縛る）。戻り値と記録に載せる参照表の情報は
sha256・単位・ファイル名の末尾だけで、パスと中身は載せない。

## 封筒

規則のメッセージ・`details` は封筒へ写す前に ms-data-parser が整える（300 字まで・絶対パスはファイル名の末尾へ・`details` は JSON のスカラ値だけ）。**規則側も、メッセージに参照表の行を載せないこと。**

エラーは他の mzTab-M 系ツールと同じ形（`{"error": {"code", "message", "details"?}}`）。

| code | 意味 |
|---|---|
| `quant_rules_unavailable` | `[quant] rules` が未設定、フォルダ・`__init__.py` が無い、または import に失敗した。message に設定の書き方か理由が入る |
| `quant_rules_incompatible` | `API_VERSION` が 1 ではない（`details.api_version` / `supported`） |
| `QUANT_RULES_REJECTED` | 規則が入力を拒んだ（`details.rule_code`） |
| `QUANT_BAD_REQUEST` | 呼び出しの引数が不正（参照表が未指定、`concentrations` に量の指定が欠ける、など） |
| `QUANT_REFERENCE_NOT_FOUND` | 参照表・付加体の採用表のファイルが無い |

## テスト

ms-data-parser のテストは、`tests/quant_fake_rules.py` が書き出す小さな偽のプラグイン
（4 関数を表で決めた単純な規則で実装したもの。本物の写しではない）を `[quant] rules` に向けて使う。
本物の規則との同値確認は `tests/test_quant_parity.py` が行い、pytest の
`--quant-rules=<フォルダ>` を渡したときだけ走る（既定では skip）。
