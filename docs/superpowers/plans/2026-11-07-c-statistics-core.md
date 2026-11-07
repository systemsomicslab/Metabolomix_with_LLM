# C 統計エンジンの中核（配線なし）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** spec C の §4〜§12（前処理・0 と欠測・moderated t / F・効果量と判定・多重性・交絡・PCA・組成・v1 互換）を、
ツールに配線しない純粋な関数とテストとして `metabolomix/analysis/stats/` に作る。

**Architecture:** `DatasetState` を `StudyData`（値・`detected_mask`・試料の表・定量の状態）へ写す adapter を入口にし、
前処理（§4）が名前付きの `PreparedMatrix`（log2・個体単位）を作り、比較のエンジン `compare`（§5〜§9）・PCA（§10）・
組成（§11）がそれを読む。analysis-spec.v2（新しい版）の正規化と実行もここに置く。analysis-spec.v1 の実行系
（`analysis_spec.py` / `spec_statistics.py` / `recipe_matrix.py`）は版の意味を固定したまま触らず、黄金値のテストで縛る（§12）。

**Tech Stack:** Python 3.14、numpy、scipy（`special.digamma` / `polygamma`、`stats`）、statsmodels（LOWESS）、scikit-learn（PCA）。
依存は足さない。moderated t の参照値は R の limma で一度だけ作る。

**Spec:** `docs/superpowers/specs/2026-11-06-statistics-engine-design.md`（全体は `2026-11-06-overhaul-overview-design.md`）

## Global Constraints

- 触らない: `metabolomix/tools/`・`server.py`・`USAGE.md`・`metabolomix/arf/`・`metabolomix/analysis/analysis_spec.py`・
  `spec_statistics.py`・`recipe_matrix.py`・`differential.py`・`statistics_v2.py`・`preprocessing.py`（B-2 と衝突させない。読んで再利用は可）。
- 既存の `metabolomix/analysis/*.py` を変えるのは Task 8 の `pca.py` への引数 1 つの追加だけ（後方互換。HISTRY に書く）。
- `metabolomix/analysis/stats/` は MCP にもセッションにも依存しない（`metabolomix.core.session_state`・`mcp_core`・`metabolomix.tools.*` を import しない）。
- 結果に時刻・乱数の ID を入れない。同じ入力と設定なら同じ `matrix_id`・同じ数字（決定論）。
- 0 以下・非有限は欠測。pseudocount も clip も使わない。補完した値で検定しない。
- 既定（spec の値をそのまま）: q の閾値 0.05、同等の幅 ±20%（|log2FC| < log2(1.2) ≈ 0.263）、区間 95%、QC RSD 0.30、
  D-ratio 0.5、量の検定に要る実検出の割合 0.5、各群 2 個体未満の特徴は検定しない。
- n<2 の群は黙って落とさず `STATISTIC_GROUP_TOO_SMALL` で止める（v2 の既存の規則）。
- テストは main ツリーの `.venv`: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。
  コミットは pre-commit が全テストを回すので、バックグラウンドで実行し出力をファイルへ落とす。`git stash` は使わない。
- エラーは `metabolomix.core.atomic_io.DomainError(code, message, details)`。新しい code は各タスクに書いたものだけ。
- 文書・docstring は日本語（リポジトリの流儀）。

## Review Focus

1. **`detected_mask` が無いデータセット**（mzTab-M だけ）: 出方の検定は走らず `detection_available=false`、§5 の印は
   付けない（付けられない理由を書く）、量の検定は走る。落ちない。→ Task 7 `test_no_detected_mask_runs_quantity_only`
2. **片方の群で全欠測の特徴**（出方だけが違う）: 量は `not_tested`、出方の Fisher と q は出る、判定は `inconclusive`。
   → Task 7 `test_feature_absent_in_one_group_is_presence_only`
3. **1 匹 3 注入 × 各群 1 匹**: 個体でまとめると各群 n=1 → `STATISTIC_GROUP_TOO_SMALL`（n=3 と数えない）。
   → Task 7 `test_replicates_collapse_before_group_size_check`
4. **QC が 1 本も無い**: ドリフト補正・QC RSD・D-ratio は理由つきで飛ばし、PQN の参照は試料の中央値。
   → Task 6 `test_no_qc_skips_qc_steps_with_reasons`
5. **分散 0 の特徴・検定できる特徴が 1 つだけ**: 事前分布の推定（`fit_f_dist`）が NaN・inf・例外を出さない
   （limma と同じく `max(x, 1e-5·median)`、n=1 は `df_prior=0`）。→ Task 2 `test_fit_f_dist_degenerate_inputs`

---

## §15「実装者に任せること」の判断

| 項目 | 判断 | 理由 |
|---|---|---|
| エンジンの形 | `metabolomix/analysis/stats/` パッケージ。比較の入口は `stats.compare.compare` の 1 つ（moderated t・welch・student・mannwhitney・F をすべてここから）。 | spec §15「入口は 1 つ」。 |
| v1 の扱い | analysis-spec.v1 は今の実行系（`analysis_spec.run_spec`）を**版の固定した意味**として残し、黄金値で縛る。v1 を新エンジンで再実装しない。対話の v1（`differential.two_group_test`・`statistics_v2` の welch・`pathway/compare.py` の Welch）を消してエンジンへ寄せるのは配線段。 | §12「版ごとの意味を固定」。再実装は同じ数字を保証しにくい。 |
| welch / student / mannwhitney | エンジンの選択肢として `spec_statistics.two_group(..., transform="none")` を log2 済みの値に呼ぶ（検定の実装を 1 つに保つ）。 | 重複を作らない。 |
| LOESS の幅 | 交差検証はしない。既定 `span=0.75`（matrix-recipe.v2 と同じ）、局所線形、`it=0`、spec で変えられる。バッチごとに当てる。 | バッチあたりの QC は 5〜20 本で、交差検証の幅の選択が不安定。決定論も保てる。 |
| ドリフト補正の可否 | 全バッチが「注入順あり・有限正値の QC が `min_qc_per_batch`（既定 4）以上・QC が試料の注入順を挟む」を満たすときだけ行う。1 つでも欠ければ**全体を**理由つきで飛ばす。 | 補正したバッチとしないバッチを混ぜない。 |
| バッチ間の水準合わせ | 各バッチを自分の QC 中央値へ補正したあと、特徴ごとに `全 QC の中央値 / そのバッチの QC 中央値` を掛ける。 | spec §4 の 2。 |
| 検出率の足切りの既定 | `by="group"`、`min_rate=0.5`: どれか 1 つの群で生体試料の半数以上が実検出なら残す。`detected_mask` が無ければ飛ばし理由を書く。 | 出方だけが違う特徴（§5）を足切りで消さない。 |
| ブランク比の既定 | `min_ratio=3.0`（`preprocessing.blank_filter` の既定と同じ）。ブランクが無ければ飛ばす。 | 既存と同じ値。 |
| 正規化の既定 | `auto`: 内部標準で換算済み（`QuantState.normalized`）なら `none`、そうでなければ PQN（参照は QC の中央値、QC が無ければ生体試料の中央値）。 | spec §4 の 4。 |
| 反復のまとめ | 同じ `subject_id` の注入を log2 の上で平均（有限値だけ）。個体の「実検出」は 1 注入でも実検出なら真。`subject_id` が無い試料は各注入が 1 個体（`subject_source="injection"` を結果に書く）。 | spec §4 の 6・§6。 |
| 反復注入らしい試料名 | 名前が `^(.*\d)[_\-.](?:inj|rep|tech|r)?\d{1,2}$`（茎が数字で終わる 2 段の番号）か `(?i)[_\-.](inj|rep|tech)\d+$` に合い、同じ茎が 2 つ以上あるとき注意に出す。まとめない。 | `HC_1`・`HC_2` のような普通の個体名を誤検出しない。 |
| §5 の印の数え方 | 対比の 2 群それぞれで「実検出の個体数 / 群の個体数 < `min_detected_fraction`」なら `quantity_judgeable=false`。 | spec §5。 |
| 検定しない特徴と事前分布 | 各群 2 個体未満の特徴は当てはめから外し、事前分布（`s2_prior`・`df_prior`）の推定にも入れない。 | 判定できない特徴の分散で事前分布を動かさない。R の参照 fixture はこの種の特徴を含めない。 |
| 判定の優先 | ①検定していない・§5 の印 → `inconclusive`、②q ≤ 閾値 → `changed`、③区間が ±幅の内 → `unchanged`、④それ以外 → `inconclusive`。区間が無い（算術平均比・mannwhitney）ときは `unchanged` にしない。 | spec §7。 |
| 出方の BH | 対比ごとに、Fisher の p を持つ特徴（`detected_mask` あり）全体で BH。量の BH とは別。 | spec §5・§8。 |
| 部分的な交絡 | 生体試料の群 × バッチの表で、バッチが 2 つ以上あり (a) 1 群しか含まないバッチがある、または (b) Cramér の V ≥ 0.5 なら注意。完全交絡（どのバッチも 1 群・どの群も 1 バッチ）は `allow_confounded` なしで `STATISTIC_CONFOUNDED`。 | spec §9。 |
| 注入順の偏り | 注入順が全生体試料にあり、各群 2 本以上なら Kruskal–Wallis。p < 0.05 で注意（群ごとの注入順の範囲を添える）。 | spec §9（私の判断の項）。 |
| PCA の補完 | 既定 `half_min`（log2 の上では特徴ごとの生体個体の最小値 − 1）。`knn` は `recipe_matrix.knn_impute`。補完したセル数を結果に書く。 | spec §4 の 7・§10。 |
| CLR の 0 | クラスの合計は、その個体で欠測の分子種を 0 として足す。合計が 0 以下の個体が 1 つでもあるクラスは CLR から**除き**、除いたクラスと理由を書く（置き換えはしない）。 | 補完した値で検定しない（反目標）。 |
| 組成の個体単位 | `PreparedMatrix` の log2（個体単位）を `2**x` に戻して足す（個体内の注入は幾何平均でまとめた後）。 | 前処理を 1 本にする。 |
| `matrix_id` | `"mx_" + canonical_hash({dataset_hash, 試料の表の hash, 正規化済みの前処理設定})[:16]`。 | 決定論。 |
| spec v2 の欄の名前 | Task 10 の表のとおり。 | spec §15。 |

---

## ファイル構成

```
metabolomix/analysis/stats/
  __init__.py        公開名の __all__ だけ（中身は各モジュール）
  data.py            StudyData・SampleInfo・QuantState・from_dataset・quant_state_of・cell_states・replicate_name_warnings（§5・§6）
  limma.py           lm_fit・contrasts_fit・fit_f_dist・trigamma_inverse・ebayes（§6 の moderated t / F と §7 の区間）
  design.py          Design・build_design（群・2 因子・共変量・ブロック）と対比（§6）
  confounding.py     confounding_report（§9）
  preprocess.py      PreparedMatrix・NeedsInput・normalize_preprocess_config・preprocess・impute_for_pca（§4）
  compare.py         compare（§5〜§9 をまとめる比較のエンジン。唯一の入口）
  pca.py             pca（§10）
  composition.py     composition（§11）
  spec_v2.py         analysis-spec.v2 の normalize_spec_v2・run_spec_v2（§8 の primary・§12 の新しい版）
metabolomix/analysis/pca.py   run_pca_fit_subset に scaling="autoscale"|"pareto" を足す（Task 8 だけ）
tests/stats_reference/
  analysis_spec_v1_golden.json     Task 1 の黄金値
  limma_input.tsv                  Task 2 の入力（特徴 × 試料）
  limma_design.tsv                 Task 2 の試料の表
  make_limma_expected.R            期待値の作り方（一度だけ実行）
  limma_expected.json              R の出力
tests/test_stats_v1_compat.py  test_stats_limma.py  test_stats_design.py  test_stats_data.py
tests/test_stats_confounding.py  test_stats_preprocess.py  test_stats_compare.py  test_stats_pca.py
tests/test_stats_composition.py  test_stats_spec_v2.py
```

---

### Task 1: analysis-spec.v1 の黄金値で数字を固定する（§12）

**Files:**
- Create: `tests/test_stats_v1_compat.py`, `tests/stats_reference/analysis_spec_v1_golden.json`

**Interfaces:**
- Consumes: `analysis_spec.normalize_spec`・`run_spec`（変えない）
- Produces: 配線段が v1 の実行経路を差し替えても数字が動かないことを縛る黄金値

- [ ] **Step 1: 合成入力を作る関数と黄金値の比較テストを書く**

テスト自身が `tmp_path` に入力を作る（`tests/test_analysis_spec.py` の `_write` の流儀）。2 つの入力を使う:
(a) 外部 TSV（試料 12: 群 HC×4・RA×4・QC×4、特徴 8。0・空欄・gap 的な小さな値を含む。注入順と共変量 `sex`/`age` の TSV）、
(b) `SimpleNamespace` の疑似データセット（`feature_matrix`・`feature_ids`・`sample_names`・`sample_ids`・
`sample_metadata_rows`（role/group/include を含む）・`detected_mask`・`feature_candidates`（一部の特徴だけ）・
`feature_annotations={}`・`quantification_measure`・`source_files`、`merge_split=false`）。
解析は次をすべて含む: `student`・`welch`・`mannwhitney`・`ols`（共変量 2 つ、欠測 1）、`fold_change` の 2 定義、
`bh_population` の `matrix`・`reported`（`features.feature_ids`）・`annotated`（(b) だけ）、recipe の 5 段すべて
（`detection_filter`・`impute`（half_min と knn）・`normalize`（pqn/tic/median）・`drift_correct`・`qc_rsd_filter`）、`aggregate` の sum と ratio。

```python
CASES = {"tsv_student_recipe": ..., "tsv_ols": ..., "dataset_annotated": ..., ...}  # 名前 → spec を作る関数

@pytest.mark.parametrize("case", sorted(CASES))
def test_analysis_spec_v1_numbers_are_frozen(case, tmp_path):
    result = _run_case(case, tmp_path)
    expected = json.loads(GOLDEN.read_text(encoding="utf-8"))[case]
    _assert_same(_numbers(result), expected)   # 数は math.isclose(rel_tol=1e-9, abs_tol=1e-12)、それ以外は ==

def test_golden_covers_every_v1_option():
    # KINDS・FC_DEFINITIONS・bh_population の 3 種・recipe の STEPS がどこかのケースで使われている
```

`_numbers(result)` は `rule_hash`・`analyses[*].rows` の全欄・`matrix.steps`・`bh_population_n` を取り出す（パスを含む欄は除く）。

- [ ] **Step 2: 黄金値を作る**

`UPDATE_V1_GOLDEN=1` のときだけ `_numbers` を JSON に書き出す分岐をテストに置き、一度実行して JSON を作る。
Run: `UPDATE_V1_GOLDEN=1 <venv python> -m pytest tests/test_stats_v1_compat.py -q` → JSON ができる。

- [ ] **Step 3: 通常実行で通ることを確かめる**

Run: `<venv python> -m pytest tests/test_stats_v1_compat.py -q` → PASS。試しに `spec_statistics.bh` の 1 行を壊すと FAIL することを手元で確かめ、戻す。

- [ ] **Step 4: Commit** — `test(stats): analysis-spec.v1 の数字を黄金値で固定する（C §12）`

---

### Task 2: limma 型の moderated t / F（§6・§7 の区間）と R の参照値

**Files:**
- Create: `metabolomix/analysis/stats/__init__.py`, `metabolomix/analysis/stats/limma.py`, `tests/test_stats_limma.py`,
  `tests/stats_reference/{limma_input.tsv,limma_design.tsv,make_limma_expected.R,limma_expected.json}`

**Interfaces:**
- Produces:
  - `@dataclass LinearFit: coefficients (G×p), stdev_unscaled (G×p), sigma (G,), df_residual (G,), cov_coefficients (p×p), column_names: list[str]`（行は特徴。NaN は推定できない係数）
  - `lm_fit(y: np.ndarray (G×n), design: np.ndarray (n×p), column_names: list[str]) -> LinearFit` — 特徴ごとに NaN の列（試料）を除いて最小二乗（limma `lmFit` の欠測の扱い）。`cov_coefficients` は欠測なしの計画行列から。
  - `contrasts_fit(fit: LinearFit, contrasts: np.ndarray (p×k), names: list[str]) -> LinearFit` — limma `contrasts.fit` と同じ（特徴ごとの `stdev_unscaled` は欠測を除いた計画から計算し直す）。
  - `trigamma_inverse(x: np.ndarray) -> np.ndarray` — limma `trigammaInverse`（Newton 法）。
  - `fit_f_dist(s2: np.ndarray, df1: np.ndarray) -> tuple[float, float]` — `(s2_prior, df_prior)`。limma `fitFDist`（`covariate=NULL`）と同じ: `x=max(x, 1e-5·median)`、`e=log x − ψ(df1/2) + log(df1/2)`、`evar = Σ(e−ē)²/(n−1) − mean ψ'(df1/2)`、`evar>0` なら `df_prior = 2·trigamma_inverse(evar)`、さもなくば `inf`。n=1 は `(s2[0], 0.0)`。
  - `@dataclass Moderated: t (G×k), p (G×k), coefficients (G×k), ci_low (G×k), ci_high (G×k), df_total (G,), s2_post (G,), s2_prior: float, df_prior: float, f (G,), f_p (G,), contrast_names: list[str]`
  - `ebayes(fit: LinearFit, *, ci_level: float = 0.95) -> Moderated` — `trend=False, robust=False`。`df_total = min(df_residual + df_prior, Σdf_residual)`。区間は `coef ± t_{df_total, (1+level)/2} · √s2_post · stdev_unscaled`（limma `topTable(confint=)`）。F は limma `classifyTestsF(fstat.only=TRUE)`（`cov_coefficients` から作る対比の相関で直交化）、`f_p` は `F(k, df_total)` の上側。

- [ ] **Step 1: R の参照値を作る**（一度だけ。下の「R の参照値について」を参照）

`limma_input.tsv`: 特徴 30 × 試料 12（3 群 A/B/C × 4。log2 の値。乱数の種を固定して作り、TSV に固定値で保存。
NaN を 6 特徴に 1〜2 個（各群 2 個体以上は残す）、分散がほぼ 0 の特徴 1 つ）。`limma_design.tsv`: `sample, group, batch, subject`
（`subject` は B と C が同じ 4 個体＝対応のある設計）。`make_limma_expected.R` は 3 つの設計を当てて JSON に書く:
1. 群の cell means（`~0+group`）、対比 B−A・C−A。`topTable(coef=, confint=0.95, number=Inf, sort.by="none")` の
   `logFC, CI.L, CI.R, t, P.Value` と `fit$df.total, fit$s2.prior, fit$df.prior, fit$F, fit$F.p.value`。
2. 群 + バッチの共変量（`~0+group+batch`）、対比 B−A。
3. 対応のある設計（B と C の 8 試料、`~0+group+subject`）、対比 C−B。

R のバージョン・limma のバージョンを JSON の `_meta` に書く。スクリプトの冒頭に実行方法を書く。

- [ ] **Step 2: 失敗するテストを書く**

```python
@pytest.mark.parametrize("design_id", ["cell_means", "batch_covariate", "paired_block"])
def test_moderated_t_matches_limma(design_id):
    exp = EXPECTED[design_id]
    mod = ebayes(contrasts_fit(lm_fit(Y, X, cols), C, names))
    assert mod.s2_prior == approx(exp["s2_prior"], rel=1e-6); assert mod.df_prior == approx(exp["df_prior"], rel=1e-6)
    for each contrast: t, p, coefficients, ci_low, ci_high, df_total → approx(rel=1e-6, abs=1e-10)
    if design_id == "cell_means": f, f_p → approx(rel=1e-6)

def test_fit_f_dist_degenerate_inputs():
    assert fit_f_dist(np.array([0.3]), np.array([4.0])) == (0.3, 0.0)
    s2p, d0 = fit_f_dist(np.zeros(5), np.full(5, 4.0)); assert math.isfinite(s2p) and s2p > 0
    s2p, d0 = fit_f_dist(np.full(5, 0.2), np.full(5, 4.0)); assert s2p == approx(0.2) and d0 == math.inf
    mod = ebayes(lm_fit(...全特徴が同じ分散...)); assert np.all(np.isfinite(mod.p))   # df_prior=inf でも

def test_trigamma_inverse_roundtrip():
    x = np.array([0.05, 0.5, 2.0, 50.0]); assert polygamma(1, trigamma_inverse(x)) == approx(x, rel=1e-8)
```

- [ ] **Step 3: `limma.py` を実装して通す** — Run: `pytest tests/test_stats_limma.py -q` → PASS
- [ ] **Step 4: Commit** — `feat(stats): limma 型の moderated t / F と区間（R の limma と突き合わせ）`

**R の参照値について:** 手元の R 4.6.0 に limma が無い（`r/` の renv にも無い）。利用者の許可を得て、Bioconductor から
`BiocManager` と `limma`（と `jsonlite`）を**一時のライブラリ**（スクラッチ。`r/` の renv にも R 本体のライブラリにも入れない）へ入れ、
`make_limma_expected.R` を一度だけ実行して `limma_expected.json` を作る。許可が無いときは、limma の公開文書の数値例は使わず
`statsmodels` の OLS（特徴ごとの係数・SE・残差の自由度）と、Smyth (2004) の式を手で書いた別実装（テスト内）の 2 つと
突き合わせ、HISTRY と JSON の `_meta` に「R とは突き合わせていない」と明記する。

---

### Task 3: 計画行列と対比（§6）

**Files:**
- Create: `metabolomix/analysis/stats/design.py`, `tests/test_stats_design.py`

**Interfaces:**
- Consumes: なし（`units: list[UnitInfo]` の欄 `group`・`factors`・`subject_id`・`batch` と共変量の値だけを読む。`UnitInfo` は Task 6）
- Produces:
  - `@dataclass Design: matrix (n×p), columns: list[str], cells: list[str], contrasts: dict[str, np.ndarray (p,)], contrast_info: dict[str, dict], omnibus: np.ndarray (p×m)`
  - `build_design(units, *, factors: list[str], comparisons: list[dict], covariates: list[dict] = (), block: str | None = None) -> Design`
    - 因子の組み合わせごとの cell means（切片なし）。1 因子なら群がそのまま cell。共変量は `{"name", "type": "numeric"|"categorical", "reference"}`（v1 の共変量の定義と同じ形）、カテゴリは基準以外の水準のダミー。`block` は個体のダミー（最初の個体を基準）。
    - `comparisons` の形: `{"id", "factor", "reference", "test"}`（2 因子のときは `"at": {他の因子: 水準}` を必須にした単純効果）か `{"id", "kind": "interaction", "factors": [A, B], "levels": {A: [ref, test], B: [ref, test]}}`。
    - `omnibus` は「どれかの cell が違うか」の対比（各 cell − 最初の cell）。moderated F に使う。
    - 計画行列が階数不足（個体が群の中に入れ子のブロックなど）なら `DomainError("DESIGN_RANK_DEFICIENT", ..., {"columns": ...})`。存在しない水準は `DomainError("DESIGN_INVALID")`。

- [ ] **Step 1: テストを書く**

```python
def test_two_groups_contrast_is_test_minus_reference(): Design の contrasts["RA_vs_HC"] が cell RA に +1・HC に −1
def test_two_by_two_simple_effect_and_interaction(): interaction の対比が (a2b2 − a1b2) − (a2b1 − a1b1) の係数を持つ
def test_simple_effect_requires_at_for_two_factors(): "at" 無し → DESIGN_INVALID
def test_block_nested_in_group_is_rank_deficient(): 各個体が 1 群だけ → DESIGN_RANK_DEFICIENT
def test_paired_block_is_full_rank(): 同じ個体が両群 → 列 = cell 2 + 個体 (k−1)
def test_categorical_covariate_dummies_and_reference(): batch の基準以外の水準だけ列になる
def test_omnibus_has_cells_minus_one_columns()
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_design.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): 計画行列と対比（群・2 因子・共変量・個体のブロック）`

---

### Task 4: DatasetState からの入口と 0・gap-fill・欠測の 3 状態（§5・§6）

**Files:**
- Create: `metabolomix/analysis/stats/data.py`, `tests/test_stats_data.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) SampleInfo: sample_id: str, role: str, group: str | None, batch: str | None, injection_order: float | None, subject_id: str | None, include: bool, attributes: dict[str, str]`（`attributes` は因子・共変量に使う追加の欄）
  - `@dataclass(frozen=True) QuantState: normalized: bool, method: str | None, unit: str | None`
  - `@dataclass StudyData: values (n_samples×n_features, 生の値), detected (同形の bool | None), sample_ids, feature_ids, samples: list[SampleInfo], quant: QuantState | None, dataset_hash: str | None, feature_classes: dict[str, str], annotated_features: frozenset[str]`
  - `from_dataset(ds) -> StudyData` — 値は `ds.feature_matrix.T`、`detected` は `ds.detected_mask.T`（None はそのまま）。試料の表の出どころの優先: `ds.sample_metadata_rows`（`subject_id`、無ければ `biological_sample_id`、ほかの未知の欄は `attributes`）→ `ds.msdial_samples`（group・role・batch・injection_order）→ `ds.roles` と `ds.sample_meta`（`run_order`）。`include=false` の試料も残す（選択は前処理の 1 段目）。`dataset_hash` は `result_state.dataset_fingerprint(ds)`。`feature_classes` は `feature_metadata[fid]["ontology"]`（空は除く）、`annotated_features` は `feature_candidates` が空でない特徴。
  - `quant_state_of(ds) -> QuantState | None` — **B-3 との受け口（2026-11-07 に B-3 と合意）**。B-3 は `DatasetState.__init__` に `self.quant_state: dict | None = None` を足し、
    `{"schema": "lipidmix.quant-state.v1", "normalized": bool, "method": "internal_standard" | None, "unit": str | None, "record_id", "steps", "skipped"?}` を置く。
    ここは `normalized`・`method`・`unit` だけを読み、ほかの鍵は無視する。属性が無い・None なら None。`normalized=False` は換算していない扱い。
    B-3 の正規化後の `detected_mask` は「対象と内部標準の両方が検出」のセルだけ真（gap-fill 由来は偽）。
  - `DETECTED, GAP_FILLED, MISSING, PRESENT_UNKNOWN = 1, 2, 0, 3`
  - `cell_states(values, detected) -> np.ndarray[int8]` — 有限正値かつ `detected` 真 → DETECTED、有限正値かつ偽 → GAP_FILLED、それ以外 → MISSING、`detected` が None なら有限正値は PRESENT_UNKNOWN。
  - `replicate_name_warnings(sample_ids: list[str]) -> list[dict]` — `{"code": "possible_technical_replicates", "stem": ..., "samples": [...]}`（規則は §15 の表）。

- [ ] **Step 1: テストを書く**

```python
def test_from_dataset_reads_metadata_rows_and_subject_id(): subject_id 優先、biological_sample_id へ落ちる、未知の欄は attributes
def test_from_dataset_falls_back_to_msdial_samples_then_roles()
def test_quant_state_port_absent_and_present(): 属性が無ければ None、dict なら QuantState
def test_cell_states_three_way_and_unknown(): 0・負・NaN は MISSING、gap-fill は GAP_FILLED、mask None は PRESENT_UNKNOWN
def test_replicate_name_warnings_flags_two_level_numbers_only():
    assert replicate_name_warnings(["HC_1", "HC_2", "RA_1"]) == []
    w = replicate_name_warnings(["M01_1", "M01_2", "M02_1", "S3_inj1", "S3_inj2"]); stems == {"M01", "S3"}
def test_from_dataset_does_not_mutate_dataset(): feature_matrix・detected_mask の配列が書き換わらない
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_data.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): DatasetState からの入口と検出・gap-fill・欠測の 3 状態`

---

### Task 5: 交絡（§9）

**Files:**
- Create: `metabolomix/analysis/stats/confounding.py`, `tests/test_stats_confounding.py`

**Interfaces:**
- Consumes: `UnitInfo`（Task 6。欄 `group`・`batch`・`injection_order` だけ。テストでは `SimpleNamespace` で足りる）
- Produces: `confounding_report(units, groups: list[str]) -> dict`:
  `{"batch_table": {group: {batch: n}}, "batch_status": "not_evaluable"|"single_batch"|"none"|"partial"|"complete", "cramers_v": float|None, "injection_order": {"status": "not_evaluable"|"ok"|"biased", "kruskal_p": float|None, "range_by_group": {g: [min, max]}}, "warnings": [{"code", "message"}]}`。
  バッチが 1 つでも欠ける生体個体があれば `not_evaluable`（「交絡なし」と書かない）。判定の規則は §15 の表。
  `assert_not_confounded(report, *, allow_confounded: bool) -> None` — `complete` で `allow_confounded` が偽なら `DomainError("STATISTIC_CONFOUNDED", ..., {"batch_table": ...})`。

- [ ] **Step 1: テストを書く**

```python
def test_complete_confounding_blocks_unless_allowed()
def test_partial_bias_single_group_batch_warns(): batch_status == "partial" と warning
def test_cramers_v_threshold(): 偏りが強い表 → partial、均等 → none、V の値を確かめる（手計算の値）
def test_missing_batch_is_not_evaluable_not_none()
def test_injection_order_grouped_runs_warn(): 群ごとに固めて測った注入順 → biased、交互 → ok
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_confounding.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): 群 × バッチの表・完全と部分の交絡・注入順の偏り`

---

### Task 6: 前処理（§4）

**Files:**
- Create: `metabolomix/analysis/stats/preprocess.py`, `tests/test_stats_preprocess.py`

**Interfaces:**
- Consumes: `StudyData`・`SampleInfo`・`QuantState`・`cell_states`（Task 4）
- Produces:
  - `PREPROCESS_DEFAULTS = {"drift": {"mode": "auto", "span": 0.75, "min_qc_per_batch": 4}, "filters": {"blank_ratio": {"min": 3.0}, "detection": {"by": "group", "min_rate": 0.5}, "qc_rsd": {"max": 0.30}, "d_ratio": {"max": 0.5}}, "normalize": {"method": "auto"}, "collapse_replicates": True}`
    （`drift.mode`: `auto`|`on`|`off`。各足切りは `null` で無効。`normalize.method`: `auto`|`none`|`tic`|`median`|`pqn`）
  - `normalize_preprocess_config(raw: dict | None) -> dict` — 既定を埋め、未知のキーは `DomainError("PREPROCESS_CONFIG_INVALID")`。
  - `@dataclass UnitInfo: unit_id: str, group: str | None, batch: str | None, subject_id: str | None, injection_order: float | None, n_injections: int, sample_ids: list[str], attributes: dict[str, str]`
  - `@dataclass PreparedMatrix: matrix_id: str, log2 (n_units×n_features, NaN=欠測), units: list[UnitInfo], feature_ids: list[str], detected_units (bool | None), qc_log2 (n_qc×n_features), qc_ids: list[str], steps: list[dict], config: dict, quant: QuantState | None, dataset_hash: str | None, subject_source: str, warnings: list[dict], feature_classes: dict, annotated_features: frozenset`
  - `@dataclass NeedsInput: field: str, reason: str, options: list[dict]`（`to_envelope()` → `{"status": "needs_input", "field", "reason", "options"}`）
  - `preprocess(study: StudyData, config: dict) -> PreparedMatrix | NeedsInput`
  - `impute_for_pca(log2: np.ndarray, *, method: str = "half_min", reference_rows: np.ndarray | None = None) -> tuple[np.ndarray, dict]` — `half_min` は `reference_rows`（既定は全行）の特徴ごとの最小 − 1、`knn` は `recipe_matrix.knn_impute(k=10, rowmax=0.5, colmax=0.8)` を log2 の上で。報告に `{"method", "imputed_cells"}`。

`preprocess` の段（`steps` に順に記録。各段 `{"step", "status": "applied"|"skipped", "reason"?, ...}`）:
1. 試料の選択: `include=false` を外す。役割は `sample`（生体）・`qc`・`blank`。それ以外の役割は外して記録。
2. ドリフト補正: `mode=auto` かつ `quant.normalized` かつ QC あり → `NeedsInput(field="preprocess.drift.mode", reason=…, options=[{"value": "on", ...}, {"value": "off", ...}])` を返す。条件（§15 の表）を満たせばバッチごとの LOESS（statsmodels `lowess(frac=span, it=0, delta=0, xvals=注入順)`）→ バッチの水準合わせ。満たさなければ `skipped` と理由（`no_injection_order`・`too_few_qc:<batch>`・`qc_not_interspersed:<batch>`）。`mode=on` で条件を満たさなければ `DomainError("PREPROCESS_PRECONDITION")`。
3. 品質の足切り: ブランク比・検出率（`detected` が無ければ skipped）・QC RSD（QC<2 なら skipped）・D-ratio（QC の SD / 生体試料の SD、どちらか <2 なら skipped）を**ドリフト後の値で同時に**計算し、論理積で一度に当てる。落とした特徴を足切りごとに記録。
4. 正規化: 係数は残った特徴だけから。`auto` の解決（§15 の表）と使った参照を記録。係数が非有限・0 以下の試料は 1 にして記録（`recipe_matrix._normalize` と同じ扱い）。
5. 変換: 有限正値だけ log2。欠測にしたセルを理由別に数える（`non_positive`・`non_finite`）。
6. 反復のまとめ: §15 の規則。`subject_source` と、同じ個体の注入が別の群・バッチを持つときは `DomainError("SUBJECT_INCONSISTENT")`。
7. 補完はしない（`impute_for_pca` は別関数）。
`matrix_id` は §15 の表のとおり。`warnings` に `replicate_name_warnings` を入れる（`subject_id` が 1 つも無いときだけ）。

- [ ] **Step 1: テストを書く**

```python
def test_step_order_is_select_drift_filter_normalize_log_collapse(): steps の step 名の並び
def test_normalization_factors_use_only_features_kept_by_filters(): 足切りで落ちる外れ特徴を混ぜても係数が変わらない
def test_drift_is_corrected_per_batch_then_levels_aligned(): バッチ 2 つに逆向きの線形ドリフト → 補正後の QC の傾きが ~0、バッチの QC 中央値が揃う
def test_drift_skipped_when_qc_not_interspersed_in_a_batch(): 理由 qc_not_interspersed:<batch>
def test_is_normalized_with_qc_and_auto_drift_returns_needs_input()
def test_is_normalized_defaults_to_no_normalization(); test_default_normalization_is_pqn_with_qc_reference()
def test_no_qc_skips_qc_steps_with_reasons()          # Review Focus 4: PQN の参照は sample
def test_filters_masks_computed_after_drift_and_applied_together()
def test_d_ratio_filter_default_half()
def test_detection_filter_by_group_keeps_presence_only_features()
def test_log2_marks_non_positive_missing_without_pseudocount(): 0 → NaN、理由の数
def test_replicates_averaged_on_log_scale_by_subject(): 2 注入 [4, 16] → log2 平均 3、detected_units は「どれか」
def test_matrix_id_is_deterministic_and_config_sensitive()
def test_impute_for_pca_half_min_and_report()
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_preprocess.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): 前処理（ドリフトのバッチ別 LOESS・同時の足切りと D-ratio・残った特徴での正規化・個体でのまとめ）`

---

### Task 7: 比較のエンジン `compare`（§5〜§9）

**Files:**
- Create: `metabolomix/analysis/stats/compare.py`, `tests/test_stats_compare.py`

**Interfaces:**
- Consumes: `PreparedMatrix`・`UnitInfo`（Task 6）、`build_design`（Task 3）、`lm_fit`・`contrasts_fit`・`ebayes`（Task 2）、
  `confounding_report`・`assert_not_confounded`（Task 5）、`spec_statistics.two_group`・`bh`・`detection_fisher_p`（既存。変えない）
- Produces:
  - `TESTS = ("moderated_t", "welch", "student", "mannwhitney")`、`EFFECT_SIZES = ("log2_mean_difference", "arithmetic_mean_ratio")`、`BH_POPULATIONS = ("matrix", "annotated", "reported")`
  - `compare(prepared: PreparedMatrix, *, factors: list[str], comparisons: list[dict], covariates: list[dict] = (), block: str | None = None, test: str = "moderated_t", effect_size: str = "log2_mean_difference", bh_population: str = "matrix", reported_feature_ids: list[str] | None = None, q_threshold: float = 0.05, equivalence_fc: float = 1.2, min_detected_fraction: float = 0.5, allow_confounded: bool = False) -> dict`
  - 戻り値（丸めない）: `{"schema": "stats-compare.v1", "matrix_id", "test", "effect_size_definition", "ci_level": 0.95, "verdict_rule": {"q_threshold", "equivalence_fc", "equivalence_log2"}, "bh_population": {"name", "n_by_contrast", "q_interpretation": "discovery_fdr"|"reproduction_only"}, "presence": {"detection_available": bool, "bh_n_by_contrast", "reason"?}, "prior": {"s2_prior", "df_prior"} | None, "design": {"columns", "factors", "covariates", "block"}, "group_n": {group: 個体数}, "subject_source", "confounding": {...}, "warnings": [...], "rows": [...], "omnibus": [...] | None}`
  - `rows` は特徴 × 対比ごと: `feature_id, contrast_id, reference, test, status ("tested"|"not_tested"), reason, n_reference, n_test, detected_reference, detected_test, quantity_judgeable, log2fc, ci_low, ci_high, ci_unavailable_reason, statistic, df, p, q, presence_p, presence_q, verdict, verdict_reason`
  - `omnibus`（cell が 3 つ以上か対比が 2 つ以上で、`moderated_t` のとき）: `feature_id, f, f_p, f_q`（BH は omnibus の中で別に）

規則:
- 群の個体数（対比に出る各 cell）が 2 未満なら `STATISTIC_GROUP_TOO_SMALL`（個体でまとめた後で数える）。
- 特徴ごとに、モデルのどれかの cell で有限値の個体が 2 未満 → その特徴は全対比で `not_tested`（`insufficient_units`）、当てはめと事前分布から外す。
- `welch`・`student`・`mannwhitney` は 1 因子・共変量なし・ブロックなしだけ（さもなくば `DomainError("STATISTIC_SPECIFICATION_INVALID")`）。各対比で `spec_statistics.two_group(ref, test, kind=…, transform="none")` を log2 の値に呼ぶ。区間は student（等分散）と welch（Welch–Satterthwaite）で出し、mannwhitney は `ci_unavailable_reason="rank_test"`。
- `log2fc` は `log2_mean_difference`（moderated t ではモデルの対比の係数、他は群の log2 平均の差）。`arithmetic_mean_ratio` は `log2(mean(2**test)/mean(2**ref))`、区間なし（`ci_unavailable_reason="arithmetic_mean_ratio"`）。正なら test が高い。
- 出方: `detected_units` があれば対比ごとに `detection_fisher_p`、`presence_q` は対比ごとの BH（量とは別）。無ければ `presence.detection_available=false`、`quantity_judgeable=None`、理由 `detected_mask_unavailable`。
- `bh_population`: `matrix` は検定した全特徴、`annotated` は `prepared.annotated_features` との積、`reported` は `reported_feature_ids`（必須）との積で `q_interpretation="reproduction_only"`。母集団の外の q は None。母集団の大きさを対比ごとに書く。
- 判定は §15 の表の優先順。`verdict_reason` は `not_tested`・`quantity_not_judgeable`・`q_below_threshold`・`ci_within_equivalence`・`no_ci`・`neither`。
- 交絡: `confounding_report(prepared.units, 対比に出る群)` を必ず載せ、`assert_not_confounded` を検定の前に呼ぶ。
- 大きさの閾値（|log2FC| ≥ 1）は判定に使わない（この関数には引数を置かない）。

- [ ] **Step 1: テストを書く**

```python
def test_moderated_t_two_groups_matches_limma_reference(): Task 2 の cell_means の入力を PreparedMatrix に包み、B_vs_A の t・p・CI が limma_expected と一致
def test_verdicts_changed_unchanged_inconclusive(): 大きな差 → changed、n が大きく差 ~0 → unchanged、n 小さく差 ~0 → inconclusive（「有意でない」を unchanged にしない）
def test_quantity_flag_forces_inconclusive_even_if_significant()
def test_feature_absent_in_one_group_is_presence_only()       # Review Focus 2
def test_no_detected_mask_runs_quantity_only()                # Review Focus 1
def test_replicates_collapse_before_group_size_check()        # Review Focus 3
def test_presence_bh_is_separate_from_quantity_bh(): presence_q の母集団の大きさが別に数えられる
def test_bh_population_reported_marks_reproduction_only(); test_bh_population_annotated_uses_annotated_features()
def test_welch_option_uses_spec_statistics_numbers(): p が spec_statistics.two_group(kind="welch") と一致、区間あり
def test_mannwhitney_and_arithmetic_ratio_have_no_ci_and_never_unchanged()
def test_three_groups_omnibus_f_matches_limma(): F・F.p.value が limma_expected と一致、f_q あり
def test_two_by_two_interaction_contrast_runs()
def test_batch_covariate_matches_limma_reference()
def test_complete_confounding_raises_unless_allowed(); test_result_always_has_batch_table()
def test_result_states_definitions(): effect_size_definition・test・bh_population.name・matrix_id・ci_level が必ずある
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_compare.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): 比較のエンジン（moderated t / F・区間・判定・出方と量の分離・BH の母集団・交絡）`

---

### Task 8: PCA（§10）

**Files:**
- Modify: `metabolomix/analysis/pca.py`（`run_pca_fit_subset` に `scaling: str = "autoscale"` を足す。`"pareto"` は SD の平方根で割る。既定の挙動は変えない）
- Create: `metabolomix/analysis/stats/pca.py`, `tests/test_stats_pca.py`

**Interfaces:**
- Consumes: `PreparedMatrix`・`impute_for_pca`（Task 6）、`analysis.pca.run_pca_fit_subset`
- Produces: `pca(prepared: PreparedMatrix, *, scaling: str = "autoscale", n_components: int | None = 5, impute: str = "half_min") -> dict`
  → `{"schema": "stats-pca.v1", "matrix_id", "scaling", "log2": True, "imputation": {"method", "imputed_cells", "used": bool}, "n_fit": 生体の個体数, "unit_ids", "qc_ids", "scores": {"units": [[...]], "qc": [[...]]}, "explained_variance_ratio", "loadings", "kept_features"}`。
  軸は生体の個体だけで決め（`fit_mask`）、QC は同じ軸へ投影。補完の基準行は生体の個体（QC は同じ値で埋める）。

- [ ] **Step 1: テストを書く**

```python
def test_run_pca_fit_subset_default_unchanged(): 既定の結果が今までと同じ（scaling 引数なしと "autoscale" が一致）
def test_pareto_scaling_divides_by_sqrt_sd()
def test_axes_fit_on_biological_units_and_qc_projected(): QC を外れ値にしても軸が変わらない
def test_imputation_recorded_when_used(); test_no_imputation_when_complete()
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_pca.py tests/test_analysis_pca.py tests/test_pca_scores.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): 前処理後の行列での PCA（生体で軸・QC を投影・Pareto）`

---

### Task 9: 組成（§11）

**Files:**
- Create: `metabolomix/analysis/stats/composition.py`, `tests/test_stats_composition.py`

**Interfaces:**
- Consumes: `PreparedMatrix`（Task 6）、`compare`（Task 7）
- Produces: `composition(prepared: PreparedMatrix, *, assay_kind: str, factors, comparisons, **compare_kwargs) -> dict`
  → `{"schema": "stats-composition.v1", "class_label_kind": "lipid_class"|"msdial_ontology", "classes": {class: [feature_ids]}, "unclassified_features": n, "zero_handling": {"rule": "exclude_class_with_nonpositive_total", "excluded_classes": {class: n_units}}, "clr": <compare の結果>, "absolute": <compare の結果>, "unit": prepared.quant.unit | None}`。
  クラスは `prepared.feature_classes`（MS-DIAL の Ontology）。CLR は `log2(total) − mean_classes(log2(total))`（個体ごと、除いたクラスを除いて）。
  CLR と絶対量の比較はどちらも `compare` を `bh_population="matrix"` で（特徴＝クラス、`detected_units=None`）。`matrix_id` は元の id に `:clr`・`:class_total` を付ける。

- [ ] **Step 1: テストを書く**

```python
def test_clr_rows_sum_to_zero_per_unit()
def test_class_with_zero_total_is_excluded_and_reported()
def test_label_kind_follows_assay_kind(): lipid → lipid_class、metabolite/unknown → msdial_ontology
def test_composition_shift_detected_by_clr_not_by_absolute(): 全体量だけ 2 倍の群 → 絶対量は changed、CLR は changed にならない
def test_absolute_analysis_reports_unit_from_quant_state()
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_composition.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): クラスの組成（合計の CLR と絶対量を並べる）`

---

### Task 10: analysis-spec.v2（§8 の primary・§12 の新しい版）

**Files:**
- Create: `metabolomix/analysis/stats/spec_v2.py`, `tests/test_stats_spec_v2.py`

**Interfaces:**
- Consumes: `from_dataset`（Task 4）、`normalize_preprocess_config`・`preprocess`・`NeedsInput`（Task 6）、`compare`（Task 7）、`pca`（Task 8）、`composition`（Task 9）、`atomic_io.canonical_hash`
- Produces: `ANALYSIS_SCHEMA_V2 = "analysis-spec.v2"`、`RESULT_SCHEMA_V2 = "analysis-result.v2"`、
  `normalize_spec_v2(spec: dict) -> dict`、`rule_hash_v2(normalized: dict) -> str`、
  `run_spec_v2(normalized: dict, ds, *, assay_kind: str = "unknown") -> dict`（`NeedsInput` のときは `{"status": "needs_input", ...}` をそのまま返す）

spec v2 の欄（未知のキーは `ANALYSIS_SPEC_INVALID`。既定は埋めて返す）:

| 欄 | 形・既定 |
|---|---|
| `schema` | `"analysis-spec.v2"`（必須） |
| `input` | `{"source": "dataset"}`（既定。v2 の中核はデータセットだけ。TSV の入力は配線段で要るか決める） |
| `preprocess` | Task 6 の設定（既定は `PREPROCESS_DEFAULTS`） |
| `design` | `{"factors": ["group"], "covariates": [], "block": null}` |
| `comparisons` | Task 3 の形のリスト（必須、1 つ以上、`id` は英数字と `_ . -`、重複不可） |
| `analyses` | リスト（必須）。各 `{"analysis_id", "kind": "compare"|"composition"|"pca", ...}`。`compare`/`composition` は `test`（既定 `moderated_t`）・`effect_size`（既定 `log2_mean_difference`）・`bh_population`（既定 `matrix`）・`reported_feature_ids`（`reported` のとき必須）・`verdict`（`{"q": 0.05, "equivalence_fc": 1.2}`）・`min_detected_fraction`（0.5）・`allow_confounded`（false）。`pca` は `scaling`（`autoscale`）・`n_components`（5）・`impute`（`half_min`） |
| `primary` | `compare`/`composition` の解析が 2 つ以上なら必須（その `analysis_id`）。1 つならそれが既定 |

結果: `{"schema": "analysis-result.v2", "spec": 正規化した spec, "rule_hash", "input_fingerprint": {"dataset_hash", "sample_table_hash"}, "matrix_id", "preprocess": {"steps", "warnings", "subject_source"}, "primary", "analyses": [全部。選んだものだけにしない]}`。`bh_population="annotated"` を `spec` で宣言していなければ使わせない（既定 `matrix`）は、この正規化が保証する。

- [ ] **Step 1: テストを書く**

```python
def test_defaults_filled_and_unknown_keys_rejected()
def test_primary_required_when_several_analyses(); test_all_analyses_reported_not_only_primary()
def test_reported_requires_feature_ids()
def test_same_spec_same_dataset_same_result(): 2 回実行して JSON が同じ（決定論）
def test_v1_schema_rejected_here(): schema="analysis-spec.v1" → ANALYSIS_SPEC_INVALID（v1 は analysis_spec.run_spec の担当）
def test_needs_input_passthrough_for_is_normalized_dataset_with_qc()
def test_end_to_end_on_fake_dataset(): SimpleNamespace の ds（detected_mask・metadata rows・subject_id 付き）で compare + pca + composition が通り、定義の欄がそろう
```

- [ ] **Step 2: 実装して通す** — Run: `pytest tests/test_stats_spec_v2.py -q` → PASS
- [ ] **Step 3: Commit** — `feat(stats): analysis-spec.v2 の正規化と実行（primary・全解析を載せる）`

---

### Task 11: 記録と全体の確認

**Files:**
- Modify: `metabolomix/analysis/stats/__init__.py`（`__all__` を確定）
- 追記（main ツリー・追跡外・末尾だけ）: `C:/Users/yuu18/Metabolomix_with_LLM/docs/HISTRY.md`、`C:/Users/yuu18/Metabolomix_with_LLM/docs/task.md`
- 共有ノート: `<side-session-notes>/c-core-handover.md`（B-2・B-3 向け: `quant_state_of` の受け口の名前、配線段の一覧）

- [ ] **Step 1: 全テスト** — Run: `<venv python> -m pytest tests -q`（バックグラウンド、出力をファイルへ）→ 失敗 0。
- [ ] **Step 2: HISTRY に追記**: 作ったもの、§15 の判断（この計画書の表へのリンク）、`pca.py` の変更、R の参照値の作り方とバージョン（または突き合わせなかったこと）。
- [ ] **Step 3: task.md の末尾に配線段を追記**（下の「配線段」）。
- [ ] **Step 4: Commit**（バックグラウンド）— `docs(stats): C 中核の __all__ と記録`

---

## 配線段（この計画の範囲外。B-2 のマージ後の別段）

- §14 のツール: `dataset_preprocess`（`preprocess` → 行列の ID を `analysis_matrices` に登録）、`dataset_compare`（単純な引数 → spec v2 → `run_spec_v2`、戻り値は判定ごとの件数・上位・注意だけ、全量は `result_state.register_result`）、`dataset_pca`、`dataset_analysis_run`（`schema` で v1 → `analysis_spec.run_spec`、v2 → `run_spec_v2` に振り分け）、`pathway_project`（`pathway/compare.py` の Welch をエンジンへ。BH はパスウェイ × スコープ種別ごとのまま）。
- `needs_input` の封筒を MCP の戻り値に載せる（`mcp_errors` の形との合わせ方）。
- 対話の v1 の削除: `differential.two_group_test`・`statistics_v2` の welch・`multigroup` の ANOVA の呼び出し元をエンジンへ寄せて消す（§15「並べて残さない」）。analysis-spec.v1 の実行系は §12 のため残す。
- §13: `export_contract` にメタ行（効果量の定義・検定・判定の幅）を足す（`CONTRACT_VERSION` は上げない）。根拠バンドル v3 の manifest へも（massbank-context の `experiment/bundle.py` が未知の欄を許すか先に確かめる）。
- 前処理と解析の記録を解析フォルダの `reports/` へ書く（全体 §5）。
- B-3 の正規化済みの状態は `ds.quant_state`（`lipidmix.quant-state.v1`。B-3 と合意済み）。B-3 のマージ後に、B-3 の実物の ds で `quant_state_of` を確かめる。
- `subject_id` の書き込み（`dataset_set_sample_metadata`、B-2）。
- 文書: `docs/output_format/reanalysis.md`（v2 の結果の欄）、`USAGE.md`、CLAUDE.md の構成の節（`analysis/stats/`）、`docs/workflow/`、vault の流れ図。
