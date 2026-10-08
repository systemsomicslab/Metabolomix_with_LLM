"""DatasetState: mzTab-M を読んだ後の正準モデル。spec §11 参照。

session.arf の ARF 専用構造とは完全に独立している。session.dataset スロットへ格納する。
"""
from __future__ import annotations

import hashlib
import re
import uuid
from pathlib import Path

import numpy as np

from metabolomix.mztab import identity as mztab_identity
from metabolomix.mztab.identity import derive_inchikey
from metabolomix.mztab.reader import extract_abundance_matrix
from metabolomix.mztab.validator import detect_quantification_measure, validate_mztab

# assay メタデータ MTD 行のキーパターン。`assay[N]-<suffix>`（ms_run_ref 等）と、
# 表示名を運ぶ素の `assay[N]` を区別する。
_ASSAY_BARE_RE = re.compile(r"^assay\[(\d+)\]$")
_ASSAY_SUFFIX_RE = re.compile(r"^assay\[(\d+)\]-(.+)$")
# abundance 列名から assay 番号を取り出す（reader.py の _ABUNDANCE_RE と同じ規則）。
_ABUNDANCE_ASSAY_RE = re.compile(r"abundance_assay\[(\d+)\]", re.IGNORECASE)
# assay[N]-custom[M] が運ぶ CV term（MS-DIAL の実出力形式）。
_CV_TERM_RE = re.compile(r"^\[([^\]]*)\]$")
_ACCESSION_INJECTION_ORDER = "MS:4000089"
_ACCESSION_BATCH = "MS:4000088"
# `MTD id_confidence_measure[N]` の宣言行パターン。N は固定順が既定だが、
# `manualAssigned` があると上流（MztabFormatExport.cs `SetIdConfidenceMeasure`）が
# 9本目を追加するため、列を N の位置決め打ちで読んではいけない——ここで宣言行から
# 名前を引き、対応表でキーへ変換する。
_MEASURE_DECLARATION_RE = re.compile(r"^id_confidence_measure\[(\d+)\]$")

#: 宣言名 → スネークケースキーの対応表。[1] の総合スコアは best_id_confidence_value /
#: best_id_confidence_measure で既に持っているので、ここでは "total_score" として
#: 識別だけしておき、_extract_confidence_measures 側で除外する。
_CONFIDENCE_MEASURE_NAMES: dict[str, str] = {
    "MS-DIAL algorithm matching score": "total_score",
    "Retention time similarity": "retention_time_similarity",
    "Retention index similarity": "retention_index_similarity",
    "m/z similarity": "mz_similarity",
    "Simple dot product": "simple_dot_product",
    "Weighted dot product": "weighted_dot_product",
    "Reverse dot product": "reverse_dot_product",
    "Matched peaks count": "matched_peaks_count",
    "Matched peaks percentage": "matched_peaks_percentage",
    "CCS similarity": "ccs_similarity",
}


def _parse_cv_term(value: str | None) -> tuple[str | None, str | None]:
    """`[CV,ACCESSION,name,value]` から accession と値を取り出す。"""
    m = _CV_TERM_RE.match((value or "").strip())
    if not m:
        return None, None
    parts = [part.strip() for part in m.group(1).split(",")]
    if len(parts) < 4:
        return None, None
    return parts[1], parts[3]


def _parse_measure_declaration_name(value: str | None) -> str | None:
    """`MTD id_confidence_measure[N]` の値 `[,, X, ]` から宣言名 X を取り出す。

    `_parse_cv_term` と構文は同じ `[a, b, c, d]` 形だが、拾う位置が違う
    （こちらは name が 3 番目=index 2。`_parse_cv_term` は accession/value 用）ので
    専用に持つ。
    """
    m = _CV_TERM_RE.match((value or "").strip())
    if not m:
        return None
    parts = [part.strip() for part in m.group(1).split(",")]
    if len(parts) < 3:
        return None
    return parts[2] or None


def _snake_case(name: str) -> str:
    """未知の measure 宣言名をスネークケース化する（対応表に無い名前を落とさないため）。

    英数字以外の連続をアンダースコアに畳み、小文字化するだけの単純な変換。
    ``"CCS similarity"`` → ``"ccs_similarity"`` のような既知の対応表エントリと
    同じ形に揃う程度で十分（厳密な自然言語処理はしない）。
    """
    return re.sub(r"[^0-9a-zA-Z]+", "_", name.strip()).strip("_").lower()


def _build_confidence_measure_key_map(metadata: dict) -> dict[str, str]:
    """`MTD id_confidence_measure[N]` の宣言群から `{N(str): スネークケースキー}` を作る。

    宣言が無い（`id_confidence_measure[N]` 行が無い、または名前を取り出せない）
    index は対応表に含めない——値があっても引けないので無視される。
    """
    key_map: dict[str, str] = {}
    for key, value in metadata.items():
        m = _MEASURE_DECLARATION_RE.match(key)
        if not m:
            continue
        name = _parse_measure_declaration_name(value)
        if not name:
            continue
        key_map[m.group(1)] = _CONFIDENCE_MEASURE_NAMES.get(name, _snake_case(name))
    return key_map


def _extract_confidence_measures(evidence: dict, measure_key_map: dict[str, str]) -> dict[str, float]:
    """SME 行（`evidence`）から `id_confidence_measure[2..]` のサブスコアを回収する。

    `[1]`（total_score）は `best_id_confidence_value` / `best_id_confidence_measure`
    が既に持っているので、ここでは含めない。値が `null`（パーサが None に正規化
    済み）の列も含めない——「取得したが空」と「無い」を区別する必要が無いぶんの
    単純化（この dict 自体が「取れたものだけ」を表す）。
    """
    if not evidence or not measure_key_map:
        return {}
    measures: dict[str, float] = {}
    for index, key in measure_key_map.items():
        if key == "total_score":
            continue
        value = _to_float(evidence.get(f"id_confidence_measure[{index}]"))
        if value is None:
            continue
        measures[key] = value
    return measures


class DatasetState:
    """mzTab-M 読み込み後の正準モデル。多変量解析の共通入口。"""

    def __init__(self):
        self.source_format: str = "mztab"
        self.source_files: dict[str, str] = {}          # path -> sha256
        self.quantification_measure: str | None = None  # peak_height | peak_area_above_zero
        self.quantification_confidence: str | None = None
        self.feature_matrix: np.ndarray | None = None   # shape (n_features, n_samples)
        self.sample_names: list[str] = []                # assay 表示名（無ければ abundance 列名にフォールバック）
        self.sample_assay_ids: list[str] = []            # sample_names と同じ列順の assay ID
        self.feature_ids: list[str] = []                # SMF_ID 列
        self.assay_metadata: dict = {}                  # assay_id -> MTD 情報
        self.feature_metadata: dict = {}                # smf_id -> {name, mz, rt, inchikey, ...}
        # feature_candidates: smf_id -> [SME候補, ...]（rank昇順）。
        # feature_metadata が持つのは rank 最上位の1件だけで、Task 7 の
        # feature binding は adduct / charge / ライブラリ識別子 / スコアで
        # 候補を選び直す必要がある——最上位だけを残すと、profile が要求する
        # adduct と違う候補が付いた特徴を「同定なし」としか言えなくなる。
        self.feature_candidates: dict = {}
        # feature_annotations: smf_id -> SML 由来のラベル。
        # **証拠ではない。** MS-DIAL は Text DB 由来の同定を SME に書かない
        # （MztabFormatExport.cs `ShouldWriteSmeLine` が
        # `IsTextDbBasedRepresentative` を除外する）ので、Text DB 運用では
        # SML が唯一の同定情報源になる。
        # feature_metadata / feature_candidates とは**別スロット**に隔離する——
        # feature_bindings._check_identity は候補ゼロのとき
        # feature_metadata["inchikey"] を見て matched を返すため、ここへ混ぜると
        # MS1 注釈だけで authentic_standard_match が通ってしまう。
        self.feature_annotations: dict = {}
        self.sml_rows: list[dict] | None = None
        self.sme_rows: list[dict] | None = None
        self.validation_result: dict = {}
        self.inchikey_coverage: dict = {}
        # evidence sidecar（spec §10.1）。mzTab-M は非ゼロ値が実測か gap-fill かを
        # 区別しないので、隣接 `.arf` から (特徴 × サンプル) の検出状態を補う。
        # 42,840 セル規模になるため dict ではなく bool 行列で持ち、戻り値には出さない。
        # feature_matrix と同じ (n_features, n_samples) の並び。取り込めなかった場合は
        # None のまま残し、feature_qc に理由を書く（0 件と混同させない）。
        self.detected_mask = None
        self.feature_qc: dict = {}

        # --- ジョブ由来フィールド（dataset_load(job_path=...) で設定） ---
        # job_path: analysis-job.json の絶対パス。mzTab-M を直接指定した場合は None。
        self.job_path: str | None = None
        # artifact_paths: role -> [絶対パス, ...] のマップ。
        # Console 出力の ARF / DCL / EIC へのルックアップに使う。
        self.artifact_paths: dict[str, list[str]] = {}

        # --- 解析状態フィールド（dataset_preprocess / dataset_pca / dataset_differential が設定） ---
        # pp_matrix: 前処理済み行列。shape は (n_samples, n_features) — feature_matrix の転置。
        # ARF 側の session.arf.feature_matrix / pp_sample_names / pp_feature_names と
        # 同じ役割で、スロットだけが独立している。
        self.pp_matrix = None
        self.pp_sample_names: list[str] = []
        self.pp_feature_names: list[str] = []
        # roles: {sample_name: "sample"|"qc"|"blank"}。preprocess() と差次的解析の群構成で使う。
        self.roles: dict[str, str] = {}
        # sample_meta: {sample_name: {role, batch, batch_source, run_order}}。
        # 交絡判定（群⟂バッチ）とドリフト補正の材料。
        self.sample_meta: dict = {}
        self.preprocessing_recipe: dict = {}
        # last_pca / last_differential: 直近結果の全量。戻り値には要約だけを載せ、
        # 全量はここに置く（CLAUDE.md の戻り値肥大禁止）。
        # 現時点で読むのは dataset_export_differential のみ。図の保存ツール
        # （save_figure の pca / volcano）は session.arf 側を見ており、
        # DatasetState 経路には未対応（次フェーズ）。
        self.last_pca = None
        self.last_differential = None

        # --- 同一性と来歴（metabolomix/analysis/result_state.py が読み書きする） ---
        # dataset_id: この DatasetState 実体の ID。結果の provenance が指す先で、
        # 「別のデータセットで計算した結果」を持ち込ませないための鍵になる。
        self.dataset_id: str = f"ds_{uuid.uuid4().hex}"
        # metadata_revision: 実験情報（サンプルメタデータ）を適用するたびに増える。
        self.metadata_revision: int = 0
        # sample_metadata_rows: 明示的に与えられた実験情報シートの行（Task 9）。
        # 無いうちは None で、役割・バッチはサンプル名からの推定に頼る。
        self.sample_metadata_rows = None
        # sample_ids: sample_metadata_rowsのsample_id列を、sample_namesと同じ順で
        # 保持する安定ID。apply_metadataだけが書く。表示名（sample_names）は
        # 上書きしない——比較定義・出力はsample_idsを、表示・群選択の見た目は
        # 引き続きsample_namesを使い分けられるようにするための別スロット。
        self.sample_ids: list[str] = []
        # preprocess_id / preprocess_metadata_hash: 現在の前処理済み行列がどの計算から
        # 出たか。派生結果（PCA・差次的解析）はこの ID を parent_ids に持つ。
        self.preprocess_id: str | None = None
        self.preprocess_metadata_hash: str | None = None
        # results: result_id -> 結果全量。図・エクスポートが「どの結果か」を
        # 名指しで選べるようにする（last_* は現在の既定を指すだけの別名）。
        self.results: dict = {}
        # analysis_matrices: matrix_id -> analysis-matrix.v1（Task 8）。
        # v2 は recipe ごとに行列を持つので、`pp_matrix` のような単一スロットへは
        # 置けない——2本目の recipe が1本目を黙って上書きする。統計は
        # matrix_id で名指しして参照する。
        self.analysis_matrices: dict = {}

        # --- 出所の信用度（metabolomix/mztab/loading.py が設定する） ---
        # source_verification: verified（終了証跡と hash で裏取り済み）/
        # legacy_unverified（ジョブ経由だが証跡が無い・一致しない）/
        # direct_unverified（.mzTab を直接読んだ）。pipeline は verified だけを受ける。
        self.source_verification: str = "direct_unverified"
        # exploratory_only: 完了していない実行の出力。探索はできるが、2 群比較と
        # 差次的エクスポートは拒否する（欠けた検体を「その群に無い」と読み違えるため）。
        self.exploratory_only: bool = False
        # assay_sources: {abundance 列名: 実行時に予定した raw のパス}。
        # 完了ゲート（console/validation.map_assays）と同じ対応関係。
        self.assay_sources: dict = {}


def _sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def build_dataset_state(
    parse_result: dict,
    filename: str,
    source_path: str | Path,
) -> DatasetState:
    """parse_mztab() の戻り値から DatasetState を構築する。"""
    ds = DatasetState()

    # ファイルハッシュ
    p = Path(source_path)
    if p.is_file():
        ds.source_files[str(p)] = _sha256(p)

    # バリデーション
    ds.validation_result = validate_mztab(parse_result)

    # 定量種別
    measure, confidence = detect_quantification_measure(parse_result, filename)
    ds.quantification_measure = measure
    ds.quantification_confidence = confidence

    # abundance 行列
    matrix, sample_names, feature_ids = extract_abundance_matrix(parse_result)
    ds.feature_matrix = matrix
    ds.sample_names = sample_names
    ds.feature_ids = feature_ids

    # SMF メタデータ + InChIKey 導出
    #
    # mzTab-M 2.0.0-M では **構造・名称は SME セクションにしか無い**。SMF が持つのは
    # SMF_ID / SME_ID_REFS / exp_mass_to_charge / retention_time_in_seconds /
    # abundance_assay[N] で、database_identifier・smiles・inchi・chemical_name は
    # SME 専用の列である。SMF 行からこれらを読もうとすると常に None になり、
    # 「この測定には同定が無い」と誤読される（実データで InChIKey 0/714 になっていた）。
    smf_rows = parse_result["sections"].get("SMF", {}).get("rows", [])
    sme_by_id = _index_sme_rows(parse_result["sections"].get("SME", {}).get("rows", []))
    by_source: dict[str, int] = {"database_identifier": 0, "inchi_derived": 0, "smiles_derived": 0, "none": 0}
    derivable_but_missing = 0
    # `id_confidence_measure[N]` の宣言 → スネークケースキー。SME 行が持つ個別
    # スコアを後で feature_annotations（SML 由来）へ合流させるため、fid ごとに
    # 拾っておく（feature_annotations は SML 行だけから作られるので、SME の
    # サブスコアはここで一旦拾って後段でマージする）。
    measure_key_map = _build_confidence_measure_key_map(parse_result.get("metadata", {}))
    confidence_measures_by_fid: dict[str, dict[str, float]] = {}
    confidence_value_by_fid: dict[str, float] = {}
    for row in smf_rows:
        fid = row.get("SMF_ID", "")
        candidates = _all_candidates(row.get("SME_ID_REFS"), sme_by_id)
        ds.feature_candidates[fid] = candidates
        evidence = _best_evidence(row.get("SME_ID_REFS"), sme_by_id)
        ik, src = derive_inchikey(
            evidence.get("database_identifier"),
            evidence.get("inchi"),
            evidence.get("smiles"),
        )
        ds.feature_metadata[fid] = {
            "name": evidence.get("chemical_name"),
            "mz": _to_float(row.get("exp_mass_to_charge")),
            # mzTab-M は RT を**秒**で持つ。ARF 経路は分で持ち、両者は同じ
            # エクスポート契約の rt 列を共有するので、ここで分へ揃える。
            "rt": _seconds_to_minutes(_to_float(row.get("retention_time_in_seconds"))),
            "inchikey": ik,
            "inchikey_source": src,
            "smiles": evidence.get("smiles"),
            "inchi": evidence.get("inchi"),
        }
        by_source[src] = by_source.get(src, 0) + 1
        # RDKit があれば拾えたはずの取りこぼし。構造が無い行は RDKit の有無に
        # 関係なく InChIKey を持てないので、ここには数えない。
        if ik is None and (evidence.get("smiles") or evidence.get("inchi")):
            derivable_but_missing += 1

        measures = _extract_confidence_measures(evidence, measure_key_map)
        if measures:
            confidence_measures_by_fid[fid] = measures
        total_value = _to_float(evidence.get("best_id_confidence_value")) if evidence else None
        if total_value is not None:
            confidence_value_by_fid[fid] = total_value

    with_ik = sum(v for k, v in by_source.items() if k != "none")
    has_rdkit = mztab_identity.rdkit_available()
    ds.inchikey_coverage = {
        "total_features": len(smf_rows),
        "with_inchikey": with_ik,
        "by_source": by_source,
        "rdkit_available": has_rdkit,
    }

    # SML / SME 行
    ds.sml_rows = parse_result["sections"].get("SML", {}).get("rows")
    ds.sme_rows = parse_result["sections"].get("SME", {}).get("rows")

    # SML 由来のラベル（証拠ではない）。証拠スロットの構築が終わってから作る。
    ds.feature_annotations, annotation_warnings = _build_feature_annotations(
        ds.sml_rows, set(ds.feature_ids))
    if annotation_warnings:
        existing = ds.validation_result.setdefault("warnings", [])
        ds.validation_result["warnings"] = [*existing, *annotation_warnings]

    # SME の個別スコアを SML 由来の注釈へ合流させる。**既存の読み取りは変えない**
    # ——SML 行が自前で best_id_confidence_value を持っていればそれを優先し、
    # 欠けているときだけ（実 mzTab で稀ではない: Text DB 経由でない同定は
    # スコアが SME 側にしかない）SME 側の total を補う。曖昧（ambiguous）な
    # 注釈は「どれが正しいか決められない」ので触らない。
    for fid, annotation in ds.feature_annotations.items():
        if annotation.get("ambiguous"):
            continue
        measures = confidence_measures_by_fid.get(fid)
        if measures:
            annotation["confidence_measures"] = measures
        if annotation.get("confidence_value") is None:
            fallback = confidence_value_by_fid.get(fid)
            if fallback is not None:
                annotation["confidence_value"] = fallback

    # 同定の出所内訳。MS1 注釈を MS/MS 裏付けと取り違えないために分けて数える。
    sme_named = sum(1 for m in ds.feature_metadata.values() if m.get("name"))
    sml_only = sum(
        1 for fid, a in ds.feature_annotations.items()
        if a.get("name") and not (ds.feature_metadata.get(fid) or {}).get("name"))
    ds.inchikey_coverage["identified_by"] = {
        "sme": sme_named,
        "sml_only": sml_only,
        "none": len(ds.feature_ids) - sme_named - sml_only,
    }
    # SME が InChIKey を出せなかった特徴を SML が補ったぶんを足す。
    # by_source は全特徴で 1 回ずつ数える不変条件を保つため、'none' から移す。
    for fid, a in ds.feature_annotations.items():
        if not a.get("inchikey"):
            continue
        if (ds.feature_metadata.get(fid) or {}).get("inchikey"):
            continue
        src = a.get("inchikey_source") or "none"
        ds.inchikey_coverage["by_source"]["none"] -= 1
        ds.inchikey_coverage["by_source"][src] = (
            ds.inchikey_coverage["by_source"].get(src, 0) + 1)
        ds.inchikey_coverage["with_inchikey"] += 1

    # assay メタデータ
    # `assay[N]-<suffix>` 形式（ms_run_ref 等）に加え、素の `assay[N]` 行も拾う。
    # 素の行が MS-DIAL の表示名（例: 20220901_RAW_control_0h_1_NEG）を持つ唯一の場所で、
    # これを取りこぼすと sample_names が abundance_assay[N] という不透明な列識別子の
    # ままになり、群選択（dataset_differential）と QC/blank ロール検出
    # （detect_sample_roles はサンプル名のトークンを見る）が機能しなくなる。
    meta = parse_result.get("metadata", {})
    for k, v in meta.items():
        m = _ASSAY_SUFFIX_RE.match(k)
        if m:
            aid = f"assay[{m.group(1)}]"
            suffix = m.group(2)
            ds.assay_metadata.setdefault(aid, {})[suffix] = v
            if suffix.startswith("custom["):
                accession, term_value = _parse_cv_term(v)
                if accession == _ACCESSION_INJECTION_ORDER:
                    try:
                        ds.assay_metadata[aid]["run_order"] = int(str(term_value).strip())
                    except (TypeError, ValueError):
                        ds.assay_metadata[aid]["run_order"] = None
                elif accession == _ACCESSION_BATCH:
                    ds.assay_metadata[aid]["batch"] = term_value
            continue
        m = _ASSAY_BARE_RE.match(k)
        if m:
            aid = f"assay[{m.group(1)}]"
            # "name" キーで保持する。dataset_status など将来の呼び出し元は
            # assay_metadata[aid]["name"] を見れば表示名に到達できる。
            ds.assay_metadata.setdefault(aid, {})["name"] = v

    # abundance 列 → assay 表示名の解決。feature_matrix の列順（assay 番号昇順）は
    # extract_abundance_matrix が既に確定させているので、ここでは並べ替えず
    # 1:1 で置き換えるだけにする（列順を変えると全サンプルが黙って誤ラベルされる）。
    ds.sample_names, ds.sample_assay_ids, name_warnings = _resolve_sample_names(
        ds.sample_names, ds.assay_metadata)
    # **RDKit が実際に取り逃がした件数がある場合だけ**警告する。無条件に出すと、
    # database_identifier が InChIKey を持つ実データ（実測 162/271）でも
    # 「InChIKey は 0 件になる」と断言してしまい、同じ payload の
    # inchikey_coverage と矛盾する。警告の先頭スロットは実問題のために空けておく。
    if not has_rdkit and derivable_but_missing:
        name_warnings.insert(0, (
            f"RDKit が利用できないため、SMILES / InChI を持つ {derivable_but_missing} 件の"
            f"特徴量で InChIKey を導出できていません（InChIKey 付きは {with_ik}/"
            f"{len(smf_rows)} 件）。この {derivable_but_missing} 件は『同定が無い』のではなく"
            "『導出できていない』状態で、dataset_export_differential の書き出し対象から"
            "外れます（InChIKey が 0 件なら書き出し自体を拒否します）。"))
    if name_warnings:
        # **先頭に差す**。パーサ層の良性 warning（末尾空列の除去など）は実データで
        # 数百件になり得るので、後ろに append すると件数を絞って表示する
        # dataset_load の要約に載る余地がなくなる。表示されない warning は
        # 無いのと同じで、群選択が黙って誤るのを止められない。
        existing = ds.validation_result.setdefault("warnings", [])
        ds.validation_result["warnings"] = [*name_warnings, *existing]

    return ds


def _resolve_sample_names(abundance_cols: list[str], assay_metadata: dict) -> tuple[list[str], list[str], list[str]]:
    """abundance_assay[N] 列名を assay[N] の表示名へ解決する。

    表示名が無い assay（既存フィクスチャは全てこれに該当。現実のファイルでも
    起こり得る）は列識別子のままフォールバックする——意味のある名前が無いより、
    一意で追跡可能な旧識別子を残すほうが安全。

    表示名が複数 assay で重複する不正ファイルは、置き換え自体は行いつつ
    warning を返す（例外にはしない。読み込み自体を止めるほどではなく、
    群選択が曖昧になり得ることだけ呼び出し元に伝えれば足りる）。
    """
    resolved: list[str] = []
    assay_ids: list[str] = []
    cols_by_name: dict[str, list[str]] = {}
    for col in abundance_cols:
        m = _ABUNDANCE_ASSAY_RE.search(col)
        name = None
        aid = ""
        if m:
            aid = f"assay[{m.group(1)}]"
            name = (assay_metadata.get(aid) or {}).get("name")
        resolved_name = name if name else col
        resolved.append(resolved_name)
        assay_ids.append(aid)
        cols_by_name.setdefault(resolved_name, []).append(col)

    warnings = [
        f"assay 表示名が重複しています（{name!r}）: {', '.join(cols)}。"
        "sample_names での群選択が意図しないアッセイを指す恐れがあります。"
        for name, cols in cols_by_name.items() if len(cols) > 1
    ]
    return resolved, assay_ids, warnings

def _seconds_to_minutes(value: float | None) -> float | None:
    """mzTab-M の秒表記 RT を分へ揃える。"""
    return None if value is None else value / 60.0


def _index_sme_rows(sme_rows) -> dict:
    """SME 行を SME_ID で引けるようにする。"""
    return {str(r.get("SME_ID")): r for r in (sme_rows or []) if r.get("SME_ID") is not None}


def _sme_rank(row: dict) -> tuple[int, int]:
    """rank の昇順キー。rank が無い/数値でない証拠は最後に回す。

    mzTab-M の rank は 1 が最上位。同一特徴に複数の候補が付くのは常態なので、
    どれを採るかを暗黙にしない。
    """
    raw = row.get("rank")
    try:
        return (0, int(str(raw).strip()))
    except (TypeError, ValueError):
        return (1, 0)


def _best_evidence(refs, sme_by_id: dict) -> dict:
    """SME_ID_REFS が指す SME 行のうち、rank が最上位のものを返す。

    参照が無い（SMF_ID だけあって同定が付かなかった）特徴では空 dict を返す。
    呼び出し側はこれを「同定を取得していない」として扱う。
    """
    if not refs:
        return {}
    candidates = []
    for ref in str(refs).split("|"):
        row = sme_by_id.get(ref.strip())
        if row is not None:
            candidates.append(row)
    if not candidates:
        return {}
    return min(candidates, key=_sme_rank)


#: SME 行からそのまま写す識別・スコア列。mzTab-M 2.0.0-M の列名で、
#: MS-DIAL の出力で実際に埋まるものだけを並べる（無い列は None のまま残す
#: ——「取得していない」と「無い」を区別するため、キー自体は落とさない）。
_SME_TEXT_FIELDS = {
    "chemical_name": "chemical_name",
    "database_identifier": "database_identifier",
    "inchi": "inchi",
    "smiles": "smiles",
    "adduct": "adduct_ion",
    "identification_method": "identification_method",
    "confidence_measure": "best_id_confidence_measure",
    "spectra_ref": "spectra_ref",
}


def _all_candidates(refs, sme_by_id: dict) -> list[dict]:
    """SME_ID_REFS が指す全候補を rank 昇順で返す（同定が無ければ空リスト）。

    `_best_evidence` は最上位1件しか返さない。binding（Task 7）は profile が
    要求する adduct / charge / ライブラリ識別子 / スコア閾値で選び直すので、
    候補を捨てずに持つ。順序は rank 昇順（`_sme_rank` と同じ規則）で、
    rank が無い候補は末尾に回す。
    """
    if not refs:
        return []
    rows = []
    for ref in str(refs).split("|"):
        row = sme_by_id.get(ref.strip())
        if row is not None:
            rows.append(row)
    candidates = []
    for row in sorted(rows, key=_sme_rank):
        rank_key = _sme_rank(row)
        candidate = {
            "sme_id": str(row.get("SME_ID")),
            "rank": rank_key[1] if rank_key[0] == 0 else None,
            "charge": _to_int(row.get("charge")),
            "exp_mass_to_charge": _to_float(row.get("exp_mass_to_charge")),
            "theoretical_mass_to_charge": _to_float(
                row.get("theoretical_mass_to_charge")),
            "confidence_value": _to_float(row.get("best_id_confidence_value")),
        }
        for key, column in _SME_TEXT_FIELDS.items():
            candidate[key] = row.get(column)
        candidates.append(candidate)
    return candidates


#: SML 行からそのまま写す列（mzTab-M 2.0.0-M の列名）。
#: `inchi` は含めない——MS-DIAL は常に "null" を書く（MztabFormatExport.cs:393）。
_SML_TEXT_FIELDS = {
    "name": "chemical_name",
    "database_identifier": "database_identifier",
    "chemical_formula": "chemical_formula",
    "smiles": "smiles",
    "adduct": "adduct_ions",
    "reliability": "reliability",
    "confidence_measure": "best_id_confidence_measure",
}


def _build_feature_annotations(sml_rows, known_feature_ids: set) -> tuple[dict, list[str]]:
    """SML 行を feature_id（SMF_ID）ごとのラベルへ畳む。

    戻り値は `(annotations, warnings)`。警告は**種類**で集約する（行ごとに積むと
    実データで数百件になり、表示制限で重要な警告が埋もれる）。
    """
    by_feature: dict[str, list[dict]] = {}
    unknown_refs = 0
    for row in sml_rows or []:
        if not (row.get("chemical_name") or row.get("database_identifier")):
            continue
        for ref in str(row.get("SMF_ID_REFS") or "").split("|"):
            fid = ref.strip()
            if not fid:
                continue
            if fid not in known_feature_ids:
                unknown_refs += 1
                continue
            entry = {"sml_id": str(row.get("SML_ID")), "ambiguous": False}
            for key, column in _SML_TEXT_FIELDS.items():
                entry[key] = row.get(column)
            entry["confidence_value"] = _to_float(row.get("best_id_confidence_value"))
            # inchi は保持しないが、導出の材料としては渡す（他実装が書く余地を残す）。
            ik, src = derive_inchikey(row.get("database_identifier"),
                                      row.get("inchi"), row.get("smiles"))
            entry["inchikey"] = ik
            entry["inchikey_source"] = src
            by_feature.setdefault(fid, []).append(entry)

    annotations: dict[str, dict] = {}
    ambiguous = 0
    for fid, entries in by_feature.items():
        if len(entries) == 1:
            annotations[fid] = entries[0]
            continue
        ambiguous += 1
        annotations[fid] = {"ambiguous": True,
                            "sml_ids": [e["sml_id"] for e in entries],
                            "name": None}

    warnings: list[str] = []
    if ambiguous:
        warnings.append(
            f"1 つの特徴に複数の SML 注釈が当たっています（{ambiguous} 件）。"
            "どれが正しいか決められないため名前を付けていません。")
    if unknown_refs:
        warnings.append(
            f"SML の SMF_ID_REFS が存在しない特徴を指しています（{unknown_refs} 件）。")
    return annotations, warnings


def _to_int(v) -> int | None:
    try:
        return int(str(v).strip())
    except (TypeError, ValueError, AttributeError):
        return None


def _to_float(v: str | None) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None
