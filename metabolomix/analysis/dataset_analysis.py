"""DatasetState と analysis/ の純関数群をつなぐアダプタ（MCP 非依存）。

DatasetState.feature_matrix は (n_features, n_samples)。analysis/ の関数は
すべて (n_samples, n_features) を期待するため、ここで転置する。
呼び出し側はこの転置を意識しなくてよい。

処理順は metabolomix/arf/tools.py の arf_preprocess と同一に保つ:
  blank_filter → normalize → drift_correct → qc_rsd_filter → impute
  → drop_samples_by_role(blank)
ARF 経路と DatasetState 経路で違う数字が出ないことが、この層の存在意義。

依存は metabolomix.analysis.* のみ。arf/ mztab/ tools/ を import しない
（DatasetState は引数として受け取るだけで型 import もしない）。
"""
from __future__ import annotations

import re

import numpy as np

from metabolomix.analysis import differential, export_contract, preprocessing
from metabolomix.analysis.pca import run_pca

_DATE_RE = re.compile(r"(\d{8})")

# 差次的解析の群に混ぜてはいけないロール。
_NON_SAMPLE_ROLES = ("qc", "blank")

# ここで結果 dict に刻む版とラベルは、dataset_export_differential が検証する値と
# 同一でなければならない。ローカルに複製すると、export_contract.CONTRACT_VERSION を
# 上げた瞬間ここだけ古い値のままになり、エクスポートが「契約非互換」で永久に拒否され
# 続ける（dataset_differential を再実行しても同じ古い値を刻むだけなので、クライアント
# のリプレイでは直らない）。単一情報源として export_contract を直接参照する。


class PreconditionError(Exception):
    """前提が満たされないことを、理由付きで呼び出し側へ返す。

    MCP ツール層がこの例外を捕らえ、`kind` を見て封筒を選ぶ:
      kind="missing_state"  → missing_state() エンベロープ（リプレイで回復可能）
      kind="bad_request"    → 引数エラー。リプレイしても直らない
    None を返して理由を捨てると、上位が一律 missing_state に変換して
    クライアントを無限リプレイに落とす。

    kind="missing_state" のときは `state` に欠けている状態の識別子を入れる。
    ツール層はこれをそのまま missing_state() の第1引数に使う。
    **メッセージ文面から状態を推測させない**——文面を直した瞬間に振り分けが
    静かに壊れる結合になる。
    """

    def __init__(self, kind: str, message: str, details: dict | None = None,
                 state: str | None = None):
        super().__init__(message)
        self.kind = kind
        self.message = message
        self.details = details or {}
        if kind == "missing_state" and not state:
            raise ValueError("kind='missing_state' には state が必須です。")
        self.state = state


def _pp_inputs_from_explicit_metadata(sample_names: list[str], rows: list[dict]) -> tuple[dict, dict]:
    """Task 9 で適用済みの実験情報シートを、前処理入力（roles/sample_meta）へ変換する。

    `ds.sample_metadata_rows` は `sample_manifest.resolve_metadata` が
    `ds.sample_names` と同じ順に揃えてから `apply_metadata` が書き込むので、
    ここでは位置で対応させる（sample_id では結合しない——表示名を上書きしない
    設計と対称に、サンプル名の並びだけを唯一の対応軸にする）。

    ここでは各行の role/batch/injection_order を roles/sample_meta へ写すだけで、
    `row["include"]` は見ない——include=false を matrix/sample_names から落とす
    判断は呼び出し側（`build_dataset_pp_inputs`）が行毎に持つ（spec §8.1
    「include=falseの試料はこれらの評価前に除外する」）。ここで返す roles/
    sample_meta は除外対象も含めた全件分のまま返す。数値評価に渡す入力を絞る
    ことと、役割の写像（監査・表示用）を保持することは別の判断だからである。
    """
    roles: dict[str, str] = {}
    sample_meta: dict = {}
    for name, row in zip(sample_names, rows):
        # role はここでは "sample"/"qc"/"blank" へ絞り込まず、シートの値を
        # そのまま渡す（空欄だけ "sample" を既定にする）。明示 role="unknown"
        # （sample_manifest が正当な列挙値として保持する）はそのまま
        # preprocessing.preprocess() へ渡り、"blank" のように弾かれることも
        # "qc"/"sample" として扱われることもない——preprocessing.py は
        # role を `== "qc"` / `== "sample"` の完全一致と `drop_roles=("blank",)`
        # でしか見ないため、"unknown" は前処理を素通りして解析行列に残る。
        # role="unknown" を比較（compare_dataset）から除外するルールは Task 12
        # （resolve_comparison の比較ガード、spec §7.4）が持つ。ここで先回りして
        # 弾くと、その判定基準を2箇所に複製することになる。
        # 一方 include=false は「数値評価（前処理）からの除外」（spec §8.1）
        # ——role とは別の軸なので、ここでは弾かず、呼び出し側
        # （build_dataset_pp_inputs）が matrix/sample_names を絞る。
        role = row.get("role") or "sample"
        roles[name] = role
        provenance = row.get("provenance") or {}
        batch_prov = provenance.get("batch") or {}
        order_prov = provenance.get("injection_order") or {}
        sample_meta[name] = {
            "role": role,
            "batch": row.get("batch"),
            "batch_source": batch_prov.get("source"),
            "run_order": row.get("injection_order"),
            "run_order_source": order_prov.get("source"),
        }
    return roles, sample_meta


def _included_sample_indices(ds) -> list[int] | None:
    """明示メタデータで include=true の試料が占める ds.sample_names 上の位置。

    None は「絞り込む必要が無い」——明示メタデータが無い、行数が合わない
    （`build_dataset_pp_inputs` が別途 PreconditionError にする）、あるいは
    全件 include=true の場合。

    include=false の除外は行列（`build_dataset_pp_inputs`）と検出マスク
    （`_apply_detection_filter`）の**両方**に効かせる必要があるので、位置の
    計算はここ 1 か所に置く。片方だけに効かせると、除外した試料の未検出が
    残ったまま検出率が計算され、min_detection_rate=1.0 が正当な特徴量を
    削り落とす。
    """
    rows = getattr(ds, "sample_metadata_rows", None)
    if not rows:
        return None
    sample_names = list(getattr(ds, "sample_names", []) or [])
    if len(rows) != len(sample_names):
        return None
    keep = [i for i, row in enumerate(rows) if row.get("include", True)]
    return keep if len(keep) != len(rows) else None


def build_dataset_pp_inputs(ds):
    """DatasetState から preprocessing.preprocess() の引数を組む。

    Returns:
        matrix        : (n_samples, n_features) — feature_matrix の転置。
                        明示メタデータ適用済みなら include=false の行は含まない
                        （spec §8.1「include=falseの試料はこれらの評価前に除外する」。
                        下記 sample_names と同じ行だけを残す）。
        sample_names  : list[str] — 同上の理由で include=false のサンプル名は
                        含まない。
        feature_names : list[str]（SMF_ID）
        roles         : {sample_name: "sample"|"qc"|"blank"} — preprocess の第3引数。
                        **ただし** 明示メタデータ適用済み（`ds.sample_metadata_rows`
                        あり）の場合は、シートが持つ正当な列挙値 "unknown" もここへ
                        素通りする（"sample"/"qc"/"blank" の3値へ絞り込まない）。
                        preprocessing.py は role を `== "qc"` / `== "sample"` の
                        完全一致と `drop_roles=("blank",)` でしか見ないため、
                        "unknown" は「サンプルでもQCでもblankでもない」まま前処理を
                        通過する——現状はここでは弾かない。role="unknown" を比較
                        （2群比較）から除外するのは Task 12（`resolve_comparison`、
                        spec §7.4）の責務であり、本関数はそれを先取りしない。
                        **この辞書自体は include=false のサンプルも含めた全件分**
                        （数値評価の入力を絞ることと、監査・表示用の写像を保持
                        することは別の判断——`sample_names`/`matrix` には現れない
                        キーが残ることがある）。
        sample_meta   : {sample_name: {role, batch, batch_source, run_order,
                         run_order_source}}
                        — 交絡判定（群⟂バッチ）とドリフト補正の材料。
                        roles と同じく include=false 分も含めた全件分を返す。
    """
    if ds.feature_matrix is None or not ds.sample_names or not ds.feature_ids:
        raise PreconditionError(
            "missing_state",
            "DatasetState に定量行列がありません。dataset_load を先に実行してください。",
            state="dataset",
        )

    matrix = np.asarray(ds.feature_matrix, dtype=float).T.copy()
    sample_names = list(ds.sample_names)
    feature_names = list(ds.feature_ids)

    # 実験情報シート（Task 9 apply_metadata）が適用済みなら、そちらを優先する。
    # 未適用（sample_metadata_rows が None）のデータセットは、この下の
    # トークン判定 + mzTab MTD 推定という既存経路をそのまま使い続ける
    # （明示メタデータが無いデータセットの挙動を変えない）。
    explicit_rows = getattr(ds, "sample_metadata_rows", None)
    if explicit_rows:
        if len(explicit_rows) != len(sample_names):
            raise PreconditionError(
                "bad_request",
                "適用済みの実験情報シートの行数がサンプル数と一致しません。"
                "dataset_load をやり直してから metadata を再適用してください。",
                {"rows": len(explicit_rows), "samples": len(sample_names)},
            )
        roles, sample_meta = _pp_inputs_from_explicit_metadata(sample_names, explicit_rows)

        # include=false は数値評価（ブランク評価・正規化・ドリフト補正・RSD
        # フィルタ等）へ一切参加させない（spec §8.1「include=falseの試料はこれら
        # の評価前に除外する」）。Task 11 の resolve_policy が既に include=true
        # だけを集計してレシピの可否を決めているので（preprocess_policy.py
        # 「include=true 対象の集計」）、実行側もここで同じ集合に揃える——
        # 揃えないと、resolve_policy が「除外前提で安全」と判断した include=false
        # のQCが、実際にはPQN参照・ドリフト補正参照として使われてしまう
        # （controller裁定R15・レビュー Finding 1）。
        # roles/sample_meta の写像自体は削らない——除外理由の監査・表示
        # （dataset_status 等）は全件分の役割を必要とするため。
        keep_idx = _included_sample_indices(ds)
        if keep_idx is not None:
            matrix = matrix[keep_idx, :]
            sample_names = [sample_names[i] for i in keep_idx]
        return matrix, sample_names, feature_names, roles, sample_meta

    # class_ids は mzTab-M に対応物が無いので渡さない（既定 None）。
    roles = preprocessing.detect_sample_roles(sample_names)

    # 注入順とバッチは mzTab の MTD assay[N]-custom[...] が運ぶ
    # （MS:4000089 injection sequence label / MS:4000088 batch label）。
    # 実データの Console 出力・GUI 出力の両方に存在することを確認済み。
    assay_ids = list(getattr(ds, "sample_assay_ids", []) or [])
    assay_meta = getattr(ds, "assay_metadata", {}) or {}

    def _assay_field(index: int, key: str):
        aid = assay_ids[index] if index < len(assay_ids) else None
        return (assay_meta.get(aid) or {}).get(key) if aid else None

    mztab_batches = [_assay_field(i, "batch") for i in range(len(sample_names))]
    # MS-DIAL のバッチラベルは CSV インポートで指定しない限り全件 "1" になる。
    # 定数のラベルは情報を持たないので、その場合はファイル名日付の推定に戻す。
    # 定数を採ると、日付で分かれていた交絡が検出できなくなる。
    use_mztab_batch = len({b for b in mztab_batches if b is not None}) > 1

    sample_meta: dict = {}
    for i, name in enumerate(sample_names):
        m = _DATE_RE.search(name)
        if use_mztab_batch and mztab_batches[i] is not None:
            batch, batch_source = mztab_batches[i], "mztab_batch_label"
        elif m:
            batch, batch_source = m.group(1), "filename_date"
        else:
            batch, batch_source = None, None

        run_order = _assay_field(i, "run_order")
        sample_meta[name] = {
            "role": roles.get(name, "sample"),
            "batch": batch,
            "batch_source": batch_source,
            "run_order": run_order,
            "run_order_source": "mztab_injection_sequence" if run_order is not None else None,
        }
    return matrix, sample_names, feature_names, roles, sample_meta


def run_dataset_preprocess(ds, recipe: dict):
    """DatasetState の feature_matrix に前処理を適用する。

    recipe は preprocessing.preprocess() と同じキー（normalize / blank_min_fold /
    drift_correct / max_qc_rsd / impute）。未知キーは preprocess が無視する。

    Returns: (pp_matrix, pp_sample_names, pp_feature_names, roles, sample_meta, report)
    """
    matrix, sample_names, feature_names, roles, sample_meta = build_dataset_pp_inputs(ds)

    # 検出率フィルタは前処理より前に掛ける。正規化・補完のあとでは gap-fill セルが
    # 実測値と区別できなくなり、「何を根拠に残したか」が言えなくなる。
    matrix, feature_names, detection_report = _apply_detection_filter(
        ds, matrix, feature_names, recipe.get("min_detection_rate"),
        sample_indices=_included_sample_indices(ds))

    run_order = {n: sample_meta[n]["run_order"] for n in sample_names}

    try:
        matrix, kept_idx, report = preprocessing.preprocess(
            matrix, sample_names, roles, run_order, recipe,
        )
    except ValueError as exc:
        # normalize の未知メソッド等。引数由来なのでリプレイでは直らない。
        raise PreconditionError(
            "bad_request", f"前処理レシピが不正です: {exc}", {"recipe": recipe},
        ) from exc

    pp_feature_names = [feature_names[i] for i in kept_idx]
    report.update(detection_report)

    # ブランクは blank_filter の参照として使い終えたので解析行列から外す。
    # 残すと総強度が桁違いに低い行が PCA の PC1 を支配する（arf_preprocess と同じ理由）。
    # QC は残す——QC クラスタの締まり具合を PCA で見るのは品質確認の定番手段。
    matrix, pp_sample_names, dropped = preprocessing.drop_samples_by_role(
        matrix, sample_names, roles, drop_roles=("blank",),
    )
    report["excluded_from_matrix"] = dropped
    if dropped.get("blank"):
        report.setdefault("caveats", []).append(
            f"ブランク {len(dropped['blank'])} 件（{', '.join(dropped['blank'])}）は背景除去に"
            "使用後、解析行列（PCA/差次的解析）から除外しました。QC は PCA での品質確認の"
            "ため残しています。"
        )

    # プールQC が複数バッチに分かれているかの簡易警告（arf_preprocess:327,364 と同じ材料）。
    qc_batches = {sample_meta[n]["batch"] for n in sample_names
                  if sample_meta[n]["role"] == "qc"}
    if len(qc_batches) > 1:
        report.setdefault("caveats", []).append(
            "プールQC が複数バッチ/層に分かれています。全体一律のドリフト補正は近似です。"
        )

    # 層別プール QC の警告（arf_preprocess と同じ材料）
    qc_strata = preprocessing.detect_qc_strata(sample_names, roles)
    if len(qc_strata) > 1:
        labels = ", ".join(sorted(s for s in qc_strata if s))
        report.setdefault("caveats", []).append(
            f"プールQC が層別（{len(qc_strata)} サブグループ"
            f"{f': {labels}' if labels else ''}）と検出されました。"
            "全 QC を1系列として扱うドリフト補正/RSD フィルタは近似です。"
        )
    has_run_order = any(sample_meta[n]["run_order"] is not None for n in sample_names)
    if not has_run_order:
        report.setdefault("caveats", []).append(
            "この mzTab-M は注入順（MTD assay[N]-custom[...] の "
            "MS:4000089 injection sequence label）を持たないため、QC ドリフト補正は"
            "実施できません。注入順に依存する品質評価が必要なら ARF 経路"
            "（arf_preprocess）を使ってください。")
    elif recipe.get("drift_correct"):
        report.setdefault("caveats", []).append(
            "注入順は mzTab-M の injection sequence label から取得しました。"
            "MS-DIAL は CSV インポートで実注入順を与えない場合、**ファイル読み込み順**を"
            "そのまま注入順として書き出します。QC の挿入位置（drift_correct の "
            "qc_interspersion）で妥当性を確認してください。")

    return matrix, pp_sample_names, pp_feature_names, roles, sample_meta, report


#: 主成分と注入順の相関がこの絶対値を超えたら、分析ドリフトの疑いを caveat に出す。
_DRIFT_CORRELATION_WARN = 0.5


def run_dataset_pca(ds, n_components: int = 5, log_transform: bool = False) -> dict:
    """DatasetState の pp_matrix に PCA を実行する。"""
    matrix = _require_pp_matrix(ds)
    n_samples, n_features = matrix.shape
    if min(n_samples, n_features) < 2:
        raise PreconditionError(
            "bad_request",
            f"PCA には 2 以上のサンプルと特徴量が必要です"
            f"（現在: サンプル={n_samples}, 特徴量={n_features}）。"
            "前処理のフィルタ閾値が厳しすぎる可能性があります。",
            {"n_samples": n_samples, "n_features": n_features},
        )

    pca = run_pca(matrix, n_components=n_components, log_transform=log_transform)
    components = np.asarray(pca["components"], dtype=float)
    n_pc = components.shape[1]

    scores = []
    for i, name in enumerate(ds.pp_sample_names):
        row = {"name": name, "role": ds.roles.get(name, "sample")}
        for pc in range(n_pc):
            row[f"PC{pc + 1}"] = round(float(components[i, pc]), 4)
        scores.append(row)

    # QC が無いバッチでは qc_drift_correct も QC-RSD フィルタも動かない。
    # そこで注入順との相関を出す。補正はできなくても「その主成分が分析ドリフトを
    # 写しているか」は言えるので、生物学として読む前の歯止めになる。
    from metabolomix.analysis.preprocessing import run_order_correlation
    sample_meta = getattr(ds, "sample_meta", None) or {}
    run_order = {n: (sample_meta.get(n) or {}).get("run_order")
                 for n in ds.pp_sample_names}
    order_corr = run_order_correlation(components, ds.pp_sample_names, run_order)

    caveats: list[str] = []
    drifting = [(i + 1, r) for i, r in enumerate(order_corr)
                if r is not None and abs(r) >= _DRIFT_CORRELATION_WARN]
    if drifting:
        detail = "、".join(f"PC{pc} (r={r:+.2f})" for pc, r in drifting)
        caveats.append(
            f"{detail} が注入順と強く相関しています。分析ドリフトを写している"
            "可能性があるため、群差として読む前に確認してください。"
            "QC 試料があれば dataset_preprocess(drift_correct=True) で補正できます。")

    return {
        "explained_variance_ratio": [round(float(v), 4)
                                     for v in pca["explained_variance_ratio"]],
        "scores": scores,
        "run_order_correlation": order_corr,
        "caveats": caveats,
        "n_samples": n_samples,
        "n_features": n_features,
        "log_transform": log_transform,
        # loadings は特徴量数 × 主成分数で巨大になる。要約には載せず、
        # 呼び出し側がセッションに保持する分にだけ含める。
        "loadings": pca["loadings"],
        "singular_values": pca["singular_values"],
    }


def run_dataset_differential(
    ds,
    group_a_samples: list[str],
    group_b_samples: list[str],
    *,
    q_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    log_transform: bool = True,
    group_a_label: str = "group_a",
    group_b_label: str = "group_b",
) -> dict:
    """前処理済み DatasetState で 2 群比較を実行する。

    group_a_samples / group_b_samples は pp_sample_names に含まれるサンプル名。
    未知の名前・非 sample ロール・両群への重複指定はすべて caveat で名指しする
    （黙って落とすと群サイズが縮んだことに気付けない）。
    """
    matrix = _require_pp_matrix(ds)
    available = list(ds.pp_sample_names)
    caveats: list[str] = []

    idx_a, names_a = _resolve_group(group_a_samples, available, ds.roles,
                                    group_a_label, caveats)
    idx_b, names_b = _resolve_group(group_b_samples, available, ds.roles,
                                    group_b_label, caveats)

    # 正規化未適用の警告（arf_differential:846-848 と同じ判定）。DatasetState 側は
    # run_dataset_preprocess が recipe 自体を ds へ書き戻さないため、呼び出し側が
    # 事前に ds.preprocessing_recipe へ設定しておく契約。未設定/空は "none" 扱い。
    recipe = ds.preprocessing_recipe or {}
    if recipe.get("normalize", "none") == "none":
        caveats.append(
            "正規化が未適用のため log2FC は測定量差を含み得ます"
            "（dataset_preprocess の normalize を検討）。")

    overlap = sorted(set(names_a) & set(names_b))
    if overlap:
        raise PreconditionError(
            "bad_request",
            f"両群に同じサンプルが指定されています: {', '.join(overlap)}。"
            "群定義を見直してください。",
            {"overlap": overlap},
        )
    if len(idx_a) < 2 or len(idx_b) < 2:
        raise PreconditionError(
            "bad_request",
            f"群サイズ不足（{group_a_label}={len(idx_a)}, {group_b_label}={len(idx_b)}）: "
            "各群 n>=2 が必要です。群名の誤り、または前処理での試料脱落の可能性があります。",
            {"n_a": len(idx_a), "n_b": len(idx_b),
             "available_samples": available, "caveats": caveats},
        )

    # two_group_test はラベル列で群を切る（サンプル名リストではない）。
    # 両群のどちらにも属さない行は検定対象から外すため、行を抜いてラベルを組む。
    keep_idx = idx_a + idx_b
    sub_matrix = matrix[keep_idx, :]
    group_labels = [group_a_label] * len(idx_a) + [group_b_label] * len(idx_b)

    # 交絡（群⟂バッチ）判定は arf_differential:872-879 と同じく「実際に比較した2群」
    # （プール後）に対して行う。プール前の細粒度ラベルで判定すると偽の交絡警告が出る
    # （arf_differential:869-871 のコメントと同じ理由）。
    # names_a/names_b は idx_a/idx_b と同じ順で構築されているため、
    # keep_idx（idx_a + idx_b）と対応するバッチ列を同じ並びで組む。
    names_ordered = names_a + names_b
    batch_labels = [(ds.sample_meta.get(n) or {}).get("batch") for n in names_ordered]
    conf = differential.check_confounding(group_labels, batch_labels)
    # batch_source は比較対象に限らず全サンプルから拾う（arf_differential:858 と同じ）。
    batch_source = next(
        (m.get("batch_source") for m in (ds.sample_meta or {}).values()
         if m.get("batch_source")),
        None,
    )
    src_note = ("（バッチはファイル名の日付から推定。実バッチ設計と異なる場合あり）"
                if batch_source == "filename_date" else "")
    if conf["confounded"]:
        caveats.append("交絡: " + conf["detail"] + src_note)
    elif not conf.get("assessable", True):
        caveats.append("交絡評価不可: " + conf["detail"] + src_note)

    results = differential.two_group_test(
        sub_matrix, ds.pp_feature_names, group_labels,
        group_a_label, group_b_label, log_transform=log_transform,
    )
    results = differential.add_fdr(results)
    summary = differential.summarize_two_group(results, q_threshold, log2fc_threshold)
    volcano = differential.volcano_data(results, q_threshold, log2fc_threshold)

    n_tested = summary["n_tested"]
    if n_tested == 0:
        caveats.append(
            "検定可能な特徴が0件（全特徴で p=NaN）。群が空・分散0・または正規化で試料が"
            "NaN化した可能性があります。『有意0件』を『群間差なし』と解釈しないでください。")
    elif n_tested < 0.2 * len(ds.pp_feature_names):
        caveats.append(
            f"検定できた特徴は {n_tested}/{len(ds.pp_feature_names)} 件のみ（多くが p=NaN）。"
            "群内 n 不足・分散0・欠損が多い可能性があります（前処理の見直しを検討）。")
    if min(len(idx_a), len(idx_b)) < 4:
        caveats.append(
            f"小n（{group_a_label}={len(idx_a)}, {group_b_label}={len(idx_b)}）につき"
            "検出力が限られます。")
    caveats.append(
        f"log2FC の向き: 正なら {group_b_label} が高い（上昇）、負なら {group_a_label} が高い（低下）。")

    return {
        "kind": "two_group",
        "a": group_a_label,
        "b": group_b_label,
        "samples_a": names_a,
        "samples_b": names_b,
        "n_a": len(idx_a),
        "n_b": len(idx_b),
        "q_threshold": q_threshold,
        "log2fc_threshold": log2fc_threshold,
        "log_transform": log_transform,
        "contract_version": export_contract.CONTRACT_VERSION,
        "log2fc_sign": export_contract.LOG2FC_SIGN,
        "summary": summary,
        "caveats": caveats,
        # 全量。呼び出し側はこれをセッションに保持し、戻り値には載せない。
        "results": results,
        "volcano": volcano,
    }


# ---------- 内部ヘルパ ----------


def _restrict_mask_to_samples(ds, mask, sample_indices):
    """検出マスク (特徴 × サンプル) の列を、解析に残す試料だけへ絞る。

    `sample_indices` が None、マスクが 2 次元でない、列数が `ds.sample_names` と
    合わない（＝並びの対応を保証できない）場合は絞らずそのまま返す。列の対応が
    言えないマスクを位置で切ると、別の試料の検出状態を数えることになる。
    """
    arr = np.asarray(mask)
    if sample_indices is None or arr.ndim != 2:
        return arr
    if arr.shape[1] != len(list(getattr(ds, "sample_names", []) or [])):
        return arr
    return arr[:, list(sample_indices)]


def _apply_detection_filter(ds, matrix, feature_names, min_detection_rate, *,
                           sample_indices=None):
    """実検出率（gap-fill を除く）で特徴量を足切りし、検出状況を報告する。

    `matrix` は (サンプル × 特徴)、`ds.detected_mask` は (特徴 × サンプル) の並び。
    足切りは列（特徴）に対して行う。

    `sample_indices`（`_included_sample_indices`）を渡すと、マスクの列も
    その試料だけへ絞ってから検出率を数える。**絞らないと、解析から除外した
    試料（include=false）の未検出が検出率を押し下げる**——`matrix` 側は
    `build_dataset_pp_inputs` が既に include=true だけにしているのに、分母だけ
    全試料のままになり、`min_detection_rate=1.0` が「残す試料では全件検出
    されている」正当な特徴量を削り落とす。

    検出状態が無いのに閾値を渡された場合は例外にする。0 扱いで通すと「gap-fill
    だけの特徴を実測として数えた行列」が黙って下流に流れるため。
    """
    mask = getattr(ds, "detected_mask", None)
    requested = min_detection_rate is not None and min_detection_rate > 0.0

    if mask is None:
        if requested:
            raise PreconditionError(
                "bad_request",
                "この DatasetState は検出状態（gap-fill の区別）を持たないため "
                "min_detection_rate を適用できません。mzTab-M の非ゼロ値は実測ピークと "
                "gap-fill 補間値の区別を持たず、隣接する .arf からも取り込めませんでした"
                "（理由は dataset_status の detection を参照）。検出率での足切りが必要なら "
                "ARF 経路（arf_parser / arf_preprocess）を使ってください。",
                {"min_detection_rate": min_detection_rate,
                 "detection": (ds.feature_qc or {}).get("reason")},
            )
        return matrix, feature_names, {}

    mask = _restrict_mask_to_samples(ds, mask, sample_indices)
    rates = preprocessing.detection_rates(mask)
    n_cells = int(mask.size)
    n_detected = int(mask.sum())
    report = {"detection": {
        "source": (ds.feature_qc or {}).get("source"),
        "n_cells": n_cells,
        "n_detected": n_detected,
        "gap_filled_rate": round(1.0 - n_detected / n_cells, 4) if n_cells else None,
        # 何試料分を数えたか。include=false を除いた場合、全試料数とは一致しない。
        "n_samples": int(np.asarray(mask).shape[1]) if np.asarray(mask).ndim == 2 else None,
    }}
    if not requested:
        return matrix, feature_names, report

    keep = [i for i, rate in enumerate(rates) if rate >= min_detection_rate]
    removed = len(feature_names) - len(keep)
    report["detection_filter"] = {
        "min_detection_rate": min_detection_rate,
        "features_before": len(feature_names),
        "features_removed": removed,
    }
    if not keep:
        raise PreconditionError(
            "bad_request",
            f"min_detection_rate={min_detection_rate} を満たす特徴量が 0 件です"
            f"（{len(feature_names)} 件すべて除去）。閾値を下げてください。",
            report["detection_filter"],
        )
    return matrix[:, keep], [feature_names[i] for i in keep], report

def _require_pp_matrix(ds):
    if ds.pp_matrix is None:
        raise PreconditionError(
            "missing_state",
            "前処理済み行列がありません。dataset_preprocess を先に実行してください。",
            state="dataset_preprocessed",
        )
    return np.asarray(ds.pp_matrix, dtype=float)


def _resolve_group(requested, available, roles, label, caveats):
    """指定サンプル名を pp_sample_names の位置に解決し、落ちた分を caveat に残す。

    requested は先勝ちで重複除去する。ARF 側は class-ID/因子トークンから群を
    解決するため重複の起きようがないが、この経路は呼び出し側が生のサンプル名
    リストを渡す API なので、group_a=["S1", "S1"] のような取り違えが起こり得る。
    重複したまま通すと同一行を2回数えた群内分散ゼロの検定になり、有意性が
    水増しされたうえで誰にも気付かれない。
    """
    seen_requested: set[str] = set()
    duplicates: list[str] = []
    deduped_requested: list[str] = []
    for name in requested:
        if name in seen_requested:
            duplicates.append(name)
            continue
        seen_requested.add(name)
        deduped_requested.append(name)
    if duplicates:
        caveats.append(
            f"{label} に同じサンプル名が重複して指定されたため 1 回に丸めました: "
            f"{', '.join(duplicates)}。重複したまま検定すると同一行を二重に数え、"
            "群内分散を過小評価して有意性を水増しします。")

    index_of = {name: i for i, name in enumerate(available)}
    idx: list[int] = []
    names: list[str] = []
    unknown: list[str] = []
    non_sample: list[str] = []
    for name in deduped_requested:
        if name not in index_of:
            unknown.append(name)
            continue
        role = roles.get(name, "sample")
        if role in _NON_SAMPLE_ROLES:
            non_sample.append(f"{name}({role})")
            continue
        idx.append(index_of[name])
        names.append(name)
    if unknown:
        caveats.append(
            f"{label} に指定されたサンプルのうち {len(unknown)} 件は前処理済み行列に"
            f"存在しないため除外しました: {', '.join(unknown)}。"
            "名前の誤り、または前処理で脱落した可能性があります。")
    if non_sample:
        caveats.append(
            f"{label} から QC/ブランクを除外しました: {', '.join(non_sample)}。"
            "群に混ぜると比較が壊れるため、生体試料のみで検定します。")
    return idx, names
