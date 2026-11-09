"""DatasetState の試料の表を、役割・バッチ・注入順の写像に組むアダプタ（MCP 非依存）。

`build_dataset_pp_inputs` は DatasetState の feature_matrix を (n_samples, n_features) に転置し、
解析に残す試料（`included_sample_indices`）と、役割・バッチ・注入順（明示メタデータか、試料名と
mzTab-M の MTD から推定したもの）を返す。根拠バンドルの `samples.tsv`（`evidence_export`）が読む。
前処理と検定は統計エンジン（`analysis/stats/`）が持つ。

依存は metabolomix.analysis.* のみ。arf/ mztab/ tools/ を import しない
（DatasetState は引数として受け取るだけで型 import もしない）。
"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis import preprocessing, sample_manifest


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


def included_sample_indices(ds) -> list[int] | None:
    """解析に残す試料が占める ds.sample_names 上の位置。

    外す試料は 2 種類ある: 明示メタデータで include=false の試料と、`dataset_exclude` で
    除外した試料（`ds.excluded_samples`。名前で持つ）。

    None は「絞る必要が無い」——どちらも無い、または全件を残す場合。明示メタデータの行数が
    試料数と合わない場合（`build_dataset_pp_inputs` が別途 PreconditionError にする）は
    include の側を見ず、除外だけを適用する。

    行列（`build_dataset_pp_inputs`）・特徴検索の検出数・図の試料の並びが同じ集合を使うよう、
    位置の計算はここ 1 か所に置く。
    """
    sample_names = list(getattr(ds, "sample_names", []) or [])
    rows = getattr(ds, "sample_metadata_rows", None)
    use_rows = bool(rows) and len(rows) == len(sample_names)
    excluded = getattr(ds, "excluded_samples", None) or {}
    keep = [i for i, name in enumerate(sample_names)
            if name not in excluded and (not use_rows or rows[i].get("include", True))]
    return keep if len(keep) != len(sample_names) else None


def build_dataset_pp_inputs(ds):
    """DatasetState から、解析に残す試料の行列と役割・バッチ・注入順の写像を組む。

    Returns:
        matrix        : (n_samples, n_features) — feature_matrix の転置。
                        明示メタデータ適用済みなら include=false の行は含まない
                        （spec §8.1「include=falseの試料はこれらの評価前に除外する」。
                        下記 sample_names と同じ行だけを残す）。
        sample_names  : list[str] — 同上の理由で include=false のサンプル名と、
                        `dataset_exclude` で除外した試料（`ds.excluded_samples`）は
                        含まない。
        feature_names : list[str]（SMF_ID）
        roles         : {sample_name: "sample"|"qc"|"blank"}。
                        **ただし** 明示メタデータ適用済み（`ds.sample_metadata_rows`
                        あり）の場合は、シートが持つ正当な列挙値 "unknown" もここへ
                        素通りする（"sample"/"qc"/"blank" の3値へ絞り込まない）。
                        **この辞書自体は include=false のサンプルも含めた全件分**
                        （数値評価の入力を絞ることと、監査・表示用の写像を保持
                        することは別の判断——`sample_names`/`matrix` には現れない
                        キーが残ることがある）。
        sample_meta   : {sample_name: {role, batch, batch_source, run_order,
                         run_order_source}}
                        — 根拠バンドルの `samples.tsv` の役割と注入順の材料。
                        roles と同じく include=false 分も含めた全件分を返す。
    """
    if ds.feature_matrix is None or not ds.sample_names or not ds.feature_ids:
        raise PreconditionError(
            "missing_state",
            "DatasetState に定量行列がありません。load_dataset を先に実行してください。",
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
                "load_dataset をやり直してから metadata を再適用してください。",
                {"rows": len(explicit_rows), "samples": len(sample_names)},
            )
        roles, sample_meta = sample_manifest.pp_inputs_from_rows(sample_names, explicit_rows)

        # include=false の試料は行列から外す（spec §8.1「include=falseの試料はこれら
        # の評価前に除外する」）。
        # roles/sample_meta の写像自体は削らない——除外理由の監査・表示
        # （dataset_status 等）は全件分の役割を必要とするため。
        keep_idx = included_sample_indices(ds)
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
    use_mztab_batch = sample_manifest.batch_labels_informative(mztab_batches)

    sample_meta: dict = {}
    for i, name in enumerate(sample_names):
        date = sample_manifest.filename_date_batch(name)
        if use_mztab_batch and mztab_batches[i] is not None:
            batch, batch_source = mztab_batches[i], "mztab_batch_label"
        elif date:
            batch, batch_source = date, "filename_date"
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
    # `dataset_exclude` で除外した試料は、明示メタデータの経路と同じく行列と名前から外す
    # （roles / sample_meta の写像は全件分のまま）。
    keep_idx = included_sample_indices(ds)
    if keep_idx is not None:
        matrix = matrix[keep_idx, :]
        sample_names = [sample_names[i] for i in keep_idx]
    return matrix, sample_names, feature_names, roles, sample_meta
