"""定量の工程が作る派生データセットと、その新しいハンドルへの登録（Task 10）。

正規化・統合・選別は元の DatasetState を書き換えず、複製に値・検出マスク・欄を差し替えた
派生の DatasetState を作り、新しいハンドルとして開いて直近にする。元のハンドルは残る。

`quant_state` は統計側（C）との受け渡し。C は `getattr(ds, "quant_state", None)` を読み、
`normalized` / `method` / `unit` だけを使う。それ以外のキーは自由（`**extra`）。
"""
from __future__ import annotations

import copy
import uuid

import numpy as np

from metabolomix.core import handles, session_state
from metabolomix.mztab.polarity import dataset_polarity

__all__ = ["QUANT_STATE_SCHEMA", "check_capacity", "derive", "quant_state", "register"]

QUANT_STATE_SCHEMA = "lipidmix.quant-state.v1"


def quant_state(*, normalized: bool, method: str | None, unit: str | None, record_id: str,
                steps: list[str], skipped: list[str], **extra) -> dict:
    state = {"schema": QUANT_STATE_SCHEMA, "normalized": bool(normalized), "method": method,
             "unit": unit, "record_id": record_id, "steps": list(steps),
             "skipped": list(skipped)}
    state.update(extra)
    return state


def _sample_positions(parent_names: list[str], sample_names, sample_index) -> list[int]:
    """新しい各試料の親での位置。index があればそれ、無ければ名前で引く。"""
    if sample_index is not None:
        idx = [int(i) for i in sample_index]
        bad = [i for i in idx if not 0 <= i < len(parent_names)]
        if bad:
            raise ValueError(f"sample_index が親の範囲外: {bad}")
        return idx
    if sample_names is None:
        return list(range(len(parent_names)))
    position = {name: i for i, name in enumerate(parent_names)}
    missing = [n for n in sample_names if n not in position]
    if missing:
        raise ValueError(f"親に無い試料名（renamed なら sample_index を渡す）: {missing}")
    return [position[n] for n in sample_names]


def derive(ds, *, matrix: np.ndarray, detected: np.ndarray,
           feature_ids: list[str] | None = None, feature_metadata: dict | None = None,
           sample_names: list[str] | None = None, sample_index: list[int] | None = None,
           state: dict):
    """`ds` の複製に値・マスク・欄を差し替えた派生データセット（`ds` は動かさない）。

    古い行列から出た結果（results / analysis_matrices / pp_* / last_* / preprocess_* /
    preprocessing_recipe / handoff_report）は空に戻す。指紋は保存せず都度計算する
    （`result_state.dataset_fingerprint`）ので、古い指紋が残ることはない。

    試料: `sample_index` は新しい各試料の親での位置。試料名だけを渡せば名前で引く（親に無い
    名前は ValueError）。index を渡したときの `sample_names` は新しい名前を決めるだけ（統合は
    試料を照合した ID へ改名するので、親に無い名前もありうる）。`sample_names` に並走する
    リスト（sample_assay_ids / sample_ids / sample_raw_paths / msdial_samples）は位置で、
    名前で引く欄（roles / sample_meta / per_sample_mz の内側）は新しい名前へ、親の試料から
    引き継ぐ。

    特徴: `feature_ids` を渡すと、ID で引く欄（feature_metadata（明示しないとき）/
    feature_candidates / feature_annotations / per_sample_mz / feature_tags）をその ID に絞る。
    親の母集団の要約（sml_rows / sme_rows は None、inchikey_coverage は空、feature_qc は
    セル数・検出数・欠測率だけ新しい行列から数え直す）は持ち越さない。`curation_flags` は
    アラインメント単位のフラグ集合なので触らない。**親に無い新しい ID（統合後の
    `positive:<id>` など）には副次データを付けない。呼び出し側が埋める。**
    """
    new = copy.deepcopy(ds)
    parent_names = list(ds.sample_names)
    idx = _sample_positions(parent_names, sample_names, sample_index)
    if len(set(idx)) != len(idx):
        raise ValueError(f"同じ親の試料を 2 回選べない（重複する位置）: {idx}")
    new_names = list(sample_names) if sample_names is not None else [parent_names[i] for i in idx]
    if len(new_names) != len(idx):
        raise ValueError(f"sample_names（{len(new_names)}）と sample_index（{len(idx)}）の長さが違う")
    if feature_ids is not None:
        new.feature_ids = list(feature_ids)
    if feature_metadata is not None:
        new.feature_metadata = copy.deepcopy(feature_metadata)
    new.sample_names = new_names
    matrix = np.array(matrix, dtype=float)
    detected = np.array(detected, dtype=bool)
    shape = (len(new.feature_ids), len(new_names))
    if matrix.shape != shape or detected.shape != shape:
        raise ValueError(f"行列の形が合わない: matrix={matrix.shape} detected={detected.shape} "
                         f"期待={shape}（特徴 × 試料）")

    n_parent = len(parent_names)
    for attr in ("sample_assay_ids", "sample_ids", "sample_raw_paths", "msdial_samples"):
        values = getattr(ds, attr, None)
        if values is not None and len(values) == n_parent:
            setattr(new, attr, [copy.deepcopy(values[i]) for i in idx])
    if isinstance(new.msdial_samples, list) and len(new.msdial_samples) == len(new_names):
        for entry, name in zip(new.msdial_samples, new_names):
            if isinstance(entry, dict):
                entry["name"] = name
    for attr in ("roles", "sample_meta"):
        old = getattr(ds, attr) or {}
        setattr(new, attr, {name: copy.deepcopy(old[parent_names[i]])
                            for name, i in zip(new_names, idx) if parent_names[i] in old})
    if ds.per_sample_mz is not None:
        new.per_sample_mz = {
            fid: {name: inner[parent_names[i]] for name, i in zip(new_names, idx)
                  if parent_names[i] in inner}
            for fid, inner in ds.per_sample_mz.items()}

    if feature_ids is not None:
        keep = set(new.feature_ids)
        if feature_metadata is None:
            new.feature_metadata = {k: v for k, v in new.feature_metadata.items() if k in keep}
        for attr in ("feature_candidates", "feature_annotations", "feature_tags"):
            setattr(new, attr, {k: v for k, v in getattr(new, attr).items() if k in keep})
        if new.per_sample_mz is not None:
            new.per_sample_mz = {k: v for k, v in new.per_sample_mz.items() if k in keep}
        new.sml_rows = None
        new.sme_rows = None
        new.inchikey_coverage = {}
    if new.feature_qc and "n_cells" in new.feature_qc:
        n_cells, n_det = int(detected.size), int(detected.sum())
        new.feature_qc = {**new.feature_qc, "n_cells": n_cells, "n_detected": n_det,
                          "gap_filled_rate": round(1.0 - n_det / n_cells, 4) if n_cells else None}

    new.feature_matrix = matrix
    new.detected_mask = detected
    new.dataset_id = f"ds_{uuid.uuid4().hex}"
    new.results = {}
    new.analysis_matrices = {}
    new.pp_matrix = None
    new.pp_sample_names = []
    new.pp_feature_names = []
    new.last_pca = None
    new.last_differential = None
    new.preprocess_id = None
    new.preprocess_metadata_hash = None
    new.preprocessing_recipe = {}
    new.handoff_report = None
    new.quant_state = copy.deepcopy(state)
    return new


def check_capacity() -> None:
    """記録を書く**前に**、新しいハンドルを 1 本足せるかを確かめる（上限なら
    `DATASET_HANDLE_LIMIT`）。書いてから登録で断ると、どのハンドルにも繋がらない記録が残る。
    上限の意味は登録簿のまま（`register` も同じ検査をもう一度する）。"""
    handles.registry.check_capacity()


def register(holder, ds, *, step: str, parent_handle_ids: list[str],
             polarity: str | None = None) -> str:
    """派生データセットを新しいハンドルとして登録し、直近にして handle_id を返す。

    上限なら `DATASET_HANDLE_LIMIT`（古いハンドルは閉じず、直近も動かさない）。ライブラリは
    新しい LibraryState（親の store を 2 か所で閉じないため。必要なら呼び出し側が
    `library_load` し直す）。分光の参照先（.arf2 / .dcl）は親が 1 本のときだけ引き継ぐ。
    `polarity` は呼び出し側が決めた極性（引数・ハンドル・データセットの順に解いたもの）。
    `polarity` を省くと、データセット由来 → 親のハンドルの極性の順。
    """
    single = len(parent_handle_ids) == 1
    if len(parent_handle_ids) > 1:
        polarity = None  # 親が複数（統合）なら極性は決まらない。親のハンドルの極性へ落とさない
    elif polarity is None:
        polarity, _source, _notes = dataset_polarity(ds)
        if polarity is None:
            polarity = getattr(holder, "polarity", None)
    handle = handles.DatasetHandle(
        ds, polarity=polarity, polarity_source="quant",
        source={"kind": "quant", "step": step, "parents": list(parent_handle_ids)},
        arf2_path=getattr(holder, "arf2_path", None) if single else None,
        dcl_path=getattr(holder, "dcl_path", None) if single else None)
    handles.registry.add(handle)
    session_state.session.set_current(handle.handle_id)
    return handle.handle_id
