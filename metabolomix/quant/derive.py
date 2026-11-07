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

__all__ = ["QUANT_STATE_SCHEMA", "derive", "quant_state", "register"]

QUANT_STATE_SCHEMA = "lipidmix.quant-state.v1"


def quant_state(*, normalized: bool, method: str | None, unit: str | None, record_id: str,
                steps: list[str], skipped: list[str], **extra) -> dict:
    state = {"schema": QUANT_STATE_SCHEMA, "normalized": bool(normalized), "method": method,
             "unit": unit, "record_id": record_id, "steps": list(steps),
             "skipped": list(skipped)}
    state.update(extra)
    return state


def derive(ds, *, matrix: np.ndarray, detected: np.ndarray,
           feature_ids: list[str] | None = None, feature_metadata: dict | None = None,
           sample_names: list[str] | None = None, state: dict):
    """`ds` の複製に値・マスク・欄を差し替えた派生データセット（`ds` は動かさない）。

    古い行列から出た結果（results / analysis_matrices / pp_* / last_* / preprocess_* /
    preprocessing_recipe / handoff_report）は空に戻す。指紋は保存せず都度計算する
    （`result_state.dataset_fingerprint`）ので、古い指紋が残ることはない。
    """
    new = copy.deepcopy(ds)
    if feature_ids is not None:
        new.feature_ids = list(feature_ids)
    if feature_metadata is not None:
        new.feature_metadata = copy.deepcopy(feature_metadata)
    if sample_names is not None:
        new.sample_names = list(sample_names)
    matrix, detected = np.asarray(matrix), np.asarray(detected)
    shape = (len(new.feature_ids), len(new.sample_names))
    if matrix.shape != shape or detected.shape != shape:
        raise ValueError(f"行列の形が合わない: matrix={matrix.shape} detected={detected.shape} "
                         f"期待={shape}（特徴 × 試料）")
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
    new.quant_state = state
    return new


def register(holder, ds, *, step: str, parent_handle_ids: list[str]) -> str:
    """派生データセットを新しいハンドルとして登録し、直近にして handle_id を返す。

    上限なら `DATASET_HANDLE_LIMIT`（古いハンドルは閉じず、直近も動かさない）。ライブラリは
    新しい LibraryState（親の store を 2 か所で閉じないため。必要なら呼び出し側が
    `library_load` し直す）。分光の参照先（.arf2 / .dcl）は親が 1 本のときだけ引き継ぐ。
    """
    polarity, _source, _notes = dataset_polarity(ds)
    if polarity is None:
        polarity = getattr(holder, "polarity", None)
    single = len(parent_handle_ids) == 1
    handle = handles.DatasetHandle(
        ds, polarity=polarity, polarity_source="quant",
        source={"kind": "quant", "step": step, "parents": list(parent_handle_ids)},
        arf2_path=getattr(holder, "arf2_path", None) if single else None,
        dcl_path=getattr(holder, "dcl_path", None) if single else None)
    handles.registry.add(handle)
    session = session_state.session
    with session._current_lock:
        session.current_handle_id = handle.handle_id
    return handle.handle_id
