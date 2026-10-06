"""版管理した解析行列（spec §8.2 `analysis-matrix.v1`）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

統計は「直近の前処理結果」ではなく **matrix result ID** を参照する。同じ
データセットから recipe 違いの行列が2本出るのが前提で（未補正 height と
内部標準比、filter の強弱）、`ds.pp_matrix` のような単一スロットへ上書きすると、
図と表が別 recipe の数字を並べてもどちらも「最新」を名乗れてしまう。
だから `make_matrix` は **ds を一切書き換えない**。

## 2段階に分ける理由

`make_matrix`（補完前）と `finalize_matrix`（QC後 filter と補完）を分ける。
補完を先にやると、QC で落ちる予定の feature の値が補完の材料に混ざり、
「落とすはずだったもの」が他のセルへ数値として溶け込む。両者は別の result
（別 ID）で、`parent_id` で繋ぐ。

## 列は消さない

検出率 filter で解析対象から外れた feature も列としては残す。内部標準は解析
対象外でも support として必要で、消すと比の分母を後から再計算できない。
「解析対象か」は `eligibility_mask`、「値が有るか」は値そのもの、という別の軸に
する。

## 保存

npz（`allow_pickle=False`）＋ JSON meta を原子的に書く。pickle を許すと、読み込みが
任意コード実行になる。load 時に軸順・値hash・mask形状を確認し、1つでも食い違えば
拒否する——黙って別の数字を統計へ渡さない。
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import numpy as np

from metabolomix.analysis import internal_standards, preprocessing
from metabolomix.analysis.result_state import array_fingerprint, dataset_fingerprint
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "MATRIX_SCHEMA",
    "finalize_matrix",
    "load_matrix",
    "make_matrix",
    "save_matrix",
]

MATRIX_SCHEMA = "analysis-matrix.v1"

_RECIPE_INVALID = "MATRIX_RECIPE_INVALID"
_INTEGRITY = "MATRIX_INTEGRITY_MISMATCH"
_NOT_FOUND = "MATRIX_NOT_FOUND"

#: `preprocessing.qc_drift_correct` が要求する QC 本数。ここを下回るときは
#: 別手法へ切り替えず停止する（spec §9）。
_MIN_QC_FOR_DRIFT = 4

#: npz に入れる配列。値と mask はここに列挙したものだけで、meta（JSON）には
#: 配列を入れない——2か所に同じ数字を持つと、片方だけ書き換えられる。
_ARRAY_KEYS = ("values", "eligibility_mask", "locked_mask", "imputed_mask",
               "detected_mask")


# ---------- 検出 mask ----------

def _detected_mask(ds, evidence: dict, assay_ids: list[str],
                   feature_ids: list[str]) -> tuple[np.ndarray | None, list[str]]:
    """(assay × feature) の検出 mask を証拠から組む。

    証拠が使えないときは `None` を返し、理由を caveat に残す——「検出状態が
    分からない」を「全部未検出」や「全部検出」に畳まない（前者は検出率を
    語ってはいけない状態で、後者は嘘になる）。
    """
    if (evidence or {}).get("availability"):
        mask = np.ones((len(assay_ids), len(feature_ids)), dtype=bool)
        assay_index = {a: i for i, a in enumerate(assay_ids)}
        feature_index = {f: i for i, f in enumerate(feature_ids)}
        for row in evidence.get("rows") or []:
            i = assay_index.get(row.get("assay_id"))
            j = feature_index.get(row.get("feature_id"))
            if i is None or j is None:
                continue
            mask[i, j] = row.get("detection_status") == "detected"
        return mask, []

    existing = getattr(ds, "detected_mask", None)
    if existing is not None:
        # DatasetState 側は (feature × assay)。解析行列は (assay × feature)。
        return np.asarray(existing, dtype=bool).T, []
    return None, ["detection_unknown"]


# ---------- recipe の適用 ----------

def _check_drift_prerequisites(ds) -> tuple[dict, dict, list[int]]:
    rows = getattr(ds, "sample_metadata_rows", None) or []
    selected = [i for i, r in enumerate(rows) if r.get("include", True)]
    active = [rows[i] for i in selected]
    qc = [r for r in active if r.get("role") == "qc"]
    batches = {r.get("batch") for r in active}
    pools = {r.get("qc_pool") for r in qc}
    orders = [r.get("injection_order") for r in active]
    if (len(qc) < _MIN_QC_FOR_DRIFT or len(batches) != 1 or not all(batches)
            or len(pools) != 1 or not all(pools)
            or any(o is None for o in orders) or len(set(orders)) != len(orders)):
        raise DomainError("QC_PREREQUISITE_MISSING",
                          "ドリフト補正にはincluded QC4本以上、単一の確認済みbatch/pool、全注入の一意なorderが必要です。",
                          {"n_qc": len(qc)})
    roles = {r["sample_id"]: r.get("role") for r in active}
    run_order = {r["sample_id"]: r.get("injection_order") for r in active}
    return roles, run_order, selected


def _apply_recipe(values, recipe: dict, ds, assay_ids: list[str]) -> tuple:
    """normalize と drift_correct を適用し、補正履歴を返す（補完はしない）。"""
    history: list[dict] = []
    caveats: list[str] = []
    base = recipe.get("base")
    method = recipe.get("normalize", "none")

    if base == internal_standards.RATIO_UNIT and method != "none":
        # 内部標準比の上に TIC/median/PQN を重ねると、どの量の比なのかが
        # 言えなくなる（二重正規化）。初期版では拒否する。
        raise DomainError(
            _RECIPE_INVALID,
            f"内部標準比には全体正規化を重ねられません（normalize={method!r}）。"
            "二重正規化を暗黙に実行しません。",
            {"base": base, "normalize": method})

    rows = getattr(ds, "sample_metadata_rows", None) or []
    if rows and len(rows) != len(assay_ids):
        raise DomainError(_RECIPE_INVALID, "metadataとassay軸の件数が一致しません。", {})
    if method != "none":
        if not rows:
            values, factors, report = preprocessing.normalize(values, method)
        else:
            included = [i for i, r in enumerate(rows)
                        if r.get("include", True) and r.get("role") in ("sample", "qc")]
            if not included:
                raise DomainError(_RECIPE_INVALID, "正規化の参照試料がありません。", {})
            # 参照の推定だけを対象試料に限定し、全assay軸は保持する。
            if method == "pqn":
                qc = [i for i in included if rows[i].get("role") == "qc"]
                reference_rows = qc or [i for i in included if rows[i].get("role") == "sample"]
                reference = np.nanmedian(values[reference_rows], axis=0)
                factors = np.nanmedian(values / np.where(reference == 0, np.nan, reference), axis=1)
                report = {"method": method, "pqn_reference": "qc_median" if qc else "sample_median"}
            elif method in ("tic", "median"):
                factors = (np.nansum(values, axis=1) if method == "tic"
                           else np.nanmedian(values, axis=1))
                valid_ref = factors[included]
                valid_ref = valid_ref[np.isfinite(valid_ref) & (valid_ref != 0)]
                reference = np.nanmedian(valid_ref) if valid_ref.size else 1.
                factors = factors / reference
                reference_rows = included
                report = {"method": method}
            else:
                raise DomainError(_RECIPE_INVALID, "未知の正規化です。", {"normalize": method})
            degenerate = ~np.isfinite(factors) | (factors == 0)
            values = values / np.where(degenerate, 1., factors)[:, None]
            report.update(reference_assay_ids=[assay_ids[i] for i in reference_rows],
                          unscaled_samples=int(degenerate.sum()))
            if degenerate.any():
                report["caveat"] = "正規化係数が不正な試料は未正規化のまま保持しました。"
        history.append({"step": "normalize", **report})
        if report.get("caveat"):
            caveats.append(report["caveat"])

    if recipe.get("drift_correct"):
        roles, run_order, selected = _check_drift_prerequisites(ds)
        names = [rows[i]["sample_id"] for i in selected]
        corrected, report = preprocessing.qc_drift_correct(
            values[selected], roles, names, run_order, min_qc=_MIN_QC_FOR_DRIFT)
        if report.get("status") != "applied":
            raise DomainError("QC_PREREQUISITE_MISSING", "ドリフト補正の前提を満たしません。", report)
        values = values.copy()
        values[selected] = corrected
        history.append({"step": "drift_correct", "assay_ids": [assay_ids[i] for i in selected], **report})

    return values, history, caveats


def _detection_eligibility(eligibility, detected, recipe: dict,
                           feature_ids: list[str]) -> tuple[np.ndarray, list[dict]]:
    """検出率 filter を **eligibility にだけ**効かせる（列は消さない）。"""
    spec = recipe.get("filter")
    if not spec:
        return eligibility, []
    min_rate = spec.get("min_detection_rate")
    if min_rate is None:
        return eligibility, []
    if detected is None:
        # 検出状態が不明なまま「検出率で落とす」判断はできない。落とさずに
        # 理由を残す（gap-fill 後の非ゼロで代用しない。spec §8）。
        return eligibility, [{"step": "detection_filter", "status": "not_evaluable",
                              "reason": "detection_unknown"}]
    # 率の定義と閾値判定は preprocessing 側が唯一の実装。ここは (assay × feature)
    # なので、(feature × サンプル) を取る契約へ転置して渡す。
    keep, report = preprocessing.detection_rate_filter(detected.T, min_rate)
    return eligibility & keep, [{
        "step": "detection_filter", "status": "applied", **report,
        "removed_from_eligibility": sorted(
            fid for fid, k in zip(feature_ids, keep) if not k)}]


# ---------- 公開関数 ----------

def _matrix_id(payload: dict) -> str:
    return f"mat_{canonical_hash(payload)[:32]}"


def make_matrix(ds, recipe: dict, bindings: dict, evidence: dict,
                eligibility) -> dict:
    """recipe 1件ぶんの解析行列を作る（補完はしない・ds は書き換えない）。

    `eligibility` は「この時点で解析対象にしてよい feature」の bool 配列。
    検出率 filter はこれを狭めるだけで、値の列は削らない。
    """
    assay_ids = list(getattr(ds, "sample_assay_ids", []) or [])
    feature_ids = list(getattr(ds, "feature_ids", []) or [])
    source = np.asarray(getattr(ds, "feature_matrix", None), dtype=float)
    if source.ndim != 2:
        raise DomainError(
            _RECIPE_INVALID, "定量行列がありません。",
            {"dataset_id": getattr(ds, "dataset_id", None)})
    # DatasetState は (feature × assay)。解析行列は pp_matrix と同じ
    # (assay × feature) に揃える。
    values = source.T.copy()

    detected, caveats = _detected_mask(ds, evidence, assay_ids, feature_ids)
    pairs = list(((bindings or {}).get("internal_standard_map") or {})
                 .get("pairs") or [])

    base = recipe.get("base", "peak_height")
    units = [internal_standards.RAW_UNIT] * len(feature_ids)
    locked = np.zeros(values.shape, dtype=bool)
    missing_reasons: dict[str, dict] = {}
    support: list[str] = []

    if base == internal_standards.RATIO_UNIT:
        applied = internal_standards.apply_internal_standards(
            values, feature_ids, pairs, detected)
        values = applied["values"]
        units = applied["units"]
        locked = applied["locked_mask"]
        support = applied["support_feature_ids"]
        missing_reasons = {
            fid: {assay_ids[int(row)]: reason for row, reason in rows.items()}
            for fid, rows in applied["missing_reasons"].items()
        }
    elif base != "peak_height":
        raise DomainError(
            _RECIPE_INVALID, f"未知のbaseです: {base!r}",
            {"base": base, "supported": ["peak_height",
                                         internal_standards.RATIO_UNIT]})

    values, history, recipe_caveats = _apply_recipe(values, recipe, ds, assay_ids)
    caveats = [*caveats, *recipe_caveats]
    if recipe.get("normalize", "none") != "none":
        units = ["normalized_height"] * len(feature_ids)

    eligibility = np.asarray(eligibility, dtype=bool).copy()
    metadata_rows = getattr(ds, "sample_metadata_rows", None) or []
    filter_detected = detected
    if metadata_rows and detected is not None:
        selected = [i for i, row in enumerate(metadata_rows)
                    if row.get("include", True) and row.get("role") == "sample"]
        filter_detected = detected[selected] if selected else None
    eligibility, filter_history = _detection_eligibility(
        eligibility, filter_detected, recipe, feature_ids)
    history.extend(filter_history)

    identity = {
        "dataset": dataset_fingerprint(ds),
        "recipe": recipe,
        "binding": {"rule_hash": (bindings or {}).get("rule_hash"),
                    "dataset_hash": (bindings or {}).get("dataset_hash"),
                    "pairs": pairs},
        "evidence": {"availability": bool((evidence or {}).get("availability")),
                     "source_artifact_hash": (evidence or {}).get(
                         "source_artifact_hash")},
        "metadata": {"assay_ids": assay_ids, "feature_ids": feature_ids,
                     "units": units,
                     "metadata_revision": getattr(ds, "metadata_revision", 0),
                     "rows": metadata_rows},
        "eligibility": eligibility.tolist(),
        "stage": "preprocessed",
    }

    return {
        "schema": MATRIX_SCHEMA,
        "matrix_id": _matrix_id(identity),
        "parent_id": None,
        "stage": "preprocessed",
        "dataset_id": getattr(ds, "dataset_id", None),
        "recipe_id": recipe.get("recipe_id"),
        "recipe": dict(recipe),
        "recipe_hash": canonical_hash(recipe),
        "binding_rule_hash": (bindings or {}).get("rule_hash"),
        "assay_ids": assay_ids,
        "feature_ids": feature_ids,
        "values": values,
        "units": units,
        "base": base,
        "eligibility_mask": eligibility,
        "detected_mask": detected,
        "imputed_mask": np.zeros(values.shape, dtype=bool),
        "locked_mask": locked,
        "missing_reasons": missing_reasons,
        "support_feature_ids": support,
        "correction_history": history,
        "caveats": caveats,
    }


def finalize_matrix(matrix: dict, eligibility, impute: str) -> dict:
    """QC 後 filter を適用し、必要なら補完して**別の result** を作る。

    補完するのは eligibility の立っている列だけ。対象外の列まで埋めると、
    落としたはずの feature が「値のある列」として残り、次の工程が拾う。

    `locked`（分母不正）のセルは補完後も欠損に戻す。ここを埋めると「分母が
    無かった」という事実が消え、平均・分散に他試料の値が混ざる。
    """
    values = np.array(matrix["values"], dtype=float, copy=True)
    locked = np.asarray(matrix["locked_mask"], dtype=bool)
    eligibility = (np.asarray(matrix["eligibility_mask"], dtype=bool)
                   & np.asarray(eligibility, dtype=bool))

    imputed = np.zeros(values.shape, dtype=bool)
    history = list(matrix.get("correction_history") or [])
    if impute and impute != "none" and eligibility.any():
        columns = np.flatnonzero(eligibility)
        subset = values[:, columns]
        filled, report = preprocessing.impute(subset, method=impute)
        was_missing = np.isnan(subset)
        values[:, columns] = filled
        imputed[:, columns] = was_missing & ~np.isnan(filled)
        history.append({"step": "impute", **report,
                        "columns": [matrix["feature_ids"][c] for c in columns]})

    # 補完後に locked を復元する。補完器は locked を知らない（知らせると
    # 補完器ごとに同じ規則を書く羽目になる）ので、ここで一度だけ戻す。
    values[locked] = np.nan
    imputed[locked] = False

    identity = {
        "parent": matrix["matrix_id"],
        "eligibility": eligibility.tolist(),
        "impute": impute,
        "stage": "finalized",
    }
    result = dict(matrix)
    result.update({
        "matrix_id": _matrix_id(identity),
        "parent_id": matrix["matrix_id"],
        "stage": "finalized",
        "values": values,
        "eligibility_mask": eligibility,
        "imputed_mask": imputed,
        "locked_mask": locked,
        "correction_history": history,
    })
    return result


# ---------- 保存と読み戻し ----------

def _atomic_write_bytes(path: Path, payload: bytes) -> None:
    """同じ親ディレクトリの一時ファイルへ書き、fsync してから置換する。

    `core/atomic_io.py` は stdlib だけの leaf（別プロセスの worker が最初に
    import する）なので、numpy を必要とする npz 保存はここに置く。
    """
    parent = path.parent
    parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
                dir=parent, prefix=f".{path.name}.", suffix=".tmp",
                delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
            tmp_file.write(payload)
            tmp_file.flush()
            os.fsync(tmp_file.fileno())
        os.replace(tmp_path, path)
        tmp_path = None
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink(missing_ok=True)


def _array_hashes(matrix: dict) -> dict:
    return {key: array_fingerprint(matrix.get(key)) for key in _ARRAY_KEYS}


def save_matrix(result: dict, directory: Path) -> dict:
    """行列を npz + JSON meta で保存し、load に必要な参照を返す。"""
    import io

    directory = Path(directory)
    matrix_id = result["matrix_id"]
    arrays_file = f"{matrix_id}.npz"
    meta_file = f"{matrix_id}.json"

    arrays = {key: np.asarray(result[key]) for key in _ARRAY_KEYS
              if result.get(key) is not None}
    buffer = io.BytesIO()
    np.savez(buffer, **arrays)
    _atomic_write_bytes(directory / arrays_file, buffer.getvalue())

    meta = {key: result[key] for key in (
        "schema", "matrix_id", "parent_id", "stage", "dataset_id", "recipe_id",
        "recipe", "recipe_hash", "binding_rule_hash", "assay_ids", "feature_ids",
        "units", "base", "missing_reasons", "support_feature_ids",
        "correction_history", "caveats")}
    meta["array_hashes"] = _array_hashes(result)
    meta["arrays_file"] = arrays_file
    _atomic_write_bytes(
        directory / meta_file,
        json.dumps(meta, ensure_ascii=False, separators=(",", ":"),
                   allow_nan=False).encode("utf-8"))

    return {"matrix_id": matrix_id, "meta_file": meta_file,
            "arrays_file": arrays_file,
            "meta_hash": canonical_hash(meta)}


def load_matrix(reference: dict, directory: Path) -> dict:
    """保存した行列を読み、軸順・値hash・mask形状を確認してから返す。"""
    directory = Path(directory)
    meta_path = directory / reference["meta_file"]
    arrays_path = directory / reference["arrays_file"]
    if not meta_path.is_file() or not arrays_path.is_file():
        raise DomainError(
            _NOT_FOUND,
            "保存済みの解析行列が見つかりません。",
            {"matrix_id": reference.get("matrix_id"),
             "meta_path": str(meta_path), "arrays_path": str(arrays_path)})

    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta.get("matrix_id") != reference.get("matrix_id"):
        raise DomainError(
            _INTEGRITY, "matrix_id が参照と一致しません。",
            {"expected": reference.get("matrix_id"), "found": meta.get("matrix_id")})
    if canonical_hash(meta) != reference.get("meta_hash"):
        # 軸順（assay_ids / feature_ids）の書き換えもここで捕まる。列順が
        # 変わった行列を通すと、全 feature のラベルが黙って1つずれる。
        raise DomainError(
            _INTEGRITY, "meta が保存時から変更されています（軸順・単位を含む）。",
            {"matrix_id": reference.get("matrix_id")})

    with np.load(arrays_path, allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}

    expected = meta["array_hashes"]
    for key in _ARRAY_KEYS:
        actual = array_fingerprint(arrays.get(key))
        if actual != expected.get(key):
            raise DomainError(
                _INTEGRITY,
                f"保存した配列が変更されています: {key}",
                {"matrix_id": reference.get("matrix_id"), "array": key})

    n_assays, n_features = len(meta["assay_ids"]), len(meta["feature_ids"])
    if arrays["values"].shape != (n_assays, n_features):
        raise DomainError(
            _INTEGRITY, "値の形が軸の長さと一致しません。",
            {"shape": list(arrays["values"].shape),
             "axes": [n_assays, n_features]})

    result = dict(meta)
    result.update(arrays)
    result.setdefault("detected_mask", None)
    return result
