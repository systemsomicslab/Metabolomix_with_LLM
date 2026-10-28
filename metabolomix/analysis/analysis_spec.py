"""analysis-spec.v1: 論文の方法（行列の作り方・検定・補正・閾値）を宣言した spec を実行する
（spec 2026-10-22 §4 の 3・4・5・10・11・12）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない（DatasetState は引数で
受ける）。行列は `recipe_matrix`、検定は `spec_statistics` が唯一の実装で、ここは「どの
試料・どの群・どの特徴に、どの順で当てるか」を spec から決めるだけ。

**判定は返さない。** 返すのは数字（統計量・p・q・FC・群平均）と向きと、spec に書いた閾値を
機械的に当てた `meets_thresholds` だけ。再現したかどうかは再解析の側が決める。

**同じ spec・同じ入力なら同じ結果。** 時刻や乱数の ID を入れない。`rule_hash` は正規化した
spec からパスの欄を除いて作り、入力のファイルは `input_fingerprint` に sha256 で入れる。
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

from metabolomix.analysis import alignment_qc, recipe_matrix, spec_statistics
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["ANALYSIS_SCHEMA", "RESULT_SCHEMA", "normalize_spec", "rule_hash", "run_spec"]

ANALYSIS_SCHEMA = "analysis-spec.v1"
RESULT_SCHEMA = "analysis-result.v1"
_INVALID = "ANALYSIS_SPEC_INVALID"
_PATH_KEYS = frozenset({"path", "matrix_path", "alignment_qc"})
_SLUG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")
_MISSING_TOKENS = frozenset({"", "na", "nan", "null", "none"})


def _fail(message: str, **details) -> DomainError:
    return DomainError(_INVALID, message, details)


def _exact_keys(value, allowed: set, required: set, label: str) -> dict:
    if not isinstance(value, dict):
        raise _fail(f"{label} は辞書で指定する", field=label)
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise _fail(f"{label} の未知のキー: {unknown}", field=label, unknown=unknown)
    lacking = sorted(required - set(value))
    if lacking:
        raise _fail(f"{label} に必須のキーが無い: {lacking}", field=label, missing=lacking)
    return value


def _string(value, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _fail(f"{label} は空でない文字列で指定する", field=label)
    return value.strip()


def _choice(value, choices, label: str):
    if value not in choices:
        raise _fail(f"{label} は {' / '.join(map(str, choices))}: {value!r}", field=label)
    return value


def _optional_positive(value, label: str):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
        raise _fail(f"{label} は正の数か null: {value!r}", field=label)
    return float(value)


# ---------- 正規化 ----------

def _normalize_tsv(raw, label: str, columns: dict[str, str | None]) -> dict:
    """{"source": "tsv", "path": ..., <列名の欄>...}。columns は欄 → 既定。"""
    raw = _exact_keys(raw, {"source", "path", *columns}, {"source", "path"}, label)
    out = {"source": "tsv", "path": _string(raw["path"], f"{label}.path")}
    for key, default in columns.items():
        value = raw.get(key, default)
        out[key] = None if value is None else _string(value, f"{label}.{key}")
    return out


def _normalize_input(raw) -> dict:
    raw = raw or {"source": "dataset"}
    source = raw.get("source") if isinstance(raw, dict) else None
    if source == "dataset":
        _exact_keys(raw, {"source"}, {"source"}, "input")
        return {"source": "dataset"}
    if source != "tsv":
        raise _fail("input.source は dataset / tsv", field="input.source")
    raw = _exact_keys(raw, {"source", "matrix_path", "orientation", "id_column"},
                      {"source", "matrix_path"}, "input")
    orientation = _choice(raw.get("orientation", "features_by_samples"),
                          ("features_by_samples", "samples_by_features"), "input.orientation")
    default_id = "feature_id" if orientation == "features_by_samples" else "sample"
    return {"source": "tsv", "matrix_path": _string(raw["matrix_path"], "input.matrix_path"),
            "orientation": orientation,
            "id_column": _string(raw.get("id_column", default_id), "input.id_column")}


def _normalize_groups(raw, input_source: str) -> dict:
    raw = raw or ({"source": "metadata"} if input_source == "dataset" else None)
    if raw is None:
        raise _fail("外部の行列（input.source=tsv）には groups の TSV が要る", field="groups")
    if raw.get("source") == "metadata":
        if input_source != "dataset":
            raise _fail("groups.source=metadata はデータセットの入力でだけ使える", field="groups")
        _exact_keys(raw, {"source"}, {"source"}, "groups")
        return {"source": "metadata"}
    if raw.get("source") != "tsv":
        raise _fail("groups.source は metadata / tsv", field="groups.source")
    return _normalize_tsv(raw, "groups", {"sample_column": "sample", "group_column": "group",
                                          "role_column": None})


def _normalize_order(raw, input_source: str):
    if raw is None:
        return None
    if raw.get("source") == "metadata":
        if input_source != "dataset":
            raise _fail("injection_order.source=metadata はデータセットの入力でだけ使える",
                        field="injection_order")
        _exact_keys(raw, {"source"}, {"source"}, "injection_order")
        return {"source": "metadata"}
    if raw.get("source") != "tsv":
        raise _fail("injection_order.source は metadata / tsv", field="injection_order.source")
    return _normalize_tsv(raw, "injection_order", {"sample_column": "sample",
                                                   "order_column": "order"})


def _normalize_covariates(raw):
    if raw is None:
        return None
    raw = _exact_keys(raw, {"path", "sample_column", "columns", "missing"},
                      {"path", "columns"}, "covariates")
    columns = raw["columns"]
    if not isinstance(columns, list) or not columns:
        raise _fail("covariates.columns は 1 つ以上の列の定義", field="covariates.columns")
    out_columns, seen = [], set()
    for i, column in enumerate(columns):
        label = f"covariates.columns[{i}]"
        kind = column.get("type") if isinstance(column, dict) else None
        if kind == "numeric":
            column = _exact_keys(column, {"name", "type"}, {"name", "type"}, label)
            item = {"name": _string(column["name"], f"{label}.name"), "type": "numeric"}
        elif kind == "categorical":
            column = _exact_keys(column, {"name", "type", "reference"},
                                 {"name", "type", "reference"}, label)
            item = {"name": _string(column["name"], f"{label}.name"), "type": "categorical",
                    "reference": _string(column["reference"], f"{label}.reference")}
        else:
            raise _fail(f"{label}.type は numeric / categorical", field=f"{label}.type")
        if item["name"] in seen:
            raise _fail(f"covariates の列が重複している: {item['name']}", field=label)
        seen.add(item["name"])
        out_columns.append(item)
    return {"path": _string(raw["path"], "covariates.path"),
            "sample_column": _string(raw.get("sample_column", "sample"), "covariates.sample_column"),
            "missing": _choice(raw.get("missing", "drop"), ("drop", "error"), "covariates.missing"),
            "columns": out_columns}


def _normalize_features(raw) -> dict:
    """報告する特徴（`all` か `feature_ids`）と、割れたスポットのまとめ方。

    `merge_split`（既定 true）は、割れた 2 スポットを試料ごとの大きい方の値で 1 つにして
    から検定する（`alignment_qc` の `merge_groups` だけを対象にする）。`alignment_qc` は
    `dataset_alignment_qc` が書いた JSON のパスで、無ければその場で `find_split_pairs` を回す。
    """
    raw = raw if raw is not None else {"all": True}
    if not isinstance(raw, dict):
        raise _fail('features は {"all": true} か {"feature_ids": [...]}（と merge_split・'
                    "alignment_qc）", field="features")
    unknown = sorted(set(raw) - {"all", "feature_ids", "merge_split", "alignment_qc"})
    if unknown:
        raise _fail(f"features の未知のキー: {unknown}", field="features", unknown=unknown)
    if "all" in raw and "feature_ids" in raw:
        raise _fail('features は {"all": true} か {"feature_ids": [...]} のどちらか一方',
                    field="features")
    merge = raw.get("merge_split", True)
    if not isinstance(merge, bool):
        raise _fail("features.merge_split は true / false", field="features.merge_split")
    if "alignment_qc" in raw and not merge:
        raise _fail("features.alignment_qc は merge_split=true のときだけ使える（効かない設定を"
                    "書かない）", field="features.alignment_qc")
    if "feature_ids" in raw:
        ids = raw["feature_ids"]
        if not isinstance(ids, list) or not ids or not all(isinstance(i, str) and i for i in ids):
            raise _fail("features.feature_ids は空でない文字列のリスト", field="features.feature_ids")
        if len(set(ids)) != len(ids):
            raise _fail("features.feature_ids に重複がある", field="features.feature_ids")
        out = {"feature_ids": list(ids)}
    else:
        if raw.get("all", True) is not True:
            raise _fail('features は {"all": true} か {"feature_ids": [...]}', field="features")
        out = {"all": True}
    out["merge_split"] = merge
    if "alignment_qc" in raw:
        out["alignment_qc"] = _string(raw["alignment_qc"], "features.alignment_qc")
    return out


def _normalize_aggregate(raw, input_source: str):
    if raw is None:
        return None
    kind = raw.get("kind") if isinstance(raw, dict) else None
    if kind == "ratio":
        raw = _exact_keys(raw, {"kind", "label", "numerator", "denominator"},
                          {"kind", "label", "numerator", "denominator"}, "aggregate")
        out = {"kind": "ratio", "label": _string(raw["label"], "aggregate.label"),
               "numerator": _string(raw["numerator"], "aggregate.numerator"),
               "denominator": _string(raw["denominator"], "aggregate.denominator")}
        if out["numerator"] == out["denominator"]:
            raise _fail("aggregate の分子と分母が同じ特徴", field="aggregate")
        return out
    if kind != "sum":
        raise _fail("aggregate.kind は sum / ratio", field="aggregate.kind")
    raw = _exact_keys(raw, {"kind", "label", "feature_ids", "select"}, {"kind", "label"},
                      "aggregate")
    if ("feature_ids" in raw) == ("select" in raw):
        raise _fail("aggregate（sum）は feature_ids か select のどちらか一方", field="aggregate")
    out = {"kind": "sum", "label": _string(raw["label"], "aggregate.label")}
    if "feature_ids" in raw:
        ids = raw["feature_ids"]
        if not isinstance(ids, list) or not ids or not all(isinstance(i, str) and i for i in ids):
            raise _fail("aggregate.feature_ids は空でない文字列のリスト", field="aggregate")
        out["feature_ids"] = sorted(set(ids), key=ids.index)
        return out
    if input_source != "dataset":
        raise _fail("aggregate.select は同定を持つデータセットの入力でだけ使える", field="aggregate")
    select = _exact_keys(raw["select"], {"name_regex", "evidence"}, {"name_regex"},
                         "aggregate.select")
    try:
        re.compile(select["name_regex"])
    except (re.error, TypeError) as exc:
        raise _fail(f"aggregate.select.name_regex を読めない: {exc}", field="aggregate.select")
    evidence = select.get("evidence", ["sme"])
    if (not isinstance(evidence, list) or not evidence
            or not set(evidence) <= {"sme", "sml"}):
        raise _fail("aggregate.select.evidence は sme / sml のリスト", field="aggregate.select")
    out["select"] = {"name_regex": select["name_regex"],
                     "evidence": [e for e in ("sme", "sml") if e in evidence]}
    return out


def _normalize_analysis(raw, index: int, covariate_names: set) -> dict:
    label = f"analyses[{index}]"
    raw = _exact_keys(raw, {"analysis_id", "matrix", "test", "fold_change", "correction",
                            "bh_population", "thresholds"}, {"analysis_id", "test"}, label)
    analysis_id = _string(raw["analysis_id"], f"{label}.analysis_id")
    if not _SLUG_RE.match(analysis_id):
        raise _fail(f"{label}.analysis_id は英数字と _ . - だけ: {analysis_id!r}", field=label)
    test = _exact_keys(raw["test"], {"kind", "transform", "covariates"}, {"kind"}, f"{label}.test")
    kind = _choice(test["kind"], spec_statistics.KINDS, f"{label}.test.kind")
    used = test.get("covariates", [])
    if not isinstance(used, list) or not all(isinstance(c, str) for c in used):
        raise _fail(f"{label}.test.covariates は列名のリスト", field=label)
    if kind != "ols" and used:
        raise _fail(f"{label}: 共変量は test.kind=ols でだけ使える", field=label)
    unknown = sorted(set(used) - covariate_names)
    if unknown:
        raise _fail(f"{label}: covariates.columns に無い共変量 {unknown}", field=label)
    thresholds = _exact_keys(raw.get("thresholds", {}), {"p", "q", "fc"}, set(),
                             f"{label}.thresholds")
    correction = _choice(raw.get("correction", "bh"), ("bh", "none"), f"{label}.correction")
    return {
        "analysis_id": analysis_id,
        "matrix": recipe_matrix.normalize_recipe(raw.get("matrix")),
        "test": {"kind": kind,
                 "transform": _choice(test.get("transform", "log2"), spec_statistics.TRANSFORMS,
                                      f"{label}.test.transform"),
                 "covariates": list(used)},
        "fold_change": _choice(raw.get("fold_change", "arithmetic_mean_ratio"),
                               spec_statistics.FC_DEFINITIONS, f"{label}.fold_change"),
        "correction": correction,
        "bh_population": _choice(raw.get("bh_population", "matrix"),
                                 ("matrix", "reported", "annotated"), f"{label}.bh_population"),
        "thresholds": {"p": _optional_positive(thresholds.get("p"), f"{label}.thresholds.p"),
                       "q": _optional_positive(thresholds.get("q"), f"{label}.thresholds.q"),
                       "fc": _optional_positive(thresholds.get("fc"), f"{label}.thresholds.fc")},
    }


def normalize_spec(spec) -> dict:
    """spec を検証して既定を埋める。未知のキーは拒否する。"""
    spec = _exact_keys(spec, {"schema", "input", "groups", "comparison", "injection_order",
                              "covariates", "features", "aggregate", "analyses"},
                       {"comparison", "analyses"}, "spec")
    if spec.get("schema", ANALYSIS_SCHEMA) != ANALYSIS_SCHEMA:
        raise _fail(f"schema は {ANALYSIS_SCHEMA}", field="schema")
    input_ = _normalize_input(spec.get("input"))
    comparison = _exact_keys(spec["comparison"], {"reference", "test"}, {"reference", "test"},
                             "comparison")
    reference = _string(comparison["reference"], "comparison.reference")
    test = _string(comparison["test"], "comparison.test")
    if reference == test:
        raise _fail("comparison の reference と test が同じ群", field="comparison")
    covariates = _normalize_covariates(spec.get("covariates"))
    names = {c["name"] for c in (covariates or {}).get("columns", [])}
    analyses = spec["analyses"]
    if not isinstance(analyses, list) or not analyses:
        raise _fail("analyses は 1 つ以上の解析のリスト", field="analyses")
    normalized_analyses = [_normalize_analysis(a, i, names) for i, a in enumerate(analyses)]
    ids = [a["analysis_id"] for a in normalized_analyses]
    if len(set(ids)) != len(ids):
        raise _fail(f"analysis_id が重複している: {ids}", field="analyses")
    aggregate = _normalize_aggregate(spec.get("aggregate"), input_["source"])
    features = _normalize_features(spec.get("features"))
    if aggregate is not None and "feature_ids" in features:
        raise _fail("aggregate を使うときは features.feature_ids を書かない（報告するのは集約した 1 本）",
                    field="features")
    if input_["source"] != "dataset":
        # 外部の行列にはアライメントが無い: merge_split は何もしない（既定の true のまま通す）。
        if "alignment_qc" in features:
            raise _fail("features.alignment_qc はデータセットの入力でだけ使える", field="features")
        if any(a["bh_population"] == "annotated" for a in normalized_analyses):
            raise _fail("bh_population=annotated（SME のある特徴）はデータセットの入力でだけ使える",
                        field="analyses")
    return {
        "schema": ANALYSIS_SCHEMA,
        "input": input_,
        "groups": _normalize_groups(spec.get("groups"), input_["source"]),
        "comparison": {"reference": reference, "test": test},
        "injection_order": _normalize_order(spec.get("injection_order"), input_["source"]),
        "covariates": covariates,
        "features": features,
        "aggregate": aggregate,
        "analyses": normalized_analyses,
    }


def _strip_paths(value):
    if isinstance(value, dict):
        return {k: _strip_paths(v) for k, v in value.items() if k not in _PATH_KEYS}
    if isinstance(value, list):
        return [_strip_paths(v) for v in value]
    return value


def merge_groups_hash(groups) -> str:
    """まとめた割れの群（特徴 ID の並びの並び）の hash。順によらない。"""
    return canonical_hash(sorted(sorted(str(f) for f in group) for group in groups))


def rule_hash(normalized: dict, merged_groups=None) -> str:
    """正規化した spec の hash（パスの欄を除く——置き場所が違っても同じ規則なら同じ）。

    データセットの入力で `merge_split` が真なら、実際にまとめた割れの群の hash も入れる
    （`merged_groups`。同じ spec でも割れの組が変われば別の規則）。外部の行列にはアライメントが
    無く `merge_split` は何もしないので、hash にも入れない。
    """
    material = _strip_paths(normalized)
    if normalized["input"]["source"] != "dataset":
        material["features"] = {k: v for k, v in material["features"].items()
                                if k != "merge_split"}
    elif normalized["features"]["merge_split"] and merged_groups is not None:
        material["merged_groups_hash"] = merge_groups_hash(merged_groups)
    return canonical_hash(material)


# ---------- 入力 ----------

def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path_text: str, label: str) -> tuple[list[str], list[dict], str]:
    path = Path(path_text)
    if not path.is_file():
        raise DomainError("ANALYSIS_INPUT_NOT_FOUND", f"{label} の TSV が見つからない: {path}",
                          {"field": label, "path": str(path)})
    with open(path, encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    rows = list(reader)
    return list(reader.fieldnames or []), rows, _sha256(path)


def _require_column(columns: list[str], column: str, label: str) -> None:
    if column not in columns:
        raise DomainError("ANALYSIS_INPUT_INVALID", f"{label} に列 {column!r} が無い（列: {columns}）",
                          {"field": label, "column": column, "columns": columns})


def _keyed(rows: list[dict], column: str, label: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    duplicated = []
    for row in rows:
        key = (row.get(column) or "").strip()
        if key in out:
            duplicated.append(key)
        out[key] = row
    if duplicated:
        raise DomainError("SAMPLE_MISMATCH", f"{label} に同じ試料が複数ある: {sorted(set(duplicated))}",
                          {"field": label, "duplicated": sorted(set(duplicated))})
    return out


def _float_or_nan(text) -> float:
    text = (text or "").strip()
    if text.lower() in _MISSING_TOKENS:
        return math.nan
    return float(text)


def _external_matrix(spec_input: dict) -> tuple[np.ndarray, list[str], list[str], str]:
    columns, rows, digest = _read_tsv(spec_input["matrix_path"], "input")
    id_column = spec_input["id_column"]
    _require_column(columns, id_column, "input")
    others = [c for c in columns if c != id_column]
    # csv.DictReader は同じ名前の列の最後の値しか残さない。重複した列名を通すと、
    # 1 つの試料（特徴）が 2 回数えられて片方の値が黙って消える。
    duplicated = sorted({c for c in others if others.count(c) > 1})
    if duplicated or not all(c.strip() for c in others):
        raise DomainError("ANALYSIS_INPUT_INVALID",
                          f"input の列名が重複しているか空です（重複: {duplicated}）",
                          {"field": "input", "duplicated": duplicated})
    ids = [(row.get(id_column) or "").strip() for row in rows]
    if not all(ids):
        raise DomainError("ANALYSIS_INPUT_INVALID", "input に ID が空の行がある",
                          {"field": "input"})
    if len(set(ids)) != len(ids):
        raise DomainError("ANALYSIS_INPUT_INVALID", "input の行の ID が重複している",
                          {"field": "input",
                           "duplicated": sorted({i for i in ids if ids.count(i) > 1})})
    try:
        table = np.array([[_float_or_nan(row.get(c)) for c in others] for row in rows],
                         dtype=float)
    except ValueError as exc:
        raise DomainError("ANALYSIS_INPUT_INVALID", f"input に数として読めない値がある: {exc}",
                          {"field": "input"}) from exc
    if spec_input["orientation"] == "features_by_samples":
        return table.T, others, ids, digest
    return table, ids, others, digest


def _check_cover(wanted: list[str], available, label: str, *, allow_extra: bool) -> None:
    lacking = [s for s in wanted if s not in available]
    extra = [] if allow_extra else sorted(set(available) - set(wanted))
    if lacking or extra:
        raise DomainError(
            "SAMPLE_MISMATCH",
            f"{label} の試料が行列と合わない（無い: {lacking}、行列に無い: {extra}）。"
            "試料名は完全一致で照合する。",
            {"field": label, "missing": lacking, "extra": extra})


def _load(normalized: dict, ds) -> dict:
    """行列・試料の役割と群・注入順・共変量を揃える。"""
    files: dict[str, str] = {}
    detected = None
    metadata_rows = None
    if normalized["input"]["source"] == "dataset":
        if ds is None:
            raise DomainError("DATASET_REQUIRED", "input.source=dataset にはデータセットが要る", {})
        rows = list(getattr(ds, "sample_metadata_rows", None) or [])
        names = list(getattr(ds, "sample_ids", None) or []) if rows else []
        names = names or list(ds.sample_names)
        values = np.asarray(ds.feature_matrix, dtype=float).T
        feature_ids = list(ds.feature_ids)
        mask = getattr(ds, "detected_mask", None)
        detected = None if mask is None else np.asarray(mask, dtype=bool).T
        if rows:
            include = [i for i, row in enumerate(rows) if row.get("include", True)]
            values, names = values[include], [names[i] for i in include]
            detected = None if detected is None else detected[include]
            metadata_rows = {names[k]: rows[i] for k, i in enumerate(include)}
    else:
        values, names, feature_ids, files["matrix"] = _external_matrix(normalized["input"])

    roles: dict[str, str] = {}
    groups: dict[str, str | None] = {}
    spec_groups = normalized["groups"]
    if spec_groups["source"] == "metadata":
        if metadata_rows is None:
            raise DomainError("GROUPS_UNAVAILABLE",
                              "このデータセットには実験情報シート（群）が無い。groups を TSV で与えるか、"
                              "pipeline_path で読む（試料対応表ごと載る）。", {})
        for name in names:
            roles[name] = metadata_rows[name].get("role") or "sample"
            groups[name] = metadata_rows[name].get("group") or None
    else:
        columns, rows_, files["groups"] = _read_tsv(spec_groups["path"], "groups")
        for key in ("sample_column", "group_column", "role_column"):
            if spec_groups[key] is not None:
                _require_column(columns, spec_groups[key], "groups")
        keyed = _keyed(rows_, spec_groups["sample_column"], "groups")
        _check_cover(names, keyed, "groups", allow_extra=False)
        for name in names:
            row = keyed[name]
            groups[name] = (row.get(spec_groups["group_column"]) or "").strip() or None
            if spec_groups["role_column"] is not None:
                roles[name] = (row.get(spec_groups["role_column"]) or "").strip() or "sample"
            elif metadata_rows is not None:
                roles[name] = metadata_rows[name].get("role") or "sample"
            else:
                roles[name] = "sample"

    order = None
    spec_order = normalized["injection_order"]
    if spec_order is not None:
        if spec_order["source"] == "metadata":
            order = {n: (metadata_rows or {}).get(n, {}).get("injection_order") for n in names}
        else:
            columns, rows_, files["injection_order"] = _read_tsv(spec_order["path"],
                                                                 "injection_order")
            _require_column(columns, spec_order["sample_column"], "injection_order")
            _require_column(columns, spec_order["order_column"], "injection_order")
            keyed = _keyed(rows_, spec_order["sample_column"], "injection_order")
            _check_cover(names, keyed, "injection_order", allow_extra=True)
            try:
                order = {n: float(keyed[n][spec_order["order_column"]]) for n in names}
            except ValueError as exc:
                raise DomainError("ANALYSIS_INPUT_INVALID", f"注入順を数として読めない: {exc}",
                                  {"field": "injection_order"}) from exc

    covariates = None
    if normalized["covariates"] is not None:
        spec_cov = normalized["covariates"]
        columns, rows_, files["covariates"] = _read_tsv(spec_cov["path"], "covariates")
        _require_column(columns, spec_cov["sample_column"], "covariates")
        for column in spec_cov["columns"]:
            _require_column(columns, column["name"], "covariates")
        covariates = _keyed(rows_, spec_cov["sample_column"], "covariates")

    return {"values": values, "samples": names, "feature_ids": feature_ids, "roles": roles,
            "groups": groups, "detected": detected, "order": order, "covariates": covariates,
            "files": files}


# ---------- 割れたスポットのまとめ ----------

def _read_alignment_qc(path_text: str) -> tuple[dict, str]:
    path = Path(path_text)
    if not path.is_file():
        raise DomainError("ANALYSIS_INPUT_NOT_FOUND",
                          f"features.alignment_qc の JSON が見つからない: {path}",
                          {"field": "features.alignment_qc", "path": str(path)})
    try:
        qc = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        raise DomainError("ANALYSIS_INPUT_INVALID",
                          f"features.alignment_qc を JSON として読めない: {exc}",
                          {"field": "features.alignment_qc"}) from exc
    groups = qc.get("merge_groups") if isinstance(qc, dict) else None
    if (not isinstance(qc, dict) or qc.get("schema") != alignment_qc.SCHEMA
            or not isinstance(groups, list)):
        raise DomainError("ANALYSIS_INPUT_INVALID",
                          f"features.alignment_qc は {alignment_qc.SCHEMA}（dataset_alignment_qc の "
                          "output_path の JSON）で、merge_groups を持つこと",
                          {"field": "features.alignment_qc"})
    return qc, _sha256(path)


def _merge_plan(normalized: dict, ds, data: dict) -> dict:
    """まとめる割れの群を決める（`merge_groups` の mergeable な組だけ。まとめの判断は
    `alignment_qc` が持ち、ここでは足さない・広げない）。

    外部の行列にはアライメントが無いので何もしない。`groups` は、行列に 2 つ以上の特徴が
    ある群（実際に使う群）。
    """
    plan = {"enabled": False, "source": None, "groups": [], "groups_hash": None,
            "n_groups_given": 0}
    features = normalized["features"]
    if normalized["input"]["source"] != "dataset" or not features["merge_split"]:
        return plan
    path = features.get("alignment_qc")
    if path is None:
        qc, source = alignment_qc.find_split_pairs(ds), "find_split_pairs"
    else:
        qc, data["files"]["alignment_qc"] = _read_alignment_qc(path)
        source = "alignment_qc"
    given = [[str(f) for f in group["features"]] for group in qc.get("merge_groups") or []
             if group.get("mergeable")]
    present = {str(f) for f in data["feature_ids"]}
    used = sorted(sorted(group) for group in given if sum(f in present for f in group) >= 2)
    seen: set[str] = set()
    for group in used:
        if seen.intersection(group):
            raise DomainError("ANALYSIS_INPUT_INVALID",
                              f"割れの群が特徴を共有している: {sorted(seen.intersection(group))}",
                              {"field": "features.alignment_qc"})
        seen.update(group)
    return {"enabled": True, "source": source, "groups": used,
            "groups_hash": merge_groups_hash(used), "n_groups_given": len(given)}


def _fold_members(members: list[str], rename: dict[str, str]) -> list[str]:
    """まとめて消えた特徴を、残った ID に置き換える（重複は 1 つ）。"""
    out: list[str] = []
    for fid in members:
        fid = rename.get(fid, fid)
        if fid not in out:
            out.append(fid)
    return out


# ---------- 集約・共変量 ----------

def _selected_by_name(ds, select: dict) -> list[str]:
    pattern = re.compile(select["name_regex"])
    chosen = []
    for fid in ds.feature_ids:
        names = []
        if "sme" in select["evidence"]:
            candidates = (ds.feature_candidates or {}).get(fid) or []
            if candidates:
                names.append(candidates[0].get("chemical_name"))
        if "sml" in select["evidence"]:
            names.append(((ds.feature_annotations or {}).get(fid) or {}).get("name"))
        if any(isinstance(n, str) and pattern.search(n) for n in names):
            chosen.append(fid)
    return chosen


def _aggregate(matrix: dict, aggregate: dict, members: list[str]) -> tuple[np.ndarray, dict]:
    index = {fid: j for j, fid in enumerate(matrix["feature_ids"])}
    present = [fid for fid in members if fid in index]
    dropped = [fid for fid in members if fid not in index]
    if aggregate["kind"] == "sum":
        if not present:
            column = np.full(len(matrix["sample_ids"]), np.nan)
        else:
            block = matrix["values"][:, [index[f] for f in present]]
            column = np.where(np.all(np.isfinite(block), axis=1), block.sum(axis=1), np.nan)
        return column, {"members": members, "summed": present, "dropped": dropped}
    if dropped:
        column = np.full(len(matrix["sample_ids"]), np.nan)
    else:
        num = matrix["values"][:, index[aggregate["numerator"]]]
        den = matrix["values"][:, index[aggregate["denominator"]]]
        with np.errstate(all="ignore"):
            column = np.where(np.isfinite(num) & np.isfinite(den) & (den > 0), num / den, np.nan)
    return column, {"numerator": aggregate["numerator"],
                    "denominator": aggregate["denominator"], "dropped": dropped}


def _design(samples: list[str], used: list[str], spec_cov: dict | None, table) -> tuple:
    """共変量の計画行列（符号化済み）と、欠測で除いた試料。"""
    if not used:
        return np.zeros((len(samples), 0)), []
    by_name = {c["name"]: c for c in spec_cov["columns"]}
    _check_cover(samples, table, "covariates", allow_extra=True)
    columns, excluded = [], set()
    for name in used:
        definition = by_name[name]
        raw = [(table[s].get(name) or "").strip() for s in samples]
        missing = [s for s, v in zip(samples, raw) if v.lower() in _MISSING_TOKENS]
        if missing and spec_cov["missing"] == "error":
            raise DomainError("COVARIATE_MISSING", f"共変量 {name} が欠測の試料: {missing}",
                              {"covariate": name, "samples": missing})
        excluded.update(missing)
        if definition["type"] == "numeric":
            try:
                columns.append([math.nan if s in missing else float(v)
                                for s, v in zip(samples, raw)])
            except ValueError as exc:
                raise DomainError("ANALYSIS_INPUT_INVALID",
                                  f"共変量 {name} を数として読めない: {exc}",
                                  {"covariate": name}) from exc
        else:
            levels = sorted({v for s, v in zip(samples, raw) if s not in missing})
            if definition["reference"] not in levels:
                raise DomainError("ANALYSIS_INPUT_INVALID",
                                  f"共変量 {name} の基準の水準 {definition['reference']!r} が"
                                  f"データに無い（水準: {levels}）",
                                  {"covariate": name, "levels": levels})
            for level in levels:
                if level == definition["reference"]:
                    continue
                columns.append([math.nan if s in missing else float(v == level)
                                for s, v in zip(samples, raw)])
    return np.array(columns, dtype=float).T, sorted(excluded)


# ---------- 実行 ----------

def _direction(row: dict) -> str | None:
    if row["status"] != "tested":
        return None
    effect = row["beta_log2"] if row["beta_log2"] is not None else row["log2fc_test_over_ref"]
    if effect is None or effect == 0:
        return None
    return "test_higher" if effect > 0 else "ref_higher"


def _meets(row: dict, thresholds: dict):
    if row["status"] != "tested" or all(v is None for v in thresholds.values()):
        return None
    checks = []
    if thresholds["p"] is not None:
        checks.append(row["p"] is not None and row["p"] < thresholds["p"])
    if thresholds["q"] is not None:
        checks.append(row["q"] is not None and row["q"] < thresholds["q"])
    if thresholds["fc"] is not None:
        log2fc = row["log2fc_test_over_ref"]
        checks.append(log2fc is not None and abs(log2fc) >= math.log2(thresholds["fc"]))
    return all(checks)


def _run_analysis(analysis: dict, data: dict, normalized: dict, ds, reported: list[str] | None,
                  merge: dict):
    matrix = recipe_matrix.build_matrix(
        data["values"], data["samples"], data["feature_ids"], data["roles"], analysis["matrix"],
        detected=data["detected"], injection_order=data["order"])
    # 割れた 2 スポットを 1 列にする（recipe の後・集約と検定の前）。recipe の篩で片方が
    # 落ちた群は列に 2 つ残らないので、まとめない（`skipped`）。
    merged_from: dict[str, list[str]] = {}
    rename: dict[str, str] = {}
    merge_info = None
    if merge["enabled"]:
        merged = alignment_qc.merge_matrix_columns(matrix["values"], matrix["feature_ids"],
                                                   merge["groups"])
        merged_from = merged.merged_from
        rename = {f: kept for kept, group in merged_from.items() for f in group if f != kept}
        matrix = {**matrix, "values": merged.matrix, "feature_ids": merged.feature_ids}
        merged_ids = {f for group in merged_from.values() for f in group}
        merge_info = {"merged_from": merged_from,
                      "skipped": [g for g in merge["groups"] if not merged_ids.intersection(g)]}
    origin = {fid: merged_from.get(fid, [fid]) for fid in matrix["feature_ids"]}
    aggregate_info = None
    aggregate = normalized["aggregate"]
    if aggregate is not None:
        members = (aggregate.get("feature_ids") or
                   ([aggregate["numerator"], aggregate["denominator"]] if aggregate["kind"] == "ratio"
                    else _selected_by_name(ds, aggregate["select"])))
        folded = _fold_members(members, rename)
        if aggregate["kind"] == "ratio":
            if len(folded) < 2:
                raise _fail("aggregate の分子と分母が、まとめた後に同じ特徴になる",
                            field="aggregate", members=members)
            aggregate = {**aggregate, "numerator": folded[0], "denominator": folded[1]}
        column, aggregate_info = _aggregate(matrix, aggregate, folded)
        if merge_info is not None:
            # 選んだ特徴のうち割れの相手を持つものは、まとめた 1 つとして 1 回だけ数える。
            aggregate_info["members"] = list(members)
            aggregate_info["merged_from"] = {f: merged_from[f] for f in folded if f in merged_from}
        used = (aggregate_info["summed"] if aggregate["kind"] == "sum"
                else [aggregate["numerator"], aggregate["denominator"]])
        origin = {aggregate["label"]: [o for f in used if f in origin for o in origin[f]]}
        matrix = {**matrix, "values": column[:, None], "feature_ids": [aggregate["label"]]}
        merged_from = {}
        reported = [aggregate["label"]]
    elif reported is not None:
        reported = _fold_members(reported, rename)

    reference, test = normalized["comparison"]["reference"], normalized["comparison"]["test"]
    index = {s: i for i, s in enumerate(matrix["sample_ids"])}
    members = {label: [s for s in matrix["sample_ids"]
                       if data["roles"].get(s, "sample") == "sample" and data["groups"].get(s) == label]
               for label in (reference, test)}
    too_small = {label: len(m) for label, m in members.items() if len(m) < 2}
    if too_small:
        raise DomainError("STATISTIC_GROUP_TOO_SMALL",
                          f"n<2 の群がある: {too_small}（群の名前・role・include を確かめる）",
                          {"group_n": {k: len(v) for k, v in members.items()}})
    rows_ref = [index[s] for s in members[reference]]
    rows_test = [index[s] for s in members[test]]

    kind = analysis["test"]["kind"]
    transform = analysis["test"]["transform"]
    design, excluded = None, []
    if kind == "ols":
        samples = members[reference] + members[test]
        design, excluded = _design(samples, analysis["test"]["covariates"], normalized["covariates"],
                                   data["covariates"])

    results = {}
    for j, fid in enumerate(matrix["feature_ids"]):
        ref_raw = matrix["values"][rows_ref, j]
        test_raw = matrix["values"][rows_test, j]
        if kind == "ols":
            y = np.concatenate([ref_raw, test_raw])
            keep = np.array([s not in excluded for s in members[reference] + members[test]])
            group = np.array([0.0] * len(rows_ref) + [1.0] * len(rows_test))
            outcome = spec_statistics.ols_group(np.where(keep, y, np.nan), group, design,
                                                transform=transform)
        else:
            outcome = spec_statistics.two_group(ref_raw, test_raw, kind=kind, transform=transform)
            outcome["beta_log2"] = None
        fc = spec_statistics.fold_change(ref_raw, test_raw, analysis["fold_change"])
        finite_ref, finite_test = ref_raw[np.isfinite(ref_raw)], test_raw[np.isfinite(test_raw)]
        results[fid] = {
            "feature_id": fid, "status": outcome["status"], "reason": outcome["reason"],
            "group_ref": reference, "group_test": test,
            "n_ref": outcome["n_ref"], "n_test": outcome["n_test"],
            "mean_ref": float(finite_ref.mean()) if finite_ref.size else None,
            "mean_test": float(finite_test.mean()) if finite_test.size else None,
            "fc_test_over_ref": fc,
            "log2fc_test_over_ref": math.log2(fc) if fc else None,
            "statistic": outcome["statistic"], "p": outcome["p"], "q": None,
            "beta_log2": outcome["beta_log2"],
        }

    wanted = matrix["feature_ids"] if reported is None else reported
    if analysis["bh_population"] == "matrix":
        population = matrix["feature_ids"]
    elif analysis["bh_population"] == "annotated":
        # SME（`ds.feature_candidates`）のある特徴だけ。まとめた特徴は、まとめた元のどちらかに
        # SME があれば入る。母集団の外の特徴の q は null。
        candidates = getattr(ds, "feature_candidates", None) or {}
        population = [f for f in matrix["feature_ids"]
                      if any(candidates.get(o) for o in origin.get(f, [f]))]
    else:
        population = [f for f in wanted if f in results]
    if analysis["correction"] == "bh":
        qs = spec_statistics.bh([results[f]["p"] for f in population])
        for fid, q in zip(population, qs):
            results[fid]["q"] = q
    bh_n = sum(1 for f in population if results[f]["p"] is not None)

    rows = []
    for fid in wanted:
        row = results.get(fid)
        if row is None:
            row = {"feature_id": fid, "status": "not_in_matrix", "reason": "removed_by_matrix_recipe",
                   "group_ref": reference, "group_test": test, "n_ref": None, "n_test": None,
                   "mean_ref": None, "mean_test": None, "fc_test_over_ref": None,
                   "log2fc_test_over_ref": None, "statistic": None, "p": None, "q": None,
                   "beta_log2": None}
        row = dict(row)
        row["merged_from"] = merged_from.get(fid)
        row["direction"] = _direction(row)
        row["meets_thresholds"] = _meets(row, analysis["thresholds"])
        rows.append(row)
    out = {
        "analysis_id": analysis["analysis_id"],
        "matrix": {"n_features_in": matrix["n_features_in"],
                   "n_features_out": matrix["n_features_out"], "steps": matrix["steps"]},
        "test": analysis["test"],
        "bh_population_n": bh_n if analysis["correction"] == "bh" else None,
        "thresholds": analysis["thresholds"],
        "rows": rows,
    }
    if kind == "ols":
        out["covariate_excluded_samples"] = excluded
    if aggregate_info is not None:
        out["aggregate"] = aggregate_info
    if merge_info is not None:
        out["merge_split"] = merge_info
    return out


def run_spec(normalized: dict, ds=None) -> dict:
    """正規化済みの spec を実行して analysis-result.v1 を返す（丸めない）。"""
    data = _load(normalized, ds)
    merge = _merge_plan(normalized, ds, data)
    reported = normalized["features"].get("feature_ids")
    if reported is not None:
        unknown = [f for f in reported if f not in set(data["feature_ids"])]
        if unknown:
            raise DomainError("FEATURE_NOT_FOUND", f"行列に無い feature_id: {unknown}",
                              {"feature_ids": unknown})
    aggregate = normalized["aggregate"]
    if aggregate is not None and "select" not in aggregate:
        members = aggregate.get("feature_ids") or [aggregate["numerator"], aggregate["denominator"]]
        unknown = [f for f in members if f not in set(data["feature_ids"])]
        if unknown:
            raise DomainError("FEATURE_NOT_FOUND", f"行列に無い feature_id: {unknown}",
                              {"feature_ids": unknown})
    if aggregate is not None and "select" in aggregate and ds is None:
        raise DomainError("DATASET_REQUIRED", "aggregate.select にはデータセットが要る", {})

    from metabolomix.analysis.result_state import dataset_fingerprint

    analyses = [_run_analysis(a, data, normalized, ds, reported, merge)
                for a in normalized["analyses"]]
    reference, test = normalized["comparison"]["reference"], normalized["comparison"]["test"]
    return {
        "schema": RESULT_SCHEMA,
        "spec": copy.deepcopy(normalized),
        "rule_hash": rule_hash(normalized, merge["groups"] if merge["enabled"] else None),
        "merge_split": merge,
        "input_fingerprint": {
            "dataset_hash": dataset_fingerprint(ds) if normalized["input"]["source"] == "dataset"
            else None,
            "files": dict(sorted(data["files"].items())),
        },
        "comparison": f"{test}_vs_{reference}",
        "group_ref": reference,
        "group_test": test,
        "analyses": analyses,
    }
