"""pipeline-request.v3 の厳密な検証・解決・更新契約（spec 2026-11-06 §5・plan「要求」）。

pipeline は「生データ → MS-DIAL Console → 検証済みの解析フォルダ」までを受け持つ。要求が
語るのは上流の条件だけで、下流（前処理・比較・PCA）は対話の鎖（`load_dataset` →
`dataset_preprocess` / `dataset_compare` / `dataset_pca`）が受け持つ。

**形は 2 つ受け、保存は 1 つ（v3 の正規形）にする。**

- 1 極性の形（今までの形）: `polarity` `method_file` `lbm_file` `msp_file` `sample_manifest`
  をトップに置く。正規化で `units: [{"raw_dir": ".", …}]` の 1 単位になる。
- 2 極性の形: `units: [{"polarity", "raw_dir", "method_file"?, "lbm_file"?, "msp_file"?,
  "sample_manifest"?, "keep_extension"?}]`（1〜2 件。`polarity` と `raw_dir` は必須、
  `raw_dir` は `dataset_root` 相対、極性もフォルダも重ならない）。`units` があるときに
  トップへ単位のキーを置くと `PIPELINE_REQUEST_INVALID`。

共有のキー（`library_mode` `measure` `keep_extension` `timeout_s` `save_project`
`output_root` `alignment_reference`）はトップに置き、意味は今までのまま。

**未指定と明示 null の区別**は今までどおり。`sample_manifest` だけが明示 null を受け、
それ以外の省略可能なキーに明示 null を渡すと `PIPELINE_REQUEST_INVALID`。出所は
`value_sources`（共有のキーはキーごと、単位のキーは `value_sources["units"][i]` に
キーごと。`explicit` / `request_file` / `default` / `explicit_update` /
`method_declaration`）。

**下流のキー**（`target` `preprocess` `comparisons`）は消えた。明示した要求にあれば
`PIPELINE_DOWNSTREAM_REMOVED`（`details.replacements` に代わりのツール）で止め、
`analysis-request.json` にあれば読み飛ばして警告 `REQUEST_FILE_DOWNSTREAM_IGNORED` を
`intake_warnings`（内部キー）に残す（判断 5: 置いてある古いファイルで Console まで
止めるのは過剰）。

解決順序は「明示値 > `analysis-request.json` > 既定値」。1 極性の形の単位のキーは
キーごとに重ねる（今までと同じ）。`units` は層ごとに丸ごと採る（2 極性の単位をキーごとに
混ぜると、どの極性に効く値かが曖昧になるため）。

`schema` は省略・`"pipeline-request.v1"`・`"pipeline-request.v3"` を受け、保存するのは
v3。`"pipeline-request.v2"`（実データ 0 本で廃止）は `PIPELINE_REQUEST_INVALID`
（`details.removed_schema`）。

resume で変えられるのは試料の対応表（`sample_manifest`）だけ（1 極性の run はトップの
`sample_manifest`、2 極性は `units: {<極性>: {"sample_manifest": …}}`）。上流の条件を
変えるなら新しい run（`NEW_PIPELINE_REQUIRED`）。
"""
from __future__ import annotations

import copy
import json
from pathlib import Path, PurePosixPath, PureWindowsPath

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "DOWNSTREAM_KEYS",
    "DOWNSTREAM_REPLACEMENTS",
    "REQUEST_FILE_NAME",
    "SCHEMA",
    "UPDATABLE",
    "merge_updates",
    "pin_polarity",
    "read_request_file",
    "request_fingerprint",
    "resolve_request",
    "unit_index",
    "unit_request",
    "unit_source_root",
    "validate_request",
]

SCHEMA = "pipeline-request.v3"

#: 受ける要求の版（省略は v1 と同じ扱い）。保存するのは常に `SCHEMA`。
_ACCEPTED_SCHEMAS = frozenset({"pipeline-request.v1", SCHEMA})

#: 廃止した要求の版（実データ 0 本。spec 2026-11-06 §7）。
_REMOVED_SCHEMAS = frozenset({"pipeline-request.v2"})

#: 元フォルダ直下から探す任意の要求ファイル（spec §7.1）。
REQUEST_FILE_NAME = "analysis-request.json"

#: resume で変えられるトップのキー（試料の対応表だけ）。
UPDATABLE = frozenset({"sample_manifest", "units"})

#: 消えた下流のキー（判断 5）。
DOWNSTREAM_KEYS = frozenset({"target", "preprocess", "comparisons"})

#: 下流のキーの代わり（対話の鎖。`details.replacements` に載せる順）。
DOWNSTREAM_REPLACEMENTS = ("load_dataset", "dataset_preprocess", "dataset_compare",
                           "dataset_pca")

#: 共有のキー（単位をまたいで 1 つ）。
_SHARED_KEYS = frozenset({
    "library_mode", "measure", "keep_extension", "timeout_s", "save_project",
    "output_root", "alignment_reference",
})

#: 1 極性の形でトップに置く単位のキー。
_UNIT_KEYS = frozenset({"polarity", "method_file", "lbm_file", "msp_file", "sample_manifest"})

#: `units[i]` が受けるキー。
_UNIT_ENTRY_KEYS = _UNIT_KEYS | {"raw_dir", "keep_extension"}
_UNIT_REQUIRED_KEYS = frozenset({"polarity", "raw_dir"})

#: 要求（明示・ファイル）が受けるトップのキー。
_TOP_LEVEL_KEYS = frozenset({"schema", "units"}) | _SHARED_KEYS | _UNIT_KEYS

#: 内部でだけ使うキー。外部入力には現れてはいけない。
_INTERNAL_KEYS = frozenset({"value_sources", "intake_warnings"})

#: 保存する正規形のトップのキー。
_NORMALIZED_KEYS = frozenset({"schema", "units"}) | _SHARED_KEYS

#: `library_mode="msp_then_lbm"`（v8）は研究室 MSP を優先し、`[library] lbm` の LBM で補う。
_LIBRARY_MODE_VALUES = frozenset({"auto", "msp_only", "msp_then_lbm"})
#: MS-DIAL の `Alignment reference file ID`（v8）。
_ALIGNMENT_REFERENCE_VALUES = frozenset({"median_injection", "method"})
_POLARITY_VALUES = frozenset({"positive", "negative"})
_MEASURE_VALUES = frozenset({"peak_height"})

#: 内容 hash から外す既定値（既定値のキーが増えただけで hash を動かさない）。
_FINGERPRINT_OMITTED_DEFAULTS = {"library_mode": "auto",
                                 "alignment_reference": "median_injection"}

_DEFAULT_TIMEOUT_S = 21600  # 既存 Console の既定に合わせる

_SHARED_DEFAULTS = {
    "library_mode": "auto", "measure": "peak_height", "keep_extension": None,
    "timeout_s": _DEFAULT_TIMEOUT_S, "save_project": True, "output_root": None,
    "alignment_reference": "median_injection",
}

#: 明示 null を拒むキー（null は「無効化」を意味する設定だけの語彙。`sample_manifest` だけが受ける）。
_NULL_REJECTED_KEYS = frozenset({
    "method_file", "lbm_file", "polarity", "keep_extension", "output_root",
    "library_mode", "msp_file", "alignment_reference", "units", "raw_dir",
})

_MAX_UNITS = 2


def _fail(message: str, **details) -> None:
    raise DomainError("PIPELINE_REQUEST_INVALID", message, details)


def _downstream_removed(keys, *, where: str) -> DomainError:
    keys = sorted(keys)
    return DomainError(
        "PIPELINE_DOWNSTREAM_REMOVED",
        f"{where}の {keys} は pipeline の下流の工程として消えました。pipeline は Console で解析フォルダを"
        "作るところまでを受け持ち、前処理・比較・PCA は load_dataset(pipeline_path=…, polarity=…) の後に"
        " dataset_preprocess / dataset_compare / dataset_pca で行います。",
        {"removed_keys": keys, "replacements": list(DOWNSTREAM_REPLACEMENTS)})


def _reject_disallowed_explicit_null(source: dict, *, where: str = "") -> None:
    """`source` に**ある**キーのうち、明示 null を許さないものを拒む（未指定は素通り）。"""
    for key in _NULL_REJECTED_KEYS:
        if key in source and source[key] is None:
            _fail(f"{where}{key}に明示nullは許可されません（無効化を意味しない設定のため、"
                  "未指定にしてください）。", **{key: None})


def _is_positive_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def _validate_optional_str(value, field: str) -> None:
    if value is not None and (not isinstance(value, str) or value == ""):
        _fail(f"{field}は空でない文字列またはnullである必要があります: {value!r}", **{field: value})


def _check_schema(schema) -> None:
    if schema is None or schema in _ACCEPTED_SCHEMAS:
        return
    if schema in _REMOVED_SCHEMAS:
        _fail(f"{schema!r} は廃止しました（{SCHEMA!r} の形で送ってください）。",
              removed_schema=schema)
    _fail(f"schemaは{sorted(_ACCEPTED_SCHEMAS)}のいずれかか省略です: {schema!r}", schema=schema)


def _normalize_raw_dir(value, index: int) -> str:
    """`raw_dir` を `dataset_root` 相対の正規の文字列（`/` 区切り）にする。外へ出る指定は拒む。"""
    if not isinstance(value, str) or not value.strip():
        _fail(f"units[{index}].raw_dir は空でない文字列である必要があります: {value!r}",
              index=index, raw_dir=value)
    text = value.strip().replace("\\", "/")
    if PureWindowsPath(value).is_absolute() or PurePosixPath(text).is_absolute() \
            or PureWindowsPath(value).drive:
        _fail(f"units[{index}].raw_dir は dataset_root 相対で指定してください: {value!r}",
              index=index, raw_dir=value)
    parts = [part for part in PurePosixPath(text).parts if part not in ("", ".")]
    if any(part == ".." for part in parts):
        _fail(f"units[{index}].raw_dir が dataset_root の外を指しています: {value!r}",
              index=index, raw_dir=value)
    return "/".join(parts) if parts else "."


def _unit_entry(source: dict, raw_dir: str) -> dict:
    return {
        "polarity": source.get("polarity"),
        "raw_dir": raw_dir,
        "method_file": source.get("method_file"),
        "lbm_file": source.get("lbm_file"),
        "msp_file": source.get("msp_file"),
        "sample_manifest": source.get("sample_manifest"),
        "keep_extension": source.get("keep_extension"),
    }


def _units_from_list(units, *, origin: str) -> tuple[list, list]:
    """2 極性の形の `units` を検査して (単位の一覧, 単位ごとの出所) を返す。"""
    if not isinstance(units, list) or not 1 <= len(units) <= _MAX_UNITS:
        _fail(f"unitsは1〜{_MAX_UNITS}件の配列である必要があります。", value=units)
    normalized, sources = [], []
    for index, entry in enumerate(units):
        if not isinstance(entry, dict):
            _fail(f"units[{index}]はオブジェクトである必要があります。", index=index)
        downstream = set(entry) & DOWNSTREAM_KEYS
        if downstream:
            raise _downstream_removed(downstream, where=f"units[{index}]")
        unknown = set(entry) - _UNIT_ENTRY_KEYS
        if unknown:
            _fail(f"units[{index}]に未知のキーがあります: {sorted(unknown)}",
                  index=index, unknown_keys=sorted(unknown))
        missing = _UNIT_REQUIRED_KEYS - set(entry)
        if missing:
            _fail(f"units[{index}]に必須のキーが足りません: {sorted(missing)}",
                  index=index, missing_keys=sorted(missing))
        _reject_disallowed_explicit_null(entry, where=f"units[{index}]の")
        normalized.append(_unit_entry(entry, _normalize_raw_dir(entry["raw_dir"], index)))
        sources.append({key: (origin if key in entry else "default")
                        for key in sorted(_UNIT_ENTRY_KEYS)})
    return normalized, sources


def _validate_units(data: dict) -> None:
    units = data.get("units")
    if not isinstance(units, list) or not 1 <= len(units) <= _MAX_UNITS:
        _fail(f"unitsは1〜{_MAX_UNITS}件の配列である必要があります。", value=units)
    library_mode = data.get("library_mode", "auto")
    seen_polarity: dict = {}
    seen_raw_dir: dict = {}
    for index, unit in enumerate(units):
        if not isinstance(unit, dict) or set(unit) != _UNIT_ENTRY_KEYS:
            _fail(f"units[{index}]の形が正規形ではありません。", index=index)
        polarity = unit["polarity"]
        if polarity is not None and polarity not in _POLARITY_VALUES:
            _fail(f"polarityが不正です: {polarity!r}", polarity=polarity, index=index)
        if polarity is None and len(units) > 1:
            _fail(f"units[{index}].polarity は必須です（2 単位の要求）。", index=index)
        for field in ("method_file", "lbm_file", "msp_file", "sample_manifest", "keep_extension"):
            _validate_optional_str(unit[field], field)
        if library_mode == "msp_only" and unit["lbm_file"]:
            _fail("library_mode='msp_only'はLBMを使わないので、lbm_fileと同時には指定できません。",
                  library_mode=library_mode)
        if library_mode == "msp_then_lbm" and unit["lbm_file"]:
            _fail("library_mode='msp_then_lbm'のLBMは設定ファイルの[library] lbmだけから決めるので、"
                  "lbm_fileと同時には指定できません。", library_mode=library_mode)
        if library_mode == "auto" and unit["msp_file"]:
            _fail("msp_fileはlibrary_mode='msp_only'か'msp_then_lbm'のときだけ使えます"
                  "（autoではMSPを解決しない）。", library_mode=library_mode)
        if polarity is not None:
            if polarity in seen_polarity:
                _fail(f"units の極性が重なっています: {polarity!r}",
                      polarity=polarity, indexes=[seen_polarity[polarity], index])
            seen_polarity[polarity] = index
        # 同じフォルダだけでなく、片方がもう片方の中にある（`.` は何でも含む）のも拒む
        # ——外側の単位の Console が内側の単位の raw まで読み、同じ raw を 2 度数える。
        parts = tuple(p for p in unit["raw_dir"].lower().split("/") if p not in ("", "."))
        for other_parts, (other_index, other_raw_dir) in seen_raw_dir.items():
            shorter = min(len(parts), len(other_parts))
            if parts[:shorter] == other_parts[:shorter]:
                _fail(f"units の raw_dir が重なっています（同じか、片方がもう片方の中）: "
                      f"{other_raw_dir!r} と {unit['raw_dir']!r}",
                      raw_dir=unit["raw_dir"], other_raw_dir=other_raw_dir,
                      indexes=[other_index, index])
        seen_raw_dir[parts] = (index, unit["raw_dir"])


def _validate_fields(data: dict) -> None:
    """正規形（`_NORMALIZED_KEYS` が全部そろった形）の型・列挙値・数値域を検査する。"""
    if data.get("schema") != SCHEMA:
        _fail(f"保存する要求のschemaは{SCHEMA!r}のみです。", schema=data.get("schema"))
    for field in ("output_root", "keep_extension"):
        _validate_optional_str(data.get(field), field)
    library_mode = data.get("library_mode")
    if library_mode not in _LIBRARY_MODE_VALUES:
        _fail(f"library_modeが不正です: {library_mode!r}", library_mode=library_mode)
    alignment_reference = data.get("alignment_reference")
    if alignment_reference not in _ALIGNMENT_REFERENCE_VALUES:
        _fail(f"alignment_referenceが不正です: {alignment_reference!r}",
              alignment_reference=alignment_reference)
    measure = data.get("measure")
    if measure not in _MEASURE_VALUES:
        _fail(f"measureが不正です: {measure!r}", measure=measure)
    if not _is_positive_int(data.get("timeout_s")):
        _fail(f"timeout_sは正の整数である必要があります（boolは不可）: {data.get('timeout_s')!r}",
              timeout_s=data.get("timeout_s"))
    if not isinstance(data.get("save_project"), bool):
        _fail(f"save_projectはboolのみ許可されます: {data.get('save_project')!r}",
              save_project=data.get("save_project"))
    _validate_units(data)


# ---------- 公開API ----------

def validate_request(data: dict, *, internal: bool = False) -> dict:
    """正規形の要求を検証して返す（`internal=True` は内部キーを持つ保存済みの要求）。"""
    allowed = _NORMALIZED_KEYS | _INTERNAL_KEYS if internal else _NORMALIZED_KEYS
    unknown = set(data) - allowed
    if unknown:
        _fail(f"未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))
    _validate_fields(data)
    return data


def read_request_file(source_root: Path) -> tuple[dict, list[dict]]:
    """元フォルダ直下の `analysis-request.json` を読む（無ければ `({}, [])`）。

    戻り値は (値, 警告)。下流のキー（`target` `preprocess` `comparisons`）は読み飛ばして
    `REQUEST_FILE_DOWNSTREAM_IGNORED` を警告に積む（reanalysis-study の古い要求ファイルで
    Console まで止めない。判断 5）。壊れている・オブジェクトでない・未知のキー・
    許さない null は `PIPELINE_REQUEST_INVALID`（置いてあるのに黙って無視しない）。
    """
    path = Path(source_root) / REQUEST_FILE_NAME
    if not path.is_file():
        return {}, []
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        _fail(f"{REQUEST_FILE_NAME}を読めません: {exc}", path=str(path))
    try:
        data = json.loads(raw)
    except ValueError as exc:
        _fail(f"{REQUEST_FILE_NAME}がJSONとして壊れています: {exc}", path=str(path))
    if not isinstance(data, dict):
        _fail(f"{REQUEST_FILE_NAME}はオブジェクトである必要があります。", path=str(path), value=data)
    warnings: list[dict] = []
    ignored = sorted(set(data) & DOWNSTREAM_KEYS)
    if ignored:
        data = {key: value for key, value in data.items() if key not in DOWNSTREAM_KEYS}
        warnings.append({
            "code": "REQUEST_FILE_DOWNSTREAM_IGNORED",
            "stage_id": None,
            "message": (f"{REQUEST_FILE_NAME} の {ignored} は pipeline の下流の工程として消えたので"
                        "読み飛ばしました。前処理・比較・PCA は load_dataset の後に dataset_preprocess /"
                        " dataset_compare / dataset_pca で行います。"),
            "details": {"ignored_keys": ignored, "replacements": list(DOWNSTREAM_REPLACEMENTS)},
        })
    unknown = set(data) - _TOP_LEVEL_KEYS
    if unknown:
        _fail(f"{REQUEST_FILE_NAME}に未知のキーがあります: {sorted(unknown)}",
              path=str(path), unknown_keys=sorted(unknown))
    _reject_disallowed_explicit_null(data, where=f"{REQUEST_FILE_NAME}の")
    _check_schema(data.get("schema"))
    return data, warnings


def _layer_units(explicit: dict, from_file: dict) -> tuple[list, list]:
    """単位の層を決める。`units` は層ごとに丸ごと、1 極性の形はキーごとに重ねる。"""
    for layer, name in ((explicit, "request"), (from_file, REQUEST_FILE_NAME)):
        if "units" in layer and set(layer) & _UNIT_KEYS:
            _fail(f"{name} に units とトップの単位のキー {sorted(set(layer) & _UNIT_KEYS)} が"
                  "同居しています（2 極性の形では単位のキーを units の中へ置いてください）。",
                  unit_keys=sorted(set(layer) & _UNIT_KEYS))
    if "units" in explicit:
        return _units_from_list(explicit["units"], origin="explicit")
    if "units" in from_file:
        if set(explicit) & _UNIT_KEYS:
            _fail(f"{REQUEST_FILE_NAME} は units（2 極性の形）を持つので、明示の単位のキー"
                  f"{sorted(set(explicit) & _UNIT_KEYS)} をどの極性に効かせるか決められません。"
                  "units ごと明示してください。", unit_keys=sorted(set(explicit) & _UNIT_KEYS))
        return _units_from_list(from_file["units"], origin="request_file")

    def _pick(key):
        if key in explicit:
            return explicit[key], "explicit"
        if key in from_file:
            return from_file[key], "request_file"
        return None, "default"

    picked = {key: _pick(key) for key in _UNIT_KEYS}
    unit = _unit_entry({key: value for key, (value, _src) in picked.items()}, ".")
    sources = {key: picked[key][1] if key in picked else "default"
               for key in sorted(_UNIT_ENTRY_KEYS)}
    return [unit], [sources]


def resolve_request(source_root: Path, explicit: dict | None = None) -> dict:
    """要求を解決し、v3 の正規形（`units` を必ず持つ）に `value_sources` を付けて返す。

    `analysis-request.json` の下流のキーを読み飛ばしたときは `intake_warnings` に警告が入る
    （受付が新しい run の `warnings` へ写す）。
    """
    source_root = Path(source_root)
    explicit = explicit if explicit is not None else {}
    if not isinstance(explicit, dict):
        _fail("requestはオブジェクトである必要があります。", value=explicit)
    downstream = set(explicit) & DOWNSTREAM_KEYS
    if downstream:
        raise _downstream_removed(downstream, where="request")
    unknown = set(explicit) - _TOP_LEVEL_KEYS
    if unknown:
        _fail(f"未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))
    _reject_disallowed_explicit_null(explicit)
    _check_schema(explicit.get("schema"))

    from_file, warnings = read_request_file(source_root)
    units, unit_sources = _layer_units(explicit, from_file)
    shadowed = sorted(set(from_file) & _UNIT_KEYS) if "units" in explicit else []
    if shadowed:
        warnings.append({
            "code": "REQUEST_FILE_UNIT_KEYS_IGNORED",
            "stage_id": None,
            "message": (f"明示の units が単位を丸ごと決めるので、{REQUEST_FILE_NAME} のトップの"
                        f"単位のキー {shadowed} は使いませんでした。効かせたいなら units の"
                        "各単位の中へ書いてください。"),
            "details": {"ignored_keys": shadowed},
        })

    data: dict = {"schema": SCHEMA}
    value_sources: dict = {}
    for key in sorted(_SHARED_KEYS):
        if key in explicit:
            data[key], value_sources[key] = explicit[key], "explicit"
        elif key in from_file:
            data[key], value_sources[key] = from_file[key], "request_file"
        else:
            data[key], value_sources[key] = copy.deepcopy(_SHARED_DEFAULTS[key]), "default"
    data["units"] = units
    value_sources["units"] = unit_sources

    validated = validate_request(data)
    validated["value_sources"] = value_sources
    validated["intake_warnings"] = warnings
    return validated


def unit_index(request: dict, polarity: str) -> int:
    """`polarity` の単位の添字。無ければ `PIPELINE_REQUEST_INVALID`。"""
    for index, unit in enumerate(request.get("units") or []):
        if unit.get("polarity") == polarity:
            return index
    _fail(f"この要求に極性 {polarity!r} の単位はありません。",
          polarity=polarity,
          available=[u.get("polarity") for u in request.get("units") or []])


def pin_polarity(request: dict, index: int, polarity: str) -> dict:
    """1 極性の形で省いた極性を、計画が選んだメソッドの宣言で埋めた写しを返す。

    run の単位・工程・置き場所は極性で名前が付くので、受付が計画を作った直後に埋める。
    既に極性があるなら何もしない（明示が勝つ）。
    """
    out = copy.deepcopy(request)
    unit = out["units"][index]
    if unit.get("polarity") is None:
        unit["polarity"] = polarity
        sources = out.setdefault("value_sources", {}).setdefault("units", [])
        while len(sources) <= index:
            sources.append({})
        sources[index]["polarity"] = "method_declaration"
    return validate_request(out, internal=True)


def unit_source_root(dataset_root: Path, unit: dict) -> Path:
    """単位の生データのフォルダ（`dataset_root / raw_dir`）。"""
    raw_dir = unit.get("raw_dir") or "."
    root = Path(dataset_root)
    return root if raw_dir == "." else root / raw_dir


def unit_request(request: dict, index: int, dataset_root: Path) -> dict:
    """単位 `index` の、入力計画（`inputs.inspect_inputs`）が読む平らな要求を返す。

    共有のキーと単位のキーを 1 つにし、単位の `keep_extension` が共有のものより勝つ。
    `raw_dir` が `.` でない単位では、相対で書かれたパス（`method_file` `lbm_file`
    `msp_file` `sample_manifest`）を `dataset_root` 基準の絶対パスにする——計画は単位の
    フォルダを起点に読むので、そのままだと別のフォルダ基準に化ける。
    """
    unit = request["units"][index]
    flat = {key: request.get(key) for key in _SHARED_KEYS}
    for key in _UNIT_KEYS:
        flat[key] = unit.get(key)
    if unit.get("keep_extension"):
        flat["keep_extension"] = unit["keep_extension"]
    flat["raw_dir"] = unit.get("raw_dir") or "."
    if flat["raw_dir"] != ".":
        base = Path(dataset_root)
        for key in ("method_file", "lbm_file", "msp_file", "sample_manifest"):
            value = flat.get(key)
            if value and not Path(value).expanduser().is_absolute():
                flat[key] = str(base / value)
    return flat


def merge_updates(request: dict, updates: dict) -> dict:
    """resume の入力訂正（試料の対応表だけ）を反映し、再検証した要求を返す。

    下流のキーは `PIPELINE_DOWNSTREAM_REMOVED`、それ以外の上流の条件（メソッド・ライブラリ・
    極性・フォルダ…）は `NEW_PIPELINE_REQUIRED`（上流の条件が変わるなら別の run）。
    """
    if not isinstance(updates, dict):
        _fail("updatesはオブジェクトである必要があります。", value=updates)
    downstream = set(updates) & DOWNSTREAM_KEYS
    if downstream:
        raise _downstream_removed(downstream, where="updates")
    if set(updates) - UPDATABLE:
        raise DomainError(
            "NEW_PIPELINE_REQUIRED", "上流条件の変更には新しい解析が必要です",
            {"keys": sorted(set(updates) - UPDATABLE), "updatable": sorted(UPDATABLE)})
    if "sample_manifest" in updates and "units" in updates:
        _fail("sample_manifest と units は同時に指定できません（2 極性の run は units で極性ごとに）。")

    out = copy.deepcopy(request)
    units = out["units"]
    changed: dict[int, object] = {}
    if "sample_manifest" in updates:
        if len(units) != 1:
            _fail("2 極性の run ではトップの sample_manifest を使えません。"
                  "units: {<極性>: {\"sample_manifest\": …}} で極性ごとに指定してください。",
                  polarities=[u.get("polarity") for u in units])
        changed[0] = updates["sample_manifest"]
    if "units" in updates:
        unit_updates = updates["units"]
        if not isinstance(unit_updates, dict) or not unit_updates:
            _fail("updates.units は {<極性>: {\"sample_manifest\": …}} の形である必要があります。",
                  value=unit_updates)
        for polarity, fields in unit_updates.items():
            index = unit_index(out, polarity)
            if not isinstance(fields, dict):
                _fail(f"updates.units.{polarity} はオブジェクトである必要があります。", value=fields)
            others = set(fields) - {"sample_manifest"}
            if set(fields) & DOWNSTREAM_KEYS:
                raise _downstream_removed(set(fields) & DOWNSTREAM_KEYS,
                                          where=f"updates.units.{polarity}")
            if others:
                raise DomainError(
                    "NEW_PIPELINE_REQUIRED", "上流条件の変更には新しい解析が必要です",
                    {"keys": sorted(others), "updatable": ["sample_manifest"],
                     "polarity": polarity})
            if "sample_manifest" in fields:
                changed[index] = fields["sample_manifest"]
    for index, value in changed.items():
        _validate_optional_str(value, "sample_manifest")
        units[index]["sample_manifest"] = copy.deepcopy(value)
        sources = out.setdefault("value_sources", {}).setdefault("units", [])
        while len(sources) <= index:
            sources.append({})
        sources[index]["sample_manifest"] = "explicit_update"
    return validate_request(out, internal=True)


def request_fingerprint(request: dict) -> str:
    """要求の内容 hash（`value_sources` のような由来情報は外す）。

    受付の冪等性（`store.find_or_create_run`）と再開の同一性はこの hash で決まる。v3 の
    正規形の上で計算するので、旧 v1 の run とは同じ要求でも hash が違う（再利用されない）。
    既定値のままの `library_mode` / `alignment_reference` は落として数える。
    """
    content = {key: request[key] for key in _NORMALIZED_KEYS if key in request}
    for key, default in _FINGERPRINT_OMITTED_DEFAULTS.items():
        if key in content and content[key] == default:
            del content[key]
    return canonical_hash(content)
