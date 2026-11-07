"""POS と NEG の統合（Task 12）。

正規化済み（または正規化を飛ばした）2 つのデータセットを、特徴を縦に連結した 1 つの派生
データセットにする。試料は規則プラグインの `match_samples` で共通 ID に照合し、食い違えば
黙って部分集合で統合せず止まる（`QUANT_MERGE_SAMPLE_MISMATCH`。対応表で進められる）。

特徴 ID は `positive:<id>` / `negative:<id>`（msdial-workflow と同じ書式）。`derive` は親に無い
ID の副次データを落とすので、特徴ごとの副次データ（候補・注釈・タグ・試料ごとの m/z）は
ここで両側の ID に付け替えて載せる。統合前のキュレーションのフラグは極性ごとに
`quant_state["curation_by_polarity"]` へ移す（統合後の ID と ID 空間が違うため）。
"""
from __future__ import annotations

import copy
import csv
import hashlib
from pathlib import Path

import numpy as np

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import missing_state
from metabolomix.core.serialization import json_payload
from metabolomix.quant import derive as qderive
from metabolomix.quant.inputs import quant_error
from metabolomix.quant.normalize import analysis_folder
from metabolomix.quant.plugin import call, load_rules, rules_settings
from metabolomix.quant.records import merged_record_folder, write_record
from metabolomix.quant.table import resolve_quant_polarity

__all__ = ["MAPPING_COLUMNS", "merge", "merge_tool", "read_mapping"]

MAPPING_COLUMNS = ("mode", "source_sample", "sample_id", "group", "sample_type")
#: 規則プラグイン（msdial-workflow の config.MODES・対応表の `mode` 列）が使うモード名。
MODES = {"positive": "positive", "negative": "negative"}
_MODE_ALIASES = {"pos": "positive", "positive": "positive", "neg": "negative",
                 "negative": "negative"}
#: ハンドルのロックを待つ秒数（超えたら QUANT_HANDLE_BUSY）。
LOCK_TIMEOUT = 30.0
_MAX_LISTED = 20
_HOW_MAPPING = ("mapping_csv（列: mode, source_sample, sample_id, group, sample_type。"
                "mode は positive / negative）")


def read_mapping(path: Path) -> tuple[list[dict], str]:
    """対応表 CSV の (行, sha256)。列が足りなければ QUANT_BAD_REQUEST。"""
    path = Path(path)
    if not path.is_file():
        raise DomainError("QUANT_MAPPING_NOT_FOUND", f"対応表が見つかりません: {path.name}",
                          {"label": path.name})
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise DomainError("QUANT_MAPPING_NOT_FOUND",
                          f"対応表を読めません: {path.name}（{type(exc).__name__}）",
                          {"label": path.name}) from None
    reader = csv.DictReader(data.decode("utf-8-sig").splitlines())
    missing = [c for c in MAPPING_COLUMNS if c not in (reader.fieldnames or [])]
    if missing:
        raise DomainError("QUANT_BAD_REQUEST",
                          f"対応表の列が足りません: {', '.join(missing)}"
                          f"（必要な列: {', '.join(MAPPING_COLUMNS)}）",
                          {"field": "mapping_csv", "missing": missing})
    rows = [{k: (v or "").strip() for k, v in row.items() if k is not None} for row in reader]
    for row in rows:  # pos / NEG などの別名を、規則が比べる正式名へ
        row["mode"] = _MODE_ALIASES.get(row["mode"].lower(), row["mode"])
    return rows, hashlib.sha256(data).hexdigest()


def _check_polarities(pos_ds, neg_ds, polarities) -> None:
    if polarities is None:
        polarities = (resolve_quant_polarity(pos_ds, None, None),
                      resolve_quant_polarity(neg_ds, None, None))
    if tuple(polarities) != ("positive", "negative"):
        raise DomainError("QUANT_MERGE_POLARITY_INVALID",
                          "positive_handle は陽イオン、negative_handle は陰イオンのデータセットを"
                          "渡してください（同じ極性の統合・取り違えはしません）。",
                          {"positive": polarities[0], "negative": polarities[1]})


def _check_states(pos_ds, neg_ds) -> tuple[dict, dict]:
    states = []
    for label, ds in (("positive", pos_ds), ("negative", neg_ds)):
        state = getattr(ds, "quant_state", None)
        if not state:
            raise DomainError(
                "QUANT_NORMALIZE_STEP_REQUIRED",
                f"{label} のデータセットは正規化の工程を通っていません。統合の前に "
                "dataset_quant_apply で正規化してください（正規化しないなら skip=True で"
                "飛ばした記録を作ります）。",
                {"polarity": label})
        states.append(state)
    pos, neg = states
    diffs = {k: [pos.get(k), neg.get(k)] for k in ("normalized", "method", "unit")
             if pos.get(k) != neg.get(k)}
    if diffs:
        raise DomainError("QUANT_MERGE_STATE_MISMATCH",
                          "正規化の状態（normalized・method・unit）が 2 つのデータセットで違うため、"
                          "統合しません。両方に同じ正規化を掛けるか、両方を飛ばしてください。",
                          {"differs": diffs})
    return pos, neg


def _source_names(plugin, names, mismatched, result, mapping_rows) -> dict:
    """食い違った共通 ID を、利用者が書く元の試料名に戻す（戻せなければ ID のまま）。"""
    out = {}
    for mode, bad in mismatched.items():
        bad_set = set(bad)
        rows = (result.get("rows_by_mode") or {}).get(mode)
        found = []
        for name in names.get(mode, []):
            if mapping_rows is not None and rows is not None and name in rows:
                sid = rows[name].get("sample_id")
            else:
                try:
                    single = call(plugin, "match_samples", names_by_mode={mode: [name]})
                    sid = single["ids_by_mode"][mode][0]
                except Exception:  # noqa: BLE001 規則が単独の照合を拒むなら ID のまま返す
                    found = None
                    break
            if sid in bad_set:
                found.append(name)
        out[mode] = found if found else list(bad)
    return out


def _match(plugin, pos_ds, neg_ds, mapping_rows):
    names = {MODES["positive"]: list(pos_ds.sample_names),
             MODES["negative"]: list(neg_ds.sample_names)}
    result = call(plugin, "match_samples", names_by_mode=names, mapping_rows=mapping_rows)
    mismatched = {m: list(v) for m, v in (result.get("mismatched") or {}).items() if v}
    if mismatched:
        raise DomainError("QUANT_MERGE_SAMPLE_MISMATCH",
                          "試料名が 2 つのデータセットで揃いません。mapping_csv で対応づけてください。",
                          {"rule": result.get("rule"), "only_in_ids": mismatched,
                           "only_in": _source_names(plugin, names, mismatched, result,
                                                    mapping_rows)})
    ids = result["ids_by_mode"]
    pos_ids, neg_ids = list(ids[MODES["positive"]]), list(ids[MODES["negative"]])
    for label, got, src in (("positive", pos_ids, pos_ds), ("negative", neg_ids, neg_ds)):
        if len(got) != len(src.sample_names) or len(set(got)) != len(got):
            raise DomainError("QUANT_MERGE_SAMPLE_AMBIGUOUS",
                              f"{label} の試料を共通 ID へ 1 対 1 に対応づけられません。"
                              "mapping_csv で対応づけてください。",
                              {"rule": result.get("rule"), "polarity": label})
    if set(pos_ids) != set(neg_ids):  # 規則が mismatched に出さなかった場合の保険
        raise DomainError("QUANT_MERGE_SAMPLE_MISMATCH", "試料名が 2 つのデータセットで揃いません。",
                          {"rule": result.get("rule"),
                           "only_in": {"pos": sorted(set(pos_ids) - set(neg_ids)),
                                       "neg": sorted(set(neg_ids) - set(pos_ids))}})
    return result.get("rule"), pos_ids, neg_ids


def _group_role(ds, position: int, name: str) -> dict:
    out = {}
    samples = getattr(ds, "msdial_samples", None)
    if samples and position < len(samples) and isinstance(samples[position], dict):
        out["group"] = samples[position].get("group")
        out["role"] = samples[position].get("role")
    role = (getattr(ds, "roles", None) or {}).get(name)
    if role is not None:
        out["role_map"] = role
    return out


def _check_groups(pos_ds, neg_ds, pos_ids, neg_pos) -> None:
    conflicts = []
    for i, sid in enumerate(pos_ids):
        j = neg_pos[sid]
        a = _group_role(pos_ds, i, pos_ds.sample_names[i])
        b = _group_role(neg_ds, j, neg_ds.sample_names[j])
        for field in ("group", "role", "role_map"):
            if field in a and field in b and a[field] != b[field]:
                conflicts.append({"sample_id": sid, "field": field,
                                  "positive": a[field], "negative": b[field]})
    if conflicts:
        raise DomainError("QUANT_MERGE_SAMPLE_CONFLICT",
                          "同じ試料の群・役割が陽イオンと陰イオンで食い違います。実験情報をそろえて"
                          "から統合してください。",
                          {"n_conflicts": len(conflicts), "conflicts": conflicts[:_MAX_LISTED]})


def _is_ids(state: dict, prefix: str) -> list[str]:
    return [f"{prefix}:{i}" for i in state.get("internal_standard_ids") or []]


def _union(a, b) -> list:
    return list(a) + [x for x in b if x not in a]


def _rename_inner(per_sample, names, ids, prefix) -> dict:
    rename = dict(zip(names, ids))
    return {f"{prefix}:{fid}": {rename[k]: v for k, v in inner.items() if k in rename}
            for fid, inner in (per_sample or {}).items()}


def merge(pos_ds, neg_ds, *, mapping_rows: list[dict] | None = None,
          pos_handle_id: str | None = None, neg_handle_id: str | None = None,
          mapping_input: tuple[str, str] | None = None, folder: Path | None = None,
          polarities: tuple[str, str] | None = None):
    """統合して記録を書き、派生データセットを返す。`(merged, info)`。

    `mapping_input` は (対応表のファイル名, sha256)（記録の入力に `mapping:<名前>` で載る）。
    正規化の状態（normalized・method・unit）が 2 つで違えば統合しない（method の違いも）。
    対応表の group / sample_type 列は規則が検証するだけで適用しない（試料の実験情報は POS の
    ものを使う）。
    `folder` を省くと 2 つの解析フォルダの共通の親（無理なら QUANT_RECORD_LOCATION_REQUIRED）。
    `polarities` は (positive 側, negative 側) の極性（省くとデータセットから決める）。
    記録は検査がすべて通ってから書く（止まるときは何も書かない）。
    """
    plugin = load_rules()
    _check_polarities(pos_ds, neg_ds, polarities)
    pos_state, neg_state = _check_states(pos_ds, neg_ds)
    rule, pos_ids, neg_ids = _match(plugin, pos_ds, neg_ds, mapping_rows)
    neg_pos = {sid: j for j, sid in enumerate(neg_ids)}
    _check_groups(pos_ds, neg_ds, pos_ids, neg_pos)
    if folder is None:
        folder = merged_record_folder(analysis_folder(pos_ds), analysis_folder(neg_ds))
    folder = Path(folder)

    order = [neg_pos[sid] for sid in pos_ids]  # 共通 ID（POS の順）ごとの NEG の列
    matrix = np.vstack([np.asarray(pos_ds.feature_matrix, dtype=float),
                        np.asarray(neg_ds.feature_matrix, dtype=float)[:, order]])
    detected = np.vstack([np.asarray(pos_ds.detected_mask, dtype=bool),
                          np.asarray(neg_ds.detected_mask, dtype=bool)[:, order]])
    fids, metadata = [], {}
    for prefix, ds, handle_id in (("positive", pos_ds, pos_handle_id),
                                  ("negative", neg_ds, neg_handle_id)):
        files = sorted(Path(p).name for p in (ds.source_files or {}))
        for fid in ds.feature_ids:
            new_id = f"{prefix}:{fid}"
            fids.append(new_id)
            meta = copy.deepcopy((ds.feature_metadata or {}).get(fid, {}))
            meta["origin"] = {"polarity": prefix, "handle_id": handle_id, "feature_id": fid,
                              "source_files": files}
            metadata[new_id] = meta

    inputs = {**pos_ds.source_files, **neg_ds.source_files}
    if mapping_input is not None:
        inputs[f"mapping:{mapping_input[0]}"] = mapping_input[1]
    parents = [pos_state.get("record_id"), neg_state.get("record_id")]
    n_pos, n_neg = len(pos_ds.feature_ids), len(neg_ds.feature_ids)
    counts = {"positive": n_pos, "negative": n_neg, "merged": n_pos + n_neg}
    record = write_record(
        folder, "merge", inputs=inputs, parent_ids=[p for p in parents if p],
        settings={"rules": rules_settings(plugin), "match_rule": rule,
                  "mapping": mapping_input[0] if mapping_input else None,
                  "handles": {"positive": pos_handle_id, "negative": neg_handle_id}},
        summary={"n_features": counts, "n_samples": len(pos_ids), "unit": pos_state.get("unit"),
                 "normalized": bool(pos_state.get("normalized"))})
    state = qderive.quant_state(
        normalized=bool(pos_state.get("normalized")), method=pos_state.get("method"),
        unit=pos_state.get("unit"), record_id=record["record_id"],
        steps=_union(pos_state.get("steps") or [], neg_state.get("steps") or []) + ["merge"],
        skipped=_union(pos_state.get("skipped") or [], neg_state.get("skipped") or []),
        record_folder=str(folder),
        curation_by_polarity={"positive": copy.deepcopy(pos_ds.curation_flags),
                              "negative": copy.deepcopy(neg_ds.curation_flags)},
        internal_standard_ids=_is_ids(pos_state, "positive") + _is_ids(neg_state, "negative"))
    new = qderive.derive(pos_ds, matrix=matrix, detected=detected, feature_ids=fids,
                         feature_metadata=metadata, sample_names=pos_ids,
                         sample_index=list(range(len(pos_ids))), state=state)

    # derive は親に無い ID の副次データを落とすので、両側の ID に付け替えて載せる。
    for attr in ("feature_candidates", "feature_annotations", "feature_tags"):
        merged_attr = {}
        for prefix, ds in (("positive", pos_ds), ("negative", neg_ds)):
            merged_attr.update({f"{prefix}:{k}": copy.deepcopy(v)
                                for k, v in (getattr(ds, attr) or {}).items()})
        setattr(new, attr, merged_attr)
    if pos_ds.per_sample_mz is None and neg_ds.per_sample_mz is None:
        new.per_sample_mz = None
    else:
        new.per_sample_mz = {
            **_rename_inner(pos_ds.per_sample_mz, pos_ds.sample_names, pos_ids, "positive"),
            **_rename_inner(neg_ds.per_sample_mz, neg_ds.sample_names, neg_ids, "negative")}
    # 片側（POS）から引き継いだ欄のうち、統合後に意味が変わるものを直す。
    new.source_files = dict(inputs)
    new.artifact_sha256 = {**(pos_ds.artifact_sha256 or {}), **(neg_ds.artifact_sha256 or {})}
    new.artifact_paths = {}
    new.job_path = None
    new.ion_mode_majority = None  # 陽陰の両方を持つので極性は決まらない
    new.curation_flags = None  # 極性ごとに quant_state["curation_by_polarity"] にある
    new.sample_metadata_rows = None
    new.library_record = None
    info = {"record_id": record["record_id"], "rule": rule, "n_features": counts,
            "n_samples": len(pos_ids), "unit": state["unit"], "normalized": state["normalized"]}
    return new, info


def _needs_sample_mapping(exc: DomainError) -> str:
    return json_payload({"status": "needs_input", "reason": "sample_mismatch",
                         "rule": (exc.details or {}).get("rule"),
                         "only_in": (exc.details or {}).get("only_in"),
                         "message": exc.message, "how_to_answer": _HOW_MAPPING})


def merge_tool(pos_holder, neg_holder, mapping_csv: str | None, record_folder: str | None) -> str:
    """dataset_merge_polarities の実体（ロック・読み取り・登録・封筒）。"""
    holders = []
    for holder in (pos_holder, neg_holder):
        if getattr(holder, "dataset", None) is None:
            return missing_state("dataset", ["dataset_load", "dataset_handle_open"],
                                 "データセットがありません。dataset_handle_open で開いたハンドルを"
                                 "渡してください。")
        holders.append(holder)
    folder = None
    if record_folder:
        folder = Path(record_folder).expanduser()
        if not folder.is_dir():
            return quant_error("QUANT_RECORD_FOLDER_INVALID",
                               f"record_folder は既存のフォルダを指してください: {folder.name}",
                               {"label": folder.name})
    # 2 つのハンドルのロックを handle_id の順に取る（逆順の呼び出しと噛み合っても止まらない）
    ordered = sorted({id(h): h for h in holders}.values(), key=lambda h: h.handle_id)
    acquired = []
    try:
        for handle in ordered:
            if not handle.lock.acquire(timeout=LOCK_TIMEOUT):
                return quant_error("QUANT_HANDLE_BUSY",
                                   "ハンドルが別の処理で使われていて待ちきれませんでした。"
                                   "少し待ってからやり直してください。",
                                   {"dataset_handle": handle.handle_id})
            acquired.append(handle)
        pos_ds, neg_ds = pos_holder.dataset, neg_holder.dataset
        try:
            mapping_rows, mapping_input = None, None
            if mapping_csv:
                mapping_rows, sha = read_mapping(Path(mapping_csv).expanduser())
                mapping_input = (Path(mapping_csv).name, sha)
            polarities = (resolve_quant_polarity(pos_ds, pos_holder, None),
                          resolve_quant_polarity(neg_ds, neg_holder, None))
            merged, info = merge(
                pos_ds, neg_ds, mapping_rows=mapping_rows, pos_handle_id=pos_holder.handle_id,
                neg_handle_id=neg_holder.handle_id, mapping_input=mapping_input, folder=folder,
                polarities=polarities)
            handle_id = qderive.register(
                pos_holder, merged, step="merge",
                parent_handle_ids=[pos_holder.handle_id, neg_holder.handle_id])
        except DomainError as exc:
            if exc.code == "QUANT_MERGE_SAMPLE_MISMATCH":
                return _needs_sample_mapping(exc)
            if exc.code == "QUANT_RECORD_LOCATION_REQUIRED":
                return json_payload({
                    "status": "needs_input", "reason": "record_folder_required",
                    "message": exc.message,
                    "how_to_answer": "record_folder（記録を置く既存のフォルダ）"})
            return quant_error(exc.code, exc.message, exc.details or None)
    finally:
        for handle in reversed(acquired):
            handle.lock.release()
    return json_payload({"status": "merged", "dataset_handle": handle_id, **info})

