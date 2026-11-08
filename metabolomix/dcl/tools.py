"""DCL（`.dcl` / MSDecResult＝デコンボリューション済み MS/MS）ツール群。

MS/MS は同定確度の最強証拠だが、MS-DIAL は実スペクトルを `.pai2` ではなく同名の
`.dcl` に置く（`.pai2` の `has_msms` は取得参照の有無を示すだけ）。本モジュールは
その実体を MCP 面へ出す入口（`dcl_find_msms` / `dcl_spot_spectrum` / `dcl_diagnostic_ions`）。

deps: mcp_core / session_state / tools._holder / dcl_reader。tools_* / server は
import しない。
"""
import os
from pathlib import Path

from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.dcl.reader import deserialize_dcl, get_msms_by_precursor

__all__ = ["dcl_find_msms", "dcl_spot_spectrum", "dcl_diagnostic_ions"]

# LLM へ返すスペクトルは主要フラグメントに絞る（全ピークは容易に数百本になる）。
_DEFAULT_TOP_PEAKS = 10


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dcl_find_msms(precursor_mz: float, file_path: str | None = None,
                  rt: float | None = None, mz_tol: float = 0.01,
                  rt_tol: float = 0.2, top_n_peaks: int | None = None,
                  dataset_handle: str | None = None) -> str:
    """precursor m/z（任意で RT）に一致する MS/MS を `.dcl` から引く。

    アノテーションの裏取り（期待フラグメントが実際に出ているか）に使う。RT を渡すと
    同一 m/z の別溶出ピークを分離できる。該当なしは `status="not_found"`。
    file_path: 読む `.dcl`。省けば直近のデータセット（または dataset_handle）の `.dcl`
        （データフォルダは探さない。データセットが無ければ missing_state）。
    """
    from metabolomix.tools._holder import alignment_paths, resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    if file_path is not None and not os.path.isfile(file_path):
        return json_payload({"status": "error", "message":
                             f"指定された .dcl（{os.path.basename(file_path)}）がありません。"})
    _arf2, resolved, lacking = alignment_paths(holder, dataset_handle, None, file_path,
                                               need=("dcl",))
    if lacking:
        return json_payload(lacking)
    limit = _DEFAULT_TOP_PEAKS if top_n_peaks is None else (top_n_peaks or None)
    try:
        results = deserialize_dcl(resolved, include_spectrum=True, top_n_peaks=limit)
    except Exception as exc:  # noqa: BLE001
        return json_payload({"status": "error", "message": f"DCL 解析に失敗しました: {exc}"})

    hits = get_msms_by_precursor(results, precursor_mz, tol=mz_tol, rt=rt, rt_tol=rt_tol)
    if not hits:
        return json_payload({
            "status": "not_found",
            "message": (
                f"precursor m/z={precursor_mz}（±{mz_tol}"
                + (f", RT={rt}±{rt_tol}" if rt is not None else "")
                + f"）に一致する MS/MS は {Path(resolved).name} にありません。"
                "MS/MS 未取得のピークである可能性を、同定確度の判断に反映すること。"
            ),
        })

    return json_payload({
        "status": "success",
        "file": Path(resolved).name,
        "query": {"precursor_mz": precursor_mz, "mz_tol": mz_tol,
                  "rt": rt, "rt_tol": rt_tol},
        "matches": [
            {
                "dcl_index": r["dcl_index"],
                "precursor_mz": round(r["precursor_mz"], 4),
                "rt": round(r["rt"], 2),
                "n_msms_peaks": r["n_msms_peaks"],
                "signal_to_noise": round(r["signal_to_noise"], 2),
                "msms_spectrum": [[round(m, 4), i] for m, i in r["msms_spectrum"]],
            }
            for r in hits
        ],
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dcl_spot_spectrum(spot_id: int, arf2_path: str | None = None, dcl_file: str | None = None,
                      dataset_handle: str | None = None) -> str:
    """アラインメントのスポット ID で代表 MS/MS（全ピーク）を引く（再解析で使う）。

    `dcl_find_msms` は precursor m/z で引くので、同じ m/z の別スポットを拾いうる。ここは
    `library_match_feature(spot_id=...)` と同じ規則で、`.arf2` の precursor・RT と照らした
    アラインメント `.dcl` の `dcl[spot_id]` を返す（索引の崩れは `dcl_precursor_mismatch`、
    ID が無ければ `spot_not_found`、MS/MS が無ければ `not_found`）。mzTab-M の SMF_ID は
    スポット ID と同じ。
    dataset_handle: `load_dataset` が返したハンドル（省略時は直近）。`arf2_path` / `dcl_file` を省けば
        そのデータセットの `.arf2` / `.dcl` を使う（データフォルダは探さない）。
    """
    from metabolomix.library.tools import spot_msms_query
    from metabolomix.tools._holder import alignment_paths, resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    arf2_path, dcl_file, lacking = alignment_paths(holder, dataset_handle, arf2_path, dcl_file)
    if lacking:
        return json_payload(lacking)
    spot = spot_msms_query(int(spot_id), arf2_path=arf2_path, dcl_file=dcl_file, rt=None,
                           ion_mode=None)
    if spot["status"] != "ok":
        return json_payload(spot)
    return json_payload({
        "status": "ok",
        "spot_id": int(spot_id),
        "precursor_mz": round(float(spot["precursor_mz"]), 4),
        "rt": None if spot["rt"] is None else round(float(spot["rt"]), 3),
        "ion_mode": spot["ion_mode"],
        "n_peaks": len(spot["measured"]),
        "peaks": [[round(float(mz), 4), float(f"{float(i):.6g}")] for mz, i in spot["measured"]],
    })


#: dcl_diagnostic_ions が 1 回に引くスポットの上限。
_MAX_DIAGNOSTIC_SPOTS = 200
_DIAGNOSTIC_COLUMNS = ["spot_id", "label", "kind", "target_mz", "found", "observed_mz",
                       "relative_intensity", "status"]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dcl_diagnostic_ions(spot_ids: list[int], ions: list[dict], arf2_path: str | None = None,
                        dcl_file: str | None = None, dataset_handle: str | None = None) -> str:
    """スポットの代表 MS/MS に、診断イオン・中性損失があるかを照会する（再解析で使う）。

    spot_ids: アラインメントのスポット ID（mzTab-M の SMF_ID）。1 回 200 件まで。
    ions: 各要素は `{"mz": 184.0733}`（診断イオン）か `{"neutral_loss": 141.0191}`（前駆体からの
        中性損失）に、任意で `label`・`tol`（Da、既定 0.01）・`min_relative_intensity`（診断イオン
        だけ。最大ピークに対する比、既定 0）。
    TSV（`spot_id label kind target_mz found observed_mz relative_intensity status`）。
    `status` が `ok` 以外の行（MS/MS が無い・索引の崩れ）は照会していない（`found` は空）——
    「イオンが無かった」と読まない。
    dataset_handle: そのハンドルの `.arf2` / `.dcl` を使う。省略時は直近のデータセット。
    """
    from metabolomix.analysis import msms_rules
    from metabolomix.arf2.reader import format_spots_as_table
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.core.mcp_errors import mztab_error
    from metabolomix.library.tools import spot_msms_query
    from metabolomix.tools._holder import alignment_paths, resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    if not isinstance(spot_ids, list) or not spot_ids or len(spot_ids) > _MAX_DIAGNOSTIC_SPOTS:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"spot_ids は 1〜{_MAX_DIAGNOSTIC_SPOTS} 件のリスト", {"field": "spot_ids"})
    if not isinstance(ions, list) or not ions:
        return mztab_error("DATASET_BAD_REQUEST", "ions は 1 件以上のリスト", {"field": "ions"})
    try:
        normalized = []
        for i, ion in enumerate(ions):
            if not isinstance(ion, dict) or ("mz" in ion) == ("neutral_loss" in ion):
                raise DomainError("MSMS_RULES_INVALID",
                                  f"ions[{i}] は mz か neutral_loss のどちらか一方を持つ辞書",
                                  {"field": f"ions[{i}]"})
            if "mz" in ion:
                normalized += msms_rules.normalize_ions([ion], kind="ion", label=f"ions[{i}]")
            else:
                renamed = {("delta" if k == "neutral_loss" else k): v for k, v in ion.items()}
                normalized += msms_rules.normalize_ions([renamed], kind="loss", label=f"ions[{i}]")
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    arf2_path, dcl_file, lacking = alignment_paths(holder, dataset_handle, arf2_path, dcl_file)
    if lacking and "error" in lacking:  # データセットが無い（missing_state）。行には並べない
        return json_payload(lacking)
    rows = []
    for spot_id in spot_ids:
        spot = lacking or spot_msms_query(int(spot_id), arf2_path=arf2_path, dcl_file=dcl_file,
                                          rt=None, ion_mode=None)
        if spot["status"] != "ok":
            rows += [{"spot_id": int(spot_id), "label": ion["label"],
                      "kind": "ion" if "mz" in ion else "neutral_loss",
                      "target_mz": None, "found": None, "observed_mz": None,
                      "relative_intensity": None, "status": spot["status"]} for ion in normalized]
            continue
        for hit in msms_rules.find_ions(spot["measured"], normalized,
                                        precursor_mz=spot["precursor_mz"]):
            rows.append({"spot_id": int(spot_id), **hit, "status": "ok"})
    cells = []
    for row in rows:
        cell = {}
        for column in _DIAGNOSTIC_COLUMNS:
            value = row.get(column)
            if isinstance(value, bool):
                value = "true" if value else "false"
            elif isinstance(value, float):
                value = round(value, 4)
            cell[column] = value
        cells.append(cell)
    return json_payload({"status": "success", "n_spots": len(spot_ids),
                         "ions": format_spots_as_table(cells, columns=_DIAGNOSTIC_COLUMNS)})
