"""DCL（`.dcl` / MSDecResult＝デコンボリューション済み MS/MS）ツール群。

MS/MS は同定確度の最強証拠だが、MS-DIAL は実スペクトルを `.pai2` ではなく同名の
`.dcl` に置く（`.pai2` の `has_msms` は取得参照の有無を示すだけ）。本モジュールは
その実体を MCP 面へ出し、`verify_peak_annotation` の判断材料に載せるための入口。

deps: mcp_core / session_state / path_resolvers / dcl_reader。tools_* / server は
import しない。
"""
from pathlib import Path

from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.core.path_resolvers import resolve_dcl_file_path
from metabolomix.dcl.reader import deserialize_dcl, summarize_dcl, get_msms_by_precursor

__all__ = ["dcl_parser", "dcl_find_msms"]

# LLM へ返すスペクトルは主要フラグメントに絞る（全ピークは容易に数百本になる）。
_DEFAULT_TOP_PEAKS = 10


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dcl_parser(file_path: str | None = None, top_n_peaks: int | None = None,
               preview: int = 5) -> str:
    """1つの `.dcl`（MSDecResult）を解析し、MS/MS 在庫の要約と先頭数件を返す。

    `.dcl` は測定ファイル1つにつき1個で、`dcl_index` が同名 `.pai2` のピーク順に
    対応する。返すのは取得率・ピーク本数分布・precursor m/z / RT 範囲の要約と、
    先頭 `preview` 件の主要フラグメント。個々のピークの同定確度に効かせるには
    `pai2_parser`（同名 `.dcl` を自動で付与する）→ `verify_peak_annotation` を使う。

    - `top_n_peaks`: 各スペクトルを強度上位 N 本に間引く（既定 10。全件は 0 を指定）。
    - `preview`: 要約に添える MSDecResult の件数。
    """
    resolved = resolve_dcl_file_path(file_path)
    if not resolved:
        return "データディレクトリに .dcl ファイルが見つかりませんでした。"

    limit = _DEFAULT_TOP_PEAKS if top_n_peaks is None else (top_n_peaks or None)
    try:
        results = deserialize_dcl(resolved, include_spectrum=True, top_n_peaks=limit)
    except Exception as exc:  # noqa: BLE001 - 破損ファイルは文言で返す（MCP が扱いやすい）
        return f"[ERROR] DCL 解析に失敗しました: {exc}"

    payload = {
        "status": "success",
        "file": Path(resolved).name,
        "summary": summarize_dcl(results),
        "preview": [
            {
                "dcl_index": r["dcl_index"],
                "precursor_mz": round(r["precursor_mz"], 4),
                "rt": round(r["rt"], 2),
                "n_msms_peaks": r["n_msms_peaks"],
                "signal_to_noise": round(r["signal_to_noise"], 2),
                "top_fragments": [[round(m, 4), i] for m, i in r["msms_spectrum"]],
            }
            for r in results[:max(0, preview)]
        ],
    }
    text = (
        f"### DCL 解析完了: {Path(resolved).name}\n"
        + json_payload(payload)
    )
    return session_state.session.maybe_prepend_caveat(text, topic="dcl")


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dcl_find_msms(precursor_mz: float, file_path: str | None = None,
                  rt: float | None = None, mz_tol: float = 0.01,
                  rt_tol: float = 0.2, top_n_peaks: int | None = None) -> str:
    """precursor m/z（任意で RT）に一致する MS/MS を `.dcl` から引く。

    アノテーションの裏取り（期待フラグメントが実際に出ているか）に使う。RT を渡すと
    同一 m/z の別溶出ピークを分離できる。該当なしは `status="not_found"`。
    """
    resolved = resolve_dcl_file_path(file_path)
    if not resolved:
        return json_payload({"status": "error", "message": "データディレクトリに .dcl ファイルが見つかりませんでした。"})
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
