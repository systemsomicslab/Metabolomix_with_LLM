"""データセットのハンドルを開閉する MCP ツール（spec 2026-10-22 §4 の 1）。

**再解析で使う**（並行のエージェントが、別々のデータセット・極性を互いに干渉せずに扱う）。
対話の使い方（`dataset_load` → `dataset_*`）は変わらない——ハンドルを渡さない呼び出しは
対話のセッションの直近のハンドル（`dataset_load` で読んだもの）に効く。

読み込んだデータセットはすべてハンドルの登録簿（`core.handles.registry`）に住む
（spec 2026-11-06 §5）。`dataset_load` も同じ登録簿にハンドルを足して直近にする。
`dataset_handle_open` は `dataset_load` と同じ入口（`dataset.opening.open_dataset`）を呼ぶが、
直近は動かさない。同じ入力を開き直したときは読み込み済みの基を使い回し、ハンドルには複製を
渡す（`open_dataset` が持つ）。

ライブラリは極性で設定（`[library] msp_positive` / `msp_negative`。環境変数
`MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）だけから引き（`library_mode=msp_only` と同じ規則）、
読み取り専用で開く。戻り値にはファイル名と sha256 だけを出し、パスは出さない。
"""
from __future__ import annotations

import re
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import handles, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.dataset import opening

__all__ = ["dataset_handle_open", "dataset_handle_list", "dataset_handle_close"]

_LIBRARY_MODES = ("auto", "none")
_ARF_SUFFIX_RE = re.compile(r"_(?:PeakProperties|DriftSopts|DriftSpots)\.arf$", re.IGNORECASE)


def alignment_files_for(opened, arf2_path: str | None = None) -> tuple[str | None, str | None]:
    """`open_dataset` が開いたデータセットのスポットの照合に使う `.arf2` と `.dcl`。

    `dataset_handle_open` と `dataset_load`（直近のハンドル）が同じ規則で決め、ハンドルに
    `arf2_path` / `dcl_path` として持たせる。ハンドル ID を明示した呼び出しはこの値を使う。
    **ID を省いた呼び出しは `_holder.alignment_paths` とデータフォルダの探索で決めるので、
    同じファイルになるとは限らない**（データフォルダに別のバッチがあれば食い違いうる）。
    解析フォルダで読んだなら（明示が無ければ）その選択、そうでなければ `_alignment_files` が探す。
    """
    selection = opened.selection
    if selection is not None and not arf2_path:
        return str(selection.arf2), (str(selection.dcl) if selection.dcl else None)
    return _alignment_files(opened.dataset, arf2_path, opened.kind, opened.path)


def _alignment_files(ds, arf2_path: str | None, kind: str | None = None,
                     path: str | None = None) -> tuple[str | None, str | None]:
    """スポットの照合に使う `.arf2` と `.dcl`。決まらなければ None。

    探す順: 明示の `arf2_path` → 検出状態を取り込んだ `.arf` の同じバッチ → mzTab の隣
    → （`kind == "pipeline_path"` のとき）run の `input/AlignResult-*.arf2`。
    """
    from metabolomix.curation import apply as curation_apply
    from metabolomix.msdial import alignment_spot

    arf2 = Path(arf2_path) if arf2_path else None
    if arf2 is None:
        arf = ((getattr(ds, "feature_qc", None) or {}).get("join") or {}).get("path")
        if arf:
            candidate = Path(arf).with_name(_ARF_SUFFIX_RE.sub("", Path(arf).name) + ".arf2")
            if candidate.is_file():
                arf2 = candidate
    if arf2 is None:
        arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p, feature_mz=curation_apply.feature_mz_of(ds))
                                  for p in (getattr(ds, "source_files", None) or {}))), None)
    if arf2 is None and kind == "pipeline_path" and path:
        arf2 = opening.pipeline_arf2(path)
    if arf2 is None or not arf2.is_file():
        return None, None
    dcl = alignment_spot.alignment_sibling(arf2, ".dcl")
    return str(arf2), (str(dcl) if dcl else None)


def _open_library_for(ds) -> tuple[object | None, str | None]:
    """(store, note)。データセットが覚えたライブラリ（`ds.library_record`）を開く。

    記憶が `known` で `verify_record` が通ればその実体を読み取り専用で開く。`mismatch` /
    `unknown` / `lbm_only` は開かず理由を返す。**極性の設定（`[library] msp_positive` /
    `msp_negative`）には落ちない**——データセットと違う極性のライブラリを掴まないため。
    理由にはファイル名だけを書く（研究室のライブラリの置き場所を戻り値に出さない）。
    """
    from metabolomix.dataset.library_memory import unknown_record, verify_record
    from metabolomix.library.store import open_store

    record = verify_record(getattr(ds, "library_record", None) or unknown_record())
    status, file = record.get("status"), record.get("file")
    if status == "mismatch":
        return None, (f"記憶したライブラリ（{file}）が不一致（内容が変わったか、ファイルが"
                      "無い）なので、読んでいない。library_load の file_path で明示できる。")
    if status == "lbm_only":
        return None, ("データセットの同定は LBM（脂質のルールベース）だけで、MS/MS の照合に使う"
                      "MSP のライブラリが無いので、読んでいない（library_* は not_evaluable を返す）。"
                      "別のライブラリで照合するなら library_load の file_path で明示する。")
    if status != "known" or not record.get("path"):
        candidates = record.get("candidates") or []
        hint = f"（候補: {', '.join(candidates)}）" if candidates else ""
        return None, ("データセットが同定に使ったライブラリを特定できないので、読んでいない"
                      f"{hint}。極性の設定からは引かない（`library_load(dataset_handle=…)` を `file_path` なしで"
                      "呼ぶと、記憶が無いときだけ極性の設定から引く）。`file_path` で明示もできる。")
    path = Path(record["path"])
    try:
        store = open_store(path, read_only=True)  # 同じ構築先の同時構築は open_store が直列にする
    except Exception as exc:  # noqa: BLE001 - 壊れたライブラリでもハンドルは開く
        # 例外の文言は置き場所（フルパス）を含みうるので、種類だけを返す。
        return None, (f"ライブラリ（{path.name}）を開けなかったので、読んでいない"
                      f"（{type(exc).__name__}）。")
    store.library_file_name = path.name
    return store, None


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_handle_open(mztab_path: str | None = None, job_path: str | None = None,
                        pipeline_path: str | None = None, allow_incomplete: bool = False,
                        library: str = "auto", arf2_path: str | None = None,
                        directory: str | None = None, source: str = "auto",
                        measure: str = "peak_height") -> str:
    """データセットを**ハンドル**として開く（再解析で使う。対話のセッションは変えない）。

    入口と `source`・`measure`・`allow_incomplete` の意味は `dataset_load` と同じ
    （`directory` / `mztab_path` / `job_path` / `pipeline_path` のどれか 1 つ。`job_path` /
    `pipeline_path` の既定は解析フォルダで、mzTab-M で読むなら `source="mztab"`。pipeline の
    run を解析フォルダで読むと、run から引き継ぐのは解決済みの試料対応表だけで、解析行列・
    差次的結果は `source="mztab"` のときだけ載る）。`measure` が効くのは解析フォルダの
    アダプタだけで、mzTab-M のアダプタでは定量種別がファイル（job / pipeline ならジョブの
    宣言）で決まるので、既定値以外は `DATASET_BAD_REQUEST`。`allow_incomplete` も bool
    以外は `DATASET_BAD_REQUEST`。返った `dataset_handle` を
    `dataset_*`・`library_*`・再解析のツールの `dataset_handle` に渡すと、そのハンドルの
    データセットとライブラリで動く。ハンドルを渡した呼び出しはサーバの別スレッドで回り、
    同じハンドルへの呼び出しは 1 本ずつ、別のハンドルどうしは並行に進む。

    library: `"auto"`（既定）はデータセットが同定に使ったライブラリの記憶
        （`library_record`）を読み取り専用で開く。実体の sha256 を測り直し、変わっている・
        無いとき（`mismatch`）、記憶が決められないとき（`unknown`）、LBM だけのとき
        （`lbm_only`）は開かず、`library_note` に理由。**極性の設定
        （`[library] msp_positive` / `msp_negative`）には落ちない**。別のライブラリは
        `library_load` の `file_path` で明示する。`"none"` は読まない。
    arf2_path: スポットの照合（spot_id = SMF_ID）に使う `.arf2`。省略時、解析フォルダで
        読んだならその `.arf2` と `.dcl`。mzTab-M で読んだなら検出状態を取り込んだ `.arf` の
        同じバッチ、無ければ mzTab-M の隣、`pipeline_path` なら run の
        `input/AlignResult-*.arf2`（複数あれば名前の新しいもの）を探す。

    同じ入力（読み込みが読んだ・探したファイルの size・更新時刻が同じ）を開き直すと
    読み込み済みのものを使い回す
    （`reused_base`）。それでもハンドルどうしの状態は別々。開けるのは
    `LIPIDMIX_MAX_DATASET_HANDLES`（既定 8）本まで。使い終えたら `dataset_handle_close`。
    """
    if library not in _LIBRARY_MODES:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"library は {' / '.join(_LIBRARY_MODES)} のどちらかです: {library!r}",
                           {"field": "library"})

    from metabolomix.dataset.library_memory import public_view
    from metabolomix.mztab.polarity import dataset_polarity

    try:
        handles.registry.check_capacity()  # 長い読み込みの前に断る（add でもう一度確かめる）
        opened = opening.open_dataset(directory=directory, mztab_path=mztab_path,
                                      job_path=job_path, pipeline_path=pipeline_path,
                                      source=source, measure=measure,
                                      allow_incomplete=allow_incomplete)
        ds = opened.dataset
        polarity, polarity_source, polarity_notes = dataset_polarity(ds)
        selection = opened.selection
        arf2, dcl = alignment_files_for(opened, arf2_path)
        handle = handles.DatasetHandle(
            ds, polarity=polarity, polarity_source=polarity_source,
            source={"kind": opened.kind, "path": opened.path,
                    "allow_incomplete": bool(allow_incomplete),
                    "adapter": opened.adapter, "measure": measure},
            arf2_path=arf2, dcl_path=dcl)
        store, note = (None, "library=\"none\" なので読んでいない。") if library == "none" \
            else _open_library_for(ds)
        handle.library.store = store
        # 対話のセッションの source_path はフルパスだが、ハンドルはファイル名だけを持つ
        # （研究室のライブラリの置き場所を戻り値・記録に出さない）。
        handle.library.source_path = store.library_file_name if store is not None else None
        try:
            handles.registry.add(handle)
        except DomainError:
            handle.close()
            raise
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    library_info = None
    if store is not None:
        summary = store.summary()
        library_info = {"file": handle.library.source_path,
                        "sha256": summary["source_sha256"],
                        "record_count": summary["record_count"]}
    detection = getattr(ds, "detected_mask", None) is not None
    notes = [*polarity_notes, *opened.notes]
    return json_payload({
        "status": "success",
        "dataset_handle": handle.handle_id,
        "adapter": opened.adapter,
        "polarity": polarity,
        "polarity_source": polarity_source,
        "source_format": ds.source_format,
        "n_features": len(ds.feature_ids),
        "n_samples": len(ds.sample_names),
        "detection_available": detection,
        "arf": selection.arf.name if selection is not None else None,
        "not_chosen": list(selection.not_chosen) if selection is not None else None,
        "library": library_info,
        **({"library_note": note} if note else {}),
        "library_record": public_view(ds.library_record),
        "arf2": Path(arf2).name if arf2 else None,
        "dcl": Path(dcl).name if dcl else None,
        "reused_base": opened.reused_base,
        "open_handles": len(handles.registry.list()),
        "max_handles": handles.max_handles(),
        **({"notes": notes} if notes else {}),
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_handle_list() -> str:
    """開いているハンドルの一覧（TSV: `dataset_handle current polarity source_kind source
    n_features n_samples library_sha256`）。古く開いた順。`current` は対話のセッションの直近の
    ハンドル（`dataset_handle` を省いた呼び出しが使うもの。`dataset_load` で読んだものは
    `source_kind` が `session`）。"""
    from metabolomix.arf2.reader import format_spots_as_table

    current = session_state.session.current_handle_id
    rows = []
    for handle in handles.registry.list():
        store = handle.library.store
        rows.append({
            "dataset_handle": handle.handle_id,
            "current": handle.handle_id == current,
            "polarity": handle.polarity,
            "source_kind": handle.source.get("kind"),
            "source": Path(handle.source.get("path") or "").name,
            "n_features": len(handle.dataset.feature_ids),
            "n_samples": len(handle.dataset.sample_names),
            "library_sha256": store.summary()["source_sha256"][:16] if store is not None else None,
        })
    columns = ["dataset_handle", "current", "polarity", "source_kind", "source", "n_features",
               "n_samples", "library_sha256"]
    return json_payload({
        "n_open": len(rows),
        "max_handles": handles.max_handles(),
        "handles": format_spots_as_table(rows, columns=columns) if rows else "",
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=True,
                                      idempotentHint=True),
          structured_output=False)
def dataset_handle_close(dataset_handle: str) -> str:
    """ハンドルを閉じる（データセットの複製とライブラリの接続を手放す）。

    対話のセッションの直近のハンドル（`dataset_handle_list` の `current`）を閉じると、
    `dataset_handle` を省いた呼び出しはデータセットが無い状態（`missing_state`）になる。
    開いている別のハンドルへは落ちない（もう一度 `dataset_load` で読む）。
    """
    try:
        handles.registry.close(dataset_handle)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    return json_payload({"status": "success", "closed": dataset_handle,
                         "open_handles": len(handles.registry.list())})
