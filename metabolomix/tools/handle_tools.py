"""データセットの入口 `load_dataset` と、開いているハンドルの一覧・クローズ。

読み込んだデータセットはすべてハンドルの登録簿（`core.handles.registry`）に住み、`load_dataset`
は新しいハンドルを足して**いつも直近にする**（spec 2026-11-06 §5・B-2b 決めたこと 2）。
`dataset_handle` を省いたツールの呼び出しは直近のハンドルに効く。再解析のエージェントは、返った
`dataset_handle` をツールに渡して使い、直近を当てにしない。同じ入力を開き直したときは読み込み
済みの基を使い回し、ハンドルには複製を渡す（`dataset.opening.open_dataset` が持つ）。

ライブラリはデータセットが同定に使ったものの記憶（`library_record`）から、読み取り専用で開く。
極性の設定（`[library] msp_positive` / `msp_negative`）には落ちない。戻り値にはファイル名と
sha256 だけを出し、パスは出さない。
"""
from __future__ import annotations

import re
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import handles, session_state, user_config
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.dataset import opening

__all__ = ["load_dataset", "dataset_list", "dataset_close"]

_LIBRARY_MODES = ("auto", "none")
_ARF_SUFFIX_RE = re.compile(r"_(?:PeakProperties|DriftSopts|DriftSpots)\.arf$", re.IGNORECASE)


def alignment_files_for(opened, arf2_path: str | None = None) -> tuple[str | None, str | None]:
    """`open_dataset` が開いたデータセットのスポットの照合に使う `.arf2` と `.dcl`。

    `load_dataset` が決め、ハンドルに `arf2_path` / `dcl_path` として持たせる。ハンドル ID を
    明示した呼び出しはこの値を使う。
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
        lookup = curation_apply.arf2_lookup_kwargs(ds)
        arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p, **lookup)
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
    except user_config.ConfigInvalidError as exc:  # キャッシュの置き場所（[paths] library_cache）が読めない
        return None, (f"設定ファイルが読めないので、ライブラリ（{path.name}）を開いていない"
                      f"（CONFIG_INVALID: {exc.message}）。")
    except Exception as exc:  # noqa: BLE001 - 壊れたライブラリでもハンドルは開く
        # 例外の文言は置き場所（フルパス）を含みうるので、種類だけを返す。
        return None, (f"ライブラリ（{path.name}）を開けなかったので、読んでいない"
                      f"（{type(exc).__name__}）。")
    store.library_file_name = path.name
    return store, None


def _max_handles_or_none() -> int | None:
    """要約に載せる上限。開いた直後に設定ファイルが壊れていたら None（開いたハンドルは返す）。"""
    try:
        return handles.max_handles()
    except DomainError:
        return None


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def load_dataset(directory: str | None = None, mztab_path: str | None = None,
                 job_path: str | None = None, pipeline_path: str | None = None,
                 allow_incomplete: bool = False, source: str = "auto",
                 measure: str = "peak_height", arf2_path: str | None = None,
                 library: str = "auto") -> str:
    """データセットを読み込み、新しいハンドル（`dh_…`）として開いて**直近にする**。入口はこの 1 本。

    入口（**どれか 1 つ**。全省略・2 つ以上の指定はエラー）:
      - directory: MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata`）。フォルダの中の
        特定の `.arf` / `.arf2` を渡してもよい。同じフォルダに `DriftSpots.arf` と
        `PeakProperties.arf` があれば `PeakProperties.arf`、複数のバッチがあれば最新を選び、
        選ばなかった `.arf` を戻り値に並べる。
      - mztab_path: .mzTab ファイルへの絶対パス。
      - job_path: analysis-job.json へのパス。ジョブの宣言（polarity・measure）で正準
        mzTab-M を一意に選ぶ。
      - pipeline_path: `pipeline_run` が返した pipeline-run.json（またはその親ディレクトリ）。

    source: どのアダプタで読むか。`"auto"`（既定）/ `"folder"` / `"mztab"`。
        - `directory` は解析フォルダ（`"mztab"` はエラー）、`mztab_path` は mzTab-M
          （`"folder"` はエラー）。
        - `job_path` / `pipeline_path` は既定（`"auto"` / `"folder"`）で、そのジョブの正準
          mzTab-M と同じアラインメントの**解析フォルダ**を読む。候補（mzTab-M の隣の `.arf`。
          pipeline はどれも合わなければ run の `input/AlignResult-*.arf2` と同じバッチも）を
          mzTab-M と特徴数・m/z・試料名で照合し、合うものがちょうど 1 本ならそれを読む。
          無ければ `ALIGNMENT_NOT_FOUND`、2 本以上なら `ALIGNMENT_AMBIGUOUS`（どちらも
          候補のファイル名を `details.candidates` に載せる）で止まり、mzTab-M へ黙って
          切り替えない。`"mztab"` を渡すと mzTab-M を読む。
        - **pipeline の run を解析フォルダで読んだときに run から引き継ぐのは、解決済みの
          試料対応表だけ**。run の解析行列・feature binding・差次的結果は
          `source="mztab"` のときだけ載せる——その run の mzTab に加えて、解決済みの
          試料対応表・feature binding・解析行列を**同じ run からまとめて**載せる
          （別 run の行列を混ぜると群の対応が黙ってずれるため）。載せた行列は
          `dataset_statistic` へ `matrix_id` で渡す。run の差次的結果も、データと群の指紋が
          今のデータセットと一致するものだけを戻す（result_id は run の TSV と同じ。
          直近の結果にはしないので result_id で名指しして使う）。戻せなかったものは理由つきで
          `summary` に並ぶ。
    measure: 解析フォルダのアダプタが値に使う `.arf` の欄（`"peak_height"`（既定）/
        `"peak_area_above_zero"`）。**この引数が効くのは解析フォルダのアダプタだけ**。
        mzTab-M のアダプタでは定量種別がファイル（job / pipeline ならジョブの宣言）で
        決まるので、既定値以外を渡すと `DATASET_BAD_REQUEST`。
    allow_incomplete: **完了していない実行**（partial / failed / running）の出力を
        読み込む（既定 False。bool 以外は `DATASET_BAD_REQUEST`——`"false"` という
        文字列を真と読まないため）。中断時点の生成物なので、既定では拒否する
        （`INCOMPLETE_ANALYSIS_JOB`）。True で読んだデータセットは探索専用になり、
        2 群比較（dataset_differential）と差次的エクスポートは拒否される——欠けた検体を
        「その群には無い」と読み違えたまま結論が出てしまうため。
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

    成功すると、いつも**新しいハンドルを直近にする**。`dataset_handle` を省いた `dataset_*` /
    `library_*` はこの直近のハンドルに効く。再解析のエージェントは返った `dataset_handle` を
    各ツールに渡して使い、直近を当てにしない（渡した呼び出しはサーバの別スレッドで回り、同じ
    ハンドルへの呼び出しは 1 本ずつ、別のハンドルどうしは並行に進む）。前のハンドルは閉じずに
    残る（`dataset_list` / `dataset_close`）。開いているハンドルが上限（`lipidmix.local.toml` の
    `[server] max_dataset_handles`、既定 8）なら `DATASET_HANDLE_LIMIT` で止め、古いハンドルを
    黙って閉じない（直近も動かさない）。同じ入力（読み込みが読んだ・探したファイルの size・
    更新時刻が同じ）を開き直すと読み込み済みのものを使い回す（`reused_base`）。それでも
    ハンドルどうしの状態は別々。

    戻り値: `dataset_handle`・読んだアダプタ・選んだ `.arf`・ライブラリの記録（ファイル名と
    sha256。置き場所は出さない）に、`summary`（特徴・試料の数、検出状態、同定の内訳、読み込みの
    注記）、`sample_counts`（群と役割ごとの試料数）、`next`（次の一手の例）。`.arf2` の概観や
    前処理前の PCA は返さない。
    """
    if library not in _LIBRARY_MODES:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"library は {' / '.join(_LIBRARY_MODES)} のどちらかです: {library!r}",
                           {"field": "library"})

    from metabolomix.dataset.library_memory import public_view
    from metabolomix.mztab.polarity import dataset_polarity
    from metabolomix.core import version
    from metabolomix.tools import mztab_tools

    try:
        handles.registry.check_capacity()  # 長い読み込みの前に断る（adopt_handle でもう一度確かめる）
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
        # ライブラリはファイル名だけを持つ（研究室のライブラリの置き場所を戻り値・記録に出さない）。
        handle.library.source_path = store.library_file_name if store is not None else None
        try:
            # 直近にする。`library_load` を先に呼んでいて、このデータセットのライブラリが
            # 決まらなかった（`"auto"`）ときは、セッションのライブラリをこのハンドルへ移す。
            session_state.session.adopt_handle(handle, move_library=(library == "auto"))
        except DomainError:
            handle.close()
            raise
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    library_info = None
    if handle.library.store is not None:
        summary = handle.library.store.summary()
        library_info = {"file": handle.library.source_path,
                        "sha256": summary["source_sha256"],
                        "record_count": summary["record_count"]}
    detection = getattr(ds, "detected_mask", None) is not None
    notes = [*polarity_notes, *opened.notes]
    return json_payload({
        "status": "success",
        # 配布先のクローンが origin/main より遅れているときだけ載せる（最新なら鍵ごと出さない）。
        **({"update_available": update} if (update := version.update_status()) else {}),
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
        "max_handles": _max_handles_or_none(),
        "summary": mztab_tools.load_summary_lines(ds, opened),
        "sample_counts": mztab_tools.sample_group_role_counts(ds),
        **({"notes": notes} if notes else {}),
        "next": list(_NEXT_STEPS),
    })


#: 読み込みの後の次の一手の例（ツール名は公開面のもの）。
_NEXT_STEPS = (
    "dataset_status() — 試料の群・役割と検出状態を確かめる（差次的解析の群名はここに出る）",
    "curation_review() — 注釈の一覧確認（解析フォルダで読んだデータセットのとき）",
    "dataset_preprocess() → dataset_pca() — 前処理して主成分分析で全体像を見る",
)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_list() -> str:
    """開いているハンドルの一覧（TSV: `dataset_handle current polarity source_kind source
    n_features n_samples library_sha256`）。古く開いた順。`current` は対話のセッションの直近の
    ハンドル（`dataset_handle` を省いた呼び出しが使うもの。`load_dataset` が開いたものが直近）。"""

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
    try:
        limit = handles.max_handles()
    except DomainError as exc:  # 設定ファイルが読めない（CONFIG_INVALID）
        return mztab_error(exc.code, exc.message, exc.details or None)
    return json_payload({
        "n_open": len(rows),
        "max_handles": limit,
        "handles": format_spots_as_table(rows, columns=columns) if rows else "",
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=True,
                                      idempotentHint=True),
          structured_output=False)
def dataset_close(dataset_handle: str) -> str:
    """ハンドルを閉じる（データセットの複製とライブラリの接続を手放す）。

    対話のセッションの直近のハンドル（`dataset_list` の `current`）を閉じると、
    `dataset_handle` を省いた呼び出しはデータセットが無い状態（`missing_state`）になる。
    開いている別のハンドルへは落ちない（もう一度 `load_dataset` で読む）。
    """
    try:
        handles.registry.close(dataset_handle)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    return json_payload({"status": "success", "closed": dataset_handle,
                         "open_handles": len(handles.registry.list())})
