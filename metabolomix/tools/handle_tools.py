"""データセットのハンドルを開閉する MCP ツール（spec 2026-10-22 §4 の 1）。

**再解析で使う**（並行のエージェントが、別々のデータセット・極性を互いに干渉せずに扱う）。
対話の使い方（`dataset_load` → `dataset_*`）は変わらない——ハンドルを渡さない呼び出しは
今までどおり対話のセッションに効く。

`dataset_handle_open` は `dataset_load` と同じ読み込み関数を呼ぶが、結果を対話のセッションに
置かず、ハンドルの登録簿（`core.handles.registry`）に置く。同じ入力を開き直したときは
読み込み済みの基を使い回し、ハンドルには複製を渡す。

ライブラリは極性で設定（`[library] msp_positive` / `msp_negative`。環境変数
`MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）だけから引き（`library_mode=msp_only` と同じ規則）、
読み取り専用で開く。戻り値にはファイル名と sha256 だけを出し、パスは出さない。
"""
from __future__ import annotations

import copy
import re
import threading
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import handles, user_config
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.serialization import json_payload

__all__ = ["dataset_handle_open", "dataset_handle_list", "dataset_handle_close"]

_LIBRARY_MODES = ("auto", "none")
#: ライブラリの store の構築（初回だけ数分）を 1 本ずつにする。
_BUILD_LOCK = threading.Lock()
_ARF_SUFFIX_RE = re.compile(r"_(?:PeakProperties|DriftSopts|DriftSpots)\.arf$", re.IGNORECASE)


def _signature(ds, input_path: Path) -> tuple:
    """入力の指紋（読んだファイルの size と mtime）。中身が変われば基を使い回さない。"""
    paths = sorted({str(input_path), *(getattr(ds, "source_files", None) or {})})
    out = []
    for path in paths:
        try:
            stat = Path(path).stat()
            out.append((path, stat.st_size, stat.st_mtime_ns))
        except OSError:
            out.append((path, None, None))
    return tuple(out)


def _load(kind: str, path: str, allow_incomplete: bool):
    if kind == "pipeline_path":
        from metabolomix.pipeline import session_handoff
        return session_handoff.load_run_dataset(path, allow_incomplete=allow_incomplete)
    from metabolomix.mztab.loading import load_dataset_state
    return load_dataset_state(
        mztab_path=Path(path) if kind == "mztab_path" else None,
        job_path=Path(path) if kind == "job_path" else None,
        allow_incomplete=allow_incomplete)


def _pipeline_arf2(pipeline_path: str) -> Path | None:
    """pipeline の run の `input/AlignResult-*.arf2`。無ければ None。

    run では mzTab が `console/attempt-*/msdial/` の下、`.arf`・`.arf2`・`.dcl` が
    `input/` にあるので、mzTab の隣を探す既定の経路では見つからない。`AlignResult-*.arf2`
    が複数あるときは名前が最も新しい（降順で先頭の）ものを選ぶ（決定的にするため。
    `AlignResult-<日時など>` は名前の順が新しさの順になる）。
    """
    candidates = sorted((Path(pipeline_path) / "input").glob("AlignResult-*.arf2"),
                        key=lambda p: p.name, reverse=True)
    return next((p for p in candidates if p.is_file()), None)


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
        arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p)
                                  for p in (getattr(ds, "source_files", None) or {}))), None)
    if arf2 is None and kind == "pipeline_path" and path:
        arf2 = _pipeline_arf2(path)
    if arf2 is None or not arf2.is_file():
        return None, None
    dcl = alignment_spot.alignment_sibling(arf2, ".dcl")
    return str(arf2), (str(dcl) if dcl else None)


def _open_library(polarity: str | None) -> tuple[object | None, str | None]:
    """(store, note)。設定に無い・極性が決まらないときは store なしで理由を返す。"""
    if polarity is None:
        return None, ("データセットの極性が決まらないので、ライブラリを読んでいない"
                      "（逆の極性のライブラリを掴まないため）。")
    key = user_config.MSP_SETTING_KEYS[polarity]
    setting = user_config.get_setting(key)
    if setting is None:
        return None, (f"[library] {key.split('.', 1)[1]} が設定されていないので、ライブラリを"
                      "読んでいない（library_* はこのハンドルでは missing_state を返す）。")
    from metabolomix.library.store import open_store

    path = Path(setting.value)
    if not path.is_file():
        return None, (f"{user_config.setting_label(setting)} のライブラリ（{path.name}）が"
                      "見つからないので、読んでいない。")
    try:
        with _BUILD_LOCK:
            store = open_store(path, read_only=True)
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
                        library: str = "auto", arf2_path: str | None = None) -> str:
    """データセットを**ハンドル**として開く（再解析で使う。対話のセッションは変えない）。

    引数の意味は `dataset_load` と同じ（3 つのパスのうち 1 つ）。返った `dataset_handle` を
    `dataset_*`・`library_*`・再解析のツールの `dataset_handle` に渡すと、そのハンドルの
    データセットとライブラリで動く。ハンドルを渡した呼び出しはサーバの別スレッドで回り、
    同じハンドルへの呼び出しは 1 本ずつ、別のハンドルどうしは並行に進む。

    library: `"auto"`（既定）はデータセットの極性の参照ライブラリを設定
        （`[library] msp_positive` / `msp_negative`）だけから引き、読み取り専用で開く。
        設定に無い・極性が決まらないときはライブラリなしで開き、`library_note` に理由。
        `"none"` は読まない。
    arf2_path: スポットの照合（spot_id = SMF_ID）に使う `.arf2`。省略時は検出状態を取り込んだ
        `.arf` の同じバッチ、無ければ mzTab-M の隣、`pipeline_path` なら run の
        `input/AlignResult-*.arf2`（複数あれば名前の新しいもの）を探す。

    同じ入力（ファイルの size・更新時刻が同じ）を開き直すと読み込み済みのものを使い回す
    （`reused_base`）。それでもハンドルどうしの状態は別々。開けるのは
    `LIPIDMIX_MAX_DATASET_HANDLES`（既定 8）本まで。使い終えたら `dataset_handle_close`。
    """
    given = {k: v for k, v in (("mztab_path", mztab_path), ("job_path", job_path),
                               ("pipeline_path", pipeline_path)) if v}
    if len(given) != 1:
        return mztab_error("DATASET_BAD_REQUEST",
                           "mztab_path / job_path / pipeline_path のどれか 1 つを指定してください。",
                           {"given": sorted(given)})
    if library not in _LIBRARY_MODES:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"library は {' / '.join(_LIBRARY_MODES)} のどちらかです: {library!r}",
                           {"field": "library"})
    kind, path = next(iter(given.items()))
    input_path = Path(path).resolve()
    key = canonical_hash({"kind": kind, "path": str(input_path),
                          "allow_incomplete": bool(allow_incomplete)})

    from metabolomix.mztab.polarity import dataset_polarity

    try:
        cached = handles.registry.cached_base(key)
        reused = cached is not None and _signature(cached[0], input_path) == cached[1]
        if reused:
            base = cached[0]
        else:
            base = _load(kind, path, bool(allow_incomplete))
            handles.registry.put_base(key, base, _signature(base, input_path))
        ds = copy.deepcopy(base)
        polarity, polarity_source, polarity_notes = dataset_polarity(ds)
        arf2, dcl = _alignment_files(ds, arf2_path, kind, path)
        handle = handles.DatasetHandle(
            ds, polarity=polarity, polarity_source=polarity_source,
            source={"kind": kind, "path": str(input_path),
                    "allow_incomplete": bool(allow_incomplete)},
            arf2_path=arf2, dcl_path=dcl)
        store, note = (None, "library=\"none\" なので読んでいない。") if library == "none" \
            else _open_library(polarity)
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
    return json_payload({
        "status": "success",
        "dataset_handle": handle.handle_id,
        "polarity": polarity,
        "polarity_source": polarity_source,
        "source_format": ds.source_format,
        "n_features": len(ds.feature_ids),
        "n_samples": len(ds.sample_names),
        "detection_available": detection,
        "library": library_info,
        **({"library_note": note} if note else {}),
        "arf2": Path(arf2).name if arf2 else None,
        "dcl": Path(dcl).name if dcl else None,
        "reused_base": reused,
        "open_handles": len(handles.registry.list()),
        "max_handles": handles.max_handles(),
        **({"notes": polarity_notes} if polarity_notes else {}),
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_handle_list() -> str:
    """開いているハンドルの一覧（TSV: `dataset_handle polarity source_kind source n_features
    n_samples library_sha256`）。古く開いた順。"""
    from metabolomix.arf2.reader import format_spots_as_table

    rows = []
    for handle in handles.registry.list():
        store = handle.library.store
        rows.append({
            "dataset_handle": handle.handle_id,
            "polarity": handle.polarity,
            "source_kind": handle.source.get("kind"),
            "source": Path(handle.source.get("path") or "").name,
            "n_features": len(handle.dataset.feature_ids),
            "n_samples": len(handle.dataset.sample_names),
            "library_sha256": store.summary()["source_sha256"][:16] if store is not None else None,
        })
    columns = ["dataset_handle", "polarity", "source_kind", "source", "n_features",
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
    """ハンドルを閉じる（データセットの複製とライブラリの接続を手放す）。対話のセッションは変えない。"""
    try:
        handles.registry.close(dataset_handle)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    return json_payload({"status": "success", "closed": dataset_handle,
                         "open_handles": len(handles.registry.list())})
