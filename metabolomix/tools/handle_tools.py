"""データセットの入口 `load_dataset` と、開いているハンドルの一覧・クローズ。

読み込んだデータセットはすべてハンドルの登録簿（`core.handles.registry`）に住み、`load_dataset`
は新しいハンドルを足して**いつも直近にする**（spec 2026-11-06 §5・B-2b 決めたこと 2）。
`dataset_handle` を省いたツールの呼び出しは直近のハンドルに効く。再解析のエージェントは、返った
`dataset_handle` をツールに渡して使い、直近を当てにしない。同じ入力を開き直したときは読み込み
済みの基を使い回し、ハンドルには複製を渡す（`dataset.opening.open_dataset` が持つ）。

ライブラリはデータセットが同定に使ったものの記憶（`library_record`）から、読み取り専用で開く。
極性の設定（`[library] msp_positive` / `msp_negative`）には落ちない。戻り値にはファイル名と
sha256 だけを出し、パスは出さない。

読み込みは受付（`core.pending_loads`）のスレッドで回す。`wait_seconds` の中で終われば戻り値を
そのまま返し、終わらなければ受付番号（`ld_…`）を返して `dataset_load_status` で引かせる
（デスクトップアプリの中継が 36〜60 秒ほどで呼び出しを打ち切るため）。本体は `_load_now`。
長い道具（`handle_dispatch.LONG_RUNNING`）の受付番号（`job_…`）は `job_status` で引く。
"""
from __future__ import annotations

import functools
import re
import threading
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import (handle_dispatch, handles, mcp_core, pending_loads, session_state,
                              user_config)
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error, needs_input
from metabolomix.core.serialization import json_payload, json_safe
from metabolomix.dataset import opening

__all__ = ["load_dataset", "dataset_load_status", "job_status", "dataset_list", "dataset_close"]

_LIBRARY_MODES = ("auto", "none")
_ARF_SUFFIX_RE = re.compile(r"_(?:PeakProperties|DriftSopts|DriftSpots)\.arf$", re.IGNORECASE)


def alignment_files_for(opened, arf2_path: str | None = None) -> tuple[str | None, str | None]:
    """`open_dataset` が開いたデータセットのスポットの照合に使う `.arf2` と `.dcl`。

    `load_dataset` が決め、ハンドルに `arf2_path` / `dcl_path` として持たせる。ハンドル ID を
    明示した呼び出しはこの値を使う。
    ID を省いた呼び出しも、直近のハンドルのこの値を使う（`_holder.alignment_paths`。
    データフォルダは探さない）。
    解析フォルダで読んだなら（明示が無ければ）その選択、そうでなければ `_alignment_files` が探す。
    """
    selection = opened.selection
    if selection is not None and not arf2_path:
        return str(selection.arf2), (str(selection.dcl) if selection.dcl else None)
    return _alignment_files(opened.dataset, arf2_path, opened.kind, opened.path, opened.polarity)


def _alignment_files(ds, arf2_path: str | None, kind: str | None = None,
                     path: str | None = None,
                     polarity: str | None = None) -> tuple[str | None, str | None]:
    """スポットの照合に使う `.arf2` と `.dcl`。決まらなければ None。

    探す順: 明示の `arf2_path` → 検出状態を取り込んだ `.arf` の同じバッチ → mzTab の隣
    → （`kind == "pipeline_path"` のとき）受領が照合した単位のバッチ、無ければ単位の
    `input/*.arf2` の名前が最も新しいもの（`opening.pipeline_arf2`）。
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
        arf2 = opening.pipeline_arf2(path, polarity)
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
                      f"{hint}。極性の設定からは引かない（`library_load` も `file_path` なしでは"
                      "開けない）。`file_path` で明示できる。")
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
                 polarity: str | None = None,
                 allow_incomplete: bool = False, source: str = "auto",
                 measure: str = "peak_height", arf2_path: str | None = None,
                 library: str = "auto", overview: bool = False,
                 wait_seconds: float = 20.0) -> str:
    """データセットを読み込み、新しいハンドル（`dh_…`）として開いて直近にする。

    入口（どれか 1 つ。全省略・2 つ以上は誤り）:
      - directory: MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata`）か、その中の `.arf` / `.arf2`。
        `DriftSpots.arf` と `PeakProperties.arf` があれば `PeakProperties.arf`、複数のバッチがあれば最新を選び、
        選ばなかった `.arf` を戻り値に並べる。
      - mztab_path: .mzTab ファイルの絶対パス。
      - job_path: analysis-job.json。ジョブの宣言（polarity・measure）で正準の mzTab-M を 1 つに選ぶ。
      - pipeline_path: `pipeline_run` が返した pipeline-run.json（かその親フォルダ）。受領
        （`results/receipt.json`）が指す極性ごとの解析フォルダを開く。

    polarity: `pipeline_path` のときだけ。開く極性（`"negative"` / `"positive"`）。1 極性の run では省ける。
        2 極性の run で省くと `needs_input`（検証を通った極性から選ぶ）。検証を通っていない単位は
        `PIPELINE_UNIT_NOT_VALIDATED`（`allow_incomplete=True` なら探索用に読む）。ほかの入口で渡すと
        `PIPELINE_POLARITY_NOT_APPLICABLE`。2 極性の run はハンドルが 2 つになる。
    source: `"auto"`（既定）/ `"folder"` / `"mztab"`。
        - `directory` は解析フォルダ（`"mztab"` は誤り）、`mztab_path` は mzTab-M（`"folder"` は誤り）。
        - `job_path` / `pipeline_path` は既定（`"auto"` / `"folder"`）で、正準 mzTab-M と同じアラインメントの
          解析フォルダを読む。候補（mzTab-M の隣の `.arf`。pipeline はどれも合わなければ `<極性>/input/*.arf2`
          のバッチを新しい順に）を mzTab-M と特徴数・m/z・試料名で照合し、ちょうど 1 本ならそれ。0 本は
          `ALIGNMENT_NOT_FOUND`、2 本以上は `ALIGNMENT_AMBIGUOUS`（どちらも `details.candidates` にファイル名）で
          止まり、mzTab-M へ切り替えない。`"mztab"` なら mzTab-M を読む。
        - pipeline の run から引き継ぐのは試料の対応表だけ（v3 は受領の `units.<極性>.samples.rows`、v1 は記録の
          `inputs.manifest`。どちらのアダプタでも同じ）。run の差次的結果は戻さない（比較は `dataset_compare`）。
          受領が無い・照合できない v3 の run は記録の単位から解析フォルダを開き、対応表が載らない旨を注記する。
    measure: 解析フォルダのアダプタが値に使う `.arf` の欄（`"peak_height"`（既定）/ `"peak_area_above_zero"`）。
        mzTab-M のアダプタでは定量種別がファイル（job / pipeline はジョブの宣言）で決まり、既定以外は
        `DATASET_BAD_REQUEST`。
    allow_incomplete: 完了していない実行（partial / failed / running）の出力を読む（既定 False は
        `INCOMPLETE_ANALYSIS_JOB`。bool 以外は `DATASET_BAD_REQUEST`）。True で読んだデータセットは探索専用
        （検証を通っていない pipeline の単位も、`directory` で run の中のそうした解析フォルダを開いたときも同じ印。
        断らない）。比較も書き出しもできるが、差次的 TSV のメタ行（`# exploratory = true`）・根拠バンドルの
        manifest（`exploratory`）・図・最初の概観の比較の事実に印と理由が載る（欠けた検体を「その群に無い」と
        読み違えうるため）。
    library: `"auto"`（既定）はデータセットが同定に使ったライブラリの記憶（`library_record`）を読み取り専用で
        開く。sha256 を測り直して変わった・無い（`mismatch`）、決められない（`unknown`）、LBM だけ（`lbm_only`）
        なら開かず、理由を `library_note` に書く。極性の設定（`[library] msp_positive` / `msp_negative`）には
        落ちない。別のライブラリは `library_load(file_path=…)`。`"none"` は読まない。
    arf2_path: スポットの照合（spot_id = SMF_ID）に使う `.arf2`。省くと、解析フォルダならその `.arf2` と `.dcl`、
        mzTab-M なら検出状態を取り込んだ `.arf` と同じバッチ、無ければ mzTab-M の隣、`pipeline_path` なら受領が
        照合した `<極性>/input/` の `.arf2`、受領に無ければ `<極性>/input/*.arf2`（複数なら名前の新しいもの）。
    overview: True なら最後に最初の概観（`dataset_overview` と同じ）を計算して `overview` に載せる（既定 False。
        bool 以外は `DATASET_BAD_REQUEST`）。概観が失敗しても読み込みは成功し、`overview` に
        `status:"unavailable"` と理由が載る。
    wait_seconds: 待つ秒数（既定 20。0〜30。範囲の外は丸めて `wait_seconds_note`。NaN は `DATASET_BAD_REQUEST`。
        真偽値と数に読める文字列は MCP の引数検査が数に変え（true は 1 秒）、読めない値は封筒でない引数検査の
        誤り）。終わらなければ読み込みを続け、`{"status": "loading", "load_id": "ld_…",
        "next": "dataset_load_status(load_id=…)"}` を返す。0 は「終わっていなければすぐ受付番号」。
        ハンドルの登録と直近の切り替えは読み込みが終わったときに 1 回だけ（受付の後に失敗すれば登録しない。
        その間に別のハンドルが直近になっていたら、登録だけして直近にせずライブラリも移さず、`current: false` と
        `note`）。引数の誤りと上限（`DATASET_HANDLE_LIMIT`。読み込み中の受付も 1 本と数える）は受け付ける前に返す。

    成功すると新しいハンドルを直近にする（`dataset_handle` を省いた `dataset_*` / `library_*` が使う）。並行の
    エージェントは返った `dataset_handle` を渡す。前のハンドルは閉じずに残る（`dataset_list` / `dataset_close`）。
    開いているハンドルが上限（`lipidmix.local.toml` の `[server] max_dataset_handles`、既定 8）なら
    `DATASET_HANDLE_LIMIT`（古いものを閉じず、直近も動かさない）。同じ入力（読んだ・探したファイルの size と
    更新時刻が同じ）は読み込み済みを使い回す（`reused_base`。ハンドルの状態は別々）。
    `assay_kind` はハンドルごと（解析フォルダ配下 `reports/` の目的ファイルから読み戻す。無ければ `unknown` で、
    `record_objective(assay_kind=...)` で確定）。

    戻り値: `dataset_handle`・`assay_kind`・読んだアダプタ・選んだ `.arf`・ライブラリ（ファイル名と sha256）・
    `caveat`（意味論の要約。未読でその種別にまだ出していないときだけ。種別の確定後の最初の 1 回は
    `record_objective` が運ぶ）・`summary`・`sample_counts`・`overview`・`next`。pipeline の run の単位なら
    `processing`（Console に渡した実効メソッドの `key_params`・宣言の無い `undeclared`・上書きの `overrides`。
    変えて処理し直すのは `pipeline_run` の `request.method_overrides`）。概観の比較と PCA は `result_id`
    で掘れるが直近の結果にはしない。欄の意味、ハンドルの並行・待ち、種別の引き継ぎは
    `docs_read(topic="reanalysis", section="データセットのハンドル")`。
    """
    bad = _argument_error(overview, library)
    if bad is not None:
        return bad
    try:
        wait, wait_note = pending_loads.clamp_wait(wait_seconds)
        # 上限なら受け付けない。読み込み中の受付も 1 本と数える（登録で必ず断られる読み込みを
        # 何本も回さない）。登録（adopt_handle）でもう一度確かめる。
        reservation = handles.registry.reserve()
        ticket = _CurrentTicket(session_state.session.current_generation)
        try:
            # 受付番号を返すまでは、省いた dataset_load_status の候補に並べない（ここで結果を直に待つので）。
            pending = pending_loads.registry.submit(functools.partial(
                _load_reserved, reservation, ticket, directory=directory, mztab_path=mztab_path,
                job_path=job_path, pipeline_path=pipeline_path, polarity=polarity,
                allow_incomplete=allow_incomplete, source=source, measure=measure,
                arf2_path=arf2_path, library=library, overview=overview),
                tool="load_dataset",
                label=_input_label(directory, mztab_path, job_path, pipeline_path, polarity),
                announced=False)
        except BaseException:
            reservation.release()
            raise
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    # 受付そのものを待つ（受付簿から捨てられても結果を失わない）。
    payload = pending.wait(wait)
    if payload is not None:
        pending_loads.registry.discard(pending.load_id)  # 受付番号を誰にも渡していないので覚えない
        return pending_loads.with_wait_note(payload, wait_note)
    ticket.detach()  # ここから後に終わる読み込みは、受け付けた後の直近の動きを見る
    pending_loads.registry.announce(pending.load_id)
    return json_payload({"status": "loading", "load_id": pending.load_id,
                         **({"wait_seconds_note": wait_note} if wait_note else {}),
                         "next": f'dataset_load_status(load_id="{pending.load_id}", '
                                 f'wait_seconds={MAX_WAIT_SECONDS:g})'})


#: `wait_seconds` の上限（秒）。正準は `pending_loads.MAX_WAIT_SECONDS`。
MAX_WAIT_SECONDS = pending_loads.MAX_WAIT_SECONDS


def _input_name(path) -> str | None:
    """入力の名前だけ（パスは出さない）。`pipeline-run.json` はそれを置いた run のフォルダの名前。"""
    if not path:
        return None
    path = Path(str(path))
    return (path.parent.name if path.name.lower() == "pipeline-run.json" else path.name) or str(path)


def _input_label(directory, mztab_path, job_path, pipeline_path, polarity) -> str | None:
    """受付の一覧（`dataset_load_status` の選択肢）に見せる短い説明。入力の名前だけ（パスは出さない）。"""
    name = _input_name(directory or mztab_path or job_path or pipeline_path)
    if name is None:
        return None
    return f"{name} ({polarity})" if polarity else name


#: 受け付けた後に別のハンドルが直近になったため、直近にしなかったときの `note`。
_NOT_CURRENT_NOTE = ("受け付けた後に別のハンドルが直近になったので、このハンドルは直近にしていない。"
                     "dataset_handle を渡して使う。")


class _CurrentTicket:
    """受け付けた読み込みが、終わったときに直近になってよいか。

    受け付けたときの直近の世代（`session.current_generation`）を覚える。`load_dataset` が受付番号を
    返した（`detach`）後に終わった読み込みは、その間に別のハンドルが直近になっていれば（世代が
    進んでいれば）直近にしない。`detach` の前に登録まで済んだ読み込み（期限の中で終わったもの）は
    今までどおり直近にする。判定（`adopt_handle` が `_current_lock` の下で呼ぶ）と `detach` は
    このロックで順序を決める。
    """

    def __init__(self, generation: int):
        self._lock = threading.Lock()
        self._generation = generation
        self._detached = False

    def detach(self) -> None:
        with self._lock:
            self._detached = True

    def may_become_current(self, generation: int) -> bool:
        with self._lock:
            return not self._detached or generation == self._generation


def _load_reserved(reservation, ticket, **kwargs) -> str:
    """受付のスレッドで本体を回す。登録で使わなかった枠（失敗したとき）は必ず返す。"""
    try:
        return _load_now(ticket=ticket, reservation=reservation, **kwargs)
    finally:
        reservation.release()  # 登録で使った枠なら何もしない


def _argument_error(overview, library) -> str | None:
    """すぐに分かる引数の誤り（受け付ける前に返す）。無ければ None。"""
    if not isinstance(overview, bool):
        return mztab_error("DATASET_BAD_REQUEST", f"overview は true / false です: {overview!r}",
                           {"field": "overview"})
    if library not in _LIBRARY_MODES:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"library は {' / '.join(_LIBRARY_MODES)} のどちらかです: {library!r}",
                           {"field": "library"})
    return None


def _recall_assay_kind(handle) -> tuple[str | None, str | None]:
    """解析フォルダの目的ファイルに記録した種別をハンドルへ読み戻す。(読み戻した種別, 注記)。

    記録の置き場所が決まらない・設定ファイルが読めない・読めないファイルは、読み戻さずに `unknown`
    のまま続ける（読み込みは止めない）。食い違いで読み戻さなかったときは (None, 理由)。
    """
    from metabolomix.tools import objective

    try:
        kind, note = objective.recall_assay_kind(
            handle.dataset, mcp_core._report_dir_candidates(handle))
    except (DomainError, OSError):
        return None, None
    if kind is not None:
        handle.assay_kind = kind
    return kind, note


class _SessionKindCarry:
    """データセットを開く前に記録したセッションの種別を、直近になるハンドルへ当てる
    （`AnalysisSession.adopt_handle` の `take_session_assay_kind`。セッションのロックの下で呼ばれる）。

    目的ファイルから読み戻せた種別と同じなら何も足さない。読み戻した種別と食い違えば、どちらも採らずに
    `unknown` にする（目的ファイルどうしの食い違いと同じ扱い）。目的ファイルどうしが食い違って読み戻さな
    かったときも当てない。目的ファイルに記録が無ければ当てる。どの場合も `note` に 1 文を残す（同じなら無し）。

    当てた種別は、登録の後に `persist` が目的ファイルへ残す（その種別を書いた目的にこのハンドルの解析フォルダを
    書く。`objective.persist_carried_kind`）——残さないと、同じ解析フォルダを開き直したときに読み戻せない。
    """

    def __init__(self, handle, recalled: str | None, recall_note: str | None):
        self._handle = handle
        self._recalled = recalled
        self._recall_conflict = recalled is None and recall_note is not None
        self._carried: tuple[str, str | None] | None = None
        self.note: str | None = None

    def __call__(self, kind: str, analysis_id: str | None = None) -> None:
        if self._recalled == kind:
            return
        if self._recalled is not None:
            self._handle.assay_kind = "unknown"
            self.note = (f"目的ファイルの assay_kind={self._recalled} と、データセットを開く前に記録した "
                         f"assay_kind={kind} が食い違うので、どちらも採らず unknown にした"
                         "（record_objective(assay_kind=...) で確定する）。")
        elif self._recall_conflict:
            self.note = (f"データセットを開く前に記録した assay_kind={kind} も、目的ファイルが食い違うので"
                         "当てていない（unknown のまま）。")
        else:
            self._handle.assay_kind = kind
            self._carried = (kind, analysis_id)
            self.note = (f"assay_kind={kind} を、データセットを開く前の record_objective の記録から"
                         "このハンドルへ引き継いだ。")

    def persist(self) -> None:
        """当てた種別を目的ファイルに残す（登録の後に、セッションのロックの外で呼ぶ）。

        書く間はハンドルのロックを持つ（同じハンドルへの `record_objective` と重ならない）。登録の後なので
        ほかの呼び出しが先にロックを取りうる——待つのは `handle_dispatch.HANDLE_LOCK_WAIT` 秒まで。
        残せなかったときは `note` にそう書く（開き直すと unknown に戻るので、記録し直すよう促す）。
        """
        if self._carried is None:
            return
        from metabolomix.tools import objective

        kind, analysis_id = self._carried
        redo = ("同じ解析フォルダを開き直すと unknown に戻るので、record_objective(assay_kind=...) で"
                "記録し直す。")
        if analysis_id is None:
            self.note += f"目的ファイルには残していない（どの目的の種別か分からない）。{redo}"
            return
        if not self._handle.lock.acquire(timeout=handle_dispatch.HANDLE_LOCK_WAIT):
            self.note += f"目的ファイルに残せなかった（ハンドルが別の呼び出しで塞がっていた）。{redo}"
            return
        try:
            name, placed = objective.persist_carried_kind(self._handle, analysis_id, kind)
        except (DomainError, OSError, ValueError) as exc:
            self.note += f"目的ファイルに残せなかった（{type(exc).__name__}）。{redo}"
            return
        finally:
            self._handle.lock.release()
        self.note += (f"目的ファイル（{name}）にもこの解析フォルダの種別として残した"
                      "（開き直すと読み戻す）。" + (f"{placed}。" if placed else ""))


def _load_now(*, directory: str | None = None, mztab_path: str | None = None,
              job_path: str | None = None, pipeline_path: str | None = None,
              polarity: str | None = None, allow_incomplete: bool = False,
              source: str = "auto", measure: str = "peak_height",
              arf2_path: str | None = None, library: str = "auto",
              overview: bool = False,
              ticket: _CurrentTicket | None = None,
              reservation: handles.Reservation | None = None) -> str:
    """`load_dataset` の本体（受付のスレッドで回る）。戻り値は `load_dataset` の JSON の文字列
    （成功・誤りの封筒・`needs_input`）。成功ならハンドルを登録して直近にする（ここで 1 回だけ）。

    `ticket` があれば直近にするかはその判定に従い、直近にしなかったら戻り値に `current: false` と
    `note` を足す（ほかの欄は変えない）。`reservation` は受付で取った枠で、登録と同時に消える。"""
    bad = _argument_error(overview, library)
    if bad is not None:
        return bad

    from metabolomix.dataset.library_memory import public_view
    from metabolomix.mztab.polarity import dataset_polarity
    from metabolomix.core import version
    from metabolomix.tools import mztab_tools

    try:
        handles.registry.check_capacity()  # 長い読み込みの前に断る（adopt_handle でもう一度確かめる）
        opened = opening.open_dataset(directory=directory, mztab_path=mztab_path,
                                      job_path=job_path, pipeline_path=pipeline_path,
                                      polarity=polarity, source=source, measure=measure,
                                      allow_incomplete=allow_incomplete)
        ds = opened.dataset
        polarity, polarity_source, polarity_notes = dataset_polarity(ds)
        selection = opened.selection
        arf2, dcl = alignment_files_for(opened, arf2_path)
        handle = handles.DatasetHandle(
            ds, polarity=polarity, polarity_source=polarity_source,
            source={"kind": opened.kind, "path": opened.path,
                    **({"polarity": opened.polarity} if opened.polarity else {}),
                    "allow_incomplete": bool(allow_incomplete),
                    "adapter": opened.adapter, "measure": measure},
            arf2_path=arf2, dcl_path=dcl)
        # 種別はハンドル（データセット）ごと。目的ファイルに記録があれば読み戻す（登録の前に）。
        recalled, kind_note = _recall_assay_kind(handle)
        # ハンドルが無いうちに記録したセッションの種別は、直近になるこのハンドルへ移す（adopt_handle）。
        carry = _SessionKindCarry(handle, recalled, kind_note)
        store, note = (None, "library=\"none\" なので読んでいない。") if library == "none" \
            else _open_library_for(ds)
        handle.library.store = store
        # ライブラリはファイル名だけを持つ（研究室のライブラリの置き場所を戻り値・記録に出さない）。
        handle.library.source_path = store.library_file_name if store is not None else None
        try:
            # 直近にする。`library_load` を先に呼んでいて、このデータセットのライブラリが
            # 決まらなかった（`"auto"`）ときは、セッションのライブラリをこのハンドルへ移す。
            # 受け付けた後に別のハンドルが直近になっていたら（`ticket`）、登録だけして移さない。
            made_current = session_state.session.adopt_handle(
                handle, move_library=(library == "auto"), reservation=reservation,
                may_become_current=ticket.may_become_current if ticket is not None else None,
                take_session_assay_kind=carry)
        except DomainError:
            handle.close()
            raise
        carry.persist()   # 引き継いだ種別を目的ファイルへ（登録の後。セッションのロックの外）
    except DomainError as exc:
        if exc.code == "PIPELINE_POLARITY_REQUIRED":
            return needs_input("polarity", exc.message, (exc.details or {})["options"],
                               argument="polarity", retry_tools=["load_dataset"])
        return mztab_error(exc.code, exc.message, exc.details or None)

    library_info = None
    if handle.library.store is not None:
        summary = handle.library.store.summary()
        # `library_load`（ハンドルなし）が先に読んだライブラリは source_path がフルパスなので、名前だけを返す。
        library_info = {"file": Path(handle.library.source_path).name,
                        "sha256": summary["source_sha256"],
                        "record_count": summary["record_count"]}
    detection = getattr(ds, "detected_mask", None) is not None
    notes = [*polarity_notes, *opened.notes,
             *(note for note in (kind_note, carry.note) if note)]
    # 意味論ダイジェスト（未読・種別ごとに 1 回・`caveat_mode=off` で止まる）。最上位の文字列として
    # 先頭に置く——結論より前に解釈の前提が目に入るように。無ければ鍵ごと出さない。
    # 種別はこのハンドルのもの（直近にしなかった読み込みでも、そのデータセットの規則を出す）。
    caveat = session_state.session.take_caveat(handle.assay_kind)
    first_overview = _first_overview(handle) if overview else None
    payload = {
        "status": "success",
        **({"caveat": caveat} if caveat else {}),
        # 配布先のクローンが origin/main より遅れているときだけ載せる（最新なら鍵ごと出さない）。
        **({"update_available": update} if (update := version.update_status()) else {}),
        "dataset_handle": handle.handle_id,
        **({} if made_current else {"current": False, "note": _NOT_CURRENT_NOTE}),
        "assay_kind": handle.assay_kind,
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
        # pipeline の run の単位なら、Console に渡した実効メソッドの主なキー（ピーク高・足切り・
        # 上書き）。変えて処理し直すのは pipeline_run の request.method_overrides。
        **({"processing": opened.method_params} if opened.method_params else {}),
        "dcl": Path(dcl).name if dcl else None,
        "reused_base": opened.reused_base,
        "open_handles": len(handles.registry.list()),
        "max_handles": _max_handles_or_none(),
        "summary": mztab_tools.load_summary_lines(ds, opened),
        "sample_counts": mztab_tools.sample_group_role_counts(ds),
        **({"overview": first_overview} if first_overview is not None else {}),
        **({"notes": notes} if notes else {}),
        "next": [step for step in _NEXT_STEPS
                 if first_overview is None or not step.startswith("dataset_overview()")],
    }
    return json_payload(payload)


def _first_overview(handle) -> dict:
    """最初の概観（`dataset_overview` と同じ組み立て。`counts` は summary と重なるので省く）。

    `load_dataset` は `handle_dispatch` のロックの外で動くので、新しいハンドルのロックの下で計算する
    （直後の別の呼び出しと重ならない）。どんな失敗でも読み込みは止めず、理由を返す。
    """
    from metabolomix.tools.overview_tools import build_overview_payload

    try:
        with handle.lock:
            return build_overview_payload(handle, embedded=True)
    except DomainError as exc:
        error = {"code": exc.code, "message": exc.message}
        if exc.details:
            error["details"] = json_safe(exc.details)
    except Exception as exc:  # noqa: BLE001 - 概観の失敗で読み込みを失敗にしない
        error = {"code": type(exc).__name__, "message": str(exc)[:300]}
    return {"status": "unavailable", "error": error,
            "next": "dataset_overview() を呼び直す／dataset_preprocess() で前処理の段の理由を見る"}


#: 読み込みの後の次の一手の例（ツール名は公開面のもの）。
_NEXT_STEPS = (
    "dataset_overview() — 最初の概観（品質 → 全体 → クラス → 分子種）。各段の next から掘る（鎖を"
    "進めた後は出し直す。2 群の向きは仮なので dataset_compare(reference_group=..., test_group=...) で決める）",
    "dataset_status() — 試料の群・役割と検出状態を確かめる（差次的解析の群名はここに出る）",
    "curation_review() — 注釈の一覧確認（解析フォルダで読んだデータセットのとき）",
)


def _pick_pending(book, given: str | None, *, field: str, tool: str,
                  start_tools: list[str], state: str,
                  foreign=None) -> tuple[str | None, str | None]:
    """受付番号を決める。(番号, 封筒)。省いたら、まだ結果を渡していない受付から選ぶ。

    並ぶのは受付番号を呼び出し側へ渡した受付だけ（`PendingCalls.unclaimed`）。ちょうど 1 本ならそれ、
    2 本以上なら `needs_input`（選べる受付つき）、無ければ `missing_state`。`foreign(entry)` が理由の文を
    返す 1 本（直近でないハンドルの仕事など。別のエージェントのものかもしれない）は黙って選ばず、
    その 1 本を選択肢にした `needs_input` にする（選ぶと結果を渡した印が付き、持ち主の省いた問い合わせ
    から消えるため）。
    """
    if given is not None:
        return given, None
    open_entries = book.unclaimed()
    reason = None
    if len(open_entries) == 1:
        reason = foreign(open_entries[0]) if foreign is not None else None
        if reason is None:
            return open_entries[0]["id"], None
    if not open_entries:
        return None, missing_state(
            state, start_tools,
            f"結果を渡していない受付がありません（{field} を渡すか、"
            f"{' / '.join(start_tools)} を先に呼ぶ）。")
    options = [{"value": e["id"],
                "label": " ".join(str(part) for part in (
                    e["tool"], e["label"], e["state"], f"{e['elapsed_s']}s") if part)}
               for e in open_entries]
    if reason is None:
        reason = f"結果を渡していない受付が {len(options)} 本あります。"
    return None, needs_input(field, f"{reason}どれを引くか {field} で選ぶ。", options,
                             argument=field, retry_tools=[tool])


def _poll(book, entry_id: str, wait_seconds, running, done_fields: dict | None = None) -> str:
    """受付を引く（最長 `wait_seconds` 秒待つ）。終わっていれば戻り値（渡した印を付ける。JSON の dict なら
    `done_fields` を先頭へ足す）、走っていれば `running`（dict か、待った後に dict を組む関数）に
    `elapsed_s` を足したもの。"""
    try:
        wait, wait_note = pending_loads.clamp_wait(wait_seconds)
        if wait > 0:
            # 受付そのものを待つ（`PendingCalls.wait`）。終わればその戻り値をそのまま返す。
            book.wait(entry_id, wait)
        state, payload, elapsed = book.snapshot(entry_id)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    if state == "done":
        book.claim(entry_id)
        payload = pending_loads.with_fields(payload, done_fields or {})
        return pending_loads.with_wait_note(payload, wait_note)
    return json_payload({**(running() if callable(running) else running), "elapsed_s": round(elapsed, 1),
                         **({"wait_seconds_note": wait_note} if wait_note else {})})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_load_status(load_id: str | None = None, wait_seconds: float = 0.0) -> str:
    """`load_dataset` が返した受付番号（`ld_…`）の読み込みを引く。引いても待っても状態は変わらない。

    load_id: 省くと、まだ結果を渡していない読み込みが 1 本ならそれ、2 本以上なら `needs_input`（受付の一覧つき）、
        無ければ `missing_state`（`load_dataset`）。
    wait_seconds: 終わるのを待つ秒数（既定 0 はすぐ返す。0〜30。範囲の外は丸めて `wait_seconds_note`。NaN は
        `DATASET_BAD_REQUEST`。真偽値と数に読める文字列は MCP の引数検査が数に変え、読めない値は封筒でない
        引数検査の誤り）。

    読み込み中は `{"status": "loading", "load_id", "elapsed_s"}`（受け付けてからの秒）。終われば `load_dataset`
    の戻り値そのもの: 成功ならそのハンドルが直近になっている（受付の後に別のハンドルが直近になっていたら直近に
    せず `current: false` と `note`）、失敗なら同じ誤りの封筒（ハンドルは登録されず直近も動かない）。終わった
    受付は 1 時間覚え、何度引いても同じ文字列を返す。知らない番号（1 時間を過ぎた・サーバが再起動した）は
    `DATASET_LOAD_NOT_FOUND`。
    """
    load_id, envelope = _pick_pending(
        pending_loads.registry, load_id, field="load_id", tool="dataset_load_status",
        start_tools=["load_dataset"], state="load")
    if envelope is not None:
        return envelope
    return _poll(pending_loads.registry, load_id, wait_seconds,
                 {"status": "loading", "load_id": load_id})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def job_status(job_id: str | None = None, wait_seconds: float = 0.0) -> str:
    """長い道具（`pathway_project` / `dataset_homolog_check`）が返した受付番号（`job_…`）を引く。

    その道具は自分の `wait_seconds`（既定 20 秒）で終わらなければ、仕事を続けたまま
    `{"status": "running", "job_id": "job_…", "tool", "dataset_handle", "next"}` を返す。仕事はそのハンドルの
    ロックの下で続き、同じハンドルへのほかの呼び出しは最長 25 秒待って何もせず `HANDLE_BUSY`（`details.job_id`
    がこの仕事）を返すので、ここで終わりを待ってから呼び直す。

    job_id: 省くと、まだ結果を渡していない受付が 1 本で、直近のハンドル（かハンドルなし）の仕事ならそれ。
        直近でないハンドルの仕事 1 本や 2 本以上は `needs_input`（label にハンドル）、無ければ `missing_state`。
        受け取った番号を渡すのが確実。
    wait_seconds: 終わるのを待つ秒数（既定 0。0〜30。範囲の外は丸めて `wait_seconds_note`）。

    走っている間は `{"status": "running", "job_id", "tool", "dataset_handle", "elapsed_s"}`（別の長い仕事の
    ロック待ちなら `queued_behind`）。終われば道具の戻り値（成功か誤りの封筒。JSON のオブジェクトなら先頭に
    `job_id` と `dataset_handle` を足す）。終わった受付は 1 時間覚え、何度引いても同じ戻り値を返す。知らない
    番号は `JOB_NOT_FOUND`。
    """
    current = session_state.session.current_handle_id

    def foreign(entry: dict) -> str | None:
        label = entry.get("label")
        if label is None or label == current:
            return None
        return (f"結果を渡していない受付は 1 本ですが、直近でないハンドル {label} の仕事です"
                "（別のエージェントのものかもしれない。受け取った job_id を渡すのが確実）。")

    job_id, envelope = _pick_pending(
        pending_loads.jobs, job_id, field="job_id", tool="job_status",
        start_tools=["pathway_project", "dataset_homolog_check"], state="job", foreign=foreign)
    if envelope is not None:
        return envelope
    tool = pending_loads.jobs.tool_of(job_id)
    handle_id = pending_loads.jobs.label_of(job_id)

    def running() -> dict:   # 待った後に組む（順番待ちは待つ間に変わりうる）
        return {"status": "running", "job_id": job_id, **({"tool": tool} if tool else {}),
                **({"dataset_handle": handle_id} if handle_id else {}),
                **handle_dispatch.queued_behind(handle_id, job_id)}

    return _poll(pending_loads.jobs, job_id, wait_seconds, running,
                 done_fields={"job_id": job_id, "dataset_handle": handle_id})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_list() -> str:
    """開いているハンドルの一覧（TSV: `dataset_handle current polarity assay_kind source_kind source
    n_features n_samples library_sha256`。古く開いた順）。`current` は直近のハンドル（`dataset_handle` を省いた
    呼び出しが使う。`load_dataset` が開いたものが直近）。`assay_kind` はハンドルごとの種別
    （`record_objective(assay_kind=...)` で確定）。"""

    from metabolomix.arf2.reader import format_spots_as_table

    current = session_state.session.current_handle_id
    rows = []
    for handle in handles.registry.list():
        store = handle.library.store
        rows.append({
            "dataset_handle": handle.handle_id,
            "current": handle.handle_id == current,
            "polarity": handle.polarity,
            "assay_kind": handle.assay_kind,
            "source_kind": handle.source.get("kind"),
            "source": _input_name(handle.source.get("path")) or "",
            "n_features": len(handle.dataset.feature_ids),
            "n_samples": len(handle.dataset.sample_names),
            "library_sha256": store.summary()["source_sha256"][:16] if store is not None else None,
        })
    columns = ["dataset_handle", "current", "polarity", "assay_kind", "source_kind", "source",
               "n_features",
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

    直近のハンドル（`dataset_list` の `current`）を閉じると、`dataset_handle` を省いた呼び出しは
    `missing_state` になる（開いている別のハンドルへは落ちない。`load_dataset` で開き直す）。
    そのハンドルで長い道具（`pathway_project` / `dataset_homolog_check`）が走っていれば、最長 25 秒待って
    何も閉じずに `HANDLE_BUSY`（`details.job_id` / `tool`）を返す（仕事は止めず、後で黙って閉じることもない）。
    `job_status(job_id)` で終わりを待ってから閉じ直す。

    戻り値: `closed`・`source`（入口の種類・入力の名前（`pipeline-run.json` は run のフォルダ名）・極性。パスは
    出さない）・`was_current`・`open_handles`（残りの数）。
    """
    handle = handles.registry.peek(dataset_handle)
    was_current = dataset_handle == session_state.session.current_handle_id
    try:
        handles.registry.close(dataset_handle)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    source = getattr(handle, "source", None) or {}
    closed_source = {key: value for key, value in (
        ("kind", source.get("kind")),
        ("name", _input_name(source.get("path"))),
        ("polarity", getattr(handle, "polarity", None))) if value}
    return json_payload({"status": "success", "closed": dataset_handle,
                         **({"source": closed_source} if closed_source else {}),
                         "was_current": was_current,
                         "open_handles": len(handles.registry.list())})
