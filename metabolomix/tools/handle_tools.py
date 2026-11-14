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
        run の受領（`results/receipt.json`）が指す、極性ごとの解析フォルダを開く。

    polarity: `pipeline_path` のときだけ使う。開く極性の単位（`"negative"` / `"positive"`）。
        単位が 1 つの run では省ける。2 極性の run で省くと `needs_input`（選べる極性つき。
        検証を通った単位だけ）。検証を通っていない単位は `PIPELINE_UNIT_NOT_VALIDATED`
        （`allow_incomplete=True` なら探索用に読む）。他の入口で渡すと
        `PIPELINE_POLARITY_NOT_APPLICABLE`。2 極性の run はハンドルが 2 つになる。

    source: どのアダプタで読むか。`"auto"`（既定）/ `"folder"` / `"mztab"`。
        - `directory` は解析フォルダ（`"mztab"` はエラー）、`mztab_path` は mzTab-M
          （`"folder"` はエラー）。
        - `job_path` / `pipeline_path` は既定（`"auto"` / `"folder"`）で、そのジョブの正準
          mzTab-M と同じアラインメントの**解析フォルダ**を読む。候補（mzTab-M の隣の `.arf`。
          pipeline はどれも合わなければ単位の解析フォルダ `<極性>/input/*.arf2` のバッチを
          新しい順に 1 本ずつ）を mzTab-M と特徴数・m/z・試料名で照合し、合うものがちょうど
          1 本ならそれを読む。無ければ `ALIGNMENT_NOT_FOUND`、2 本以上なら
          `ALIGNMENT_AMBIGUOUS`（どちらも候補のファイル名を `details.candidates` に載せる）で
          止まり、mzTab-M へ黙って切り替えない。`"mztab"` を渡すと mzTab-M を読む。
        - **pipeline の run から引き継ぐのは試料の対応表だけ**（v3 の run は受領
          `results/receipt.json` の `units.<極性>.samples.rows`、旧 v1 の run は記録の
          `inputs.manifest`。どちらのアダプタでも同じ）。run が計算した差次的結果は戻さない
          （比較は読み込んだあと `dataset_compare` で回す）。受領が無い・照合できない v3 の
          run は記録の単位から解析フォルダを開き、対応表は載らない旨を注記する。
    measure: 解析フォルダのアダプタが値に使う `.arf` の欄（`"peak_height"`（既定）/
        `"peak_area_above_zero"`）。**この引数が効くのは解析フォルダのアダプタだけ**。
        mzTab-M のアダプタでは定量種別がファイル（job / pipeline ならジョブの宣言）で
        決まるので、既定値以外を渡すと `DATASET_BAD_REQUEST`。
    allow_incomplete: **完了していない実行**（partial / failed / running）の出力を
        読み込む（既定 False。bool 以外は `DATASET_BAD_REQUEST`——`"false"` という
        文字列を真と読まないため）。中断時点の生成物なので、既定では拒否する
        （`INCOMPLETE_ANALYSIS_JOB`）。True で読んだデータセットは探索専用になる（検証を
        通っていない pipeline の単位も同じ。`directory` で run の中の解析フォルダを開いたときも、その
        単位が検証を通っていなければ同じ印を付ける。断らない）。比較も書き出しもできるが、欠けた検体を
        「その群には無い」と読み違えうるので、差次的 TSV のメタ行（`# exploratory = true`）・
        根拠バンドルの manifest（`exploratory`）・図・最初の概観（dataset_overview）の比較を使う段の事実に
        印と理由が載る。
    library: `"auto"`（既定）はデータセットが同定に使ったライブラリの記憶
        （`library_record`）を読み取り専用で開く。実体の sha256 を測り直し、変わっている・
        無いとき（`mismatch`）、記憶が決められないとき（`unknown`）、LBM だけのとき
        （`lbm_only`）は開かず、`library_note` に理由。**極性の設定
        （`[library] msp_positive` / `msp_negative`）には落ちない**。別のライブラリは
        `library_load` の `file_path` で明示する。`"none"` は読まない。
    arf2_path: スポットの照合（spot_id = SMF_ID）に使う `.arf2`。省略時、解析フォルダで
        読んだならその `.arf2` と `.dcl`。mzTab-M で読んだなら検出状態を取り込んだ `.arf` の
        同じバッチ、無ければ mzTab-M の隣、`pipeline_path` なら受領が照合した単位の
        `<極性>/input/` の `.arf2`、受領に無ければ単位の `<極性>/input/*.arf2`（複数あれば
        名前の新しいもの）を探す。

    overview: True なら読み込みの最後に最初の概観（`dataset_overview` と同じ。品質 → 全体 → クラス →
        分子種）を計算して戻り値の `overview` 欄に載せる（既定 False。読み込みの戻り値を軽く保ち、
        概観は `dataset_overview()` で出す）。bool 以外は `DATASET_BAD_REQUEST`。概観が失敗しても
        読み込みは成功し、`overview` に `status:"unavailable"` と理由が載る。

    wait_seconds: 読み込みを待つ秒数（既定 20。0〜30 秒。範囲の外の数は丸めて `wait_seconds_note` を
        付ける。NaN は `DATASET_BAD_REQUEST`。数でない値は MCP の引数検査が先に扱う: 真偽値と数に読める
        文字列は数に変えて通り（true は 1 秒）、数に読めない値は封筒でない引数検査のエラー）。
        この中で終われば下の戻り値をそのまま返す。終わらなければ読み込みは続け、受付番号を
        `{"status": "loading", "load_id": "ld_…", "next": "dataset_load_status(load_id=…)"}` で返す。
        `dataset_load_status()` を呼ぶと、読み込み中は `status: "loading"`、終われば
        この関数の戻り値そのもの（成功か誤りの封筒）が返る（その `wait_seconds` を渡すと、最長
        30 秒まで終わりを待ってから返す）。ハンドルの登録と直近の切り替えは
        読み込みが終わったときに 1 回だけ起き、受付の後に失敗したらハンドルは登録しない。
        受付番号を返した後に別のハンドルが直近になっていたら、終わった読み込みはハンドルを
        登録するだけで直近にせず（ライブラリも移さない）、戻り値に `current: false` と `note` が
        付く（そのハンドルは `dataset_handle` を渡して使う）。0 は「既に終わっていなければすぐ
        受付番号を返す」。引数の誤りと上限（`DATASET_HANDLE_LIMIT`。読み込み中の受付も 1 本と
        数える）は受け付ける前に返す。

    アッセイ種別（`assay_kind`）はハンドルごとに持つ。解析フォルダ配下 `reports/` の目的ファイル
    （`record_objective` が書いたもの）に `lipid` / `metabolite` が記録されていれば、それを読み戻して
    このハンドルの種別にし（`notes` に注記）、無ければ `unknown`（`record_objective(assay_kind=...)` で確定）。
    データセットを開く前に（`[paths] reports_dir` の目的へ）記録した種別は、そのあと最初に直近になる
    ハンドルへ引き継ぐ（目的ファイルの種別と食い違えば unknown。どちらも `notes` に注記）。

    成功すると、いつも**新しいハンドルを直近にする**。`dataset_handle` を省いた `dataset_*` /
    `library_*` はこの直近のハンドルに効く。再解析のエージェントは返った `dataset_handle` を
    各ツールに渡して使い、直近を当てにしない（渡した呼び出しはサーバの別スレッドで回り、同じ
    ハンドルへの呼び出しは 1 本ずつ、別のハンドルどうしは並行に進む）。前のハンドルは閉じずに
    残る（`dataset_list` / `dataset_close`）。開いているハンドルが上限（`lipidmix.local.toml` の
    `[server] max_dataset_handles`、既定 8）なら `DATASET_HANDLE_LIMIT` で止め、古いハンドルを
    黙って閉じない（直近も動かさない）。同じ入力（読み込みが読んだ・探したファイルの size・
    更新時刻が同じ）を開き直すと読み込み済みのものを使い回す（`reused_base`）。それでも
    ハンドルどうしの状態は別々。

    戻り値: `dataset_handle`・`assay_kind`（このハンドルの種別）・読んだアダプタ・選んだ `.arf`・
    ライブラリの記録（ファイル名と sha256。置き場所は出さない）に、`caveat`（意味論の要約。未読で、その
    種別にまだ出していないときだけ先頭に載る。`assay_kind` が確定したあとの最初の 1 回は
    `record_objective` の戻り値が運ぶ）、`summary`（特徴・試料の数、検出状態、同定の内訳、読み込みの
    注記）、`sample_counts`（群と役割ごとの試料数。群は統計に入る試料だけを数え、実験情報シートの include=false
    の試料は `excluded_by_sheet` に分ける）、`overview`（最初の概観。`overview=True` のとき）、
    `next`（次の一手の例）。概観の比較と PCA は `result_id` で掘れるが直近の結果にはしない。
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
    """`load_dataset` が返した受付番号（`ld_…`）の読み込みを引く（最長 `wait_seconds` 秒待つ）。

    load_id: 引く受付番号。省くと、まだ結果を渡していない読み込みから選ぶ（ちょうど 1 本ならそれ、
        2 本以上なら `needs_input`（選べる受付つき）、無ければ `missing_state`（`load_dataset`））。
    wait_seconds: 読み込みが終わるのを待つ秒数（既定 0 は待たずにすぐ返す。0〜30 秒。範囲の外の数は
        丸めて `wait_seconds_note` を付ける。NaN は `DATASET_BAD_REQUEST`。数でない値は MCP の引数検査が先に
        扱う: 真偽値と数に読める文字列は数に変えて通り、数に読めない値は封筒でない引数検査のエラー）。
        その中で終われば、すぐ下の戻り値を返す。呼び直しの回数を減らす。

    読み込み中は `{"status": "loading", "load_id", "elapsed_s"}`（`elapsed_s` は受け付けてからの秒）。
    終われば `load_dataset` の戻り値そのもの: 成功なら `dataset_handle` つきの要約（そのハンドルが
    直近になっている。受け付けた後に別のハンドルが直近になっていたら直近にせず、`current: false` と
    `note` が付く）、失敗なら `load_dataset` と同じ誤りの封筒（ハンドルは登録されず、直近も
    動かない）。引いても待っても状態は変わらない——ハンドルの登録と直近の切り替えは、読み込みが
    終わったときに 1 回だけ起きている。終わった受付は 1 時間覚えていて、番号を渡せば何度引いても同じ
    文字列を返す。知らない受付番号（終わってから 1 時間を過ぎた・サーバが再起動した）は
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

    その道具は自分の `wait_seconds`（既定 20 秒）の中で終わらなければ、仕事を続けたまま
    `{"status": "running", "job_id": "job_…", "tool", "dataset_handle", "next"}` を返す（デスクトップアプリの
    中継が 36〜60 秒ほどで呼び出しを打ち切るため）。仕事はそのハンドルのロックの下で続き、その間、同じ
    ハンドルへのほかの呼び出しは最長 25 秒待ってから何もせずに `HANDLE_BUSY`（`details.job_id` がこの仕事）を
    返す——この関数で終わりを待ってから呼び直す。

    job_id: 引く受付番号（受け取った番号を渡すのが確実）。省くと、受付番号を渡し済みで、まだ結果を渡して
        いない受付から選ぶ: ちょうど 1 本で、それが直近のハンドル（または ハンドルなし）の仕事ならそれ。
        直近でないハンドルの仕事 1 本（別のエージェントのものかもしれない）や 2 本以上なら `needs_input`
        （選べる受付つき。label にハンドル）、無ければ `missing_state`。
    wait_seconds: 終わるのを待つ秒数（既定 0。0〜30 秒。範囲の外の数は丸めて `wait_seconds_note`）。

    走っている間は `{"status": "running", "job_id", "tool", "dataset_handle", "elapsed_s"}`（別の長い仕事が
    ロックを持っていて順番を待っていれば `queued_behind`）。終われば、その道具の戻り値（成功か
    誤りの封筒。JSON のオブジェクトなら、その先頭に `job_id` と `dataset_handle`（その仕事のハンドル）を足したもの）。終わった受付は
    1 時間覚えていて、番号を渡せば何度引いても同じ戻り値を返す。知らない受付番号は `JOB_NOT_FOUND`。
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
    n_features n_samples library_sha256`）。古く開いた順。`current` は対話のセッションの直近の
    ハンドル（`dataset_handle` を省いた呼び出しが使うもの。`load_dataset` が開いたものが直近）。
    `assay_kind` はそのハンドルの種別（ハンドルごと。`record_objective(assay_kind=...)` で確定）。"""

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

    対話のセッションの直近のハンドル（`dataset_list` の `current`）を閉じると、
    `dataset_handle` を省いた呼び出しはデータセットが無い状態（`missing_state`）になる。
    開いている別のハンドルへは落ちない（もう一度 `load_dataset` で読む）。

    閉じるのもそのハンドルのロックの下で行う。長い道具（`pathway_project` / `dataset_homolog_check`）が
    そのハンドルで走っている間は、最長 25 秒待ってから**何も閉じずに** `HANDLE_BUSY`（`details.job_id` と
    `tool` がその仕事）を返す——仕事を止めはしない。`job_status(job_id)` で仕事の終わり（か結果）を待ってから
    閉じ直す。待ちきれなかった閉じ方が、仕事の後で黙って効くことはない。

    戻り値: `closed`（閉じたハンドル）・`source`（何を閉じたか: 入口の種類・入力の名前（`pipeline-run.json`
    はその run のフォルダの名前）・極性。パスは出さない。開き直すときの手がかり）・`was_current`（直近だったか）・
    `open_handles`（残りの数）。
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
