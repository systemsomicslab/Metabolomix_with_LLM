"""参照ライブラリ（MS/MS スペクトル照合）の MCP ツール群。

下の層（`metabolomix.library.dbs` / `.msp` / `.store`、`metabolomix.analysis.spectral_match`、
`metabolomix.plots.mirror`）はすべて完成済み。ここはそれらを MCP の公開面へ繋ぐだけの層で、
`dcl/tools.py` と同じ流儀（`missing_state` 封筒・`json_payload`・TSV 一覧・
`structured_output=False`）に揃える。

deps: mcp_core / session_state / mcp_errors / path_resolvers / library.{dbs,msp,store} /
analysis.spectral_match / plots.{mirror,render} / dcl.reader / arf2.reader
（`format_spots_as_table` の借用のみ）。server は import しない。
"""
from __future__ import annotations

from pathlib import Path

from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.analysis.spectral_match import fragment_matches_excluding_precursor, match_spectrum
from metabolomix.arf2.reader import format_spots_as_table
from metabolomix.core import mcp_errors, user_config
from metabolomix.core.mcp_core import mcp
from metabolomix.core.path_resolvers import (
    LIBRARY_SETTING_KEYS,
    LibraryPathError,
    resolve_arf2_file_path,
    resolve_dcl_file_path,
    resolve_library_path,
)
from metabolomix.msdial import alignment_spot
from metabolomix.core.serialization import json_payload, round_floats
from metabolomix.dcl.reader import deserialize_dcl, get_msms_by_precursor
from metabolomix.library.defaults import DEFAULT_MS2_TOL as _DEFAULT_MS2_TOL
from metabolomix.library.defaults import DEFAULT_MZ_TOL as _DEFAULT_MZ_TOL
from metabolomix.library.defaults import DEFAULT_RT_TOL as _DEFAULT_RT_TOL
from metabolomix.library import matching
from metabolomix.library.store import open_store
from metabolomix.plots import mirror as mirror_plot
from metabolomix.plots import render as plot_render

__all__ = ["library_load", "library_match_feature", "library_plot_mirror", "library_search"]

# library_load が要約に添える化合物クラス分布の上位件数。
_TOP_COMPOUND_CLASSES = 10

# library_match_feature の候補の並び順。MS-DIAL の総合スコア
# （`spectral_match.total_score` = 上流 `GetTotalScore`）を採用している。
# 以前は weighted_dot_product 単独だったが、同名・同 precursor の別レコードで
# 実用上の最良候補が 1 位に来ないことがあった（実データ 120 feature の
# top-1 一致 88.3% → 94.2%。HISTRY 2026-09-22(1)）。戻り値の "ranked_by" と
# docstring の両方がこの定数を指す——コード内コメントだけにしておくと、
# rank 列が何の降順かを呼び出し側（MCP クライアント）が知る手段が無くなる。
_RANK_KEY = "total_score"

# 候補一覧 TSV の列（先頭に列名を 1 回だけ出す）。スペクトル座標・alignment は
# 意図的に含めない（座標は session.library.last_match にだけ持つ）。
_CANDIDATE_TABLE_COLUMNS = [
    "rank", "name", "precursor_mz", "ion_mode", "adduct", "rt", "formula",
    "ontology", "compound_class", "total_score",
    "simple_dot_product", "weighted_dot_product",
    "reverse_dot_product", "matched_peaks_percentage", "matched_peaks_count",
    "entropy_similarity",
]
#: `detail=True` のときだけ末尾に足す列（既存の列の並びは変えない）。
_DETAIL_COLUMNS = ["inchikey", "fragment_matches_excl_precursor"]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def library_load(file_path: str | None = None, rebuild: bool = False,
                 ion_mode: str | None = None, dataset_handle: str | None = None) -> str:
    """参照ライブラリを読み込み、照合用の SQLite store を構築（または既存キャッシュを
    再利用）する。

    解決順: `file_path` の明示 → 直近のデータセット（または `dataset_handle`）が覚えた、同定に
    使ったライブラリ（`library_record`。実体の sha256 を測り直して通れば開く）→ `ion_mode`（`"positive"` / `"negative"`）に対応する
    設定（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の
    `[library] msp_positive` / `msp_negative`）→ データディレクトリの
    `*_Loaded.msp2.dbs` → 設定済みの極性（両方あれば `ion_mode` を求める）→
    データディレクトリの `*.msp`。候補が 1 つに決まらなければ `code` 付きの
    エラー（`MSP_AMBIGUOUS` など）を返し、黙って選ばない。一度に保持する
    ライブラリは 1 つだけ——測定の極性に合わせて `ion_mode` を指定する。

    記憶が `mismatch`（内容が変わった・ファイルが無い）や `lbm_only`（LBM だけで MSP を使って
    いない）のときは、上の設定の順には**落ちず**理由を返す（`file_path` で明示すれば開ける）。
    `file_path` が記憶と違うライブラリなら、そのまま開き、戻り値の `differs_from_record` に
    記憶したファイル名と開いたファイル名を書く。記憶から開いたときは `library_source` が
    `"dataset_record"`。

    初回の構築は大きな `.msp`（1 GB 級）で数分かかる。MCP のタイムアウトに当たる
    なら、手元で `python -m metabolomix.library.store --ion-mode <極性>` を 1 度
    走らせておけば、以後はキャッシュを開くだけで済む。

    `library_match_feature` の前提。別ライブラリへ切り替えると直近の照合結果
    （`library_plot_mirror` が読む座標）は破棄する——古い照合を新ライブラリの
    結果と取り違えないため。

    `rebuild=True` はキャッシュを無視して構築し直す（元ファイルが壊れている
    疑いがあるときなど）。通常は不要——store は元ファイルの sha256 をキーに
    キャッシュされるため、内容が変わらない限り再構築しない。

    dataset_handle: `dataset_handle_open` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッション。ハンドルのライブラリ・照合結果だけを読み書きする。
        ハンドルへは読み取り専用で開く。
    """
    from metabolomix.tools._holder import resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    from_record = False
    differs = None
    record = getattr(holder.dataset, "library_record", None)
    if record is not None and record.get("status") != "unknown":
        from metabolomix.dataset.library_memory import verify_record

        record = verify_record(record)
        status, remembered = record["status"], record.get("file")
        if file_path is None:
            if status == "known":
                file_path, from_record = record["path"], True
            elif status == "mismatch":
                return json_payload({
                    "status": "error", "code": "LIBRARY_RECORD_MISMATCH",
                    "message": (f"データセットが覚えたライブラリ（{remembered}）が不一致です"
                                "（内容が変わったか、ファイルが無い）。設定の順には落ちません。"
                                "使うライブラリを file_path で明示してください。")})
            else:  # lbm_only
                return json_payload({
                    "status": "error", "code": "LIBRARY_RECORD_LBM_ONLY",
                    "message": ("このデータセットの同定は LBM（脂質のルールベース）だけで、MS/MS の"
                                "照合に使う MSP のライブラリがありません。設定の順には落ちません。"
                                "別のライブラリで照合するなら file_path で明示してください。")})
        elif status in ("known", "mismatch") and record.get("path"):
            if Path(file_path).expanduser().resolve() != Path(record["path"]).resolve():
                differs = {"remembered": remembered, "opened": Path(file_path).name}
    if dataset_handle is not None and file_path is None:
        # ハンドルは極性の設定だけから引く（dataset_handle_open と同じ規則）。対話の
        # データフォルダの .dbs / .msp を探さず、データセットと逆の極性も読まない。
        mode = (ion_mode or holder.polarity or "").strip().lower() or None
        if ion_mode and holder.polarity and mode != holder.polarity:
            return json_payload({
                "status": "error", "code": "LIBRARY_POLARITY_MISMATCH",
                "message": (f"このハンドルのデータセットは {holder.polarity} です"
                            f"（ion_mode={ion_mode!r}）。逆の極性のライブラリは読みません。")})
        key = user_config.MSP_SETTING_KEYS.get(mode or "")
        setting = user_config.get_setting(key) if key else None
        if setting is None:
            reason = ("データセットの極性が決まらないので ion_mode を指定してください。"
                      if key is None else f"{key} が設定されていません。")
            return json_payload({
                "status": "error", "code": "LIBRARY_NOT_CONFIGURED",
                "message": ("ハンドルのライブラリは極性の設定（[library] msp_positive / "
                            "msp_negative）から引きます。" + reason)})
        file_path = setting.value
    try:
        resolved = resolve_library_path(file_path, ion_mode=ion_mode)
    except LibraryPathError as exc:
        payload = {"status": "error", "code": exc.code, "message": exc.message}
        if exc.details:
            payload["details"] = exc.details
        return json_payload(payload)
    if not resolved:
        # resolve_library_path は大文字・空白を許して正規化するので、ここでも同じに揃える。
        mode = (ion_mode or "positive").strip().lower()
        key = LIBRARY_SETTING_KEYS.get(mode, "library.msp_positive")
        hint = user_config.missing_hint(key, "研究室の参照ライブラリ（.msp）")
        if ion_mode is None:
            hint += " 負イオンは [library] msp_negative です。"
        return json_payload({
            "status": "error",
            "message": "参照ライブラリが見つかりませんでした（データディレクトリにも "
                       "*_Loaded.msp2.dbs / *.msp がありません）。file_path で指定するか、" + hint,
            "details": user_config.describe_missing(key),
        })

    try:
        store_obj = open_store(resolved, rebuild=rebuild, read_only=dataset_handle is not None)
    except Exception as exc:  # noqa: BLE001 - 壊れたライブラリは文言で返す（MCP が扱いやすい）
        # OSError の文言は置き場所（フルパス）を含む。研究室ライブラリの置き場所は
        # 戻り値（LLM の文脈）に出さないので、ディレクトリ部分を伏せる。
        # OS が返す文言は `\` 区切り・`/` 区切りのどちらでも出うる（環境変数に `/` で
        # 書かれた値はそのまま `/` で残る）ので、両方の表記を伏せる。
        parent = Path(resolved).parent
        detail = f"{type(exc).__name__}: {exc}"
        for shown in (str(parent), parent.as_posix()):
            detail = detail.replace(shown, "…")
        return json_payload({"status": "error", "message": f"参照ライブラリの読み込みに失敗しました: {detail}"})

    # 古い store の sqlite3 接続を閉じてから差し替える（Minor 8: 閉じずに上書きすると
    # library_load を繰り返すたびに接続が漏れる。LibraryStore.close() は最初から
    # 存在していたが、ここで呼んでいなかった）。
    old_store = holder.library.store
    if old_store is not None:
        old_store.close()

    # 別ライブラリへの切り替え起こりうるので、古い照合結果を新ライブラリのものと
    # 取り違えないよう破棄する（ArfState.reset_analysis と同じ考え方）。
    holder.library.store = store_obj
    # ハンドルにはファイル名だけを持たせる（研究室のライブラリの置き場所を記録に残さない）。
    holder.library.source_path = resolved if dataset_handle is None else Path(resolved).name
    holder.library.last_match = None

    path = Path(resolved)
    summary = store_obj.summary()
    search_params = summary.get("search_params")
    skipped = summary.get("skipped_no_precursor_mz", 0)

    payload = {
        "status": "success",
        "file": path.name,
        "source_sha256": summary["source_sha256"],
        "record_count": summary["record_count"],
        "ion_modes": summary["ion_modes"],
        "compound_classes": store_obj.compound_class_counts(_TOP_COMPOUND_CLASSES),
        "search_params": search_params,
        "skipped_no_precursor_mz": skipped,
        "records_without_ion_mode": summary["records_without_ion_mode"],
        "non_utf8_lines": summary["non_utf8_lines"],
    }
    notes = []
    if from_record:
        payload["library_source"] = "dataset_record"
    if differs is not None:
        payload["differs_from_record"] = differs
        notes.append(
            f"データセットが同定に使ったライブラリ（{differs['remembered']}）とは違う"
            f"（{differs['opened']}）を開きました。")
    if search_params is None:
        notes.append(
            f"`.msp` には照合の許容幅（search_params）が同梱されないため、"
            f"library_match_feature は既定値 mz_tol={_DEFAULT_MZ_TOL} / "
            f"ms2_tol={_DEFAULT_MS2_TOL} を使用します（明示的に指定すれば上書きできます）。"
        )
    if skipped:
        # Important 3: PRECURSORMZ が無い/パース不能なレコードは黙って捨てず、
        # 件数を表に出す（黙って捨てると「record_count が元ファイルと合わない」で悩む）。
        notes.append(
            f"precursor m/z が無い（または解釈できなかった）レコードを {skipped} 件、"
            f"読み飛ばしました（record_count には含まれません）。"
        )
    if summary["records_without_ion_mode"]:
        notes.append(
            f"極性（IONMODE）の無いレコードが {summary['records_without_ion_mode']} 件あります。"
            f"library_match_feature の ion_mode 絞り込みはこれらを極性不明として候補に残します"
            f"——極性ごとに分かれたライブラリなら、測定の極性に合ったファイルを読んでいるか確認してください。"
        )
    if summary["non_utf8_lines"]:
        notes.append(
            f"UTF-8 で読めない行が {summary['non_utf8_lines']} 行あり、cp932（無理なら latin-1）で"
            f"読みました。化合物名が化けている可能性があります。"
        )
    if notes:
        payload["note"] = " ".join(notes)
    return json_payload(round_floats(payload))


def _measured_spectrum(precursor_mz: float, *, rt: float | None = None,
                       dcl_file: str | None = None, mz_tol: float = _DEFAULT_MZ_TOL,
                       rt_tol: float = _DEFAULT_RT_TOL) -> list[list[float]] | None:
    """`.dcl` から測定 MS/MS を引く。**全ピーク**（`top_n_peaks` で間引かない——
    間引くと採点の数値が変わる。Task 1 で確認済み）。

    `.dcl` が無い、または該当 precursor の MS/MS が無ければ `None`。

    **`mz_tol` / `rt_tol` の既定はここだけで完結させる**（呼び出し元の
    `library_match_feature` は library 側の許容幅——`.dbs` の `search_params`
    由来なら実値、RT は上流既定 100.0 相当——をここへ絶対に流し込まないこと。
    最終レビュー Important 7: 以前は library 側で解決した `resolved_rt_tol`
    をここにも渡していたため、その `.dbs` を読ませると `.dcl` 側の RT 窓が
    ±100 分になり、precursor が近い無関係なピークまで拾っていた）。

    `.dcl` の同じ precursor に複数ヒットしうる（同一 m/z の別溶出ピーク）ため、
    `rt` が渡されたときは **RT 距離が最も近いヒットを選ぶ**（最終レビュー
    Important 7: 以前は `hits[0]`——`.dcl` ファイル内の並び順の先頭——を無条件に
    採っており、別のピークの測定スペクトルを黙って採点しうる不具合があった）。
    """
    resolved = resolve_dcl_file_path(dcl_file)
    if not resolved:
        return None
    try:
        results = deserialize_dcl(resolved, include_spectrum=True, top_n_peaks=None)
    except Exception:  # noqa: BLE001 - 壊れた .dcl は「見つからない」と同じ扱いにする
        return None
    hits = get_msms_by_precursor(results, precursor_mz, tol=mz_tol, rt=rt, rt_tol=rt_tol)
    if not hits:
        return None
    if rt is not None:
        hits = sorted(hits, key=lambda hit: abs(hit["rt"] - rt))
    return hits[0]["msms_spectrum"]


_NO_MSMS_MESSAGE = (
    "MS/MS 未取得（.dcl が無い、または該当の記録が無い）のであって、"
    "期待フラグメントが合わなかったのではありません。"
)


def _spot_query(spot_id: int, *, arf2_path: str | None, dcl_file: str | None,
                rt: float | None, ion_mode: str | None) -> dict:
    """アラインメントのスポットから照合の入力を組む。失敗は `status` 付きの封筒。

    precursor m/z・RT・極性は ARF2、測定スペクトルはアラインメント `.dcl` の
    `dcl[spot_id]`（precursor を m/z で探さない——同じ m/z の別スポットを拾わないため）。
    `not_found` には `reason` が付く: `no_arf2` / `no_dcl`（ファイルが無い＝スペクトルの有無は
    分からない）と `no_spectrum`（照会して無かった）。`status` の値は変えない。
    """
    arf2 = resolve_arf2_file_path(arf2_path)
    if not arf2:
        return {"status": "not_found", "reason": "no_arf2",
                "message": "spot_id の解決に要る .arf2 が見つかりませんでした（arf2_path で指定できます）。"}
    record = alignment_spot.spot_record(arf2, spot_id)
    if record is None or record.get("MassCenter") is None:
        return {"status": "spot_not_found",
                "message": f"spot_id={spot_id} は {Path(arf2).name} にありません。"}
    precursor_mz = float(record["MassCenter"])
    dcl = Path(dcl_file) if dcl_file else alignment_spot.alignment_sibling(arf2, ".dcl")
    if dcl is None or not dcl.is_file():
        return {"status": "not_found", "reason": "no_dcl",
                "message": f"spot_id={spot_id} のアラインメント .dcl が見つかりません。" + _NO_MSMS_MESSAGE}
    msms = alignment_spot.alignment_msms(dcl, spot_id, expected_mz=precursor_mz,
                                         expected_rt=record.get("RT"))
    if msms["status"] == "precursor_mismatch":
        return {"status": "dcl_precursor_mismatch",
                "message": (
                    f"{dcl.name} の {spot_id} 番（precursor m/z {msms['precursor_mz']:.4f}・"
                    f"RT {msms['rt']:.2f}）が ARF2（{precursor_mz:.4f}・{float(record.get('RT') or 0):.2f}）"
                    f"と ±{alignment_spot.PRECURSOR_TOLERANCE} Da / ±{alignment_spot.RT_TOLERANCE} 分"
                    "で合いません（索引対応の崩れ。照合しません）。"
                )}
    if msms["status"] != "ok":
        return {"status": "not_found", "reason": "no_spectrum",
                "message": f"spot_id={spot_id} に一致する MS/MS が見つかりませんでした。" + _NO_MSMS_MESSAGE}
    polarity = ion_mode or str(record.get("IonMode") or "").lower() or None
    return {
        "status": "ok",
        "precursor_mz": precursor_mz,
        "rt": rt if rt is not None else record.get("RT"),
        "ion_mode": polarity,
        "measured": msms["spectrum"],
    }


#: スポットの代表 MS/MS を引く口（`dcl_spot_spectrum`・再解析の証拠の束が共有する）。
spot_msms_query = _spot_query


def lbm_only_not_evaluable(holder) -> dict | None:
    """ライブラリが開かれておらず、データセットの同定が LBM だけだったときの `not_evaluable`。

    MSP のスペクトルが無いので、照合は「合わなかった」のではなく評価できない。
    `missing_state`（先に library_load を）とは違い、library_load でも開けない状態を指す。
    """
    record = getattr(holder.dataset, "library_record", None)
    if not record or record.get("status") != "lbm_only":
        return None
    return {"status": "not_evaluable",
            "reason": ("このデータセットの同定は LBM（脂質のルールベース）だけで、MS/MS の照合に"
                       "使う MSP のライブラリがありません。評価できない（照合が合わなかった"
                       "のではない）。別のライブラリで照合するなら library_load の file_path で"
                       "明示する。")}


def identified_by_lbm(holder, spot_id) -> dict:
    """`{"identified_by_lbm": True}`、または `{}`（False は書かない）。

    データセットの同定が MSP と LBM の両方を使い（`lbm_also`）、照合した特徴の代表の注釈器が
    LBM 由来と分かるときだけ。照合結果が LBM の同定と食い違って見えるのを読み手が取り違えない
    ための印で、注釈器 ID そのもの（ライブラリの置き場所を含みうる）は返さない。
    """
    from metabolomix.dataset.library_memory import is_lbm_annotator

    ds = holder.dataset
    record = getattr(ds, "library_record", None)
    if spot_id is None or not record or not record.get("lbm_also"):
        return {}
    meta = (getattr(ds, "feature_metadata", None) or {}).get(str(spot_id)) or {}
    return {"identified_by_lbm": True} if is_lbm_annotator(meta.get("annotator_id")) else {}


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def library_match_feature(
    precursor_mz: float | None = None,
    rt: float | None = None,
    ion_mode: str | None = None,
    dcl_file: str | None = None,
    mz_tol: float | None = None,
    ms2_tol: float | None = None,
    rt_tol: float | None = None,
    top_n: int = 5,
    spot_id: int | None = None,
    arf2_path: str | None = None,
    detail: bool = False,
    precursor_exclusion_tol: float | None = None,
    dataset_handle: str | None = None,
) -> str:
    """測定 MS/MS（`.dcl`）を参照ライブラリ（`library_load` 済みの store）と照合し、
    上位候補を返す。

    照合の対象は `precursor_mz` か `spot_id` の**どちらか一方**で指定する。
    `spot_id`（MasterAlignmentID。`eic_search` 等で引ける）を渡すと、precursor m/z・
    RT・極性を同じバッチの `.arf2`（`arf2_path`、省略時は現在のデータフォルダ）から取り、
    測定スペクトルはアラインメント `.dcl`（`dcl_file`、省略時は `.arf2` と同じバッチ）の
    そのスポットの代表 MS/MS を使う。precursor が ARF2 と ±0.01 で合わなければ
    照合せず `status="dcl_precursor_mismatch"`、ID が無ければ `spot_not_found`。
    `rt` / `ion_mode` を明示すれば ARF2 の値より優先する。

    `ion_mode` は `"positive"` / `"negative"`（大小は問わない——`store.candidates()`
    が `COLLATE NOCASE` で比較するため、PAI2 由来の `IonMode.Positive.name`
    のような `"Positive"` 表記もそのまま渡してよい）。

    `mz_tol` / `ms2_tol` / `rt_tol` は明示指定が無ければ store の `search_params`
    （`.dbs` 由来なら実値）を使い、それも無ければ既定値
    （`_DEFAULT_MZ_TOL` / `_DEFAULT_MS2_TOL` / `_DEFAULT_RT_TOL`、現在の値は
    0.01 / 0.025 / 0.2）を使う。**これらは参照ライブラリの候補検索（`store.candidates()`
    の precursor m/z 窓・RT 窓）にだけ効く。** `.dcl` から測定 MS/MS を引く窓は
    `_measured_spectrum` 内部の固定既定（`dcl/reader.py` と同じ `mz_tol=0.01`/
    `rt_tol=0.2`）を常に使い、ここの引数や store の `search_params` の影響を受けない
    （最終レビュー Important 7: `.dbs` の `RtTolerance` 既定 100.0 のような
    ライブラリ用の値を `.dcl` 側の RT 窓に流用すると、同じ precursor m/z の
    別ピークの測定スペクトルを黙って拾ってしまうため、2 つの窓を分離した）。
    足切り（`relative_amp_cutoff` / `absolute_amp_cutoff`）と質量範囲
    （`mass_range_begin` / `mass_range_end`）も同じ `search_params` から採って
    採点（`match_spectrum`）に渡す——検証 CLI（`scripts/verify_spectral_match.py`）
    と同じ集合（最終レビュー Important 2）。

    `.dcl` に MS/MS が無い（未取得）場合は `status="not_found"` を返す——
    「候補と合わなかった」のではなく「照合する測定スペクトルがそもそも無い」ことを
    区別するため、`dcl_find_msms` と同じ文言の方針に揃える。

    候補は **`total_score` の降順**で返す（`rank` 列・`ranked_by` フィールドが
    同じ基準を指す）。MS-DIAL の総合スコア（上流 `GetTotalScore`）で、
    `rt + precursor_mz + (weighted+simple+reverse)/3 + matched_peaks_percentage`
    の**正規化しない和**——1 を超える。組み立て方は戻り値の `scoring` が明示する
    （RT 項が入ったかは `.dbs` の `IsUseTimeForAnnotationScoring` 次第で、
    数値だけを見ても分からないため）。同点は候補順を温存する（安定ソート）。

    **上流の順位付けと完全に同じではない。** `MsScanMatchResultContainer.ResultOrder`
    は `(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority,
    TotalScore)` の辞書順で、`TotalScore` は最後のタイブレークにすぎない。
    支配項の `IsReferenceMatched` は脂質クラス固有の判定
    （`Lipidomics/`、上流 66,932 行）に依存するため移植していない。実データでは
    この判定が無いと全候補が `True` になってゲートとして働かないので、
    **`total_score` 単独比較が現実的な最良の近似**という位置づけ。
    帰結として、上位に化学的にありえない候補が残ることがある
    （`docs/output_format/library.md` §14.9）。

    `top_n` は上位何件を返すかの上限（既定 5）。負値を渡すと全件を返す
    （明示的な仕様ではなく Python のスライス挙動に由来する副作用的な動作
    ——件数を絞りたくないだけなら大きな正の値を渡すこと）。

    候補一覧は TSV（列名 1 回）で返す。**スペクトル座標と alignment は戻り値に
    含めない**（LLM の文脈を食うため）。座標は `session.library.last_match` に
    持ち、`library_plot_mirror` がそこから読む。

    `detail=True`（opt-in。再解析で使う）のときだけ、TSV の末尾に `inchikey` と
    `fragment_matches_excl_precursor`（前駆体 ± `precursor_exclusion_tol` を除いた一致の
    本数。省略時の幅は採点に使った `ms2_tol`）の列を足し、payload に
    `precursor_exclusion_tol`・`measured_peaks`（測定スペクトル）・`alignments`
    （候補ごとの一致したピーク `[m/z, 測定, 参照]`）を足す——直前の照合の状態
    （`last_match`）を読まずに、戻り値だけで照合の中身を確かめられる。
    **`matched_peaks_count` は前駆体の一致も 1 本に数える**（意味は変えていない）。

    dataset_handle: `dataset_handle_open` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッション。ハンドルのライブラリ・照合結果だけを読み書きする。
        `spot_id` の解決は、`arf2_path` / `dcl_file` を省けばハンドルの `.arf2` / `.dcl` を使う。
    """
    from metabolomix.tools._holder import alignment_paths, resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    store_obj = holder.library.store
    if store_obj is None:
        lbm_only = lbm_only_not_evaluable(holder)
        if lbm_only is not None:
            return json_payload(lbm_only)
        return mcp_errors.missing_state(
            "library_store", ["library_load"],
            "先に library_load を実行してください（参照ライブラリが読み込まれていません）。",
        )
    if (precursor_mz is None) == (spot_id is None):
        return json_payload({
            "status": "invalid_arguments",
            "message": "precursor_mz か spot_id のどちらか一方を指定してください。",
        })
    arf2_path, dcl_file, lacking = alignment_paths(
        holder, dataset_handle, arf2_path, dcl_file,
        need=("arf2", "dcl") if spot_id is not None else ("dcl",))
    if lacking:
        return json_payload(lacking)
    measured = None
    if spot_id is not None:
        spot = _spot_query(spot_id, arf2_path=arf2_path, dcl_file=dcl_file,
                           rt=rt, ion_mode=ion_mode)
        if spot["status"] != "ok":
            return json_payload(spot)
        precursor_mz, rt, ion_mode, measured = (
            spot["precursor_mz"], spot["rt"], spot["ion_mode"], spot["measured"])

    if precursor_exclusion_tol is not None and not precursor_exclusion_tol > 0:
        return json_payload({
            "status": "invalid_arguments",
            "message": f"precursor_exclusion_tol は正の数で指定してください: {precursor_exclusion_tol!r}",
        })

    # `.dcl` 側の窓は library 側の resolved_mz_tol/resolved_rt_tol を渡さない
    # （Important 7 — 上のクラス docstring参照）。_measured_spectrum 自身の
    # 固定既定（dcl/reader.py と同じ mz_tol=0.01/rt_tol=0.2）を使う。
    if measured is None:
        measured = _measured_spectrum(precursor_mz, rt=rt, dcl_file=dcl_file)
    if not measured:
        return json_payload({
            "status": "not_found",
            "message": (
                f"precursor m/z={precursor_mz}"
                + (f", RT={rt}" if rt is not None else "")
                + " に一致する MS/MS が見つかりませんでした。"
                "MS/MS 未取得（.dcl が無い、または該当 precursor の記録が無い）のであって、"
                "期待フラグメントが合わなかったのではありません。"
            ),
        })

    # 許容幅の解決（明示 → store の search_params → 既定）・候補の検索・採点・並べ替えは
    # library/matching.py が唯一の実装。採点の関数はこの module の名前を渡す
    # （既存のテストがここの match_spectrum を差し替えて順位付けを確かめている）。
    matched = matching.score_candidates(
        store_obj, precursor_mz=precursor_mz, rt=rt, ion_mode=ion_mode, measured=measured,
        mz_tol=mz_tol, ms2_tol=ms2_tol, rt_tol=rt_tol, match_fn=match_spectrum)
    resolved_mz_tol = matched["query"]["mz_tol"]
    resolved_ms2_tol = matched["query"]["ms2_tol"]
    if matched["status"] == "no_candidates":
        return json_payload({
            "status": "no_candidates",
            "message": (
                # 許容幅は書式を指定して埋める。`.dbs` 由来の値は C# の 32bit float を
                # 64bit へ広げたもので `0.009999999776482582` になり、しかも
                # `round_floats()` は payload 構造の中の float しか辿らないので
                # **文字列へ焼いた数値には届かない**（成功時の query.mz_tol は
                # 構造の中の float なので丸まる——だからここだけ素の値が出ていた）。
                f"precursor m/z={precursor_mz}±{resolved_mz_tol:.4g}"
                + (f", ion_mode={ion_mode}" if ion_mode else "")
                + " に該当する参照レコードが見つかりませんでした。"
            ),
        })

    # 足切りで採点から外れたピーク（図が「描かれてはいるが採点に入っていない」ピークを
    # 区別するのに使う）と、総合スコアの降順（同点は候補順）に並べた候補。
    unscored_mz = matched["unscored_mz"]
    scored = matched["scored"]
    candidates = scored
    top = scored[:top_n] if top_n >= 0 else scored

    holder.library.last_match = {
        "spot_id": spot_id,
        "precursor_mz": precursor_mz,
        "rt": rt,
        "ion_mode": ion_mode,
        "ms2_tol": resolved_ms2_tol,
        "measured": measured,
        "unscored_mz": unscored_mz,
        "candidates": [
            {
                "name": record["name"],
                "precursor_mz": record["precursor_mz"],
                "ion_mode": record["ion_mode"],
                "adduct": record["adduct"],
                "rt": record["rt"],
                "formula": record["formula"],
                "inchikey": record["inchikey"],
                "smiles": record["smiles"],
                "compound_class": record["compound_class"],
                "ontology": record["ontology"],
                "library_id": record["library_id"],
                "record_index": record["record_index"],
                "spectrum": record["spectrum"],
                "scores": {k: v for k, v in result.items() if k != "alignment"},
                "alignment": result["alignment"],
            }
            for record, result in top
        ],
    }

    rows = [
        {
            "rank": i + 1,
            "name": record["name"],
            "precursor_mz": record["precursor_mz"],
            "ion_mode": record["ion_mode"],
            "adduct": record["adduct"],
            "rt": record["rt"],
            "formula": record["formula"],
            "ontology": record["ontology"],
            "compound_class": record["compound_class"],
            "total_score": result["total_score"],
            "simple_dot_product": result["simple_dot_product"],
            "weighted_dot_product": result["weighted_dot_product"],
            "reverse_dot_product": result["reverse_dot_product"],
            "matched_peaks_percentage": result["matched_peaks_percentage"],
            "matched_peaks_count": result["matched_peaks_count"],
            "entropy_similarity": result["entropy_similarity"],
        }
        for i, (record, result) in enumerate(top)
    ]
    columns = list(_CANDIDATE_TABLE_COLUMNS)
    exclusion_tol = precursor_exclusion_tol or resolved_ms2_tol
    if detail:
        columns += _DETAIL_COLUMNS
        for row, (record, result) in zip(rows, top):
            row["inchikey"] = record["inchikey"]
            row["fragment_matches_excl_precursor"] = fragment_matches_excluding_precursor(
                result["alignment"], precursor_mz, exclusion_tol)
    table = format_spots_as_table(rows, columns=columns)

    payload = {
        "status": "success",
        **identified_by_lbm(holder, spot_id),
        "query": {"spot_id": spot_id, **matched["query"]},
        # 総合スコアの組み立て方。RT 項の有無は戻り値の数値だけでは分からないので、
        # 候補ごとではなく 1 回だけ載せる（内訳の 3 列は session 側にだけ持つ）。
        "scoring": matched["scoring"],
        "measured_peak_count": len(measured),
        "n_candidates": len(candidates),
        "top_n": len(top),
        "ranked_by": _RANK_KEY,
        "candidates_table": table,
    }
    if detail:
        payload["precursor_exclusion_tol"] = exclusion_tol
        payload["measured_peaks"] = [[float(mz), float(intensity)] for mz, intensity in measured]
        payload["alignments"] = [
            {"rank": i + 1,
             "matched": [[row["mz"], row["measured"], row["reference"]]
                         for row in result["alignment"] if row.get("matched")]}
            for i, (_record, result) in enumerate(top)
        ]
    return json_payload(round_floats(payload))


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def library_plot_mirror(rank: int = 1, output: str | None = None,
                        scale: str = mirror_plot.RELATIVE,
                        label_policy: str = mirror_plot.AUTO,
                        dataset_handle: str | None = None) -> list | str:
    """直近の `library_match_feature` 結果から対向プロット（mirror plot）を描く。

    上段が測定、下段が参照（`rank` 位の候補）。`output="image"`（既定）は
    サーバ側で描画した PNG を返し、`output="payload"` は座標を丸めた
    `lipidmix.mirror.v2` JSON を返す（Use-LLLM 等、自前で描くクライアント向け。
    v2 で採点外ピーク `unscored_mz` の層を追加）。

    測定側の一致色分けには許容幅（`ms2_tol`）が要る——`library_match_feature`
    が使った値（store の `search_params` か既定値）を自動で引き継ぐ。

    **`.dbs` の強度足切りで採点に入らなかった測定ピークは灰色で薄く描く**
    （`last_match["unscored_mz"]` を引き継ぐ）。採点は `normalize_measured` が
    足切りした後のスペクトルに対して行うのに、図は足切り前の生ピークを描くため、
    区別が無いと「描かれているのに採点されていない」ピークが黙って混ざる。
    消さずに描くのは、消すと「MS/MS が取れていない」と読めてしまうため。
    この層はラベル枠を取らず、凡例と caption の件数は**該当ピークがあるときだけ**
    出る（足切りが 0 の run では図も caption も従来と変わらない）。

    `scale` は縦軸の写し方で `"relative"`（既定）/ `"sqrt"` / `"log10"`。
    **precursor がベースピークのスペクトル**（脂質の [M-H]- など）は `relative`
    だと診断イオンが相対数 % に潰れて読めない——そのときに `"sqrt"` を使う。
    上流 MS-DIAL も同じ問題を軸の切り替えで解いている
    （`ObservableMsSpectrum.CreateAxisPropertySelectors2` の Relative / Absolute /
    Log10 / Sqrt）。`Absolute` は用意しない——対向プロットは単位の違う 2 つの
    スペクトルを上下に並べるので、生の強度で並べても比較にならない。

    `label_policy` は m/z ラベルの衝突回避の方式。どちらも**強度降順に走査して
    既に置いたラベルと重なるものを飛ばす**（上流 `Annotator.OnRender` と同じ
    貪欲法）が、重なりの見方が違う:

    - `"auto"`（既定）— 水平と垂直の両方が近いときだけ飛ばす。箱は
      ラベル自身の文字列を実測する。対向プロットは同じ側でもピークの高さが
      大きく違うので、2 次元で見るほうが読めるラベルを多く残せる。
    - `"msdial"` — 上流に忠実。**水平距離だけ**を見て、箱は代表文字列
      `"1000.00000"` 1 つで全ラベル共通にする。上流の MS2 ビューは
      `Overlap="Horizontal, Direct"` だが、この合成は OR で、全ラベルが同じ箱を
      使うため `Direct` は `Horizontal` の部分集合になり実質 `Horizontal` 単独に
      縮退している。縦にどれだけ離れていても m/z が近ければ飛ばすので、
      `"auto"` よりラベルは少なくなる。

    どちらも上下は別枠で数える（上流は上下で別々の `Annotator` を持つため、
    側をまたぐ衝突は起きない）。

    dataset_handle: `dataset_handle_open` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッション。ハンドルのライブラリ・照合結果だけを読み書きする。
    """
    from metabolomix.tools._holder import resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    last_match = holder.library.last_match
    if not last_match or not last_match.get("candidates"):
        return mcp_errors.missing_state(
            "library_match", ["library_match_feature"],
            "先に library_match_feature を実行してください（直近の照合結果がありません）。",
        )

    candidates = last_match["candidates"]
    if rank < 1 or rank > len(candidates):
        return json_payload({
            "status": "error",
            "message": f"rank は 1〜{len(candidates)} の範囲で指定してください（受け取った値: {rank!r}）。",
        })

    try:
        mode = plot_render.resolve_plot_output(output)
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})

    candidate = candidates[rank - 1]
    payload = mirror_plot.build_mirror_payload(
        last_match["measured"], candidate["spectrum"], candidate["alignment"],
        title=f"{candidate['name']} (rank {rank})",
        ms2_tol=last_match.get("ms2_tol"),
        unscored_mz=last_match.get("unscored_mz"),
    )

    if mode == plot_render.PAYLOAD:
        return json_payload(round_floats(payload))

    try:
        png = mirror_plot.render_mirror(payload, scale=scale, label_policy=label_policy)
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})
    caption = (
        f"Mirror: measured vs {candidate['name']} "
        f"(rank {rank}, precursor m/z={candidate['precursor_mz']}, "
        f"y-axis={scale}, labels={label_policy})."
    )
    # 足切りで採点から外れたピークがあるときだけ件数を添える。図の灰色の層は
    # 見落としやすく、「描かれているのに採点されていない」は数値の読み違いに直結する。
    if payload["unscored_peak_count"]:
        total_peaks = payload["scored_peak_count"] + payload["unscored_peak_count"]
        caption += (
            f" Scored {payload['scored_peak_count']}/{total_peaks} measured peaks; "
            f"{payload['unscored_peak_count']} below cutoff (grey, not scored)."
        )
    return [caption, Image(data=png, format="png")]


#: library_search の一覧の列。
_SEARCH_COLUMNS = ["name", "precursor_mz", "ion_mode", "adduct", "rt", "formula", "inchikey",
                   "ontology", "compound_class", "n_peaks"]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def library_search(name: str | None = None, formula: str | None = None,
                   inchikey: str | None = None, ion_mode: str | None = None,
                   limit: int = 50, dataset_handle: str | None = None) -> str:
    """読み込んだ参照ライブラリを名前・分子式・InChIKey で引く（読み取りだけ。再解析で使う）。

    name: 大小を問わない部分一致。formula: 完全一致（例 `C42H82NO8P`）。
    inchikey: 先頭 14 文字（構造の骨格のブロック）で照合する（全体を渡してよい）。
    ion_mode: `positive` / `negative`（極性不明のレコードは残す。`library_match_feature` と同じ）。
    条件は与えたものの AND で、1 つも無ければ `invalid_arguments`。`limit`（既定 50）件まで
    TSV（列名 1 回）で返し、打ち切った件数を `truncated` に書く。スペクトルは返さない
    （ピーク数 `n_peaks` だけ）。
    dataset_handle: `dataset_handle_open` が返したハンドル（そのハンドルのライブラリを引く）。
        省略時は対話のセッションの `library_load` のライブラリ。
    """
    from metabolomix.tools._holder import resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    store_obj = holder.library.store
    if store_obj is None:
        return mcp_errors.missing_state(
            "library_store", ["library_load", "dataset_handle_open"],
            "参照ライブラリが読み込まれていません（library_load、またはライブラリつきの"
            "dataset_handle_open）。")
    if name is None and formula is None and inchikey is None:
        return json_payload({"status": "invalid_arguments",
                             "message": "name / formula / inchikey のどれか 1 つは指定してください。"})
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        return json_payload({"status": "invalid_arguments",
                             "message": f"limit は 1 以上の整数で指定してください: {limit!r}"})
    block = inchikey.strip()[:14].upper() if inchikey else None
    records, truncated = store_obj.search(
        name=name.strip() if name else None, formula=formula.strip() if formula else None,
        inchikey_block=block, ion_mode=ion_mode, limit=limit)
    rows = [{**{k: record.get(k) for k in _SEARCH_COLUMNS if k != "n_peaks"},
             "n_peaks": len(record.get("spectrum") or [])} for record in records]
    return json_payload(round_floats({
        "status": "success",
        "query": {"name": name, "formula": formula, "inchikey_block": block, "ion_mode": ion_mode,
                  "limit": limit},
        "n_hits": len(rows),
        "truncated": truncated,
        "records": format_spots_as_table(rows, columns=_SEARCH_COLUMNS) if rows else "",
    }))
