"""参照ライブラリ（MS/MS スペクトル照合）の MCP ツール群。

下の層（`lipidmix.library.dbs` / `.msp` / `.store`、`lipidmix.analysis.spectral_match`、
`lipidmix.plots.mirror`）はすべて完成済み。ここはそれらを MCP の公開面へ繋ぐだけの層で、
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

from lipidmix.analysis.spectral_match import cutoff_mask, match_spectrum, total_score
from lipidmix.arf2.reader import format_spots_as_table
from lipidmix.core import mcp_errors, session_state, user_config
from lipidmix.core.mcp_core import mcp
from lipidmix.core.path_resolvers import (
    LIBRARY_SETTING_KEYS,
    LibraryPathError,
    resolve_dcl_file_path,
    resolve_library_path,
)
from lipidmix.core.serialization import json_payload, round_floats
from lipidmix.dcl.reader import deserialize_dcl, get_msms_by_precursor
from lipidmix.library.defaults import DEFAULT_MS2_TOL as _DEFAULT_MS2_TOL
from lipidmix.library.defaults import DEFAULT_MZ_TOL as _DEFAULT_MZ_TOL
from lipidmix.library.defaults import DEFAULT_RT_TOL as _DEFAULT_RT_TOL
from lipidmix.library.defaults import pick_tol as _pick_tol
from lipidmix.library.store import open_store
from lipidmix.plots import mirror as mirror_plot
from lipidmix.plots import render as plot_render

__all__ = ["library_load", "library_match_feature", "library_plot_mirror"]

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


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def library_load(file_path: str | None = None, rebuild: bool = False,
                 ion_mode: str | None = None) -> str:
    """参照ライブラリを読み込み、照合用の SQLite store を構築（または既存キャッシュを
    再利用）する。

    解決順: `file_path` の明示 → `ion_mode`（`"positive"` / `"negative"`）に対応する
    設定（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の
    `[library] msp_positive` / `msp_negative`）→ データディレクトリの
    `*_Loaded.msp2.dbs` → 設定済みの極性（両方あれば `ion_mode` を求める）→
    データディレクトリの `*.msp`。候補が 1 つに決まらなければ `code` 付きの
    エラー（`MSP_AMBIGUOUS` など）を返し、黙って選ばない。一度に保持する
    ライブラリは 1 つだけ——測定の極性に合わせて `ion_mode` を指定する。

    初回の構築は大きな `.msp`（1 GB 級）で数分かかる。MCP のタイムアウトに当たる
    なら、手元で `python -m lipidmix.library.store --ion-mode <極性>` を 1 度
    走らせておけば、以後はキャッシュを開くだけで済む。

    `library_match_feature` の前提。別ライブラリへ切り替えると直近の照合結果
    （`library_plot_mirror` が読む座標）は破棄する——古い照合を新ライブラリの
    結果と取り違えないため。

    `rebuild=True` はキャッシュを無視して構築し直す（元ファイルが壊れている
    疑いがあるときなど）。通常は不要——store は元ファイルの sha256 をキーに
    キャッシュされるため、内容が変わらない限り再構築しない。
    """
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
        store_obj = open_store(resolved, rebuild=rebuild)
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
    old_store = session_state.session.library.store
    if old_store is not None:
        old_store.close()

    # 別ライブラリへの切り替え起こりうるので、古い照合結果を新ライブラリのものと
    # 取り違えないよう破棄する（ArfState.reset_analysis と同じ考え方）。
    session_state.session.library.store = store_obj
    session_state.session.library.source_path = resolved
    session_state.session.library.last_match = None

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


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def library_match_feature(
    precursor_mz: float,
    rt: float | None = None,
    ion_mode: str | None = None,
    dcl_file: str | None = None,
    mz_tol: float | None = None,
    ms2_tol: float | None = None,
    rt_tol: float | None = None,
    top_n: int = 5,
) -> str:
    """測定 MS/MS（`.dcl`）を参照ライブラリ（`library_load` 済みの store）と照合し、
    上位候補を返す。

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
    """
    store_obj = session_state.session.library.store
    if store_obj is None:
        return mcp_errors.missing_state(
            "library_store", ["library_load"],
            "先に library_load を実行してください（参照ライブラリが読み込まれていません）。",
        )

    search_params = store_obj.summary().get("search_params") or {}
    resolved_mz_tol = _pick_tol(mz_tol, search_params, "ms1_tolerance", _DEFAULT_MZ_TOL)
    resolved_ms2_tol = _pick_tol(ms2_tol, search_params, "ms2_tolerance", _DEFAULT_MS2_TOL)
    resolved_rt_tol = _pick_tol(rt_tol, search_params, "rt_tolerance", _DEFAULT_RT_TOL)
    # match_spectrum の採点前処理パラメータ。検証 CLI（scripts/verify_spectral_match.py）
    # と同じ集合を search_params から採る（Important 2）。既定は
    # MsRefSearchParameterBase の Key 0/1/7/8（`dbs.py` の `_SEARCH_PARAM_KEYS`）。
    mass_begin = _pick_tol(None, search_params, "mass_range_begin", 0.0)
    mass_end = _pick_tol(None, search_params, "mass_range_end", 2000.0)
    relative_amp_cutoff = _pick_tol(None, search_params, "relative_amp_cutoff", 0.0)
    absolute_amp_cutoff = _pick_tol(None, search_params, "absolute_amp_cutoff", 0.0)

    # `.dcl` 側の窓は library 側の resolved_mz_tol/resolved_rt_tol を渡さない
    # （Important 7 — 上のクラス docstring参照）。_measured_spectrum 自身の
    # 固定既定（dcl/reader.py と同じ mz_tol=0.01/rt_tol=0.2）を使う。
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

    candidates = store_obj.candidates(
        precursor_mz, mz_tol=resolved_mz_tol, ion_mode=ion_mode, rt=rt, rt_tol=resolved_rt_tol,
    )
    if not candidates:
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

    # 足切りで採点から外れたピーク。候補ごとに変わらない（測定スペクトルと cutoff
    # だけで決まる）ので、候補ループの外で 1 回だけ求める。図はこの一覧を使って
    # 「描かれてはいるが採点に入っていない」ピークを区別する。
    kept = cutoff_mask(measured, relative_amp_cutoff=relative_amp_cutoff,
                       absolute_amp_cutoff=absolute_amp_cutoff)
    unscored_mz = [float(peak[0]) for peak, keep in zip(measured, kept) if not keep]

    # RT 項を入れるかは `.dbs` の `IsUseTimeForAnnotationScoring`（Key 16）次第。
    # `.msp` のようにフラグを持たないライブラリは上流の既定（False）に倣う。
    use_rt_scoring = bool(search_params.get("use_time_for_annotation_scoring", False))

    scored = []
    for record in candidates:
        result = match_spectrum(
            measured, record["spectrum"], ms2_tol=resolved_ms2_tol,
            mass_begin=mass_begin, mass_end=mass_end,
            relative_amp_cutoff=relative_amp_cutoff, absolute_amp_cutoff=absolute_amp_cutoff,
        )
        result.update(total_score(
            result,
            precursor_mz=precursor_mz, reference_precursor_mz=record.get("precursor_mz"),
            ms1_tol=resolved_mz_tol,
            rt=rt, reference_rt=record.get("rt"),
            rt_tol=resolved_rt_tol, use_rt=use_rt_scoring,
        ))
        scored.append((record, result))
    # ランク付けの基準は _RANK_KEY（total_score、MS-DIAL の総合スコア）。
    # 同点はそのまま候補順で温存する（安定ソート）。docstring と payload["ranked_by"]
    # がこの基準を呼び出し側へ明示する（コード内コメントだけに留めない）。
    scored.sort(key=lambda item: item[1][_RANK_KEY], reverse=True)
    top = scored[:top_n] if top_n >= 0 else scored

    session_state.session.library.last_match = {
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
    table = format_spots_as_table(rows, columns=_CANDIDATE_TABLE_COLUMNS)

    payload = {
        "status": "success",
        "query": {
            "precursor_mz": precursor_mz, "rt": rt, "ion_mode": ion_mode,
            "mz_tol": resolved_mz_tol, "ms2_tol": resolved_ms2_tol, "rt_tol": resolved_rt_tol,
            "mass_begin": mass_begin, "mass_end": mass_end,
            "relative_amp_cutoff": relative_amp_cutoff, "absolute_amp_cutoff": absolute_amp_cutoff,
        },
        # 総合スコアの組み立て方。RT 項の有無は戻り値の数値だけでは分からないので、
        # 候補ごとではなく 1 回だけ載せる（内訳の 3 列は session 側にだけ持つ）。
        "scoring": {
            "rule": "MS-DIAL GetTotalScore: rt + precursor_mz + (weighted+simple+reverse)/3 "
                    "+ matched_peaks_percentage（正規化しない和。各項は > 0 のときだけ加算）",
            "use_rt": use_rt_scoring,
            "ms1_tol": resolved_mz_tol,
            "rt_tol": resolved_rt_tol if use_rt_scoring else None,
        },
        "measured_peak_count": len(measured),
        "n_candidates": len(candidates),
        "top_n": len(top),
        "ranked_by": _RANK_KEY,
        "candidates_table": table,
    }
    return json_payload(round_floats(payload))


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def library_plot_mirror(rank: int = 1, output: str | None = None,
                        scale: str = mirror_plot.RELATIVE,
                        label_policy: str = mirror_plot.AUTO) -> list | str:
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
    """
    last_match = session_state.session.library.last_match
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
