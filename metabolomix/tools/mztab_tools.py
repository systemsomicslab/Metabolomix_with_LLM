"""mzTab-M 入口 MCP ツール: dataset_load, dataset_status。spec §12 参照。"""
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import handles, session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.mztab.reader import parse_mztab
from metabolomix.mztab import evidence as mztab_evidence
from metabolomix.mztab.dataset_state import build_dataset_state
# 生成物パスの解決は loading.py が正準。console_cleanup など既存の呼び出しが
# ここを参照しているため、名前だけ再エクスポートする（実装は 1 か所）。
from metabolomix.mztab.loading import artifact_abs_path as _artifact_abs_path

__all__ = ["dataset_load", "dataset_status"]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def dataset_load(mztab_path: str | None = None, job_path: str | None = None,
                 pipeline_path: str | None = None,
                 allow_incomplete: bool = False, directory: str | None = None,
                 source: str = "auto", measure: str = "peak_height") -> str:
    """データセットを読み込んで正準状態（DatasetState）を作る。

    入口（**どれか 1 つ**。全省略・2 つ以上の指定はエラー）:
      - directory: MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata`）。フォルダの中の
        特定の `.arf` / `.arf2` を渡してもよい。同じフォルダに `DriftSpots.arf` と
        `PeakProperties.arf` があれば `PeakProperties.arf`、複数のバッチがあれば最新を選び、
        選ばなかった `.arf` を要約に並べる。
      - mztab_path: .mzTab ファイルへの絶対パスを直接指定する場合。
      - job_path: analysis-job.json へのパスを指定する場合。ジョブの宣言
        （polarity・measure）で正準 mzTab-M を一意に選ぶ。
      - pipeline_path: `pipeline_run` が返した pipeline-run.json（またはその親
        ディレクトリ）。

    source: どのアダプタで読むか。`"auto"`（既定）/ `"folder"` / `"mztab"`。
        - `directory` は解析フォルダ（`"mztab"` はエラー）、`mztab_path` は mzTab-M
          （`"folder"` はエラー）。
        - `job_path` / `pipeline_path` は既定（`"auto"` / `"folder"`）で、そのジョブの正準
          mzTab-M と同じアラインメントの**解析フォルダ**を読む。候補（mzTab-M の隣の `.arf`。
          pipeline はどれも合わなければ run の `input/AlignResult-*.arf2` と同じバッチも）を
          mzTab-M と特徴数・m/z・試料名で照合し、合うものがちょうど 1 本ならそれを読む。
          無ければ `ALIGNMENT_NOT_FOUND`、2 本以上なら `ALIGNMENT_AMBIGUOUS`（どちらも
          候補のファイル名を `details.candidates` に載せる）で止まり、mzTab-M へ黙って
          切り替えない。`"mztab"` を渡すと今までどおり mzTab-M を読む。
        - **pipeline の run を解析フォルダで読んだときに run から引き継ぐのは、解決済みの
          試料対応表だけ**。run の解析行列・feature binding・差次的結果は
          `source="mztab"` のときだけ載せる——その run の mzTab に加えて、解決済みの
          試料対応表・feature binding・解析行列を**同じ run からまとめて**載せる
          （行列だけを別に載せる入口は置いていない。別 run の行列を混ぜると群の対応が
          黙ってずれるため）。載せた行列は `dataset_statistic` へ `matrix_id` で渡す。
          run の差次的結果も、データと群の指紋が今のデータセットと一致するものだけを戻す
          （result_id は run の TSV と同じ。直近の結果にはしないので result_id で名指しして
          使う）。戻せなかったものは理由つきで要約に並ぶ。
    measure: 解析フォルダのアダプタが値に使う `.arf` の欄（`"peak_height"`（既定）/
        `"peak_area_above_zero"`）。**この引数が効くのは解析フォルダのアダプタだけ**。
        mzTab-M のアダプタでは定量種別がファイル（job / pipeline ならジョブの宣言）で
        決まるので、既定値以外を渡すと `DATASET_BAD_REQUEST`。

    allow_incomplete: **完了していない実行**（partial / failed / running）の出力を
        読み込みます（既定 False。bool 以外は `DATASET_BAD_REQUEST`——`"false"` という
        文字列を真と読まないため。どちらのアダプタでも同じ門を通る）。中断時点の生成物
        なので、既定では拒否します（`INCOMPLETE_ANALYSIS_JOB`）。True で読んだ
        データセットは探索専用になり、2 群比較（dataset_differential）と差次的
        エクスポートは拒否されます——欠けた検体を「その群には無い」と読み違えたまま
        結論が出てしまうためです。

    成功すると DatasetState を新しいデータセットのハンドル（`dh_…`）として登録簿に入れ、
    対話のセッションの直近にする（要約の `dataset_handle`）。`dataset_handle` を省いた
    `dataset_*` / `library_*` はこの直近のハンドルに効き、dataset_status で内容を確認できる。
    前に読んだデータセットのハンドルは閉じずに残る（`dataset_handle_list` / `dataset_handle_close`）。
    開いているハンドルが上限（`LIPIDMIX_MAX_DATASET_HANDLES`、既定 8）なら
    `DATASET_HANDLE_LIMIT` で止め、古いハンドルを黙って閉じない（直近も動かさない）。
    要約には読んだアダプタ・選んだ `.arf`・
    同定に使ったライブラリの記録（ファイル名と sha256。置き場所は出さない）が並ぶ。
    同じ入力を読み直すと読み込み済みのものの複製を使う（読み込みが読んだ・探したファイル
    ——mzTab-M の隣に後から置かれた `.arf` なども——の size・更新時刻が同じとき）。
    既存の session.arf / .arf2 / .pai2 / .eic は変更しない。
    """
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.dataset.opening import open_dataset

    try:
        # 長い読み込みの前に上限を確かめる。読み込みの間に別の呼び出しが足しうるので、
        # 直近にする代入（registry.add）でももう一度確かめる。
        handles.registry.check_capacity()
        opened = open_dataset(directory=directory, mztab_path=mztab_path, job_path=job_path,
                              pipeline_path=pipeline_path, source=source, measure=measure,
                              allow_incomplete=allow_incomplete)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    # 成功したときだけ直近を差し替える。途中で失敗した読み込み（上限を含む）が
    # 現在のデータセットを壊さない。
    ds = opened.dataset
    try:
        session_state.session.dataset = ds
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    handle_id = session_state.session.current_handle_id
    # ハンドルに持たせる `.arf2` / `.dcl`（`dataset_handle_open` と同じ規則）。ID を明示した
    # 呼び出しはこれを使う。ID を省いた呼び出しは `_holder.alignment_paths` とデータフォルダの
    # 探索で決めるので、同じファイルになるとは限らない。
    from metabolomix.tools.handle_tools import alignment_files_for

    handle = handles.registry.peek(handle_id)
    if handle is not None:
        handle.arf2_path, handle.dcl_path = alignment_files_for(opened)

    if opened.selection is not None:
        lines = [_folder_summary_text(ds, opened.selection.arf.name)]
    else:
        filename = Path(ds.source_files and list(ds.source_files)[-1]
                        or (mztab_path or "")).name or "mzTab-M"
        lines = [_summary_text(ds, filename)]
    lines.append(f"- dataset_handle: {handle_id}（直近のデータセット。省略の呼び出しが使う）")
    lines.append(f"- 出所の裏取り: {ds.source_verification}"
                 + ("（探索専用: 2 群比較と差次的エクスポートは拒否されます）"
                    if ds.exploratory_only else ""))
    if ds.analysis_matrices:
        lines.append("- 解析行列（analysis-matrix.v1）: "
                     + ", ".join(sorted(ds.analysis_matrices)))
    if ds.artifact_paths:
        lines.append("- 付随アーティファクト: "
                     + str({role: len(paths)
                            for role, paths in ds.artifact_paths.items()}))
    lines.extend(_handoff_lines(ds.handoff_report))
    lines.extend(_opening_lines(opened))
    return "\n".join(lines)


def _opening_lines(opened) -> list[str]:
    """どのアダプタで何を読んだか・同定に使ったライブラリの記録（置き場所は出さない）。"""
    from metabolomix.dataset.library_memory import public_view

    lines = [f"- アダプタ: {opened.adapter}"]
    if opened.selection is not None:
        lines.append(f"- 読んだ .arf: {opened.selection.arf.name}")
        if opened.selection.not_chosen:
            lines.append("- 選ばなかった .arf: " + ", ".join(opened.selection.not_chosen))
    lines.append("- 同定に使ったライブラリ（library_record）: "
                 + json_payload(public_view(opened.dataset.library_record)))
    lines.extend(f"- 注記: {note}" for note in opened.notes)
    return lines


#: 戻せなかった理由ごとの案内（`session_handoff.restore_differentials` の reason）。
#: `pipeline_resume` には完了済みの差次的工程だけを再計算させる手順が無い
#: （要求の変更か `rerun_upstream=true` でしか差し戻らない）ので、旧 run は
#: session で計算し直すよう案内する。
_HANDOFF_REASON_NOTES = {
    "legacy_result_without_dataset_hash":
        "指紋を持たない旧 run の結果。session で dataset_preprocess → "
        "dataset_differential を実行して計算し直してください。",
    "dataset_mismatch":
        "run の計算時とデータ（mzTab）が一致しません。今読み込んだデータで使うには、"
        "session で dataset_preprocess → dataset_differential を実行して計算し直してください。",
    "group_mismatch":
        "run の計算時と群割当が一致しません。実験情報シート（または "
        "dataset_set_sample_metadata で与えた群）を確かめるか、今の群のまま session で "
        "dataset_preprocess → dataset_differential を実行して計算し直してください。",
}


def _handoff_lines(report: dict | None) -> list[str]:
    """pipeline run から引き継いだ差次的結果の要約（TSV は列名 1 回）。"""
    if not report:
        return []
    lines: list[str] = []
    restored = report.get("restored") or []
    rejected = report.get("rejected") or []
    if restored:
        lines.append("- 引き継いだ差次的結果（result_id で名指しして使う。"
                     "直近の結果 last_differential にはしない）:")
        lines.append("comparison_id\tresult_id\tgroup_a\tgroup_b")
        lines.extend(f"{r['comparison_id']}\t{r['result_id']}\t{r['group_a']}\t{r['group_b']}"
                     for r in restored)
        if any(r.get("evidence_attached_after") for r in restored):
            lines.append("  - evidence_attached_after: run の計算は検出状態（.arf の gap-fill）を"
                         "取り込む前のもの（結果に evidence_attached_after: true）。"
                         "検出率による絞り込みは、その計算では掛かっていない。")
    if rejected:
        lines.append("- 引き継がなかった差次的結果:")
        lines.append("comparison_id\tresult_id\treason")
        lines.extend(f"{r['comparison_id']}\t{r['result_id']}\t{r['reason']}"
                     for r in rejected)
        for reason in sorted({r["reason"] for r in rejected}):
            lines.append(f"  - {reason}: {_HANDOFF_REASON_NOTES.get(reason, '')}")
    return lines


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_status(dataset_handle: str | None = None) -> str:
    """現在の DatasetState の概要を返す。

    `samples` に name/role の TSV が入る。dataset_differential の
    group_a / group_b はここに出る名前をそのまま使う。

    dataset_load を先に実行しておくこと。
    未実行の場合は missing_state エンベロープを返す。

    dataset_handle: `dataset_handle_open` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
    """
    from metabolomix.tools._holder import resolve_holder

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state(
            "dataset",
            ["dataset_load"],
            "DatasetState がありません。先に dataset_load を実行してください。"
            "console_run 完了後は dataset_load(job_path=<analysis-job.json へのパス>) を推奨します。",
        )
    from metabolomix.core import version
    update = version.update_status()
    payload = {
        # 更新後に MCP サーバを再起動し忘れると古いプロセスが黙って動き続ける。
        # 「おかしい」と思ったときに 1 回で分かるよう、ここに刻む。
        "server_version": version.server_version(),
        # 最新なら鍵ごと出さない（戻り値は LLM の文脈をそのまま食う）。
        **({"update_available": update} if update else {}),
        "source_format": ds.source_format,
        "source_files": list(ds.source_files.keys()),
        "quantification_measure": ds.quantification_measure,
        "quantification_confidence": ds.quantification_confidence,
        "n_features": len(ds.feature_ids),
        "n_samples": len(ds.sample_names),
        "validation": {
            "ok": ds.validation_result.get("ok"),
            "n_errors": len(ds.validation_result.get("errors", [])),
            "n_warnings": len(ds.validation_result.get("warnings", [])),
        },
        "inchikey_coverage": ds.inchikey_coverage,
        # 群指定（dataset_differential）に必要なサンプル名。件数だけ返していた頃は、
        # 名前を知る手段が「わざと群サイズ不足のエラーを起こして details を読む」
        # しか無かった。行が並ぶ一覧なので TSV（列名 1 回）で返す。
        "samples": _samples_tsv(ds),
        # 検出状態の有無。「無い」と「全部未検出」は解釈が正反対なので、
        # 率だけを返して available を省くことはしない。
        "detection": _detection_summary(ds),
    }
    if ds.job_path:
        payload["job_path"] = ds.job_path
    if ds.artifact_paths:
        payload["artifact_roles"] = {role: len(paths) for role, paths in ds.artifact_paths.items()}
    return json_payload(payload)


# ---------- 内部ヘルパ ----------


def _detection_summary(ds) -> dict:
    """検出状態（gap-fill）の要約。行列そのものは返さない（42,840 セル規模）。"""
    qc = ds.feature_qc or {}
    if ds.detected_mask is None:
        return {"available": False,
                "reason": qc.get("reason", "no_candidate"),
                "note": ("mzTab-M の非ゼロ値には gap-fill による補間値が含まれ得ます。"
                         "検出状態が取り込めていないため、検出率・欠測率は語れません。")}
    return {"available": True,
            "source": qc.get("source"),
            "n_cells": qc.get("n_cells"),
            "n_detected": qc.get("n_detected"),
            "gap_filled_rate": qc.get("gap_filled_rate")}

def _samples_tsv(ds) -> str:
    """サンプル名と役割を TSV で返す（列名 1 回 + 1 行 1 サンプル）。

    前処理済みなら ds.roles を、まだなら同じ判定関数（detect_sample_roles）を
    その場で適用する。ARF 経路（arf_list_sample_roles）と同じ判定を使うので、
    どちらの経路でも同じ試料が同じ役割になる。
    """
    from metabolomix.analysis import preprocessing
    names = list(ds.sample_names)
    roles = ds.roles or preprocessing.detect_sample_roles(names)
    lines = ["name\trole"]
    lines += [f"{n}\t{roles.get(n, 'sample')}" for n in names]
    return "\n".join(lines)


def _describe(entries) -> list[dict]:
    return [{"path": e.path, "polarity": e.polarity, "measure": e.measure}
            for e in entries]


def _validate_or_error(parse_result: dict, filename: str) -> dict | str:
    """バリデーションを実行し、失敗ならエラーエンベロープ文字列を返す。成功なら validation dict。"""
    from metabolomix.mztab.validator import validate_mztab
    result = validate_mztab(parse_result)
    if not result["ok"]:
        return mztab_error(
            "MZTAB_STRUCTURE_INVALID",
            f"mzTab-M 構造不正 ({filename}): " + "; ".join(result["errors"]),
            {"errors": result["errors"], "warnings": result["warnings"]},
        )
    return result


def _detection_line(ds) -> str:
    qc = ds.feature_qc or {}
    if ds.detected_mask is None:
        return ("- 検出状態: 取り込めていません"
                f"（{qc.get('reason', 'no_candidate')}）。非ゼロ値に gap-fill 補間が"
                "混在し得るため検出率・欠測率は語れません")
    rate = qc.get("gap_filled_rate")
    percent = f"{rate * 100:.1f}%" if rate is not None else "不明"
    return (f"- 検出状態: {qc.get('source')} 由来 — 実測 {qc.get('n_detected')}/"
            f"{qc.get('n_cells')} セル（gap-fill {percent}）")


def _identification_line(ds) -> str:
    """同定の出所内訳。MS1 注釈と MS/MS 裏付けを入口で区別する。

    MCP サーバ指示の「`analytical_checks.msms.band` の PASS と FLAG_ONLY を
    同一視するな」と同じ趣旨。SML 由来の注釈は m/z 照合だけで付く。
    """
    by = (ds.inchikey_coverage or {}).get("identified_by") or {}
    return (f"- 同定: MS/MS 証拠あり {by.get('sme', 0)} 件 / "
            f"MS1 注釈のみ {by.get('sml_only', 0)} 件 / "
            f"同定なし {by.get('none', 0)} 件"
            "（MS1 注釈は m/z 照合のみ。MS/MS 裏付けと同一視しないこと）")


def _summary_text(ds, filename: str) -> str:
    lines = [
        f"## dataset_load 完了: {filename}",
        f"- source_format: {ds.source_format}",
        f"- 特徴量数: {len(ds.feature_ids)}",
        f"- サンプル数: {len(ds.sample_names)}",
        f"- 定量種別: {ds.quantification_measure or '不明'} ({ds.quantification_confidence})",
        f"- 検証: {'OK' if ds.validation_result['ok'] else 'NG — ' + '; '.join(ds.validation_result['errors'])}",
        f"- InChIKey 付き: {ds.inchikey_coverage.get('with_inchikey', 0)}/{ds.inchikey_coverage.get('total_features', 0)} 件",
        _identification_line(ds),
        # 検出状態は入口で言う。実データでは 70% のセルが gap-fill だったので、
        # これを知らずに非ゼロを検出と数えると検出率を 3 倍以上に過大評価する。
        _detection_line(ds),
    ]
    if ds.validation_result.get("warnings"):
        lines.append(f"- 警告 {len(ds.validation_result['warnings'])} 件: " +
                     "; ".join(ds.validation_result["warnings"][:3]))
    return "\n".join(lines)


def _folder_summary_text(ds, filename: str) -> str:
    """解析フォルダのアダプタで読んだときの要約（mzTab-M の検証・InChIKey 表の行は無い）。"""
    n_features = len(ds.feature_ids)
    with_inchikey = sum(1 for meta in ds.feature_metadata.values() if meta.get("inchikey"))
    lines = [
        f"## dataset_load 完了: {filename}",
        f"- source_format: {ds.source_format}",
        f"- 特徴量数: {n_features}",
        f"- サンプル数: {len(ds.sample_names)}",
        f"- 定量種別: {ds.quantification_measure}（.arf の試料ごとの値）",
        f"- InChIKey 付き: {with_inchikey}/{n_features} 件",
        _folder_identification_line(ds),
        _folder_detection_line(ds),
    ]
    if ds.load_notes:
        lines.append(f"- 読み込みの注記 {len(ds.load_notes)} 件: " + "; ".join(ds.load_notes[:3]))
    return "\n".join(lines)


def _folder_identification_line(ds) -> str:
    """同定の出所内訳（mzTab-M の SME / SML / なしと同じ 3 つに写したもの）。"""
    sme = {fid for fid, candidates in ds.feature_candidates.items() if candidates}
    sml_only = {fid for fid, annotation in ds.feature_annotations.items()
                if fid not in sme and (annotation or {}).get("name")}
    return (f"- 同定: MS/MS 証拠あり {len(sme)} 件 / "
            f"MS1 注釈のみ {len(sml_only)} 件 / "
            f"同定なし {len(ds.feature_ids) - len(sme) - len(sml_only)} 件"
            "（MS1 注釈は m/z 照合のみ。MS/MS 裏付けと同一視しないこと）")


def _folder_detection_line(ds) -> str:
    mask = ds.detected_mask
    if mask is None or not mask.size:
        return "- 検出状態: 無し（.arf にセルが無い）"
    n_cells, n_detected = int(mask.size), int(mask.sum())
    return (f"- 検出状態: .arf の gap-fill の印 — 実測 {n_detected}/{n_cells} セル"
            f"（gap-fill・欠測 {(1.0 - n_detected / n_cells) * 100:.1f}%）")
