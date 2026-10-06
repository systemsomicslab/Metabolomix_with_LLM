"""mzTab-M 入口 MCP ツール: dataset_load, dataset_status。spec §12 参照。"""
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import session_state
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
                 allow_incomplete: bool = False) -> str:
    """mzTab-M 2.0 ファイルを読み込んで正準状態（DatasetState）を作る。

    呼び出し方:
      - mztab_path: .mzTab ファイルへの絶対パスを直接指定する場合。
      - job_path: analysis-job.json へのパスを指定する場合。ジョブの宣言
        （polarity・measure）で正準エントリを一意に選ぶ。console_run 完了後は
        こちらを推奨。
      - pipeline_path: `pipeline_run` が返した pipeline-run.json（またはその親
        ディレクトリ）を指定する場合。その run の mzTab に加えて、**解決済みの
        試料対応表・feature binding・解析行列も同じ run からまとめて**載せる。
        行列だけを別に載せる入口は置いていない——別 run の行列を混ぜると群の
        対応が黙ってずれるため。載せた行列は `dataset_statistic` へ
        `matrix_id` で渡す。

    3つとも省略するとエラー、2つ以上指定するとエラー。

    allow_incomplete: **完了していない実行**（partial / failed / running）の出力を
        読み込みます（既定 False）。中断時点の生成物なので、既定では拒否します
        （`INCOMPLETE_ANALYSIS_JOB`）。True で読んだデータセットは探索専用になり、
        2 群比較（dataset_differential）と差次的エクスポートは拒否されます——
        欠けた検体を「その群には無い」と読み違えたまま結論が出てしまうためです。

    成功すると session.dataset に DatasetState が格納され、
    dataset_status で内容を確認できる。
    既存の session.arf / .arf2 / .pai2 / .eic は変更しない。
    """
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.mztab.loading import load_dataset_state
    from metabolomix.pipeline import session_handoff

    try:
        if pipeline_path:
            if mztab_path or job_path:
                raise DomainError(
                    "DATASET_BAD_REQUEST",
                    "pipeline_path は mztab_path / job_path と同時に指定できません"
                    "（run の成果物はまとめて載せます）。",
                    {"pipeline_path": pipeline_path})
            ds = session_handoff.load_run_dataset(
                pipeline_path, allow_incomplete=allow_incomplete)
        else:
            ds = load_dataset_state(
                mztab_path=Path(mztab_path) if mztab_path else None,
                job_path=Path(job_path) if job_path else None,
                allow_incomplete=allow_incomplete)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    # 成功したときだけ session を差し替える。途中で失敗した読み込みが
    # 現在のデータセットを壊さない。
    session_state.session.dataset = ds

    filename = Path(ds.source_files and list(ds.source_files)[-1]
                    or (mztab_path or "")).name or "mzTab-M"
    lines = [_summary_text(ds, filename)]
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
    return "\n".join(lines)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_status() -> str:
    """現在の DatasetState の概要を返す。

    `samples` に name/role の TSV が入る。dataset_differential の
    group_a / group_b はここに出る名前をそのまま使う。

    dataset_load を先に実行しておくこと。
    未実行の場合は missing_state エンベロープを返す。
    """
    ds = session_state.session.dataset
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
