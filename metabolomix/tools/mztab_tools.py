"""データセットの状態 MCP ツール: dataset_status と、`load_dataset` の要約の部品。spec §12 参照。"""
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

__all__ = ["dataset_status"]


def load_summary_lines(ds, opened) -> list[str]:
    """`load_dataset` の戻り値の `summary`（特徴・試料の数、定量、同定の内訳、検出状態、読み込みの注記）。

    解析フォルダのアダプタと mzTab-M のアダプタで書き分ける。出所の裏取り・探索専用の印・
    解析行列・付随アーティファクト・pipeline の run から引き継いだ差次的結果もここに載る。
    """
    if opened.selection is not None:
        lines = _folder_summary_text(ds, opened.selection.arf.name).splitlines()
    else:
        filename = Path(ds.source_files and list(ds.source_files)[-1]
                        or opened.path or "").name or "mzTab-M"
        lines = _summary_text(ds, filename).splitlines()
    lines.append(f"- 出所の裏取り: {ds.source_verification}"
                 + ("（探索専用: 2 群比較と差次的エクスポートは拒否されます）"
                    if ds.exploratory_only else ""))
    if ds.artifact_paths:
        lines.append("- 付随アーティファクト: "
                     + str({role: len(paths)
                            for role, paths in ds.artifact_paths.items()}))
    lines.extend(_handoff_lines(ds.handoff_report))
    return lines


def _metadata_rows(ds) -> list[dict]:
    """試料の対応表の行（`ds.sample_names` の順）。

    適用済みの `ds.sample_metadata_rows`（解決済みの行）があればそれ、無ければ読み込みが決めた
    既定の行（`sample_manifest.resolve_metadata`）。
    """
    from metabolomix.analysis import sample_manifest

    return getattr(ds, "sample_metadata_rows", None) or sample_manifest.resolve_metadata(ds, None)


def sample_group_role_counts(ds) -> dict:
    """試料の群（Class ID）ごと・役割ごとの件数。`load_dataset` の戻り値と `dataset_status` の要約が共有する。

    群と役割は試料の対応表（`_metadata_rows`）から数える。群の宣言が無い試料は `n_ungrouped`。
    """
    rows = _metadata_rows(ds)
    groups: dict[str, int] = {}
    roles: dict[str, int] = {}
    n_ungrouped = 0
    for row in rows:
        group = row.get("group")
        if group:
            groups[group] = groups.get(group, 0) + 1
        else:
            n_ungrouped += 1
        role = row.get("role") or "sample"
        roles[role] = roles.get(role, 0) + 1
    return {"groups": groups, "roles": roles, "n_ungrouped": n_ungrouped}


#: 戻せなかった理由ごとの案内（`session_handoff.restore_differentials` の reason）。
#: `pipeline_resume` には完了済みの差次的工程だけを再計算させる手順が無い
#: （要求の変更か `rerun_upstream=true` でしか差し戻らない）ので、旧 run は
#: session で計算し直すよう案内する。
_HANDOFF_REASON_NOTES = {
    "legacy_result_without_dataset_hash":
        "指紋を持たない旧 run の結果。session で dataset_preprocess → "
        "dataset_compare を実行して計算し直してください。",
    "dataset_mismatch":
        "run の計算時とデータ（mzTab）が一致しません。今読み込んだデータで使うには、"
        "session で dataset_preprocess → dataset_compare を実行して計算し直してください。",
    "group_mismatch":
        "run の計算時と群割当が一致しません。実験情報シート（または "
        "dataset_set_sample_metadata で与えた群）を確かめるか、今の群のまま session で "
        "dataset_preprocess → dataset_compare を実行して計算し直してください。",
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
def dataset_status(dataset_handle: str | None = None, view: str = "summary",
                   query: str | None = None) -> str:
    """現在の DatasetState の状態を返す。`view` で何を見るかを選ぶ。

    view:
      - `"summary"`（既定）: 概要。`samples` に name/role の TSV が入る。群（Class ID）ごとの件数 `groups` に
        dataset_compare の reference_group / test_group の群名が並ぶ（そのまま使う）・
        役割ごとの件数 `roles`・群の無い試料の数 `n_ungrouped`・特徴のタグごとの件数 `tags`・
        除外の件数 `excluded`（`dataset_exclude`）を含む。
      - `"samples"`: 試料の TSV（name・group・role・batch・injection_order・subject_id・excluded。
        `excluded` は除外の理由、除外していなければ空）。`query` で絞れる。
      - `"tags"`: タグごとの特徴 ID（各タグ先頭 50 件まで。`n_features` は全件数）。`query` に
        タグ名を渡すとそのタグだけ。
    query: `view="samples"` では因子トークン。`_` 区切りで全トークンを含む試料に一致（順不同の AND）し、
        `,` で区切った複数の指定は OR。トークンは試料名・群（Class ID）から取る（日付や拡張子の
        処理タイムスタンプは除く）。QC・ブランクも落とさない（role 列で見分ける）。例 `"ILG_6h"`、
        `"ILG_6h,control_6h"`。1 件も一致しなければ使えるトークンの語彙を添えて拒む。

    load_dataset を先に実行しておくこと。
    未実行の場合は missing_state エンベロープを返す。

    dataset_handle: `load_dataset` が返したハンドル（再解析の並行のエージェント向け）。
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
            ["load_dataset"],
            "DatasetState がありません。先に load_dataset を実行してください。"
            "console_run 完了後は load_dataset(job_path=<analysis-job.json へのパス>) を推奨します。",
        )
    if view not in VIEWS:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"view は {list(VIEWS)} のどれかです: {view!r}", {"view": view})
    if query is not None and view == "summary":
        return mztab_error("DATASET_BAD_REQUEST",
                           "query は view=samples（因子トークン）か view=tags（タグ名）で使います。",
                           {"view": view})
    if view == "samples":
        return _samples_view(ds, query)
    if view == "tags":
        return _tags_view(ds, query)

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
        "validation": _validation_summary(ds),
        "inchikey_coverage": ds.inchikey_coverage,
        # 群指定（dataset_compare）に必要な試料名と群名。件数だけ返していた頃は、
        # 名前を知る手段が「わざと群サイズ不足のエラーを起こして details を読む」
        # しか無かった。行が並ぶ一覧なので TSV（列名 1 回）で返す。
        "samples": _samples_tsv(ds),
        # 検出状態の有無。「無い」と「全部未検出」は解釈が正反対なので、
        # 率だけを返して available を省くことはしない。
        "detection": _detection_summary(ds),
        # 群（Class ID）・役割ごとの件数、特徴のタグごとの件数、除外の件数。
        # 試料ごと・特徴ごとの中身は view="samples" / view="tags"。
        **sample_group_role_counts(ds),
        "tags": _tag_counts(ds),
        "excluded": {"samples": len(ds.excluded_samples), "features": len(ds.excluded_features)},
    }
    if ds.load_notes:
        # 読み込み時の注記（ジョブの未完了・未検証の警告ほか）。load_dataset の要約を見ていない
        # クライアントが、状態を組み直すときにも気づけるよう件数と先頭だけ出す。
        payload["load_notes"] = {"n": len(ds.load_notes), "first": ds.load_notes[:3]}
    if ds.job_path:
        payload["job_path"] = ds.job_path
    if ds.artifact_paths:
        payload["artifact_roles"] = {role: len(paths) for role, paths in ds.artifact_paths.items()}
    return json_payload(payload)


# ---------- 内部ヘルパ ----------

VIEWS = ("summary", "samples", "tags")

#: `view="tags"` がタグごとに返す特徴 ID の上限。
_TAG_FEATURE_CAP = 50


def _tag_counts(ds) -> dict[str, int]:
    """タグ -> 特徴の件数（`ds.feature_tags`: 特徴 ID -> タグのリスト）。"""
    counts: dict[str, int] = {}
    for tags in (ds.feature_tags or {}).values():
        for tag in tags:
            counts[tag] = counts.get(tag, 0) + 1
    return counts


def _tags_view(ds, query: str | None) -> str:
    wanted = (query or "").strip().casefold()
    by_tag: dict[str, list[str]] = {}
    for feature_id in ds.feature_ids:
        for tag in (ds.feature_tags or {}).get(feature_id, []):
            if not wanted or tag.casefold() == wanted:
                by_tag.setdefault(tag, []).append(feature_id)
    return json_payload({
        "view": "tags",
        "tags": {tag: {"n_features": len(ids), "feature_ids": ids[:_TAG_FEATURE_CAP],
                       "truncated": len(ids) > _TAG_FEATURE_CAP}
                 for tag, ids in by_tag.items()},
    })


def _tsv_cell(value) -> str:
    """TSV の 1 セル。None は空、タブ・改行は空白に潰す（行を壊さない）。"""
    return "" if value is None else str(value).replace("\t", " ").replace("\n", " ")


def _excluded_cell(ds, name: str, row: dict) -> str:
    """除外の理由。`dataset_exclude` の除外はその理由、実験情報シートの include=false はその旨。"""
    reason = ds.excluded_samples.get(name)
    if reason is not None:
        return _tsv_cell(reason)
    return "include=false（実験情報シート）" if not row.get("include", True) else ""


def _samples_view(ds, query: str | None) -> str:
    """試料の TSV。`query` は因子トークン（`msdial/sample_factors.py` の文法。旧 ARF 経路の試料検索と同じ）。"""
    from metabolomix.msdial import sample_factors

    names = list(ds.sample_names)
    rows = _metadata_rows(ds)
    by_name = dict(zip(names, rows))
    facets = sample_factors.build_sample_facets(
        names, None,
        {name: {"role": row.get("role"), "group": row.get("group")}
         for name, row in by_name.items()})

    specs = [part.strip() for part in (query or "").split(",") if part.strip()]
    selected = names
    if specs:
        try:
            # 役割では落とさない（include_roles=None）。QC・ブランクは role 列で見分ける。
            matches, _ = sample_factors.expand_sample_specs(specs, facets, include_roles=None)
        except ValueError as exc:
            return mztab_error("DATASET_BAD_REQUEST", str(exc), {
                "query": query,
                "token_vocabulary": sample_factors.token_vocabulary(facets)})
        hit = {name for found in matches.values() for name in found}
        selected = [name for name in names if name in hit]

    lines = ["name\tgroup\trole\tbatch\tinjection_order\tsubject_id\texcluded"]
    for name in selected:
        row = by_name[name]
        subject = (ds.sample_meta.get(name) or {}).get("subject_id", row.get("subject_id"))
        lines.append("\t".join([
            _tsv_cell(name), _tsv_cell(row.get("group")), _tsv_cell(row.get("role")),
            _tsv_cell(row.get("batch")), _tsv_cell(row.get("injection_order")),
            _tsv_cell(subject), _excluded_cell(ds, name, row)]))
    return json_payload({"view": "samples", "query": query, "n_total": len(names),
                         "n_matched": len(selected), "samples": "\n".join(lines)})



def _validation_summary(ds) -> dict:
    """mzTab-M の構造検査の要約。解析フォルダ（`.arf` / `.arf2`）は mzTab-M ではないので当てはまらない。"""
    if ds.source_format == "msdial_folder":
        return {"applicable": False,
                "note": "解析フォルダのデータセットに mzTab-M の構造検査は無い（読み込みの注記は load_notes）"}
    return {
        "ok": ds.validation_result.get("ok"),
        "n_errors": len(ds.validation_result.get("errors", [])),
        "n_warnings": len(ds.validation_result.get("warnings", [])),
    }


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
    その場で適用する。旧 ARF 経路の役割判定と同じ判定を使うので、
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
        f"## load_dataset 完了: {filename}",
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
        f"## load_dataset 完了: {filename}",
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
