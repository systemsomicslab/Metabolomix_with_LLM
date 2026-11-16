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
# 生成物パスの解決は loading.py が正準。既存の呼び出しが
# ここを参照しているため、名前だけ再エクスポートする（実装は 1 か所）。
from metabolomix.mztab.loading import artifact_abs_path as _artifact_abs_path

__all__ = ["dataset_status"]


def load_summary_lines(ds, opened) -> list[str]:
    """`load_dataset` の戻り値の `summary`（特徴・試料の数、定量、同定の内訳、検出状態、読み込みの注記）。

    解析フォルダのアダプタと mzTab-M のアダプタで書き分ける。出所の裏取り・探索専用の印・
    付随アーティファクトもここに載る。
    """
    if opened.selection is not None:
        lines = _folder_summary_text(ds, opened.selection.arf.name).splitlines()
    else:
        filename = Path(ds.source_files and list(ds.source_files)[-1]
                        or opened.path or "").name or "mzTab-M"
        lines = _summary_text(ds, filename).splitlines()
    lines.append(f"- 出所の裏取り: {ds.source_verification}"
                 + ("（探索専用: 書き出し・図・概観の比較の段に探索専用の印が付きます）"
                    if ds.exploratory_only else ""))
    if ds.artifact_paths:
        lines.append("- 付随アーティファクト: "
                     + str({role: len(paths)
                            for role, paths in ds.artifact_paths.items()}))
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
    `groups` と `n_ungrouped` は統計に入る試料（`included_sample_indices` と同じ集合）だけを数え、外す
    試料は外した理由ごとに分ける（どちらも群 → 件数。群が無ければ `(群なし)`。1 件も無ければ鍵ごと出さない）:
    実験情報シートの `include=false`（`pipeline_sample_sheet` の群 `unassigned` など。`include` の読み方は
    統計エンジンの `_include` と同じ）は `excluded_by_sheet`、`dataset_exclude` で外した試料
    （`ds.excluded_samples`）は `excluded_by_dataset_exclude`（シートでも外した試料はシートの側だけに数える）。
    `groups` に並ぶ群名は dataset_compare にそのまま渡せる。`roles` は全試料を数える。
    """
    rows = _metadata_rows(ds)
    names = list(ds.sample_names)
    dropped = getattr(ds, "excluded_samples", None) or {}
    groups: dict[str, int] = {}
    roles: dict[str, int] = {}
    by_sheet: dict[str, int] = {}
    by_exclude: dict[str, int] = {}
    n_ungrouped = 0
    for index, row in enumerate(rows):
        role = row.get("role") or "sample"
        roles[role] = roles.get(role, 0) + 1
        group = row.get("group")
        name = names[index] if index < len(names) else None
        if not bool(row.get("include", True)):
            by_sheet[group or "(群なし)"] = by_sheet.get(group or "(群なし)", 0) + 1
        elif name is not None and name in dropped:
            by_exclude[group or "(群なし)"] = by_exclude.get(group or "(群なし)", 0) + 1
        elif group:
            groups[group] = groups.get(group, 0) + 1
        else:
            n_ungrouped += 1
    out = {"groups": groups, "roles": roles, "n_ungrouped": n_ungrouped}
    if by_sheet:
        out["excluded_by_sheet"] = by_sheet
    if by_exclude:
        out["excluded_by_dataset_exclude"] = by_exclude
    return out


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_status(dataset_handle: str | None = None, view: str = "summary",
                   query: str | None = None) -> str:
    """現在の DatasetState の状態を返す。`view` で何を見るかを選ぶ。データセットが無ければ missing_state（先に load_dataset）。

    view:
      - `"summary"`（既定）: 概要。`samples` に name/role の TSV。群（Class ID）ごとの件数 `groups`（dataset_compare の
        reference_group / test_group にそのまま使う群名）・役割ごとの件数 `roles`・群の無い試料の数 `n_ungrouped`・特徴のタグごとの件数
        `tags`・除外の件数 `excluded`。`groups` と `n_ungrouped` は統計に入る試料だけを数え、外した試料は群 → 件数で
        `excluded_by_sheet`（シートの include=false）と `excluded_by_dataset_exclude` に分ける。
      - `"samples"`: 試料の TSV（name・group・role・batch・batch_source・injection_order・injection_order_source・subject_id・
        excluded）。`*_source` は値の出所（sheet / mddata / mzTab / raw_acquisition_time / filename。後の 2 つはシートの宣言
        `# injection_order_source = ...`。シートの空欄を .mddata・mzTab-M で埋めた値もその出所）。`query` で絞れる。
      - `"tags"`: タグごとの特徴 ID（各タグ先頭 50 件まで。`n_features` は全件数）。`query` にタグ名を渡すとそのタグだけ。
      - `"results"`: このハンドルの統計の行列 `matrices` と結果 `results` の TSV 2 本（`current` は今の除外・試料の表で作ったものか。
        古い行列 ID はツールが missing_state で断る）。掘るときは `result_id`（と `contrast_id`）を渡す。
      各欄の意味は `docs_read(topic="mztab")` の 13.5b。
    query: `view="samples"` では因子トークン。`_` 区切りで全トークンを含む試料に一致（順不同の AND）し、`,` で区切った複数の指定は
        OR（例 `"ILG_6h"`、`"ILG_6h,control_6h"`）。トークンは試料名・群（Class ID）から取る（日付や拡張子の処理タイムスタンプは除く）。
        QC・ブランクも落とさない。1 件も一致しなければ使えるトークンの語彙を添えて拒む。`view="summary"` / `"results"` では効かない
        ので、使わずに `notes` に書く（断らない）。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
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
            "pipeline の完了後は load_dataset(pipeline_path=...) を推奨します。",
        )
    if view not in VIEWS:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"view は {list(VIEWS)} のどれかです: {view!r}", {"view": view})
    # query は samples（因子トークン）と tags（タグ名）でだけ効く。ほかの view では使わずに notes に書く
    notes = ([f"query は view=samples（因子トークン）か view=tags（タグ名）でだけ効きます。view={view} では"
              "使わずに返しました。"] if query is not None and view in ("summary", "results") else [])
    if view == "samples":
        return _samples_view(ds, query)
    if view == "tags":
        return _tags_view(ds, query)
    if view == "results":
        return _results_view(ds, notes)

    from metabolomix.core import version
    update = version.update_status()
    payload = {
        # 更新後に MCP サーバを再起動し忘れると古いプロセスが黙って動き続ける。
        # 「おかしい」と思ったときに 1 回で分かるよう、ここに刻む。
        "server_version": version.server_version(),
        # 最新なら鍵ごと出さない（戻り値は LLM の文脈をそのまま食う）。
        **({"update_available": update} if update else {}),
        "source_format": ds.source_format,
        # 解釈規則の種別（ハンドルごと。`record_objective(assay_kind=...)` で確定）。
        "assay_kind": getattr(holder, "assay_kind", None) or "unknown",
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
    if notes:
        payload["notes"] = notes
    return json_payload(payload)


# ---------- 内部ヘルパ ----------

VIEWS = ("summary", "samples", "tags", "results")

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


def _source_cell(row: dict, field: str) -> str:
    """値の出所（sheet / mddata / mzTab / raw_acquisition_time / filename）。値が無いバッチは試料名の 8 桁の候補
    （`candidate:<値>`）。"""
    from metabolomix.analysis import sample_manifest

    prov = (row.get("provenance") or {}).get(field) or {}
    if row.get(field) is None:
        candidate = prov.get("candidate") or {}
        return _tsv_cell(f"candidate:{candidate['value']}") if candidate.get("value") else ""
    return _tsv_cell(sample_manifest.order_provenance_label(prov))


def _samples_view(ds, query: str | None) -> str:
    """試料の TSV。`query` は因子トークン（`msdial/sample_factors.py` の文法。旧 ARF 経路の試料検索と同じ）。"""
    from metabolomix.analysis import sample_manifest
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

    # 足した因子の列（sample_manifest の因子の列。名前の順）を後ろに並べる。無ければ列を足さない。
    factor_columns = sorted(sample_manifest._factor_columns_of(rows))
    lines = ["\t".join(["name", "group", "role", "batch", "batch_source", "injection_order",
                       "injection_order_source", "subject_id", "excluded",
                       *map(_tsv_cell, factor_columns)])]
    for name in selected:
        row = by_name[name]
        subject = (ds.sample_meta.get(name) or {}).get("subject_id", row.get("subject_id"))
        lines.append("\t".join([
            _tsv_cell(name), _tsv_cell(row.get("group")), _tsv_cell(row.get("role")),
            _tsv_cell(row.get("batch")), _source_cell(row, "batch"),
            _tsv_cell(row.get("injection_order")), _source_cell(row, "injection_order"),
            _tsv_cell(subject), _excluded_cell(ds, name, row),
            *(_tsv_cell(row.get(column)) for column in factor_columns)]))
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

    役割は試料の表（`_metadata_rows`。`view="samples"`・`sample_counts`・統計エンジンと同じ正準）から
    読む。試料の表で役割を直した（`dataset_set_sample_metadata`）なら、その役割がここにも出る。
    """
    names = list(ds.sample_names)
    rows = _metadata_rows(ds)
    lines = ["name\trole"]
    lines += [f"{_tsv_cell(n)}\t{_tsv_cell(row.get('role') or 'sample')}"
              for n, row in zip(names, rows)]
    return "\n".join(lines)


def _matrix_recipe(prepared) -> tuple[str, str, str]:
    """行列の作り方の要約（ドリフト補正の状態・効いた足切り・正規化）。"""
    steps = {step.get("step"): step for step in (getattr(prepared, "steps", None) or [])}
    filters = (steps.get("filter") or {}).get("filters") or {}
    applied = [name for name, value in filters.items()
               if isinstance(value, dict) and value.get("status") == "applied"]
    return ((steps.get("drift") or {}).get("status") or "",
            "/".join(applied),
            (steps.get("normalize") or {}).get("method") or "")


def _results_view(ds, notes: list[str] | None = None) -> str:
    """`view="results"`: ハンドルの統計の行列と結果の一覧（TSV 2 本）。"""
    from metabolomix.analysis import stats_results
    from metabolomix.analysis.result_state import dataset_fingerprint, is_current, is_stats_current
    from metabolomix.arf2.reader import format_spots_as_table

    fingerprint = dataset_fingerprint(ds)
    revision = getattr(ds, "metadata_revision", 0)
    matrix_rows = []
    for matrix_id, entry in (getattr(ds, "stat_matrices", None) or {}).items():
        prepared = entry.get("prepared")
        drift, filters, normalize = _matrix_recipe(prepared)
        exclusions = entry.get("exclusions") or {}
        matrix_rows.append({
            "matrix_id": matrix_id,
            "current": (entry.get("dataset_hash") == fingerprint
                        and entry.get("metadata_revision") == revision),
            "drift": drift, "filters": filters, "normalize": normalize,
            "n_features": len(getattr(prepared, "feature_ids", None) or []),
            "n_units": len(getattr(prepared, "units", None) or []),
            "excluded_samples": len(exclusions.get("samples") or []),
            "excluded_features": len(exclusions.get("features") or []),
            "group_overrides": len(entry.get("group_overrides") or {}),
        })
    last = {record.get("result_id")
            for record in (getattr(ds, "last_differential", None), getattr(ds, "last_pca", None))
            if isinstance(record, dict)}
    result_rows = []
    for result_id, record in (getattr(ds, "results", None) or {}).items():
        if not isinstance(record, dict):
            continue
        prov = record.get("provenance") or {}
        current = (is_stats_current(ds, record, dataset_hash=fingerprint)
                   if prov.get("dataset_hash") is not None
                   else is_current(ds, record, dataset_hash=fingerprint))
        result_rows.append({
            "result_id": result_id,
            "kind": record.get("kind") or prov.get("kind"),
            "matrix_id": record.get("matrix_id"),
            "contrasts": ",".join(stats_results.contrast_ids(record)),
            "origin": record.get("origin"),
            "current": current,
            "last": result_id in last,
        })
    matrix_columns = ["matrix_id", "current", "drift", "filters", "normalize", "n_features",
                      "n_units", "excluded_samples", "excluded_features", "group_overrides"]
    result_columns = ["result_id", "kind", "matrix_id", "contrasts", "origin", "current", "last"]
    return json_payload({
        "view": "results",
        "n_matrices": len(matrix_rows),
        "matrices": format_spots_as_table(matrix_rows, columns=matrix_columns) if matrix_rows else "",
        "n_results": len(result_rows),
        "results": format_spots_as_table(result_rows, columns=result_columns) if result_rows else "",
        "note": ("current=false の行列・結果は除外・試料の表の更新・データの入れ替えの前に作ったもの"
                 "（作り直す）。結果を掘るときは result_id（と contrast_id）を渡す。"),
        **({"notes": list(notes)} if notes else {}),
    })


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
    # 読み込みの注記（pipeline の run の試料の対応表を結べず当てなかった旨など）。解析フォルダの要約と同じ行。
    if ds.load_notes:
        lines.append(f"- 読み込みの注記 {len(ds.load_notes)} 件: " + "; ".join(ds.load_notes[:3]))
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
