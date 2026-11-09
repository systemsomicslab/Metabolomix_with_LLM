"""最初の概観 `dataset_overview`（spec 2026-11-06-dataset-overview-design、計画 2026-11-09-a-overview）。

品質 → 全体 → クラス → 分子種の 4 段で「群の間で何が変わり、何が変わっていないか」を粗く示す。
`load_dataset` は読み込みの最後に `build_overview_payload(embedded=True)` を呼び、鎖を進めた後は
`dataset_overview` で同じ形を出し直す。統計はすべて C のエンジン（既定の前処理・PCA・spec v2 の比較・組成）を通す。

概観の比較と PCA は `ds.results` に登録する（`result_id` で掘れる）が、直近（`last_differential` / `last_pca`）は
動かさない（判断 #3）。概観そのものは `ds.results[ovw_…]` と `reports/stats/`（step=overview）に残す。
"""
from __future__ import annotations

from mcp.types import ToolAnnotations

from metabolomix.analysis import overview, stats_records, stats_results
from metabolomix.analysis.stats import pca as stats_pca
from metabolomix.analysis.stats import quality, spec_v2
from metabolomix.analysis.stats.composition import composition
from metabolomix.analysis.stats.confounding import confounding_report
from metabolomix.analysis.stats.data import from_dataset
from metabolomix.analysis.stats.preprocess import NeedsInput, normalize_preprocess_config
from metabolomix.core import mcp_core, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload, json_safe, round_floats
from metabolomix.tools._holder import resolve_holder
from metabolomix.tools._matrix_store import register_default

__all__ = ["dataset_overview"]

_TOP_N_MAX = 20
_Q, _EQ_FC, _LOG2FC = 0.05, 1.2, 1.0
_MSDIAL_SOURCE = "msdial_folder"          # quant.table.MSDIAL_SOURCE と同じ（定量の鎖が使える入口）
_NOTE = ("考察は事実 ID（例 {oid}#S1）を挙げて、事実と別の欄に書く。サーバは考察を書かない。"
         "判定は changed / unchanged / inconclusive だけ（有意でないことは変わっていないことではない）。")
_DRIFT_NOTE = ("ドリフト補正の既定が未決（内部標準で換算済みかつ QC あり）のため、概観は補正なしで計算した"
               "（決めるには dataset_preprocess(drift_mode=...)）")


def _curation_judgments(holder, ds) -> int | None:
    """`.arf2` に当たるキュレーションの判断の件数（無い・読めなければ None）。"""
    from metabolomix.curation.apply import flags_for_arf2
    from metabolomix.dataset.artifacts import artifact_path

    path = getattr(holder, "arf2_path", None) or artifact_path(ds, "arf2")
    if not path:
        return None
    try:
        return int(flags_for_arf2(path)["n"])
    except Exception:  # noqa: BLE001 - 判断の件数は状態の表示だけ。読めなくても概観は止めない
        return None


def _recipe(prepared) -> str:
    steps = {s.get("step"): s for s in prepared.steps}
    drift = steps.get("drift", {})
    filters = (steps.get("filter") or {}).get("filters") or {}
    applied = [k for k, v in filters.items() if v.get("status") == "applied"]
    norm = (steps.get("normalize") or {}).get("method")
    return (f"前処理 {prepared.matrix_id}（drift {drift.get('status', '?')}・足切り {'/'.join(applied) or 'なし'}"
            f"・正規化 {norm}・特徴 {len(prepared.feature_ids)}・個体 {len(prepared.units)}）")


def _state_line(holder, ds, chain: dict, prepared, assay_kind: str) -> str:
    source = getattr(holder, "source", None) or {}
    done = [name for name in ("normalize", "merge", "select") if chain[name]]
    cur = chain["curation"]
    return ("・".join([
        f"入口 {source.get('kind') or ds.source_format}",
        f"極性 {getattr(holder, 'polarity', None) or '不明'}",
        f"assay_kind {assay_kind}",
        ("暫定（選別の前: " + ", ".join(chain["pending"]) + " が済んでいない）") if chain["provisional"]
        else "鎖は選別まで済み",
        f"キュレーションの判断 {'不明' if cur is None else cur} 件",
        f"済んだ工程 {', '.join(done) or 'なし'}",
    ]) + "。" + _recipe(prepared))


def _counts(ds, study, prepared) -> dict:
    roles: dict[str, int] = {}
    for s in study.samples:
        if s.include:
            roles[s.role or "unknown"] = roles.get(s.role or "unknown", 0) + 1
    levels = {"msms": 0, "ms1": 0, "none": 0}
    for fid in ds.feature_ids:
        levels[overview.id_level(ds, fid)] += 1
    return {"samples_by_role": roles, "n_features": len(ds.feature_ids), "id_levels": levels}


def _compare(ds, study, prepared, config, inputs, comps):
    spec = spec_v2.normalize_spec_v2(stats_results.compare_spec(
        factors=["group"], comparisons=comps, covariates=[], block=None, method="moderated_t",
        effect_size="log2_mean_difference", bh_population="matrix", q_threshold=_Q,
        equivalence_fc=_EQ_FC, allow_confounded=False, preprocess=config, inputs=inputs))
    run = spec_v2.run_spec_v2(spec, study, prepared=prepared)
    return stats_results.register_compare(ds, run, log2fc_threshold=_LOG2FC, units=prepared.units,
                                          make_current=False, origin="overview")


def _write_record(holder, entry, oid, matrix_id, settings, summary) -> dict:
    try:
        reports_dir = mcp_core._resolve_report_dir(holder)
        written = stats_records.write_stats_record(
            reports_dir, "overview",
            inputs={"dataset_hash": entry.get("dataset_hash"), "matrix_id": matrix_id, "result_id": oid},
            parent_ids=[entry["record_id"]] if entry.get("record_id") else [],
            settings=settings, summary=summary)
    except DomainError as exc:
        return {"written": False, "reason": f"{exc.code}: {exc.message}"}
    except OSError as exc:
        return {"written": False, "reason": f"OSError: {exc}"}
    return {"written": True, "record_id": written["record_id"], "path": written["path"]}


def _texts(section: dict) -> dict:
    """戻り値の段: 事実は 1 文字列ずつ（構造は記録と ds.results に残す）。"""
    out = {k: v for k, v in section.items() if k != "facts" and v not in (None, "")}
    out["facts"] = [x["text"] for x in section.get("facts", [])]
    return out


def build_overview_payload(holder, *, top_n: int = 5, embedded: bool = False) -> dict:
    """概観を計算して登録し、戻り値の dict を返す。`DomainError` はそのまま送出する（呼び出し側が封筒にする）。

    `embedded=True`（`load_dataset` の中）は `counts` を省く（`load_dataset` の summary と重なるため）。
    """
    ds = holder.dataset
    assay_kind = getattr(session_state.session, "assay_kind", "unknown") or "unknown"

    # 1. 既定の前処理（ドリフト補正の既定が未決なら補正なしで計算する。§7）
    _, inputs, config = stats_results.matrix_inputs(ds, None)
    study = from_dataset(ds, **inputs)
    prepared = stats_results.prepare_default(study, config)
    drift_note = None
    if isinstance(prepared, NeedsInput):
        config = normalize_preprocess_config({**config, "drift": {**config["drift"], "mode": "off"}})
        prepared = stats_results.prepare_default(study, config)
        drift_note = _DRIFT_NOTE
    entry, _record = register_default(holder, ds, study, prepared, config, inputs)

    # 2. 品質
    plan = overview.group_plan(ds, prepared)
    conf = confounding_report(prepared.units, plan["levels"])
    rsd = quality.qc_rsd_summary(study)
    failed = quality.failed_qc(study)
    drift = quality.drift_signs(study)
    replicates = [w for w in prepared.warnings if w.get("code") == "possible_technical_replicates"]
    verdict, reasons = quality.quality_verdict(rsd=rsd, failed=failed, drift=drift, confounding=conf,
                                               replicate_warnings=replicates)
    n_blank = sum(1 for s in study.samples if s.include and s.role == "blank")

    # 3. 全体（PCA は直近にしない）
    pca_rec, pca_result = None, None
    labels = stats_results.unit_group_labels(ds, prepared)
    try:
        pca_result = stats_pca.pca(prepared)
        pca_rec = stats_results.register_pca(ds, pca_result, labels, make_current=False)
    except DomainError:
        pca_result = None

    # 4. 比較（群が 2 つ以上・探索専用でない・完全交絡でないとき。直近にしない）
    comps, pairs_omitted = overview.group_comparisons(plan)
    record, omitted = None, None
    if plan["mode"] == "none":
        omitted = plan["reason"]
    elif getattr(ds, "exploratory_only", False):
        omitted = "探索専用のデータセット（allow_incomplete=True で開いた未完了の実行）なので比較しない"
    elif conf.get("batch_status") == "complete":
        omitted = "群とバッチが完全交絡なので比較しない（品質の段を見る）"
    else:
        try:
            record = _compare(ds, study, prepared, config, inputs, comps)
        except DomainError as exc:
            omitted = f"{exc.code}: {exc.message}"

    q_note = None if verdict != "poor" else "品質 Q1 が poor のため参考"
    q_sec = overview.quality_facts(rsd=rsd, failed=failed, detection=quality.detection_summary(study),
                                   drift=drift, confounding=conf, replicate_warnings=replicates,
                                   verdict=verdict, reasons=reasons, n_qc=len(prepared.qc_ids),
                                   n_blank=n_blank, drift_note=drift_note,
                                   pca_id=pca_rec["result_id"] if pca_rec else None)
    g_sec = overview.global_facts(pca_result=pca_result, labels=labels, record=record, plan=plan,
                                  omitted=omitted, quality_note=q_note)

    # 5. クラス（付加体の選別の後だけ。§4.3）
    chain = overview.chain_state(ds, curation_judgments=_curation_judgments(holder, ds))
    if getattr(ds, "source_format", None) != _MSDIAL_SOURCE:
        gate = "解析フォルダ以外（mzTab-M など）から開いたので定量の鎖（付加体の選別）が使えない"
        gate_next = 'dataset_analysis_run（analysis-spec.v2 の kind="composition"）で組成を明示に計算する'
    elif not chain["select"]:
        gate = "付加体の選別の前（内部標準の行と同じ分子種の複数の付加体がクラスの合計に入り、見せかけの unchanged が出る）"
        gate_next = ("dataset_quant_preview → dataset_quant_apply（換算しないなら skip=True）→ "
                     "dataset_select_features → dataset_overview()")
    elif record is None:
        gate, gate_next = f"比較を回していない（{omitted}）", g_sec["next"]
    else:
        gate, gate_next = None, ""
    comp_result = None
    if gate is None:
        try:
            comp_result = composition(prepared, assay_kind=assay_kind, factors=["group"], comparisons=comps,
                                      test="moderated_t", effect_size="log2_mean_difference",
                                      bh_population="matrix", q_threshold=_Q, equivalence_fc=_EQ_FC)
        except DomainError as exc:
            gate, gate_next = f"{exc.code}: {exc.message}", g_sec["next"]
    c_sec = overview.class_facts(composition_result=comp_result, gate=gate, gate_next=gate_next,
                                 plan=plan, quality_note=q_note)

    # 6. 分子種
    s_sec = overview.species_facts(ds, record=record, plan=plan, top_n=top_n, pairs_omitted=pairs_omitted,
                                   omitted=omitted, quality_note=q_note)

    sections = {"quality": q_sec, "global": g_sec, "classes": c_sec, "species": s_sec}
    ids = {"compare": record["result_id"] if record else None,
           "pca": pca_rec["result_id"] if pca_rec else None}
    oid = overview.overview_id(dataset_hash=prepared.dataset_hash,
                               metadata_revision=getattr(ds, "metadata_revision", 0),
                               matrix_id=prepared.matrix_id, compare_id=ids["compare"], pca_id=ids["pca"],
                               extra={"top_n": top_n, "drift_off_undecided": drift_note is not None,
                                      "assay_kind": assay_kind, "chain": chain})
    settings = {"thresholds": quality.QUALITY_THRESHOLDS, "top_n": top_n, "plan": plan,
                "verdict_rule": {"q": _Q, "equivalence_fc": _EQ_FC}, "log2fc_threshold": _LOG2FC,
                "result_ids": ids, "chain": chain, "assay_kind": assay_kind}
    reg = overview.register_overview(ds, result_id=oid, matrix_id=prepared.matrix_id, sections=sections,
                                     parents=[prepared.matrix_id, ids["compare"], ids["pca"]],
                                     settings=json_safe(settings))
    written = _write_record(holder, entry, oid, prepared.matrix_id, json_safe(settings),
                            json_safe(round_floats(reg["facts"])))
    payload = {
        "status": "success",
        "overview_id": oid,
        "provisional": chain["provisional"],
        "state": _state_line(holder, ds, chain, prepared, assay_kind),
        **({} if embedded else {"counts": _counts(ds, study, prepared)}),
        "quality": _texts(q_sec),
        "global": _texts(g_sec),
        "classes": _texts(c_sec),
        "species": _texts(s_sec),
        "result_ids": {k: v for k, v in ids.items() if v},
        "record": written,
        "note": _NOTE.format(oid=oid),
    }
    return json_safe(payload)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False),
          structured_output=False)
def dataset_overview(top_n: int = 5, dataset_handle: str | None = None) -> str:
    """最初の概観: 品質 → 全体 → クラス → 分子種の 4 段で「群の間で何が変わり、何が変わっていないか」を粗く示す。

    `load_dataset` が読み込みの最後に同じものを返す（`overview` 欄）。鎖を進めた後（キュレーション・換算・統合・
    選別・除外・試料の表の更新）に呼び直すと、同じ形で出し直す。図は入れない（各段の `next` の dataset_plot で描く）。

    - 計算は C の既定の前処理（ドリフト補正の既定が未決なら補正なしで計算し、そう書く）・PCA（QC は投影）・
      moderated t（2 群は 1 対比、3〜6 群は全対と moderated F、7 群以上は基準との対比と F）。判定は
      changed / unchanged / inconclusive だけ（q≤0.05、同等 ±1.2 倍）。
    - 2 群の向きは**仮**（試料の表で先に現れた群を基準。名前から対照を推測しない）。向きを決めるのは
      `dataset_compare(reference_group=..., test_group=...)`。
    - 品質の判定（`ok` / `caution` / `poor` / `not_evaluable`）を先に置く。`poor` でも差の段は「参考」と添えて出すが、
      群とバッチが完全交絡なら比較を回さない。群が 1 つ・探索専用のデータセットも比較しない。
    - クラスの段は付加体の選別（`dataset_select_features`）の後だけ出す（それまでは理由と次の一手）。
    - 事実には段ごとの ID（`Q1` `G1` `C1` `S1`）が付き、外からは `<overview_id>#Q1` で引く。考察は事実 ID を挙げて
      事実と別に書く（サーバは考察を書かない）。
    - 比較と PCA は `result_ids` で掘れる（volcano・特徴検索・書き出し）が、直近の結果にはしない。
      概観は `reports/stats/` に記録を 1 件追記する（書けなくても概観は成功し `record.written=false`）。

    top_n: 分子種の段の上位の件数（既定 5、0〜20）。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。
    """
    if isinstance(top_n, bool) or not isinstance(top_n, int) or not 0 <= top_n <= _TOP_N_MAX:
        return mztab_error("DATASET_BAD_REQUEST", f"top_n は 0〜{_TOP_N_MAX} の整数です。", {"top_n": top_n})
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    if holder.dataset is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    try:
        return json_payload(build_overview_payload(holder, top_n=top_n))
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, json_safe(exc.details) if exc.details else None)
