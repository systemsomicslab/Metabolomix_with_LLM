"""再解析の spec で呼ぶ MCP ツール（spec 2026-10-22 §4）。

**再解析で使う**（判定の数字を、その場で書くコードではなく版のあるテスト済みの実装から
出すため）。数値の実装は `metabolomix/analysis/` の純ロジックにあり、ここは spec を
正規化して呼び、戻り値を小さくするだけ。どのツールも `dataset_handle` を受け、省略時は
対話のセッションのデータセットを読む（書き換えない）。

大きな結果は `output_path` に JSON（キーの順を固定・改行なし）で書き、要約だけを返す。
同じ spec・同じ入力なら同じファイルになる。
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Literal

from mcp.types import ToolAnnotations

from metabolomix.core.atomic_io import DomainError, atomic_write_bytes
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload, round_floats
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_alignment_qc", "dataset_analysis_run", "dataset_claim_evidence",
           "dataset_feature_relations", "dataset_homolog_check", "dataset_mass_accuracy",
           "dataset_novel_candidates", "library_msms_qualify"]

#: output_path なしで行を返す上限（解析ごと）。超えれば output_path を求める。
MAX_INLINE_ROWS = 50

ROW_COLUMNS = ("feature_id", "status", "reason", "group_ref", "group_test", "n_ref", "n_test",
               "mean_ref", "mean_test", "fc_test_over_ref", "log2fc_test_over_ref",
               "statistic", "p", "q", "beta_log2", "direction", "meets_thresholds",
               "merged_from")
_SIGNIFICANT = frozenset({"mean_ref", "mean_test", "statistic", "p", "q", "beta_log2"})
_DECIMAL4 = frozenset({"fc_test_over_ref", "log2fc_test_over_ref"})


def _cell(column: str, value):
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return ";".join(str(v) for v in value)
    if isinstance(value, float) and math.isfinite(value):
        if column in _SIGNIFICANT:
            return str(float(f"{value:.4g}"))
        if column in _DECIMAL4:
            rounded = round(value, 4)
            return str(0.0 if rounded == 0 else rounded)
        # そのほかの float は有効 6 桁（既定 repr の 17 桁を返さない）
        return str(float(f"{value:.6g}"))
    return str(value)


def rows_table(rows: list[dict], columns=ROW_COLUMNS) -> str:
    from metabolomix.arf2.reader import format_spots_as_table

    cells = [{c: _cell(c, row.get(c)) for c in columns} for row in rows]
    return format_spots_as_table(cells, columns=list(columns)) if cells else ""


def _finite(value):
    """JSON に書けない非有限の float を null にする（集合は並べた list にする）。"""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {k: _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(_finite(v) for v in value)
    return value


def write_result(output_path: str, payload: dict) -> str:
    """結果の JSON を書く（キーの順を固定・改行なし・非有限は null）。書いた絶対パスを返す。"""
    path = Path(output_path)
    data = json.dumps(_finite(payload), ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)
    atomic_write_bytes(path, data.encode("utf-8"))
    return str(path.resolve())


def _dataset_or_error(holder, needed: bool):
    """(ds, error)。データセットが要るのに無ければ missing_state の封筒。"""
    ds = getattr(holder, "dataset", None)
    if needed and ds is None:
        return None, missing_state(
            "dataset", ["dataset_load", "dataset_handle_open"],
            "データセットがありません。dataset_handle_open で開いたハンドルを渡すか、"
            "dataset_load を実行してください。")
    return ds, None


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=True),
          structured_output=False)
def dataset_analysis_run(spec: dict, dataset_handle: str | None = None,
                         output_path: str | None = None) -> str:
    """論文の方法を写した spec（analysis-spec.v1）で、行列を作って 2 群を検定する（再解析で使う）。

    spec の形と欄の意味は `lipidmix://docs/output-format/reanalysis` §2。要点:
    `comparison`（`reference` / `test` の群の名前。比較名は `<test>_vs_<ref>`）、
    `analyses`（名前つきの解析を並べる——感度解析は 1 つの spec に並べる。各解析は
    `matrix`（matrix-recipe.v2: 段の順序・検出率・KNN・PQN・LOESS・RSD）と
    `test`（`student` / `welch` / `mannwhitney` / `ols`）・`fold_change`・`correction`・
    `bh_population`・`thresholds`）、入力は `input`（データセットか外部の行列 TSV）、
    群は `groups`（実験情報シートか TSV）、`covariates`（OLS の共変量の TSV）、
    `aggregate`（クラスの和・2 特徴の比を 1 本の特徴として検定）、`features`（報告する特徴と、
    `merge_split`——割れたスポットを試料ごとに大きい方の値の 1 列にまとめてから検定する。既定 true、
    行の `merged_from`、外部の行列では何もしない）。`bh_population` は `matrix` / `reported` /
    `annotated`（SME のある特徴だけ）。

    返すのは数字と向き（`direction`）と閾値を機械的に当てた `meets_thresholds` だけで、
    再現したかの判定は返さない。結果には正規化した spec・`rule_hash`・`input_fingerprint`
    が入り、同じ spec・同じ入力なら同じ JSON になる。

    output_path: 結果の JSON をここへ書き、要約だけを返す。省略時は解析ごとの行が
        50 以下のときだけ行を TSV で返す（超えれば `ANALYSIS_OUTPUT_REQUIRED`）。
    対話のセッション・ハンドルの状態は変えない（読むだけ）。
    """
    from metabolomix.analysis import analysis_spec

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    try:
        normalized = analysis_spec.normalize_spec(spec)
        ds, error = _dataset_or_error(holder, normalized["input"]["source"] == "dataset"
                                      or normalized["aggregate"] is not None
                                      and "select" in (normalized["aggregate"] or {}))
        if error:
            return error
        result = analysis_spec.run_spec(normalized, ds)
        too_many = [a["analysis_id"] for a in result["analyses"]
                    if len(a["rows"]) > MAX_INLINE_ROWS]
        if output_path is None and too_many:
            raise DomainError(
                "ANALYSIS_OUTPUT_REQUIRED",
                f"報告する行が {MAX_INLINE_ROWS} を超える解析がある（{too_many}）。output_path を"
                "指定して JSON に書くか、features.feature_ids で絞ってください。",
                {"analyses": too_many})
        written = write_result(output_path, result) if output_path else None
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    summary = []
    for analysis in result["analyses"]:
        rows = analysis["rows"]
        item = {
            "analysis_id": analysis["analysis_id"],
            "n_features_in": analysis["matrix"]["n_features_in"],
            "n_features_out": analysis["matrix"]["n_features_out"],
            "bh_population_n": analysis["bh_population_n"],
            "n_rows": len(rows),
            "n_tested": sum(1 for r in rows if r["status"] == "tested"),
            "n_meets_thresholds": sum(1 for r in rows if r["meets_thresholds"]),
        }
        if "aggregate" in analysis:
            item["aggregate"] = analysis["aggregate"]
        if "covariate_excluded_samples" in analysis:
            item["covariate_excluded_samples"] = analysis["covariate_excluded_samples"]
        if written is None:
            item["rows"] = rows_table(rows)
        summary.append(item)
    return json_payload({
        "status": "success",
        "schema": result["schema"],
        "rule_hash": result["rule_hash"],
        "merge_split": {"enabled": result["merge_split"]["enabled"],
                        "source": result["merge_split"]["source"],
                        "n_groups_used": len(result["merge_split"]["groups"]),
                        "n_groups_skipped": len(result["merge_split"]["skipped"]),
                        "groups_hash": result["merge_split"]["groups_hash"]},
        "input_fingerprint": result["input_fingerprint"],
        "comparison": result["comparison"],
        "group_ref": result["group_ref"],
        "group_test": result["group_test"],
        **({"output_path": written} if written else {}),
        "analyses": summary,
        "note": "判定（再現したか）は返さない。direction は group_test から見た向き。",
    })


# ---------- 共通 ----------

def _polarity(holder, ds) -> str | None:
    """ハンドルなら開いたときに決めた極性、対話のセッションならその場で決める。"""
    if hasattr(holder, "handle_id"):
        return holder.polarity
    from metabolomix.mztab.polarity import dataset_polarity
    return dataset_polarity(ds)[0]


def _spot_query_for(holder, dataset_handle):
    """ハンドルに `.arf2` / `.dcl` が無ければ、対話のフォルダを探さずに断る spot の照会。"""
    from metabolomix.library import tools as library_tools
    from metabolomix.tools._holder import alignment_paths

    arf2_path, dcl_path, lacking = alignment_paths(holder, dataset_handle, None, None)
    if lacking:
        return (lambda *_a, **_k: dict(lacking)), arf2_path, dcl_path
    return library_tools.spot_msms_query, arf2_path, dcl_path


def _library_sha(store) -> str | None:
    return store.summary()["source_sha256"] if store is not None else None


_QUALIFY_COLUMNS = ("rank", "name", "formula", "inchikey", "adduct", "total_score",
                    "matched_peaks_count", "fragment_matches_excl_precursor")


# ---------- MS/MS の合否 ----------

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def library_msms_qualify(rules: dict, spot_id: int | None = None,
                         precursor_mz: float | None = None, rt: float | None = None,
                         ion_mode: str | None = None, arf2_path: str | None = None,
                         dcl_file: str | None = None, top_n: int = 10,
                         dataset_handle: str | None = None) -> str:
    """MS/MS の照合に、spec で書いた合否の規則（msms-rules.v1）を当てる（再解析で使う）。

    照合は `library_match_feature` と同じ（`spot_id` か `precursor_mz` のどちらか一方）。
    規則の形と検査の意味は `lipidmix://docs/output-format/reanalysis` §4（`target` は
    `inchikey` か、InChIKey が一意に決まらない脂質の種の `lipid_species`）。設定した検査だけを
    `checks` に並べ、`verdict` は `qualified` / `not_qualified` / `not_evaluable`。
    `rule_hash` は正規化した規則から作る。MS/MS が無ければ `not_evaluable`（「断片が
    合わなかった」ではない）。「明らかに別物か」は判断しない（再解析の規則で決める）。
    dataset_handle: ハンドルのライブラリ・`.arf2` / `.dcl` を使う。省略時は対話のセッション。
    """
    from metabolomix.analysis import msms_rules
    from metabolomix.analysis.spectral_match import fragment_matches_excluding_precursor
    from metabolomix.library import matching
    from metabolomix.library import tools as library_tools

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    store = holder.library.store
    if store is None:
        return missing_state("library_store", ["library_load", "dataset_handle_open"],
                             "参照ライブラリが読み込まれていません。")
    if (spot_id is None) == (precursor_mz is None):
        return mztab_error("DATASET_BAD_REQUEST", "spot_id か precursor_mz のどちらか一方を指定する",
                           {"field": "spot_id"})
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        return mztab_error("DATASET_BAD_REQUEST", "top_n は 1 以上の整数", {"field": "top_n"})
    try:
        normalized = msms_rules.normalize_rules(rules)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    from metabolomix.tools._holder import alignment_paths

    arf2_path, dcl_file, lacking = alignment_paths(
        holder, dataset_handle, arf2_path, dcl_file,
        need=("arf2", "dcl") if spot_id is not None else ("dcl",))
    if lacking:
        return json_payload({**lacking, "verdict": "not_evaluable", "checks": []})
    if spot_id is not None:
        spot = library_tools.spot_msms_query(int(spot_id), arf2_path=arf2_path, dcl_file=dcl_file,
                                             rt=rt, ion_mode=ion_mode)
        measured = spot.get("measured") if spot["status"] == "ok" else None
        if spot["status"] == "ok":
            precursor_mz, rt, ion_mode = spot["precursor_mz"], spot["rt"], spot["ion_mode"]
    else:
        measured = library_tools._measured_spectrum(precursor_mz, rt=rt, dcl_file=dcl_file)
        spot = {"status": "ok" if measured else "not_found"}
    base = {"rule_hash": msms_rules.rules_hash(normalized), "rules": normalized,
            "query": {"spot_id": spot_id, "precursor_mz": precursor_mz, "rt": rt,
                      "ion_mode": ion_mode}}
    if not measured:
        return json_payload(round_floats({
            "status": spot["status"], "verdict": "not_evaluable", "checks": [], **base,
            "message": spot.get("message") or "MS/MS が無い（未取得）。断片が合わなかったのではない。"}))
    matched = matching.score_candidates(store, precursor_mz=precursor_mz, rt=rt,
                                        ion_mode=ion_mode, measured=measured)
    ms2_tol = matched["query"]["ms2_tol"]
    if matched["status"] == "no_candidates":
        return json_payload(round_floats({"status": "no_candidates", "verdict": "not_evaluable",
                                          "checks": [], **base}))
    outcome = msms_rules.qualify(matched["scored"], measured, precursor_mz, normalized,
                                 ms2_tol=ms2_tol)
    tol = normalized["precursor_exclusion_tol"] or ms2_tol
    rows = [{"rank": rank, "name": record.get("name"), "formula": record.get("formula"),
             "inchikey": record.get("inchikey"), "adduct": record.get("adduct"),
             "total_score": result["total_score"],
             "matched_peaks_count": result["matched_peaks_count"],
             "fragment_matches_excl_precursor": fragment_matches_excluding_precursor(
                 result["alignment"], precursor_mz, tol)}
            for rank, (record, result) in enumerate(matched["scored"][:top_n], start=1)]
    return json_payload(round_floats({
        "status": "success", **outcome, **base,
        "n_candidates": matched["n_candidates"],
        "candidates": rows_table(rows, _QUALIFY_COLUMNS),
    }))


# ---------- 同族列 ----------

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_homolog_check(spec: dict, dataset_handle: str | None = None) -> str:
    """同族列の RT の単調性を検査する（homolog-check.v1。再解析で使う）。

    `members`（`feature_id`・`carbons`・`double_bonds` を手で並べる）か `auto`
    （`lipid_class` の SME rank 1 の名前を GOSLIN で種のレベルに読んで集める。`evidence`
    に `sml` を足すと SML の名前も見る。`adduct` で SME の付加体を絞る）のどちらか一方と、
    `rule`（`carbon_direction` 既定 increasing・`double_bond_direction` 既定 decreasing・
    `max_violations` 既定 0）。同じ（炭素数, 二重結合数）が複数あれば RT の中央値で並べ
    `duplicates` に挙げる。戻り値は `verdict`（pass / fail）・違反の TSV・`rule_hash`。
    dataset_handle: ハンドルのデータセット。省略時は対話のセッション。
    """
    from metabolomix.analysis import homolog

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds, error = _dataset_or_error(holder, True)
    if error:
        return error
    try:
        normalized = homolog.normalize_spec(spec)
        if normalized["members"] is not None:
            members = normalized["members"]
            known = set(ds.feature_ids)
            unknown = [m["feature_id"] for m in members if m["feature_id"] not in known]
            if unknown:
                raise DomainError("FEATURE_NOT_FOUND", f"特徴が無い: {unknown}",
                                  {"feature_ids": unknown})
        else:
            features = []
            for fid in ds.feature_ids:
                candidates = (ds.feature_candidates or {}).get(fid) or []
                top = candidates[0] if candidates else {}
                annotation = (ds.feature_annotations or {}).get(fid) or {}
                features.append({"feature_id": fid, "adduct": top.get("adduct"),
                                 "names": {"sme": top.get("chemical_name"),
                                           "sml": annotation.get("name")}})
            members = homolog.members_from_names(features, normalized["auto"])
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    rt_of = {fid: (ds.feature_metadata.get(fid) or {}).get("rt") for fid in ds.feature_ids}
    outcome = homolog.check(members, rt_of, normalized["rule"])
    violations = [{"axis": v["axis"], "low": v["low"], "high": v["high"], "rt_low": v["rt_low"],
                   "rt_high": v["rt_high"], "features_low": ";".join(v["features_low"]),
                   "features_high": ";".join(v["features_high"])} for v in outcome["violations"]]
    return json_payload(round_floats({
        "status": "success",
        "verdict": outcome["verdict"],
        "n_members": len(outcome["members"]),
        "n_violations": outcome["n_violations"],
        "violations": rows_table(violations, ("axis", "low", "high", "rt_low", "rt_high",
                                              "features_low", "features_high")),
        "members": rows_table(outcome["members"], ("feature_id", "carbons", "double_bonds", "rt",
                                                   "status")),
        "duplicates": outcome["duplicates"],
        "rule": normalized["rule"],
        "rule_hash": homolog.rule_hash(normalized),
    }, digits=4))


# ---------- 特徴どうしの関係 ----------

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_feature_relations(feature_id: str, rt_window: float = 0.05, mz_tol: float = 0.01,
                              min_r: float = 0.7, adducts: list[str] | None = None,
                              dataset_handle: str | None = None) -> str:
    """mzTab-M の特徴 1 つについて、同時溶出の特徴と付加体・同位体・源内断片の関係を返す
    （再解析で使う。関係の規則はアラインメントのキュレーションの関係と同じ）。

    rt_window: 同時溶出とみなす RT の差（分、既定 0.05）。mz_tol: 関係の m/z の許容（Da）。
    min_r: 試料間の強度の相関（log1p の Pearson）がこれ以上なら裏付けあり。
    adducts: 付加体の候補（省略時はデータセットの極性の既定の層を平らにしたもの）。
    `as_ion_of`（この特徴が相手のイオン）と `ions_of_this`（相手がこの特徴のイオン）を TSV で
    返す。相手の付加体は SME rank 1 の付加体で、同定の無い相手には付加体の関係は出ない
    （同位体・源内断片だけ）。MS-DIAL のリンクは使わない（裏付けは相関だけ）。
    dataset_handle: ハンドルのデータセット。省略時は対話のセッション。
    """
    from metabolomix.analysis import claim_evidence, feature_search

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds, error = _dataset_or_error(holder, True)
    if error:
        return error
    polarity = _polarity(holder, ds)
    if adducts is None:
        adducts = [a for tier in feature_search.DEFAULT_ADDUCT_TIERS.get(polarity or "", [])
                   for a in tier]
    try:
        out = claim_evidence.feature_relations(ds, str(feature_id), rt_window=rt_window,
                                               mz_tol=mz_tol, min_r=min_r, adducts=adducts)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    relation_columns = ("of", "of_name", "relation", "dmz_mda", "drt", "r", "isotope_ratio",
                        "expected_ratio", "isotope_consistent", "strong")
    ion_columns = ("ion", "relation", "dmz_mda", "drt", "r", "isotope_consistent", "strong")
    return json_payload(round_floats({
        "status": "success",
        "feature_id": out["feature_id"],
        "query": {"rt_window": rt_window, "mz_tol": mz_tol, "min_r": min_r, "adducts": adducts},
        "n_coeluting": len(out["coeluting"]),
        "coeluting": rows_table(out["coeluting"], ("feature_id", "mz", "rt", "drt", "r", "name")),
        "as_ion_of": rows_table(out["as_ion_of"], relation_columns),
        "ions_of_this": rows_table(out["ions_of_this"], ion_columns),
    }, digits=4))


# ---------- 証拠の束 ----------

#: output_path なしで候補を返す上限。
MAX_INLINE_CANDIDATES = 10


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=True),
          structured_output=False)
def dataset_claim_evidence(spec: dict, dataset_handle: str | None = None,
                           output_path: str | None = None) -> str:
    """1 主張（目標の物質）の証拠の束（claim-evidence.v1）を作る（再解析で使う）。

    spec: `claim_id`・`target`（`dataset_feature_search` の目標と同じ形。`target_id` は
    `claim_id`。脂質の種の主張は `inchikey` の代わりに `lipid_species`）・`scan`（`ppm` と `rt_tol_min` の水準のリスト、任意で `adduct_tiers`——
    全組で走査する）・`msms`（`top_n` と msms-rules.v1 の `rules`。target は書かない）・
    `relations`（`rt_window` / `mz_tol` / `min_r`）・`max_candidates`（既定 50）・`merge_split`
    （既定 true。割れたスポットの組を判定で 1 つの候補（残す ID）と数え、窓の外の相手も候補に
    `merge_partner` で足す。false は v7 と同じ規則・同じ `rule_hash`）。
    束に入るのは、規則で選んだ既定の特徴（`default_selection`。`dataset_feature_search`
    の目標照合と同じ規則で、決まらなければ `needs_input`）、走査ごとの判定、候補ごとの
    照合（spot_id = SMF_ID）・合否・関係。**特徴は決めない**（実行役が決め、統計は
    `dataset_analysis_run` を別に呼ぶ）。`rule_hash`・`input_fingerprint`（データの指紋・
    ライブラリの sha256）を含み、同じ spec・同じ入力なら同じ JSON になる。
    output_path: JSON をここへ書き、要約だけを返す（候補が 10 を超えるときは必須）。
    dataset_handle: ハンドルのデータセット・ライブラリ・`.arf2` / `.dcl`。省略時は対話のセッション。
    """
    from metabolomix.analysis import claim_evidence
    from metabolomix.library import tools as library_tools

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds, error = _dataset_or_error(holder, True)
    if error:
        return error
    polarity = _polarity(holder, ds)
    spot_query, arf2_path, dcl_path = _spot_query_for(holder, dataset_handle)
    store = holder.library.store
    try:
        normalized = claim_evidence.normalize_claim_spec(spec, dataset_polarity=polarity)
        bundle = claim_evidence.build(ds, normalized, store=store,
                                      spot_query=spot_query,
                                      arf2_path=arf2_path, dcl_path=dcl_path,
                                      dataset_polarity=polarity,
                                      library_sha256=_library_sha(store))
        if output_path is None and len(bundle["candidates"]) > MAX_INLINE_CANDIDATES:
            raise DomainError(
                "ANALYSIS_OUTPUT_REQUIRED",
                f"候補が {MAX_INLINE_CANDIDATES} を超える（{len(bundle['candidates'])} 件）。"
                "output_path を指定して JSON に書いてください。",
                {"n_candidates": len(bundle["candidates"])})
        written = write_result(output_path, bundle) if output_path else None
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    if written is None:
        return json_payload(round_floats({"status": "success", **bundle}))
    return json_payload(round_floats({
        "status": "success", "output_path": written, "claim_id": bundle["claim_id"],
        "rule_hash": bundle["rule_hash"], "input_fingerprint": bundle["input_fingerprint"],
        "default_selection": bundle["default_selection"],
        "n_candidates": len(bundle["candidates"]), "truncated": bundle["truncated"],
        "scans": [{k: s[k] for k in ("ppm", "rt_tol_min", "status", "selected_feature_id")}
                  for s in bundle["scans"]],
    }))


# ---------- 質量精度・新しい知見の候補 ----------

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_mass_accuracy(band_edges: list[float] | None = None, rank: int = 1,
                          subset: Literal["rank1", "msms_matched"] = "rank1",
                          dataset_handle: str | None = None) -> str:
    """SME（MS/MS の照合で付いた同定）の質量誤差（ppm）の分布を m/z 帯ごとに返す（再解析で使う）。

    band_edges: 帯の境界（理論 m/z。既定 [200] → `<200` と `≥200`）。rank: 見る SME 候補の順位
    （既定 1）。帯ごとに `n median q1 q3 within_5ppm within_10ppm`（四分位は
    `statistics.quantiles(n=4)`）。許容幅は決めない（再解析の規則で決める）。
    subset: `rank1`（既定。rank の SME すべて）か `msms_matched`（そのうち、特徴の MS/MS が参照
    スペクトルに照合されたもの。SML の `opt_global_spectrum_matched` が真。LBM 由来も研究室 MSP
    由来も数える。この旗は rank 1 の同定を表すので `rank=1` のときだけ。列の無い mzTab では
    `subset_available: false` と `note` を付けて n=0 を返す——手順の後退（帯の n 不足 → 極性全体 →
    rank1）を止めない）。
    `injection_trend`: 試料ごとの実測 m/z（隣接 `.arf` の検出されたセル）があれば、試料ごとの
    誤差の中央値 `per_injection`（`sample` `injection_order` `median_ppm` `n`）と注入順との
    `spearman_rho`（5 試料未満、または注入順の出所が `filename`＝採取時刻が読めないときは null。
    理由は `rho_unavailable_reason`）・`n_samples`・`injection_order_source`。無ければ
    `{"status": "unavailable"}`。
    dataset_handle: ハンドルのデータセット。省略時は対話のセッション。
    """
    from metabolomix.analysis import mass_accuracy

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds, error = _dataset_or_error(holder, True)
    if error:
        return error
    edges = [200.0] if band_edges is None else band_edges
    if (not isinstance(edges, list) or not all(isinstance(e, (int, float)) and not isinstance(e, bool)
                                                and e > 0 for e in edges)):
        return mztab_error("DATASET_BAD_REQUEST", "band_edges は正の数のリスト",
                           {"field": "band_edges"})
    if isinstance(rank, bool) or not isinstance(rank, int) or rank < 1:
        return mztab_error("DATASET_BAD_REQUEST", "rank は 1 以上の整数", {"field": "rank"})
    if subset not in mass_accuracy.SUBSETS:
        return mztab_error("DATASET_BAD_REQUEST", f"subset は {list(mass_accuracy.SUBSETS)}",
                           {"field": "subset"})
    if subset == "msms_matched" and rank != 1:
        return mztab_error("DATASET_BAD_REQUEST",
                           "subset=msms_matched は rank 1 の同定だけ（SML の照合の旗は rank 1 の同定を"
                           "表す）", {"field": "rank"})
    available = subset != "msms_matched" or mass_accuracy.has_match_information(ds)
    rows = mass_accuracy.sme_ppm(ds, rank=rank, subset=subset) if available else []
    summary = mass_accuracy.band_summary(rows, [float(e) for e in edges])
    extra = {} if available else {
        "subset_available": False,
        "note": "この mzTab の SML に opt_global_spectrum_matched が無く、MS/MS で照合された同定を"
                "選べない（n=0）。subset=rank1 を使う"}
    return json_payload(round_floats({"status": "success", "rank": rank, "subset": subset,
                                      **extra,
                                      "band_edges": sorted(float(e) for e in edges),
                                      **summary,
                                      "injection_trend": mass_accuracy.injection_trend(ds, rows)},
                                     digits=4))


def _novel_result_stats(ds, result: dict, source: dict, *,
                        dataset_hash: str | None = None) -> tuple[dict, dict]:
    """2 群比較の結果からの統計（特徴 ID → log2fc・p・q）と、割れのまとめの記録。

    `merge_split` が真でまとめる組が 1 つでもあれば、すべての特徴の統計をまとめた行列で検定し
    直した値にする（q はまとめた後の全特徴の BH。`dataset_feature_search` の `merge_split` と
    同じ規則・同じ実装）。組が無ければ結果の値のまま。
    """
    from metabolomix.analysis import feature_search

    def _rows(found):
        return {row["feature"]: {"log2fc_test_over_ref": row.get("log2fc"),
                                 "p": row.get("p"), "q": row.get("q")}
                for row in found.get("results") or []}

    stats_rows = _rows(result)
    enabled = bool(source.get("merge_split"))
    merge = {"enabled": enabled, "source": None, "groups": []}
    if enabled:
        plan = feature_search.merge_plan(ds, alignment_qc_path=source.get("alignment_qc"),
                                         field="statistics.alignment_qc",
                                         dataset_hash=dataset_hash)
        merged = feature_search.merged_results(ds, [result], plan["groups"], ds.feature_ids)
        for redone in merged["results"].values():
            stats_rows = _rows(redone)
        merge.update(source=plan["source"], groups=list(merged["touched"].values()))
    return stats_rows, merge


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=True),
          structured_output=False)
def dataset_novel_candidates(spec: dict, dataset_handle: str | None = None,
                             output_path: str | None = None) -> str:
    """新しい知見の候補の一覧（novel-candidates.v1。再解析で使う）。

    SME（rank 1）の同定があり、質量誤差が m/z 帯ごとの許容幅（`ppm.default` と
    `ppm.bands: [{max_mz, ppm}]`）に入り、q が `q_max` 未満の特徴から、`exclude`
    （`feature_ids`・`inchikey_blocks`）を除いた一覧。統計は `statistics` が
    `{"result_id"}`（`dataset_differential` の 2 群比較）か `{"analysis_spec", "analysis_id"}`
    （`dataset_analysis_run` と同じ spec のその解析）。q の昇順。`rule_hash` を返す。
    論文の同定リストとの突き合わせはしない（再解析の側）。
    割れたアライメントスポットの組（`dataset_alignment_qc` の `merge_groups`）は 1 候補に数え
    （`feature_id` は残す小さい ID、`merged_from` に 2 つの ID、同定と質量誤差は SME を持つ方の
    スポット＝`sme_feature_id`。両方にあれば残す ID）、統計はまとめた行列で検定し直した値。
    `{"result_id"}` では `statistics.merge_split`（既定 true）・`statistics.alignment_qc`
    （`dataset_alignment_qc` の JSON。省略時はその場で探す）、analysis_spec では spec の
    `features.merge_split` で決まる。
    spec の `labels`（`{path, name, feature_column="feature_id", label_column="level",
    order=["A","B","C","D"]}`）は外部の印の TSV。行に `labels.<name>` の列（特徴ごとの最も強い
    水準、割れの組はどちらかの最も強い水準）と水準ごとの件数 `labels.counts`（全体と上位 20）を
    足すだけで、並び・件数・除外・`rule_hash` は変えない。
    output_path: JSON をここへ書き、要約だけを返す（行が 50 を超えるときは必須）。
    dataset_handle: ハンドルのデータセット。省略時は対話のセッション。
    """
    from metabolomix.analysis import analysis_spec, mass_accuracy
    from metabolomix.analysis.feature_search import select_results
    from metabolomix.analysis.result_state import dataset_fingerprint

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds, error = _dataset_or_error(holder, True)
    if error:
        return error
    try:
        normalized = mass_accuracy.normalize_novel_spec(spec)
        source = normalized["statistics"]
        if "result_id" in source:
            # データの指紋は 1 回だけ計算し、結果の照合と割れの照合で使い回す
            dataset_hash = dataset_fingerprint(ds)
            chosen = select_results(ds, [source["result_id"]], dataset_hash=dataset_hash)
            if not chosen["results"]:
                return missing_state(
                    "dataset_differential_result", ["dataset_differential"],
                    f"result_id={source['result_id']} の今も有効な 2 群比較がありません"
                    f"（無い: {chosen['missing']} / 古い: {chosen['stale']}）。")
            result = chosen["results"][0]
            group_ref, group_test = result.get("a"), result.get("b")
            stats_rows, merge = _novel_result_stats(ds, result, source,
                                                    dataset_hash=dataset_hash)
        else:
            run = analysis_spec.run_spec(source["analysis_spec"], ds)
            analysis = next(a for a in run["analyses"] if a["analysis_id"] == source["analysis_id"])
            group_ref, group_test = run["group_ref"], run["group_test"]
            stats_rows = {row["feature_id"]: row for row in analysis["rows"]
                          if row["status"] == "tested"}
            merge = {key: run["merge_split"][key] for key in ("enabled", "source", "groups")}
        # まとめないとき（merge_split=false）は行も戻り値も v7 と同じ形
        rows = mass_accuracy.novel_candidates(
            ds, normalized, stats_rows, group_ref=group_ref, group_test=group_test,
            merged_from={g[0]: g for g in merge["groups"]} if merge["enabled"] else None)
        labels_info = None
        if "labels" in normalized:
            spec_labels = normalized["labels"]
            read = mass_accuracy.read_labels(spec_labels)
            counts = mass_accuracy.attach_labels(rows, read["by_feature"], name=spec_labels["name"],
                                                 order=spec_labels["order"])
            known = set(ds.feature_ids)
            labels_info = {"name": spec_labels["name"], "sha256": read["sha256"],
                           "order": spec_labels["order"], "n_rows": read["n_rows"],
                           "n_unmatched": sum(1 for f in read["by_feature"] if f not in known),
                           "counts": counts}
        # まとめた組の hash は、まとめるとき（組が 0 でも）だけ規則に入れる
        groups_hash = (analysis_spec.merge_groups_hash(merge["groups"]) if merge["enabled"]
                       else None)
        merge_info = ({"enabled": True, "source": merge["source"],
                       "n_groups": len(merge["groups"]),
                       "groups_hash": groups_hash} if merge["enabled"] else None)
        result_payload = {"schema": mass_accuracy.NOVEL_SCHEMA,
                          "rule": normalized,
                          "rule_hash": mass_accuracy.rule_hash(normalized, groups_hash),
                          "group_ref": group_ref, "group_test": group_test,
                          **({"merge_split": merge_info} if merge_info else {}),
                          **({"labels": labels_info} if labels_info else {}),
                          "rows": rows}
        if output_path is None and len(rows) > MAX_INLINE_ROWS:
            raise DomainError("ANALYSIS_OUTPUT_REQUIRED",
                              f"行が {MAX_INLINE_ROWS} を超える（{len(rows)} 件）。output_path を"
                              "指定してください。", {"n_rows": len(rows)})
        written = write_result(output_path, result_payload) if output_path else None
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    columns = ("feature_id", "name", "inchikey", "adduct", "mz", "rt", "ppm", "group_ref",
               "group_test", "log2fc_test_over_ref", "p", "q", "direction")
    if merge_info:
        columns += ("merged_from",)
    if labels_info:
        columns += (f"labels.{labels_info['name']}",)
    return json_payload({
        "status": "success",
        "rule_hash": result_payload["rule_hash"],
        "group_ref": group_ref, "group_test": group_test,
        **({"merge_split": merge_info} if merge_info else {}),
        **({"labels": labels_info} if labels_info else {}),
        "n_rows": len(rows),
        **({"output_path": written} if written else
           {"rows": rows_table([{**r, "mz": round(r["mz"], 4), "ppm": round(r["ppm"], 2),
                                 "rt": None if r["rt"] is None else round(r["rt"], 3)}
                                for r in rows], columns)}),
    })


# ---------- アライメントの品質 ----------

#: 戻り値に group_aligned の組を TSV で並べる上限（チームの監査に回す組。超えればファイルだけ）。
_GROUP_ALIGNED_COLUMNS = ("pair_id", "feature_a", "feature_b", "mz", "rt_a", "rt_b", "delta_rt",
                          "adduct", "name_a", "name_b")


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=True),
          structured_output=False)
def dataset_alignment_qc(output_path: str | None = None, ppm: float = 15.0,
                         max_rt_gap: float = 1.0, max_interleaved_rt_gap: float = 0.15,
                         min_block_size: int = 2, dataset_handle: str | None = None) -> str:
    """1 つの物質が 2 つ以上のアライメントスポットに割れた組を探す（alignment-qc.v1。再解析で使う）。

    候補は同じ付加体・|Δm/z| ≤ ppm・0.05 < ΔRT ≤ max_rt_gap 分の 2 スポットで、試料ごとの
    受け持ち（検出状態があれば検出、無ければ値の大きい方）が補い合う組を割れとする。型は
    `injection_block`（片方が注入順で連続した塊だけを持つ。両方が min_block_size（既定 2）
    試料以上を受け持つときだけ。群と重なっても優先）・
    `group_aligned`（受け持ちが群と一致。まとめない）・`interleaved`。`mergeable` は注釈が
    両立する `injection_block` と、ΔRT ≤ max_interleaved_rt_gap（既定 0.15 分）の `interleaved`。
    そのうち互いに |ΔRT| が最も近い相手どうしの 2 スポット（`merge_groups`）だけを、統計の側が
    試料ごとに大きい方の値で 1 つにまとめる（ほかは `merge_reason: not_mutual_nearest` で印だけ）。注入順はハンドルの実験情報シート、無ければファイル名の順
    （`injection_order_source`）。欄の意味は `lipidmix://docs/output-format/reanalysis`。

    output_path: 全体（組・連結成分・試料ごとの受け持ち）を JSON でここへ書く
        （例 `results/alignment_qc_pos.json`）。戻り値は要約と `group_aligned` の組だけ。
    dataset_handle: ハンドルのデータセット。省略時は対話のセッション。状態は変えない。
    """
    from metabolomix.analysis import alignment_qc

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds, error = _dataset_or_error(holder, True)
    if error:
        return error
    try:
        qc = round_floats(alignment_qc.find_split_pairs(
            ds, ppm=ppm, max_rt_gap=max_rt_gap, max_interleaved_rt_gap=max_interleaved_rt_gap,
            min_block_size=min_block_size),
            digits=6)
        written = write_result(output_path, qc) if output_path else None
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    flagged = [p for p in qc["pairs"] if p["type"] == "group_aligned"]
    payload = {"status": "success", "schema": qc["schema"],
               "injection_order_source": qc["injection_order_source"],
               "signal": qc["params"]["signal"],
               "summary": round_floats(qc["summary"], digits=4)}
    if written:
        payload["output_path"] = written
    if flagged and len(flagged) <= MAX_INLINE_ROWS:
        payload["group_aligned"] = rows_table(
            [{"pair_id": p["pair_id"], "feature_a": p["features"][0],
              "feature_b": p["features"][1], "mz": round(p["mz"][0], 4),
              "rt_a": round(p["rt"][0], 3), "rt_b": round(p["rt"][1], 3),
              "delta_rt": round(p["delta_rt"], 3), "adduct": p["adduct"],
              "name_a": p["names"][0], "name_b": p["names"][1]} for p in flagged],
            _GROUP_ALIGNED_COLUMNS)
    return json_payload(payload)
