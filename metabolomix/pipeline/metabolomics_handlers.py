"""v2（LC–MSメタボロミクス）の工程handler（spec §6.2・§8〜§11）。

**session を一切 import しない。** ここは独立 worker プロセスが回す層で、
進行状況は `pipeline-run.json` と、そのプロセスだけが持つ `runtime` dict に住む
（`metabolomix/pipeline/` の規約。`tests/test_pipeline_engine.py` の AST テストが
固定している）。

## v1 と同じ stage 名でも中身が違う

`preprocess` / `export` / `report` は v1 と v2 で別の計算をする。1つの辞書へ
後勝ちで入れると、片方の pipeline が黙ってもう片方の計算を走らせる——しかも
成果物の形が違うので、気付くのは最後のレポートを読んだ人になる。
`by_schema` で要求の schema を見て選ぶ（`service.build_handlers` が1か所で組む）。

## 工程の繋ぎ方

    load_assay_evidence → resolve_feature_bindings → qc_raw
      → preprocess → qc_processed → statistics:* → export:* → report

- `qc_raw` が **filter 前に** QC 評価集合を固定する。以後の QC はその集合を
  再利用し、filter で縮まない（spec §9）。
- `preprocess` は recipe ごとに行列を作り、**補完しない**。`ds.pp_matrix` へは
  書かない——recipe が2本あるので、単一スロットは必ずどちらかを失う。
- `qc_processed` が処理後 QC を評価してから `finalize_matrix` を呼び、
  「QC 結果」と「最終行列」を**1つの stage outcome** で確定する。片方だけが
  記録に残る瞬間を作らない。
- 統計は最終行列だけを参照する（`matrix_id` を結果に書く）。

## needs_input は止めたまま

binding が未解決なら `FEATURE_BINDING_UNRESOLVED` で止める。同じ未解決条件で
自動再試行しない——同じ答えが返るだけで、記録だけが増える。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from metabolomix.analysis import assay_evidence, assay_qc, feature_bindings
from metabolomix.analysis import feature_export, matrix_state, statistics_v2
from metabolomix.analysis.preprocess_policy import explicit_policy
from metabolomix.analysis.result_state import dataset_fingerprint
from metabolomix.analysis.sample_manifest import validate_standard_assays
from metabolomix.console import profile_adapter, profiles
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import evidence as mztab_evidence
from metabolomix.pipeline import report as report_mod
from metabolomix.pipeline import stage_plan, store

__all__ = ["build_handlers", "by_schema", "restore_results"]

#: `restore_results` が復元する成果物（output_name → runtime のキー）。
_RESTORABLE = {
    "assay_evidence": "assay_evidence",
    "feature_bindings": "bindings",
    "qc_population": "qc_population",
    "qc": "qc",
}


def by_schema(v1_handler, v2_handler):
    """要求の schema で v1/v2 の handler を選ぶ（stage 名が重なるキー用）。"""
    def dispatch(context: dict) -> dict:
        handler = (v2_handler if stage_plan.is_v2_request(context.get("request"))
                   else v1_handler)
        return handler(context)
    dispatch.__name__ = getattr(v2_handler, "__name__", "dispatch")
    dispatch.__doc__ = v2_handler.__doc__
    return dispatch


# ---------- 共通部品 ----------

def _ok(result_refs=None, warnings=None, **extra) -> dict:
    return {"status": "succeeded", "result_refs": list(result_refs or []),
            "warnings": list(warnings or []), "error": None, **extra}


def _needs_input(code: str, message: str, details: dict) -> dict:
    return {"status": "needs_input", "result_refs": [], "warnings": [],
            "error": {"code": code, "message": message, "details": details}}


def _profile(context: dict) -> dict:
    """要求が指す profile を読む（runtime に1度だけ載せる）。"""
    runtime = context["runtime"]
    if "profile" not in runtime:
        request = context["request"]
        runtime["profile"] = profiles.load_profile(
            Path(request["profile_file"]),
            request.get("execution_purpose", "routine"))
    return runtime["profile"]


def _adapter(context: dict) -> dict:
    runtime = context["runtime"]
    if "adapter" not in runtime:
        profile = _profile(context)
        version = str(profile["software"]["msdial_version"])
        major = version.split(".", 1)[0]
        runtime["adapter"] = profile_adapter.adapter_capabilities(f"msdial{major}")
    return runtime["adapter"]


def _persist(context: dict, output_name: str, kind: str, result_id: str,
             data: dict, parent_ids=None) -> dict:
    return report_mod.persist_result(context["pipeline_root"], {
        "output_name": output_name, "kind": kind, "result_id": result_id,
        "parent_ids": list(parent_ids or []),
        "data": data,
        "request_revision": (context.get("request_meta") or {}).get("revision"),
    })


#: 行列の配列。要約 JSON には入れず、`save_matrix` が npz へ書く。
_MATRIX_ARRAY_KEYS = ("values", "detected_mask", "imputed_mask", "locked_mask")

#: pipeline_root 相対の行列置き場。`persist_result` が書く
#: `results/<result_id>.json` と衝突しないよう1段下げる（`save_matrix` も
#: `<matrix_id>.json` を書くため、同じ場所だと片方が片方を潰す）。
_MATRIX_SUBDIR = ("results", "matrices")


def _persist_matrix(context: dict, matrix: dict, parent_ids=None) -> dict:
    """行列を npz+meta で保存し、要約 JSON の成果物 ref を返す。

    要約に値は入れない（run record と戻り値を肥大させない）。数値は
    `matrix_state.save_matrix` が別ファイルへ書き、その参照だけを
    `storage` として要約に載せる。完全性はこれで一本に繋がる:
    run record（追記専用）の `storage.meta_hash` → meta の `array_hashes`
    → npz。どこを書き換えても `load_matrix` が拒否する。
    """
    directory = Path(context["pipeline_root"]).joinpath(*_MATRIX_SUBDIR)
    directory.mkdir(parents=True, exist_ok=True)
    storage = matrix_state.save_matrix(matrix, directory)
    summary = _jsonable({k: v for k, v in matrix.items()
                         if k not in _MATRIX_ARRAY_KEYS})
    return _persist(context, "matrix", "analysis_matrix", matrix["matrix_id"],
                    {**summary, "storage": storage}, parent_ids=parent_ids)


def _jsonable(value):
    """numpy 配列を持つ結果を JSON 保存できる形へ落とす（mask は list にする）。"""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _recipes(context: dict) -> dict:
    """profile の matrix_recipes に、要求の recipe 別 override を重ねる。"""
    profile = _profile(context)
    overrides = (context["request"].get("preprocess") or {})
    recipes = {}
    for recipe_id, recipe in (profile.get("matrix_recipes") or {}).items():
        merged = dict(recipe)
        merged.update(overrides.get(recipe_id) or {})
        merged["recipe_id"] = recipe_id
        recipes[recipe_id] = merged
    return recipes


def _qc_policy(context: dict) -> list[dict]:
    """profile の qc_policy（metric キーの dict）を metric 入りの list にする。"""
    profile = _profile(context)
    return [{"metric": metric, **entry}
            for metric, entry in (profile.get("qc_policy") or {}).items()]


def _standard_assay_ids(context: dict) -> dict:
    """要求の `standard_assays`（sample_id）を assay ID へ解決する。

    解決は manifest 経由（`sample_id` → 行位置 → `ds.sample_assay_ids`）。
    role=standard・include=true の確認は `sample_manifest` が唯一の実装。
    """
    metadata = context["runtime"]["metadata"]
    ds = context["runtime"]["dataset"]
    requested = context["request"].get("standard_assays") or {}
    validate_standard_assays(metadata, requested)

    position = {row["sample_id"]: i for i, row in enumerate(metadata)}
    assay_ids = list(getattr(ds, "sample_assay_ids", []) or [])
    resolved: dict[str, list[str]] = {}
    for target_id, sample_ids in requested.items():
        ids = []
        for sample_id in sample_ids:
            index = position.get(sample_id)
            if index is None or index >= len(assay_ids):
                raise DomainError(
                    "STANDARD_ASSAY_INVALID",
                    f"standard_assaysのsample_idを解決できません: {sample_id!r}",
                    {"target_id": target_id, "sample_id": sample_id})
            ids.append(assay_ids[index])
        resolved[target_id] = ids
    return resolved


# ---------- prepare_inputs / resolve_metadata への v2 の上乗せ ----------

def _assert_raw_unchanged(context: dict, plan: dict) -> None:
    """前回固定した profile snapshot の raw hash と、今回の再解決結果を突き合わせる。

    v2 は `resolve_profile_inputs` で全 raw を hash しており、その値は profile
    成果物の snapshot に入っている。突き合わせなければ、入力が差し替わっても
    気付かないまま解析が進む。v1 の `verify_inputs` が raw の内容 hash を検証
    しないのは仕様上の決定（「rawの全量ハッシュは初期版の必須条件にしない」）
    なので、そちらの契約は変えず、**自分が既に持っている hash** をここで使う。

    前回の記録が無い初回は何もしない（比べる相手がいない）。
    """
    previous = store.read_result_data(
        Path(context["pipeline_root"]), context.get("results") or [], "profile")
    if not previous:
        return
    before = (previous[-1].get("snapshot") or {}).get("raw_files") or {}
    baseline = {f["relative_path"]: f.get("sha256")
                for f in (before.get("files") or []) if f.get("sha256")}
    if not baseline:
        return
    current = {f["relative_path"]: f.get("sha256")
               for f in ((plan.get("raw_files") or {}).get("files") or [])}

    changed = sorted(rel for rel, digest in baseline.items()
                     if rel in current and current[rel] != digest)
    missing = sorted(rel for rel in baseline if rel not in current)
    added = sorted(rel for rel in current if rel not in baseline)
    if not (changed or missing or added):
        return
    raise DomainError(
        "INPUT_CHANGED",
        "固定済みの生データが変化しています（内容・件数のいずれか）。"
        "同じ run では入力を差し替えられません——別の run を起こしてください。",
        {"changed": changed, "missing": missing, "added": added})


def snapshot_profile_outcome(context: dict, outcome: dict) -> dict:
    """解決済み profile と実行環境 manifest を run へ固定し、ref を足す。

    v2 の必須成果物（spec §11「解決済みprofile、実行条件/依存manifest」）。
    v1 の `prepare_input`（raw の実配置）の結果へ**足す**形にするのは、
    どちらも「実行前に入力を固定する」同じ工程だから——別 stage にすると、
    片方だけ成功した状態が記録に残りうる。
    """
    profile = _profile(context)
    source_root = Path(context["identity"]["source_root"])
    profile_path = Path(context["request"]["profile_file"])
    if not profile_path.is_absolute():
        profile_path = source_root / profile_path
    plan = profiles.resolve_profile_inputs(profile, profile_path.parent, raw_root=source_root)
    _assert_raw_unchanged(context, plan)
    snapshot = profiles.snapshot_profile(plan, Path(context["pipeline_root"]))

    refs = list(outcome.get("result_refs") or [])
    refs.append(_persist(context, "profile", "profile",
                         f"res_profile_{context['pipeline_id']}",
                         _jsonable({"profile": profile, "snapshot": snapshot})))
    refs.append(_persist(context, "execution_manifest", "execution_manifest",
                         f"res_execution_manifest_{context['pipeline_id']}",
                         _jsonable({
                             "execution_environment": plan.get("execution_environment"),
                             "dependencies": plan.get("dependencies"),
                             "method": plan.get("method"),
                             "raw": plan["raw_files"]})))
    context["runtime"]["profile_plan"] = plan
    return {**outcome, "result_refs": refs}


def sample_manifest_outcome(context: dict, outcome: dict) -> dict:
    """解決済みの試料対応表を成果物として固定する（spec §11「試料対応表」）。

    v1 は `record["inputs"]["manifest"]` にだけ書く。v2 はそれを result ref に
    もする——必須成果物は hash 付きで名指しできなければ「揃っている」と言えない。
    """
    metadata = context["runtime"].get("metadata") or []
    refs = list(outcome.get("result_refs") or [])
    refs.append(_persist(context, "sample_manifest", "sample_manifest",
                         f"res_sample_manifest_{context['pipeline_id']}",
                         _jsonable({"rows": metadata, "n_rows": len(metadata)})))
    return {**outcome, "result_refs": refs}


# ---------- load_assay_evidence ----------

def _handle_load_assay_evidence(context: dict) -> dict:
    """注入ごとの RT/m/z 証拠を読む（取れなくても理由 JSON を必ず残す）。

    証拠が無いこと自体は異常ではない（別アライメントの `.arf` しか無い、等）。
    ここで止めず、`availability=False` と理由を成果物として残して先へ進む
    ——止めると、証拠を必要としない profile（`mass_rt` だけ）まで動かなくなる。
    必要かどうかを判断するのは binding と QC の側。
    """
    ds = context["runtime"]["dataset"]
    adapter = _adapter(context)

    candidates = []
    for path in (getattr(ds, "artifact_paths", {}) or {}).get("peak_matrix_source", []):
        candidates.append(Path(path))
    for path_text in (getattr(ds, "source_files", {}) or {}):
        candidates.extend(mztab_evidence.arf_candidates(
            path_text, getattr(ds, "artifact_paths", None)))

    result = None
    for candidate in candidates:
        result = assay_evidence.build_assay_evidence(ds, candidate, adapter)
        if result["availability"]:
            break
    if result is None:
        result = {
            "schema": assay_evidence.EVIDENCE_SCHEMA,
            "dataset_id": getattr(ds, "dataset_id", None),
            "availability": False,
            "reasons": [{"code": "no_candidate_artifact",
                         "message": "注入ごとの証拠を読める成果物が見つかりません。",
                         "details": {}}],
            "rows": [], "coverage": {}, "detail": {},
            "source_artifact_hash": None, "source_locator": None,
            "reader_version": None, "rt_unit": None,
        }

    context["runtime"]["assay_evidence"] = result
    ref = _persist(context, "assay_evidence", "assay_evidence",
                   f"res_assay_evidence_{context['pipeline_id']}",
                   _jsonable({k: v for k, v in result.items() if k != "rows"}
                             | {"n_rows": len(result.get("rows") or [])}))
    warnings = [{"code": "ASSAY_EVIDENCE_UNAVAILABLE",
                 "message": reason.get("message", reason.get("code"))}
                for reason in (result.get("reasons") or [])]
    return _ok([ref], warnings)


# ---------- resolve_feature_bindings ----------

def _binding_overrides(context: dict, ds) -> dict | None:
    """要求の `feature_bindings` payload を `bind_features` の overrides へ移す。

    要求側（`request_v2`）の形は `{dataset_hash, selections: {target_id:
    {feature_id, reason}}}`——「どの dataset に対する選択か」を payload 自身が
    宣言する。ここで**実際の dataset の指紋と突き合わせる**: 一致しなければ、
    前のバッチで決めた選択を今のバッチへ持ち込もうとしている。同じ feature_id が
    別の化合物を指しうるので、黙って適用しない。
    """
    payload = context["request"].get("feature_bindings")
    if not payload:
        return None
    expected = payload.get("dataset_hash")
    actual = dataset_fingerprint(ds)
    if expected != actual:
        raise DomainError(
            "FEATURE_BINDING_OVERRIDE_INVALID",
            "feature_bindingsが別のデータセットに対する選択です"
            "（dataset_hashが今回のデータセットと一致しません）。"
            "同じprofile規則でこのバッチに対応付け直してください。",
            {"expected_dataset_hash": actual, "payload_dataset_hash": expected})
    return {target_id: {**selection, "dataset_id": getattr(ds, "dataset_id", None)}
            for target_id, selection in (payload.get("selections") or {}).items()}


def _handle_resolve_feature_bindings(context: dict) -> dict:
    ds = context["runtime"]["dataset"]
    profile = _profile(context)
    evidence = context["runtime"].get("assay_evidence") or {}
    standard_assays = _standard_assay_ids(context)
    overrides = _binding_overrides(context, ds)

    result = feature_bindings.bind_features(ds, profile, evidence,
                                            standard_assays, overrides)
    context["runtime"]["bindings"] = result

    # 未解決でも結果は**必ず**保存する。全候補と不採用理由が spec §8.1 の
    # 必須成果物で、利用者はそれを見てから `feature_bindings` で選ぶ
    # ——保存しないと「何を選べるのか」が記録のどこにも無い。
    ref = _persist(context, "feature_bindings", "feature_bindings",
                   f"res_bindings_{context['pipeline_id']}", _jsonable(result))

    if result["status"] != "resolved":
        # 未解決のまま先へ進めない。自動再試行もしない——同じ入力で同じ答えが
        # 返るだけで、attempt ディレクトリと記録だけが増える。
        outcome = _needs_input(
            "FEATURE_BINDING_UNRESOLVED",
            "profileのfeature_targetsをこのバッチのfeatureへ一意に対応付けられません。"
            "候補と不採用理由を確認し、feature_bindingsで選択するかprofileを"
            "改訂してください。",
            {"unresolved": result["unresolved"],
             "bindings": {tid: {"status": b["status"], "reason": b["reason"]}
                          for tid, b in result["bindings"].items()}})
        return {**outcome, "result_refs": [ref]}

    return _ok([ref])


# ---------- qc_raw ----------

def _raw_recipe() -> dict:
    """元データ QC 用の行列 recipe。filter も補完もしない（評価集合を縮めない）。"""
    return {"recipe_id": "__raw__", "base": "peak_height", "normalize": "none",
            "drift_correct": False, "filter": None, "impute": "none"}


def _handle_qc_raw(context: dict) -> dict:
    """filter 前の行列で QC を評価し、**評価集合を固定**する（spec §9）。

    評価集合は profile の recipe に依存しない（`__raw__` は filter も補完も
    持たない固定 recipe で、全 feature × 全 assay が母集団）。だから recipe を
    変えても元データ QC の対象集合は動かず、`stage_plan` の `qc_raw_scope`
    token をここから立てる必要は無い——recipe 変更は `preprocess` 以降だけを
    無効化すれば足りる（Task 4 から持ち越した判断の決着）。
    """
    ds = context["runtime"]["dataset"]
    evidence = context["runtime"].get("assay_evidence") or {}
    bindings = context["runtime"].get("bindings") or {}
    metadata = context["runtime"]["metadata"]

    eligibility = np.ones(len(getattr(ds, "feature_ids", []) or []), dtype=bool)
    raw_matrix = matrix_state.make_matrix(ds, _raw_recipe(), bindings, evidence,
                                          eligibility)
    context["runtime"]["raw_matrix"] = raw_matrix

    population = _target_population(context)
    qc = assay_qc.evaluate_qc(raw_matrix, evidence, metadata,
                              _qc_policy(context), population)
    context["runtime"]["qc_population"] = qc["population"]
    context["runtime"]["qc_raw"] = qc

    ref = _persist(context, "qc_population", "qc_population",
                   f"res_qc_population_{context['pipeline_id']}",
                   _jsonable({"population": qc["population"],
                              "raw_batch_status": qc["batch_status"],
                              "matrix_id": raw_matrix["matrix_id"]}))
    return _ok([ref])


def _target_population(context: dict) -> dict | None:
    """標準品 metric が要る期待値（target→feature と RT/m/z）を population へ載せる。

    `assay_qc` は profile を読まない（数値層に profile を持ち込まない）ので、
    解決済みの対応と期待値をここで渡す。binding が無い target は載せない
    ——載せると「期待値ゼロで全部合格」になる。
    """
    profile = _profile(context)
    bindings = (context["runtime"].get("bindings") or {}).get("bindings") or {}
    targets = {}
    for target_id, binding in bindings.items():
        if binding.get("status") != "resolved":
            continue
        target = (profile.get("feature_targets") or {}).get(target_id) or {}
        targets[target_id] = {
            "feature_id": binding["selected_feature_id"],
            "expected_rt_min": target.get("expected_rt_min"),
            "expected_mz": target.get("expected_mz"),
        }
    if not targets:
        return None
    ds = context["runtime"]["dataset"]
    return {entry["metric"]: {
        "feature_ids": list(getattr(ds, "feature_ids", []) or []),
        "assay_ids": list(getattr(ds, "sample_assay_ids", []) or []),
        "targets": targets,
        "hash": None,
    } for entry in _qc_policy(context)}


# ---------- preprocess ----------

def _handle_preprocess(context: dict) -> dict:
    """recipe ごとに解析行列を作る（補完前・ds は書き換えない）。"""
    ds = context["runtime"]["dataset"]
    evidence = context["runtime"].get("assay_evidence") or {}
    bindings = context["runtime"].get("bindings") or {}
    eligibility = np.ones(len(getattr(ds, "feature_ids", []) or []), dtype=bool)

    matrices: dict[str, dict] = {}
    plans: dict[str, dict] = {}
    refs = []
    for recipe_id, recipe in _recipes(context).items():
        # v2 は profile が全stepを明示する。conservative-v1 の自動判断を
        # 重ねない（前提不足は make_matrix が QC_PREREQUISITE_MISSING で止める）。
        plans[recipe_id] = explicit_policy(recipe)
        matrix = matrix_state.make_matrix(ds, recipe, bindings, evidence,
                                          eligibility)
        matrices[recipe_id] = matrix
        refs.append(_persist_matrix(context, {**matrix, "policy": plans[recipe_id]}))

    context["runtime"]["matrices"] = matrices
    context["runtime"]["matrix_plans"] = plans
    warnings = [{"code": "MATRIX_CAVEAT", "message": caveat}
                for matrix in matrices.values()
                for caveat in matrix.get("caveats", [])]
    return _ok(refs, warnings)


# ---------- qc_processed ----------

def _handle_qc_processed(context: dict) -> dict:
    """処理後 QC を評価してから補完し、QC と最終行列を1つの outcome で確定する。"""
    evidence = context["runtime"].get("assay_evidence") or {}
    metadata = context["runtime"]["metadata"]
    matrices = context["runtime"]["matrices"]
    population = context["runtime"].get("qc_population")
    recipes = _recipes(context)

    qc_results: dict[str, dict] = {}
    final: dict[str, dict] = {}
    refs = []
    for recipe_id, matrix in matrices.items():
        qc = assay_qc.evaluate_qc(matrix, evidence, metadata,
                                  _qc_policy(context), population)
        qc_results[recipe_id] = qc
        # QC で落ちた feature は eligibility にだけ反映する。QC は再集計しない。
        exclude = set(qc["eligibility_suggestion"]["exclude_feature_ids"])
        eligibility = np.array([fid not in exclude
                                for fid in matrix["feature_ids"]], dtype=bool)
        finalized = matrix_state.finalize_matrix(
            matrix, eligibility, recipes[recipe_id].get("impute", "none"))
        final[recipe_id] = finalized

        refs.append(_persist(
            context, "qc", "assay_qc", f"res_qc_{recipe_id}_{context['pipeline_id']}",
            _jsonable(qc), parent_ids=[matrix["matrix_id"]]))
        refs.append(_persist_matrix(context, finalized,
                                    parent_ids=[matrix["matrix_id"]]))

    context["runtime"]["qc"] = qc_results
    context["runtime"]["final_matrices"] = final
    warnings = [{"code": "QC_WARNING", "message": message}
                for qc in qc_results.values() for message in qc["warnings"]]
    return _ok(refs, warnings)


# ---------- statistics ----------

def _specification(context: dict) -> dict:
    statistic_id = context["statistic_id"]
    for spec in (context["request"].get("statistics") or []):
        if spec.get("statistic_id") == statistic_id:
            return spec
    raise DomainError(
        "STATISTIC_SPECIFICATION_INVALID",
        f"要求に無い統計です: {statistic_id!r}",
        {"statistic_id": statistic_id})


def _target_features(context: dict) -> dict:
    return feature_bindings.resolved_targets(context["runtime"].get("bindings"))


def _handle_statistics(context: dict) -> dict:
    """指定統計を**最終行列**に対して実行する（評価不能でも結果JSONを残す）。"""
    specification = _specification(context)
    recipe_id = specification.get("matrix_recipe_id") or "default"
    final = context["runtime"].get("final_matrices") or {}
    if recipe_id not in final:
        return _needs_input(
            "MATRIX_NOT_FOUND",
            f"統計が参照するmatrix_recipe_idの行列がありません: {recipe_id!r}",
            {"statistic_id": specification.get("statistic_id"),
             "matrix_recipe_id": recipe_id, "available": sorted(final)})

    result = statistics_v2.run_statistic(
        final[recipe_id], specification, context["runtime"]["metadata"],
        _target_features(context))
    context["runtime"].setdefault("statistics", {})[
        specification["statistic_id"]] = result

    ref = _persist(context, f"statistic:{specification['statistic_id']}",
                   "statistic", f"res_statistic_{specification['statistic_id']}",
                   _jsonable(result), parent_ids=[final[recipe_id]["matrix_id"]])
    return _ok([ref])


# ---------- export ----------

def _statistic_result(context: dict, statistic_id: str) -> dict | None:
    """統計結果を runtime から、無ければ**保存済み result から hash 照合で**取る。

    `statistics:<id>` stage は成果物 JSON が完全なので、再開時に skip されうる
    （行列と違って値を丸ごと保存してある）。そのとき runtime は空なので、
    記録から読み直す——再計算して別の result ID を作らない。
    """
    cached = (context["runtime"].get("statistics") or {}).get(statistic_id)
    if cached is not None:
        return cached
    restored = _restore_output(context, f"statistic:{statistic_id}")
    if restored is not None:
        context["runtime"].setdefault("statistics", {})[statistic_id] = restored
    return restored


def _restore_output(context: dict, output_name: str):
    """`context["results"]` の最新 ref を hash 照合して読む（無ければ None）。"""
    payloads = store.read_result_data(
        Path(context["pipeline_root"]), context.get("results") or [], output_name)
    return payloads[-1] if payloads else None


def _handle_export(context: dict) -> dict:
    """全feature表と統計表を書く（統計1件につき1回。全feature表は上書き更新）。"""
    ds = context["runtime"]["dataset"]
    statistic_id = context["statistic_id"]
    result = _statistic_result(context, statistic_id)
    if result is None:
        return _needs_input(
            "STATISTIC_RESULT_MISSING",
            f"統計結果がまだありません: {statistic_id!r}",
            {"statistic_id": statistic_id})

    out_dir = Path(context["pipeline_root"]) / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    matrices = list((context["runtime"].get("final_matrices") or {}).values())
    feature_summary = feature_export.export_features(
        ds, matrices, out_dir / "feature-table.tsv")
    statistic_summary = feature_export.export_statistic(
        result, out_dir / f"statistic-{statistic_id}.tsv")

    refs = [
        report_mod.persist_result(context["pipeline_root"], {
            "output_name": "feature_table", "kind": "feature_table",
            "result_id": f"res_feature_table_{statistic_id}",
            "path": feature_summary["path"]}),
        report_mod.persist_result(context["pipeline_root"], {
            "output_name": f"statistic_table:{statistic_id}", "kind": "statistic_table",
            "result_id": f"res_statistic_table_{statistic_id}",
            "path": statistic_summary["path"]}),
    ]
    return _ok(refs)


# ---------- report ----------

def _axis_statuses(context: dict) -> dict:
    statistics = _all_statistics(context).values()
    analysis = report_mod.analysis_status(
        [{"status": "ready" if r.get("status") == "completed" else "not_evaluable"}
         for r in statistics])
    qc_results = (context["runtime"].get("qc") or {}).values()
    statuses = {qc["batch_status"] for qc in qc_results}
    qc_status = ("fail" if "fail" in statuses
                 else "not_evaluable" if ("not_evaluable" in statuses or not statuses)
                 else "pass")
    return {"qc_status": qc_status, "analysis_status": analysis}


def _all_statistics(context: dict) -> dict:
    """要求が挙げた全統計の結果（runtimeに無いものは記録から復元する）。"""
    results: dict = {}
    for spec in (context["request"].get("statistics") or []):
        statistic_id = spec.get("statistic_id")
        result = _statistic_result(context, statistic_id)
        if result is not None:
            results[statistic_id] = result
    return results


def _handle_report(context: dict) -> dict:
    """3軸（実行・QC・解析）を分けて書く。1つの成否へ畳まない（spec §11）。"""
    axes = _axis_statuses(context)
    statistics = _all_statistics(context)
    lines = [
        "# LC–MS metabolomics pipeline report (v2)", "",
        f"- qc_status: {axes['qc_status']}",
        f"- analysis_status: {axes['analysis_status']}",
        "",
        "QC と解析は独立した軸です。QC が fail でも統計が計算できる場合は"
        "診断目的で続行し、結果には QC 不合格を付記します。", "",
        "## Statistics", "",
    ]
    for statistic_id, result in sorted(statistics.items()):
        lines.append(f"- {statistic_id} ({result.get('kind')}): "
                     f"{result.get('status')}"
                     + (f" — {result.get('reason')}" if result.get("reason") else ""))
    if not statistics:
        lines.append("- (no statistic result)")

    path = Path(context["pipeline_root"]) / "quality-report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    ref = report_mod.persist_result(context["pipeline_root"], {
        "output_name": "quality_report", "kind": "quality_report",
        "result_id": f"res_report_{context['pipeline_id']}", "path": str(path)})
    return _ok([ref], record_updates={"analysis": axes})


# ---------- 復元 ----------

def restore_results(record: dict, root: Path, ds) -> dict:
    """再開時に、current な成果物を hash 照合してから runtime へ戻す。

    再計算した値で古い result ID を上書きしない（`results` は追記専用）。hash が
    合わなければ復元せずに止める——「途中まで正しい」状態で先へ進むと、どの
    工程の数字がどの入力から出たのかが言えなくなる。
    """
    root = Path(root)
    restored: dict = {}
    for output_name, runtime_key in _RESTORABLE.items():
        matches = [r for r in (record.get("results") or [])
                   if isinstance(r, dict) and r.get("output_name") == output_name]
        if not matches:
            continue
        ref = matches[-1]
        failures = store.verify_result_refs(root, [ref])
        if failures:
            raise DomainError(
                "RESULT_INTEGRITY_MISMATCH",
                f"保存済みの成果物が変更されています: {output_name}",
                {"output_name": output_name, "failures": failures})
        payload = json.loads((root / ref["relative_path"]).read_text(encoding="utf-8"))
        restored[runtime_key] = payload.get("data")
    return restored


# ---------- 組み立て ----------

def build_handlers() -> dict:
    """v2 固有の handler（stage 名が v1 と重ならないものと、重なるものの v2 側）。

    `prepare_inputs` / `execute_console` / `validate_outputs` / `load_dataset` /
    `resolve_metadata` はここには無い——v1 の実装をそのまま使う工程で、
    `service.build_handlers` が同じ関数を両方の stage 名へ割り当てる。
    """
    return {
        "load_assay_evidence": _handle_load_assay_evidence,
        "resolve_feature_bindings": _handle_resolve_feature_bindings,
        "qc_raw": _handle_qc_raw,
        "preprocess": _handle_preprocess,
        "qc_processed": _handle_qc_processed,
        "statistics": _handle_statistics,
        "export": _handle_export,
        "report": _handle_report,
    }
