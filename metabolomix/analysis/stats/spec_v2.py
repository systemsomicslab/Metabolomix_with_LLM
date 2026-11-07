"""analysis-spec.v2 の正規化と実行（統計エンジン spec §8・§12）。

宣言的な入口: 前処理を 1 回だけ回し、spec の順に compare / composition / pca を実行して
**全部**を載せて返す（選んだものだけにしない）。`primary` は報告の主役を指すだけ。
丸めはここでは行わない（配線層の仕事）。v1（analysis_spec）とは独立で、v1 の schema は拒否する。
"""
from __future__ import annotations

import copy
import math
import re

from metabolomix.analysis.stats.compare import BH_POPULATIONS, EFFECT_SIZES, TESTS, compare
from metabolomix.analysis.stats.composition import composition
from metabolomix.analysis.stats.data import from_dataset
from metabolomix.analysis.stats.pca import pca
from metabolomix.analysis.stats.preprocess import (NeedsInput, normalize_preprocess_config,
                                                   preprocess)
from metabolomix.core.atomic_io import DomainError, canonical_hash

ANALYSIS_SCHEMA_V2 = "analysis-spec.v2"
RESULT_SCHEMA_V2 = "analysis-result.v2"
_INVALID = "ANALYSIS_SPEC_INVALID"
_ID_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")
_KINDS = ("compare", "composition", "pca")
_PCA_SCALING = ("autoscale", "pareto")
_PCA_IMPUTE = ("half_min", "none")
_VERDICT_DEFAULT = {"q": 0.05, "equivalence_fc": 1.2}


def _fail(message: str, **details) -> DomainError:
    return DomainError(_INVALID, message, details)


def _keys(value, allowed: set, required: set, path: str) -> dict:
    if not isinstance(value, dict):
        raise _fail(f"{path} は辞書で指定する", field=path)
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise _fail(f"{path} の未知のキー: {unknown}", field=path, unknown=unknown)
    lacking = sorted(required - set(value))
    if lacking:
        raise _fail(f"{path} に必須のキーが無い: {lacking}", field=path, missing=lacking)
    return value


def _ident(value, path: str) -> str:
    if not isinstance(value, str) or not _ID_RE.match(value):
        raise _fail(f"{path} は英数字と _ . - の文字列", field=path, value=repr(value))
    return value


def _text(value, path: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _fail(f"{path} は空でない文字列で指定する", field=path, value=repr(value))
    return value


def _number(value, path: str, ok) -> float:
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or not ok(value)):
        raise _fail(f"{path} の値が範囲の外", field=path, value=repr(value))
    return float(value)


def _design(raw) -> dict:
    d = _keys({} if raw is None else raw, {"factors", "covariates", "block"}, set(), "design")
    factors = d.get("factors", ["group"])
    if (not isinstance(factors, list) or not factors
            or any(not isinstance(f, str) or not f for f in factors)
            or len(set(factors)) != len(factors)):
        raise _fail("design.factors は重複のない文字列の 1 つ以上のリスト", field="design.factors")
    covs = d.get("covariates", [])
    if not isinstance(covs, list):
        raise _fail("design.covariates はリスト", field="design.covariates")
    out_cov = []
    for i, c in enumerate(covs):
        path = f"design.covariates[{i}]"
        c = _keys(c, {"name", "type", "reference"}, {"name", "type"}, path)
        if c["type"] not in ("numeric", "categorical"):
            raise _fail(f"{path}.type は numeric / categorical", field=path + ".type")
        item = {"name": _text(c["name"], path + ".name"), "type": c["type"]}
        if "reference" in c:
            item["reference"] = _text(c["reference"], path + ".reference")
        out_cov.append(item)
    block = d.get("block")
    if block is not None:
        block = _text(block, "design.block")
    return {"factors": list(factors), "covariates": out_cov, "block": block}


def _comparison(raw, i: int, factors: list[str]) -> dict:
    path = f"comparisons[{i}]"
    if isinstance(raw, dict) and raw.get("kind") == "interaction":
        c = _keys(raw, {"id", "kind", "factors", "levels"}, {"id", "kind", "factors", "levels"}, path)
        fa, lv = c["factors"], c["levels"]
        if not (isinstance(fa, list) and len(fa) == 2 and all(isinstance(f, str) for f in fa)):
            raise _fail(f"{path}.factors は 2 因子のリスト", field=path + ".factors")
        if any(f not in factors for f in fa):
            raise _fail(f"{path}.factors が design.factors に無い", field=path + ".factors")
        if not isinstance(lv, dict) or set(lv) != set(fa) or any(
                not (isinstance(v, list) and len(v) == 2 and all(isinstance(x, str) for x in v))
                for v in lv.values()):
            raise _fail(f"{path}.levels は各因子の [基準, 比較]", field=path + ".levels")
        return {"id": _ident(c["id"], path + ".id"), "kind": "interaction", "factors": list(fa),
                "levels": {f: list(lv[f]) for f in fa}}
    c = _keys(raw, {"id", "factor", "reference", "test", "at"}, {"id", "factor", "reference", "test"},
              path)
    if c["factor"] not in factors:
        raise _fail(f"{path}.factor が design.factors に無い", field=path + ".factor",
                    value=repr(c["factor"]))
    out = {"id": _ident(c["id"], path + ".id"), "factor": c["factor"],
           "reference": _text(c["reference"], path + ".reference"),
           "test": _text(c["test"], path + ".test")}
    if "at" in c:
        at = c["at"]
        if not isinstance(at, dict) or any(not isinstance(k, str) or not isinstance(v, str)
                                           for k, v in at.items()):
            raise _fail(f"{path}.at は 因子 → 水準 の辞書", field=path + ".at")
        out["at"] = dict(sorted(at.items()))
    return out


def _analysis(raw, i: int) -> dict:
    path = f"analyses[{i}]"
    if not isinstance(raw, dict):
        raise _fail(f"{path} は辞書で指定する", field=path)
    kind = raw.get("kind")
    if kind not in _KINDS:
        raise _fail(f"{path}.kind は {list(_KINDS)} のどれか", field=path + ".kind", value=repr(kind))
    aid = _ident(raw.get("analysis_id"), path + ".analysis_id")
    if kind == "pca":
        a = _keys(raw, {"analysis_id", "kind", "scaling", "n_components", "impute"},
                  {"analysis_id", "kind"}, path)
        scaling, impute = a.get("scaling", "autoscale"), a.get("impute", "half_min")
        if scaling not in _PCA_SCALING:
            raise _fail(f"{path}.scaling は {list(_PCA_SCALING)}", field=path + ".scaling")
        if impute not in _PCA_IMPUTE:
            raise _fail(f"{path}.impute は {list(_PCA_IMPUTE)}", field=path + ".impute")
        n = a.get("n_components", 5)
        if isinstance(n, bool) or not isinstance(n, int) or n < 1:
            raise _fail(f"{path}.n_components は 1 以上の整数", field=path + ".n_components")
        return {"analysis_id": aid, "kind": kind, "scaling": scaling, "n_components": n,
                "impute": impute}
    a = _keys(raw, {"analysis_id", "kind", "test", "effect_size", "bh_population",
                    "reported_feature_ids", "verdict", "min_detected_fraction",
                    "allow_confounded"}, {"analysis_id", "kind"}, path)
    test, eff = a.get("test", "moderated_t"), a.get("effect_size", "log2_mean_difference")
    bh = a.get("bh_population", "matrix")
    for name, val, allowed in (("test", test, TESTS), ("effect_size", eff, EFFECT_SIZES),
                               ("bh_population", bh, BH_POPULATIONS)):
        if val not in allowed:
            raise _fail(f"{path}.{name} は {list(allowed)} のどれか", field=f"{path}.{name}",
                        value=repr(val))
    ids = a.get("reported_feature_ids")
    if bh == "reported":
        if (not isinstance(ids, list) or not ids
                or any(not isinstance(x, str) for x in ids)):
            raise _fail(f"{path}: bh_population=reported には reported_feature_ids（文字列の"
                        "空でないリスト）が要る", field=path + ".reported_feature_ids")
    elif ids is not None:
        raise _fail(f"{path}: reported_feature_ids は bh_population=reported のときだけ",
                    field=path + ".reported_feature_ids")
    verdict = _keys(a.get("verdict", {}), {"q", "equivalence_fc"}, set(), path + ".verdict")
    q = _number(verdict.get("q", _VERDICT_DEFAULT["q"]), path + ".verdict.q",
                lambda v: 0 < v <= 1)
    fc = _number(verdict.get("equivalence_fc", _VERDICT_DEFAULT["equivalence_fc"]),
                 path + ".verdict.equivalence_fc", lambda v: v > 1)
    mdf = _number(a.get("min_detected_fraction", 0.5), path + ".min_detected_fraction",
                  lambda v: 0 <= v <= 1)
    conf = a.get("allow_confounded", False)
    if not isinstance(conf, bool):
        raise _fail(f"{path}.allow_confounded は真偽値", field=path + ".allow_confounded")
    out = {"analysis_id": aid, "kind": kind, "test": test, "effect_size": eff,
           "bh_population": bh, "verdict": {"q": q, "equivalence_fc": fc},
           "min_detected_fraction": mdf, "allow_confounded": conf}
    if ids is not None:
        out["reported_feature_ids"] = list(ids)
    return out


def normalize_spec_v2(spec) -> dict:
    """spec を検証して既定を埋めた正準の dict にする。未知のキーは ANALYSIS_SPEC_INVALID。"""
    spec = _keys(spec, {"schema", "input", "preprocess", "design", "comparisons", "analyses",
                        "primary"}, {"schema", "comparisons", "analyses"}, "spec")
    if spec["schema"] != ANALYSIS_SCHEMA_V2:
        raise _fail(f"schema は {ANALYSIS_SCHEMA_V2}（v1 は analysis_spec.run_spec の担当）",
                    field="schema", value=repr(spec["schema"]))
    inp = _keys(spec.get("input", {"source": "dataset"}), {"source"}, set(), "input")
    if inp.get("source", "dataset") != "dataset":
        raise _fail("input.source は dataset だけ", field="input.source")
    try:
        pre = normalize_preprocess_config(spec.get("preprocess"))
    except DomainError as e:
        raise _fail(f"preprocess が不正: {e.message}", field="preprocess",
                    details=e.details) from e
    design = _design(spec.get("design"))
    comps = spec["comparisons"]
    if not isinstance(comps, list) or not comps:
        raise _fail("comparisons は 1 つ以上のリスト", field="comparisons")
    comparisons = [_comparison(c, i, design["factors"]) for i, c in enumerate(comps)]
    cids = [c["id"] for c in comparisons]
    if len(set(cids)) != len(cids):
        raise _fail(f"comparisons の id が重複している: {cids}", field="comparisons")
    alist = spec["analyses"]
    if not isinstance(alist, list) or not alist:
        raise _fail("analyses は 1 つ以上のリスト", field="analyses")
    analyses = [_analysis(a, i) for i, a in enumerate(alist)]
    aids = [a["analysis_id"] for a in analyses]
    if len(set(aids)) != len(aids):
        raise _fail(f"analysis_id が重複している: {aids}", field="analyses")
    stat_ids = [a["analysis_id"] for a in analyses if a["kind"] != "pca"]
    primary = spec.get("primary")
    if primary is None:
        if len(stat_ids) >= 2:
            raise _fail("compare / composition の解析が 2 つ以上あるので primary が要る",
                        field="primary", candidates=stat_ids)
        primary = stat_ids[0] if stat_ids else None
    elif primary not in stat_ids:
        raise _fail("primary は compare / composition の analysis_id", field="primary",
                    value=repr(primary), candidates=stat_ids)
    return {"schema": ANALYSIS_SCHEMA_V2, "input": {"source": "dataset"}, "preprocess": pre,
            "design": design, "comparisons": comparisons, "analyses": analyses,
            "primary": primary}


def rule_hash_v2(normalized: dict) -> str:
    """正規化済み spec の hash（実行のたびに変わる値を含まない）。"""
    return canonical_hash(normalized)


def _sample_table_hash(study) -> str:
    rows = [{"id": s.sample_id, "role": s.role, "group": s.group, "batch": s.batch,
             "injection_order": s.injection_order, "subject_id": s.subject_id,
             "include": s.include, "attributes": dict(sorted(s.attributes.items()))}
            for s in study.samples]
    return canonical_hash(rows)


def run_spec_v2(normalized: dict, ds, *, assay_kind: str = "unknown") -> dict:
    """正規化済み spec を実行する。前処理は 1 回で、全解析が同じ行列を読む。

    前処理が利用者への確認を要するときは `{"status": "needs_input", ...}` を返す。
    """
    spec = copy.deepcopy(normalized)
    rule_hash = rule_hash_v2(spec)
    study = from_dataset(ds)
    prepared = preprocess(study, spec["preprocess"])
    if isinstance(prepared, NeedsInput):
        return {**prepared.to_envelope(), "spec": spec, "rule_hash": rule_hash}
    design = spec["design"]
    common = dict(factors=design["factors"], comparisons=spec["comparisons"])
    results = []
    for a in spec["analyses"]:
        kind = a["kind"]
        if kind == "pca":
            res = pca(prepared, scaling=a["scaling"], n_components=a["n_components"],
                      impute=a["impute"])
        else:
            kw = dict(common, covariates=design["covariates"], block=design["block"],
                      test=a["test"], effect_size=a["effect_size"],
                      bh_population=a["bh_population"],
                      reported_feature_ids=a.get("reported_feature_ids"),
                      q_threshold=a["verdict"]["q"], equivalence_fc=a["verdict"]["equivalence_fc"],
                      min_detected_fraction=a["min_detected_fraction"],
                      allow_confounded=a["allow_confounded"])
            res = (compare(prepared, **kw) if kind == "compare"
                   else composition(prepared, assay_kind=assay_kind, **kw))
        results.append({"analysis_id": a["analysis_id"], "kind": kind, "result": res})
    return {"schema": RESULT_SCHEMA_V2, "spec": spec, "rule_hash": rule_hash,
            "input_fingerprint": {"dataset_hash": prepared.dataset_hash,
                                  "sample_table_hash": _sample_table_hash(study)},
            "matrix_id": prepared.matrix_id,
            "preprocess": {"steps": prepared.steps, "warnings": prepared.warnings,
                           "subject_source": prepared.subject_source},
            "primary": spec["primary"], "analyses": results}
