"""m/z 帯ごとの質量精度と、新しい知見の候補の一覧（spec 2026-10-22 §4 の 9(c)・11）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

- 質量精度: SME（MS/MS の照合で付いた同定）の実測 m/z と理論 m/z の差（ppm）の分布。
  許容幅の決定はしない（再解析の規則で決める）。四分位は `statistics.quantiles(n=4)`
  （exclusive 法）。`subset="msms_matched"` は、その特徴の MS/MS が参照スペクトルに照合された
  同定（SML の `opt_global_spectrum_matched` が真。研究室 MSP と LBM は区別しない）に絞る。
  注入順の傾向（`injection_trend`）は、試料ごとの実測 m/z（`ds.per_sample_mz`）の誤差の中央値と
  注入順の Spearman 相関。
- 新しい知見の候補: SME の同定・m/z 帯ごとの許容幅・q の閾値で絞り、除く対象を外した
  特徴の一覧。統計は 2 群比較の結果（`result_id`）か analysis-spec の解析から引く。
  割れたアライメントスポットの組（`alignment_qc` の `merge_groups`）は 1 つの物質として
  1 候補に数える（spec 2026-10-28 §4.3 (c)。`novel_candidates` の `merged_from`）。
"""
from __future__ import annotations

import math
import re
import statistics

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["NOVEL_SCHEMA", "SUBSETS", "band_summary", "has_match_information", "injection_trend",
           "normalize_novel_spec", "novel_candidates", "sme_ppm", "summarize"]

NOVEL_SCHEMA = "novel-candidates.v1"
#: `sme_ppm` の絞り込み。rank1 = rank の SME すべて、msms_matched = そのうち MS/MS が照合されたもの。
SUBSETS = ("rank1", "msms_matched")
#: 注入順との Spearman 相関を出す最小の試料数（これ未満は rho を None にする）。
MIN_INJECTIONS_FOR_RHO = 5
_INVALID = "NOVEL_SPEC_INVALID"


def _number(value) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    return float(value)


def _spectrum_matched(ds, fid) -> bool | None:
    return ((getattr(ds, "feature_annotations", None) or {}).get(fid) or {}).get("spectrum_matched")


def has_match_information(ds) -> bool:
    """SML の `opt_global_spectrum_matched` を読めた特徴が 1 つでもあるか（列が無い mzTab は偽）。"""
    return any(_spectrum_matched(ds, fid) is not None for fid in ds.feature_ids)


def sme_ppm(ds, rank: int = 1, subset: str = "rank1") -> list[dict]:
    """特徴ごとに、その rank の SME 候補の (実測 − 理論)/理論 × 1e6。並びは ds.feature_ids。

    `subset="msms_matched"` は、その特徴の SML の `opt_global_spectrum_matched` が真のものだけ
    （偽・不明の特徴は入れない。LBM 由来と MSP 由来は区別しない）。
    """
    if subset not in SUBSETS:
        raise ValueError(f"subset は {SUBSETS}")
    out = []
    for fid in ds.feature_ids:
        if subset == "msms_matched" and _spectrum_matched(ds, fid) is not True:
            continue
        for candidate in (ds.feature_candidates or {}).get(fid) or []:
            if candidate.get("rank") != rank:
                continue
            exp = _number(candidate.get("exp_mass_to_charge"))
            theo = _number(candidate.get("theoretical_mass_to_charge"))
            if exp is None or theo is None or theo <= 0:
                continue
            out.append({"feature_id": fid, "theoretical_mz": theo, "mz": exp,
                        "ppm": (exp - theo) / theo * 1e6, "name": candidate.get("chemical_name"),
                        "inchikey": candidate.get("inchikey"), "adduct": candidate.get("adduct")})
            break
    return out


def summarize(values: list[float]) -> dict:
    values = sorted(values)
    if len(values) < 2:
        return {"n": len(values)}
    q1, _, q3 = statistics.quantiles(values, n=4)
    return {"n": len(values), "median": statistics.median(values), "q1": q1, "q3": q3,
            "within_5ppm": sum(abs(v) <= 5 for v in values) / len(values),
            "within_10ppm": sum(abs(v) <= 10 for v in values) / len(values)}


def injection_trend(ds, rows: list[dict]) -> dict:
    """試料ごとの質量誤差（ppm）の中央値と、注入順との Spearman 相関。

    `rows` は `sme_ppm` の行（選んだ subset のものだけ）。全試料を数える（ブランク・QC が
    あれば含む。使える特徴が 1 つも無い試料は載せない）。試料の誤差は、その試料で**検出された**（gap-fill でない）セルの
    実測 m/z（`ds.per_sample_mz`）と理論 m/z の差。注入順は割れの組の判定と同じ規則
    （`alignment_qc.sample_injection_order`。出所を `injection_order_source` に書く）。
    試料ごとの m/z が無ければ `{"status": "unavailable"}`。rho は、注入順の出所が
    `filename`（採取時刻が読めない）か、データのある試料が `MIN_INJECTIONS_FOR_RHO` 未満なら
    None で、理由を `rho_unavailable_reason`（`filename_order` / `too_few_samples`）に書く。
    """
    per_sample = getattr(ds, "per_sample_mz", None)
    if per_sample is None:
        return {"status": "unavailable"}
    from metabolomix.analysis.alignment_qc import order_label, sample_injection_order

    names = [str(n) for n in ds.sample_names]
    order, source = sample_injection_order(ds)
    errors: dict[str, list[float]] = {name: [] for name in names}
    for row in rows:
        theo = row["theoretical_mz"]
        for name, mz in (per_sample.get(row["feature_id"]) or {}).items():
            if name in errors:
                errors[name].append((mz - theo) / theo * 1e6)
    per_injection = [{"sample": names[i], "injection_order": order_label(order[i]),
                      "median_ppm": statistics.median(errors[names[i]]), "n": len(errors[names[i]])}
                     for i in sorted(range(len(names)), key=lambda i: (order[i], i))
                     if errors[names[i]]]
    rho = None
    reason = None
    if source == "filename":
        # 採取時刻（注入順）が読めずファイル名の順に落ちたときは、順位相関を出さない
        # （ファイル名の順は注入順とは限らない。設計 4.3(c)）。試料ごとの中央値は残す
        reason = "filename_order"
    elif len(per_injection) < MIN_INJECTIONS_FOR_RHO:
        reason = "too_few_samples"
    else:
        from scipy.stats import spearmanr

        value = float(spearmanr([p["injection_order"] for p in per_injection],
                                [p["median_ppm"] for p in per_injection]).statistic)
        rho = value if math.isfinite(value) else None      # 片方が定数なら相関は定義できない
    return {"per_injection": per_injection, "spearman_rho": rho,
            **({"rho_unavailable_reason": reason} if reason else {}),
            "n_samples": len(per_injection), "injection_order_source": source}


def band_summary(rows: list[dict], band_edges: list[float]) -> dict:
    edges = sorted(band_edges)
    bounds = list(zip([None, *edges], [*edges, None]))
    bands = []
    for low, high in bounds:
        chosen = [r["ppm"] for r in rows
                  if (low is None or r["theoretical_mz"] >= low)
                  and (high is None or r["theoretical_mz"] < high)]
        bands.append({"min_mz": low, "max_mz": high, **summarize(chosen)})
    return {"overall": summarize([r["ppm"] for r in rows]), "bands": bands}


# ---------- 新しい知見の候補 ----------

def normalize_novel_spec(spec) -> dict:
    if not isinstance(spec, dict) or set(spec) - {"schema", "statistics", "ppm", "q_max", "exclude",
                                                         "labels", "ppm_offset"}:
        raise DomainError(_INVALID, "spec は {statistics, ppm, q_max, exclude, labels, ppm_offset}",
                          {})
    stats_spec = spec.get("statistics")
    if not isinstance(stats_spec, dict):
        raise DomainError(_INVALID, "statistics は {result_id} か {analysis_spec, analysis_id}",
                          {"field": "statistics"})
    if "result_id" in stats_spec and not set(stats_spec) - {"result_id", "merge_split",
                                                            "alignment_qc"}:
        # 割れたスポットのまとめ（既定 true）。analysis_spec の側は features.merge_split で決める。
        merge = stats_spec.get("merge_split", True)
        if not isinstance(merge, bool):
            raise DomainError(_INVALID, "statistics.merge_split は true / false",
                              {"field": "statistics.merge_split"})
        statistics_out = {"result_id": str(stats_spec["result_id"])}
        if merge:
            # 真（既定）のときだけ規則に載せる。false は v7 と同じ規則・同じ rule_hash になる
            statistics_out["merge_split"] = True
        if "alignment_qc" in stats_spec:
            if not merge or not isinstance(stats_spec["alignment_qc"], str):
                raise DomainError(_INVALID, "statistics.alignment_qc は merge_split=true のときの "
                                  "JSON のパス", {"field": "statistics.alignment_qc"})
            statistics_out["alignment_qc"] = stats_spec["alignment_qc"]
    elif set(stats_spec) == {"analysis_spec", "analysis_id"}:
        from metabolomix.analysis import analysis_spec
        statistics_out = {"analysis_spec": analysis_spec.normalize_spec(stats_spec["analysis_spec"]),
                          "analysis_id": str(stats_spec["analysis_id"])}
        ids = [a["analysis_id"] for a in statistics_out["analysis_spec"]["analyses"]]
        if statistics_out["analysis_id"] not in ids:
            raise DomainError(_INVALID, f"analysis_id が spec に無い（{ids}）",
                              {"field": "statistics.analysis_id"})
    else:
        raise DomainError(_INVALID, "statistics は {result_id} か {analysis_spec, analysis_id}",
                          {"field": "statistics"})
    ppm = spec.get("ppm", {"default": 10.0, "bands": []})
    if (not isinstance(ppm, dict) or set(ppm) - {"default", "bands"}
            or _number(ppm.get("default", 10.0)) is None):
        raise DomainError(_INVALID, "ppm は {default, bands: [{max_mz, ppm}]}", {"field": "ppm"})
    bands = []
    for i, band in enumerate(ppm.get("bands") or []):
        if (not isinstance(band, dict) or set(band) != {"max_mz", "ppm"}
                or _number(band["max_mz"]) is None or _number(band["ppm"]) is None):
            raise DomainError(_INVALID, f"ppm.bands[{i}] は {{max_mz, ppm}}", {"field": "ppm"})
        bands.append({"max_mz": float(band["max_mz"]), "ppm": float(band["ppm"])})
    q_max = _number(spec.get("q_max", 0.05))
    if q_max is None or not 0 < q_max <= 1:
        raise DomainError(_INVALID, "q_max は 0 < q_max ≤ 1", {"field": "q_max"})
    exclude = spec.get("exclude") or {}
    if not isinstance(exclude, dict) or set(exclude) - {"feature_ids", "inchikey_blocks"}:
        raise DomainError(_INVALID, "exclude は {feature_ids, inchikey_blocks}", {"field": "exclude"})
    normalized = {"schema": NOVEL_SCHEMA, "statistics": statistics_out,
                  "ppm": {"default": float(ppm.get("default", 10.0)),
                          "bands": sorted(bands, key=lambda b: b["max_mz"])},
                  "q_max": q_max,
                  "exclude": {"feature_ids": sorted(set(map(str, exclude.get("feature_ids") or []))),
                              "inchikey_blocks": sorted({str(b).strip()[:14].upper()
                                                         for b in exclude.get("inchikey_blocks") or []})}}
    if "labels" in spec:
        normalized["labels"] = _normalize_labels(spec["labels"])
    if spec.get("ppm_offset") is not None:
        # 渡したときだけ規則に載せる（省けば v8 と同じ spec・同じ rule_hash）
        normalized["ppm_offset"] = normalize_ppm_offset(spec["ppm_offset"], code=_INVALID)
    return normalized


def normalize_ppm_offset(raw, *, code: str, field: str = "ppm_offset") -> dict | None:
    """質量の窓の中心のずれ `{"default": 0.0, "bands": [{"max_mz", "offset"}]}` を正規化する。

    中心は「実測 − 理論」の ppm の符号つきの値。帯の選び方は ppm の帯と同じ（理論 m/z が
    `max_mz` 未満の最初の帯、無ければ `default`）。None なら None。
    """
    if raw is None:
        return None
    detail = {"field": field}
    if (not isinstance(raw, dict) or set(raw) - {"default", "bands"}
            or _number(raw.get("default", 0.0)) is None):
        raise DomainError(code, f"{field} は {{default, bands: [{{max_mz, offset}}]}}", detail)
    bands = []
    for i, band in enumerate(raw.get("bands") or []):
        if (not isinstance(band, dict) or set(band) != {"max_mz", "offset"}
                or _number(band["max_mz"]) is None or _number(band["offset"]) is None):
            raise DomainError(code, f"{field}.bands[{i}] は {{max_mz, offset}}", detail)
        bands.append({"max_mz": float(band["max_mz"]), "offset": float(band["offset"])})
    return {"default": float(raw.get("default", 0.0)),
            "bands": sorted(bands, key=lambda b: b["max_mz"])}


def offset_for(mz: float, offset_spec: dict | None) -> float:
    """理論 m/z の帯の中心のずれ（ppm）。帯の境目の m/z は次の帯（`_tolerance` と同じ `<`）。"""
    if offset_spec is None:
        return 0.0
    for band in offset_spec["bands"]:
        if mz < band["max_mz"]:
            return band["offset"]
    return offset_spec["default"]


#: 外部の印の既定（強い順）
LABEL_ORDER = ["A", "B", "C", "D"]
_LABEL_NAME = re.compile(r"^[A-Za-z0-9_]+$")


def _normalize_labels(labels) -> dict:
    """印の spec。`path` と `name` は必須。印は選び方を変えないので `rule_hash` には入れない。"""
    allowed = {"path", "name", "feature_column", "label_column", "order"}
    if (not isinstance(labels, dict) or set(labels) - allowed
            or not isinstance(labels.get("path"), str) or not labels["path"]
            or not isinstance(labels.get("name"), str)
            or not _LABEL_NAME.match(labels["name"])):
        raise DomainError(_INVALID, "labels は {path, name（英数字と _）, feature_column, "
                          "label_column, order}", {"field": "labels"})
    order = labels.get("order", LABEL_ORDER)
    columns = [labels.get("feature_column", "feature_id"), labels.get("label_column", "level")]
    if (not isinstance(order, list) or not order
            or not all(isinstance(x, str) and x for x in order) or len(set(order)) != len(order)
            or not all(isinstance(c, str) and c for c in columns)):
        raise DomainError(_INVALID, "labels.order は重ならない水準の並び（強い順）、列名は文字列",
                          {"field": "labels"})
    return {"path": labels["path"], "name": labels["name"], "feature_column": columns[0],
            "label_column": columns[1], "order": list(order)}


def read_labels(labels: dict) -> dict:
    """印の TSV を読む。特徴ごとに最も強い水準を残す。空の水準は無視、`order` に無い水準は拒む。"""
    from metabolomix.analysis import analysis_spec

    columns, rows, sha256 = analysis_spec._read_tsv(labels["path"], "labels")
    analysis_spec._require_column(columns, labels["feature_column"], "labels")
    analysis_spec._require_column(columns, labels["label_column"], "labels")
    rank = {level: i for i, level in enumerate(labels["order"])}
    by_feature: dict[str, str] = {}
    unknown: set[str] = set()
    for row in rows:
        level = (row.get(labels["label_column"]) or "").strip()
        feature = (row.get(labels["feature_column"]) or "").strip()
        if not level or not feature:
            continue
        if level not in rank:
            unknown.add(level)
        elif feature not in by_feature or rank[level] < rank[by_feature[feature]]:
            by_feature[feature] = level
    if unknown:
        raise DomainError(_INVALID, f"labels に order に無い水準がある: {sorted(unknown)}",
                          {"field": "labels", "unknown_levels": sorted(unknown)})
    return {"by_feature": by_feature, "sha256": sha256, "n_rows": len(rows)}


def attach_labels(rows: list[dict], by_feature: dict[str, str], *, name: str, order: list[str],
                  top_n: int = 20) -> dict:
    """各行に `labels.<name>`（水準か None）を足し、水準ごとの件数を返す。並びは変えない。
    割れの組（`merged_from`）は群のどれかの最も強い水準。"""
    rank = {level: i for i, level in enumerate(order)}
    column = f"labels.{name}"
    for row in rows:
        group = row.get("merged_from") or [row["feature_id"]]
        found = [by_feature[f] for f in group if f in by_feature]
        row[column] = min(found, key=rank.__getitem__) if found else None

    def count(part: list[dict]) -> dict:
        out = {level: 0 for level in order}
        out["none"] = 0
        for row in part:
            out[row[column] or "none"] += 1
        return out

    return {"all": count(rows), "top20": count(rows[:top_n])}


def rule_hash(normalized: dict, merge_groups_hash: str | None = None) -> str:
    """規則の hash（パスの欄を除く）。まとめた割れの組の hash（`merge_groups_hash`）を渡したときだけ
    規則に足す（同じ spec でも組が変われば別の規則。`analysis_spec.rule_hash` と同じ扱い）。
    まとめない（None）なら v7 と同じ hash。"""
    from metabolomix.analysis.analysis_spec import _strip_paths

    material = _strip_paths(normalized)
    material.pop("labels", None)          # 印は選び方を変えない
    if merge_groups_hash is not None:
        material["merged_groups_hash"] = merge_groups_hash
    return canonical_hash(material)


def _tolerance(mz: float, ppm_spec: dict) -> float:
    for band in ppm_spec["bands"]:
        if mz < band["max_mz"]:
            return band["ppm"]
    return ppm_spec["default"]


def novel_candidates(ds, normalized: dict, stats_rows: dict, *, group_ref, group_test,
                     merged_from: dict[str, list[str]] | None = None) -> list[dict]:
    """`stats_rows` は feature_id → {log2fc_test_over_ref, p, q}（検定できた行）。統計エンジンの結果の行は
    `verdict` も持ち、q < `q_max` に加えて `verdict == "changed"` の行だけを候補にする（R16: 検出が足りず
    量を判定できない `inconclusive` の行は q が小さくても通さない）。v1 の行（`verdict` 無し）は q だけ。

    `merged_from` は割れたスポットをまとめた群（残す ID → 群。残す ID が先頭）。None（まとめない）
    なら行は v7 と同じ形（`merged_from`・`sme_feature_id` を足さない）。群は 1 候補に
    数え、`feature_id` は残す ID・統計はその ID の行（まとめた行列での検定）。同定と質量誤差
    （`name` `inchikey` `adduct` `mz` `ppm` `rt`）は SME（rank 1）を持つ方のスポット（両方に
    あれば残す ID）から取り、その ID を `sme_feature_id` に書く。`exclude.feature_ids` は群の
    どちらかの ID に当たれば群ごと外す（同じ物質なので）。
    """
    excluded_ids = set(normalized["exclude"]["feature_ids"])
    excluded_blocks = set(normalized["exclude"]["inchikey_blocks"])
    merging = merged_from is not None
    merged_from = merged_from or {}
    absorbed = {fid for group in merged_from.values() for fid in group[1:]}
    sme = {row["feature_id"]: row for row in sme_ppm(ds, rank=1)}
    out = []
    for fid in ds.feature_ids:
        if fid in absorbed:
            continue
        group = merged_from.get(fid) or [fid]
        row = next((sme[f] for f in group if f in sme), None)
        if row is None:
            continue
        stats = stats_rows.get(fid)
        if stats is None or stats.get("q") is None or not stats["q"] < normalized["q_max"]:
            continue
        if stats.get("verdict", "changed") != "changed":
            continue                    # 統計エンジンの行は判定が changed のときだけ（R16。v1 の行は判定を持たない）
        center = offset_for(row["theoretical_mz"], normalized.get("ppm_offset"))
        if abs(row["ppm"] - center) > _tolerance(row["theoretical_mz"], normalized["ppm"]):
            continue
        block = (row["inchikey"] or "").strip()[:14].upper()
        if excluded_ids.intersection(group) or (block and block in excluded_blocks):
            continue
        log2fc = stats.get("log2fc_test_over_ref")
        source = row["feature_id"]
        out.append({**row, "feature_id": fid,
                    "rt": (ds.feature_metadata.get(source) or {}).get("rt"),
                    "group_ref": group_ref, "group_test": group_test,
                    "log2fc_test_over_ref": log2fc, "p": stats.get("p"), "q": stats["q"],
                    "direction": None if not log2fc else
                    ("test_higher" if log2fc > 0 else "ref_higher")})
        if merging:
            out[-1].update(sme_feature_id=source,
                           merged_from=list(group) if len(group) > 1 else None)
    out.sort(key=lambda r: (r["q"], r["feature_id"]))
    return out
