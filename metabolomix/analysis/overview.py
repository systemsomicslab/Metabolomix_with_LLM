"""最初の概観の事実の組み立て（spec 2026-11-06-dataset-overview-design、計画 2026-11-09-a-overview）。

品質 → 全体 → クラス → 分子種の 4 段それぞれについて、事実（`{"id", "text", "data"}`）・判定・次の一手 1 行を作る。
`text` は戻り値に載る 1 文（条件付きの数字）、`data` は記録に残す構造。事実の ID は段ごとの通し番号
（`Q1` `G1` `C1` `S1`）で、外からは `ovw_<hash>#Q1` で引く。

計算はしない（統計は C のエンジンが済ませた結果を読む）。例外として、比較の記録は `as_two_group` を通さず
`engine.rows` / `engine.omnibus` を直に読む——概観は 3 群以上も扱い、**件数**を数えるだけで 2 群の形を要らないため。
判定は行の `verdict`（changed / unchanged / inconclusive）だけを使い、q や log2fc から「変わっていない」を作らない。

`tools.*`・`core.session_state`・`core.mcp_core` を import しない（`stats_results` と同じ位置づけ）。
"""
from __future__ import annotations

import json
import math

from metabolomix.analysis.result_state import dataset_fingerprint
from metabolomix.analysis.stats.quality import QUALITY_THRESHOLDS
from metabolomix.analysis.stats_results import contrast_id_for, sample_ids_by_name
from metabolomix.core.atomic_io import canonical_hash

__all__ = ["OVERVIEW_SCHEMA", "MAX_PAIRWISE_GROUPS", "group_plan", "group_comparisons", "chain_state",
           "id_level", "quality_facts", "global_facts", "class_facts", "species_facts",
           "overview_id", "register_overview"]

OVERVIEW_SCHEMA = "lipidmix.overview.v1"
#: 全対を比べる群の数の上限（6 群で 15 対）。超えたら基準との対比と moderated F だけ。
MAX_PAIRWISE_GROUPS = 6
_CLASS_ROWS = 5            # クラスの表に載せる変わった／変わっていないクラスの各上限


# --- 書式 ---------------------------------------------------------------------------------

def _g(x, digits: int = 3) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "NA"
    return format(float(x), f".{digits}g")


def _pct(x) -> str:
    return "NA" if x is None else f"{100 * float(x):.0f}%"


def _j(value) -> str:
    """次の一手の引数の書き方（JSON。日本語の群名もそのまま）。"""
    return json.dumps(value, ensure_ascii=False)


def _cell(text) -> str:
    return str(text if text is not None else "").replace("\t", " ").replace("\n", " ")


def _tsv(columns: list[str], rows: list[list]) -> str:
    return "\n".join(["\t".join(columns), *("\t".join(_cell(v) for v in r) for r in rows)]) if rows else ""


class _Facts:
    """段の事実を通し番号で積む。"""

    def __init__(self, prefix: str):
        self.prefix, self.items = prefix, []

    def add(self, text: str, **data) -> str:
        fid = f"{self.prefix}{len(self.items) + 1}"
        self.items.append({"id": fid, "text": f"{fid} {text}", "data": data})
        return fid


# --- 群の計画 ------------------------------------------------------------------------------

def group_plan(ds, prepared) -> dict:
    """比べる群と仮の向き。群は行列の個体（QC・ブランクは入らない）から、順序は試料の表で先に現れた順。

    `mode`: `none`（群が 1 つ以下）/ `two` / `multi`。`two` の基準は先に現れた群（`direction="provisional"`）。
    """
    position = {sid: k for k, (name, sid) in enumerate(sample_ids_by_name(ds).items())}
    order = {name: k for k, name in enumerate(ds.sample_names)}
    first: dict[str, int] = {}
    for unit in prepared.units:
        if unit.group is None:
            continue
        idx = min((position.get(s, order.get(s, 10 ** 9)) for s in unit.sample_ids), default=10 ** 9)
        first[unit.group] = min(first.get(unit.group, idx), idx)
    levels = sorted(first, key=lambda g: (first[g], g))
    if len(levels) < 2:
        return {"mode": "none", "levels": levels, "reference": None, "direction": None,
                "reason": "群が 1 つ以下（Class ID が空か、すべて同じ）"}
    return {"mode": "two" if len(levels) == 2 else "multi", "levels": levels, "reference": levels[0],
            "direction": "provisional",
            "direction_rule": "試料の表で先に現れた群を基準にした（仮。名前から対照を推測していない）"}


def group_comparisons(plan: dict) -> tuple[list[dict], bool]:
    """spec v2 の comparisons と、全対を省いたか。`two` は 1 対比、`multi` は群が 6 つまで全対
    （基準は先に現れた方）、7 群以上は基準との対比だけ（moderated F は同じ当てはめから出る）。"""
    levels = plan["levels"]
    taken: set[str] = set()

    def comp(ti, ri):
        return {"id": contrast_id_for(levels[ti], levels[ri], taken, test_pos=ti, ref_pos=ri),
                "factor": "group", "reference": levels[ri], "test": levels[ti]}

    if plan["mode"] == "none":
        return [], False
    if plan["mode"] == "two" or len(levels) > MAX_PAIRWISE_GROUPS:
        return [comp(i, 0) for i in range(1, len(levels))], plan["mode"] == "multi"
    return [comp(j, i) for i in range(len(levels)) for j in range(i + 1, len(levels))], False


# --- 状態 ----------------------------------------------------------------------------------

def chain_state(ds, *, curation_judgments: int | None) -> dict:
    """鎖のどの工程が済んでいるか（B §6）。選別まで済んでいなければ暫定。"""
    steps = list(((getattr(ds, "quant_state", None) or {}).get("steps")) or [])
    state = {"curation": curation_judgments,
             "normalize": "normalize" in steps, "merge": "merge" in steps, "select": "select" in steps}
    pending = [name for name in ("normalize", "select") if not state[name]]
    state["provisional"] = not state["select"]
    state["pending"] = pending
    return state


def id_level(ds, fid: str) -> str:
    """同定の証拠の水準: `msms`（MS/MS の候補あり）/ `ms1`（MS1 注釈のみ）/ `none`（load_dataset の要約と同じ 3 分類）。"""
    if (getattr(ds, "feature_candidates", None) or {}).get(fid):
        return "msms"
    if ((getattr(ds, "feature_annotations", None) or {}).get(fid) or {}).get("name"):
        return "ms1"
    return "none"


def _name(ds, fid: str) -> str:
    return str(((getattr(ds, "feature_metadata", None) or {}).get(fid) or {}).get("name") or "")


# --- 品質 ----------------------------------------------------------------------------------

def quality_facts(*, rsd: dict, failed: dict, detection: dict, drift: list[dict], confounding: dict,
                  replicate_warnings: list, verdict: str, reasons: list[str], n_qc: int, n_blank: int,
                  drift_note: str | None, pca_id: str | None) -> dict:
    """品質の段。判定（ok / caution / poor / not_evaluable）と理由のコードを最初の事実に置く。"""
    f = _Facts("Q")
    f.add(f"判定 {verdict}" + (f"（{', '.join(reasons)}）" if reasons else "")
          + f"。QC {n_qc} 本・ブランク {n_blank} 本"
          + ("（ブランクが無いのでブランクの足切りは行っていない）" if not n_blank else ""),
          verdict=verdict, reasons=reasons, n_qc=n_qc, n_blank=n_blank, thresholds=QUALITY_THRESHOLDS)
    if rsd.get("status") == "ok":
        f.add(f"QC の RSD（前処理前）: 中央値 {_g(rsd['median'], 2)}・RSD≤30% の特徴 {_pct(rsd['frac_le_030'])}"
              f"（{rsd['n_computable']} 特徴）", **rsd)
    else:
        f.add(f"QC の RSD は評価できない（{rsd.get('reason')}。QC {rsd.get('n_qc', 0)} 本）", **rsd)
    if failed.get("failed"):
        f.add(f"失敗した QC 注入 {len(failed['failed'])} 本（総量が QC の中央値の 0.2 倍未満）: "
              + ", ".join(failed["failed"][:5]), **failed)
    if detection.get("status") == "ok":
        f.add(f"生体試料のセル: 実検出 {_pct(detection['detected_frac'])}・gap-fill {_pct(detection['gap_filled_frac'])}"
              f"・欠測 {_pct(detection['missing_frac'])}", **detection)
    else:
        f.add(f"実検出と gap-fill の別は評価できない（{detection.get('reason')}）", **detection)
    flagged = [d for d in drift if d.get("flag")]
    if flagged:
        f.add("ドリフトの兆候: " + "; ".join(f"バッチ {d['batch']} で QC の総量と注入順の ρ={_g(d['rho'], 2)}"
                                           f"（QC {d['n_qc']} 本）" for d in flagged[:3]), batches=drift)
    elif drift:
        f.add(f"ドリフトの兆候なし（QC が 5 本以上の {len(drift)} バッチで |ρ|<0.6 か p≥0.05）", batches=drift)
    status = confounding.get("batch_status")
    inj = (confounding.get("injection_order") or {}).get("status")
    conf_text = {"complete": "群とバッチが完全交絡（群差とバッチ差を区別できない。比較は回していない）",
                 "partial": f"群とバッチが部分交絡（Cramér の V={_g(confounding.get('cramers_v'), 2)}）",
                 "none": "群とバッチの交絡なし", "single_batch": "バッチは 1 つ",
                 "not_evaluable": "群とバッチの交絡は評価できない（バッチが不明な試料がある）"}.get(status, str(status))
    f.add(conf_text + ("・注入順が群で偏っている" if inj == "biased" else ""),
          batch_status=status, cramers_v=confounding.get("cramers_v"), injection_order=inj,
          batch_table=confounding.get("batch_table"))
    if replicate_warnings:
        f.add(f"反復注入らしい試料名が {len(replicate_warnings)} 組あるのに subject_id が無い（独立な試料として数えている）",
              stems=[w.get("stem") for w in replicate_warnings])
    if drift_note:
        f.add(drift_note, drift="off_undecided")
    nxt = [f'dataset_plot(kind="pca", result_id="{pca_id}", color_by="batch")'] if pca_id else []
    if failed.get("failed"):
        nxt.append(f"dataset_exclude(samples={_j(failed['failed'][:5])})")
    if drift_note:
        nxt.append('dataset_preprocess(drift_mode="on" か "off")')
    if replicate_warnings:
        nxt.append("dataset_set_sample_metadata（subject_id を入れて反復をまとめる）")
    return {"verdict": verdict, "facts": f.items, "next": " / ".join(nxt) or "なし"}


# --- 全体 ----------------------------------------------------------------------------------

def _eta_squared(scores, labels: list[str]) -> float | None:
    """一元配置の平方和比（群で説明される分散の割合）。"""
    import numpy as np

    x = np.asarray(scores, dtype=float)
    if x.size < 2 or len(set(labels)) < 2:
        return None
    total = float(((x - x.mean()) ** 2).sum())
    if total <= 0:
        return None
    between = sum(len(idx) * (x[idx].mean() - x.mean()) ** 2
                  for idx in ([i for i, g in enumerate(labels) if g == lvl] for lvl in set(labels)))
    return float(between / total)


def _omnibus_counts(record: dict) -> tuple[int, int]:
    q_thr = record["engine"]["verdict_rule"]["q_threshold"]
    tested = [o for o in (record["engine"].get("omnibus") or []) if o.get("f_p") is not None]
    return len(tested), sum(1 for o in tested if o.get("f_q") is not None and o["f_q"] <= q_thr)


def _verdict_counts(rows: list[dict]) -> dict:
    out = {"changed": 0, "up": 0, "down": 0, "unchanged": 0, "inconclusive": 0}
    for r in rows:
        out[r["verdict"]] += 1
        if r["verdict"] == "changed" and r.get("log2fc") is not None:
            out["up" if r["log2fc"] > 0 else "down"] += 1
    return out


def global_facts(*, pca_result: dict | None, labels: dict, record: dict | None, plan: dict,
                 omitted: str | None, quality_note: str | None) -> dict:
    """全体の段: PCA の寄与率・群の重心・群で説明される分散（η²）と、比較の全体の件数。"""
    f = _Facts("G")
    if pca_result is not None:
        evr = list(pca_result.get("explained_variance_ratio") or [])
        units = pca_result["unit_ids"]
        scores = pca_result["scores"]["units"]
        groups = [labels.get(u, "(none)") for u in units]
        n_pc = min(2, len(evr))
        eta = [_eta_squared([row[k] for row in scores], groups) for k in range(n_pc)]
        f.add("PCA（前処理後・log2・autoscale。QC は投影）: "
              + " / ".join(f"PC{k + 1} {_pct(evr[k])}" for k in range(n_pc))
              + "。群で説明される分散 η² " + " / ".join(f"PC{k + 1} {_g(eta[k], 2)}" for k in range(n_pc))
              + f"（個体 {len(units)}・QC {len(pca_result.get('qc_ids') or [])}）",
              explained_variance_ratio=evr[:n_pc], eta_squared=eta, n_units=len(units))
        cents = {}
        for g in dict.fromkeys(groups):
            idx = [i for i, x in enumerate(groups) if x == g]
            cents[g] = [sum(scores[i][k] for i in idx) / len(idx) for k in range(n_pc)]
        f.add("群の重心 (PC1, PC2): " + " / ".join(
            f"{g} ({', '.join(_g(v, 2) for v in c)})" for g, c in list(cents.items())[:8]), centroids=cents)
    if record is None:
        f.add(f"比較は回していない: {omitted}", omitted=omitted)
    elif plan["mode"] == "two":
        rows = record["engine"]["rows"]
        c = _verdict_counts(rows)
        r0 = rows[0] if rows else {}
        f.add(f"{r0.get('test')} vs {r0.get('reference')}（仮の向き）: changed {c['changed']}"
              f"（上昇 {c['up']}・下降 {c['down']}）/ unchanged {c['unchanged']} / inconclusive {c['inconclusive']}"
              + (f"。{quality_note}" if quality_note else ""), **c)
    else:
        n_tested, n_sig = _omnibus_counts(record)
        f.add(f"moderated F（どれかの群で差）: q≤{_g(record['engine']['verdict_rule']['q_threshold'])} の特徴 "
              f"{n_sig} / 検定 {n_tested}" + (f"。{quality_note}" if quality_note else ""),
              n_tested=n_tested, n_q_below=n_sig)
    if plan["mode"] == "two":
        nxt = (f'dataset_compare(reference_group="{plan["levels"][0]}", test_group="{plan["levels"][1]}")'
               "（向きを決めて比べ直す。逆なら入れ替える）")
    elif plan["mode"] == "multi":
        nxt = f"dataset_compare(groups={_j(plan['levels'])})（先頭が基準。対照を決めて並べ替える）"
    else:
        nxt = "dataset_set_sample_metadata（群＝Class ID を入れる）→ dataset_overview()"
    return {"facts": f.items, "next": nxt}


# --- クラス --------------------------------------------------------------------------------

def class_facts(*, composition_result: dict | None, gate: str | None, gate_next: str,
                plan: dict, quality_note: str | None) -> dict:
    """クラスの段。門（選別の前・定量の鎖が使えない入口）を通らなければ `omitted`。"""
    if gate is not None:
        return {"omitted": gate, "facts": [], "next": gate_next}
    comp = composition_result or {}
    absolute = comp.get("absolute")
    if absolute is None:
        return {"omitted": "クラスの付いた特徴が無い", "facts": [], "next": gate_next}
    f = _Facts("C")
    label = "脂質クラス" if comp.get("class_label_kind") == "lipid_class" else "MS-DIAL のオントロジー"
    classes = comp.get("classes") or {}
    clr = comp.get("clr")
    q_thr = absolute["verdict_rule"]["q_threshold"]
    if plan["mode"] == "two":
        cid = absolute["rows"][0]["contrast_id"] if absolute["rows"] else None
        tot = {r["feature_id"]: r for r in absolute["rows"] if r["contrast_id"] == cid}
        clr_rows = {r["feature_id"]: r for r in (clr or {}).get("rows", []) if r["contrast_id"] == cid}
        c = _verdict_counts(list(tot.values()))
        f.add(f"{label}の合計（{len(tot)} クラス、{cid}・仮の向き）: changed {c['changed']}"
              f"（上昇 {c['up']}・下降 {c['down']}）/ unchanged {c['unchanged']} / inconclusive {c['inconclusive']}"
              + (f"。{quality_note}" if quality_note else ""), contrast_id=cid, **c)
        if clr is not None:
            cc = _verdict_counts(list(clr_rows.values()))
            f.add(f"組成（CLR）: changed {cc['changed']} / unchanged {cc['unchanged']} / inconclusive {cc['inconclusive']}"
                  "（総量だけが動く変化は合計にだけ、取り分の変化は CLR に出る）", **cc)
        changed = sorted((r for r in tot.values() if r["verdict"] == "changed"),
                         key=lambda r: (r["q"] if r["q"] is not None else 1.0, r["feature_id"]))
        same = sorted((r for r in tot.values() if r["verdict"] == "unchanged"),
                      key=lambda r: (abs(r["log2fc"] or 0.0), r["feature_id"]))
        rows = []
        for r in [*changed[:_CLASS_ROWS], *same[:_CLASS_ROWS]]:
            cr = clr_rows.get(r["feature_id"]) or {}
            rows.append([r["feature_id"], len(classes.get(r["feature_id"], [])), _g(r["log2fc"]),
                         _g(r["ci_low"]), _g(r["ci_high"]), r["verdict"], _g(cr.get("log2fc")),
                         cr.get("verdict", "")])
        table = _tsv(["class", "n_species", "log2fc_total", "ci_low", "ci_high", "verdict_total",
                      "clr_log2fc", "verdict_clr"], rows)
        items = [r["feature_id"] for r in changed[:3]] or [r["feature_id"] for r in same[:3]]
    else:
        n_tested, n_sig = _omnibus_counts({"engine": absolute})
        f.add(f"{label}の合計: moderated F で q≤{_g(q_thr)} のクラス {n_sig} / 検定 {n_tested}"
              + (f"。{quality_note}" if quality_note else ""), n_tested=n_tested, n_q_below=n_sig)
        om = sorted((o for o in absolute.get("omnibus") or [] if o.get("f_q") is not None),
                    key=lambda o: (o["f_q"], o["feature_id"]))
        rows = [[o["feature_id"], len(classes.get(o["feature_id"], [])), _g(o["f_q"])] for o in om[:2 * _CLASS_ROWS]]
        table = _tsv(["class", "n_species", "f_q_total"], rows)
        items = [o["feature_id"] for o in om[:3]]
    for w in comp.get("warnings") or []:
        if w.get("code") == "class_total_includes_missing_members":
            f.add("欠測のメンバーを 0 として合計したクラスがある（欠測の多い群で低めに出る）",
                  classes=w.get("classes", [])[:20])
    nxt = f'dataset_plot(kind="group_intensity", items={_j(items)})' if items else gate_next
    return {"facts": f.items, "table": table, "next": nxt}


# --- 分子種 --------------------------------------------------------------------------------

def species_facts(ds, *, record: dict | None, plan: dict, top_n: int, pairs_omitted: bool,
                  omitted: str | None, quality_note: str | None) -> dict:
    """分子種の段: 判定ごとの件数（上昇・下降）、出方だけが違う特徴、変わった上位。"""
    if record is None:
        return {"omitted": omitted, "facts": [], "next": "dataset_overview() を条件を直してから呼び直す"}
    f = _Facts("S")
    engine = record["engine"]
    rid = record["result_id"]
    rows = engine["rows"]
    q_thr = engine["verdict_rule"]["q_threshold"]
    cond = (f"{engine['test']}、q≤{_g(q_thr)}、同等 ±{_g(engine['verdict_rule']['equivalence_fc'])} 倍、"
            f"BH={engine['bh_population']['name']}")
    det = bool(engine["presence"]["detection_available"])
    if plan["mode"] == "two":
        cid = rows[0]["contrast_id"] if rows else None
        c = _verdict_counts(rows)
        f.add(f"{cid}（仮の向き）: changed {c['changed']}（上昇 {c['up']}・下降 {c['down']}）/ unchanged {c['unchanged']}"
              f" / inconclusive {c['inconclusive']}（{cond}）" + (f"。{quality_note}" if quality_note else ""),
              contrast_id=cid, condition=cond, **c)
        if det:
            pres = [r["feature_id"] for r in rows if r.get("presence_q") is not None
                    and r["presence_q"] <= q_thr and r["verdict"] == "inconclusive"]
            f.add(f"出方（検出の有無）だけが違う特徴 {len(pres)}（量は判定できない）", n=len(pres),
                  feature_ids=pres[:50])
        else:
            f.add("出方（検出の有無）の違いは評価できない（検出の別が無い）", n=None)
        top = sorted((r for r in rows if r["verdict"] == "changed" and r["q"] is not None),
                     key=lambda r: (r["q"], r["feature_id"]))[:top_n]
        table = _tsv(["contrast_id", "feature_id", "name", "id_level", "log2fc", "ci_low", "ci_high", "q"],
                     [[r["contrast_id"], r["feature_id"], _name(ds, r["feature_id"]), id_level(ds, r["feature_id"]),
                       _g(r["log2fc"]), _g(r["ci_low"]), _g(r["ci_high"]), _g(r["q"])] for r in top])
        ids = [r["feature_id"] for r in top]
        nxt = (f'dataset_plot(kind="volcano", result_id="{rid}", contrast_id="{cid}")'
               + (f' / dataset_feature_search(result_ids=["{rid}#{cid}"], feature_ids={_j(ids[:5])})' if ids else ""))
        return {"facts": f.items, "top": table, "next": nxt}

    n_tested, n_sig = _omnibus_counts(record)
    sig = {o["feature_id"]: o for o in engine.get("omnibus") or []
           if o.get("f_q") is not None and o["f_q"] <= q_thr}
    f.add(f"どれかの群で差（moderated F）: {n_sig} / 検定 {n_tested}（{cond}）"
          + (f"。{quality_note}" if quality_note else ""), n_tested=n_tested, n_q_below=n_sig, condition=cond)
    pair_rows, per = [], {}
    for r in rows:
        if r["feature_id"] in sig:
            per.setdefault(r["contrast_id"], []).append(r)
    for comp_rows in per.values():
        c = _verdict_counts(comp_rows)
        r0 = comp_rows[0]
        pair_rows.append([f"{r0['test']} vs {r0['reference']}", c["changed"], c["up"], c["down"]])
    pairs = _tsv(["pair", "n_changed", "n_up", "n_down"], pair_rows)
    f.add("F で差のある特徴について、対ごとの changed の件数は pairs の表"
          + ("（群が 7 つ以上なので基準との対比だけ。全対は dataset_compare で）" if pairs_omitted else ""),
          pairs={row[0]: row[1:] for row in pair_rows}, pairs_omitted=pairs_omitted)
    top = sorted(sig.values(), key=lambda o: (o["f_q"], o["feature_id"]))[:top_n]
    table = _tsv(["feature_id", "name", "id_level", "f_q"],
                 [[o["feature_id"], _name(ds, o["feature_id"]), id_level(ds, o["feature_id"]), _g(o["f_q"])]
                  for o in top])
    ids = [o["feature_id"] for o in top]
    nxt = (f'dataset_feature_search(result_ids=["{rid}"], feature_ids={_j(ids[:5])})' if ids
           else f"dataset_compare(groups={_j(plan['levels'])})")
    return {"facts": f.items, "top": table, "pairs": pairs, "next": nxt}


# --- 登録 ----------------------------------------------------------------------------------

def overview_id(*, dataset_hash, metadata_revision, matrix_id, compare_id, pca_id, extra=None) -> str:
    return "ovw_" + canonical_hash({"dataset_hash": dataset_hash, "metadata_revision": metadata_revision,
                                    "matrix_id": matrix_id, "compare_id": compare_id, "pca_id": pca_id,
                                    "extra": extra})[:12]


def register_overview(ds, *, result_id: str, matrix_id: str, sections: dict, parents: list[str],
                      settings: dict) -> dict:
    """概観を `ds.results[result_id]`（`kind:"overview"`）に置く。事実は段ごとの `data`（構造）で持つ。"""
    record = {
        "kind": "overview", "schema": OVERVIEW_SCHEMA, "result_id": result_id, "matrix_id": matrix_id,
        "settings": settings,
        "facts": {name: [{"id": x["id"], "text": x["text"], "data": x["data"]} for x in sec.get("facts", [])]
                  for name, sec in sections.items()},
        "provenance": {"result_id": result_id, "kind": "overview", "dataset_hash": dataset_fingerprint(ds),
                       "metadata_revision": getattr(ds, "metadata_revision", 0),
                       "parent_ids": [p for p in parents if p]},
    }
    ds.results[result_id] = record
    return record
