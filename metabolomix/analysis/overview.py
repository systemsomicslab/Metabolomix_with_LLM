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
from metabolomix.analysis.stats.data import SAMPLE_TABLE_POLARITY
from metabolomix.analysis.stats.preprocess import REMERGE_LOSS, runs_remedy
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

def _polarity_of(row: dict) -> str:
    """統合したデータセットのドリフトの兆候の行の極性（`polarity` があれば「<極性> の」）。"""
    return f"{row['polarity']} の" if row.get("polarity") else ""


def _unsure_polarities(drift: list[dict]) -> list[str]:
    """ρ に使った QC の注入順がファイルの順かもしれない行の極性（統合したデータセットの行だけ。出た順）。"""
    return list(dict.fromkeys(str(d["polarity"]) for d in drift
                              if d.get("injection_order_may_be_file_order") and d.get("polarity")))


def _drift_order_caveat(drift: list[dict]) -> str:
    """ドリフトの兆候の事実に添える但し書き（ρ に使った QC の注入順がファイルの順かもしれない行が無ければ空）。

    統合したデータセットでは、そうである極性を名指しし、その極性の測定の入れ直し方（`preprocess.runs_remedy`。
    陽イオンを最後に。`dataset_preprocess` の注意と同じ文と順）を添える: 統合したハンドルへのシートは陽イオンの
    測定だけを決めるので、陰イオンの注入順は陰イオンのハンドルにシートを当てて統合し直すまで変わらない。統合で
    なければ試料の表の injection_order に入れる。"""
    unsure = [d for d in drift if d.get("injection_order_may_be_file_order")]
    if not unsure:
        return ""
    pols = _unsure_polarities(drift)
    where = ", ".join(pols)
    remedy = (runs_remedy(pols).rstrip("。") if pols
              else "実際の測定順は試料の表の injection_order に入れる")
    return (f"（ただし {where + ' の ' if where else ''}QC の注入順は MS-DIAL の既定〔ファイルの順〕かシートの宣言で"
            "ファイル名の順かもしれず、ρ はファイル名の並びとの相関で、実際の測定順でのドリフトを見ていないかもしれない。"
            + remedy + "）")


def _order_next(drift: list[dict], *, drift_unsure: bool, sheet_order_biased: bool) -> list[str]:
    """実際の測定順の入れ方（品質の段の `next`）。

    統合したハンドルへのシートは陽イオン（試料の表の極性）の測定だけを決める。陽イオン以外の極性の注入順が
    確かでなければ、その極性のハンドルにシートを当てて統合し直す。統合し直すと統合したハンドルに当てたシートと
    除外（`dataset_exclude`）は引き継がれない（統合は陽イオンのハンドルのシートと極性ごとのハンドルの除外を引き継ぐ。
    `REMERGE_LOSS`）ので、陽イオンの測定も入れるなら陽イオンのハンドルに当ててから一緒に統合し直し、入れないなら、
    統合したハンドルで入れていれば陽イオンのハンドルにも当てる。統合したハンドルに当てていた除外は、統合し直して
    できたハンドルに当て直す。`sheet_order_biased` は交絡の注入順の偏り（統合したデータセットの試料の表＝陽イオンの
    測定）。"""
    pols = _unsure_polarities(drift)
    others = [p for p in pols if p != SAMPLE_TABLE_POLARITY]
    sheet = sheet_order_biased or (drift_unsure and (not pols or SAMPLE_TABLE_POLARITY in pols))
    if others:
        handles = ([SAMPLE_TABLE_POLARITY] if sheet else []) + others
        keep = ("" if sheet else f"統合したハンドルに {SAMPLE_TABLE_POLARITY} のシートを当てていれば "
                f"{SAMPLE_TABLE_POLARITY} のハンドルにも当てておき、")
        return [f"{'・'.join(handles)} のハンドルに dataset_set_sample_metadata（その極性の raw と実際の測定順の"
                f"シート）→ dataset_merge_polarities（{REMERGE_LOSS}。{keep}統合したハンドルに当てていた"
                "除外は統合し直してできたハンドルに当て直す）"]
    return ["dataset_set_sample_metadata（injection_order に実際の測定順を入れる）"] if sheet else []


_BLOCK_GROUPS = 10                 # 技術反復らしいまとまりの事実の文に載せる群の上限（全部は data に）
_BLOCK_LAYOUT = {"consecutive": "で連続", "interleaved": "で巡回（1, 2, 3, 1, 2, 3, …）"}
_BLOCK_ORDER = {"injection_order": "注入順", "table": "試料の表の順"}


def _range(values: list[float], digits: int = 2) -> str:
    lo, hi = min(values), max(values)
    return _g(lo, digits) if _g(lo, digits) == _g(hi, digits) else f"{_g(lo, digits)}〜{_g(hi, digits)}"


def _replicate_block_fact(found: dict) -> str:
    """データから見つけた技術反復らしいまとまり（`quality.replicate_blocks` の `found`）の事実の文。"""
    groups = found["groups"]
    k = found["k"]
    shown = []
    for g in groups[:_BLOCK_GROUPS]:
        name = "(群なし)" if g["group"] is None else g["group"]
        layout = ("相関から束ねた" if g["layout"] == "data"
                  else _BLOCK_ORDER.get(g["order"], g["order"]) + _BLOCK_LAYOUT[g["layout"]])
        shown.append(f"{name}（{layout}）" + " ".join("[" + ", ".join(b) + "]" for b in g["blocks"]))
    more = f"; ほか {len(groups) - _BLOCK_GROUPS} 群は記録の data" if len(groups) > _BLOCK_GROUPS else ""
    skipped = found.get("skipped") or []
    other = ("。見ていない群（注入が少ない・割り切れない・実検出の特徴が足りない）: "
             + ", ".join("(群なし)" if s["group"] is None else str(s["group"]) for s in skipped) if skipped else "")
    return (f"subject_id が無いのに、群の中の注入が {k} 本ずつの技術反復らしいまとまりに分かれる（{len(groups)} 群。"
            "群のすべての注入で実検出の特徴の log2 の相関で、まとまりの中と間の差 "
            + _range([g["margin"] for g in groups])
            + f"・近い {k - 1} 本が仲間の割合 " + _range([100 * g["purity"] for g in groups], 3) + "%"
            + "・特徴ごとの ICC の中央値 " + _range([g["median_icc"] for g in groups]) + "）: "
            + "; ".join(shown) + more + other
            + "。いまは注入を独立な試料として数えている。データから見た候補なので、同じ個体の反復かを実験の記述で"
              "確かめてから subject_id を入れる（自動では入れていない）")


def quality_facts(*, rsd: dict, failed: dict, detection: dict, drift: list[dict], confounding: dict,
                  replicate_warnings: list, verdict: str, reasons: list[str], n_qc: int, n_blank: int,
                  drift_note: str | None, pca_id: str | None, drift_rerunnable: bool = True,
                  replicate_blocks: dict | None = None) -> dict:
    """品質の段。判定（ok / caution / poor / not_evaluable）と理由のコードを最初の事実に置く。

    `drift_rerunnable` は drift_note の補正を dataset_preprocess(drift_mode="on") で当てられるか（統合した
    データセットで、ほかの極性が補正の条件を満たさなければ偽。`next` に drift_mode="on" を勧めない）。
    交絡の注入順が「ファイルの順かもしれない」（`injection_order.may_be_file_order`）なら、偏りの事実にその但し書きを
    添え、事実のデータに `injection_order_may_be_file_order=True` を足す（偽なら足さない）。ドリフトの兆候も同じ:
    ρ に使った QC の注入順がファイルの順かもしれない行（`drift` の行の `injection_order_may_be_file_order`）があれば、
    兆候の有無の事実に但し書きとデータの印を足す。どちらかで但し書きを付けたら `next` に実際の測定順の入れ方を足す
    （統合したデータセットで陽イオン以外の極性なら、その極性のハンドルへのシートと統合し直し。`_order_next`）。
    統合し直す手は `next` の中で除外（`dataset_exclude`）と `dataset_preprocess(drift_mode="on")` より前に置く: 統合し直すと
    統合したハンドルに当てた除外は引き継がれず、前処理の行列も前のハンドルに残るので、書いた順に行えば両方が
    統合し直してできたハンドルに当たる。

    `replicate_blocks`（`quality.replicate_blocks` の戻り値）が `found` なら、群ごとのまとまりの候補を挙げた事実を足し
    （データに戻り値そのもの）、`next` に subject_id の入れ方を足す（名前からの反復の注意と同じ 1 つ。subject_id は入れない）。"""
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
        by = {p: ids for p, ids in (failed.get("by_polarity") or {}).items() if ids}
        where = ("（極性ごとの総量で: " + "; ".join(f"{p} {', '.join(ids[:5])}" for p, ids in by.items()) + "）"
                 if by else "")
        f.add(f"失敗した QC 注入 {len(failed['failed'])} 本（総量が QC の中央値の 0.2 倍未満）: "
              + ", ".join(failed["failed"][:5]) + where, **failed)
    if detection.get("status") == "ok":
        f.add(f"生体試料のセル: 実検出 {_pct(detection['detected_frac'])}・gap-fill {_pct(detection['gap_filled_frac'])}"
              f"・欠測 {_pct(detection['missing_frac'])}", **detection)
    else:
        f.add(f"実検出と gap-fill の別は評価できない（{detection.get('reason')}）", **detection)
    flagged = [d for d in drift if d.get("flag")]
    drift_unsure = _drift_order_caveat(drift)
    drift_mark = {"injection_order_may_be_file_order": True} if drift_unsure else {}
    if flagged:
        f.add("ドリフトの兆候: " + "; ".join(f"{_polarity_of(d)}バッチ {d['batch']} で QC の総量と注入順の"
                                           f" ρ={_g(d['rho'], 2)}（QC {d['n_qc']} 本）"
                                           for d in flagged[:3]) + drift_unsure, batches=drift, **drift_mark)
    elif drift:
        unit = "極性・バッチ" if any(d.get("polarity") for d in drift) else "バッチ"
        f.add(f"ドリフトの兆候なし（QC が 5 本以上の {len(drift)} {unit}で |ρ|<0.6 か p≥0.05）" + drift_unsure,
              batches=drift, **drift_mark)
    status = confounding.get("batch_status")
    inj = (confounding.get("injection_order") or {}).get("status")
    may_file = bool((confounding.get("injection_order") or {}).get("may_be_file_order"))
    conf_text = {"complete": "群とバッチが完全交絡（群差とバッチ差を区別できない。比較は回していない）",
                 "partial": f"群とバッチが部分交絡（Cramér の V={_g(confounding.get('cramers_v'), 2)}）",
                 "none": "群とバッチの交絡なし", "single_batch": "バッチは 1 つ",
                 "not_evaluable": "群とバッチの交絡は評価できない（バッチが不明な試料がある）",
                 "not_computed": "群とバッチの交絡は評価していない（前処理の行列が作れない）"}.get(status, str(status))
    inj_text = ""
    if inj == "biased":
        inj_text = "・注入順が群で偏っている" + (
            "（ただし注入順は MS-DIAL の既定〔ファイルの順〕かシートの宣言でファイル名の順かもしれず、偏りは"
            "ファイル名の並びを見ているだけかもしれない。実際の測定順は試料の表の injection_order に入れる）"
            if may_file else "")
    f.add(conf_text + inj_text,
          batch_status=status, cramers_v=confounding.get("cramers_v"), injection_order=inj,
          batch_table=confounding.get("batch_table"),
          **({"injection_order_may_be_file_order": True} if may_file else {}))
    if replicate_warnings:
        f.add(f"反復注入らしい試料名が {len(replicate_warnings)} 組あるのに subject_id が無い（独立な試料として数えている）",
              stems=[w.get("stem") for w in replicate_warnings])
    blocks_found = (replicate_blocks or {}).get("status") == "found"
    if blocks_found:
        f.add(_replicate_block_fact(replicate_blocks), replicate_blocks=replicate_blocks)
    if drift_note:
        f.add(drift_note, drift="off_undecided")
    nxt = [f'dataset_plot(kind="pca", result_id="{pca_id}", color_by="batch")'] if pca_id else []
    order_next = _order_next(drift, drift_unsure=bool(drift_unsure), sheet_order_biased=may_file and inj == "biased")
    remerge = [x for x in order_next if "dataset_merge_polarities" in x]
    nxt += remerge                      # 統合し直す手は除外・drift_mode="on" より前（統合し直すと引き継がれない）
    if failed.get("failed"):
        nxt.append(f"dataset_exclude(samples={_j(failed['failed'][:5])})")
    if drift_note and drift_rerunnable:
        nxt.append('dataset_preprocess(drift_mode="on")')
    nxt += [x for x in order_next if x not in remerge]
    if replicate_warnings or blocks_found:
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


def _exploratory_fact(f: _Facts, exploratory: str | None) -> None:
    """探索専用のデータセットなら、その段の先頭に印の事実を足す（`exploratory` は理由。無ければ何もしない）。

    比較は断らずに回す（dataset_compare・書き出しと同じ）が、未完了の実行では欠けた検体が「その群には無い」
    ように見えうるので、比較の数字を使う段（全体・クラス・分子種）ごとに同じ印を付ける。"""
    if exploratory is None:
        return
    f.add(f"探索専用のデータセット（{exploratory}）。比較は回したが、未完了の実行で欠けた検体が「その群には無い」"
          "ように見えうるので、この段の差は参考（書き出し・図にも探索専用の印が載る）",
          exploratory=True, reason=exploratory)


def global_facts(*, pca_result: dict | None, labels: dict, record: dict | None, plan: dict,
                 omitted: str | None, quality_note: str | None, exploratory: str | None = None) -> dict:
    """全体の段: PCA の寄与率・群の重心・群で説明される分散（η²）と、比較の全体の件数。

    `exploratory`（探索専用のデータセットの理由）があれば、比較を回した段の先頭に印の事実を置く。"""
    f = _Facts("G")
    if record is not None:
        _exploratory_fact(f, exploratory)
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

def _non_lipid_count(comp: dict) -> int:
    """組成が脂質クラスから外した特徴（脂質と読めない注釈。`non_lipid`）の数。外していなければ 0。"""
    return int((comp.get("non_lipid") or {}).get("n_features") or 0)


def _non_lipid_fact(f: _Facts, comp: dict) -> None:
    """脂質アッセイで、脂質と読めない注釈の特徴を非脂質として外した事実（無ければ足さない）。"""
    block = comp.get("non_lipid") or {}
    n = _non_lipid_count(comp)
    if not n:
        return
    top = list((block.get("by_ontology") or {}).items())[:5]
    f.add(f"脂質と読めない注釈の特徴 {n} 個は非脂質（{block.get('label')}）としてクラスの集計"
          "（合計・CLR・分子種の数）から外した。元のクラスの上位: " + "、".join(f"{c} {k}" for c, k in top),
          n_features=n, by_ontology=dict(top),
          examples=[e.get("name") for e in (block.get("examples") or [])[:5]])


def _relabel_fact(f: _Facts, comp: dict) -> None:
    """脂質アッセイで、名前か ClassyFire のクラスと組成式から MS-DIAL の脂質クラスへ移して数えた事実（無ければ足さない）。"""
    block = comp.get("lipid_class_relabel") or {}
    n = int(block.get("n_features") or 0)
    if not n:
        return
    top = (block.get("moves") or [])[:5]
    f.add(f"MS-DIAL の脂質クラスでないクラス（研究室 MSP の ClassyFire など）の注釈の特徴 {n} 個は、名前かクラスと組成式から"
          "MS-DIAL の脂質クラスへ移して数えた"
          f"（うち名前で読めず組成式で救ったもの {int(block.get('n_rescued') or 0)} 個。注釈そのものは変えない）。"
          "上位: " + "、".join(f"{m['from']}→{m['to']} {m['n_features']}" for m in top),
          n_features=n, n_rescued=int(block.get("n_rescued") or 0),
          moves=[{k: m.get(k) for k in ("from", "to", "route", "n_features")} for m in top])


def class_facts(*, composition_result: dict | None, gate: str | None, gate_next: str,
                plan: dict, quality_note: str | None, exploratory: str | None = None) -> dict:
    """クラスの段。門（選別の前・定量の鎖が使えない入口）を通らなければ `omitted`。

    `exploratory`（探索専用のデータセットの理由）があれば、先頭に印の事実を置く。"""
    if gate is not None:
        return {"omitted": gate, "facts": [], "next": gate_next}
    comp = composition_result or {}
    absolute = comp.get("absolute")
    if absolute is None:
        n_non_lipid = _non_lipid_count(comp)
        omitted = ("クラスの付いた特徴が無い" if not n_non_lipid else
                   f"脂質クラスの付いた特徴が無い（脂質と読めない注釈の特徴 {n_non_lipid} 個は非脂質として外した）")
        return {"omitted": omitted, "facts": [], "next": gate_next}
    f = _Facts("C")
    _exploratory_fact(f, exploratory)
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
    _relabel_fact(f, comp)
    _non_lipid_fact(f, comp)
    nxt = f'dataset_plot(kind="group_intensity", items={_j(items)})' if items else gate_next
    return {"facts": f.items, "table": table, "next": nxt}


# --- 分子種 --------------------------------------------------------------------------------

def species_facts(ds, *, record: dict | None, plan: dict, top_n: int, pairs_omitted: bool,
                  omitted: str | None, quality_note: str | None, exploratory: str | None = None,
                  report_cut: bool = False) -> dict:
    """分子種の段: 判定ごとの件数（上昇・下降）、出方だけが違う特徴、変わった上位。

    `exploratory`（探索専用のデータセットの理由）があれば、先頭に印の事実を置く。`report_cut`（求めた件数を
    上限に下げたとき）なら、上位の表を `top_n` で切ったときに `top_truncated: True` と `top_n_more`（載せなかった
    行数）を足す。"""
    if record is None:
        return {"omitted": omitted, "facts": [], "next": "dataset_overview() を条件を直してから呼び直す"}
    f = _Facts("S")
    _exploratory_fact(f, exploratory)
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
        ranked = sorted((r for r in rows if r["verdict"] == "changed" and r["q"] is not None),
                        key=lambda r: (r["q"], r["feature_id"]))
        top = ranked[:top_n]
        table = _tsv(["contrast_id", "feature_id", "name", "id_level", "log2fc", "ci_low", "ci_high", "q"],
                     [[r["contrast_id"], r["feature_id"], _name(ds, r["feature_id"]), id_level(ds, r["feature_id"]),
                       _g(r["log2fc"]), _g(r["ci_low"]), _g(r["ci_high"]), _g(r["q"])] for r in top])
        ids = [r["feature_id"] for r in top]
        nxt = (f'dataset_plot(kind="volcano", result_id="{rid}", contrast_id="{cid}")'
               + (f' / dataset_feature_search(result_ids=["{rid}#{cid}"], feature_ids={_j(ids[:5])})' if ids else ""))
        return {"facts": f.items, "top": table, "next": nxt, **_cut(len(ranked), top_n, report_cut)}

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
    ranked = sorted(sig.values(), key=lambda o: (o["f_q"], o["feature_id"]))
    top = ranked[:top_n]
    table = _tsv(["feature_id", "name", "id_level", "f_q"],
                 [[o["feature_id"], _name(ds, o["feature_id"]), id_level(ds, o["feature_id"]), _g(o["f_q"])]
                  for o in top])
    ids = [o["feature_id"] for o in top]
    nxt = (f'dataset_feature_search(result_ids=["{rid}"], feature_ids={_j(ids[:5])})' if ids
           else f"dataset_compare(groups={_j(plan['levels'])})")
    return {"facts": f.items, "top": table, "pairs": pairs, "next": nxt, **_cut(len(ranked), top_n, report_cut)}


def _cut(n_ranked: int, top_n: int, report: bool) -> dict:
    """上位の表を上限で切った印（`report` のときだけ。切っていなければ足さない）。"""
    more = n_ranked - top_n
    return {"top_truncated": True, "top_n_more": more} if report and more > 0 else {}


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
