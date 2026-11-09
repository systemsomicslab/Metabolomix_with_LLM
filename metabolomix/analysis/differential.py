"""差次的結果の volcano 点列（純ロジック層、MCP 非依存）。

統計エンジンの 2 群の結果（`stats_results.as_two_group`）の行から、volcano の点（log2fc・-log10 p・
up/down/ns）を作る。検定そのものは統計エンジン（`analysis/stats/`）が持つ。
"""

from __future__ import annotations

import math
import sys


def volcano_data(results, q_thr=0.05, log2fc_thr=1.0):
    """volcano 用の点列（log2fc, -log10 p, 有意フラグ up/down/ns）を返す。

    行が `significant` を持てば（統計エンジンの行。`stats_results.as_two_group` が判定込みで付けた印。R16）、
    色はその印と log2fc の符号だけで決める——判定できない行（`inconclusive`）は q が小さくても `ns`。
    印を持たない行は q と log2fc の閾値で決める。
    p が 0.0（アンダーフロー。有限の検定済みの値）は最小の正の float で -log10 を取る（R17。約 307.65。
    NaN＝検定不能にしない）。行の p 自体は変えない。
    """
    points = []
    for r in results:
        p = r.get("p")
        q = r.get("q", p)
        fc = r.get("log2fc")
        neg_log10_p = (-math.log10(max(p, sys.float_info.min))
                       if (p is not None and math.isfinite(p) and p >= 0) else math.nan)
        sig = "ns"
        flag = r.get("significant")
        if flag is not None:
            if flag and fc is not None and math.isfinite(fc):
                sig = "up" if fc >= 0 else "down"
        elif q is not None and math.isfinite(q) and q <= q_thr and fc is not None and math.isfinite(fc):
            if fc >= log2fc_thr:
                sig = "up"
            elif fc <= -log2fc_thr:
                sig = "down"
        points.append({"feature": r["feature"], "log2fc": fc,
                       "neg_log10_p": neg_log10_p, "sig": sig})
    return points
