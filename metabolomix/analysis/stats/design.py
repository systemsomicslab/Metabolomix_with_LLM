"""計画行列と対比（統計エンジン spec §6）。

因子の組み合わせごとの cell means（切片なし）に、共変量の列と個体のブロックの
ダミー列を足した行列を作る。`limma.lm_fit` は階数不足を ValueError で拒むので、
ここで先に同じ基準（`limma._is_full_rank`）で調べて `DomainError("DESIGN_RANK_DEFICIENT")` にする。

数値の共変量は中心化し、母標準偏差（ddof=0）で割って入れる。値の桁が大きい共変量（例: 総タンパク
5e7）が XᵀX の固有値の比を潰し、lm_fit の階数の判定に落ちるのを防ぐため。cell の対比は共変量の
平行移動・尺度の変更で変わらない（共変量の係数は標準化した尺度になるが、報告はしない）。
全個体で同じ値の共変量は cell の列と区別できないので DESIGN_RANK_DEFICIENT。

単位（units）は属性 `group` `batch` `subject_id` `attributes`（dict[str, str]）を
持つ任意のオブジェクト。因子・共変量の名前が "group" / "batch" ならその属性を、
それ以外は `attributes[name]` を読む。
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

from metabolomix.analysis.stats.limma import _is_full_rank
from metabolomix.core.atomic_io import DomainError


@dataclass
class Design:
    matrix: np.ndarray            # n×p
    columns: list[str]
    cells: list[str]              # 先頭 len(cells) 列が cell means
    contrasts: dict[str, np.ndarray]   # id -> (p,)
    contrast_info: dict[str, dict]
    omnibus: np.ndarray           # p×m（各 cell − 最初の cell）
    dropped_cells: list[str] = field(default_factory=list)  # 単位が無かった cell（複数因子では DESIGN_INVALID にするので通常は空）


def _invalid(message: str, **details) -> DomainError:
    return DomainError("DESIGN_INVALID", message, details)


def _value(unit, name: str, *, block: bool = False) -> str:
    """単位から因子・共変量・ブロックの値を読む。無ければ DESIGN_INVALID。"""
    if block and name == "subject_id":
        raw = getattr(unit, "subject_id", None)
    elif name in ("group", "batch"):
        raw = getattr(unit, name, None)
    else:
        raw = (getattr(unit, "attributes", None) or {}).get(name)
    if raw is None or str(raw) == "":
        label = (getattr(unit, "unit_id", None) or getattr(unit, "name", None)
                 or getattr(unit, "subject_id", None) or "?")
        raise _invalid(f"単位 {label} に {name} の値が無い", unit=str(label), factor=name)
    return str(raw)


def _cell_name(levels) -> str:
    return ":".join(levels)


def build_design(units, *, factors: list[str], comparisons: list[dict],
                 covariates=(), block: str | None = None) -> Design:
    units = list(units)
    if not units:
        raise _invalid("単位が 1 つも無い")
    if not factors:
        raise _invalid("因子が 1 つも無い")

    values = {f: [_value(u, f) for u in units] for f in factors}
    levels = {f: sorted(set(values[f])) for f in factors}
    unit_cell = [_cell_name(values[f][i] for f in factors) for i in range(len(units))]
    present = set(unit_cell)
    all_cells = [_cell_name(c) for c in itertools.product(*(levels[f] for f in factors))]
    cells = [c for c in all_cells if c in present]
    dropped = [c for c in all_cells if c not in present]
    if dropped and len(factors) >= 2:
        raise _invalid(f"単位が無い cell がある: {dropped}（cell means の交互作用・omnibus が解釈できない）",
                       empty_cells=dropped)
    columns = list(cells)
    cols = [np.array([1.0 if c == cell else 0.0 for c in unit_cell]) for cell in cells]

    for cov in covariates:
        name, kind = cov["name"], cov.get("type")
        raw = [_value(u, name) for u in units]
        if kind == "numeric":
            nums = []
            for text in raw:
                try:
                    x = float(text)
                except ValueError:
                    raise _invalid(f"共変量 {name} に数値でない値 {text!r}", factor=name) from None
                if not math.isfinite(x):
                    raise _invalid(f"共変量 {name} に有限でない値 {text!r}", factor=name)
                nums.append(x)
            arr = np.array(nums)
            sd = float(arr.std())
            if not (sd > 0 and math.isfinite(sd)):
                raise DomainError("DESIGN_RANK_DEFICIENT",
                                  f"数値の共変量 {name} が全個体で同じ値（cell の列と区別できない）",
                                  {"columns": columns + [name], "covariate": name})
            columns.append(name)
            cols.append((arr - arr.mean()) / sd)
        elif kind == "categorical":
            lv = sorted(set(raw))
            ref = cov.get("reference", lv[0])
            if ref not in lv:
                raise _invalid(f"共変量 {name} の基準 {ref!r} が水準に無い", factor=name, levels=lv)
            for level in lv:
                if level != ref:
                    columns.append(f"{name}:{level}")
                    cols.append(np.array([1.0 if v == level else 0.0 for v in raw]))
        else:
            raise _invalid(f"共変量 {name} の type は numeric / categorical", factor=name)

    if block is not None:
        raw = [_value(u, block, block=True) for u in units]
        subjects = sorted(set(raw))
        for s in subjects[1:]:
            columns.append(f"{block}:{s}")
            cols.append(np.array([1.0 if v == s else 0.0 for v in raw]))

    matrix = np.column_stack(cols)
    p = matrix.shape[1]
    if not _is_full_rank(matrix) or np.linalg.matrix_rank(matrix) < p:
        raise DomainError("DESIGN_RANK_DEFICIENT",
                          "計画行列が階数不足（個体が群に入れ子など）", {"columns": columns})

    index = {c: i for i, c in enumerate(columns)}

    def cell_of(assign: dict, cid: str) -> str:
        for f in factors:
            if assign.get(f) not in levels[f]:
                raise _invalid(f"比較 {cid}: 因子 {f} の水準 {assign.get(f)!r} が存在しない",
                               id=cid, factor=f, levels=levels[f])
        name = _cell_name(assign[f] for f in factors)
        if name not in cells:
            raise _invalid(f"比較 {cid}: cell {name} に単位が無い", id=cid, cell=name)
        return name

    contrasts: dict[str, np.ndarray] = {}
    info: dict[str, dict] = {}
    for comp in comparisons:
        cid = comp.get("id")
        if not cid:
            raise _invalid("比較に id が無い")
        if cid in contrasts:
            raise _invalid(f"比較 id が重複: {cid}", id=cid)
        vec = np.zeros(p)
        if comp.get("kind") == "interaction":
            fa = comp.get("factors")
            if not (isinstance(fa, (list, tuple)) and len(fa) == 2 and len(factors) == 2
                    and set(fa) == set(factors)):
                raise _invalid(f"比較 {cid}: 交互作用は 2 因子のモデルで、その 2 因子を factors に",
                               id=cid)
            lv = comp.get("levels") or {}
            a, b = fa
            try:
                (a0, a1), (b0, b1) = lv[a], lv[b]
            except (KeyError, TypeError, ValueError):
                raise _invalid(f"比較 {cid}: levels は各因子の [基準, 比較]", id=cid) from None
            if a0 == a1 or b0 == b1:
                raise _invalid(f"比較 {cid}: 交互作用の基準と比較が同じ水準", id=cid)
            ref_cells, test_cells = [], []
            for la, lb, sign in ((a1, b1, 1.0), (a0, b1, -1.0), (a1, b0, -1.0), (a0, b0, 1.0)):
                name = cell_of({a: la, b: lb}, cid)
                vec[index[name]] += sign
                (test_cells if sign > 0 else ref_cells).append(name)
            info[cid] = {"kind": "interaction", "factors": [a, b],
                         "levels": {a: [a0, a1], b: [b0, b1]},
                         "reference_cells": ref_cells, "test_cells": test_cells}
        else:
            f = comp.get("factor")
            if f not in factors:
                raise _invalid(f"比較 {cid}: 因子 {f!r} がモデルに無い", id=cid)
            at = comp.get("at")
            others = [o for o in factors if o != f]
            if others and (not isinstance(at, dict) or set(at) != set(others)):
                raise _invalid(f"比較 {cid}: 複数因子では at に他の因子の水準が必須", id=cid)
            if comp.get("reference") == comp.get("test"):
                raise _invalid(f"比較 {cid}: reference と test が同じ水準", id=cid)
            ref_name = cell_of({**(at or {}), f: comp.get("reference")}, cid)
            test_name = cell_of({**(at or {}), f: comp.get("test")}, cid)
            vec[index[test_name]] += 1.0
            vec[index[ref_name]] -= 1.0
            info[cid] = {"kind": "simple", "factor": f, "reference": comp["reference"],
                         "test": comp["test"], "reference_cells": [ref_name],
                         "test_cells": [test_name]}
            if at:
                info[cid]["at"] = dict(at)
        contrasts[cid] = vec

    k = len(cells)
    omnibus = np.zeros((p, max(k - 1, 0)))
    for j in range(1, k):
        omnibus[j, j - 1] = 1.0
        omnibus[0, j - 1] = -1.0
    return Design(matrix, columns, cells, contrasts, info, omnibus, [])
