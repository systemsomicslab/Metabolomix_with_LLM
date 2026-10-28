"""割れたアライメントスポットの検出と、事前に決めた規則でのまとめ（spec 2026-10-28 §4.3）。

MS-DIAL のアライメントは、1 つの物質を 2 つ以上のスポットに割ることがある——RT がずれた
注入（例: 最初の数注入）だけが別のスポットに行く、ピークの形の揺れで試料ごとに行き先が
入れ替わる、など。割れたまま検定すると、片方のスポットの「無い」試料に gap-fill の背景が
入り、群の差でも何でもない差が出る。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

- 候補: 同じ付加体、|Δm/z| ≤ `ppm`、`min_rt_gap` < |ΔRT| ≤ `max_rt_gap` の 2 スポット
  （m/z の昇順に並べて掃くので、全組を比べない）。
- 割れの判定: 試料ごとの受け持ち（検出状態があれば検出、無ければ 2 スポットの値の大きい方）が
  補い合うこと。検出状態があれば「両方で検出された試料」が `max_both_detected` 以下で、
  どちらでも検出されていない試料は値で補い（`_owners`）、受け持ちの決まる試料が
  `MIN_COMPLEMENTARY` 以上。無ければ「片方の値が他方の `MIN_DOMINANCE` 倍以上の試料」が
  `MIN_COMPLEMENTARY` 以上。どちらの場合も、2 スポットがそれぞれ 1 試料以上を受け持つ。
- 型: 片方の受け持ちが注入順で連続した塊 → `injection_block`（群と重なっても優先）、
  そうでなく受け持ちが群と一致 → `group_aligned`、ほか → `interleaved`。
- まとめてよい（`mergeable`）のは注釈が両立する組のうち、`injection_block` と、ΔRT が
  `max_interleaved_rt_gap`（既定 0.15 分 = MS-DIAL のアライメントの RT 許容 0.1 分の 1.5 倍。
  アライメントの窓で割れたもので、別の異性体ではない。MTBLS6965 の陰イオンでは FA 20:3 の
  異性体 401/402 が 0.20 分離れており、0.2 では別の異性体をまとめてしまう）以下の
  `interleaved` だけ。
- 実際にまとめる（`merged`）のは、mergeable の組のうち**互いに最も近い相手**どうしだけ
  （`_mutual_nearest`: A の mergeable な相手の中で |ΔRT| が最小なのが B、かつ B のそれが A。
  最小が同点なら相手なし）。まとめは常に 2 スポットで、1 つのスポットは高々 1 組にしか入らない。
  ほかの mergeable な組は `merge_reason: not_mutual_nearest` で載せるだけ。MTBLS6965 の陰イオンの
  FA 20:3 は 2 つの異性体が同じ後半の注入で一緒にずれ、401（後半の n-9）が 403（前半の n-9、
  ΔRT 0.088）とも 402（前半の n-3/n-6、ΔRT 0.199）とも組になる——互いに最も近い組だけを
  残すと 401↔403 だけがまとまる。
  まとめると判定に使うデータが変わるので、迷う組はまとめない側に倒す（コントローラの決定 2026-10-29）。
  `group_aligned` は群ごとに別の異性体が出ている可能性があるので、印を付けるだけにする。
- まとめ方は試料ごとに 2 スポットの大きい方の値（gap-fill の値は背景で、本物のピークより
  大きくならない）。`merge_matrix_columns` が唯一の実装で、統計の側（feature_search・
  analysis_spec）はこれを呼ぶ。
"""
from __future__ import annotations

import math
import statistics
from pathlib import Path
from typing import NamedTuple

import numpy as np

from metabolomix.core.atomic_io import DomainError

__all__ = ["SCHEMA", "TYPES", "MergedMatrix", "find_split_pairs", "merge_matrix_columns",
           "names_compatible", "split_partners"]

SCHEMA = "alignment-qc.v1"
TYPES = ("injection_block", "group_aligned", "interleaved")
#: 「補い合い」: 片方の値が他方の MIN_DOMINANCE 倍以上の試料（検出状態があれば、受け持ちの
#: 決まる試料）が、信号のある試料の MIN_COMPLEMENTARY 以上（`_owners`）。
MIN_DOMINANCE = 3.0
MIN_COMPLEMENTARY = 0.90
#: 「受け持ちが群と一致」: どの群でも多い方のスポットの受け持ちが GROUP_PURITY 以上で、
#: 2 スポットのどちらも少なくとも 1 つの群の多い方になる（試料 1 本の揺れで
#: `interleaved` に落としてまとめない——迷う組はまとめない側に倒す）。
GROUP_PURITY = 0.90
#: 群の中で型を 1 つにするときの順（前ほど強い）。群ごとに受け持ちが分かれる組が
#: 1 つでもあれば群全体を `group_aligned` と読む。
_TYPE_PRECEDENCE = ("group_aligned", "interleaved", "injection_block")
#: 注釈が無いとみなす名前。
_UNANNOTATED = frozenset({"", "unknown", "null", "na", "n/a"})
_BAD = "DATASET_BAD_REQUEST"


def _number(value) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    return float(value)


def _check_params(ppm, min_rt_gap, max_rt_gap, max_both_detected,
                  max_interleaved_rt_gap) -> dict:
    values = {"ppm": _number(ppm), "min_rt_gap": _number(min_rt_gap),
              "max_rt_gap": _number(max_rt_gap), "max_both_detected": _number(max_both_detected),
              "max_interleaved_rt_gap": _number(max_interleaved_rt_gap)}
    if values["ppm"] is None or values["ppm"] <= 0:
        raise DomainError(_BAD, "ppm は正の数", {"field": "ppm"})
    if values["min_rt_gap"] is None or values["min_rt_gap"] < 0:
        raise DomainError(_BAD, "min_rt_gap は 0 以上の数", {"field": "min_rt_gap"})
    if values["max_rt_gap"] is None or values["max_rt_gap"] <= values["min_rt_gap"]:
        raise DomainError(_BAD, "max_rt_gap は min_rt_gap より大きい数", {"field": "max_rt_gap"})
    if values["max_both_detected"] is None or not 0 <= values["max_both_detected"] <= 1:
        raise DomainError(_BAD, "max_both_detected は 0〜1", {"field": "max_both_detected"})
    if values["max_interleaved_rt_gap"] is None or values["max_interleaved_rt_gap"] < 0:
        raise DomainError(_BAD, "max_interleaved_rt_gap は 0 以上の数",
                          {"field": "max_interleaved_rt_gap"})
    return values


# ---------- 試料の注入順と群 ----------

def _rows_by_sample(ds, names: list[str]) -> list[dict | None]:
    """実験情報シートの行を ds.sample_names の並びへ当てる（当たらない試料は None）。

    行は `ds.sample_names` と同じ順の契約（`apply_metadata`）だが、pipeline の run は
    run が解決した行をそのまま載せる。生ファイル名の語幹・sample_id の順に照合し
    （末尾の空白は無視）、どれでも全部当たらず件数が同じなら並びで当てる。
    """
    rows = [r for r in (getattr(ds, "sample_metadata_rows", None) or []) if isinstance(r, dict)]
    if not rows:
        return [None] * len(names)
    wanted = [n.strip() for n in names]
    for key in (lambda r: Path(str(r.get("source_file") or "").replace("\\", "/")).stem.strip(),
                lambda r: str(r.get("sample_id") or "").strip()):
        keyed = {}
        for row in rows:
            keyed.setdefault(key(row), row)
        if all(n in keyed for n in wanted):
            return [keyed[n] for n in wanted]
    if len(rows) == len(names):
        return list(rows)
    return [None] * len(names)


def _injection_order(names: list[str], rows: list[dict | None],
                     given: dict | None) -> tuple[list[float], str]:
    """試料ごとの注入順の値と、その出所（manifest / metadata / filename）。

    引数で渡されればそれ（`metadata`）、無ければ実験情報シートの `injection_order`
    （出所が利用者のシートなら `manifest`、mzTab などから推した値なら `metadata`）。
    1 試料でも欠ければ、推した順を信じずファイル名（試料名）の順にする
    （pipeline の参照ファイルの選び方 `choose_alignment_reference` と同じ規約）。
    """
    if given is not None:
        order = [_number(given.get(n, given.get(n.strip()))) for n in names]
        if all(v is not None for v in order):
            return order, "metadata"
    order = [_number((row or {}).get("injection_order")) for row in rows]
    if all(v is not None for v in order) and names:
        sources = {(((row or {}).get("provenance") or {}).get("injection_order") or {})
                   .get("source") for row in rows}
        return order, "manifest" if sources == {"user_manifest"} else "metadata"
    ranked = sorted(range(len(names)), key=lambda i: (names[i], i))
    order = [0.0] * len(names)
    for rank, i in enumerate(ranked, start=1):
        order[i] = float(rank)
    return order, "filename"


def _order_label(value: float):
    return int(value) if float(value).is_integer() else value


# ---------- 特徴の注釈 ----------

def _clean_name(name) -> str | None:
    if not isinstance(name, str) or name.strip().casefold() in _UNANNOTATED:
        return None
    return name.strip()


def _feature_name(ds, fid: str) -> str | None:
    """MS/MS の照合で付いた rank 1 の名前（SME）、無ければ SML の注釈。"""
    for candidate in (getattr(ds, "feature_candidates", None) or {}).get(fid) or []:
        name = _clean_name(candidate.get("chemical_name"))
        if name:
            return name
    name = _clean_name(((getattr(ds, "feature_metadata", None) or {}).get(fid) or {}).get("name"))
    if name:
        return name
    annotation = (getattr(ds, "feature_annotations", None) or {}).get(fid) or {}
    return None if annotation.get("ambiguous") else _clean_name(annotation.get("name"))


def _feature_adducts(ds) -> dict[str, str | None]:
    """特徴ごとのスポットの付加体（SML の `adduct_ions`。注釈の無い特徴にも入っている）。

    SML が無い特徴は SME の rank 1 候補の付加体で補う。
    """
    out: dict[str, str | None] = {}
    for row in getattr(ds, "sml_rows", None) or []:
        adduct = _clean_name(row.get("adduct_ions"))
        for ref in str(row.get("SMF_ID_REFS") or "").split("|"):
            if ref.strip():
                out.setdefault(ref.strip(), adduct)
    candidates = getattr(ds, "feature_candidates", None) or {}
    for fid in ds.feature_ids:
        if out.get(fid) is None:
            first = (candidates.get(fid) or [{}])[0]
            out[fid] = _clean_name(first.get("adduct"))
    return out


def names_compatible(a, b) -> bool:
    """2 スポットの注釈が両立するか（同じ名前・同じ和組成、または片方が注釈なし）。

    和組成は `lipid_identity.lipid_species_of`（GOSLIN の種のレベル。`PC 16:0_18:1` と
    `PC 34:1` は同じ）。脂質名として読めない名前は、大文字小文字を無視した一致だけを見る。
    """
    a, b = _clean_name(a), _clean_name(b)
    if a is None or b is None or a.casefold() == b.casefold():
        return True
    from metabolomix.msdial.lipid_identity import lipid_species_of
    species_a, species_b = lipid_species_of(a), lipid_species_of(b)
    return species_a is not None and species_a == species_b


# ---------- 割れの判定と型 ----------

def _owners(va, vb, da, db, max_both: float):
    """(受け持ち, 両方で検出された割合)。割れでなければ None。

    受け持ちは試料ごとに 0（a）/ 1（b）/ -1（決まらない・どちらにも信号が無い）。

    - 検出状態が無いとき: 信号（値 > 0）のある試料で、片方の値が他方の `MIN_DOMINANCE` 倍以上の
      試料が `MIN_COMPLEMENTARY` 以上なら割れ。受け持ちは値の大きい方。
    - 検出状態があるとき: 両方で検出された試料（分母はどちらかで検出された試料）が
      `max_both` 以下。受け持ちは検出された方（両方なら値の大きい方）。**どちらでも検出されて
      いない試料は値で補う**——片方が他方の `MIN_DOMINANCE` 倍以上ならその方、そうでなければ
      決まらない。そのうえで、信号か検出のある試料のうち受け持ちの決まった試料が
      `MIN_COMPLEMENTARY` 以上なら割れ。
      検出だけで決めないのは、MS-DIAL の検出が疎な run（MTBLS6965 の陽イオンは全セルの 20%、
      特徴の中央値で 42 試料中 3）で、本物のピークの多くが gap-fill の側に入るから。検出だけ
      だと「どちらも数試料でしか検出されない無関係な 2 スポット」が両方検出 0% で割れに見える
      （陽イオンで 1560 組）。値だけだと、ΔRT が小さい割れ（陰イオンの FA 20:3 401/403、
      0.088 分）で gap-fill が隣のピークを拾い 3 倍の差が出ず、見逃す。
    """
    va, vb = np.nan_to_num(va, nan=0.0), np.nan_to_num(vb, nan=0.0)
    owner = np.full(va.shape, -1, dtype=int)
    high, low = np.maximum(va, vb), np.minimum(va, vb)
    signal = high > 0
    dominant = signal & (high >= MIN_DOMINANCE * low)
    if da is None:
        if not signal.any() or dominant.sum() / signal.sum() < MIN_COMPLEMENTARY:
            return None
        fraction = None
        owner[signal] = np.where(va[signal] >= vb[signal], 0, 1)
    else:
        either = da | db
        if not either.any():
            return None
        both = da & db
        fraction = float(both.sum() / either.sum())
        if fraction > max_both:
            return None
        owner[da & ~db] = 0
        owner[db & ~da] = 1
        owner[both] = np.where(va[both] >= vb[both], 0, 1)
        filled = ~either & dominant
        owner[filled] = np.where(va[filled] >= vb[filled], 0, 1)
        considered = signal | either
        if (owner[considered] >= 0).mean() < MIN_COMPLEMENTARY:
            return None
    if not (owner == 0).any() or not (owner == 1).any():
        return None
    return owner, fraction


def _block(owner, ranked: list[int]) -> tuple[int, list[int]] | None:
    """注入順で連続した塊だけを受け持つスポット（0/1）と、その試料の添字。

    信号の無い試料は塊を切らない（飛ばして並べる）。両方が塊なら受け持ちの少ない方。
    """
    sequence = [i for i in ranked if owner[i] >= 0]
    blocks = []
    for spot in (0, 1):
        positions = [k for k, i in enumerate(sequence) if owner[i] == spot]
        if positions[-1] - positions[0] + 1 == len(positions):
            blocks.append((len(positions), spot, [sequence[k] for k in positions]))
    if not blocks:
        return None
    _, spot, members = min(blocks)
    return spot, members


def _group_aligned(owner, groups: list[str | None]) -> bool:
    by_group: dict[str, list[int]] = {}
    for spot, group in zip(owner, groups):
        if group is not None and spot >= 0:
            by_group.setdefault(group, []).append(int(spot))
    if len(by_group) < 2:
        return False
    majority = set()
    for spots in by_group.values():
        n_b = sum(spots)
        top = 1 if n_b > len(spots) - n_b else 0
        if max(n_b, len(spots) - n_b) / len(spots) < GROUP_PURITY:
            return False
        majority.add(top)
    return majority == {0, 1}


# ---------- 公開: 割れの組を探す ----------

def find_split_pairs(ds, *, ppm: float = 15.0, min_rt_gap: float = 0.05,
                     max_rt_gap: float = 1.0, max_both_detected: float = 0.10,
                     max_interleaved_rt_gap: float = 0.15,
                     injection_order: dict[str, int] | None = None) -> dict:
    """割れたスポットの組と、その連結成分（3 スポット以上の連鎖を 1 つにする）。

    戻り値は alignment-qc.v1（`lipidmix://docs/output-format/reanalysis` §9）。
    `pairs[].features` は、`injection_block` なら [塊でない側, 塊の側]、ほかは RT の昇順。
    `delta_rt` は `rt[1] − rt[0]`（`injection_block` では塊の側のずれ。符号つき）。
    `groups[].features` の先頭は代表（受け持ちの最も多いスポット）で、
    `merge_matrix_columns` はこれを代表の ID として残す。
    `pairs[].mergeable` は注釈が両立し、型が `injection_block`、または ΔRT ≤
    `max_interleaved_rt_gap` の `interleaved`。`groups[]` はすべての組の連結成分（報告用）、
    `pairs[].merged` は互いに最も近い mergeable な組（`_mutual_nearest`）で、`merge_groups[]` は
    その組ごとの 2 スポット（まとめる単位）。mergeable だがまとめない組は `merge_reason:
    not_mutual_nearest`。`groups[].mergeable` は成分全体がその 1 つのまとめのときだけ真。
    """
    params = _check_params(ppm, min_rt_gap, max_rt_gap, max_both_detected,
                           max_interleaved_rt_gap)
    if getattr(ds, "feature_matrix", None) is None:
        raise DomainError("DATASET_REQUIRED", "定量行列が無い", {})
    matrix = np.asarray(ds.feature_matrix, dtype=float)
    feature_ids = [str(f) for f in ds.feature_ids]
    names = [str(n) for n in ds.sample_names]
    mask = getattr(ds, "detected_mask", None)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != matrix.shape:
            mask = None
    rows = _rows_by_sample(ds, names)
    order, order_source = _injection_order(names, rows, injection_order)
    ranked = sorted(range(len(names)), key=lambda i: (order[i], i))
    groups = [(row or {}).get("group") or None for row in rows]
    adducts = _feature_adducts(ds)
    metadata = getattr(ds, "feature_metadata", None) or {}

    features = []
    for index, fid in enumerate(feature_ids):
        meta = metadata.get(fid) or {}
        mz, rt = _number(meta.get("mz")), _number(meta.get("rt"))
        if mz is not None and rt is not None and mz > 0:
            features.append((mz, index, fid, rt))
    features.sort()

    pairs = []
    for k, (mz_a, index_a, fid_a, rt_a) in enumerate(features):
        limit = mz_a * (1 + params["ppm"] * 1e-6)
        for j in range(k + 1, len(features)):
            mz_b, index_b, fid_b, rt_b = features[j]
            if mz_b > limit:
                break
            gap = abs(rt_b - rt_a)
            if not params["min_rt_gap"] < gap <= params["max_rt_gap"]:
                continue
            if adducts.get(fid_a) != adducts.get(fid_b):
                continue
            assessed = _owners(matrix[index_a], matrix[index_b],
                               None if mask is None else mask[index_a],
                               None if mask is None else mask[index_b],
                               params["max_both_detected"])
            if assessed is None:
                continue
            pairs.append(_pair(ds, (index_a, fid_a, mz_a, rt_a), (index_b, fid_b, mz_b, rt_b),
                               assessed, ranked, order, groups, names, adducts,
                               params["max_interleaved_rt_gap"]))
    pairs.sort(key=lambda p: (p["mz"][0], p["features"]))
    for n, pair in enumerate(pairs, start=1):
        pair["pair_id"] = f"P{n}"
    matched = _mutual_nearest(pairs)
    for pair in pairs:
        pair["merged"] = pair["pair_id"] in matched
        pair["merge_reason"] = ("not_mutual_nearest" if pair["mergeable"] and not pair["merged"]
                                else None)
    components = _components(pairs, feature_ids, matrix, mask)
    for n, component in enumerate(components, start=1):
        component["group_id"] = f"G{n}"
        # 成分全体が 1 つのまとめのときだけ真（まとめの判断は merge_groups）。
        component["mergeable"] = (len(component["features"]) == 2
                                  and all(p["merged"] for p in pairs
                                          if p["pair_id"] in component["pair_ids"]))
    group_of = {f: c["group_id"] for c in components for f in c["features"]}
    # 互いに最も近い組は特徴を共有しないので、成分はどれも 2 スポット（= その組）。
    merges = _components([p for p in pairs if p["merged"]], feature_ids, matrix, mask)
    for n, merge in enumerate(merges, start=1):
        merge["merge_group_id"] = f"M{n}"
        merge["group_id"] = group_of[merge["features"][0]]
        merge["mergeable"] = True

    by_type = {t: sum(p["type"] == t for p in pairs) for t in TYPES}
    block_shifts = [p["delta_rt"] for p in pairs if p["type"] == "injection_block"]
    return {
        "schema": SCHEMA,
        "params": {**params, "min_dominance": MIN_DOMINANCE,
                   "min_complementary": MIN_COMPLEMENTARY, "group_purity": GROUP_PURITY,
                   "signal": "detected" if mask is not None else "values"},
        "injection_order_source": order_source,
        "pairs": pairs,
        "groups": components,
        "merge_groups": merges,
        "summary": {
            "n_pairs": len(pairs), "by_type": by_type,
            "n_features_affected": len({f for p in pairs for f in p["features"]}),
            "injection_block_median_delta_rt": (statistics.median(block_shifts)
                                                if block_shifts else None),
            "n_mergeable_pairs": sum(p["mergeable"] for p in pairs),
            "n_groups": len(components),
            "n_merged_pairs": len(merges),
        },
    }


def _pair(ds, a, b, assessed, ranked, order, groups, names, adducts,
          max_interleaved_rt_gap: float) -> dict:
    owner, fraction = assessed
    block = _block(owner, ranked)
    if block is not None:
        kind = "injection_block"
        # 塊の側を後ろに置き、delta_rt を「塊のずれ」にする。
        if block[0] == 0:
            a, b = b, a
            owner = np.where(owner >= 0, 1 - owner, owner)
    else:
        kind = "group_aligned" if _group_aligned(owner, groups) else "interleaved"
        if (b[3], b[1]) < (a[3], a[1]):
            a, b = b, a
            owner = np.where(owner >= 0, 1 - owner, owner)
    fids = [a[1], b[1]]
    pair_names = [_feature_name(ds, a[1]), _feature_name(ds, b[1])]
    compatible = names_compatible(*pair_names)
    block_out = None
    if block is not None:
        members = block[1]
        block_out = {"spot": b[1],
                     "injections": [_order_label(order[members[0]]),
                                    _order_label(order[members[-1]])]}
    return {
        "pair_id": None,
        "features": fids,
        "mz": [a[2], b[2]],
        "rt": [a[3], b[3]],
        "delta_rt": b[3] - a[3],
        "delta_ppm": (b[2] - a[2]) / a[2] * 1e6,
        "adduct": adducts.get(a[1]),
        "names": pair_names,
        "names_compatible": compatible,
        "type": kind,
        "owner": {name: (fids[int(s)] if s >= 0 else None) for name, s in zip(names, owner)},
        "n_owned": [int((owner == 0).sum()), int((owner == 1).sum())],
        "both_detected_fraction": fraction,
        "block": block_out,
        "mergeable": compatible and (
            kind == "injection_block"
            or (kind == "interleaved" and abs(b[3] - a[3]) <= max_interleaved_rt_gap + 1e-9)),
    }


#: |ΔRT| を同点とみなす幅（分。浮動小数の差だけで最も近い相手が決まらないように）。
_RT_TIE = 1e-6


def _mutual_nearest(pairs: list[dict]) -> set[str]:
    """mergeable の組のうち、互いに最も近い（|ΔRT| が最小の）相手どうしの組の pair_id。

    特徴ごとに mergeable な相手の |ΔRT| の最小を取り、最小が同点なら相手なしとする
    （どちらとまとめるか決められない組はまとめない）。
    """
    best: dict[str, tuple[float, str | None]] = {}
    for pair in pairs:
        if not pair["mergeable"]:
            continue
        gap = abs(pair["delta_rt"])
        a, b = pair["features"]
        for fid, other in ((a, b), (b, a)):
            current = best.get(fid)
            if current is None or gap < current[0] - _RT_TIE:
                best[fid] = (gap, other)
            elif abs(gap - current[0]) <= _RT_TIE:
                best[fid] = (current[0], None)
    nearest = {fid: other for fid, (_, other) in best.items() if other is not None}
    return {p["pair_id"] for p in pairs if p["mergeable"]
            and nearest.get(p["features"][0]) == p["features"][1]
            and nearest.get(p["features"][1]) == p["features"][0]}


def _components(pairs: list[dict], feature_ids: list[str], matrix, mask) -> list[dict]:
    """`pairs` を辺とする連結成分（ID と mergeable は呼ぶ側が付ける）。

    代表は検出された試料（無ければ値が最大の試料）の最も多いスポット。
    """
    parent: dict[str, str] = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for pair in pairs:
        a, b = (find(f) for f in pair["features"])
        if a != b:
            parent[max(a, b)] = min(a, b)
    position = {fid: i for i, fid in enumerate(feature_ids)}
    members: dict[str, list[str]] = {}
    for fid in sorted(parent, key=position.__getitem__):
        members.setdefault(find(fid), []).append(fid)

    components = []
    for features in sorted(members.values(), key=lambda m: position[m[0]]):
        rows = [position[f] for f in features]
        if mask is not None:
            counts = mask[rows].sum(axis=1)
        else:
            values = np.nan_to_num(matrix[rows], nan=-np.inf)
            winners = values.argmax(axis=0)
            counts = np.bincount(winners, minlength=len(rows))
        representative = features[int(np.argmax(counts))]
        ordered = [representative] + [f for f in features if f != representative]
        chosen = [p for p in pairs if p["features"][0] in features]
        types = sorted({p["type"] for p in chosen}, key=_TYPE_PRECEDENCE.index)
        components.append({"features": ordered, "pair_ids": [p["pair_id"] for p in chosen],
                           "type": types[0], "types": types})
    return components


def split_partners(qc: dict) -> dict[str, dict]:
    """特徴 ID → `{"group_id", "merge_group_id", "partners", "component", "type", "mergeable"}`。

    割れの組を 1 つでも持つ特徴だけが載る。まとめの判断は `merge_groups`（互いに最も近い
    mergeable な組）から取る: まとめる特徴は `partners` がその相手 1 つで `mergeable` が真、`type` は
    その組の型。まとめない特徴は `partners` が空・`mergeable` が偽・`type` は全体の成分の型。
    `component` は全体の成分（すべての組）のほかの特徴で、報告用。
    型は混ざれば `group_aligned` > `interleaved` > `injection_block` の順で強い方。
    """
    out = {}
    for group in qc.get("groups") or []:
        for fid in group["features"]:
            out[fid] = {"group_id": group["group_id"], "merge_group_id": None, "partners": [],
                        "component": [f for f in group["features"] if f != fid],
                        "type": group["type"], "mergeable": False}
    for merge in qc.get("merge_groups") or []:
        for fid in merge["features"]:
            out[fid].update({"merge_group_id": merge["merge_group_id"],
                             "partners": [f for f in merge["features"] if f != fid],
                             "type": merge["type"], "mergeable": merge["mergeable"]})
    return out


# ---------- 公開: まとめる ----------

class MergedMatrix(NamedTuple):
    matrix: np.ndarray                 # (試料, 特徴)。まとめた列は代表の位置に 1 列
    feature_ids: list[str]
    merged_from: dict[str, list[str]]  # 代表 → まとめた特徴（代表が先頭）
    detected: np.ndarray | None        # 渡されたときだけ。まとめた列は「どれかで検出」


def merge_matrix_columns(matrix, feature_ids, groups, *, detected=None) -> MergedMatrix:
    """割れの群ごとに列を 1 本へまとめる（試料ごとに大きい方の値。NaN は無視）。

    matrix は (試料, 特徴)（`analysis_spec` / `recipe_matrix` の並び）。`groups` は特徴 ID の
    並びの並びで、先頭を代表とする（`find_split_pairs` の `merge_groups[].features` をそのまま
    渡せる。まとめてよいかの判断は `find_split_pairs` の側で済んでいる）。行列に無い特徴は
    飛ばし、残りが 1 つ以下ならまとめない。群どうしが特徴を共有していれば ValueError。
    """
    values = np.asarray(matrix, dtype=float)
    ids = [str(f) for f in feature_ids]
    position = {fid: i for i, fid in enumerate(ids)}
    mask = None if detected is None else np.asarray(detected, dtype=bool)
    seen: set[str] = set()
    plans = []
    for group in groups:
        present = [str(f) for f in group if str(f) in position]
        if seen.intersection(present):
            raise ValueError(f"割れの群が特徴を共有している: {sorted(seen.intersection(present))}")
        seen.update(present)
        if len(present) >= 2:
            plans.append(present)
    drop: set[int] = set()
    merged_from: dict[str, list[str]] = {}
    values = values.copy()
    mask = None if mask is None else mask.copy()
    for present in plans:
        columns = [position[f] for f in present]
        target = columns[0]
        stacked = values[:, columns]
        with np.errstate(all="ignore"):
            merged = np.fmax.reduce(stacked, axis=1)
        values[:, target] = merged
        if mask is not None:
            mask[:, target] = mask[:, columns].any(axis=1)
        drop.update(columns[1:])
        merged_from[present[0]] = present
    keep = [i for i in range(len(ids)) if i not in drop]
    return MergedMatrix(matrix=values[:, keep], feature_ids=[ids[i] for i in keep],
                        merged_from=merged_from,
                        detected=None if mask is None else mask[:, keep])
