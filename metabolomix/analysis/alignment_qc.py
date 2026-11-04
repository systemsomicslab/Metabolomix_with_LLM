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
- 残すのは受け持ち（`n_owned`）の多い側、同数なら小さい ID（`_keep_first`）。割れの群の並びは
  **残す ID が先頭**（`merge_groups[].features`・`usable_merge_groups` の戻り・`merged_from`）。
"""
from __future__ import annotations

import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import NamedTuple

import numpy as np

from metabolomix.core.atomic_io import DomainError

__all__ = ["SCHEMA", "TYPES", "MergedMatrix", "find_split_pairs", "merge_matrix_columns",
           "id_order", "names_compatible", "order_label", "read_qc_file",
           "sample_injection_order", "split_partners", "usable_merge_groups"]

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
                  max_interleaved_rt_gap, min_block_size) -> dict:
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
    if isinstance(min_block_size, bool) or not isinstance(min_block_size, int) or min_block_size < 1:
        raise DomainError(_BAD, "min_block_size は 1 以上の整数", {"field": "min_block_size"})
    values["min_block_size"] = min_block_size
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

    引数で渡されればそれ（`metadata`。1 試料でも欠ける・数でなければ DomainError——
    呼び出し側が明示した順を黙って捨てない）、無ければ実験情報シートの `injection_order`
    （出所が利用者のシートなら `manifest`、mzTab などから推した値なら `metadata`）。
    シートの値が 1 試料でも欠ければ、推した順を信じずファイル名（試料名）の順にする
    （pipeline の参照ファイルの選び方 `choose_alignment_reference` と同じ規約）。
    """
    if given is not None:
        if not isinstance(given, dict):
            raise DomainError(_BAD, "injection_order は {試料名: 注入順}",
                              {"field": "injection_order"})
        order = [_number(given.get(n, given.get(n.strip()))) for n in names]
        missing = [n for n, v in zip(names, order) if v is None]
        if missing:
            raise DomainError(_BAD, f"injection_order に無い・数でない試料がある（{len(missing)} 件）",
                              {"field": "injection_order", "missing": missing})
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


def sample_injection_order(ds, given: dict | None = None) -> tuple[list[float], str]:
    """`ds.sample_names` の並びの注入順の値と、その出所（manifest / metadata / filename）。

    割れの組の判定（`find_split_pairs`）と同じ規則で、質量精度の注入順の傾向
    （`mass_accuracy.injection_trend`）が共有する。
    """
    names = [str(n) for n in ds.sample_names]
    return _injection_order(names, _rows_by_sample(ds, names), given)


def order_label(value: float):
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


def _block(owner, ranked: list[int], min_block_size: int) -> tuple[int, list[int]] | None:
    """注入順で連続した塊だけを受け持つスポット（0/1）と、その試料の添字。

    信号の無い試料は塊を切らない（飛ばして並べる）。両方が塊なら受け持ちの少ない方。
    どちらかの受け持ちが `min_block_size` 未満なら塊にしない（1 試料だけを受け持つ
    スポットは、どこにあっても「連続」に見えるため。コントローラの決定）。
    """
    if min((owner == 0).sum(), (owner == 1).sum()) < min_block_size:
        return None
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
                     max_interleaved_rt_gap: float = 0.15, min_block_size: int = 2,
                     injection_order: dict[str, int] | None = None,
                     dataset_hash: str | None = None) -> dict:
    """割れたスポットの組と、その連結成分（3 スポット以上の連鎖を 1 つにする）。

    戻り値は alignment-qc.v1（`lipidmix://docs/output-format/reanalysis` §9）。
    `pairs[].features` は、`injection_block` なら [塊でない側, 塊の側]、ほかは RT の昇順。
    `delta_rt` は `rt[1] − rt[0]`（`injection_block` では塊の側のずれ。符号つき）。
    `groups[].features` は特徴 ID の小さい順（すべて整数なら数の順）、`merge_groups[].features` は
    残す ID が先頭（`keep`。受け持ちの多い側、同数なら小さい ID。`keep_basis` が `n_owned` /
    `smaller_id_tie`、`n_owned` が特徴 ID ごとの受け持ちの数）で、`merge_matrix_columns` はその先頭を
    残す。`injection_block` には両スポットが
    `min_block_size` 試料以上を受け持つことを求める。付加体が分からないスポットは組にしない。
    `pairs[].mergeable` は注釈が両立し、型が `injection_block`、または ΔRT ≤
    `max_interleaved_rt_gap` の `interleaved`。`groups[]` はすべての組の連結成分（報告用）、
    `pairs[].merged` は互いに最も近い mergeable な組（`_mutual_nearest`）で、`merge_groups[]` は
    その組ごとの 2 スポット（まとめる単位）。mergeable だがまとめない組は `merge_reason:
    not_mutual_nearest`。`groups[].mergeable` は成分全体がその 1 つのまとめのときだけ真。
    `dataset_hash` は計算済みの `dataset_fingerprint(ds)`（任意。行列全体の hash を計算し直さない）。
    """
    params = _check_params(ppm, min_rt_gap, max_rt_gap, max_both_detected,
                           max_interleaved_rt_gap, min_block_size)
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
            # 付加体が分からないスポットは組にしない（None どうしを「同じ」と読まない）。
            if adducts.get(fid_a) is None or adducts.get(fid_a) != adducts.get(fid_b):
                continue
            assessed = _owners(matrix[index_a], matrix[index_b],
                               None if mask is None else mask[index_a],
                               None if mask is None else mask[index_b],
                               params["max_both_detected"])
            if assessed is None:
                continue
            pairs.append(_pair(ds, (index_a, fid_a, mz_a, rt_a), (index_b, fid_b, mz_b, rt_b),
                               assessed, ranked, order, groups, names, adducts,
                               params["max_interleaved_rt_gap"], params["min_block_size"]))
    pairs.sort(key=lambda p: (p["mz"][0], p["features"]))
    for n, pair in enumerate(pairs, start=1):
        pair["pair_id"] = f"P{n}"
    matched = _mutual_nearest(pairs)
    for pair in pairs:
        pair["merged"] = pair["pair_id"] in matched
        pair["merge_reason"] = ("not_mutual_nearest" if pair["mergeable"] and not pair["merged"]
                                else None)
    components = _components(pairs, feature_ids)
    for n, component in enumerate(components, start=1):
        component["group_id"] = f"G{n}"
        # 成分全体が 1 つのまとめのときだけ真（まとめの判断は merge_groups）。
        component["mergeable"] = (len(component["features"]) == 2
                                  and all(p["merged"] for p in pairs
                                          if p["pair_id"] in component["pair_ids"]))
    group_of = {f: c["group_id"] for c in components for f in c["features"]}
    # 互いに最も近い組は特徴を共有しないので、成分はどれも 2 スポット（= その組）。
    merges = _components([p for p in pairs if p["merged"]], feature_ids)
    by_pair_id = {p["pair_id"]: p for p in pairs}
    for n, merge in enumerate(merges, start=1):
        merge["merge_group_id"] = f"M{n}"
        merge["group_id"] = group_of[merge["features"][0]]
        merge["mergeable"] = True
        pair = by_pair_id[merge["pair_ids"][0]]
        merge["features"], merge["keep_basis"] = _keep_first(pair)
        merge["keep"] = merge["features"][0]
        merge["n_owned"] = dict(zip(pair["features"], pair["n_owned"]))

    by_type = {t: sum(p["type"] == t for p in pairs) for t in TYPES}
    block_shifts = [p["delta_rt"] for p in pairs if p["type"] == "injection_block"]
    from metabolomix.analysis.result_state import dataset_fingerprint

    return {
        "schema": SCHEMA,
        # どのデータセット（極性・run）で作ったか。統計の側が別の run の JSON を取り込まない
        # ための照合に使う（`analysis_spec` の `features.alignment_qc`）。
        "dataset_fingerprint": dataset_hash or dataset_fingerprint(ds),
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
          max_interleaved_rt_gap: float, min_block_size: int) -> dict:
    owner, fraction = assessed
    block = _block(owner, ranked, min_block_size)
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
                     "injections": [order_label(order[members[0]]),
                                    order_label(order[members[-1]])]}
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


def _keep_first(pair: dict) -> tuple[list[str], str]:
    """組の `features` と `n_owned` から（残す ID が先頭の並び, keep_basis）。

    残すのは受け持ちの多い側（`n_owned`）、同数なら小さい ID（`smaller_id_tie`）。
    """
    owned = dict(zip(pair["features"], pair["n_owned"]))
    smaller, larger = id_order(pair["features"])
    if owned[smaller] == owned[larger]:
        return [smaller, larger], "smaller_id_tie"
    ordered = [smaller, larger] if owned[smaller] > owned[larger] else [larger, smaller]
    return ordered, "n_owned"


def id_order(ids) -> list[str]:
    """特徴 ID の小さい順（すべて整数として読めれば数の順、そうでなければ文字列の順）。"""
    ids = [str(f) for f in ids]
    if all(f.strip().isdigit() for f in ids):
        return sorted(ids, key=lambda f: (int(f), f))
    return sorted(ids)


def _components(pairs: list[dict], feature_ids: list[str]) -> list[dict]:
    """`pairs` を辺とする連結成分（ID と mergeable は呼ぶ側が付ける）。

    `features` は特徴 ID の小さい順（`id_order`）。まとめの組の残す ID は呼ぶ側が決める（`_keep_first`）。
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
        ordered = id_order(features)
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


def read_qc_file(path_text: str, ds, *, field: str,
                 dataset_hash: str | None = None) -> tuple[dict, str]:
    """`dataset_alignment_qc` が書いた JSON を読み、`(中身, sha256)` を返す。

    alignment-qc.v1 で `merge_groups` を持ち、`dataset_fingerprint` が今のデータセットと
    一致することを求める。特徴 ID は小さい整数で、別の極性・別の run にも同じ ID がある。
    別のデータセットで作った JSON を取り込むと無関係な特徴を黙ってまとめるので、指紋が
    合わない（または無い）JSON は拒む。`field` はエラーに書く欄の名前（呼ぶ側の引数・spec の欄）。
    `dataset_hash` は計算済みの `dataset_fingerprint(ds)`（任意）。
    """
    path = Path(path_text)
    if not path.is_file():
        raise DomainError("ANALYSIS_INPUT_NOT_FOUND",
                          f"{field} の JSON が見つからない: {path}",
                          {"field": field, "path": str(path)})
    try:
        raw = path.read_bytes()
        qc = json.loads(raw.decode("utf-8-sig"))
    except (OSError, ValueError) as exc:
        raise DomainError("ANALYSIS_INPUT_INVALID",
                          f"{field} を JSON として読めない: {exc}", {"field": field}) from exc
    groups = qc.get("merge_groups") if isinstance(qc, dict) else None
    if not isinstance(qc, dict) or qc.get("schema") != SCHEMA or not isinstance(groups, list):
        raise DomainError("ANALYSIS_INPUT_INVALID",
                          f"{field} は {SCHEMA}（dataset_alignment_qc の output_path の JSON）で、"
                          "merge_groups を持つこと", {"field": field})
    if dataset_hash is None:
        from metabolomix.analysis.result_state import dataset_fingerprint

        dataset_hash = dataset_fingerprint(ds)
    if qc.get("dataset_fingerprint") != dataset_hash:
        raise DomainError("ANALYSIS_INPUT_INVALID",
                          f"{field} は今のデータセットで作った JSON ではない"
                          "（dataset_fingerprint が合わない、または無い）。このデータセットで "
                          "dataset_alignment_qc を実行し直す", {"field": field})
    return qc, hashlib.sha256(raw).hexdigest()


def usable_merge_groups(qc: dict, feature_ids) -> tuple[list[list[str]], list[list[str]]]:
    """`(まとめる群, まとめない群)`。まとめるのは `merge_groups` の mergeable な群のうち、
    `feature_ids` に 2 つ以上そろうもの（まとめの判断は `find_split_pairs` が持ち、ここでは
    足さない・広げない）。そろわない群（別の run の JSON など）は 2 つ目に返す。

    群の並びは残す ID（`keep`）が先頭。`keep` の無い古い JSON は小さい ID を残す（`id_order`）。"""
    given = [_keep_leading(group) for group in qc.get("merge_groups") or []
             if group.get("mergeable")]
    present = {str(f) for f in feature_ids}
    usable = [g for g in given if sum(f in present for f in g) >= 2]
    return usable, [g for g in given if g not in usable]


def _keep_leading(group: dict) -> list[str]:
    """`merge_groups` の 1 群の特徴 ID（`keep` が先頭。無ければ `id_order`）。"""
    ordered = id_order(group["features"])
    keep = group.get("keep")
    if keep is None:
        return ordered
    return [str(keep), *(f for f in ordered if f != str(keep))]


# ---------- 公開: まとめる ----------

class MergedMatrix(NamedTuple):
    matrix: np.ndarray                 # (試料, 特徴)。まとめた列は残す ID の位置に 1 列
    feature_ids: list[str]
    merged_from: dict[str, list[str]]  # 残す ID → まとめた特徴（残す ID が先頭）
    detected: np.ndarray | None        # 渡されたときだけ。まとめた列は「どれかで検出」


def merge_matrix_columns(matrix, feature_ids, groups, *, detected=None) -> MergedMatrix:
    """割れの群ごとに列を 1 本へまとめる（試料ごとに大きい方の値。NaN は無視）。

    matrix は (試料, 特徴)（`analysis_spec` / `recipe_matrix` の並び）。`groups` は特徴 ID の
    並びの並び（`find_split_pairs` の `merge_groups[].features` をそのまま渡せる。まとめてよいかの
    判断は `find_split_pairs` の側で済んでいる）。まとめた列の ID は群の中で**行列にある最初の ID**
    （渡した並びに従う。並べ直さない）で、その列の位置に置く。行列に無い特徴は
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
