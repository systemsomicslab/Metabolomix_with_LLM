"""割れたアライメントスポットの検出と、事前に決めた規則でのまとめ（spec 2026-10-28 §4.3）。

合成の DatasetState を直に組む（mzTab を通さない）。値は「受け持つ試料では 1000 前後、
受け持たない試料では gap-fill の 20 前後」、検出状態を渡すときは受け持つ試料だけ True。

閾値の既定（両方で検出された試料 10% 以下・補い合い 3 倍が 90% 以上）は実データ
（MTBLS6965 の陽イオン・陰イオン）で変えていない。検出状態があるときは、どちらでも検出されて
いない試料を値（3 倍）で補い、受け持ちの決まる試料が 90% 以上を求める
（`test_undetected_samples_are_owned_by_value` / `test_sparse_detection_alone_is_not_a_split`）。
検出だけだと陽イオン（検出は全セルの 20%）で無関係な疎なスポットの組が 1560 組混ざり、値だけだと
陰イオンの FA 20:3 401/403（ΔRT 0.088 分）を見逃した。

まとめの規則（コントローラの決定 2026-10-29）: mergeable は注釈が両立する `injection_block`
（両スポットが `min_block_size`=2 試料以上を受け持つ）と ΔRT ≤ 0.15 分の `interleaved`。実際に
まとめるのは互いに |ΔRT| が最も近い mergeable な相手どうしの 2 スポットだけで、残す ID は小さい方。
MTBLS6965 では CE 1265/1266・1185/1186・1257/1258 と陰イオンの 401/403・441/442 をまとめ、
402 と CE 1178/1179（interleaved ΔRT 0.612）はまとめない（task-5 の報告）。群の差は見ていない。
"""
from __future__ import annotations

import numpy as np
import pytest

from metabolomix.analysis import alignment_qc as Q
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState

N = 42
NAMES = [f"S{i:02d}" for i in range(1, N + 1)]
ADDUCT = "[M+H]1+"


def _rows(groups, orders=None, source="user_manifest"):
    orders = orders or list(range(1, N + 1))
    return [{"sample_id": name, "source_file": f"C:/raw/{name}.wiff", "role": "sample",
             "group": group, "injection_order": order, "include": True,
             "provenance": {"injection_order": {"value": order, "source": source,
                                                "confidence": "confirmed"}}}
            for name, group, order in zip(NAMES, groups, orders)]


GROUPS_HALVES = ["A"] * 21 + ["B"] * 21          # 注入順で前半が A、後半が B
GROUPS_ALTERNATE = ["A", "B"] * 21               # 注入順で A と B が交互


def _make_ds(features, *, groups=GROUPS_HALVES, detected=False, rows=True, background=20.0):
    """features: (fid, mz, rt, owned: set[int] (0 始まりの試料添字), name, adduct)。

    受け持たない試料の値は gap-fill の背景 `background + (i % 3)`（0 なら 0 そのもの）。
    """
    ds = DatasetState()
    ds.sample_names = list(NAMES)
    ds.feature_ids = [f[0] for f in features]
    matrix, mask = [], []
    for k, (fid, mz, rt, owned, name, adduct) in enumerate(features):
        values = [1000.0 + 10 * i + k if i in owned else (background + (i % 3) if background else 0.0)
                  for i in range(N)]
        matrix.append(values)
        mask.append([i in owned for i in range(N)])
        ds.feature_metadata[fid] = {"name": None, "mz": mz, "rt": rt, "inchikey": None}
        ds.feature_candidates[fid] = ([{"rank": 1, "chemical_name": name, "adduct": adduct}]
                                      if name else [])
        if name:
            ds.feature_metadata[fid]["name"] = name
    ds.feature_matrix = np.asarray(matrix, dtype=float)
    ds.detected_mask = np.asarray(mask, dtype=bool) if detected else None
    ds.sml_rows = [{"SML_ID": f[0], "SMF_ID_REFS": f[0], "adduct_ions": f[5]} for f in features]
    if rows:
        ds.sample_metadata_rows = _rows(groups)
    return ds


ALL = set(range(N))
FIRST_TWO = {0, 1}


def test_injection_block_first_two_injections():
    # 最初の 2 注入だけが RT +0.6 の別スポット（MTBLS6965 の陽イオンの CE と同じ形）
    ds = _make_ds([("1265", 650.6000, 20.00, ALL - FIRST_TWO, "CE 18:2", ADDUCT),
                   ("1266", 650.6050, 20.61, FIRST_TWO, None, ADDUCT)])
    for detected in (False, True):
        ds.detected_mask = (np.asarray([[i not in FIRST_TWO for i in range(N)],
                                        [i in FIRST_TWO for i in range(N)]])
                            if detected else None)
        qc = Q.find_split_pairs(ds)
        assert qc["schema"] == "alignment-qc.v1"
        assert qc["injection_order_source"] == "manifest"
        (pair,) = qc["pairs"]
        assert pair["type"] == "injection_block"
        assert pair["features"] == ["1265", "1266"]          # 塊の側が後ろ
        assert pair["block"] == {"spot": "1266", "injections": [1, 2]}
        assert pair["delta_rt"] == pytest.approx(0.61)
        assert pair["adduct"] == ADDUCT
        assert pair["names"] == ["CE 18:2", None] and pair["names_compatible"]
        assert pair["mergeable"] is True
        assert pair["owner"]["S01"] == "1266" and pair["owner"]["S42"] == "1265"
        assert (pair["both_detected_fraction"] == 0.0) if detected else \
            (pair["both_detected_fraction"] is None)
        assert qc["summary"]["by_type"] == {"injection_block": 1, "group_aligned": 0,
                                            "interleaved": 0}
        assert qc["summary"]["n_features_affected"] == 2
        assert qc["summary"]["injection_block_median_delta_rt"] == pytest.approx(0.61)
        (group,) = qc["groups"]
        assert group["features"] == ["1265", "1266"] and group["mergeable"] is True
        (merge,) = qc["merge_groups"]
        assert merge["features"] == ["1265", "1266"] and merge["mergeable"] is True
        assert merge["group_id"] == group["group_id"]


def test_interleaved_pair():
    # 受け持ちが注入順でも群でも入れ替わる（陰イオンの FA 20:3 401/403 の形）
    owned = {i for i in range(N) if i % 3 == 0}
    ds = _make_ds([("401", 305.2464, 3.000, ALL - owned, "FA 20:3", "[M-H]1-"),
                   ("403", 305.2495, 3.088, owned, "FA 20:3", "[M-H]1-")], detected=True)
    (pair,) = Q.find_split_pairs(ds)["pairs"]
    assert pair["type"] == "interleaved" and pair["block"] is None
    assert pair["features"] == ["401", "403"]               # 塊が無いときは RT の昇順
    assert pair["names_compatible"] and pair["mergeable"] is True
    assert pair["delta_ppm"] == pytest.approx(10.2, abs=0.1)


def _interleaved_ds(delta_rt):
    owned = {i for i in range(N) if i % 3 == 0}
    return _make_ds([("1", 500.3000, 8.0, ALL - owned, None, ADDUCT),
                     ("2", 500.3010, 8.0 + delta_rt, owned, None, ADDUCT)], detected=True)


def test_interleaved_mergeable_only_within_rt_gap():
    # ΔRT ≤ 0.15 分（MS-DIAL のアライメントの RT 許容 0.1 分の 1.5 倍）だけまとめる。
    # 0.2 にしないのは、MTBLS6965 の陰イオンの FA 20:3 の異性体 401/402 が 0.20 分離れているから。
    near = Q.find_split_pairs(_interleaved_ds(0.1))
    (pair,) = near["pairs"]
    assert pair["type"] == "interleaved" and pair["mergeable"] is True
    assert near["groups"][0]["mergeable"] is True
    assert near["merge_groups"][0]["mergeable"] is True
    for delta in (0.18, 0.3):
        far = Q.find_split_pairs(_interleaved_ds(delta))
        (pair,) = far["pairs"]
        assert pair["type"] == "interleaved" and pair["mergeable"] is False, delta
        assert far["groups"][0]["mergeable"] is False and far["merge_groups"] == []
    assert far["params"]["max_interleaved_rt_gap"] == 0.15
    widened = Q.find_split_pairs(_interleaved_ds(0.3), max_interleaved_rt_gap=0.5)
    assert widened["pairs"][0]["mergeable"] is True
    # injection_block は ΔRT が大きくてもまとめる
    ds = _make_ds([("1", 650.6, 20.0, ALL - FIRST_TWO, None, ADDUCT),
                   ("2", 650.6, 20.9, FIRST_TWO, None, ADDUCT)])
    assert Q.find_split_pairs(ds)["pairs"][0]["mergeable"] is True


def test_group_aligned_not_mergeable():
    # 群が注入順で交互に並び、受け持ちが群と一致する → 注入順の塊ではない
    owned = {i for i, g in enumerate(GROUPS_ALTERNATE) if g == "B"}
    ds = _make_ds([("10", 500.3000, 8.0, ALL - owned, None, ADDUCT),
                   ("11", 500.3010, 8.3, owned, None, ADDUCT)], groups=GROUPS_ALTERNATE)
    qc = Q.find_split_pairs(ds)
    (pair,) = qc["pairs"]
    assert pair["type"] == "group_aligned" and pair["mergeable"] is False
    assert qc["groups"][0]["mergeable"] is False
    assert qc["summary"]["injection_block_median_delta_rt"] is None


def test_injection_block_wins_over_group():
    # 塊がちょうど 1 つの群と重なっても injection_block（測定順と群が重なる run）
    owned = set(range(21, N))
    ds = _make_ds([("10", 500.3000, 8.0, ALL - owned, None, ADDUCT),
                   ("11", 500.3010, 8.3, owned, None, ADDUCT)])
    (pair,) = Q.find_split_pairs(ds)["pairs"]
    assert pair["type"] == "injection_block" and pair["mergeable"] is True


def test_incompatible_names_not_mergeable():
    ds = _make_ds([("1", 650.6000, 20.0, ALL - FIRST_TWO, "CE 18:2", ADDUCT),
                   ("2", 650.6050, 20.6, FIRST_TWO, "CE 20:4", ADDUCT)])
    (pair,) = Q.find_split_pairs(ds)["pairs"]
    assert pair["type"] == "injection_block"
    assert pair["names_compatible"] is False and pair["mergeable"] is False


def test_same_sum_composition_is_compatible():
    pytest.importorskip("pygoslin")
    assert Q.names_compatible("PC 16:0_18:1", "PC 34:1")
    assert Q.names_compatible("Glucose", "glucose ")
    assert Q.names_compatible(None, "PC 34:1") and Q.names_compatible("Unknown", "PC 34:1")
    assert not Q.names_compatible("PC 34:1", "PC 34:2")
    assert not Q.names_compatible("Glucose", "Fructose")


def test_three_spot_chain_is_one_group():
    # 背景を 0 にする: どの 2 スポットの組でも、3 つ目のスポットが受け持つ試料は「両方とも
    # 信号なし」になり、補い合いの分母に入らない。背景が 2 スポットで同じ程度だと、その試料は
    # 受け持ちが決まらず組にならない（隣り合う組だけが残る。報告の懸念に記した）。
    early, middle = FIRST_TWO, set(range(2, 20))
    late = ALL - early - middle
    ds = _make_ds([("a", 400.2000, 5.0, early, None, ADDUCT),
                   ("b", 400.2010, 5.3, middle, None, ADDUCT),
                   ("c", 400.2020, 5.6, late, None, ADDUCT),
                   ("z", 900.0000, 9.0, ALL, None, ADDUCT)], detected=True, background=0.0)
    qc = Q.find_split_pairs(ds)
    assert len(qc["pairs"]) == 3
    (group,) = qc["groups"]
    assert sorted(group["features"]) == ["a", "b", "c"]
    assert group["features"] == ["a", "b", "c"]            # ID の小さい順
    assert all(p["mergeable"] for p in qc["pairs"])
    # b の最も近い相手は a と c の同点（0.3 分）→ b は相手なし。a と c の最も近い相手は b
    # なので、互いに最も近い組は無く、何もまとめない
    assert group["mergeable"] is False and qc["merge_groups"] == []
    assert {p["merge_reason"] for p in qc["pairs"]} == {"not_mutual_nearest"}
    partners = Q.split_partners(qc)
    assert set(partners) == {"a", "b", "c"}
    assert partners["a"]["partners"] == [] and partners["a"]["group_id"] == group["group_id"]
    assert partners["a"]["component"] == ["b", "c"]
    assert partners["b"]["mergeable"] is False and partners["b"]["type"] == "injection_block"


def _chain(names):
    """RT が 0.8 分ずつ離れた鎖（隣り合う組だけが ΔRT ≤ 1.0 分）。すべて injection_block。"""
    owned = [FIRST_TWO, set(range(2, 11)), set(range(11, 21)), set(range(21, N))]
    fids = "ABCD"[:len(names)]
    return _make_ds([(fid, 400.2000 + 0.001 * k, 5.0 + 0.8 * k, owned[k], name, ADDUCT)
                     for k, (fid, name) in enumerate(zip(fids, names))],
                    detected=True, background=0.0)


def test_four_spot_chain_merges_only_the_mergeable_pair():
    # A–B と C–D は名前が両立しない。B–C だけが mergeable → B と C だけをまとめる
    qc = Q.find_split_pairs(_chain(["PC 34:1", "PC 34:2", "PC 34:2", "PC 36:1"]))
    edges = {tuple(sorted(p["features"])): p["mergeable"] for p in qc["pairs"]}
    assert edges == {("A", "B"): False, ("B", "C"): True, ("C", "D"): False}
    (group,) = qc["groups"]
    assert sorted(group["features"]) == ["A", "B", "C", "D"] and group["mergeable"] is False
    (merge,) = qc["merge_groups"]
    assert merge["features"] == ["B", "C"] and merge["mergeable"] is True
    assert qc["summary"]["n_merged_pairs"] == 1
    partners = Q.split_partners(qc)
    assert partners["B"]["partners"] == ["C"] and partners["B"]["mergeable"] is True
    assert partners["A"]["partners"] == [] and partners["A"]["mergeable"] is False
    assert partners["A"]["merge_group_id"] is None and sorted(partners["A"]["component"]) == ["B", "C", "D"]
    ds = _chain(["PC 34:1", "PC 34:2", "PC 34:2", "PC 36:1"])
    merged = Q.merge_matrix_columns(ds.feature_matrix.T, ds.feature_ids,
                                    [m["features"] for m in qc["merge_groups"] if m["mergeable"]])
    assert merged.feature_ids == ["A", "B", "D"]                # 小さい ID（B）を残す
    assert merged.merged_from == {"B": ["B", "C"]}


def test_tie_for_nearest_partner_merges_nothing():
    # B の mergeable な相手 A と C がどちらも 0.8 分 → B は相手を決められず、何もまとめない
    qc = Q.find_split_pairs(_chain([None, None, None]))
    assert [p["mergeable"] for p in qc["pairs"]] == [True, True]
    assert [p["merged"] for p in qc["pairs"]] == [False, False]
    assert qc["merge_groups"] == [] and qc["summary"]["n_merged_pairs"] == 0
    assert not any(v["mergeable"] for v in Q.split_partners(qc).values())


def test_fa_20_3_like_chain_merges_only_mutual_nearest():
    # MTBLS6965 の陰イオンの FA 20:3: 402（前半の n-3/n-6）と 401（後半）は injection_block
    # ΔRT 0.2、401 と 403（前半の n-9）は interleaved ΔRT 0.088、403 と 404 は injection_block
    # ΔRT 0.38。m/z は隣どうしだけが 15 ppm に入る。互いに最も近いのは 401↔403 だけ
    late_even = {i for i in range(10, N) if i % 2 == 0}
    late_odd = {i for i in range(10, N) if i % 2 == 1}
    early = set(range(10))
    ds = _make_ds([("402", 305.2440, 3.000, early, "FA 20:3", "[M-H]1-"),
                   ("401", 305.2464, 3.200, late_even, "FA 20:3", "[M-H]1-"),
                   ("403", 305.2495, 3.288, late_odd, "FA 20:3", "[M-H]1-"),
                   ("404", 305.2522, 3.668, early, "FA 20:3", "[M-H]1-")],
                  detected=True, background=0.0)
    qc = Q.find_split_pairs(ds)
    pairs = {tuple(sorted(p["features"])): p for p in qc["pairs"]}
    assert set(pairs) == {("401", "402"), ("401", "403"), ("403", "404")}
    assert pairs[("401", "402")]["type"] == "injection_block"
    assert pairs[("401", "403")]["type"] == "interleaved"
    assert pairs[("403", "404")]["type"] == "injection_block"
    assert all(p["mergeable"] for p in pairs.values())
    assert [tuple(sorted(p["features"])) for p in qc["pairs"] if p["merged"]] == [("401", "403")]
    assert pairs[("401", "402")]["merge_reason"] == "not_mutual_nearest"
    assert pairs[("401", "403")]["merge_reason"] is None
    (merge,) = qc["merge_groups"]
    assert merge["features"] == ["401", "403"]
    partners = Q.split_partners(qc)
    assert partners["401"]["partners"] == ["403"] and partners["401"]["mergeable"] is True
    assert partners["402"]["partners"] == [] and partners["402"]["mergeable"] is False
    assert sorted(partners["402"]["component"]) == ["401", "403", "404"]


def test_far_mz_or_rt_not_paired():
    base = 500.0
    cases = [
        (base * (1 + 20e-6), 8.6, ADDUCT),      # 20 ppm 離れている
        (base, 9.6, ADDUCT),                     # ΔRT 1.6 分
        (base, 8.03, ADDUCT),                    # ΔRT 0.03 分（同じピークの揺れ）
        (base, 8.6, "[M+Na]1+"),                 # 付加体が違う
    ]
    for mz, rt, adduct in cases:
        ds = _make_ds([("1", base, 8.0, ALL - FIRST_TWO, None, ADDUCT),
                       ("2", mz, rt, FIRST_TWO, None, adduct)])
        assert Q.find_split_pairs(ds)["pairs"] == [], (mz, rt, adduct)


def test_both_detected_above_threshold_not_split():
    owned = set(range(0, 21))
    ds = _make_ds([("1", 500.0, 8.0, ALL - owned, None, ADDUCT),
                   ("2", 500.001, 8.4, owned, None, ADDUCT)], detected=True)
    assert len(Q.find_split_pairs(ds)["pairs"]) == 1
    # 両方で検出された試料が 40%（17/42）→ 割れではなく 2 つの物質
    both = set(range(21, 38))
    ds.detected_mask[1, list(both)] = True
    assert Q.find_split_pairs(ds)["pairs"] == []


def test_undetected_samples_are_owned_by_value():
    # 検出が疎な run（MTBLS6965 の陽イオンの CE 1265/1266）: 塊の側は注入 1 でしか検出されず、
    # 注入 2 は gap-fill だが値は 3 倍以上 → 値で受け持ちを補い、塊は [1, 2] のまま
    ds = _make_ds([("1265", 650.6000, 20.00, ALL - FIRST_TWO, None, ADDUCT),
                   ("1266", 650.6050, 20.61, FIRST_TWO, None, ADDUCT)], detected=True)
    ds.detected_mask[1, 1] = False
    ds.detected_mask[0, 10:30] = False
    (pair,) = Q.find_split_pairs(ds)["pairs"]
    assert pair["block"] == {"spot": "1266", "injections": [1, 2]}
    assert pair["n_owned"] == [40, 2] and pair["both_detected_fraction"] == 0.0


def test_sparse_detection_alone_is_not_a_split():
    # どちらも数試料でしか検出されず、ほかの試料では値が同程度（無関係な 2 スポット）
    ds = _make_ds([("1", 500.0, 8.0, {5}, None, ADDUCT),
                   ("2", 500.001, 8.6, {30, 31}, None, ADDUCT)], detected=True)
    assert Q.find_split_pairs(ds)["pairs"] == []


def test_without_detection_needs_threefold_dominance():
    owned = FIRST_TWO
    ds = _make_ds([("1", 500.0, 8.0, ALL - owned, None, ADDUCT),
                   ("2", 500.001, 8.6, owned, None, ADDUCT)])
    # 受け持たない側の値を本物の半分にする（gap-fill が隣のピークを拾った形）→ 補い合わない
    ds.feature_matrix[1, 2:] = ds.feature_matrix[0, 2:] / 2
    assert Q.find_split_pairs(ds)["pairs"] == []


def test_filename_order_when_metadata_lacks_injection_order():
    ds = _make_ds([("1", 650.6, 20.0, ALL - {40, 41}, None, ADDUCT),
                   ("2", 650.6, 20.6, {40, 41}, None, ADDUCT)], rows=False)
    qc = Q.find_split_pairs(ds)
    assert qc["injection_order_source"] == "filename"
    assert qc["pairs"][0]["block"] == {"spot": "2", "injections": [41, 42]}
    ds.sample_metadata_rows = _rows(GROUPS_HALVES, source="mztab")
    assert Q.find_split_pairs(ds)["injection_order_source"] == "metadata"
    # 注入順を引数で渡せば、それを使う（逆順なら塊は先頭に来る）
    order = {name: N - i for i, name in enumerate(NAMES)}
    qc = Q.find_split_pairs(ds, injection_order=order)
    assert qc["injection_order_source"] == "metadata"
    assert qc["pairs"][0]["block"] == {"spot": "2", "injections": [1, 2]}


def test_single_sample_spot_is_not_an_injection_block():
    # 片方が 1 試料（注入の途中・先頭）だけを受け持つ → 塊にせず interleaved（ΔRT 0.5 > 0.15 で
    # まとめない）
    for owned in ({20}, {0}):
        ds = _make_ds([("1", 500.3000, 8.0, ALL - owned, None, ADDUCT),
                       ("2", 500.3010, 8.5, owned, None, ADDUCT)], detected=True)
        (pair,) = Q.find_split_pairs(ds)["pairs"]
        assert pair["type"] == "interleaved" and pair["block"] is None, owned
        assert pair["mergeable"] is False
        assert Q.find_split_pairs(ds, min_block_size=1)["pairs"][0]["type"] == "injection_block"
    # 2 試料なら塊（CE の注入 1–2）
    ds = _make_ds([("1", 500.3000, 8.0, ALL - FIRST_TWO, None, ADDUCT),
                   ("2", 500.3010, 8.5, FIRST_TWO, None, ADDUCT)], detected=True)
    assert Q.find_split_pairs(ds)["pairs"][0]["type"] == "injection_block"


def test_unknown_adduct_is_not_paired():
    ds = _make_ds([("1", 650.6000, 20.0, ALL - FIRST_TWO, None, None),
                   ("2", 650.6050, 20.6, FIRST_TWO, None, None)])
    assert Q.find_split_pairs(ds)["pairs"] == []


def test_injection_order_argument_must_cover_every_sample():
    ds = _make_ds([("1", 650.6, 20.0, ALL - FIRST_TWO, None, ADDUCT),
                   ("2", 650.6, 20.6, FIRST_TWO, None, ADDUCT)])
    order = {name: i for i, name in enumerate(NAMES[:-1], start=1)}
    with pytest.raises(DomainError) as info:
        Q.find_split_pairs(ds, injection_order=order)
    assert info.value.details["missing"] == ["S42"]


def test_bad_parameters_are_domain_errors():
    ds = _make_ds([("1", 500.0, 8.0, ALL, None, ADDUCT)])
    for kwargs in ({"ppm": 0}, {"min_rt_gap": 1.0, "max_rt_gap": 0.5},
                   {"max_both_detected": 1.5}, {"max_interleaved_rt_gap": -0.1},
                   {"min_block_size": 0}, {"min_block_size": 1.5}):
        with pytest.raises(DomainError):
            Q.find_split_pairs(ds, **kwargs)


def test_merge_matrix_columns_takes_per_sample_max():
    matrix = np.array([[1.0, 50.0, 7.0, 3.0],
                       [9.0, 2.0, np.nan, 3.0],
                       [np.nan, np.nan, 1.0, 3.0]])           # (試料, 特徴)
    detected = np.array([[True, False, True, True],
                         [False, True, False, True],
                         [False, False, True, True]])
    merged = Q.merge_matrix_columns(matrix, ["a", "b", "c", "d"], [["c", "a"]],
                                    detected=detected)
    assert merged.feature_ids == ["a", "b", "d"]              # 最も小さい ID の位置に 1 列
    assert merged.merged_from == {"a": ["a", "c"]}
    np.testing.assert_array_equal(merged.matrix[:, 0], [7.0, 9.0, 1.0])
    np.testing.assert_array_equal(merged.detected[:, 0], [True, False, True])
    np.testing.assert_array_equal(merged.matrix[:, [1, 2]], matrix[:, [1, 3]])
    # 整数の ID は数の順（"9" < "10"）。渡した順によらない
    numeric = Q.merge_matrix_columns(matrix, ["10", "b", "9", "d"], [["10", "9"]])
    assert numeric.feature_ids == ["b", "9", "d"] and numeric.merged_from == {"9": ["9", "10"]}
    # 行列に無い特徴は飛ばし、1 つしか残らない群はまとめない
    same = Q.merge_matrix_columns(matrix, ["a", "b", "c", "d"], [["a", "zz"]])
    assert same.feature_ids == ["a", "b", "c", "d"] and same.merged_from == {}
    assert same.detected is None
    with pytest.raises(ValueError):
        Q.merge_matrix_columns(matrix, ["a", "b", "c", "d"], [["a", "b"], ["b", "c"]])
