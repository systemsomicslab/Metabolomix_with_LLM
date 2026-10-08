"""arf/selection.py: 項目の解決とスポット単位への展開。fixture はテストが作る。"""
from __future__ import annotations

import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state
from tests.curation_fixtures import arf2_spot_raw, arf_row, match_result, write_arf, write_arf2

STEM = "AlignmentResult_2026_01_01_00_00_00"
SAMPLES = ["blank_1", "ctrl_1", "ctrl_2", "ko_1", "ko_2", "std_1"]


def _group(heights, gap=()):
    rows = []
    for fid, (name, h) in enumerate(zip(SAMPLES, heights)):
        row = arf_row(file_id=fid, mz=700.0, rt=3.0, height=h, gap_filled=name in gap)
        row[1] = f"20261006_{name}"
        rows.append(row)
    return rows


@pytest.fixture()
def loaded(tmp_path):
    arf = write_arf(tmp_path / f"{STEM}_PeakProperties.arf", [
        _group([10, 5000, 6000, 2000, 2500, 10]),            # 0 PG 34:1（MS/MS あり）
        _group([10, 3000, 3500, 4000, 4200, 10]),            # 1 PG 35:1（MS/MS なし）
        _group([0, 900, 950, 800, 700, 0]),                  # 2 DG 34:1 [M+NH4]+
        _group([0, 300, 320, 280, 260, 0]),                  # 3 DG 34:1 [M+Na]+
    ])
    arf2 = write_arf2(tmp_path / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=0, name="PG 34:1|PG 16:0_18:1", ontology="PG", adduct="[M-H]-",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=1, name="PG 35:1|PG 16:0_19:1", ontology="PG", adduct="[M-H]-"),
        arf2_spot_raw(spot_id=2, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+NH4]+",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=3, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+Na]+",
                      matches=[match_result()]),
    ])
    session_state.session = session_state.AnalysisSession()
    with open(arf, "rb") as handle:
        session_state.session.arf.features = arf_reader.deserialize(handle)
    session_state.session.arf.current_file_path = str(arf)
    session_state.session.arf.class_index = None
    return {"arf": arf, "arf2": arf2}


def _select(loaded, **kw):
    from metabolomix.arf import selection
    kw.setdefault("groups", ["ctrl", "ko"])
    return selection.build_selection(session_state.session.arf, loaded["arf2"], **kw)


def test_build_selection_resolves_items_groups_and_msms(loaded):
    sel = _select(loaded, items=["PG", "DG"])
    assert [s["spot_id"] for s in sel.items[0]["spots"]] == [0, 1]
    assert [g["label"] for g in sel.groups] == ["ctrl", "ko"]
    assert sel.msms[0] is True and sel.msms.get(1, False) is False
    assert sel.catalog_by_id[2]["AdductType"] == "[M+NH4]+"


def test_expand_spots_orders_by_item_then_name_and_keeps_adduct_duplicates(loaded):
    from metabolomix.arf import selection
    spots, dropped = selection.expand_spots(_select(loaded, items=["DG", "PG"]))
    assert [s["spot_id"] for s in spots] == [2, 3, 0, 1]
    assert [s["adduct"] for s in spots[:2]] == ["[M+NH4]+", "[M+Na]+"]
    assert spots[0]["label"] == "DG 16:0_18:1" and spots[0]["item"] == "DG"
    assert dropped == []


def test_expand_spots_require_msms_drops_unsupported(loaded):
    from metabolomix.arf import selection
    spots, dropped = selection.expand_spots(_select(loaded, items=["PG"]), require_msms=True)
    assert [s["spot_id"] for s in spots] == [0]
    assert dropped == ["#1 PG 35:1|PG 16:0_19:1"]


def test_expand_spots_counts_a_spot_once_across_items(loaded):
    from metabolomix.arf import selection
    spots, _ = selection.expand_spots(_select(loaded, items=["PG", "PG 34:1"]))
    assert [s["spot_id"] for s in spots] == [0, 1]


def test_display_name_takes_the_last_candidate():
    from metabolomix.arf.selection import display_name
    assert display_name("PG 34:1|PG 16:0_18:1") == "PG 16:0_18:1"
    assert display_name("no MS2: PC 34:1") == "no MS2: PC 34:1"
