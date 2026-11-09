"""pipeline が作った差次的結果を、対話セッションへ引き渡す経路（spec 2026-10-19 §11）。

`load_dataset(pipeline_path=...)` は run の dataset を作り直すので、`dataset_id` は
run 側と必ず食い違う。ここで縛るのは:

- 指紋（データ・群）が一致する差次的結果だけを、**result_id を付け替えずに**
  `ds.results` へ戻す（pipeline の TSV のメタ行と同じ ID で引けるように）。
- 戻せなかったもの（旧 run・データ違い・群違い）は理由つきで報告する。

この引き渡しは pipeline が下流を持たなくなる改修（D 段）で消える。pipeline はもう
差次的結果を作らない（pipeline-run.v3）ので、実 run からの引き継ぎの検証は無い。
残るのは旧 v1 の run を読むときの規則（戻さなかった理由と `is_current` の一致）だけ。
"""
from __future__ import annotations

import copy

from metabolomix.analysis import result_state
from metabolomix.pipeline import session_handoff


# ---------- 戻さなかった理由の案内と、is_current との一致 ----------

def test_every_rejection_reason_tells_the_next_step():
    """理由だけ告げて止めると、利用者は何をすれば統計が使えるのか分からない。"""
    from metabolomix.tools.mztab_tools import _handoff_lines

    reasons = ["legacy_result_without_dataset_hash", "dataset_mismatch", "group_mismatch"]
    report = {"restored": [], "rejected": [
        {"comparison_id": f"c{i}", "result_id": f"res_{i}", "reason": reason}
        for i, reason in enumerate(reasons)]}
    notes = {line.split(":", 1)[0].strip(" -"): line
             for line in _handoff_lines(report) if line.startswith("  - ")}
    assert set(notes) == set(reasons)
    # データが違う: session で計算し直す。
    assert "dataset_preprocess" in notes["dataset_mismatch"]
    assert "dataset_compare" in notes["dataset_mismatch"]
    # 群が違う: シートか session のメタデータを確かめる、または session で計算し直す。
    assert "dataset_set_sample_metadata" in notes["group_mismatch"]
    assert "dataset_compare" in notes["group_mismatch"]


def _handed_off(ds, **overrides) -> dict:
    provenance = {"result_id": "res_x", "kind": "differential", "dataset_id": "ds_run",
                  "parent_ids": ["pp_run"],
                  "dataset_hash": result_state.dataset_fingerprint(ds),
                  "group_fingerprint": result_state.group_fingerprint(ds)}
    provenance.update(overrides)
    return {"kind": "two_group", "provenance": provenance}


def test_rejection_reasons_agree_with_is_current():
    """読み込み時に戻した結果は使えて、戻さなかった結果は使えない（2 つの規則がずれない）。"""
    ds, _ = _dataset_for_rule_agreement()
    cases = {
        None: _handed_off(ds),
        "legacy_result_without_dataset_hash": _handed_off(ds, dataset_hash=None),
        "dataset_mismatch": _handed_off(ds, dataset_hash="0" * 64),
        "group_mismatch": _handed_off(ds, group_fingerprint="1" * 64),
    }
    for expected, result in cases.items():
        provenance = result["provenance"]
        assert session_handoff._handoff_rejection(ds, provenance) == expected
        restored = copy.deepcopy(result)
        restored["provenance"]["handoff"] = {"pipeline_root": "x", "comparison_id": "c1"}
        assert result_state.is_current(ds, restored) is (expected is None), expected


def _dataset_for_rule_agreement():
    from tests.test_result_state import _dataset_with_groups

    return _dataset_with_groups("control", "treated")
