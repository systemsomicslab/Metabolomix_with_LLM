"""新しい知見の候補の外部の印（labels）。強さの順・空と未知の水準・割れの組・件数・rule_hash。"""
import pytest

from metabolomix.analysis.mass_accuracy import (attach_labels, normalize_novel_spec, read_labels,
                                                rule_hash)
from metabolomix.core.atomic_io import DomainError


def _write(tmp_path, pairs):
    path = tmp_path / "labels.tsv"
    path.write_text("feature_id\tlevel\n" + "".join(f"{f}\t{lv}\n" for f, lv in pairs),
                    encoding="utf-8")
    return path


def _spec(path):
    return {"path": str(path), "name": "xenobiotic", "feature_column": "feature_id",
            "label_column": "level", "order": ["A", "B", "C", "D"]}


def test_strongest_level_wins_for_a_feature_with_several_rows(tmp_path):
    path = _write(tmp_path, [("10", "C"), ("10", "A"), ("11", "D")])
    out = read_labels(_spec(path))
    assert out["by_feature"] == {"10": "A", "11": "D"} and out["n_rows"] == 3
    assert len(out["sha256"]) == 64


def test_unknown_level_is_refused_and_empty_level_is_ignored(tmp_path):
    assert read_labels(_spec(_write(tmp_path, [("10", "")])))["by_feature"] == {}
    with pytest.raises(DomainError) as info:
        read_labels(_spec(_write(tmp_path, [("10", "Z")])))
    assert info.value.code == "NOVEL_SPEC_INVALID" and info.value.details["unknown_levels"] == ["Z"]


def test_merged_pair_takes_the_strongest_member_label():
    rows = [{"feature_id": "5", "q": 0.01, "merged_from": ["5", "6"]}]
    attach_labels(rows, {"6": "B"}, name="xenobiotic", order=["A", "B", "C", "D"])
    assert rows[0]["labels.xenobiotic"] == "B"


def test_counts_all_and_top20_and_order_unchanged():
    rows = [{"feature_id": str(i), "q": i / 1000} for i in range(30)]
    before = [r["feature_id"] for r in rows]
    counts = attach_labels(rows, {"0": "A", "25": "B"}, name="x", order=["A", "B", "C", "D"])
    assert [r["feature_id"] for r in rows] == before
    assert counts["all"]["A"] == 1 and counts["all"]["B"] == 1 and counts["all"]["none"] == 28
    assert counts["top20"]["A"] == 1 and counts["top20"]["B"] == 0
    assert rows[1]["labels.x"] is None


def test_labels_do_not_change_rule_hash():
    base = {"statistics": {"result_id": "r"}, "q_max": 0.05}
    with_labels = {**base, "labels": {"path": "x.tsv", "name": "xenobiotic"}}
    assert rule_hash(normalize_novel_spec(base)) == rule_hash(normalize_novel_spec(with_labels))


def test_labels_spec_is_validated_and_defaulted():
    base = {"statistics": {"result_id": "r"}}
    got = normalize_novel_spec({**base, "labels": {"path": "x.tsv", "name": "xeno"}})["labels"]
    assert got["feature_column"] == "feature_id" and got["label_column"] == "level"
    assert got["order"] == ["A", "B", "C", "D"]
    for bad in ({"path": "x.tsv"}, {"name": "x"}, {"path": "x.tsv", "name": "a b"},
                {"path": "x.tsv", "name": "x", "order": ["A", "A"]},
                {"path": "x.tsv", "name": "x", "extra": 1}):
        with pytest.raises(DomainError) as info:
            normalize_novel_spec({**base, "labels": bad})
        assert info.value.code == "NOVEL_SPEC_INVALID"
