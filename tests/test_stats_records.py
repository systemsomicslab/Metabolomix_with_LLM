"""統計の工程の記録（metabolomix.analysis.stats_records）と、上書きしない JSON の書き込み
（metabolomix.core.atomic_io.write_new_json）。fixture は tmp に合成する。"""
from __future__ import annotations

import json
import math
import re

import numpy as np
import pytest

from metabolomix.analysis import stats_records
from metabolomix.core.atomic_io import DomainError, write_new_json
from metabolomix.quant import records as quant_records


def _raise(name):
    raise ValueError(f"厳密な JSON ではない定数: {name}")


def _write(directory, step="compare", **kw):
    args = dict(inputs={"dataset_hash": "abc", "matrix_id": "mx_1"}, parent_ids=["s_parent"],
                settings={"q_threshold": 0.05}, summary={"n_tested": 3})
    args.update(kw)
    return stats_records.write_stats_record(directory, step, **args)


def test_stats_record_append_only_and_strict_json(tmp_path, monkeypatch):
    reports = tmp_path / "reports"
    r = _write(reports, summary={"df_prior": math.inf, "s2": math.nan, "n": np.int64(4),
                                 "v": np.array([1.0, np.inf])})
    assert re.fullmatch(r"s_\d{8}T\d{6}Z_[0-9a-f]{8}", r["record_id"])
    path = reports / "stats" / f"{r['record_id']}.json"
    assert r["path"] == str(path) and path.is_file()
    data = json.loads(path.read_text(encoding="utf-8"), parse_constant=_raise)
    assert data["schema"] == stats_records.RECORD_SCHEMA == "lipidmix.stats-record.v1"
    assert data["step"] == "compare" and data["parent_ids"] == ["s_parent"]
    assert data["inputs"] == {"dataset_hash": "abc", "matrix_id": "mx_1"}
    assert data["settings"] == {"q_threshold": 0.05}
    assert data["summary"] == {"df_prior": "inf", "s2": None, "n": 4, "v": [1.0, "inf"]}
    assert "created_at" in data and "path" not in data   # 置き場所はファイルに書かない
    assert {k: v for k, v in r.items() if k != "path"} == data
    # 追記専用: 同じ ID の記録は上書きしない
    monkeypatch.setattr(stats_records, "_new_record_id", lambda: r["record_id"])
    with pytest.raises(DomainError) as exc:
        _write(reports, summary={"other": 1})
    assert exc.value.code == "STATS_RECORD_EXISTS"
    assert json.loads(path.read_text(encoding="utf-8"))["summary"]["df_prior"] == "inf"
    assert [p.name for p in (reports / "stats").iterdir()] == [path.name]   # 一時ファイルを残さない


def test_stats_record_steps_and_distinct_ids(tmp_path):
    ids = {_write(tmp_path, step)["record_id"] for step in ("preprocess", "compare", "pca")}
    assert len(ids) == 3
    with pytest.raises(DomainError) as exc:
        _write(tmp_path, "normalize")
    assert exc.value.code == "STATS_RECORD_STEP_INVALID"


def test_write_new_json_refuses_overwrite_and_strict(tmp_path):
    path = tmp_path / "sub" / "a.json"
    write_new_json(path, {"x": 1})                   # 親フォルダも作る
    assert json.loads(path.read_text(encoding="utf-8")) == {"x": 1}
    with pytest.raises(FileExistsError):
        write_new_json(path, {"x": 2})
    assert json.loads(path.read_text(encoding="utf-8")) == {"x": 1}
    with pytest.raises(ValueError):                  # allow_nan=False
        write_new_json(tmp_path / "nan.json", {"x": math.nan})
    assert not (tmp_path / "nan.json").exists()
    assert sorted(p.name for p in (tmp_path / "sub").iterdir()) == ["a.json"]


def test_write_new_json_without_hard_links_still_refuses(tmp_path, monkeypatch):
    import os

    def no_link(src, dst):
        raise OSError(1, "hard links not supported")

    monkeypatch.setattr(os, "link", no_link)
    path = tmp_path / "a.json"
    write_new_json(path, {"x": 1})
    with pytest.raises(FileExistsError):
        write_new_json(path, {"x": 2})
    assert json.loads(path.read_text(encoding="utf-8")) == {"x": 1}
    assert [p.name for p in tmp_path.iterdir()] == ["a.json"]


def test_quant_records_still_refuse_overwrite(tmp_path, monkeypatch):
    """write_new_json へ移した後も、定量の記録は QUANT_RECORD_EXISTS で止まる。"""
    monkeypatch.setattr(quant_records, "_new_record_id", lambda: "q_fixed")
    args = dict(inputs={"a": "x"}, parent_ids=[], settings={}, summary={})
    quant_records.write_record(tmp_path, "normalize", **args)
    with pytest.raises(DomainError) as exc:
        quant_records.write_record(tmp_path, "normalize", **args)
    assert exc.value.code == "QUANT_RECORD_EXISTS"
    path = quant_records.record_dir(tmp_path) / "q_fixed.json"
    with pytest.raises(DomainError) as exc:
        quant_records.write_audit(tmp_path, "q_fixed", {"a": 1})
        quant_records.write_audit(tmp_path, "q_fixed", {"a": 2})
    assert exc.value.code == "QUANT_RECORD_EXISTS"
    assert json.loads(path.read_text(encoding="utf-8"))["record_id"] == "q_fixed"
    # 記録の書き込みは atomic_io の 1 か所（quant に独自の排他作成を残さない）
    assert not hasattr(quant_records, "_create_exclusive_direct")
