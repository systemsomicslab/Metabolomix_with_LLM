"""`load_dataset` が最初の概観を返す（A。`overview` 欄・`overview=False`・失敗しても読み込みは成功）。"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import session_state
from metabolomix.dataset import opening
from metabolomix.tools import handle_tools, overview_tools
from tests.test_load_dataset import _folder, fresh_state  # noqa: F401  (fixture)


def _load(**kwargs) -> dict:
    return json.loads(handle_tools.load_dataset(**kwargs))


def test_load_dataset_includes_overview(tmp_path, fresh_state):
    out = _load(directory=str(_folder(tmp_path / "a")))
    assert out["status"] == "success"
    ov = out["overview"]
    assert ov["status"] == "success" and ov["overview_id"].startswith("ovw_")
    assert "counts" not in ov                         # summary と重なるので省く
    for key in ("quality", "global", "classes", "species"):
        assert "next" in ov[key]
    ds = session_state.session.dataset
    assert ds.results[ov["overview_id"]]["kind"] == "overview"
    assert ds.last_pca is None and ds.last_differential is None
    assert list(out).index("overview") > list(out).index("summary")


def test_load_dataset_overview_false_keeps_old_shape(tmp_path, fresh_state):
    out = _load(directory=str(_folder(tmp_path / "a")), overview=False)
    assert out["status"] == "success" and "overview" not in out
    assert not session_state.session.dataset.results


def test_load_dataset_overview_must_be_bool(tmp_path, fresh_state):
    out = _load(directory=str(_folder(tmp_path / "a")), overview="false")
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"


def test_load_dataset_incomplete_overview_has_no_comparisons(tmp_path, fresh_state, monkeypatch):
    real = opening.open_dataset

    def exploratory(**kwargs):
        opened = real(**kwargs)
        opened.dataset.exploratory_only = True
        return opened

    monkeypatch.setattr(opening, "open_dataset", exploratory)
    out = _load(directory=str(_folder(tmp_path / "a")))
    assert out["status"] == "success"
    assert "探索専用" in out["overview"]["species"]["omitted"]
    assert "compare" not in out["overview"]["result_ids"]


def test_load_dataset_survives_overview_failure(tmp_path, fresh_state, monkeypatch):
    def boom(holder, **kwargs):
        raise RuntimeError("synthetic failure")

    monkeypatch.setattr(overview_tools, "build_overview_payload", boom)
    out = _load(directory=str(_folder(tmp_path / "a")))
    assert out["status"] == "success" and out["dataset_handle"].startswith("dh_")
    assert out["overview"]["status"] == "unavailable"
    assert out["overview"]["error"] == {"code": "RuntimeError", "message": "synthetic failure"}
    assert "dataset_overview()" in out["overview"]["next"]
