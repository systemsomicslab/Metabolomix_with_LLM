"""`load_dataset(overview=True)` が最初の概観を返す（A。`overview` 欄・既定の `False` では載せない・
失敗しても読み込みは成功）。概観の本体は `dataset_overview` で、読み込みの `next` の最初の一手がそれ。"""
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
    out = _load(directory=str(_folder(tmp_path / "a")), overview=True)
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
    # 概観を載せたら、次の一手から dataset_overview() を外す（もう手元にある）。
    assert not any(step.startswith("dataset_overview(") for step in out["next"])


def test_load_dataset_default_has_no_overview(tmp_path, fresh_state):
    """既定（`overview=False`）は概観を計算も記録もしない。次の一手の最初が `dataset_overview()`。"""
    out = _load(directory=str(_folder(tmp_path / "a")))
    assert out["status"] == "success" and "overview" not in out
    assert not session_state.session.dataset.results
    assert out["next"][0].startswith("dataset_overview(")
    assert out["assay_kind"] == "unknown"


def test_load_dataset_overview_false_keeps_old_shape(tmp_path, fresh_state):
    out = _load(directory=str(_folder(tmp_path / "a")), overview=False)
    assert out["status"] == "success" and "overview" not in out
    assert not session_state.session.dataset.results


def test_load_dataset_overview_must_be_bool(tmp_path, fresh_state):
    out = _load(directory=str(_folder(tmp_path / "a")), overview="false")
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"


def test_load_dataset_incomplete_overview_compares_with_the_mark(tmp_path, fresh_state, monkeypatch):
    real = opening.open_dataset

    def exploratory(**kwargs):
        opened = real(**kwargs)
        opened.dataset.exploratory_only = True
        return opened

    monkeypatch.setattr(opening, "open_dataset", exploratory)
    out = _load(directory=str(_folder(tmp_path / "a")), overview=True)
    assert out["status"] == "success"
    assert "compare" in out["overview"]["result_ids"]
    assert out["overview"]["species"]["facts"][0].startswith("S1 探索専用")


def test_load_dataset_survives_overview_failure(tmp_path, fresh_state, monkeypatch):
    def boom(holder, **kwargs):
        raise RuntimeError("synthetic failure")

    monkeypatch.setattr(overview_tools, "build_overview_payload", boom)
    out = _load(directory=str(_folder(tmp_path / "a")), overview=True)
    assert out["status"] == "success" and out["dataset_handle"].startswith("dh_")
    assert out["overview"]["status"] == "unavailable"
    assert out["overview"]["error"] == {"code": "RuntimeError", "message": "synthetic failure"}
    assert "dataset_overview()" in out["overview"]["next"]
