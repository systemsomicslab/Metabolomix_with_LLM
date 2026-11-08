"""戻り値の直列化（metabolomix.core.serialization）。

`json_safe` は厳密な JSON（NaN / Infinity を含まない）にするための変換、`truncate_lists` は
長い一覧を戻り値から切るための変換。どちらも非破壊。
"""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from metabolomix.core.serialization import json_payload, json_safe, round_floats, truncate_lists


def _raise(name):
    raise ValueError(f"厳密な JSON ではない定数: {name}")


def _strict(text: str):
    return json.loads(text, parse_constant=_raise)


def test_json_safe_rejects_no_constants():
    inf, nan = math.inf, math.nan
    out = _strict(json_payload(json_safe({"a": inf, "b": nan, "c": -inf, "d": 1.5})))
    assert out == {"a": "inf", "b": None, "c": "-inf", "d": 1.5}


def test_json_safe_numpy_scalars_arrays_and_sets():
    obj = {
        "f": np.float64(2.5), "f_inf": np.float64(np.inf), "f_nan": np.float32(np.nan),
        "i": np.int64(7), "b": np.bool_(True), "arr": np.array([[1.0, np.nan], [np.inf, 2.0]]),
        "iarr": np.arange(3), "s": {"b", "a", "c"}, "fs": frozenset({3, 1, 2}), "t": (1, (2, 3)),
        np.int64(4): "numpy key",
    }
    out = json_safe(obj)
    assert type(out["f"]) is float and out["f"] == 2.5
    assert out["f_inf"] == "inf" and out["f_nan"] is None
    assert type(out["i"]) is int and out["i"] == 7
    assert out["b"] is True
    assert out["arr"] == [[1.0, None], ["inf", 2.0]]
    assert out["iarr"] == [0, 1, 2] and all(type(x) is int for x in out["iarr"])
    assert out["s"] == ["a", "b", "c"] and out["fs"] == [1, 2, 3]
    assert out["t"] == [1, [2, 3]]
    assert out[4] == "numpy key"
    # 全体が厳密な JSON になる（allow_nan=False でも通る）
    json.dumps(out, allow_nan=False)
    _strict(json_payload(out))
    # 非破壊
    assert isinstance(obj["arr"], np.ndarray) and math.isnan(obj["f_nan"])


def test_json_safe_keeps_bool_and_none_and_mixed_set():
    assert json_safe(True) is True and json_safe(None) is None and json_safe("x") == "x"
    assert json_safe(3) == 3
    # 並べられない集合も落とさずに list にする（決定的な順）
    out = json_safe({1, "a"})
    assert sorted(map(str, out)) == ["1", "a"] and json_safe({1, "a"}) == out


def test_truncate_lists_adds_n_more():
    obj = {"unit_ids": [f"u{i}" for i in range(25)], "other": list(range(30)),
           "nested": [{"samples": list("abcdefgh")}],
           "examples": ["x"] * 3}
    out = truncate_lists(obj, limit=5)
    assert out["unit_ids"] == ["u0", "u1", "u2", "u3", "u4"] and out["unit_ids_n_more"] == 20
    # 隣に置く（キーの並びで直後）
    keys = list(out)
    assert keys.index("unit_ids_n_more") == keys.index("unit_ids") + 1
    assert out["other"] == list(range(30))          # 指定外のキーは切らない
    assert out["nested"][0]["samples"] == list("abcde")
    assert out["nested"][0]["samples_n_more"] == 3
    assert out["examples"] == ["x"] * 3 and "examples_n_more" not in out   # 短いものはそのまま
    assert len(obj["unit_ids"]) == 25                 # 非破壊


def test_truncate_lists_adds_to_existing_n_more_and_custom_keys():
    # エンジンが既に [:20] で切って n_more を持つ一覧をさらに切っても、数を失わない
    obj = {"subjects": [f"s{i}" for i in range(20)], "subjects_n_more": 7}
    out = truncate_lists(obj, limit=10)
    assert len(out["subjects"]) == 10 and out["subjects_n_more"] == 17
    out = truncate_lists({"foo": [1, 2, 3], "unit_ids": [1, 2, 3]}, limit=1, keys=("foo",))
    assert out == {"foo": [1], "foo_n_more": 2, "unit_ids": [1, 2, 3]}
    with pytest.raises(ValueError):
        truncate_lists({}, limit=-1)


def test_round_floats_unchanged():
    assert round_floats({"a": [1.23456789]}, 3) == {"a": [1.235]}
