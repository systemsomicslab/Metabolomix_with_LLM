"""候補提示ツール。失敗する前に候補を列挙できることが要点。

手動実行 WebUI の「メソッドを選ぶ」画面はこれを叩く。console_plan の封筒を
待たないと候補が分からない設計では、画面を先に描けない。
"""
from __future__ import annotations

import json as _json
from pathlib import Path

from metabolomix.tools.console_tools import console_method_candidates


def _write_param(directory: Path, name: str, ion: str = "Positive") -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    p = directory / name
    p.write_text(f"Ion mode: {ion}\nTarget omics: Lipidomics\n"
                 "Minimum peak height: 1000\n", encoding="ascii")
    return p


def test_lists_a_sibling_candidate_of_the_other_polarity(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    _write_param(tmp_path / "NEG", "d_param_1.txt", ion="Negative")
    parsed = _json.loads(console_method_candidates(str(root), polarity="positive"))
    assert parsed["n_candidates"] == 1
    entry = parsed["candidates"][0]
    assert entry["origin"] == "sibling"
    assert entry["usable"] == "needs_polarity_conversion"
    assert entry["key_params"]["Ion mode"] == "Negative"


def test_reports_where_it_looked_when_nothing_is_found(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    parsed = _json.loads(console_method_candidates(str(root), polarity="positive"))
    assert parsed["n_candidates"] == 0
    assert str(root) in parsed["searched"]


def test_polarity_may_be_omitted(tmp_path):
    """polarity を渡さないと『一致するか』を判定しようがないので、usable は
    direct と誤認できる値を出さず、キーそのものを省く（needs_polarity_conversion
    と同じく key_params が None のときの省略と同じ扱い）。"""
    root = tmp_path / "POS"
    _write_param(root, "own_param_1.txt")
    parsed = _json.loads(console_method_candidates(str(root)))
    assert "usable" not in parsed["candidates"][0]


def test_truncates_the_candidate_list_beyond_the_cap(tmp_path):
    """候補が多いラボ配置（兄弟フォルダ多数）でも戻り値を無制限に太らせない。

    n_candidates は総数のまま、candidates は上限までに切り、切ったときだけ
    truncated が立つ（切っていないときは truncated キー自体が無い）。
    """
    from metabolomix.console.method_file import MAX_REPORTED_CANDIDATES
    root = tmp_path / "POS"
    root.mkdir()
    n_siblings = MAX_REPORTED_CANDIDATES + 2
    for i in range(n_siblings):
        _write_param(tmp_path / f"SIB{i}", f"p_param_{i}.txt")
    total = n_siblings  # root 自身は候補を持たないので兄弟の件数と一致
    parsed = _json.loads(console_method_candidates(str(root), polarity="positive"))
    assert parsed["n_candidates"] == total
    assert len(parsed["candidates"]) == MAX_REPORTED_CANDIDATES
    assert parsed["truncated"] is True


def test_does_not_report_truncated_when_nothing_was_cut(tmp_path):
    root = tmp_path / "POS"
    _write_param(root, "own_param_1.txt")
    parsed = _json.loads(console_method_candidates(str(root), polarity="positive"))
    assert "truncated" not in parsed


def test_missing_dataset_root_is_an_envelope_not_an_exception(tmp_path):
    parsed = _json.loads(console_method_candidates(str(tmp_path / "nope")))
    assert parsed["error"]["code"] == "DATASET_ROOT_NOT_FOUND"
