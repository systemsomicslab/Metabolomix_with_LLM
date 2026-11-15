"""再解析のツールの行の上限・効かない設定の notes・`dataset_alignment_qc` の結果の JSON と
ハンドルに持つ最後の結果・v1 の外部の表と correlate の戻り値。fixture はすべてテスト自身が作る。
"""
from __future__ import annotations

import hashlib
import json
import random

import pytest
from scipy import stats

from metabolomix.analysis import alignment_qc, mass_accuracy
from metabolomix.core import handles, session_state
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import reanalysis_tools as R
from metabolomix.tools import sample_tools
from tests.external_table_fixtures import write_mwtab
from tests.reanalysis_fixtures import SAMPLES
from tests.test_reanalysis_tools import (CE_KEY, SPLIT_CHANGED, SPLIT_FEATURES, _compare_on, _groups_tsv,
                                         _split_handle, _tsv_rows)


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()


# ---------- 結果の JSON は行の範囲で読める（B12） ----------

def test_json_lines_breaks_long_values_and_stays_parseable():
    payload = {"b": [{"pair_id": f"P{i}", "features": [str(i), str(i + 1)], "mz": [500.0, 500.001]}
                     for i in range(300)],
               "a": list(range(2000)), "small": {"x": 1}, "text": "短い"}
    text = R.result_bytes(payload).decode("utf-8")
    assert json.loads(text) == payload
    lines = text.splitlines()
    assert len(lines) > 300 and max(len(line) for line in lines) <= R.JSON_LINE_WIDTH
    assert lines[0] == "{" and lines[-1] == "}"
    assert '"small":{"x":1},' in lines                       # 短い値は 1 行に詰める
    pair_lines = [line for line in lines if line.startswith('{"features"')]
    assert len(pair_lines) == 300                            # 長いリストの dict は 1 要素 1 行
    assert R.result_bytes(payload) == R.result_bytes(dict(reversed(list(payload.items()))))
    # 短い結果は 1 行（キーは昇順）
    assert R.result_bytes({"b": 1, "a": None}) == b'{"a":null,"b":1}\n'
    # 辞書のキーが文字列でなくても json.dumps と同じ規則
    wide = {1: "x" * 600, True: "y" * 600}
    assert json.loads(R.result_bytes(wide)) == json.loads(json.dumps(wide))


def _max_line(payload) -> int:
    text = R.result_bytes(payload).decode("utf-8")
    assert json.loads(text) == json.loads(json.dumps(payload))
    return max(len(line) for line in text.splitlines())


def test_json_lines_count_the_key_and_the_trailing_comma_in_the_width():
    width = R.JSON_LINE_WIDTH
    # 詰めた値（996 字の文字列 1 つのリスト）に `"key":` が付く行
    assert _max_line({"key": ["a" * 996], "pad": "x" * 10}) <= width
    # ちょうど 1,000 字に詰めた辞書の要素に区切りの `,` が付く行
    item = {"v": "b" * (width - len('{"v":""}'))}
    assert len(json.dumps(item, separators=(",", ":"))) == width
    assert _max_line([item, item]) <= width
    # 数や文字列だけのリストを 1,000 字ちょうどまで詰めた行に `,` が付く行
    assert _max_line({"k": ["a" * 498, "b" * 497, "c" * 10, "d" * 600]}) <= width
    # 入れ子の形を変えても行は上限以下（分けられない 1 つの値が上限を超えない限り）
    rng = random.Random(7)

    def node(depth):
        kind = rng.random()
        if depth > 2 or kind < 0.35:
            return rng.choice([rng.randint(-10**6, 10**6), "s" * rng.randint(0, 900), None,
                               rng.random()])
        if kind < 0.65:
            return [node(depth + 1) for _ in range(rng.randint(0, 20))]
        return {f"k{rng.randint(0, 10**5)}" + "x" * rng.randint(0, 40): node(depth + 1)
                for _ in range(rng.randint(0, 10))}

    for _ in range(60):
        assert _max_line({"root": node(0), "nums": list(range(rng.randint(0, 2000)))}) <= width
    # 分けられない 1 つの値（上限より長い文字列）の行だけは超える
    assert _max_line({"long": "z" * (width + 10), "x": list(range(500))}) > width


def test_alignment_qc_output_is_line_broken_and_matches_the_kept_result(tmp_path):
    handle = _split_handle(tmp_path)
    target = tmp_path / "qc.json"
    out = json.loads(R.dataset_alignment_qc(output_path=str(target), dataset_handle=handle))
    assert out["kept_on_handle"] is True
    raw = target.read_bytes()
    written = json.loads(raw)
    assert written["schema"] == "alignment-qc.v1"
    if len(json.dumps(written, separators=(",", ":"))) > R.JSON_LINE_WIDTH:
        assert raw.count(b"\n") > 1                        # 1 行の JSON にしない
    assert max(len(line) for line in raw.decode("utf-8").splitlines()) <= R.JSON_LINE_WIDTH
    ds = handles.registry.get(handle).dataset
    kept, sha256 = alignment_qc.remembered(ds)
    assert kept == written and sha256 == hashlib.sha256(raw).hexdigest()


# ---------- dataset_alignment_qc の group_aligned は先頭 50 組（C7） ----------

def _synthetic_qc(n_group_aligned):
    pairs = [{"pair_id": f"P{i}", "type": "group_aligned", "features": [f"{2 * i}", f"{2 * i + 1}"],
              "mz": [500.0 + i, 500.0 + i], "rt": [1.0, 1.5], "delta_rt": 0.5,
              "adduct": "[M+H]+", "names": [None, None], "mergeable": False}
             for i in range(n_group_aligned)]
    return {"schema": "alignment-qc.v1", "dataset_fingerprint": "synthetic",
            "params": {"signal": "values"}, "injection_order_source": "filename",
            "pairs": pairs, "groups": [], "merge_groups": [],
            "summary": {"n_pairs": len(pairs), "by_type": {"group_aligned": len(pairs)}}}


def test_alignment_qc_lists_the_first_50_group_aligned_pairs_and_the_count(tmp_path, monkeypatch):
    handle = _split_handle(tmp_path)
    monkeypatch.setattr(alignment_qc, "find_split_pairs", lambda ds, **kw: _synthetic_qc(60))
    out = json.loads(R.dataset_alignment_qc(dataset_handle=handle))
    assert out["n_group_aligned"] == 60 and out["group_aligned_truncated"] is True
    rows = _tsv_rows(out["group_aligned"])
    assert len(rows) == R.MAX_INLINE_ROWS and rows[0]["pair_id"] == "P0"
    assert "output_path" in out["notes"][0] and "60" in out["notes"][0]
    monkeypatch.setattr(alignment_qc, "find_split_pairs", lambda ds, **kw: _synthetic_qc(3))
    few = json.loads(R.dataset_alignment_qc(dataset_handle=handle))
    assert few["n_group_aligned"] == 3 and few["group_aligned_truncated"] is False
    assert len(_tsv_rows(few["group_aligned"])) == 3 and "notes" not in few


def test_alignment_qc_without_group_aligned_has_no_list(tmp_path):
    out = json.loads(R.dataset_alignment_qc(dataset_handle=_split_handle(tmp_path)))
    assert "group_aligned" not in out and "n_group_aligned" not in out


# ---------- 最後の割れの QC をハンドルに持つ（C18） ----------

def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_last_alignment_qc_is_used_when_the_path_is_omitted(tmp_path):
    handle = _split_handle(tmp_path, SPLIT_CHANGED)
    diff = _compare_on(handle)
    fresh = json.loads(T.dataset_feature_search(feature_ids=["30"], dataset_handle=handle))
    assert fresh["merge_split"]["source"] == "find_split_pairs"
    assert fresh["merge_split"]["n_groups"] == 1
    claim = {"claim_id": "C9", "target": {"inchikey": CE_KEY, "expected_mz": 650.6, "adduct": "[M+H]+",
                                          "rt_min": 20.3, "rt_tol_min": 0.5}}
    claim_before = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert claim_before["default_selection"]["merged_from"] == ["30", "31"]

    # 既定でない閾値（ppm=1 では 30/31 の 6 ppm の組を拾わない）で回し、その結果が使われる
    qc_path = tmp_path / "qc_strict.json"
    out = json.loads(R.dataset_alignment_qc(ppm=1.0, output_path=str(qc_path), dataset_handle=handle))
    assert out["kept_on_handle"] is True and out["summary"]["n_pairs"] == 0
    search = json.loads(T.dataset_feature_search(feature_ids=["30"], dataset_handle=handle))
    assert search["merge_split"]["source"] == "last_alignment_qc"
    assert search["merge_split"]["alignment_qc_sha256"] == _sha(qc_path)
    assert search["merge_split"]["n_groups"] == 0
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    novel = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    assert novel["merge_split"]["source"] == "last_alignment_qc"
    assert novel["merge_split"]["n_groups"] == 0
    # パスを渡せばそのファイル（同じ中身なので同じ sha256）
    by_file = json.loads(T.dataset_feature_search(feature_ids=["30"], alignment_qc_path=str(qc_path),
                                                  dataset_handle=handle))
    assert by_file["merge_split"]["source"] == "alignment_qc"
    assert by_file["merge_split"]["alignment_qc_sha256"] == _sha(qc_path)

    # dataset_analysis_run（v1）は spec だけで再現する: 最後の結果を使わず、その場で探す
    run_spec = {"groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
                "comparison": {"reference": "HC", "test": "RA"},
                "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}}]}
    run = json.loads(R.dataset_analysis_run(run_spec, dataset_handle=handle))
    assert run["merge_split"]["source"] == "find_split_pairs"
    assert run["merge_split"]["n_groups_used"] == 1
    # dataset_claim_evidence も spec だけで再現する: 最後の結果を使わず、いつも既定の閾値の組（文書の約束。§7・§9）
    claim_after = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert claim_after["default_selection"] == claim_before["default_selection"]
    assert claim_after["rule_hash"] == claim_before["rule_hash"]
    assert claim_after["merge_groups_hash"] == claim_before["merge_groups_hash"]


def test_last_alignment_qc_is_dropped_when_exclusions_or_data_change(tmp_path):
    handle = _split_handle(tmp_path)
    ds = handles.registry.get(handle).dataset
    assert alignment_qc.remembered(ds) is None
    R.dataset_alignment_qc(dataset_handle=handle)
    assert alignment_qc.remembered(ds) is not None
    out = json.loads(sample_tools.dataset_exclude(samples=["QC2"], reason="テスト",
                                                  dataset_handle=handle))
    assert "error" not in out, out
    assert alignment_qc.remembered(ds) is None             # 除外が変われば使わない（捨てる）
    assert ds.last_alignment_qc is None
    R.dataset_alignment_qc(dataset_handle=handle)
    assert alignment_qc.remembered(ds) is not None         # 今の除外で作り直せば使う
    ds.feature_matrix = ds.feature_matrix * 2              # データが変われば使わない
    assert alignment_qc.remembered(ds) is None


def test_last_alignment_qc_lives_on_its_own_handle(tmp_path):
    first = _split_handle(tmp_path / "a")
    second = _split_handle(tmp_path / "b")
    R.dataset_alignment_qc(dataset_handle=first)
    assert alignment_qc.remembered(handles.registry.get(first).dataset) is not None
    assert alignment_qc.remembered(handles.registry.get(second).dataset) is None


# ---------- dataset_novel_candidates: 行の上限と notes（C7） ----------

def _novel_rows(n):
    return [{"feature_id": str(i), "name": f"M{i}", "inchikey": None, "adduct": "[M+H]+",
             "mz": 100.0 + i, "rt": 1.0, "ppm": 1.0, "group_ref": "HC", "group_test": "RA",
             "log2fc_test_over_ref": 1.0, "p": 0.001 * (i + 1), "q": 0.002 * (i + 1),
             "direction": "test_higher", "merged_from": None} for i in range(n)]


def test_novel_candidates_over_50_rows_are_truncated_not_refused(tmp_path, monkeypatch):
    handle = _split_handle(tmp_path, SPLIT_CHANGED)
    diff = _compare_on(handle)
    monkeypatch.setattr(mass_accuracy, "novel_candidates", lambda *a, **k: _novel_rows(60))
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    out = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    assert "error" not in out and out["n_rows"] == 60 and out["truncated"] is True
    rows = _tsv_rows(out["rows"])
    assert len(rows) == R.MAX_INLINE_ROWS and rows[0]["feature_id"] == "0" and rows[-1]["feature_id"] == "49"
    assert any("output_path" in n for n in out["notes"])
    written = tmp_path / "novel.json"
    full = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle, output_path=str(written)))
    assert full["n_rows"] == 60 and "truncated" not in full and "notes" not in full
    assert len(json.loads(written.read_text(encoding="utf-8"))["rows"]) == 60
    monkeypatch.setattr(mass_accuracy, "novel_candidates", lambda *a, **k: _novel_rows(2))
    small = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    assert small["truncated"] is False and "notes" not in small


def test_novel_candidates_ineffective_settings_are_notes(tmp_path):
    handle = _split_handle(tmp_path, SPLIT_CHANGED)
    diff = _compare_on(handle)
    base = {"statistics": {"result_id": diff["result_id"], "merge_split": False}, "q_max": 1.0}
    plain = json.loads(R.dataset_novel_candidates(base, dataset_handle=handle))
    ignored = json.loads(R.dataset_novel_candidates(
        {**base, "statistics": {**base["statistics"], "alignment_qc": "not-read.json"}},
        dataset_handle=handle))
    assert "error" not in ignored, ignored
    assert ignored["rule_hash"] == plain["rule_hash"]
    assert any("statistics.alignment_qc" in n for n in ignored["notes"])
    assert "notes" not in plain
    # 埋め込んだ analysis_spec の効かない設定は、その欄の名前を付けて notes に
    analysis = {"groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
                "comparison": {"reference": "HC", "test": "RA"},
                "analyses": [{"analysis_id": "main", "test": {"kind": "welch"},
                              "matrix": {"steps": [], "impute": {"method": "half_min"}}}]}
    from_spec = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": analysis, "analysis_id": "main"}, "q_max": 1.0},
        dataset_handle=handle))
    assert from_spec["status"] == "success", from_spec
    assert any(n.startswith("statistics.analysis_spec: analyses[0].matrix") for n in from_spec["notes"])
    # analysis_spec の隣の statistics.contrast_id は効かない（比較は spec の 1 つ）ので、断らず捨てて notes に
    with_contrast = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": analysis, "analysis_id": "main", "contrast_id": "RA_vs_HC"},
         "q_max": 1.0}, dataset_handle=handle))
    assert with_contrast["status"] == "success", with_contrast
    assert with_contrast["rule_hash"] == from_spec["rule_hash"]
    assert any(n.startswith("statistics.contrast_id") for n in with_contrast["notes"])
    # result_id の経路の merge_split / alignment_qc も analysis_spec の隣では効かない（まとめは spec の
    # features.merge_split が決める）ので、断らず捨てて欄ごとに notes に（最終検証 WK）
    for extra in ({"merge_split": True}, {"merge_split": False}, {"alignment_qc": "not-read.json"},
                  {"merge_split": False, "alignment_qc": "not-read.json", "contrast_id": "RA_vs_HC"}):
        out = json.loads(R.dataset_novel_candidates(
            {"statistics": {"analysis_spec": analysis, "analysis_id": "main", **extra}, "q_max": 1.0},
            dataset_handle=handle))
        assert out["status"] == "success", (extra, out)
        assert out["rule_hash"] == from_spec["rule_hash"]
        assert out["n_rows"] == from_spec["n_rows"]
        assert out.get("merge_split") == from_spec.get("merge_split")
        for key in extra:
            assert any(n.startswith(f"statistics.{key} は statistics.analysis_spec では効かず")
                       for n in out["notes"]), (key, out["notes"])
    # 統計の出所を 2 つ書いた形の誤りは従来どおり断る（どちらの統計か決まらない）
    both = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": analysis, "analysis_id": "main", "result_id": diff["result_id"]},
         "q_max": 1.0}, dataset_handle=handle))
    assert both["error"]["code"] == "NOVEL_SPEC_INVALID"


def test_claim_evidence_description_names_the_msms_fields():
    # ツールの説明は LLM が読む唯一の入口: msms.mz_tol_ppm と、候補の msms が返す窓・順位・MS-DIAL との一致を書く
    import asyncio

    from metabolomix.core.mcp_core import mcp

    description = {t.name: t for t in asyncio.run(mcp.list_tools())}["dataset_claim_evidence"].description
    for word in ("mz_tol_ppm", "`window`", "target_rank", "top1_is_target", "msdial_name", "msdial_rank",
                 "top_differs_from_msdial"):
        assert word in description, word


# ---------- dataset_analysis_run（v1）: notes・外部の表・correlate・寄与・群の中で一定（C7・D4・B4） ----------

def test_analysis_run_v1_returns_input_read_correlate_contributions_and_notes(tmp_path):
    values = {"F1": [10, 12, 11, 20, 22, 21, 15, 16], "F2": [5, 6, 5, 5, 6, 5, 5, 5],
              "F3": [3, 3, 3, 1, 2, 5, 0, 0]}                  # HC の 3 試料が同じ値で埋めてある
    paper = write_mwtab(tmp_path / "paper.txt", SAMPLES, list(values.items()))
    external = write_mwtab(tmp_path / "external.txt", SAMPLES,
                           [("Glucose", ["1", "3", "2", "5", "8", "6", "4", "ND"])])
    spec = {"input": {"source": "mwtab", "path": str(paper)},
            "groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"},
                          "matrix": {"steps": [], "normalize": {"method": "pqn"}}}],
            "correlate": {"external": {"source": "mwtab", "path": str(external)},
                          "pairs": [{"feature_id": "F1", "external_id": "Glucose"}]}}
    out = json.loads(R.dataset_analysis_run(spec))           # 外部の入力はデータセットが要らない
    assert out["status"] == "success", out
    assert out["input_read"]["source"] == "mwtab" and out["input_read"]["n_unreadable_cells"] == 0
    assert any(n.startswith("analyses[0].matrix:") for n in out["notes"])
    item = out["analyses"][0]
    assert item["n_constant_in_group"] == 1                   # F3 の HC がすべて 3
    rows = {r["feature_id"]: r for r in _tsv_rows(item["rows"])}
    assert rows["F3"]["constant_in_group"] == "HC" and rows["F1"]["constant_in_group"] == ""
    corr = out["correlate"]
    assert corr["stage"] == "input" and corr["n_samples_matched"] == 8 and corr["n_rows"] == 1
    assert corr["external"]["n_unreadable_cells"] == 1 and corr["truncated"] is False
    [row] = _tsv_rows(corr["rows"])
    assert row["feature_id"] == "F1" and row["external_id"] == "Glucose" and row["n"] == "7"
    # 集約（和）は寄与を TSV で返し、output_path の JSON には全量
    agg_spec = {**spec, "aggregate": {"kind": "sum", "label": "total", "feature_ids": ["F1", "F2"]}}
    agg_spec.pop("correlate")
    target = tmp_path / "agg.json"
    agg = json.loads(R.dataset_analysis_run(agg_spec))
    contributions = _tsv_rows(agg["analyses"][0]["aggregate"]["contributions"])
    assert [c["feature_id"] for c in contributions] == ["F1", "F2"]
    assert agg["analyses"][0]["aggregate"]["n_samples_for_shares"] == 6
    R.dataset_analysis_run(agg_spec, output_path=str(target))
    saved = json.loads(target.read_text(encoding="utf-8"))
    assert saved["notes"] == agg["notes"]
    assert [c["feature_id"] for c in saved["analyses"][0]["aggregate"]["contributions"]] == ["F1", "F2"]


def test_analysis_run_v1_aggregate_summary_is_capped_by_share(tmp_path):
    # 和のメンバーが 50 を超えれば、戻り値の寄与は share_median の大きい順の先頭 50 行、メンバーの一覧は
    # 先頭 50 件（件数 n_*・印 *_truncated・notes）。output_path の JSON は全メンバーのまま。
    import numpy as np

    n = 120
    values = {f"F{k:03d}": [(k + 1) * (s + 1) + 50 * (k % 3) for s in range(len(SAMPLES))]
              for k in range(n)}
    paper = write_mwtab(tmp_path / "paper.txt", SAMPLES, list(values.items()))
    spec = {"input": {"source": "mwtab", "path": str(paper)},
            "groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"}}],
            "aggregate": {"kind": "sum", "label": "total", "feature_ids": list(values)}}
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["status"] == "success", out
    agg = out["analyses"][0]["aggregate"]
    assert agg["n_members"] == agg["n_summed"] == agg["n_contributions"] == n
    assert agg["n_dropped"] == 0 and agg["dropped"] == [] and agg["dropped_truncated"] is False
    assert agg["members"] == agg["summed"] == list(values)[:R.MAX_INLINE_ROWS]
    assert agg["members_truncated"] is agg["summed_truncated"] is agg["contributions_truncated"] is True
    assert agg["contributions_order"] == "share_median_desc"
    [note] = [n for n in out["notes"] if n.startswith("集約の要約")]
    assert "['main']" in note and "寄与 contributions を share_median の大きい順の先頭 50 行" in note
    assert "members・summed を先頭 50 件" in note and "dropped" not in note
    # 独立の計算: 比べる 2 群（HC・RA の 6 試料）の試料ごとの割合（メンバー ÷ 和）の中央値、大きい順
    matrix = np.array([values[f] for f in values], dtype=float)[:, :6]
    shares = np.median(matrix / matrix.sum(axis=0), axis=1)
    expected = sorted(range(n), key=lambda k: (-shares[k], k))[:R.MAX_INLINE_ROWS]
    shown = _tsv_rows(agg["contributions"])
    assert [r["feature_id"] for r in shown] == [f"F{k:03d}" for k in expected]
    assert [float(r["share_median"]) for r in shown] == pytest.approx(
        [round(float(shares[k]), 4) for k in expected], abs=1e-12)
    # output_path の JSON は全メンバー（summed の順）。戻り値の要約は同じく切る
    target = tmp_path / "agg.json"
    written = json.loads(R.dataset_analysis_run(spec, output_path=str(target)))
    assert written["analyses"][0]["aggregate"]["contributions_truncated"] is True
    saved = json.loads(target.read_text(encoding="utf-8"))["analyses"][0]["aggregate"]
    assert [c["feature_id"] for c in saved["contributions"]] == list(values)
    assert saved["members"] == saved["summed"] == list(values)
    assert written["rule_hash"] == out["rule_hash"]


def test_analysis_run_v1_aggregate_notes_name_only_what_was_cut(tmp_path):
    # 和に入ったメンバーが 50 以下なら、members だけが 50 を超えても寄与は切らず summed の順のまま。
    # notes は切ったものだけを名指し、寄与を並べ替えたとは書かない（レビュー指摘: 切っていない寄与の
    # 先頭行を「和を支配するメンバー」と読ませない）。
    import numpy as np

    n, n_zero = 60, 20
    values = {f"F{k:03d}": ([0] * len(SAMPLES) if k < n_zero
                            else [(k + 1) * (s + 1) + 50 * (k % 3) for s in range(len(SAMPLES))])
              for k in range(n)}
    ids = list(values)
    paper = write_mwtab(tmp_path / "paper.txt", SAMPLES, list(values.items()))
    spec = {"input": {"source": "mwtab", "path": str(paper)},
            "groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": [{"analysis_id": "raw", "test": {"kind": "student"}},
                         {"analysis_id": "filtered", "test": {"kind": "student"},
                          "matrix": {"steps": ["detection_filter"],
                                     "detection_filter": {"basis": "positive",
                                                          "min_rate": {"sample": 0.5}}}}],
            "aggregate": {"kind": "sum", "label": "total", "feature_ids": ids}}
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["status"] == "success", out
    by_id = {a["analysis_id"]: a["aggregate"] for a in out["analyses"]}

    filtered = by_id["filtered"]
    assert filtered["n_members"] == n and filtered["members_truncated"] is True
    assert filtered["members"] == ids[:R.MAX_INLINE_ROWS]
    assert filtered["n_summed"] == filtered["n_contributions"] == n - n_zero
    assert filtered["summed"] == ids[n_zero:] and filtered["summed_truncated"] is False
    assert filtered["n_dropped"] == n_zero and filtered["dropped"] == ids[:n_zero]
    assert filtered["dropped_truncated"] is False
    assert filtered["contributions_truncated"] is False and "contributions_order" not in filtered
    # 寄与は summed の順の全行。値は独立の計算（HC・RA の 6 試料の割合 メンバー ÷ 和 の中央値）と合う
    matrix = np.array([values[f] for f in ids[n_zero:]], dtype=float)[:, :6]
    shares = np.median(matrix / matrix.sum(axis=0), axis=1)
    shown = _tsv_rows(filtered["contributions"])
    assert [r["feature_id"] for r in shown] == ids[n_zero:]
    assert [float(r["share_median"]) for r in shown] == pytest.approx(
        [round(float(x), 4) for x in shares], abs=1e-12)
    assert list(shares) != sorted(shares, reverse=True)  # summed の順は share の大きい順ではない

    raw = by_id["raw"]
    assert raw["n_summed"] == raw["n_contributions"] == n
    assert raw["contributions_truncated"] is raw["summed_truncated"] is raw["members_truncated"] is True
    assert raw["contributions_order"] == "share_median_desc"

    aggregate_notes = [note for note in out["notes"] if note.startswith("集約の要約")]
    assert len(aggregate_notes) == 2
    [raw_note] = [note for note in aggregate_notes if "['raw']" in note]
    [filtered_note] = [note for note in aggregate_notes if "['filtered']" in note]
    assert "寄与 contributions を share_median の大きい順の先頭 50 行" in raw_note
    assert "members・summed を先頭 50 件" in raw_note
    assert "share_median" not in filtered_note and "contributions" not in filtered_note
    assert "、members を先頭 50 件だけ載せた" in filtered_note
    assert "output_path を指定して書く JSON" in filtered_note
    # output_path を渡せば全量はその JSON（notes の案内も変わる）
    written = json.loads(R.dataset_analysis_run(spec, output_path=str(tmp_path / "agg.json")))
    assert any(note.startswith("集約の要約（解析 ['filtered']）") and "output_path の JSON" in note
               for note in written["notes"])


def test_aggregate_summary_caps_dropped_and_merged_from_and_keeps_small_ones():
    info = {"members": [f"M{i}" for i in range(60)], "summed": ["M0", "M1"],
            "dropped": [f"M{i}" for i in range(2, 60)],
            "merged_from": {f"M{i}": [f"M{i}", f"N{i}"] for i in range(55)},
            "contributions": [{"feature_id": "M0", "share_median": 0.25},
                              {"feature_id": "M1", "share_median": 0.75}],
            "n_samples_for_shares": 4}
    out, cut = R._aggregate_summary(info)
    # 切ったものだけを名指す（寄与と summed は 2 件で切らない）
    assert cut == ["members", "dropped", "merged_from"]
    assert out["n_dropped"] == 58 and len(out["dropped"]) == R.MAX_INLINE_ROWS and out["dropped_truncated"]
    assert out["n_merged_from"] == 55 and list(out["merged_from"]) == [f"M{i}" for i in range(50)]
    assert out["merged_from_truncated"] is True
    # 寄与が 50 以下なら summed の順のまま全行（並べ替えない）
    assert [r["feature_id"] for r in _tsv_rows(out["contributions"])] == ["M0", "M1"]
    assert out["contributions_truncated"] is False and "contributions_order" not in out
    assert out["n_samples_for_shares"] == 4
    small, small_cut = R._aggregate_summary({"numerator": "A", "denominator": "B", "dropped": []})
    assert small_cut == []
    assert small == {"numerator": "A", "denominator": "B", "dropped": [], "n_dropped": 0,
                     "dropped_truncated": False}


def test_analysis_run_v1_correlate_rows_keep_small_p(tmp_path):
    # 行の p は有効 4 桁（小数 4 桁に丸めると 5e-5 未満の p が 0.0 に化ける）、ρ は小数 4 桁。
    # 8 組の順位が一致: 正確な並べ替えの p は 2/8! = 4.96e-5（t 近似の 0 ではない）
    ours, theirs = list(range(1, 9)), [10, 20, 30, 40, 50, 60, 70, 80]
    paper = write_mwtab(tmp_path / "paper.txt", SAMPLES, [("F1", ours)])
    external = write_mwtab(tmp_path / "external.txt", SAMPLES, [("Glucose", theirs)])
    spec = {"input": {"source": "mwtab", "path": str(paper)},
            "groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"}}],
            "correlate": {"external": {"source": "mwtab", "path": str(external)},
                          "pairs": [{"feature_id": "F1", "external_id": "Glucose"}]}}
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["status"] == "success", out
    [row] = _tsv_rows(out["correlate"]["rows"])
    expected = stats.spearmanr(ours, theirs)
    exact = 2 / 40320
    assert exact < 5e-5
    assert float(row["p"]) == float(f"{exact:.4g}") and float(row["p"]) > 0
    assert row["p_method"] == "exact_permutation"
    assert row["rho"] == str(round(float(expected.statistic), 4))


def test_analysis_run_aggregate_with_feature_ids_is_spec_invalid(tmp_path):
    paper = write_mwtab(tmp_path / "paper.txt", SAMPLES, [("F1", list(range(8))), ("F2", list(range(8)))])
    spec = {"input": {"source": "mwtab", "path": str(paper)},
            "groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "features": {"feature_ids": ["F1"]},
            "aggregate": {"kind": "sum", "label": "total", "feature_ids": ["F1", "F2"]},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"}}]}
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["error"]["code"] == "ANALYSIS_SPEC_INVALID"
    assert out["error"]["details"]["field"] == "features"


def test_split_features_fixture_is_unchanged():
    # このファイルが借りる fixture（SPLIT_FEATURES の 30/31 の割れ）の前提
    assert [f["id"] for f in SPLIT_FEATURES] == ["30", "31", "32"]
