import json
import shutil
import subprocess

import pytest

from lipidmix.arf2.reader import load_catalog
from lipidmix.curation import evidence, flags, judge, review, viewer
from lipidmix.library import store as library_store
from tests.curation_fixtures import write_alignment_set


@pytest.fixture()
def built(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg")
    s = library_store.open_store(paths["msp"])
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    result = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                               th=judge.resolve_thresholds(None), file_ids=None, max_traces=12,
                               selection={"kind": "annotated"})
    yield paths, result
    s.close()


def test_verdicts_follow_the_evidence(built):
    _, result = built
    by_id = {s["spot_id"]: s for s in result["spots"]}
    assert by_id[0]["verdict"] == "ok"
    # spot 1: low_score（弱）+ drt_out 1.5 分（弱）+ ギザギザの EIC（弱）→ suspect
    assert by_id[1]["verdict"] == "suspect"
    assert {"low_score", "drt_out", "eic_poor"} <= set(by_id[1]["reasons"])
    assert result["counts"] == {"ok": 1, "suspect": 1, "likely_wrong": 0}


def test_summary_tsv_lists_only_flagged_spots(built):
    _, result = built
    lines = review.summary_tsv(result).splitlines()
    assert lines[0].split("\t")[:3] == ["spot_id", "name", "ontology"]
    assert [line.split("\t")[0] for line in lines[1:]] == ["1"]


def test_save_and_load_round_trip(built):
    paths, result = built
    saved = review.save_review(result)
    assert saved["json"].parent == flags.curation_dir(paths["arf2"])
    loaded = review.load_review(paths["arf2"], result["review_id"])
    assert loaded["review_id"] == result["review_id"]
    html = saved["html"].read_text(encoding="utf-8")
    assert "</script>" in html
    assert result["review_id"] in html


def test_embedded_json_cannot_break_out_of_the_script_tag(built):
    _, result = built
    result["spots"][0]["name"] = "</script><b>x</b>"
    html = viewer.render_html(result)
    assert "</script><b>" not in html


def test_existing_flags_are_attached_to_spots(built):
    paths, result = built
    store = flags.FlagStore(flags.curation_dir(paths["arf2"]))
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=result["alignment"],
                 review_id="old", source="user")
    s = library_store.open_store(paths["msp"])
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        again = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                  th=judge.resolve_thresholds(None), file_ids=None,
                                  max_traces=12, selection={"kind": "annotated"})
    finally:
        s.close()
    assert {sp["spot_id"]: sp["flag"] for sp in again["spots"]} == {0: None, 1: "wrong"}


def test_low_reference_resolution_warns_about_the_library(built, tmp_path):
    paths, _ = built
    other = tmp_path / "other.msp"
    other.write_text("NAME: X\nPRECURSORMZ: 100.0\nIONMODE: Positive\nNum Peaks: 0\n", encoding="utf-8")
    s = library_store.open_store(other)
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        result = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                   th=judge.resolve_thresholds(None), file_ids=None,
                                   max_traces=12, selection={"kind": "annotated"})
    finally:
        s.close()
    assert any("別のライブラリ" in w for w in result["warnings"])


def test_page_splits_spots(built):
    _, result = built
    first = review.page(result, 0, page_size=1)
    assert first["n_pages"] == 2 and len(first["spots"]) == 1


def test_app_template_has_no_embedded_data():
    html = viewer.render_html(None)
    assert "const EMBEDDED = null" in html


def test_template_avoids_horizontal_overflow_at_phone_width():
    # 実ファイル名にはスペースがなく、header h1 が縮められないと sticky header ごと
    # 横スクロールが出る（binding requirement: no horizontal page scroll at phone
    # width）。#trends / #grid の minmax も 320px 幅の携帯（16px ガター）では
    # 300/320px 未満に縮まないと同様に溢れる。
    html = viewer.render_html(None)
    assert "overflow-wrap:anywhere" in html
    assert "min(300px, 100%)" in html


def test_render_builds_cards_in_chunks_with_a_cancellable_generation_counter():
    # 1468 件のような大きなレビューで render() が 1 回の同期ループで全カードを
    # 組み立てると、ページが 10-15 秒応答しなくなる（実データ check）。チャンク分割し、
    # フィルタ変更などで新しい render() が割り込んだら古い方を打ち切る世代カウンタが要る。
    html = viewer.render_html(None)
    assert "renderGeneration" in html
    assert "generation !== renderGeneration" in html
    assert "requestAnimationFrame" in html
    assert "setTimeout" in html
    assert "DocumentFragment" in html


def test_draw_trend_clips_the_regression_line_to_the_plot_rectangle():
    # slope*x+intercept は RT レンジの外まで伸びうり、パネルからはみ出して見える
    # （実データ check）。プロット矩形で ctx.clip() する。
    html = viewer.render_html(None)
    start = html.index("function drawTrend")
    end = html.index("function spotCard")
    trend_source = html[start:end]
    assert "ctx.save()" in trend_source
    assert "ctx.clip()" in trend_source
    assert "ctx.restore()" in trend_source


def test_extracted_script_is_valid_javascript(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が見つからないので構文チェックを省略")
    html = viewer.render_html(None)
    start = html.index("<script>") + len("<script>")
    end = html.index("</script>")
    script_path = tmp_path / "viewer.js"
    script_path.write_text(html[start:end], encoding="utf-8")
    result = subprocess.run([node, "--check", str(script_path)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# ---------- I2: TSV の行数上限 ----------

def test_summary_tsv_keeps_the_first_max_rows_by_severity(built):
    _, result = built
    result["spots"][0]["flag"] = "suspect"      # ok だがフラグ済み → 表に載る(2 行になる)
    assert review.n_summary_rows(result) == 2
    lines = review.summary_tsv(result, max_rows=1).splitlines()
    assert len(lines) == 2                      # 見出し + 1 行
    assert lines[1].split("\t")[0] == "1"       # suspect が先(重い順)


# ---------- M1: review_id の形式 ----------

@pytest.mark.parametrize("bad", ["../x", "cr-1", "cr-20260101-000000-zzzz", "cr-20260101-000000-abcd/..", ""])
def test_load_review_rejects_malformed_review_ids(built, bad):
    paths, _ = built
    with pytest.raises(ValueError):
        review.load_review(paths["arf2"], bad)


def test_generated_review_ids_are_valid(built):
    _, result = built
    assert review.is_valid_review_id(result["review_id"])


# ---------- M5: 置換の順序と `<` のエスケープ ----------

def test_spot_names_containing_the_prefix_placeholder_survive(built):
    _, result = built
    result["spots"][0]["name"] = "__SUBMISSION_PREFIX__ <!-- x"
    html = viewer.render_html(result)
    start = html.index("const EMBEDDED = ") + len("const EMBEDDED = ")
    end = html.index(";\nconst PREFIX")
    embedded = html[start:end]
    assert "<" not in embedded
    assert json.loads(embedded)["spots"][0]["name"] == "__SUBMISSION_PREFIX__ <!-- x"


# ---------- M7: クラス別傾向の要約 ----------

def test_trend_summary_includes_outlier_counts():
    fake = {"trend": {"classes": {"PC": {"n": 9, "r2": 0.95, "n_outliers": 1, "coef": [1, 2]}}}}
    assert review.trend_summary(fake) == {"PC": {"n": 9, "r2": 0.95, "n_outliers": 1}}


# ---------- I4 / I5 / M4: ビューアのテンプレート ----------

def _script():
    html = viewer.render_html(None)
    return html[html.index("<script>"):html.index("</script>")]


def test_submission_includes_the_arf2_path():
    script = _script()
    body = script[script.index("function submission"):]
    body = body[:body.index("\n}")]
    assert "REVIEW.arf2_path" in body


def test_trend_canvas_click_jumps_to_the_spot_card():
    script = _script()
    assert "jumpToSpot" in script
    trends = script[script.index("function renderTrends"):]
    trends = trends[:trends.index("\n}")]
    assert 'addEventListener("click"' in trends
    assert "HIT_RADIUS" in script
    jump = script[script.index("function jumpToSpot"):]
    jump = jump[:jump.index("\n}")]
    assert "scrollIntoView" in jump and "highlight" in jump
    assert "resetFilters" in jump


def test_card_canvases_open_an_enlarged_dialog():
    html = viewer.render_html(None)
    assert '<dialog id="zoom"' in html
    script = _script()
    card = script[script.index("function spotCard"):]
    card = card[:card.index("\n}")]
    assert "openZoom" in card
    zoom = script[script.index("function openZoom"):]
    zoom = zoom[:zoom.index("\n}")]
    assert "drawEic" in zoom and "drawMirror" in zoom
    assert "Escape" in script
    assert 'id="zoom-close"' in html


def test_note_input_reflects_unsent_edits():
    script = _script()
    card = script[script.index("function spotCard"):]
    card = card[:card.index("\n}")]
    assert "edits.get(spot.spot_id).note" in card
    assert "spot.mz?.toFixed" not in card


def test_app_loader_shows_error_payloads_and_send_updates_spots():
    script = _script()
    assert "showError" in script
    send = script[script.index('getElementById("send").addEventListener'):]
    send = send[:send.index("\n});")]
    assert "spot.flag = " in send and "spot.flag_note = " in send
    assert send.index("spot.flag = ") < send.index("edits.clear()")


def _review(paths):
    s = library_store.open_store(paths["msp"])
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        return review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                 th=judge.resolve_thresholds(None), file_ids=None,
                                 max_traces=12, selection={"kind": "annotated"})
    finally:
        s.close()


def test_class_rule_rejection_is_likely_wrong_when_lipid_rules_ran(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    # spot 0 は既定のまま（class=True）＝このデータは脂質規則で採点されている。
    paths = write_alignment_set(tmp_path / "neg", match_overrides={1: {19: False, 20: False, 22: False}})
    result = _review(paths)
    by_id = {s["spot_id"]: s for s in result["spots"]}
    assert by_id[1]["verdict"] == "likely_wrong"
    assert by_id[1]["reasons"][0] == "class_rule_rejected"
    assert not any("脂質規則" in w for w in result["warnings"])


def test_rule_flags_are_ignored_when_no_spot_shows_lipid_rules(tmp_path, monkeypatch):
    # メタボロミクス採点器の出力を模す: どのスポットも規則フラグが全部 False。
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    off = {19: False, 20: False, 22: False}
    paths = write_alignment_set(tmp_path / "neg", match_overrides={0: off, 1: off})
    result = _review(paths)
    by_id = {s["spot_id"]: s for s in result["spots"]}
    assert by_id[1]["verdict"] == "suspect"
    assert "class_rule_rejected" not in by_id[1]["reasons"]
    assert any("脂質規則" in w for w in result["warnings"])


# --- ビューア: 自動判別 likely_wrong の赤破線枠と、クラス選択肢の下位項目（2026-09-29 ユーザー決定）---

_CLASS_FILTER_BEGIN = "// --- class filter (pure) ---"
_CLASS_FILTER_END = "// --- end class filter ---"


def _run_class_filter(tmp_path, expression: str):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が見つからないのでビューアの純関数を実行できない")
    html = viewer.render_html(None)
    source = html[html.index(_CLASS_FILTER_BEGIN):html.index(_CLASS_FILTER_END)]
    spots = [{"spot_id": 0, "ontology": "PC", "verdict": "likely_wrong"},
             {"spot_id": 1, "ontology": "PC", "verdict": "ok"},
             {"spot_id": 2, "ontology": "PG", "verdict": "suspect"},
             {"spot_id": 3, "ontology": "Cer_NS", "verdict": "likely_wrong"},
             {"spot_id": 4, "ontology": "PC", "verdict": "likely_wrong"},
             {"spot_id": 5, "ontology": "", "verdict": "likely_wrong"}]
    script = tmp_path / "class_filter.js"
    script.write_text(source + f"\nconst SPOTS = {json.dumps(spots)};\n"
                      f"console.log(JSON.stringify({expression}));\n", encoding="utf-8")
    result = subprocess.run([node, str(script)], capture_output=True, text=True, encoding="utf-8")
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_class_options_put_likely_wrong_subitems_under_each_class(tmp_path):
    options = _run_class_filter(tmp_path, "classFilterOptions(SPOTS)")
    labels = [o["label"] for o in options]
    assert labels == ["すべて",
                      "⚠ 自動判別: 間違い（全クラス 4）",
                      "Cer_NS (1)", "└ Cer_NS › 自動判別: 間違い (1)",
                      "PC (3)", "└ PC › 自動判別: 間違い (2)",
                      "PG (1)"]
    assert len({o["value"] for o in options}) == len(options)


def test_class_options_omit_the_overall_item_when_nothing_is_likely_wrong(tmp_path):
    options = _run_class_filter(
        tmp_path, "classFilterOptions(SPOTS.map(s => ({...s, verdict: 'ok'})))")
    assert [o["label"] for o in options] == ["すべて", "Cer_NS (1)", "PC (3)", "PG (1)"]


def test_class_filter_values_select_the_right_spots(tmp_path):
    picked = _run_class_filter(tmp_path, """Object.fromEntries(classFilterOptions(SPOTS).map(o =>
        [o.label, SPOTS.filter(s => matchesClassFilter(s, o.value)).map(s => s.spot_id)]))""")
    assert picked["すべて"] == [0, 1, 2, 3, 4, 5]
    assert picked["⚠ 自動判別: 間違い（全クラス 4）"] == [0, 3, 4, 5]
    assert picked["PC (3)"] == [0, 1, 4]
    assert picked["└ PC › 自動判別: 間違い (2)"] == [0, 4]
    assert picked["PG (1)"] == [2]


def test_likely_wrong_cards_get_a_dashed_red_border_that_turns_solid_when_flagged_wrong():
    html = viewer.render_html(None)
    assert ".card.verdict-likely_wrong { border:2px dashed var(--wrong); }" in html
    assert ".card.verdict-likely_wrong.flag-wrong { border-style:solid; }" in html
    start = html.index("function spotCard")
    card_source = html[start:html.index("function setFlag")]
    assert '"verdict-" + spot.verdict' in card_source


def test_summary_tsv_has_the_mz_difference_column(built):
    _, result = built
    header = review.summary_tsv(result).splitlines()[0].split("\t")
    assert header[header.index("ppm") + 1] == "dmz_mda"


def test_review_attaches_auto_notes_and_cleared_state(built):
    paths, result = built
    by_id = {s["spot_id"]: s for s in result["spots"]}
    assert by_id[0]["auto_note"] is None
    assert by_id[1]["auto_note"].startswith("自動: ")
    assert by_id[1]["flag_cleared"] is False
    store = flags.FlagStore(flags.curation_dir(paths["arf2"]))
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=result["alignment"],
                 review_id="old", source="user")
    store.append([{"spot_id": 1, "flag": "clear"}], alignment=result["alignment"],
                 review_id="old2", source="user")
    again = _review(paths)
    assert {s["spot_id"]: s["flag_cleared"] for s in again["spots"]} == {0: False, 1: True}


# --- ビューア: likely_wrong の「間違い」プリセット・根拠メモ・ミラーの m/z ラベルと縦軸目盛り
# （ユーザー決定 2026-09-29。ラベルは library_plot_mirror の _draw_labels（auto 方式）と同じ）---

def _run_block(tmp_path, block: str, prelude: str, expression: str):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が見つからないのでビューアの純関数を実行できない")
    html = viewer.render_html(None)
    begin, end = f"// --- {block} (pure) ---", f"// --- end {block} ---"
    source = html[html.index(begin):html.index(end)]
    script = tmp_path / f"{block.replace(' ', '_')}.js"
    script.write_text(source + "\n" + prelude + f"\nconsole.log(JSON.stringify({expression}));\n",
                      encoding="utf-8")
    result = subprocess.run([node, str(script)], capture_output=True, text=True, encoding="utf-8")
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


_PRESET_SPOTS = json.dumps([
    {"spot_id": 0, "verdict": "likely_wrong", "flag": None, "flag_cleared": False, "auto_note": "自動: A"},
    {"spot_id": 1, "verdict": "likely_wrong", "flag": "suspect", "flag_cleared": False, "auto_note": "自動: B"},
    {"spot_id": 2, "verdict": "likely_wrong", "flag": None, "flag_cleared": True, "auto_note": "自動: C"},
    {"spot_id": 3, "verdict": "suspect", "flag": None, "flag_cleared": False, "auto_note": "自動: D"},
    {"spot_id": 4, "verdict": "ok", "flag": None, "flag_cleared": False, "auto_note": None},
    {"spot_id": 5, "verdict": "suspect", "flag": "wrong", "flag_note": "人のメモ", "flag_cleared": False,
     "auto_note": "自動: E"}])


def test_likely_wrong_is_preset_to_wrong_unless_already_recorded_or_cleared(tmp_path):
    presets = _run_block(tmp_path, "presets", f"const SPOTS = {_PRESET_SPOTS};",
                         "SPOTS.filter(isPresetTarget).map(s => [s.spot_id, presetEdit(s)])")
    assert presets == [[0, {"flag": "wrong", "note": "自動: A", "preset": True}]]


def test_note_defaults_to_the_recorded_note_then_the_auto_note(tmp_path):
    notes = _run_block(tmp_path, "presets", f"const SPOTS = {_PRESET_SPOTS};",
                       "SPOTS.map(defaultNote)")
    assert notes == ["自動: A", "自動: B", "自動: C", "自動: D", "", "人のメモ"]


_LABEL_PRELUDE = """
const measure = text => ({w: 40, h: 8});           // 代表箱: 幅 40px・高さ 8px
const toPixel = (mz, y) => ({x: mz, y: 100 - y * 100});   // m/z をそのまま px に
"""


def test_mirror_labels_are_picked_by_height_and_skip_2d_overlaps(tmp_path):
    labels = _run_block(tmp_path, "mirror labels", _LABEL_PRELUDE,
                        "pickMirrorLabels([[100, 0.5], [110, 1.0], [300, 0.2], [120, 0.1]], toPixel, measure)"
                        ".map(l => l.text)")
    # 110 (最大) → 100 は横も縦も近い（Δx 10 < 40、Δy 50 > 8 なので縦は離れている）ので残る。
    # 120 は 110 と Δx 10・Δy 90 で縦に離れて残る。auto は横と縦の両方が近いときだけ飛ばす。
    assert labels == ["110.0000", "100.0000", "300.0000", "120.0000"]


def test_mirror_labels_skip_a_peak_that_overlaps_both_ways(tmp_path):
    labels = _run_block(tmp_path, "mirror labels", _LABEL_PRELUDE,
                        "pickMirrorLabels([[100, 1.0], [105, 0.97], [200, 0.5]], toPixel, measure)"
                        ".map(l => l.text)")
    assert labels == ["100.0000", "200.0000"]


def test_mirror_labels_stop_at_25_per_side(tmp_path):
    labels = _run_block(tmp_path, "mirror labels", _LABEL_PRELUDE,
                        "pickMirrorLabels(Array.from({length: 40}, (_, i) => [i * 100, 1 - i / 100]),"
                        " toPixel, measure).length")
    assert labels == 25


def test_mirror_draws_labels_per_side_and_relative_intensity_ticks():
    html = viewer.render_html(None)
    source = html[html.index("function drawMirror"):html.index("function drawTrend")]
    # 実測・参照で別々に選ぶ（側をまたいで干渉しない）
    assert "for (const [pts, sign] of [[measured, 1], [reference, -1]])" in source
    assert "pickMirrorLabels(pts," in source
    assert "[0, 50, 100]" in source                        # 縦軸の目盛り（相対強度 %、上下とも）


def test_preset_cards_keep_the_dashed_border_until_the_flag_is_recorded():
    html = viewer.render_html(None)
    assert '"flag-preset"' in html
    assert ".card.verdict-likely_wrong.flag-wrong { border-style:solid; }" in html
    assert ".card.flag-preset" not in html or "outline" not in html.split(".card.flag-preset")[1].split("}")[0]


def test_status_keeps_showing_pending_changes_after_render():
    # render() の完了時に countsText() で上書きするので、未送信の変更（プリセット含む）を
    # そこに含めないと「likely_wrong N 件を初期選択」が絞り込みのたびに消える（実データ check）。
    html = viewer.render_html(None)
    source = html[html.index("function countsText"):html.index("function countsText") + 400]
    assert "edits.size" in source and "presetCount" in source


def test_mirror_labels_are_kept_inside_the_plot_area(tmp_path):
    # 端のピークのラベルが縦軸の目盛りに重なった（実データ check）。中心を [xmin+w/2, xmax-w/2] に寄せる。
    labels = _run_block(tmp_path, "mirror labels", _LABEL_PRELUDE,
                        "pickMirrorLabels([[5, 1.0], [298, 0.5]], toPixel, measure, {xmin: 0, xmax: 300})"
                        ".map(l => l.x)")
    assert labels == [20, 280]


def test_assign_and_redundant_are_decisions_not_flags(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg")
    first = _review(paths)
    store = flags.FlagStore(flags.curation_dir(paths["arf2"]))
    store.append([{"spot_id": 0, "flag": "assign", "name": "PC 34:1", "level": "sum", "note": ""},
                  {"spot_id": 1, "flag": "redundant", "of": 0, "relation": "isotope_M+1", "note": ""}],
                 alignment=first["alignment"], review_id="cs-x", source="user")
    again = {s["spot_id"]: s for s in _review(paths)["spots"]}
    assert again[0]["flag"] is None and again[0]["decision"] == {"flag": "assign", "name": "PC 34:1"}
    assert again[1]["flag"] is None and again[1]["decision"] == {"flag": "redundant", "of": 0}
