import json
import shutil
import subprocess

import pytest

from metabolomix.curation import viewer

SPOT = {"spot_id": 1, "target_kind": "flagged", "preset": "R1", "measured": [[255.23, 999.0]],
        "candidates": [{"candidate_id": "L1", "name": "PG 16:0_18:1", "sum_name": "PG 34:1",
                        "mirror": {"reference": [[255.23, 999.0]], "matched_mz": [255.23],
                                   "matched_measured_mz": [255.23]}}],
        "relations": [{"candidate_id": "R1", "informational": False},
                      {"candidate_id": "R2", "informational": True}]}


def _run(tmp_path, expression, prelude=""):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が無いのでビューアの純関数を実行できない")
    html = viewer.render_suggest_html(None)
    begin, end = "// --- suggest (pure) ---", "// --- end suggest ---"
    source = html[html.index(begin):html.index(end)]
    script = tmp_path / "s.js"
    script.write_text(source + "\n" + prelude + f"\nconsole.log(JSON.stringify({expression}));\n", encoding="utf-8")
    result = subprocess.run([node, str(script)], capture_output=True, text=True, encoding="utf-8")
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_initial_choice_is_the_strong_relation_only(tmp_path):
    spot = json.dumps(SPOT)
    assert _run(tmp_path, f"initialChoice({spot})") == {"candidate": "R1", "flag": "redundant", "level": "sum"}
    no_preset = json.dumps({**SPOT, "preset": None})
    assert _run(tmp_path, f"initialChoice({no_preset})") is None


def test_choice_entries_carry_only_ids(tmp_path):
    spot = json.dumps(SPOT)
    entry = _run(tmp_path, f"choiceEntry({spot}, {{candidate: 'L1', flag: 'assign', level: 'species', note: 'x'}})")
    assert entry == {"spot_id": 1, "flag": "assign", "candidate": "L1", "level": "species", "note": "x"}
    restore = _run(tmp_path, f"choiceEntry({spot}, {{flag: 'clear', note: ''}})")
    assert restore == {"spot_id": 1, "flag": "clear", "note": ""}


def test_submission_uses_the_suggestion_id(tmp_path):
    out = _run(tmp_path, "suggestSubmission({suggestion_id: 'cs-1', arf2_path: 'C:/a.arf2', spots: []}, new Map([[1, {candidate: 'R1', flag: 'redundant'}]]))")
    assert out == {"review_id": "cs-1", "arf2_path": "C:/a.arf2",
                   "flags": [{"spot_id": 1, "flag": "redundant", "candidate": "R1", "note": ""}]}


def test_mirror_for_combines_spot_measured_with_candidate_reference(tmp_path):
    spot = json.dumps(SPOT)
    out = _run(tmp_path, f"mirrorFor({spot}, {spot}.candidates[0])")
    assert out == {"mirror": {"measured": [[255.23, 999.0]], "reference": [[255.23, 999.0]],
                              "matched_mz": [255.23]}}


def test_overlay_adds_the_partner_trace_scaled_to_the_representative(tmp_path):
    spot = json.dumps({"eic": {"samples": [
        {"file_id": 0, "representative": True, "detected": True, "points": [[10.0, 0.0], [10.1, 100.0]]},
        {"file_id": 1, "representative": False, "detected": True, "points": [[10.0, 0.0], [10.1, 50.0]]}]}})
    relation = json.dumps({"partner_eic": {"file_id": 0, "left": 9.9, "top": 10.1, "right": 10.3,
                                           "points": [[10.0, 0.0], [10.1, 1000.0]]}})
    out = _run(tmp_path, f"overlayFor({spot}, {relation}).eic.samples")
    assert len(out) == 3
    assert out[2]["partner"] is True and out[2]["representative"] is False
    assert out[2]["points"][1] == [10.1, 100.0]                     # Y の最大 1000 → X の代表の最大 100
    assert _run(tmp_path, f"overlayFor({spot}, {{partner_eic: null}}).eic.samples.length") == 2


def test_eic_draws_the_partner_in_its_own_style():
    html = viewer.render_suggest_html(None)
    source = html[html.index("function drawEic"):html.index("// --- mirror labels (pure) ---")]
    assert "s.partner" in source


def test_embedded_data_is_escaped():
    html = viewer.render_suggest_html({"suggestion_id": "cs-1", "spots": [{"name": "</script>"}]})
    assert "\\u003c/script>" in html
    assert html.count("</script>") == 1          # ページ自身の閉じタグだけ。データ側の </script> は残らない


def test_not_compared_sentinel_is_shown_as_a_dash(tmp_path):
    out = _run(tmp_path, "[scoreText(-1), scoreText(0.8567), scoreText(0), scoreText(null), scoreText(undefined), scoreText(NaN), countText(-1), countText(3), countText(0)]")
    assert out == ["–", "0.86", "0.00", "–", "–", "–", "–", "3", "0"]


def test_candidate_label_tells_records_with_the_same_sum_name_apart(tmp_path):
    same = json.dumps({"name": "PC 34:1", "sum_name": "PC 34:1", "adduct": "[M-H]-"})
    assert _run(tmp_path, f"candidateLabel({same})") == "PC 34:1 [M-H]-"
    species = json.dumps({"name": "PC 16:0_18:1", "sum_name": "PC 34:1", "adduct": "[M-H]-"})
    assert _run(tmp_path, f"candidateLabel({species})") == "PC 34:1（PC 16:0_18:1） [M-H]-"
    bare = json.dumps({"name": "X", "sum_name": None, "adduct": None})
    assert _run(tmp_path, f"candidateLabel({bare})") == "X"


def test_candidate_detail_shows_a_short_inchikey(tmp_path):
    base = {"scores": {"total_score": 0.9, "weighted_dot_product": 0.8, "matched_peaks_count": 3},
            "dmz_mda": 1.5, "soft": ["trend_outlier"], "info": []}
    with_key = json.dumps({**base, "inchikey": "ABCDEFGHIJKLMN-OPQRSTUVWX-N"})
    detail = _run(tmp_path, f"candidateDetail({with_key})")
    assert "ABCDEFGHIJKLMN" in detail and "OPQRST" not in detail and detail.endswith("trend_outlier")
    empty = json.dumps({**base, "inchikey": ""})
    assert "mDa / – / trend_outlier" in _run(tmp_path, f"candidateDetail({empty})")


def test_apply_suggest_recorded_marks_spots(tmp_path):
    out = _run(tmp_path, "(applySuggestRecorded(SPOTS, new Map([[1, {candidate: 'L1', flag: 'assign', level: 'species', note: ''}], [2, {flag: 'clear', note: ''}]])), SPOTS.map(s => s.recorded))",
               prelude="const SPOTS = [{spot_id: 1}, {spot_id: 2}, {spot_id: 3}];")
    assert out == [{"flag": "assign", "candidate": "L1", "level": "species"}, {"flag": "clear"}, None]


def test_recorded_text_names_the_decision(tmp_path):
    out = _run(tmp_path, "[recordedText({flag:'assign', candidate:'L1', level:'sum'}), recordedText({flag:'redundant', candidate:'R1'}), recordedText({flag:'clear'})]")
    assert out[0].startswith("記録済み") and "L1" in out[0]
    assert "R1" in out[1] and "元の注釈" in out[2]


def test_clear_sent_choices_keeps_choices_changed_while_sending(tmp_path):
    prelude = (
        "const A = {candidate: 'L1', flag: 'assign', level: 'sum', note: ''};\n"
        "const sent = new Map([[1, A], [2, A], [3, A], [4, A], [5, A]]);\n"
        "const cur = new Map([[1, {...A}], [2, {...A, note: 'x'}], [3, {...A, level: 'species'}],\n"
        "                     [4, {...A, candidate: 'L2'}], [6, {...A}]]);\n"   # 5 は送信中に外した、6 は新規
        "clearSentChoices(cur, sent);\n")
    out = _run(tmp_path, "[...cur.keys()]", prelude=prelude)
    assert out == [2, 3, 4, 6]


def _script():
    html = viewer.render_suggest_html(None)
    return html[html.index("<script>"):html.index("</script>")]


def test_suggest_viewer_has_submit_and_finish_wired():
    html = viewer.render_suggest_html(None)
    for element in ('id="submit"', 'id="finish"', 'id="result"'):
        assert element in html
    script = _script()
    handler = script[script.index('getElementById("submit").addEventListener'):]
    handler = handler[:handler.index("\n});")]
    assert handler.index("if (!choices.size)") < handler.index("confirm(suggestConfirmText(")
    assert handler.index("new Map(choices)") < handler.index("confirm(suggestConfirmText(")
    assert handler.index("confirm(suggestConfirmText(") < handler.index("submitClient.submit(")
    assert handler.index("disabled = true") < handler.index("submitClient.submit(") < handler.index("disabled = false")


def test_suggest_extracted_script_is_valid_javascript(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が無いので構文チェックを省略")
    path = tmp_path / "suggest.js"
    path.write_text(_script()[len("<script>"):], encoding="utf-8")
    result = subprocess.run([node, "--check", str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
