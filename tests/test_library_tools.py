"""MCP ツール。戻り値の量と missing_state 契約を固定する。"""
import json
import struct
import textwrap

import pytest

from metabolomix.core import mcp_core, session_state
from metabolomix.library import tools

_MSP = textwrap.dedent("""\
    NAME: GABA
    PRECURSORMZ: 104.0706
    IONMODE: Positive
    Num Peaks: 2
    87.0441 999
    69.0335 500
""")


@pytest.fixture(autouse=True)
def fresh_session(tmp_path, monkeypatch):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)   # 実物は pathlib.Path
    monkeypatch.setenv("LIPIDMIX_LIBRARY_CACHE_DIR", str(tmp_path / "cache"))
    (tmp_path / "lib.msp").write_text(_MSP, encoding="utf-8")
    return tmp_path


def test_matching_without_a_loaded_library_returns_missing_state():
    payload = json.loads(tools.library_match_feature(104.07))
    assert payload["error"]["code"] == "missing_state"
    assert "library_load" in payload["error"]["required_tools"]


def test_plotting_without_a_match_returns_missing_state():
    payload = json.loads(tools.library_plot_mirror())
    assert payload["error"]["code"] == "missing_state"
    assert "library_match_feature" in payload["error"]["required_tools"]


def test_loading_reports_what_was_loaded(fresh_session):
    text = tools.library_load()
    assert "lib.msp" in text
    assert session_state.session.library.store is not None


def test_the_payload_carries_no_coordinate_arrays(fresh_session, monkeypatch):
    """座標はセッションに持つ。戻り値に点列を載せない（文脈を食うため）。"""
    tools.library_load()
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])
    text = tools.library_match_feature(104.0706, ion_mode="positive")

    # 測定側の点列が戻り値に出ていないこと。座標はセッションにだけ持つ。
    assert "69.0335" not in text
    assert "87.0441" not in text
    # 一方でセッションには完全な座標がある（図保存ツールがここから読む）。
    # store のレコードは `.msp` 由来でも **m/z 昇順**（Task 4）。
    stored = session_state.session.library.last_match["candidates"][0]["spectrum"]
    assert stored == [[69.0335, 500.0], [87.0441, 999.0]]


def test_loading_does_not_disturb_the_other_slots(fresh_session):
    before = session_state.session.arf
    tools.library_load()
    assert session_state.session.arf is before


def test_missing_dcl_reports_not_found_without_confusing_it_for_no_match(fresh_session):
    """`.dcl` に MS/MS が無いときは `not_found`。「合わなかった」と読めない文言にする。"""
    tools.library_load()
    payload = json.loads(tools.library_match_feature(104.0706, ion_mode="positive"))
    assert payload["status"] == "not_found"
    assert "未取得" in payload["message"]


def test_mirror_plot_passes_ms2_tol_so_measured_side_gets_colored(fresh_session, monkeypatch):
    """Task 8 で直したバグの再発防止: ms2_tol を渡し忘れると測定側の色分けが消える。"""
    tools.library_load()
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])
    tools.library_match_feature(104.0706, ion_mode="positive")

    captured = {}
    from metabolomix.plots import mirror as mirror_plot
    real_build = mirror_plot.build_mirror_payload

    def spy(*args, **kwargs):
        captured.update(kwargs)
        return real_build(*args, **kwargs)

    monkeypatch.setattr(mirror_plot, "build_mirror_payload", spy)
    tools.library_plot_mirror(output="payload")
    assert captured.get("ms2_tol") is not None


def test_candidates_are_ranked_by_a_documented_key(fresh_session, monkeypatch):
    """rank 列が何の降順かを payload 自身が明示する（docstring だけに頼らない）。"""
    tools.library_load()
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])
    payload = json.loads(tools.library_match_feature(104.0706, ion_mode="positive"))
    assert payload["ranked_by"] == "total_score"
    assert "total_score" in payload["candidates_table"].splitlines()[0].split("	")


def test_second_load_does_not_rescan_the_source_file(fresh_session, monkeypatch):
    """SQLite store をキャッシュした目的そのもの: 2 回目の library_load は
    ライブラリファイルを再走査してはいけない（cache hit のはず）。"""
    tools.library_load()  # 1 回目でキャッシュを構築する

    from metabolomix.library import dbs as dbs_reader
    from metabolomix.library import msp as msp_reader

    def _boom(*_args, **_kwargs):
        raise AssertionError("iter_records が呼ばれた（2 回目は cache hit のはず）")

    monkeypatch.setattr(msp_reader, "iter_records", _boom)
    monkeypatch.setattr(dbs_reader, "iter_records", _boom)

    text = tools.library_load()  # 2 回目: 元ファイルへは触れないはず
    payload = json.loads(text)
    assert payload["status"] == "success"


# --------------------------------------------------------------------------
# Important 7（最終レビュー）: `.dcl` 側の測定スペクトル検索経路。
# 既存のテストは `_measured_spectrum` を全部 monkeypatch していて、
# `.dcl` 解決 → `get_msms_by_precursor` → hits の選択が一度も検証されて
# いなかった（この穴が Important 7 を通した）。ここでは実際の `.dcl` バイナリを
# 組み立てて、その経路を通しで確認する。レイアウトは
# `metabolomix/dcl/reader.py` 冒頭のドキュメントに準拠（`tests/test_parser_decode.py`
# の `dcl_bytes`/`dcl_result` と同じ組み立て方だが、fixture はテスト自身が
# 作る規約に沿ってここに複製する）。
# --------------------------------------------------------------------------
def _dcl_result(precursor_mz, rt, spectrum, scan_id=0, raw_spec_id=1,
                ion_mode=0, model_height=9000.0, sn=50.0):
    return {
        "scan_id": scan_id, "raw_spec_id": raw_spec_id, "precursor_mz": precursor_mz,
        "ion_mode": ion_mode, "rt": rt, "model_height": model_height, "sn": sn,
        "spectrum": spectrum,
    }


def _dcl_bytes(results):
    header = b"DC" + struct.pack("<i", 1) + b"\x00" + struct.pack("<i", len(results))
    table_size = 8 * len(results)
    body = b""
    pointers = []
    for result in results:
        pointers.append(len(header) + table_size + len(body))
        spectrum = result["spectrum"]
        body += struct.pack(
            "<qiididdd d",
            0, result["scan_id"], result["raw_spec_id"], result["precursor_mz"],
            result["ion_mode"], result["rt"], 0.0, 0.0, result["precursor_mz"],
        )
        body += struct.pack("<5d", 0.0, result["model_height"], 0.0, 0.0, 0.0)
        body += struct.pack("<5f", 0.0, 0.0, 0.0, result["sn"], 0.0)
        body += struct.pack("<3i", len(spectrum), 0, 0)
        for mass, intensity in spectrum:
            body += struct.pack("<2di", mass, intensity, 0)
    return header + b"".join(struct.pack("<q", p) for p in pointers) + body


def test_measured_spectrum_picks_the_dcl_hit_closest_in_rt(fresh_session):
    """Important 7 (1/2) の再発防止: 同じ precursor m/z で複数の `.dcl` ヒットが
    あるとき、ファイル内の出現順の先頭（`hits[0]`）ではなく RT が最も近いものを
    選ぶこと。以前は先頭を無条件に採っており、rt_tol の窓が広いと別ピークの
    測定スペクトルを黙って採点しうる不具合があった。"""
    far_spectrum = [[999.0, 1.0]]
    near_spectrum = [[111.0, 500.0]]
    (fresh_session / "sample.dcl").write_bytes(_dcl_bytes([
        _dcl_result(precursor_mz=760.5851, rt=3.0, spectrum=far_spectrum, scan_id=0),
        _dcl_result(precursor_mz=760.5851, rt=12.0, spectrum=near_spectrum, scan_id=1),
    ]))

    measured = tools._measured_spectrum(760.5851, rt=11.8, mz_tol=0.01, rt_tol=100.0)
    assert measured == near_spectrum


def test_library_match_feature_walks_the_real_dcl_lookup_path(fresh_session):
    """Important 7 の穴そのものへの回帰テスト: 既存テストは全部
    `_measured_spectrum` を monkeypatch していたため、`.dcl` 解決 →
    `deserialize_dcl` → `get_msms_by_precursor` → hits の選択という経路が
    1 度も通っていなかった。ここでは monkeypatch せず `library_match_feature`
    を実際の `.dcl` に対して走らせ、RT が近いヒットの測定スペクトルが
    採点に使われることまで確認する。"""
    tools.library_load()
    (fresh_session / "sample.dcl").write_bytes(_dcl_bytes([
        _dcl_result(precursor_mz=104.0706, rt=3.0, spectrum=[[999.0, 1.0]], scan_id=0),
        _dcl_result(precursor_mz=104.0706, rt=12.0,
                    spectrum=[[87.0441, 999.0], [69.0335, 500.0]], scan_id=1),
    ]))

    payload = json.loads(tools.library_match_feature(104.0706, rt=12.0, ion_mode="positive"))
    assert payload["status"] == "success"
    assert payload["measured_peak_count"] == 2  # RT=12.0 の方（RT=3.0 は 1 本だけ）


def test_dcl_lookup_window_is_decoupled_from_library_search_params(fresh_session, monkeypatch):
    """Important 7 (2/2) の再発防止: ライブラリ照合用の許容幅
    （`store.search_params` 由来、`.dbs` の `RtTolerance` 既定は実質 RT フィルタ
    無効化の 100.0）を `.dcl` 側の測定スペクトル検索へ流用しない。以前は
    `library_match_feature` が解決した `resolved_mz_tol`/`resolved_rt_tol` を
    そのまま `_measured_spectrum` に渡していたため、広いライブラリ用許容幅の
    `.dbs` を読ませると `.dcl` 側の窓まで連動して広がっていた。"""
    tools.library_load()
    real_store = session_state.session.library.store

    class _StoreWithWideSearchParams:
        """`.dbs` の実質 RT フィルタ無効化（RtTolerance=100.0）を模した store。"""
        def summary(self):
            s = dict(real_store.summary())
            s["search_params"] = {
                "ms1_tolerance": 50.0, "ms2_tolerance": 0.05, "rt_tolerance": 100.0,
            }
            return s

        def candidates(self, *args, **kwargs):
            return real_store.candidates(*args, **kwargs)

    session_state.session.library.store = _StoreWithWideSearchParams()

    captured = {}

    def spy(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        return [[87.0441, 999.0], [69.0335, 500.0]]

    monkeypatch.setattr(tools, "_measured_spectrum", spy)
    tools.library_match_feature(104.0706, ion_mode="positive")

    # 呼び出しに mz_tol/rt_tol を一切渡していない（= _measured_spectrum 自身の
    # 固定既定 0.01/0.2 を使う。ライブラリ側の 50.0/100.0 を持ち込まない）。
    assert "mz_tol" not in captured["kwargs"]
    assert "rt_tol" not in captured["kwargs"]


def test_reloading_the_library_closes_the_old_store(fresh_session):
    """Minor 8 の再発防止: `library_load` を繰り返すと sqlite3 接続が漏れる。
    `LibraryStore.close()` は最初から存在したが、古い store を上書きする前に
    呼んでいなかった。閉じた store への問い合わせは `sqlite3.ProgrammingError`
    になることを利用して、上書き前にちゃんと閉じたことを確認する。"""
    tools.library_load()
    old_store = session_state.session.library.store

    tools.library_load()  # 同じファイル・同じ sha256 でも再読み込みは起こる

    assert session_state.session.library.store is not old_store
    with pytest.raises(Exception):
        old_store.candidates(104.0706, mz_tol=0.01)  # 閉じた接続へのクエリはここで失敗するはず


# --------------------------------------------------------------------------
# 順位付け。MS-DIAL の総合スコア（`GetTotalScore`）に揃えてある。
# weighted 単独だと、同名・同 precursor の別レコードで実用上の最良候補が
# 1 位に来ないことがある（2026-09-22 の spike。実データ 120 feature で
# top-1 一致 88.3% → 94.2%）。
# --------------------------------------------------------------------------
_TWO_RECORDS = textwrap.dedent("""\
    NAME: NOISY
    PRECURSORMZ: 104.0706
    IONMODE: Positive
    Num Peaks: 2
    87.0441 999
    69.0335 500

    NAME: CLEAN
    PRECURSORMZ: 104.0706
    IONMODE: Positive
    Num Peaks: 2
    87.0441 999
    69.0335 500
""")


def _fake_scores(by_name):
    """レコード名ごとに固定スコアを返す `match_spectrum` の差し替え。

    採点そのものは `test_spectral_match.py` が縛っているので、ここでは
    「どのスコアで並べるか」だけを見る。
    """
    def fake(measured, reference, **kwargs):
        # reference の中身では名前が分からないので、呼ばれた順に割り当てる。
        scores = by_name[fake.calls]
        fake.calls += 1
        return {**scores, "entropy_similarity": 0.0, "alignment": []}
    fake.calls = 0
    return fake


def test_a_higher_weighted_score_no_longer_wins_on_its_own(fresh_session, monkeypatch):
    """実データで見つかった逆転（PS 753.5430）の再現。1 件目は weighted だけが高く、
    2 件目は他の全指標で優る。総合スコアなら 2 件目が 1 位になる。"""
    (fresh_session / "lib.msp").write_text(_TWO_RECORDS, encoding="utf-8")
    tools.library_load()
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])
    monkeypatch.setattr(tools, "match_spectrum", _fake_scores([
        {"simple_dot_product": 0.70, "weighted_dot_product": 0.99,
         "reverse_dot_product": 0.76, "matched_peaks_percentage": 0.20,
         "matched_peaks_count": 18},
        {"simple_dot_product": 0.97, "weighted_dot_product": 0.98,
         "reverse_dot_product": 0.99, "matched_peaks_percentage": 1.00,
         "matched_peaks_count": 16},
    ]))

    payload = json.loads(tools.library_match_feature(104.0706, ion_mode="positive"))
    rows = payload["candidates_table"].splitlines()
    header = rows[0].split("\t")
    assert rows[1].split("\t")[header.index("name")] == "CLEAN"
    assert rows[2].split("\t")[header.index("name")] == "NOISY"


def test_the_scoring_rules_used_are_echoed_once(fresh_session, monkeypatch):
    """RT 項を入れたかどうかは `.dbs` の `IsUseTimeForAnnotationScoring` 次第で、
    戻り値だけを見ても分からない。呼び出し側が再現できるよう 1 回だけ載せる。"""
    tools.library_load()
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])
    payload = json.loads(tools.library_match_feature(104.0706, rt=1.0, ion_mode="positive"))

    # `.msp` には注釈スコアリングのフラグが無い → 上流既定の False に倣い RT 項は入らない。
    assert payload["scoring"]["use_rt"] is False
    assert payload["scoring"]["ms1_tol"] == pytest.approx(0.01)


def test_the_candidate_table_keeps_the_subterms_out_of_the_payload(fresh_session, monkeypatch):
    """内訳（rt_similarity / mass_similarity / spectrum_score）はセッション側に持つ。
    表に出すと候補数ぶん掛け算で効いて戻り値が膨らむ。"""
    tools.library_load()
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])
    payload = json.loads(tools.library_match_feature(104.0706, ion_mode="positive"))

    header = payload["candidates_table"].splitlines()[0].split("\t")
    assert "rt_similarity" not in header
    assert "mass_similarity" not in header

    scores = session_state.session.library.last_match["candidates"][0]["scores"]
    assert "total_score" in scores
    assert "rt_similarity" in scores
    assert "mass_similarity" in scores


def test_the_no_candidates_message_does_not_leak_float32_noise(fresh_session, monkeypatch):
    """`.dbs` の許容幅は C# の 32bit float 由来で、64bit へ広げると
    `0.009999999776482582` になる。`round_floats()` は payload 構造の中の float しか
    辿らないので、**文字列へ焼いた数値には届かない**——補間箇所で書式を指定する。"""
    tools.library_load()
    store = session_state.session.library.store
    monkeypatch.setattr(store, "summary", lambda: {"search_params": {
        "ms1_tolerance": 0.009999999776482582, "ms2_tolerance": 0.02500000037252903,
        "rt_tolerance": 2.0,
    }})
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 500.0]])

    payload = json.loads(tools.library_match_feature(104.0706, ion_mode="negative"))
    assert payload["status"] == "no_candidates"
    assert "±0.01," in payload["message"] or "±0.01 " in payload["message"]
    assert "0.00999999" not in payload["message"]


# --------------------------------------------------------------------------
# 採点対象外のピーク。採点は `normalize_measured` で足切りした後のスペクトルに
# 対して行うのに、図は足切り前の生ピークを描いている。この差を図の上で
# 見分けられるようにするための受け渡し（`.dbs` の cutoff が非ゼロのときだけ現れる）。
# --------------------------------------------------------------------------
def _store_with_cutoffs(monkeypatch, **cutoffs):
    store = session_state.session.library.store
    params = {"ms1_tolerance": 0.01, "ms2_tolerance": 0.025, "rt_tolerance": 0.2}
    params.update(cutoffs)
    monkeypatch.setattr(store, "summary", lambda: {"search_params": params})


def test_peaks_dropped_by_the_cutoff_are_recorded_for_the_figure(fresh_session, monkeypatch):
    tools.library_load()
    _store_with_cutoffs(monkeypatch, relative_amp_cutoff=0.1)
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 5.0]])

    tools.library_match_feature(104.0706, ion_mode="positive")

    assert session_state.session.library.last_match["unscored_mz"] == [69.0335]


def test_nothing_is_marked_unscored_when_the_cutoffs_are_zero(fresh_session, monkeypatch):
    """既定（足切り 0）の run では従来と何も変わらないこと。"""
    tools.library_load()
    _store_with_cutoffs(monkeypatch, relative_amp_cutoff=0.0, absolute_amp_cutoff=0.0)
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 5.0]])

    tools.library_match_feature(104.0706, ion_mode="positive")

    assert session_state.session.library.last_match["unscored_mz"] == []


def test_mirror_plot_passes_the_unscored_peaks_so_the_layer_reaches_the_figure(
        fresh_session, monkeypatch):
    """`ms2_tol` と同じ轍（ツール層で渡し忘れると図から静かに消える）を踏まない。"""
    tools.library_load()
    _store_with_cutoffs(monkeypatch, relative_amp_cutoff=0.1)
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 5.0]])
    tools.library_match_feature(104.0706, ion_mode="positive")

    captured = {}
    from metabolomix.plots import mirror as mirror_plot
    real_build = mirror_plot.build_mirror_payload

    def spy(*args, **kwargs):
        captured.update(kwargs)
        return real_build(*args, **kwargs)

    monkeypatch.setattr(mirror_plot, "build_mirror_payload", spy)
    tools.library_plot_mirror(output="payload")

    assert captured.get("unscored_mz") == [69.0335]


def test_the_caption_reports_the_dropped_peaks_because_the_figure_alone_is_easy_to_miss(
        fresh_session, monkeypatch):
    tools.library_load()
    _store_with_cutoffs(monkeypatch, relative_amp_cutoff=0.1)
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 5.0]])
    tools.library_match_feature(104.0706, ion_mode="positive")

    caption = tools.library_plot_mirror()[0]

    assert "1/2" in caption and "cutoff" in caption


def test_the_caption_stays_quiet_when_every_peak_was_scored(fresh_session, monkeypatch):
    tools.library_load()
    _store_with_cutoffs(monkeypatch, relative_amp_cutoff=0.0)
    monkeypatch.setattr(tools, "_measured_spectrum",
                        lambda *a, **k: [[87.0441, 999.0], [69.0335, 5.0]])
    tools.library_match_feature(104.0706, ion_mode="positive")

    caption = tools.library_plot_mirror()[0]

    assert "cutoff" not in caption


# --------------------------------------------------------------------------
# 研究室ライブラリ（pos / neg の 2 ファイルを環境変数で指す）。
# --------------------------------------------------------------------------
_MSP_NEG_NO_MODE = textwrap.dedent("""\
    NAME: Glutamate
    PRECURSORMZ: 146.0459
    Num Peaks: 1
    102.0561 999
""")


def _lab_libraries(tmp_path, monkeypatch):
    outside = tmp_path / "outside"
    outside.mkdir()
    pos = outside / "lab_pos.msp"
    neg = outside / "lab_neg.msp"
    pos.write_text(_MSP, encoding="utf-8")
    neg.write_text(_MSP_NEG_NO_MODE, encoding="utf-8")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(pos))
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(neg))
    return pos, neg


def test_library_load_picks_the_env_var_for_the_polarity(fresh_session, monkeypatch):
    (fresh_session / "lib.msp").unlink()
    _lab_libraries(fresh_session, monkeypatch)

    payload = json.loads(tools.library_load(ion_mode="negative"))

    assert payload["status"] == "success"
    assert payload["file"] == "lab_neg.msp"
    assert "outside" not in json.dumps(payload, ensure_ascii=False)   # 置き場所は出さない


def test_library_load_asks_for_a_polarity_when_both_are_configured(fresh_session, monkeypatch):
    (fresh_session / "lib.msp").unlink()
    _lab_libraries(fresh_session, monkeypatch)

    payload = json.loads(tools.library_load())

    assert payload["status"] == "error"
    assert payload["code"] == "MSP_AMBIGUOUS"
    assert session_state.session.library.store is None


def test_library_load_reports_an_env_var_pointing_nowhere(fresh_session, monkeypatch):
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(fresh_session / "gone.msp"))
    payload = json.loads(tools.library_load(ion_mode="negative"))
    assert payload["status"] == "error"
    assert payload["code"] == "MSP_ENV_NOT_FOUND"


def test_records_without_a_polarity_are_announced_and_still_matched(fresh_session, monkeypatch):
    """IONMODE 欄の無い neg ライブラリでも ion_mode="negative" の照合で候補が出る。"""
    (fresh_session / "lib.msp").unlink()
    _lab_libraries(fresh_session, monkeypatch)

    payload = json.loads(tools.library_load(ion_mode="negative"))
    assert payload["records_without_ion_mode"] == 1
    assert "極性" in payload["note"]

    monkeypatch.setattr(tools, "_measured_spectrum", lambda *a, **k: [[102.0561, 999.0]])
    match = json.loads(tools.library_match_feature(146.0459, ion_mode="negative"))
    assert match["n_candidates"] == 1


def test_non_utf8_lines_are_announced(fresh_session):
    (fresh_session / "lib.msp").write_bytes(
        "NAME: グルタミン酸\nPRECURSORMZ: 148.06\nNum Peaks: 0\n".encode("cp932"))
    payload = json.loads(tools.library_load())
    assert payload["non_utf8_lines"] == 1
    assert "UTF-8" in payload["note"]
