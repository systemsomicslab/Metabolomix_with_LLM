"""peak_verification 純ロジックの検証。"""

import textwrap
import unittest

from metabolomix.msdial import peak_verification as pv


class FormulaMassTests(unittest.TestCase):
    def test_parse_formula_multi_element_and_multi_digit(self):
        self.assertEqual(
            pv.parse_formula("C42H82NO8P"),
            {"C": 42, "H": 82, "N": 1, "O": 8, "P": 1},
        )

    def test_parse_formula_single_atom_defaults_to_one(self):
        self.assertEqual(pv.parse_formula("CH4"), {"C": 1, "H": 4})

    def test_parse_formula_empty_raises(self):
        with self.assertRaises(ValueError):
            pv.parse_formula("")

    def test_monoisotopic_mass_pc_34_1(self):
        mass = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        self.assertAlmostEqual(mass, 759.5778, places=2)

    def test_monoisotopic_mass_unknown_element_raises(self):
        with self.assertRaises(KeyError):
            pv.monoisotopic_mass({"Xx": 1})


class AdductMassErrorTests(unittest.TestCase):
    def test_adduct_mz_protonated(self):
        neutral = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        self.assertAlmostEqual(pv.adduct_mz(neutral, "[M+H]+"), 760.5851, places=2)

    def test_adduct_mz_ammonium(self):
        self.assertAlmostEqual(pv.adduct_mz(100.0, "[M+NH4]+"), 118.03383, places=4)

    def test_adduct_mz_deprotonated(self):
        self.assertAlmostEqual(pv.adduct_mz(100.0, "[M-H]-"), 98.99272, places=4)

    def test_adduct_mz_unknown_returns_none(self):
        self.assertIsNone(pv.adduct_mz(100.0, "[M+ZZ]+"))

    def test_mass_error_pass(self):
        result = pv.mass_error_ppm(760.5851, "C42H82NO8P", "[M+H]+")
        self.assertEqual(result["band"], "PASS")
        self.assertAlmostEqual(result["ppm"], 0.0, delta=5.0)

    def test_mass_error_borderline(self):
        neutral = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        theo = pv.adduct_mz(neutral, "[M+H]+")
        observed = theo * (1 + 8e-6)  # +8 ppm
        result = pv.mass_error_ppm(observed, "C42H82NO8P", "[M+H]+")
        self.assertEqual(result["band"], "BORDERLINE")

    def test_mass_error_fail(self):
        neutral = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        theo = pv.adduct_mz(neutral, "[M+H]+")
        observed = theo * (1 + 20e-6)  # +20 ppm
        result = pv.mass_error_ppm(observed, "C42H82NO8P", "[M+H]+")
        self.assertEqual(result["band"], "FAIL")

    def test_mass_error_unknown_formula(self):
        result = pv.mass_error_ppm(760.0, "Unknown", "[M+H]+")
        self.assertEqual(result["band"], "UNKNOWN")
        self.assertIsNone(result["ppm"])

    def test_mass_error_unknown_adduct(self):
        result = pv.mass_error_ppm(760.0, "C42H82NO8P", "Unknown")
        self.assertEqual(result["band"], "UNKNOWN")


class AdductConsistencyTests(unittest.TestCase):
    def test_polarity_match_positive(self):
        result = pv.adduct_consistency("[M+H]+", "Positive", "PC")
        self.assertTrue(result["polarity_ok"])
        self.assertEqual(result["band"], "PASS")
        self.assertTrue(result["class_typical"])

    def test_polarity_mismatch(self):
        result = pv.adduct_consistency("[M+H]+", "Negative", "PC")
        self.assertFalse(result["polarity_ok"])
        self.assertEqual(result["band"], "FAIL")

    def test_class_atypical_adduct_is_advisory_not_fail(self):
        result = pv.adduct_consistency("[M+Na]+", "Positive", "PC")
        self.assertTrue(result["polarity_ok"])
        self.assertEqual(result["band"], "PASS")
        self.assertFalse(result["class_typical"])

    def test_unknown_class_gives_none_typical(self):
        result = pv.adduct_consistency("[M+H]+", "Positive", "ZZZ")
        self.assertIsNone(result["class_typical"])

    def test_unknown_adduct_is_unknown_band(self):
        result = pv.adduct_consistency("Unknown", "Positive", "PC")
        self.assertEqual(result["band"], "UNKNOWN")
        self.assertIsNone(result["polarity_ok"])

    def test_unrecognized_ion_mode_is_unknown_not_fail(self):
        result = pv.adduct_consistency("[M+H]+", "Both", "PC")
        self.assertEqual(result["band"], "UNKNOWN")
        self.assertIsNone(result["polarity_ok"])
        self.assertTrue(result["class_typical"])  # class advisory still computed


class ClassTokenTests(unittest.TestCase):
    def test_extract_class_token_prefers_ontology(self):
        self.assertEqual(pv.extract_class_token("PC 34:1", "PC"), "pc")

    def test_extract_class_token_falls_back_to_name(self):
        self.assertEqual(pv.extract_class_token("TG 52:2", ""), "tg")

    def test_extract_class_token_empty(self):
        self.assertIsNone(pv.extract_class_token("", ""))

    def test_ether_caveat_detected_for_o_prefix(self):
        caveats = pv.ether_caveats("PE O-38:5", "PE")
        self.assertTrue(caveats)
        self.assertIn("P-/O-", caveats[0])

    def test_ether_caveat_detected_for_p_prefix(self):
        self.assertTrue(pv.ether_caveats("PC P-36:4", "PC"))

    def test_no_ether_caveat_for_diacyl(self):
        self.assertEqual(pv.ether_caveats("PC 34:1", "PC"), [])

    def test_vocab_hits_match(self):
        vocab = {"pc": ["phosphatidylcholine"], "tg": ["triacylglycerol"]}
        self.assertEqual(pv.vocab_hits("pc", vocab), ["phosphatidylcholine"])

    def test_vocab_hits_none_token(self):
        self.assertEqual(pv.vocab_hits(None, {"pc": ["x"]}), [])


class TestExpandedTables(unittest.TestCase):
    def test_dimer_adduct_mz(self):
        # [2M-H]- : (2*neutral - proton) / 1
        neutral = 100.0
        mz = pv.adduct_mz(neutral, "[2M-H]-")
        self.assertAlmostEqual(mz, 2 * 100.0 - pv.PROTON_MASS, places=4)

    def test_doubly_charged_divides(self):
        neutral = 800.0
        mz = pv.adduct_mz(neutral, "[M-2H]2-")
        self.assertAlmostEqual(mz, (800.0 - 2 * pv.PROTON_MASS) / 2, places=4)

    def test_new_elements_present(self):
        for el in ("D", "F", "Br", "C13"):
            self.assertIn(el, pv.ELEMENT_MASSES)

    def test_singly_charged_unchanged(self):
        # 既存の1価挙動が数値的に不変であることを保証。
        self.assertAlmostEqual(pv.adduct_mz(500.0, "[M+H]+"), 500.0 + pv.PROTON_MASS, places=6)
        self.assertAlmostEqual(pv.adduct_mz(500.0, "[M+HCOO]-"), 500.0 + 44.99820286, places=6)


import json
import tempfile
from pathlib import Path

import server
from metabolomix.core import mcp_core
from metabolomix.core import session_state
from metabolomix.pai2.reader import IonMode


def _feat(**over):
    base = {
        "id": 1,
        "name": "PC 34:1",
        "ontology": "PC",
        "formula": "C42H82NO8P",
        "adduct": "[M+H]+",
        "m/z": 760.5851,
        "ion_mode": IonMode.Positive,
        "time": {"rt": 12.34},
        "S/N": 42.0,
    }
    base.update(over)
    return base


class VerifyPeakToolTests(unittest.TestCase):
    def setUp(self):
        self._orig_features = session_state.session.pai2.filtered_features
        self._orig_knowledge = mcp_core.KNOWLEDGE_DIR
        self._tmp = tempfile.TemporaryDirectory()
        mcp_core.KNOWLEDGE_DIR = Path(self._tmp.name) / "knowledge"
        mcp_core.KNOWLEDGE_DIR.mkdir()

    def tearDown(self):
        session_state.session.pai2.filtered_features = self._orig_features
        mcp_core.KNOWLEDGE_DIR = self._orig_knowledge
        self._tmp.cleanup()

    def test_error_when_not_loaded(self):
        session_state.session.pai2.filtered_features = None
        out = json.loads(server.verify_peak_annotation(peak_id="1"))
        self.assertEqual(out["error"]["code"], "missing_state")

    def test_error_when_no_selector(self):
        session_state.session.pai2.filtered_features = [_feat()]
        out = json.loads(server.verify_peak_annotation())
        self.assertEqual(out["status"], "error")

    def test_not_found(self):
        session_state.session.pai2.filtered_features = [_feat()]
        out = json.loads(server.verify_peak_annotation(peak_id="999"))
        self.assertEqual(out["status"], "not_found")

    def test_success_shape_and_bands(self):
        session_state.session.pai2.filtered_features = [_feat()]
        out = json.loads(server.verify_peak_annotation(peak_id="1"))
        self.assertEqual(out["status"], "success")
        self.assertEqual(out["identity"]["name"], "PC 34:1")
        self.assertEqual(out["identity"]["ion_mode"], "Positive")
        self.assertEqual(out["analytical_checks"]["mass_error"]["band"], "PASS")
        self.assertEqual(out["analytical_checks"]["adduct_consistency"]["band"], "PASS")
        self.assertIn("class_token", out["biological_plausibility"])

    def test_unknown_formula_degrades_to_unknown_band(self):
        session_state.session.pai2.filtered_features = [_feat(formula="Unknown")]
        out = json.loads(server.verify_peak_annotation(peak_id="1"))
        self.assertEqual(out["analytical_checks"]["mass_error"]["band"], "UNKNOWN")

    def test_ether_caveat_surfaced(self):
        session_state.session.pai2.filtered_features = [_feat(name="PE O-38:5", ontology="PE", id=2)]
        out = json.loads(server.verify_peak_annotation(peak_id="2"))
        self.assertTrue(out["biological_plausibility"]["caveats"])

    def test_candidate_slugs_from_knowledge(self):
        note = mcp_core.KNOWLEDGE_DIR / "pe-p-vs-pe-o-annotation.md"
        note.write_text(
            "---\ntype: knowledge\n"
            "description: PE の P-/O- 表記の混同に注意\n"
            "claim_strength: established\n"
            "tags: [annotation, plasmalogen]\n---\n\n本文\n",
            encoding="utf-8",
        )
        session_state.session.pai2.filtered_features = [_feat(name="PE 38:5", ontology="PE", id=3)]
        out = json.loads(server.verify_peak_annotation(peak_id="3"))
        self.assertIn(
            "pe-p-vs-pe-o-annotation",
            out["biological_plausibility"]["candidate_knowledge_slugs"],
        )

    def test_multiple_matches_wrapped(self):
        session_state.session.pai2.filtered_features = [_feat(id=1), _feat(id=2)]
        out = json.loads(server.verify_peak_annotation(peak_name="PC"))
        self.assertEqual(out["status"], "success")
        self.assertIn("matches", out)
        self.assertEqual(len(out["matches"]), 2)


if __name__ == "__main__":
    unittest.main()


class MsmsEvidenceTests(unittest.TestCase):
    """MS/MS の「実スペクトルがある」と「取得フラグが立っている」を峻別する。

    PAI2 の has_msms は取得参照の有無を示すだけで msms_spectrum が非空とは限らない
    （output-format core §9-8）。MSI Level 2 は MS/MS を根拠にするので、この区別を
    曖昧にしたまま確度を主張してはならない。
    """

    def test_real_spectrum_passes_and_reports_top_fragments(self):
        feat = {
            "has_msms": True,
            "n_msms_peaks": 3,
            "msms_spectrum": [[255.2, 900.0], [281.2, 1500.0], [283.2, 400.0]],
        }
        out = pv.msms_evidence(feat)
        self.assertEqual(out["band"], "PASS")
        self.assertEqual(out["source"], "spectrum")
        self.assertEqual(out["n_peaks"], 3)
        # 強度降順で上位を返す（解釈に使うのは主要フラグメント）
        self.assertEqual(out["top_fragments"][0], [281.2, 1500.0])
        self.assertIsNone(out["caveat"])

    def test_flag_without_spectrum_is_flag_only_with_caveat(self):
        out = pv.msms_evidence({"has_msms": True, "msms_spectrum": [], "n_msms_peaks": 0})
        self.assertEqual(out["band"], "FLAG_ONLY")
        self.assertEqual(out["source"], "flag")
        self.assertIsNotNone(out["caveat"])
        self.assertEqual(out["top_fragments"], [])

    def test_no_flag_and_no_spectrum_is_absent(self):
        out = pv.msms_evidence({"has_msms": False})
        self.assertEqual(out["band"], "ABSENT")
        self.assertIsNone(out["source"])
        self.assertEqual(out["n_peaks"], 0)

    def test_top_fragments_are_capped(self):
        spectrum = [[100.0 + i, float(i)] for i in range(20)]
        out = pv.msms_evidence(
            {"has_msms": True, "n_msms_peaks": 20, "msms_spectrum": spectrum}, top_n=5
        )
        self.assertEqual(len(out["top_fragments"]), 5)
        self.assertEqual(out["n_peaks"], 20)  # 元本数は保つ


def test_the_msms_band_still_has_only_three_states():
    """PASS / FLAG_ONLY / ABSENT の 3 状態は契約。照合はその内側に足す。"""
    from metabolomix.msdial.peak_verification import msms_evidence
    assert msms_evidence({"msms_spectrum": [[100.0, 999.0]]})["band"] == "PASS"
    assert msms_evidence({"has_msms": True})["band"] == "FLAG_ONLY"
    assert msms_evidence({})["band"] == "ABSENT"


def test_spectral_match_is_absent_without_a_loaded_library():
    from metabolomix.msdial.peak_verification import msms_evidence
    assert msms_evidence({"msms_spectrum": [[100.0, 999.0]]}).get("spectral_match") is None


_MSP_GABA = textwrap.dedent("""\
    NAME: GABA
    PRECURSORMZ: 104.0706
    IONMODE: Positive
    Num Peaks: 2
    87.0441 999
    69.0335 500
""")


class _FakeIonModePositive:
    """PAI2 の `pai2.reader.IonMode.Positive` の代わり——`.name` だけ持つ。"""
    name = "Positive"


def _open_gaba_store(tmp_path, monkeypatch):
    from metabolomix.library import store as library_store
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    msp_path = tmp_path / "lib.msp"
    msp_path.write_text(_MSP_GABA, encoding="utf-8")
    return library_store.open_store(msp_path)


def test_spectral_match_finds_candidates_despite_ion_mode_case_mismatch(tmp_path, monkeypatch):
    """Critical 1 の再発防止: `verify_peak_annotation` の統合経路
    （`_spectral_match_for_feature`）は feature 側の `ion_mode.name`
    （`"Positive"`、大文字始まり）をそのまま store へ渡す。store 側は
    `"positive"`（小文字、`record.ION_MODES`）。以前は SQLite の既定 BINARY
    照合で一致せず、ライブラリを読み込んでいても常に `n_candidates=0` を返した
    （`no_candidates` は『ライブラリに無い化合物だ』という実質的な主張なので、
    この不一致は誤った結論に直結する）。"""
    lib_store = _open_gaba_store(tmp_path, monkeypatch)
    try:
        feature = {
            "m/z": 104.0706,
            "ion_mode": _FakeIonModePositive(),
            "time": {"rt": None},
        }
        spectrum = [[87.0441, 999.0], [69.0335, 500.0]]
        result = pv._spectral_match_for_feature(feature, spectrum, lib_store)
        assert result["status"] == "matched"
        assert result["n_candidates"] == 1
        assert result["best_match"]["name"] == "GABA"
    finally:
        lib_store.close()


def test_spectral_match_via_msms_evidence_end_to_end(tmp_path, monkeypatch):
    """`msms_evidence` から `_spectral_match_for_feature` までを通しで確認する
    （`session.library.store` が実際に読み込まれた状態で）。"""
    from metabolomix.core import session_state

    lib_store = _open_gaba_store(tmp_path, monkeypatch)
    try:
        monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
        session_state.session.library.store = lib_store
        feature = {
            "m/z": 104.0706,
            "ion_mode": _FakeIonModePositive(),
            "time": {"rt": None},
            "has_msms": True,
            "n_msms_peaks": 2,
            "msms_spectrum": [[87.0441, 999.0], [69.0335, 500.0]],
        }
        out = pv.msms_evidence(feature)
        assert out["band"] == "PASS"
        assert out["spectral_match"]["status"] == "matched"
        assert out["spectral_match"]["n_candidates"] == 1
    finally:
        lib_store.close()


def test_scoring_forwards_search_params_amplitude_cutoffs_to_match_spectrum(monkeypatch):
    """Important 2 の再発防止: `search_params` の `mass_range_begin` /
    `mass_range_end` / `relative_amp_cutoff` / `absolute_amp_cutoff` が
    `match_spectrum` に届くこと（検証 CLI `scripts/verify_spectral_match.py`
    と同じ集合）。以前は `ms2_tol` しか渡していなかった。"""
    captured = {}

    def fake_match_spectrum(measured, reference, *, ms2_tol, **kwargs):
        captured["ms2_tol"] = ms2_tol
        captured.update(kwargs)
        return {
            "simple_dot_product": 1.0, "weighted_dot_product": 1.0,
            "reverse_dot_product": 1.0, "matched_peaks_percentage": 1.0,
            "matched_peaks_count": 1, "entropy_similarity": 1.0, "alignment": [],
        }

    import metabolomix.analysis.spectral_match as spectral_match_module
    monkeypatch.setattr(spectral_match_module, "match_spectrum", fake_match_spectrum)

    class _FakeStore:
        def summary(self):
            return {"search_params": {
                "ms1_tolerance": 0.02, "ms2_tolerance": 0.05, "rt_tolerance": 1.0,
                "mass_range_begin": 10.0, "mass_range_end": 1000.0,
                "relative_amp_cutoff": 0.05, "absolute_amp_cutoff": 100.0,
            }}

        def candidates(self, *args, **kwargs):
            return [{"name": "X", "ontology": None, "adduct": None, "spectrum": [[100.0, 10.0]]}]

    feature = {"m/z": 100.0, "ion_mode": None, "time": {}}
    pv._spectral_match_for_feature(feature, [[100.0, 10.0]], _FakeStore())

    assert captured["mass_begin"] == 10.0
    assert captured["mass_end"] == 1000.0
    assert captured["relative_amp_cutoff"] == 0.05
    assert captured["absolute_amp_cutoff"] == 100.0


def test_tolerance_lookup_keeps_an_explicit_zero_from_search_params(monkeypatch):
    """Important 5 の再発防止: `search_params.get(key)` が `0.0`（正当な値。
    例えば足切りなし）でも、`or` 判定（偽値扱い）で既定値へ差し替えてはいけない。
    以前この関数だけが `or` を使っており、`metabolomix.library.tools._pick_tol`
    （`is None` 判定）とロジックがずれていた。"""
    captured = {}

    def fake_match_spectrum(measured, reference, *, ms2_tol, **kwargs):
        captured["ms2_tol"] = ms2_tol
        return {
            "simple_dot_product": 1.0, "weighted_dot_product": 1.0,
            "reverse_dot_product": 1.0, "matched_peaks_percentage": 1.0,
            "matched_peaks_count": 1, "entropy_similarity": 1.0, "alignment": [],
        }

    import metabolomix.analysis.spectral_match as spectral_match_module
    monkeypatch.setattr(spectral_match_module, "match_spectrum", fake_match_spectrum)

    class _FakeStore:
        def summary(self):
            # ms2_tolerance=0.0 は「明示的にゼロ」——`or` だと偽値扱いで
            # 既定値 (_LIBRARY_MS2_TOL=0.025) に化ける。
            return {"search_params": {"ms2_tolerance": 0.0}}

        def candidates(self, *args, **kwargs):
            return [{"name": "X", "ontology": None, "adduct": None, "spectrum": [[100.0, 10.0]]}]

    feature = {"m/z": 100.0, "ion_mode": None, "time": {}}
    pv._spectral_match_for_feature(feature, [[100.0, 10.0]], _FakeStore())

    assert captured["ms2_tol"] == 0.0


def test_the_best_match_is_chosen_by_the_same_rule_as_library_match_feature(monkeypatch):
    """2 つのツールが違う「最良候補」を返すと、ドシエと候補一覧が食い違う。
    順位付けは `library_match_feature` と同じ総合スコア（`GetTotalScore`）に揃える。"""
    scores = {
        "NOISY": {"simple_dot_product": 0.70, "weighted_dot_product": 0.99,
                  "reverse_dot_product": 0.76, "matched_peaks_percentage": 0.20,
                  "matched_peaks_count": 18},
        "CLEAN": {"simple_dot_product": 0.97, "weighted_dot_product": 0.98,
                  "reverse_dot_product": 0.99, "matched_peaks_percentage": 1.00,
                  "matched_peaks_count": 16},
    }

    def fake_match_spectrum(measured, reference, **kwargs):
        name = "NOISY" if reference == [[100.0, 10.0]] else "CLEAN"
        return {**scores[name], "entropy_similarity": 0.0, "alignment": []}

    import metabolomix.analysis.spectral_match as spectral_match_module
    monkeypatch.setattr(spectral_match_module, "match_spectrum", fake_match_spectrum)

    class _FakeStore:
        def summary(self):
            return {"search_params": {"ms1_tolerance": 0.01, "ms2_tolerance": 0.025,
                                      "rt_tolerance": 2.0}}

        def candidates(self, *args, **kwargs):
            return [
                {"name": "NOISY", "ontology": None, "adduct": None,
                 "precursor_mz": 100.0, "rt": 1.0, "spectrum": [[100.0, 10.0]]},
                {"name": "CLEAN", "ontology": None, "adduct": None,
                 "precursor_mz": 100.0, "rt": 1.0, "spectrum": [[100.0, 20.0]]},
            ]

    feature = {"m/z": 100.0, "ion_mode": None, "time": {"rt": 1.0}}
    out = pv._spectral_match_for_feature(feature, [[100.0, 10.0]], _FakeStore())

    assert out["best_match"]["name"] == "CLEAN"
    assert "total_score" in out["best_match"]


def test_the_best_match_survives_a_candidate_without_precursor_mz_or_rt(monkeypatch):
    """store のレコードが precursor m/z / RT を欠いていても照合を止めない
    （照合の失敗で band 判定を止めない、というこの関数の約束）。"""
    def fake_match_spectrum(measured, reference, **kwargs):
        return {"simple_dot_product": 1.0, "weighted_dot_product": 1.0,
                "reverse_dot_product": 1.0, "matched_peaks_percentage": 1.0,
                "matched_peaks_count": 1, "entropy_similarity": 1.0, "alignment": []}

    import metabolomix.analysis.spectral_match as spectral_match_module
    monkeypatch.setattr(spectral_match_module, "match_spectrum", fake_match_spectrum)

    class _FakeStore:
        def summary(self):
            return {"search_params": {}}

        def candidates(self, *args, **kwargs):
            return [{"name": "X", "ontology": None, "adduct": None, "spectrum": [[100.0, 10.0]]}]

    feature = {"m/z": 100.0, "ion_mode": None, "time": {}}
    out = pv._spectral_match_for_feature(feature, [[100.0, 10.0]], _FakeStore())

    assert out["status"] == "matched"
    assert out["best_match"]["name"] == "X"
