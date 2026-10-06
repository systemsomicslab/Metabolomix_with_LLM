"""dcl_parser ツールと、MS/MS を同定検証へ統合する経路の検証。

MS/MS は同定確度の最強証拠だが、これまで .dcl は MCP 面から到達できず、
verify_peak_annotation は精密質量・アダクト・S/N だけで確度を語っていた。
ここでは (1) dcl_parser で MS/MS を読めること (2) PAI2 ピークへ索引対応で
付与されること (3) 検証ドシエが実スペクトルとフラグのみを区別することを固定する。
"""
import asyncio
import json
import tempfile
import unittest
from pathlib import Path

from metabolomix.core import mcp_core
import server
from metabolomix.core import session_state
from metabolomix.pai2 import tools as tools_pai2

from tests.dcl_fixture import build_dcl_bytes

SPECTRUM = [(255.2330, 900.0), (281.2486, 1500.0), (283.2643, 400.0)]


def _write_dcl(directory: Path, name: str, results) -> Path:
    path = directory / name
    path.write_bytes(build_dcl_bytes(results))
    return path


def _extract_json(text: str) -> dict:
    """ツール出力（見出し＋JSON＋トピック誘導1行）から JSON 本体だけを取り出す。"""
    decoder = json.JSONDecoder()
    return decoder.raw_decode(text[text.index("{"):])[0]


class DclToolRegistrationTests(unittest.TestCase):
    def test_dcl_parser_is_registered(self):
        names = {t.name for t in asyncio.run(server.mcp.list_tools())}
        self.assertIn("dcl_parser", names)


class DclParserToolTests(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_summarizes_msms_inventory(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_dcl(Path(tmp), "s.dcl", [
                {"precursor_mz": 100.0, "rt": 1.0, "spectrum": []},
                {"precursor_mz": 885.5499, "rt": 12.3, "spectrum": SPECTRUM},
            ])
            out = server.dcl_parser(file_path=str(path))
        payload = _extract_json(out)
        self.assertEqual(payload["summary"]["total_results"], 2)
        self.assertEqual(payload["summary"]["with_msms"], 1)

    def test_points_at_its_own_output_format_topic(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_dcl(Path(tmp), "s.dcl", [{"precursor_mz": 100.0, "spectrum": SPECTRUM}])
            out = server.dcl_parser(file_path=str(path))
        self.assertIn("lipidmix://docs/output-format/dcl", out)

    def test_missing_file_reports_error_not_exception(self):
        out = server.dcl_parser(file_path=str(Path("no_such_dir") / "missing.dcl"))
        self.assertIn("見つかりません", out)

    def test_resolver_picks_latest_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            _write_dcl(directory, "AlignmentResult_2024_01_01_00_00_00.dcl",
                       [{"precursor_mz": 1.0, "spectrum": []}])
            newer = _write_dcl(directory, "AlignmentResult_2026_05_15_10_13_35.dcl",
                               [{"precursor_mz": 2.0, "spectrum": SPECTRUM}])
            original = mcp_core.DATA_DIR
            mcp_core.DATA_DIR = directory
            self.addCleanup(setattr, mcp_core, "DATA_DIR", original)
            from metabolomix.core.path_resolvers import resolve_dcl_file_path
            self.assertEqual(Path(resolve_dcl_file_path()), newer)


class SiblingMsmsAttachmentTests(unittest.TestCase):
    """.pai2 のピークに、同名 .dcl の MS/MS を索引対応で付与する。"""

    def test_attaches_when_sibling_dcl_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            _write_dcl(directory, "run_A.dcl", [
                {"precursor_mz": 100.0, "spectrum": []},
                {"precursor_mz": 885.5499, "spectrum": SPECTRUM},
            ])
            features = [{"m/z": 100.0}, {"m/z": 885.5499}]
            report = tools_pai2._attach_sibling_msms(str(directory / "run_A.pai2"), features)
        self.assertEqual(report["attached"], 1)
        self.assertEqual(features[1]["n_msms_peaks"], 3)

    def test_reports_when_no_sibling_dcl(self):
        with tempfile.TemporaryDirectory() as tmp:
            features = [{"m/z": 100.0}]
            report = tools_pai2._attach_sibling_msms(str(Path(tmp) / "lonely.pai2"), features)
        self.assertEqual(report["attached"], 0)
        self.assertIsNotNone(report["caveat"])
        self.assertNotIn("msms_spectrum", features[0])

    def test_corrupt_dcl_does_not_break_pai2_parsing(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "run_B.dcl").write_bytes(b"NOPE" + b"\x00" * 40)
            features = [{"m/z": 100.0}]
            report = tools_pai2._attach_sibling_msms(str(directory / "run_B.pai2"), features)
        self.assertEqual(report["attached"], 0)
        self.assertIn("caveat", report)


class VerificationDossierMsmsTests(unittest.TestCase):
    """ドシエは「実スペクトル」と「取得フラグのみ」を区別して示す。"""

    # m/z は C42H82NO8P + [M+HCOO]- の理論値。MSI Level 2 の他の条件（名称・精密質量整合）
    # を満たさないと MS/MS の寄与を検証できないため、質量は PASS する値を使う。
    BASE = {
        "id": 1, "name": "PC 34:1", "ontology": "PC", "formula": "C42H82NO8P",
        "adduct": "[M+HCOO]-", "m/z": 804.5760, "ion_mode": "Negative",
        "time": {"rt": 5.0},
    }

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_real_spectrum_is_surfaced_as_evidence(self):
        feat = dict(self.BASE, has_msms=True, n_msms_peaks=3,
                    msms_spectrum=[list(p) for p in SPECTRUM])
        dossier = server._build_verification_dossier(feat, {})
        evidence = dossier["analytical_checks"]["msms"]
        self.assertEqual(evidence["band"], "PASS")
        self.assertEqual(evidence["n_peaks"], 3)
        self.assertEqual(evidence["top_fragments"][0][0], 281.2486)

    def test_flag_only_is_disclosed_in_dossier_and_msi_rationale(self):
        feat = dict(self.BASE, has_msms=True)
        dossier = server._build_verification_dossier(feat, {})
        self.assertEqual(dossier["analytical_checks"]["msms"]["band"], "FLAG_ONLY")
        msi = dossier["identity_normalization"]["msi"]
        self.assertIn("フラグ", msi["rationale"])

    def test_absent_msms_is_not_claimed_as_evidence(self):
        feat = dict(self.BASE, has_msms=False)
        dossier = server._build_verification_dossier(feat, {})
        self.assertEqual(dossier["analytical_checks"]["msms"]["band"], "ABSENT")
        self.assertNotEqual(dossier["identity_normalization"]["msi"]["level"], 2)

    def test_deterministic_summary_mentions_msms(self):
        feat = dict(self.BASE, has_msms=True, n_msms_peaks=3,
                    msms_spectrum=[list(p) for p in SPECTRUM])
        dossier = server._build_verification_dossier(feat, {})
        self.assertIn("msms", dossier["llm_decision"]["deterministic_summary"])


if __name__ == "__main__":
    unittest.main()
