"""差次的結果のエクスポート契約（spec §6.1）。"""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import server
from lipidmix.core import session_state

# ブリーフは `server.session_state...` と書くが、server は wildcard import
# （from lipidmix.arf.tools import *）経由で __all__ に無い session_state を
# 再エクスポートしていないため AttributeError になる（test_differential_tools.py の
# 既存テストと同じ既知の落とし穴）。直接 import した session_state を参照する。


_EXPECTED_META_KEYS = [
    "# contract_version", "# exported_at",
    "# source_arf", "# source_arf2",
    "# group_a", "# group_b", "# log2fc_sign", "# q_threshold",
    "# preprocess", "# n_features_total",
    "# msi_level は .arf2 由来の注釈確度。MS/MS の有無ではない",
]


class TestExportDifferential(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        session_state.session.arf.last_differential = {
            "kind": "two_group", "a": "gf_AIN", "b": "gf_HFD",
            "n_a": 5, "n_b": 5,
            "q_threshold": 0.05, "log2fc_threshold": 1.0,
            "contract_version": 1,
            "log2fc_sign": "positive means group_b is higher",
            "log_transform": True,
            "results": [
                {"feature": "Spot_1_height", "mean_a": 10.0, "mean_b": 40.0,
                 "log2fc": 1.9, "p": 0.001, "q": 0.01},
                {"feature": "Spot_2_height", "mean_a": 5.0, "mean_b": 5.0,
                 "log2fc": 0.0, "p": 0.9, "q": 0.95},
            ],
            "volcano": [],
        }
        self.catalog = [
            {"MasterAlignmentID": 1, "Name": "PC 34:1", "Ontology": "PC",
             "InChIKey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "SMILES": "CCO",
             "MassCenter": 760.5851, "RT": 12.34},
            {"MasterAlignmentID": 2, "Name": "Unknown", "Ontology": "",
             "InChIKey": "", "SMILES": "", "MassCenter": 100.0, "RT": 1.0},
        ]

    def test_missing_state_without_differential(self):
        session_state.session.arf.last_differential = None
        payload = json.loads(server.arf_export_differential("out.tsv"))
        self.assertEqual(payload["error"]["code"], "missing_state")
        self.assertIn("arf_differential", payload["error"]["required_tools"])

    def test_fails_when_sibling_arf2_missing(self):
        with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=None):
            payload = json.loads(server.arf_export_differential("out.tsv"))
        self.assertEqual(payload["error"]["code"], "missing_state")
        self.assertEqual(payload["error"]["state"], "sibling_arf2")

    def test_rejects_stale_log2fc_contract_before_writing(self):
        """旧符号の状態へ新契約ラベルを付けて書き出してはならない。"""
        session_state.session.arf.last_differential["log2fc_sign"] = (
            "positive means group_a is higher")
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2):
                payload = json.loads(server.arf_export_differential(str(out)))
            self.assertFalse(out.exists())
        self.assertEqual(payload["error"]["code"], "missing_state")
        self.assertEqual(payload["error"]["state"], "compatible_two_group_differential")

    def test_rejects_export_when_no_row_has_inchikey(self):
        """consumer が拒否する本文 0 行のファイルを成功扱いで残さない。"""
        self.catalog[0]["InChIKey"] = ""
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                payload = json.loads(server.arf_export_differential(str(out)))
            self.assertFalse(out.exists())
        self.assertEqual(payload["status"], "error")
        self.assertEqual(payload["n_with_inchikey"], 0)
        self.assertIn("InChIKey", payload["message"])

    def test_writes_meta_and_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                payload = json.loads(server.arf_export_differential(str(out)))
            text = out.read_text(encoding="utf-8")

        self.assertEqual(payload["status"], "success")
        self.assertEqual(payload["n_with_inchikey"], 1)
        self.assertEqual(payload["n_unannotated"], 1)
        self.assertIn("# contract_version = 1", text)
        self.assertIn("# group_a = gf_AIN", text)
        self.assertIn("# log2fc_sign = positive means group_b is higher", text)
        self.assertIn("n_unannotated = 1", text)
        header = [l for l in text.splitlines() if not l.startswith("#")][0]
        self.assertEqual(header.split("\t")[0], "spot_id")
        self.assertIn("inchikey", header.split("\t"))
        body = [l for l in text.splitlines()
                if not l.startswith("#")][1:]
        self.assertEqual(len(body), 1)
        self.assertIn("AAAAAAAAAAAAAA-BBBBBBBBBB-C", body[0])
        cells = dict(zip(header.split("\t"), body[0].split("\t")))
        self.assertEqual(cells["mz"], "760.5851")
        self.assertEqual(cells["rt"], "12.3400")
        self.assertEqual(cells["log2fc"], "1.900000")
        self.assertEqual(cells["msi_level"], "3")

    def test_significant_flag_uses_stored_thresholds(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                server.arf_export_differential(str(out))
            row = [l for l in out.read_text(encoding="utf-8").splitlines()
                   if not l.startswith("#")][1]
        # q=0.01 <= 0.05 かつ |log2fc|=1.9 >= 1.0 なので significant
        self.assertEqual(row.split("\t")[-1], "true")

    def test_nonfinite_cells_are_written_blank(self):
        """下流の massbank-context 側 float() は 'nan' 文字列も受理してしまうため、
        NaN・inf は None と同じく空欄で書き出す（コントローラ裁定）。
        小 n・分散 0・全欠損では Welch t 検定・BH 補正が普通に NaN を返す
        （arf_differential の caveats 自身が p=NaN を想定している）。
        """
        session_state.session.arf.last_differential["results"] = [
            {"feature": "Spot_1_height", "mean_a": 10.0, "mean_b": float("nan"),
             "log2fc": float("nan"), "p": float("nan"), "q": float("inf")},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                server.arf_export_differential(str(out))
            text = out.read_text(encoding="utf-8")

        header = [l for l in text.splitlines() if not l.startswith("#")][0].split("\t")
        row = [l for l in text.splitlines() if not l.startswith("#")][1].split("\t")
        cells = dict(zip(header, row))
        self.assertEqual(cells["mean_b"], "")
        self.assertEqual(cells["log2fc"], "")
        self.assertEqual(cells["p_value"], "")
        self.assertEqual(cells["q_value"], "")
        self.assertEqual(cells["mean_a"], "10")
        # NaN/inf は非有意扱い（q<=threshold の判定ができないため）
        self.assertEqual(cells["significant"], "false")

    def test_arf_export_meta_line_order_is_frozen(self):
        """メタ行の順序と文面を固定する。

        Task 5 Step 6 で export_contract.build_meta() へ寄せるとき、この順序が
        変わっていないことを保証する。行の順序も下流との契約の一部。
        """
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                server.arf_export_differential(str(out))
            text = out.read_text(encoding="utf-8")
        meta = [l for l in text.splitlines() if l.startswith("#")]
        keys = [l.split(" = ")[0].split("\t")[0] for l in meta]
        self.assertEqual(keys, _EXPECTED_META_KEYS)

    def test_log2fc_sign_meta_line_is_verbatim(self):
        """massbank-context 側の契約リーダはこの文字列と逐語一致しないと拒否する。"""
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                server.arf_export_differential(str(out))
            text = out.read_text(encoding="utf-8")
        meta_lines = [l for l in text.splitlines() if l.startswith("#")]
        self.assertIn(
            "# log2fc_sign = positive means group_b is higher",
            meta_lines,
        )

    def test_no_curation_flags_means_no_curation_meta_line(self):
        """フラグ 0 件のときの出力は現行と完全に同じ(curation メタ行が無い)。"""
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                server.arf_export_differential(str(out))
            text = out.read_text(encoding="utf-8")
        self.assertFalse(any(l.startswith("# curation") for l in text.splitlines()))

    def test_wrong_flag_drops_the_row_and_declares_curation_applied(self):
        from lipidmix.curation import flags as curation_flags

        # spot_id=1 と 2 の両方に InChIKey を付け、1 だけを wrong にする
        # (2 だけになれば「消えた」ことを行の残り方で確認できる)。setUp の
        # results には Spot_1_height / Spot_2_height が既にある。
        self.catalog[1]["InChIKey"] = "DDDDDDDDDDDDDD-EEEEEEEEEE-F"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
                [{"spot_id": 1, "flag": "wrong"}],
                alignment=curation_flags.alignment_key(arf2),
                review_id="r", source="user")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                payload = json.loads(server.arf_export_differential(str(out)))
            text = out.read_text(encoding="utf-8")
        self.assertEqual(payload["status"], "success")
        self.assertEqual(payload["n_with_inchikey"], 1)
        body = [l for l in text.splitlines() if not l.startswith("#")][1:]
        self.assertEqual(len(body), 1)
        self.assertIn("DDDDDDDDDDDDDD-EEEEEEEEEE-F", body[0])
        meta_lines = [l for l in text.splitlines() if l.startswith("# curation")]
        self.assertEqual(len(meta_lines), 1)
        self.assertTrue(meta_lines[0].startswith("# curation = applied\t"))
        self.assertIn("curation_wrong_excluded = 1", meta_lines[0])


if __name__ == "__main__":
    unittest.main()
