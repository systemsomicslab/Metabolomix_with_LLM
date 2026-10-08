import json
import unittest
from unittest import mock

import numpy as np

import server
from metabolomix.analysis import export_contract
from metabolomix.core import session_state
from metabolomix.arf import tools as tools_arf


def _run_preprocess_and_differential():
    """arf_differential の前準備（前処理後行列 + サンプルメタ）を session に積む。

    既存テスト（TestArfDifferential 内）と同じ 3+3 サンプルの下地。呼び出し側が
    続けて server.arf_differential(group_a="A", group_b="B") を呼ぶ想定。
    """
    session_state.session.arf.feature_matrix = np.array([
        [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
        [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
    ])
    names = ["a1", "a2", "a3", "b1", "b2", "b3"]
    session_state.session.arf.pp_sample_names = names
    session_state.session.arf.pp_feature_names = ["f0", "f1"]
    session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
    session_state.session.arf.sample_meta = {
        n: {"group": ("A" if n.startswith("a") else "B")} for n in names
    }


class TestArfDifferential(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_requires_matrix(self):
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(out["error"]["code"], "missing_state")

    def test_two_group_reports_significant_and_caveats(self):
        session_state.session.arf.feature_matrix = np.array([
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
        ])
        names = ["a1", "a2", "a3", "b1", "b2", "b3"]
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("A" if n.startswith("a") else "B"),
                "batch": ("d1" if n.startswith("a") else "d2")}
            for n in names
        }
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(out["status"], "success")
        self.assertGreaterEqual(out["summary"]["n_significant"], 1)
        self.assertTrue(any("バッチ" in c or "交絡" in c for c in out["caveats"]))

    def test_two_group_payload_omits_full_volcano_but_session_keeps_it(self):
        # 8000字切り詰めで summary が埋没しないよう、payload は全量 volcano を含まず
        # 要約中心にする。全量は session に残し save_figure(kind="volcano") から使える。
        session_state.session.arf.feature_matrix = np.array([
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
        ])
        names = ["a1", "a2", "a3", "b1", "b2", "b3"]
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("A" if n.startswith("a") else "B")} for n in names
        }
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertNotIn("volcano", out)
        self.assertIn("volcano_note", out)
        self.assertIn("summary", out)
        # 全量 volcano（全特徴分）は session に残る
        vol = session_state.session.arf.last_differential["volcano"]
        self.assertEqual(len(vol), 2)

    def test_discloses_zero_manual_exclusions_when_none_applied(self):
        # #6: arf_parser/load_dataset の再実行で excluded_samples/excluded_spots が
        # reset_analysis() によって黙ってゼロ化されても、差次的解析の caveats が
        # 「除外ゼロ件」を積極的に述べていれば、ユーザーは外れ値込みの統計を
        # 「除外済みの解析と同じもの」と誤解しない。
        session_state.session.arf.feature_matrix = np.array([
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
        ])
        names = ["a1", "a2", "a3", "b1", "b2", "b3"]
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("A" if n.startswith("a") else "B")} for n in names
        }
        self.assertEqual(session_state.session.arf.excluded_samples, set())
        self.assertEqual(session_state.session.arf.excluded_spots, set())
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertTrue(
            any("現在なし" in c for c in out["caveats"]),
            f"ゼロ件の手動除外が開示されていない: {out['caveats']}",
        )

    def test_discloses_active_manual_exclusions(self):
        session_state.session.arf.feature_matrix = np.array([
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
        ])
        names = ["a1", "a2", "a3", "b1", "b2", "b3"]
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("A" if n.startswith("a") else "B")} for n in names
        }
        session_state.session.arf.excluded_samples.add("a3")
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertTrue(
            any("サンプル 1 件" in c for c in out["caveats"]),
            f"除外の開示が無い: {out['caveats']}",
        )

    def test_volcano_note_points_to_structured_plot_tool(self):
        # 既存の test_two_group_* と同じ 3+3 サンプルの下地を使う（群内 n>=2 を満たし、
        # 「n 不足」caveat 経路に入らない構成）。
        session_state.session.arf.feature_matrix = np.array([
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
        ])
        names = ["a1", "a2", "a3", "b1", "b2", "b3"]
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("A" if n.startswith("a") else "B")} for n in names
        }
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertIn("arf_plot_volcano", out["volcano_note"])
        self.assertNotIn('save_figure(kind="volcano") で図示', out["volcano_note"])
        self.assertNotIn("volcano", out)
        last = session_state.session.arf.last_differential
        self.assertEqual(last["q_threshold"], 0.05)
        self.assertEqual(last["log2fc_threshold"], 1.0)
        self.assertEqual(last["n_a"], 3)
        self.assertEqual(last["n_b"], 3)


class TestArfDifferentialDegenerate(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def _prime(self, groups, matrix=None):
        n = len(groups)
        if matrix is None:
            matrix = np.arange(1, n * 2 + 1, dtype=float).reshape(n, 2)
        session_state.session.arf.feature_matrix = matrix
        names = [f"s{i}" for i in range(n)]
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {names[i]: {"group": groups[i]} for i in range(n)}

    def test_empty_requested_group_is_flagged(self):
        # 群名が1件も一致しないなら success を返してはいけない。有意0件を
        # 「群間差なし」と読ませる余地を残さず、利用可能な Class ID を示して落とす。
        self._prime(["A", "A", "A"])  # no B members at all
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(out["status"], "error")
        self.assertIn("B", out["message"])
        self.assertIn("A", out["message"])

    def test_zero_testable_features_is_flagged(self):
        # A and B each have 2 members but all values identical -> zero variance ->
        # every feature p=NaN -> n_tested 0. Must warn, not read as "no differences".
        self._prime(["A", "A", "B", "B"], matrix=np.ones((4, 2)))
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(out["summary"]["n_tested"], 0)
        self.assertTrue(any("検定" in c for c in out["caveats"]),
                        f"expected zero-tested caveat, got {out['caveats']}")


class TestSaveVolcano(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_requires_last_differential(self):
        out = server.save_figure("volcano", "A1")
        self.assertIn("error", out.lower())

    def test_writes_png(self):
        session_state.session.arf.last_differential = {
            "kind": "two_group", "a": "A", "b": "B",
            "volcano": [
                {"feature": "f0", "log2fc": 2.0, "neg_log10_p": 3.0, "sig": "up"},
                {"feature": "f1", "log2fc": -2.0, "neg_log10_p": 3.0, "sig": "down"},
                {"feature": "f2", "log2fc": 0.0, "neg_log10_p": 0.1, "sig": "ns"},
            ],
        }
        rel = server.save_figure("volcano", "A1")
        self.assertIn("volcano", rel)



class TestPooledGroupSpecs(unittest.TestCase):
    """Class ID は因子トークンの連結（24M_GF_F）。プール群比較を可能にする。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()
        # 2 週齢 x 2 菌叢、各 2 個体。加齢で f0 が上がる。
        classes = ["24M_GF", "24M_GF", "24M_SPF", "24M_SPF",
                   "9w_GF", "9w_GF", "9w_SPF", "9w_SPF"]
        names = [f"s{i}" for i in range(len(classes))]
        session_state.session.arf.feature_matrix = np.array([
            [50.0, 5.0], [52.0, 5.1], [48.0, 4.9], [51.0, 5.0],
            [10.0, 5.0], [11.0, 5.2], [9.5, 4.8], [10.5, 5.1],
        ])
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["Spot_0_height", "Spot_1_height"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": c} for n, c in zip(names, classes)}

    def test_partial_token_pools_class_ids(self):
        out = json.loads(server.arf_differential(group_a="24M", group_b="9w"))
        self.assertEqual(out["status"], "success")
        self.assertEqual(out["n_a"], 4)
        self.assertEqual(out["n_b"], 4)
        self.assertEqual(out["resolved_class_ids"]["group_a"], ["24M_GF", "24M_SPF"])
        self.assertEqual(out["summary"]["n_significant"], 1)
        self.assertTrue(any("プール群として解決" in c for c in out["caveats"]))

    def test_multi_token_spec_narrows_the_pool(self):
        out = json.loads(server.arf_differential(group_a="24M_GF", group_b="9w_GF"))
        self.assertEqual(out["n_a"], 2)
        self.assertEqual(out["resolved_class_ids"]["group_a"], ["24M_GF"])

    def test_overlapping_specs_are_rejected(self):
        # "GF" と "24M" は 24M_GF を共有する。プールが排他でないので検定してはいけない。
        out = json.loads(server.arf_differential(group_a="GF", group_b="24M"))
        self.assertEqual(out["status"], "error")
        self.assertIn("24M_GF", out["message"])

    def test_identical_specs_are_rejected(self):
        out = json.loads(server.arf_differential(group_a="GF", group_b="GF"))
        self.assertEqual(out["status"], "error")

    def test_unknown_spec_lists_available_class_ids(self):
        out = json.loads(server.arf_differential(group_a="24M", group_b="99w"))
        self.assertEqual(out["status"], "error")
        self.assertIn("99w", out["message"])
        self.assertIn("24M_GF", out["message"])


class TestTopHitAnnotation(unittest.TestCase):
    """ARF 側 Name が Unknown でも ARF2 の注釈で解釈可能にする（両者は食い違い得る）。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()
        names = ["a1", "a2", "a3", "b1", "b2", "b3"]
        session_state.session.arf.feature_matrix = np.array([
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
        ])
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["Spot_474_height", "Spot_1_height"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = {
            n: {"group": ("A" if n.startswith("a") else "B")} for n in names}

    def test_uses_arf_name_when_annotated(self):
        session_state.session.arf.features = [
            {"MasterAlignmentID": 474, "Name": "BMP 42:10"},
            {"MasterAlignmentID": 1, "Name": "Unknown"},
        ]
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        top = out["summary"]["top"][0]
        self.assertEqual(top["spot_id"], 474)
        self.assertEqual(top["name"], "BMP 42:10")
        self.assertEqual(top["name_source"], "arf")

    def test_falls_back_to_arf2_catalog(self):
        import tempfile
        from pathlib import Path

        session_state.session.arf.features = [{"MasterAlignmentID": 474, "Name": "Unknown"}]
        arf2_spots = [{"MasterAlignmentID": 474,
                       "Name": "SL 33:0;O|SL 17:0;O/16:0", "Ontology": "SL"}]
        # 実ファイルを置く。読み出しは load_catalog 経由（パス+mtime+サイズをキーに
        # 再パースを避けるキャッシュ）で、実在しないパスは stat の時点で落ちる。
        with tempfile.TemporaryDirectory() as tmp:
            sibling = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            sibling.write_bytes(b"")
            with mock.patch.object(tools_arf, "_sibling_arf2_path", return_value=sibling), \
                 mock.patch("metabolomix.arf2.reader.deserialize", return_value=arf2_spots):
                out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        top = out["summary"]["top"][0]
        self.assertEqual(top["name"], "SL 33:0;O|SL 17:0;O/16:0")
        self.assertEqual(top["ontology"], "SL")
        self.assertEqual(top["name_source"], "arf2")
        self.assertTrue(any("ARF2" in c for c in out["caveats"]))

    def test_no_sibling_arf2_leaves_name_none(self):
        session_state.session.arf.features = [{"MasterAlignmentID": 474, "Name": "Unknown"}]
        with mock.patch.object(tools_arf, "_sibling_arf2_path", return_value=None):
            out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertIsNone(out["summary"]["top"][0]["name"])



class TestSiblingArf2Resolution(unittest.TestCase):
    """MasterAlignmentID はアラインメント実行ごとに振り直される。別バッチの .arf2 を
    引くと ID 対応が黙って崩れるため、同一語幹の兄弟だけを許す。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_matches_same_alignment_stem(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            stem = "AlignmentResult_2026_07_09_17_14_19"
            arf = os.path.join(d, f"{stem}_PeakProperties.arf")
            arf2 = os.path.join(d, f"{stem}.arf2")
            open(arf, "wb").close()
            open(arf2, "wb").close()
            session_state.session.arf.current_file_path = arf
            self.assertEqual(str(tools_arf._sibling_arf2_path()), arf2)

    def test_rejects_other_batch(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            arf = os.path.join(d, "AlignmentResult_2026_07_09_17_14_19_PeakProperties.arf")
            open(arf, "wb").close()
            # 別実行の .arf2 しか無い場合は掴まない
            open(os.path.join(d, "AlignmentResult_2026_07_09_17_34_57.arf2"), "wb").close()
            session_state.session.arf.current_file_path = arf
            self.assertIsNone(tools_arf._sibling_arf2_path())

    def test_no_loaded_file(self):
        self.assertIsNone(tools_arf._sibling_arf2_path())


if __name__ == "__main__":
    unittest.main()


class TestDifferentialSampleSelection(unittest.TestCase):
    """比較対象は生体試料のみ。交絡判定は「実際に比較した2群」に対して行う。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def _prime(self, meta, matrix):
        names = list(meta)
        session_state.session.arf.feature_matrix = np.asarray(matrix, dtype=float)
        session_state.session.arf.pp_sample_names = names
        session_state.session.arf.pp_feature_names = ["f0", "f1"]
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        session_state.session.arf.sample_meta = meta

    def test_qc_sample_is_not_counted_into_a_compared_group(self):
        # QC の Class ID が group_a のトークンを含むと、素通しではプールに紛れ込む。
        meta = {
            "a1": {"group": "A", "role": "sample", "batch": "d1"},
            "a2": {"group": "A", "role": "sample", "batch": "d1"},
            "a3": {"group": "A", "role": "sample", "batch": "d1"},
            "b1": {"group": "B", "role": "sample", "batch": "d1"},
            "b2": {"group": "B", "role": "sample", "batch": "d1"},
            "b3": {"group": "B", "role": "sample", "batch": "d1"},
            "QC_1": {"group": "A", "role": "qc", "batch": "d1"},
        }
        self._prime(meta, [
            [10.0, 5.0], [11.0, 5.1], [9.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [48.0, 4.8],
            [30.0, 5.0],
        ])
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(out["status"], "success")
        self.assertEqual(out["n_a"], 3, "QC が群 A に数え込まれている")
        self.assertEqual(out["n_b"], 3)
        self.assertTrue(
            any("QC" in c for c in out["caveats"]),
            f"QC 除外の開示が無い: {out['caveats']}",
        )

    def test_confounding_is_assessed_after_pooling(self):
        # Class ID 単位では各群が単一バッチ（細粒度では交絡に見える）が、
        # プール後の A / B はどちらも d1+d2 を含むので交絡していない。
        meta = {
            "a1": {"group": "A_x", "role": "sample", "batch": "d1"},
            "a2": {"group": "A_x", "role": "sample", "batch": "d1"},
            "a3": {"group": "A_y", "role": "sample", "batch": "d2"},
            "a4": {"group": "A_y", "role": "sample", "batch": "d2"},
            "b1": {"group": "B_x", "role": "sample", "batch": "d1"},
            "b2": {"group": "B_x", "role": "sample", "batch": "d1"},
            "b3": {"group": "B_y", "role": "sample", "batch": "d2"},
            "b4": {"group": "B_y", "role": "sample", "batch": "d2"},
        }
        self._prime(meta, [
            [10.0, 5.0], [11.0, 5.1], [10.5, 4.9], [10.2, 5.0],
            [50.0, 5.0], [52.0, 5.2], [51.0, 4.8], [50.5, 5.1],
        ])
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(out["status"], "success")
        self.assertEqual({out["n_a"], out["n_b"]}, {4})
        self.assertFalse(
            any(c.startswith("交絡:") for c in out["caveats"]),
            f"プール後は交絡していないのに交絡と報告した: {out['caveats']}",
        )

    def test_confounding_still_detected_when_pooled_groups_are_confounded(self):
        meta = {
            "a1": {"group": "A_x", "role": "sample", "batch": "d1"},
            "a2": {"group": "A_x", "role": "sample", "batch": "d1"},
            "a3": {"group": "A_y", "role": "sample", "batch": "d1"},
            "b1": {"group": "B_x", "role": "sample", "batch": "d2"},
            "b2": {"group": "B_x", "role": "sample", "batch": "d2"},
            "b3": {"group": "B_y", "role": "sample", "batch": "d2"},
        }
        self._prime(meta, [
            [10.0, 5.0], [11.0, 5.1], [10.5, 4.9],
            [50.0, 5.0], [52.0, 5.2], [51.0, 4.8],
        ])
        out = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertTrue(
            any(c.startswith("交絡:") for c in out["caveats"]),
            f"群⟂バッチが完全交絡なのに報告されていない: {out['caveats']}",
        )


class TestDifferentialSignDisclosure(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_payload_and_state_disclose_sign(self):
        """符号の向きが payload と session 状態の両方に出ること。

        向きは 2026-08-31 に反転したため、古い出力と区別できないと
        解釈が静かに逆転する。
        """
        _run_preprocess_and_differential()   # 同ファイルの既存ヘルパに合わせる
        payload = json.loads(server.arf_differential(group_a="A", group_b="B"))
        self.assertEqual(payload["differential_contract_version"], 1)
        self.assertIn("group_b", payload["log2fc_sign"])
        # ブリーフは `server.session_state...` と書くが、server は wildcard import
        # （from metabolomix.arf.tools import *）経由で __all__ に無い session_state を
        # 再エクスポートしていないため AttributeError になる。本ファイルの既存テスト
        # （例: line 122 の last = session_state.session.arf.last_differential）に
        # 倣い、直接 import した session_state を参照する。
        state = session_state.session.arf.last_differential
        self.assertEqual(state["contract_version"], 1)


class TestDifferentialContractLinkage(unittest.TestCase):
    """contract_version / log2fc_sign が export_contract を単一情報源として
    参照していること（値の一致ではなく紐付け自体）を確認する。

    `== 1` のような値一致テストは、別々にリテラルを2箇所へ複製していても
    今日は両方とも 1 なので通ってしまい、意味がない。ここでは
    export_contract.CONTRACT_VERSION を書き換えて、arf_differential が
    その変更に追随することを見る。追随しなければ、CONTRACT_VERSION を
    上げた瞬間に arf_export_differential の互換性チェックへ永遠に
    引っかかる（再現ループ）。
    """
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_stored_contract_version_follows_export_contract_constant(self):
        _run_preprocess_and_differential()
        with mock.patch.object(export_contract, "CONTRACT_VERSION", 2):
            json.loads(server.arf_differential(group_a="A", group_b="B"))
        state = session_state.session.arf.last_differential
        self.assertEqual(state["contract_version"], 2)

    def test_stored_log2fc_sign_follows_export_contract_constant(self):
        _run_preprocess_and_differential()
        patched_sign = "positive means group_a is higher (patched)"
        with mock.patch.object(export_contract, "LOG2FC_SIGN", patched_sign):
            json.loads(server.arf_differential(group_a="A", group_b="B"))
        state = session_state.session.arf.last_differential
        self.assertEqual(state["log2fc_sign"], patched_sign)
