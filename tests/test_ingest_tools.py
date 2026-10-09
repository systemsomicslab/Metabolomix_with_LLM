"""文献探索ツール群の登録と knowledge_coverage のスモーク検証。"""

import asyncio
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from metabolomix.core import mcp_core
import server


class IngestToolRegistrationTests(unittest.TestCase):
    def test_new_tools_registered(self):
        names = {t.name for t in asyncio.run(server.mcp.list_tools())}
        for expected in (
            "knowledge_coverage",
            "paper_search",
            "ingest_stage",
            "ingest_review_queue",
            "ingest_promote",
            "ingest_reject",
        ):
            self.assertIn(expected, names)

    def test_inbox_resource_registered(self):
        uris = {str(r.uri) for r in asyncio.run(server.mcp.list_resources())}
        self.assertIn("lipidmix://knowledge/inbox", uris)


class SlugContainmentToolTests(unittest.TestCase):
    """ツール引数の slug は LLM 由来（材料に非信頼な抄録を含む）。置き場の外を触らせない。"""

    BAD_SLUGS = ("../../evil", "..\\evil", "sub/evil", "", "..")

    def test_ingest_reject_refuses_traversal(self):
        for bad in self.BAD_SLUGS:
            with self.subTest(slug=bad):
                out = server.ingest_reject(bad)
                self.assertIn("不正", out)

    def test_ingest_promote_refuses_traversal(self):
        for bad in self.BAD_SLUGS:
            with self.subTest(slug=bad):
                out = server.ingest_promote(bad)
                self.assertIn("不正", out)

    def test_record_objective_refuses_traversal(self):
        out = server.record_objective(
            analysis_id="../../evil", dataset="d", polarity="NEG",
            groups=["a", "b"], comparison="a vs b", sub_questions=["Q"],
        )
        self.assertIn("不正", out)


class KnowledgeCoverageSmokeTests(unittest.TestCase):
    """objective は tmp に自前で作る。

    objective は解析フォルダ配下 reports/ に住む解析記録でリポジトリには無い。
    ユーザ環境に残った記録の有無でテストの成否が変わらないよう、対象の objective は
    tmp の解析フォルダにここで作る。
    knowledge/ は追跡済みの種ノートを参照するので実ディレクトリのまま使う。
    """

    ANALYSIS_ID = "test-neg-lipidome-trt-vs-ctrl"

    def setUp(self):
        self._orig_data = mcp_core.DATA_DIR
        self._tmp = tempfile.TemporaryDirectory()
        mcp_core.DATA_DIR = Path(self._tmp.name)
        self._env = mock.patch.dict("os.environ")
        self._env.start()
        os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        server.record_objective(
            analysis_id=self.ANALYSIS_ID,
            dataset="NEG / 2_lipidome_lcms",
            polarity="NEG",
            groups=["control", "treatment"],
            comparison="群間の脂質クラスプロファイル差",
            sub_questions=[
                "どの脂質クラスに群間差が最も出るか",
                "エーテル脂質（PE P-/PC P-）に一方向の減少があるか",
                "その差は注釈の信頼度（P-/O- の取り違え）で説明できてしまわないか",
            ],
        )

    def tearDown(self):
        self._env.stop()
        mcp_core.DATA_DIR = self._orig_data
        self._tmp.cleanup()

    def test_coverage_on_sample_objective(self):
        # frontmatter の analysis_id で解決される
        out = server.knowledge_coverage(self.ANALYSIS_ID)
        # 3つの小問が状態付きで返る
        self.assertIn("どの脂質クラスに群間差", out)
        self.assertIn("エーテル脂質", out)
        self.assertTrue(any(state in out for state in ("COVERED", "WEAK", "GAP")))

    def test_missing_objective(self):
        out = server.knowledge_coverage("no-such-analysis-id")
        self.assertIn("見つかりません", out)


if __name__ == "__main__":
    unittest.main()
