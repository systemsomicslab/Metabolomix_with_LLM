"""playbook ノートの `tools:` が実在の登録済み MCP tool を指すことを保証する。

playbook は個々のツールを再ドキュメントせず「触る道具」を宣言するだけなので、
ツールの改名・削除で参照が腐る。frontmatter は実ツール名と機械照合できる
（knowledge と違い人手レビュー不要）ため、ここで自動検出する。
"""

import asyncio
import unittest

import server
from metabolomix.corpus import knowledge_store as ks


class PlaybookToolConsistencyTests(unittest.TestCase):
    def test_all_playbook_tools_are_registered(self):
        tools = asyncio.run(server.mcp.list_tools())
        registered = {t.name for t in tools}

        refs = ks.playbook_tool_references(server.PLAYBOOK_DIR)
        problems = {
            slug: [t for t in used if t not in registered]
            for slug, used in refs.items()
        }
        problems = {slug: missing for slug, missing in problems.items() if missing}

        self.assertEqual(
            problems,
            {},
            f"playbook が未登録のツールを参照しています: {problems}\n"
            f"登録済み: {sorted(registered)}",
        )


if __name__ == "__main__":
    unittest.main()
