"""再編で静かに壊れる箇所を縛る特性テスト。

BASE_DIR 系は解決先を間違えても例外を出さず、空のディレクトリを新規作成して
黙って動く。既存テストは全部緑のまま蓄積ノートだけ見えなくなるため、
リポジトリルートを指し続けることをここで固定する。
"""
import os
import unittest
from pathlib import Path

from metabolomix.core import data_config
from metabolomix.core import mcp_core
from metabolomix.core import tool_helpers

REPO_ROOT = Path(__file__).resolve().parents[1]


class TestBaseDirResolution(unittest.TestCase):
    def test_base_dir_is_repo_root(self):
        self.assertEqual(mcp_core.BASE_DIR.resolve(), REPO_ROOT)

    def test_base_dir_has_repo_markers(self):
        # metabolomix/core/ を指してしまった場合をここで落とす
        for marker in ("docs", "playbook", "reference"):
            with self.subTest(marker=marker):
                self.assertTrue(
                    (mcp_core.BASE_DIR / marker).is_dir(),
                    f"{marker}/ が見つからない: BASE_DIR={mcp_core.BASE_DIR}",
                )
        self.assertTrue((mcp_core.BASE_DIR / "server.py").is_file())

    def test_output_format_doc_resolves(self):
        self.assertTrue(mcp_core.OUTPUT_FORMAT_DOC.is_file())

    def test_data_config_root_matches_mcp_core(self):
        # 循環回避のため両者は独立に root を計算する。一致することをここで縛る。
        self.assertEqual(
            data_config.DEFAULT_DATA_DIR.resolve(),
            (mcp_core.BASE_DIR / "data").resolve(),
        )


class TestRootLayout(unittest.TestCase):
    """ルート直下の .py は server.py と check.py の 2 つだけ、という鉄則を縛る。

    CLAUDE.md が明文化している規約だが、これまで散文にしかなかった。スクラッチや
    使い捨てスクリプトがルートに置き去りにされても誰も気づかず、`server.py` が
    薄いファサードであるという前提だけが静かに崩れる。
    """

    def test_root_python_files_are_exactly_the_two(self):
        found = {p.name for p in REPO_ROOT.glob("*.py")}
        self.assertEqual(
            found, {"server.py", "check.py"},
            "ルート直下の .py が増減している。実装は metabolomix/ に置く",
        )


class TestStateDirsOutsidePackage(unittest.TestCase):
    def test_state_dirs_are_not_inside_lipidmix(self):
        package_dir = (REPO_ROOT / "metabolomix").resolve()
        for name in ("KNOWLEDGE_DIR", "PLAYBOOK_DIR", "ANALYSES_DIR"):
            with self.subTest(name=name):
                resolved = getattr(mcp_core, name).resolve()
                self.assertNotIn(
                    package_dir,
                    resolved.parents,
                    f"{name} がパッケージ内を指している: {resolved}",
                )


class TestReferenceTables(unittest.TestCase):
    """`lipid_identity.load_reference_tables` は CWD 相対で reference/ を読む。

    移動では壊れないが、テストを常にリポジトリルートから実行する前提を固定する。
    表が空でも例外にならないので、行数で縛る。
    """

    def test_identity_tables_load_rows(self):
        tables = tool_helpers._identity_tables()
        self.assertEqual(len(tables["lipidmaps"]), 16)
        self.assertEqual(len(tables["refmet"]), 16)


class TestStateDirsAreNotCreatedOnImport(unittest.TestCase):
    """`metabolomix.core.mcp_core` の import が置き場ディレクトリを作らないことを縛る。

    KNOWLEDGE_DIR / PLAYBOOK_DIR / ANALYSES_DIR はパス解決の時点で mkdir していた
    ため、`analyses/` を消しても import のたび（＝pytest を回すたび、MCP サーバを
    起動するたび）に空ディレクトリが復活していた。置き場は書き込み時に作られれば
    十分で、解決は副作用を持ってはならない。
    """

    def test_import_does_not_create_state_dirs(self):
        import subprocess
        import sys
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            env = dict(os.environ)
            env["LIPIDMIX_KNOWLEDGE_DIR"] = str(base / "knowledge")
            env["LIPIDMIX_PLAYBOOK_DIR"] = str(base / "playbook")
            env["LIPIDMIX_ANALYSES_DIR"] = str(base / "analyses")
            proc = subprocess.run(
                [sys.executable, "-c", "import metabolomix.core.mcp_core"],
                cwd=str(REPO_ROOT), env=env,
                capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            created = sorted(p.name for p in base.iterdir())
            self.assertEqual(
                created, [],
                f"import だけで置き場が作られた: {created}",
            )
