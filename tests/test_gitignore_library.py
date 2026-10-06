"""参照ライブラリとその複製が、置き場所を問わず git に載らないこと。

研究室の参照ライブラリ（.msp）は外部流出禁止の資産。本体は環境変数で
リポジトリの外を指すが、`LIPIDMIX_LIBRARY_CACHE_DIR` をリポジトリ内の
`data/` 以外へ向けたり、`.msp` を手元でリポジトリ直下に置いたりしても、
`git add -A` で巻き込まれてはいけない。誰かが `.gitignore` を緩めたらここで落ちる。
"""
import shutil
import subprocess
from pathlib import Path

import pytest

from metabolomix.library import store

ROOT = Path(__file__).resolve().parents[1]

_MUST_BE_IGNORED = [
    "lib.msp",
    "somewhere/deep/lab_pos.msp",
    "somewhere/Dataset_Loaded.msp2",
    "somewhere/Dataset_Loaded.msp2.dbs",
    "somewhere/lipids.lbm",
    "somewhere/lipids.lbm2",
    "cache/lab_pos-0123456789abcdef.sqlite",
    "anywhere/.library-cache/lab_pos-0123456789abcdef.sqlite",
    f"cache/{store.DIGEST_INDEX_NAME}",
    # 設定ファイルには研究室ライブラリのパスが入る（spec 2026-10-06）
    "lipidmix.local.toml",
]


def _git_available() -> bool:
    if shutil.which("git") is None:
        return False
    probe = subprocess.run(["git", "rev-parse", "--is-inside-work-tree"], cwd=ROOT,
                           capture_output=True, text=True)
    return probe.returncode == 0


@pytest.mark.skipif(not _git_available(), reason="git の作業ツリーでない（配備先など）")
@pytest.mark.parametrize("relpath", _MUST_BE_IGNORED)
def test_library_files_and_their_copies_are_ignored_anywhere(relpath):
    result = subprocess.run(["git", "check-ignore", "--no-index", "-q", relpath], cwd=ROOT)
    assert result.returncode == 0, f"{relpath} が .gitignore で無視されていない"


@pytest.mark.skipif(not _git_available(), reason="git の作業ツリーでない（配備先など）")
def test_no_library_file_is_already_tracked():
    """gitignore は追跡済みのファイルには効かない。既に載っていたら網は意味を持たない。"""
    tracked = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True,
                             check=True).stdout.splitlines()
    suffixes = (".msp", ".msp2", ".dbs", ".lbm", ".lbm2", ".sqlite")
    offenders = [p for p in tracked if p.lower().endswith(suffixes)
                 or p.endswith(store.DIGEST_INDEX_NAME) or "/.library-cache/" in f"/{p}"]
    assert offenders == []
