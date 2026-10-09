"""全テスト共通の隔離。

外部資産の場所（`MSDIAL_EXE` / `MSDIAL_LBM` / `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）は、
利用者が OS や `.mcp.json` に常設するか、リポジトリ直下の `lipidmix.local.toml` に書く
前提の設定なので、テストを走らせる機械にも置かれている。残したままだと
`resolve_library_path()` がテストの tmp ではなく実ライブラリを掴み、結果が機械ごとに
変わる（しかも実ライブラリは外部流出禁止の資産）。環境変数は消し、設定ファイルは
存在しないパスへ向ける。使うテストは monkeypatch.setenv で明示的に置く
（設定ファイルは tmp に書いて `LIPIDMIX_CONFIG` で指す）。
"""
import sys
from pathlib import Path

import pytest

from metabolomix.core.user_config import CONFIG_ENV, SETTINGS

#: 作られることの無いパス。autouse で tmp_path を要求すると全テストに tmp ディレクトリが
#: できるので、存在しない固定パスで済ませる。
_NO_CONFIG = Path(__file__).resolve().parent / "_no_such_dir" / "lipidmix.local.toml"


@pytest.fixture(autouse=True)
def _isolate_external_asset_settings(monkeypatch):
    for name in SETTINGS.values():
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(CONFIG_ENV, str(_NO_CONFIG))


@pytest.fixture(autouse=True)
def _stop_curation_submit_server():
    """curation_review / curation_suggest はビューアの受け口（daemon スレッドの HTTP サーバ）を
    立てる。テストをまたいで登録と待受けを残さない。未 import なら何もしない（重い import を避ける）。"""
    yield
    module = sys.modules.get("metabolomix.curation.submit_server")
    if module is not None:
        module.reset()
