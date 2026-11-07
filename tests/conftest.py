"""全テスト共通の隔離。

外部資産の場所（Console の実行体・脂質ライブラリ・研究室の参照ライブラリ）は、
リポジトリ直下の `lipidmix.local.toml` に書く前提の設定なので、テストを走らせる機械にも
置かれている。そのまま読むと `resolve_library_path()` がテストの tmp ではなく実ライブラリを
掴み、結果が機械ごとに変わる（しかも実ライブラリは外部流出禁止の資産）。設定ファイルは
存在しないパスへ向ける。使うテストは tmp に書いて `LIPIDMIX_CONFIG` で指す
（`tests/config_fixtures.py` の `use_config`）。資産の場所の環境変数は読まない
（spec 2026-11-06 §7.4）ので消す必要もない。
"""
import sys
from pathlib import Path

import pytest

from metabolomix.core.user_config import CONFIG_ENV

#: 作られることの無いパス。autouse で tmp_path を要求すると全テストに tmp ディレクトリが
#: できるので、存在しない固定パスで済ませる。
_NO_CONFIG = Path(__file__).resolve().parent / "_no_such_dir" / "lipidmix.local.toml"


@pytest.fixture(autouse=True)
def _isolate_external_asset_settings(monkeypatch):
    monkeypatch.setenv(CONFIG_ENV, str(_NO_CONFIG))


@pytest.fixture(autouse=True)
def _fresh_handle_registry(monkeypatch):
    """テストごとに新しいハンドルの登録簿にする。

    `session.dataset = ds` も登録簿にハンドルを足す（spec 2026-11-06 §5）。登録簿を
    テストをまたいで持ち越すと、前のテストのハンドルが上限（既定 8）を食い、
    `session.dataset` が前のテストのデータセットを指し続ける。
    """
    from metabolomix.core import handles

    monkeypatch.setattr(handles, "registry", handles.HandleRegistry())


@pytest.fixture(autouse=True)
def _stop_curation_submit_server():
    """curation_review / curation_suggest はビューアの受け口（daemon スレッドの HTTP サーバ）を
    立てる。テストをまたいで登録と待受けを残さない。未 import なら何もしない（重い import を避ける）。"""
    yield
    module = sys.modules.get("metabolomix.curation.submit_server")
    if module is not None:
        module.reset()


@pytest.fixture(autouse=True)
def _fresh_file_digest():
    """sha256 の記憶（パス・大きさ・更新時刻がキー）をテストごとに空にする。

    同じ tmp パスへ同じ大きさで素早く書き直すテストが、更新時刻の刻みの粗い環境で
    前の中身の値を拾わないようにする。"""
    from metabolomix.core import file_digest

    file_digest.clear()
    yield
