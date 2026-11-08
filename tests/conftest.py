"""全テスト共通の隔離。

外部資産の場所（Console の実行体・脂質ライブラリ・研究室の参照ライブラリ）は、
リポジトリ直下の `lipidmix.local.toml` に書く前提の設定なので、テストを走らせる機械にも
置かれている。そのまま読むと `resolve_library_path()` がテストの tmp ではなく実ライブラリを
掴み、結果が機械ごとに変わる（しかも実ライブラリは外部流出禁止の資産）。設定ファイルは
存在しないパスへ向ける。使うテストは tmp に書いて `LIPIDMIX_CONFIG` で指す
（`tests/config_fixtures.py` の `use_config`）。資産と置き場所の環境変数は読まない
（spec 2026-11-06 §7.4）ので消す必要もない。

向け先はテストごとの tmp の中にする。参照ライブラリのキャッシュの既定は「設定ファイルの
フォルダの data/.library-cache」なので、リポジトリの中を向けるとテストがリポジトリに
キャッシュを作る。`use_config` も `tmp_path` に書くので、そのテストのキャッシュも tmp に入る。
"""
import itertools
import sys
from pathlib import Path

import pytest

from metabolomix.core.user_config import CONFIG_ENV

def pytest_addoption(parser):
    """本物の規則プラグインとの同値確認（tests/test_quant_parity.py）だけが読む。"""
    parser.addoption("--quant-rules", action="store", default=None,
                     help="規則パッケージ msdial_rules のフォルダ。渡したときだけ同値確認を走らせる")


_REAL_RULES_SUITES = ("tests.test_quant_normalize", "tests.test_quant_merge",
                      "tests.test_quant_select", "tests.test_quant_resume")


@pytest.fixture(autouse=True)
def _real_rules_for_quant_suites(request, monkeypatch):
    """`--quant-rules` 付きなら、定量 3 スイート（と再開）の `write_fake_rules` は本物の規則の
    フォルダを返す（同じテストを本物の規則でも走らせる）。plugin / parity のテストは対象外。"""
    folder = request.config.getoption("--quant-rules")
    if folder and request.module.__name__ in _REAL_RULES_SUITES:
        from tests import quant_fake_rules

        monkeypatch.setattr(quant_fake_rules, "REAL_RULES_FOLDER", Path(folder))


_NO_CONFIG_SEQ = itertools.count()


@pytest.fixture(autouse=True)
def _isolate_external_asset_settings(monkeypatch, tmp_path_factory):
    """設定ファイルを、テストごとに別の・存在しないパスへ向ける（何も作らない）。

    autouse で `tmp_path` を要求すると全テストに tmp ディレクトリができるので、基底の tmp の
    下に番号違いのパスを組むだけにする（キャッシュを作るテストだけがフォルダを作る）。
    """
    base = tmp_path_factory.getbasetemp() / "no-config" / str(next(_NO_CONFIG_SEQ))
    monkeypatch.setenv(CONFIG_ENV, str(base / "lipidmix.local.toml"))


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
