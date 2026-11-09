"""解析セッションにパーサ別のスロットが無いことの回帰テスト。

背景: 以前は `session.arf` / `session.arf2` / `session.pai2` / `session.eic` というパーサ別の
スロットがあり、パーサ間の混線（`pai2_parser` が ARF の前処理行列を消すなど）を防ぐために
分離していた。B-2b で読み込みの入口を `load_dataset` に一本化し、データセットはハンドルの
登録簿にだけ住むことになった。パーサ別のスロットへ戻ると、同じ混線がまた起きる。
"""
import pytest

from metabolomix.core import session_state


@pytest.mark.parametrize("slot", ["arf", "arf2", "pai2", "eic", "last_loadings_plot"])
def test_session_has_no_parser_slots(slot):
    assert not hasattr(session_state.AnalysisSession(), slot)


@pytest.mark.parametrize("name", ["ArfState", "Arf2State", "Pai2State", "EicState",
                                  "_build_sample_meta", "_new_dataset_id"])
def test_session_state_module_has_no_parser_state_classes(name):
    assert not hasattr(session_state, name)


def test_session_keeps_the_handle_registry_slots():
    s = session_state.AnalysisSession()
    assert s.dataset is None
    assert s.library is not None
    assert s.curation is not None
    assert s.pathway is not None
