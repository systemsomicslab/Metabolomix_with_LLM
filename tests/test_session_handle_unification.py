"""対話のセッションとハンドルの一本化（B-1 Task 9 / spec 2026-11-06 §5）。

読み込んだデータセットはすべてハンドルの登録簿（`handles.registry`）に入り、対話のセッションは
直近のハンドル（`current_handle_id`）を指すだけになる。`dataset_handle` を省いた呼び出しも
直近のハンドルのロックの下で動き、閉じた直近のハンドルから別のハンドルへ黙って落ちない。

登録簿はテストごとに conftest の `_fresh_handle_registry` が新しくする。
fixture はすべてこのテストが tmp に作る。
"""
from __future__ import annotations

import asyncio
import json
import threading
import time

import pytest

import server
from metabolomix.core import handles, session_state
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools import handle_tools, mztab_tools
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.reanalysis_fixtures import NEG_FEATURES, write_mztab
from tests.config_fixtures import use_config


@pytest.fixture(autouse=True)
def fresh_session():
    session_state.session = session_state.AnalysisSession()
    yield
    for handle in handles.registry.list():   # ライブラリの接続を閉じる（登録簿は conftest が戻す）
        handles.registry.close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()


def _folder(path):
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
               for i in range(4)]
    spots = [{"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0]},
             {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]}]
    write_msdial_folder(path, spots=spots, samples=samples)
    return path


def _open(**kwargs) -> str:
    out = json.loads(handle_tools.dataset_handle_open(library="none", **kwargs))
    assert out["status"] == "success", out
    return out["dataset_handle"]


def _current():
    return session_state.session.current_handle_id


def _msp(path, name="Choline"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"NAME: {name}\nPRECURSORMZ: 104.107\nIONMODE: Positive\nNum Peaks: 2\n"
                    "60.0808 100\n104.107 50\n", encoding="ascii")
    return path


def _list_rows() -> list[dict]:
    table = json.loads(handle_tools.dataset_handle_list())["handles"]
    header, *lines = table.splitlines()
    columns = header.split("\t")
    return [dict(zip(columns, line.split("\t"))) for line in lines]


def _text(result) -> str:
    """`mcp.call_tool` の戻り（content の列）から本文を取り出す。"""
    return "".join(block.text for block in result)


def test_dataset_load_registers_handle_and_sets_current(tmp_path):
    out = mztab_tools.dataset_load(directory=str(_folder(tmp_path / "a")))
    hid = _current()
    assert hid and hid in out
    assert handles.registry.peek(hid).dataset is session_state.session.dataset


def test_closing_current_handle_does_not_fall_back(tmp_path):  # Review Focus 3
    mztab_tools.dataset_load(directory=str(_folder(tmp_path / "a")))
    first = _current()
    other = _open(directory=str(_folder(tmp_path / "b")))
    assert json.loads(handle_tools.dataset_handle_close(first))["status"] == "success"
    assert session_state.session.dataset is None
    assert '"missing_state"' in mztab_tools.dataset_status()
    # mcp 経由（包みを通る）でも、開いている別のハンドルへは落ちない。
    via_mcp = _text(asyncio.run(server.mcp.call_tool("dataset_status", {})))
    assert '"missing_state"' in via_mcp
    assert handles.registry.peek(other) is not None
    assert _current() is None


def test_library_follows_current_handle(tmp_path):
    sess = session_state.session
    own = sess.library
    assert isinstance(own, session_state.LibraryState)
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    handle = handles.registry.peek(_current())
    assert sess.library is handle.library and sess.library is not own
    replacement = session_state.LibraryState()
    sess.library = replacement          # 直近のハンドルがあれば、そのハンドルへ
    assert handle.library is replacement
    sess.dataset = None                 # 直近のハンドルが無ければ、セッションの LibraryState
    assert sess.library is own
    standalone = session_state.LibraryState()
    sess.library = standalone
    assert sess.library is standalone and handle.library is replacement


def test_handle_limit_refuses_instead_of_evicting(monkeypatch, tmp_path):
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 2})
    paths = [write_mztab(tmp_path / str(i) / "AlignResult-1.mzTab") for i in range(3)]
    mztab_tools.dataset_load(mztab_path=str(paths[0]))
    first = _current()
    mztab_tools.dataset_load(mztab_path=str(paths[1]))
    second = _current()
    err = json.loads(mztab_tools.dataset_load(mztab_path=str(paths[2])))["error"]
    assert err["code"] == "DATASET_HANDLE_LIMIT"
    assert {h.handle_id for h in handles.registry.list()} == {first, second}
    # 失敗した読み込みは直近のデータセットを動かさない。
    assert _current() == second
    assert session_state.session.dataset is handles.registry.peek(second).dataset


def test_interactive_call_waits_for_agent_call_on_same_handle(tmp_path):
    """省略の呼び出しと、同じハンドルを明示した呼び出しが重ならない（同じロックを取る）。"""
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    hid = _current()
    handle = handles.registry.get(hid)
    real_lock = handle.lock
    spans = []

    class SpyLock:
        def __enter__(self):
            real_lock.acquire()
            start = time.perf_counter()
            time.sleep(0.2)
            spans.append(start)
            return self

        def __exit__(self, *exc):
            spans.append(time.perf_counter())
            real_lock.release()

    handle.lock = SpyLock()

    async def main():
        return await asyncio.gather(
            server.mcp.call_tool("dataset_status", {"dataset_handle": hid}),
            server.mcp.call_tool("dataset_status", {}))

    agent, interactive = asyncio.run(main())
    # 2 本とも同じロック（SpyLock）を通った。SpyLock は本物のロックを持ったまま区間を測るので、
    # 重ならないことは作りから決まる——落ちうる主張は「省略の呼び出しもこのロックを取った」だけ。
    assert len(spans) == 4, "省略の呼び出しが直近のハンドルのロックを取っていない"
    assert '"n_features":3' in _text(agent) and '"n_features":3' in _text(interactive)


def test_interactive_call_with_current_handle_runs_off_event_loop(monkeypatch, tmp_path):
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    seen = []
    original = mztab_tools.json_payload

    def spy(*args, **kwargs):
        seen.append(threading.current_thread())
        return original(*args, **kwargs)

    monkeypatch.setattr(mztab_tools, "json_payload", spy)

    async def main():
        seen.append(threading.current_thread())
        await server.mcp.call_tool("dataset_status", {})

    asyncio.run(main())
    assert len(seen) == 2 and seen[0] is not seen[1]


def test_interactive_call_keeps_the_handle_it_locked(tmp_path):
    """呼び出しの途中で直近のデータセットが変わっても、ロックを取ったハンドルのデータを読む
    （ロックを取っていない別のハンドルの状態を触らない）。"""
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    handle = handles.registry.get(_current())
    other = DatasetState()
    other.source_format = "other"
    real_lock = handle.lock
    entered, release = threading.Event(), threading.Event()

    class GateLock:
        def __enter__(self):
            real_lock.acquire()
            entered.set()
            assert release.wait(5)
            return self

        def __exit__(self, *exc):
            real_lock.release()

    handle.lock = GateLock()

    async def main():
        task = asyncio.create_task(server.mcp.call_tool("dataset_status", {}))
        assert await asyncio.to_thread(entered.wait, 5)
        session_state.session.dataset = other   # dataset_load と同じく、イベントループ上で切り替える
        release.set()
        return await task

    out = json.loads(_text(asyncio.run(main())))
    assert out["source_format"] == "mztab" and out["n_features"] == 3
    assert session_state.session.dataset is other


def test_library_loaded_by_interactive_call_is_usable_from_another_thread(tmp_path, monkeypatch):
    """省略の `library_load` はワーカースレッドで開く。その store をイベントループ側の
    ツール（`curation_*` / `verify_peak_annotation`）や別のワーカーが読めること
    （sqlite3 の接続は既定では開いたスレッドでしか使えない）。"""
    msp = _msp(tmp_path / "lab" / "pos.msp")
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    loaded = json.loads(_text(asyncio.run(
        server.mcp.call_tool("library_load", {"file_path": str(msp)}))))
    assert loaded["status"] == "success", loaded
    store = session_state.session.library.store
    assert store is handles.registry.peek(_current()).library.store
    # このスレッド（開いたワーカーではない）から接続を使う。
    records, _ = store.search(name="Choline")
    assert [row["name"] for row in records] == ["Choline"]


def test_handle_list_marks_current(tmp_path):
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    current = _current()
    other = _open(mztab_path=str(write_mztab(tmp_path / "b" / "AlignResult-2.mzTab",
                                             features=NEG_FEATURES, sml_adduct="[M-H]1-")))
    rows = {row["dataset_handle"]: row for row in _list_rows()}
    assert rows[current]["current"] == "True" and rows[other]["current"] == "False"
    assert rows[current]["source_kind"] == "session"


def test_reassigning_the_current_dataset_does_not_open_a_handle():  # R13
    sess = session_state.session
    ds = DatasetState()
    sess.dataset = ds
    hid = sess.current_handle_id
    sess.dataset = ds
    sess.dataset = sess.dataset
    assert sess.current_handle_id == hid
    assert [h.handle_id for h in handles.registry.list()] == [hid]


# ---------- 修正 1 回目（R19・R20・Minor） ----------

def _reads_session_dataset_or_library(fn) -> bool:
    """関数**自身の**ソースが `session.dataset` / `session.library` を読むか
    （`session_state.session.…` と `getattr(session…, "dataset")` を含む。代入は数えない）。"""
    import ast
    import inspect
    import textwrap

    tree = ast.parse(textwrap.dedent(inspect.getsource(inspect.unwrap(fn))))

    def is_session(node):
        return ((isinstance(node, ast.Name) and node.id == "session")
                or (isinstance(node, ast.Attribute) and node.attr == "session"))

    for node in ast.walk(tree):
        if (isinstance(node, ast.Attribute) and node.attr in ("dataset", "library")
                and isinstance(node.ctx, ast.Load) and is_session(node.value)):
            return True
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "getattr" and len(node.args) >= 2
                and is_session(node.args[0]) and isinstance(node.args[1], ast.Constant)
                and node.args[1].value in ("dataset", "library")):
            return True
    return False


def test_tools_reading_the_session_dataset_take_the_current_handle_lock():  # R20
    """`session.dataset` / `session.library` を読むツールは、`dataset_handle` を持つか
    `handle_dispatch.SESSION_READERS` に載り、`install` に包まれている（直近のハンドルの
    ロックの下で動く）。

    見るのはツールの関数**自身の**ソースだけ。ヘルパー経由で読むもの（`pathway_project` の
    `_build_dataset`、`save_figure` / `plot_pca_loadings` の候補集め、`verify_peak_annotation` の
    `msms_evidence`）はこの検査では見つからない——推移的な読み取りは範囲外で、
    SESSION_READERS に手で載せてある。
    """
    import inspect

    from metabolomix.core import handle_dispatch

    tools = server.mcp._tool_manager._tools
    unguarded = []
    for name, tool in tools.items():
        fn = inspect.unwrap(tool.fn)
        if (_reads_session_dataset_or_library(fn)
                and "dataset_handle" not in inspect.signature(fn).parameters
                and name not in handle_dispatch.SESSION_READERS):
            unguarded.append(name)
    assert unguarded == []
    assert handle_dispatch.SESSION_READERS <= set(tools)
    assert all(tools[name].is_async for name in handle_dispatch.SESSION_READERS)


def test_library_load_does_not_close_the_store_under_curation_review(tmp_path, monkeypatch):  # R20
    """省略の `library_load`（ワーカー・直近のハンドルのロック）が古い store を閉じても、
    `curation_review` がその閉じた store を使わない（レビューも直近のハンドルのロックを待つ）。

    順序は Event で決める: library_load がロックを取って止まる → curation_review を出す。
    修正前はレビューがイベントループで古い store を掴み、library_load を進めて閉じさせてから
    使う（`Cannot operate on a closed database`）。修正後はレビューがロックを待ち、
    library_load が新しい store に替えた後で掴む。
    """
    from metabolomix.curation import flags as flag_log
    from metabolomix.curation import review as review_module
    from metabolomix.library import tools as library_tools
    from metabolomix.tools import curation_tools

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    old = json.loads(library_tools.library_load(file_path=str(_msp(tmp_path / "a" / "old.msp"))))
    assert old["status"] == "success", old
    new_msp = _msp(tmp_path / "b" / "new.msp", name="Betaine")
    handle = handles.registry.get(_current())
    real_lock = handle.lock
    load_entered, load_go, load_released, review_waiting = (threading.Event() for _ in range(4))
    entries, exits = [], []

    class GateLock:
        def __enter__(self):
            entries.append(1)
            if len(entries) == 1:          # library_load: ロックを取ったまま止まる
                real_lock.acquire()
                load_entered.set()
                assert load_go.wait(5)
            else:                          # curation_review（修正後）と後片付け
                review_waiting.set()
                real_lock.acquire()
            return self

        def __exit__(self, *exc):
            exits.append(1)
            real_lock.release()
            if len(exits) == 1:
                load_released.set()

    handle.lock = GateLock()
    used = []

    def thresholds(value):
        # store を掴んだ後。修正前（イベントループ）なら、ここで library_load を進めて閉じさせる。
        load_go.set()
        assert load_released.wait(5)
        return {}

    def run_review(arf2_path, spots, *, store, **kwargs):
        try:
            records, _ = store.search(name="")
            used.append(("ok", [row["name"] for row in records]))
        except Exception as exc:  # noqa: BLE001 - 閉じた接続の例外をそのまま記録する
            used.append(("error", repr(exc)))
        raise flag_log.FlagFileError("flags.jsonl", 1, "テストの打ち切り")

    monkeypatch.setattr(curation_tools, "resolve_arf2_file_path", lambda file_path=None: "x.arf2")
    monkeypatch.setattr(curation_tools.judge, "resolve_thresholds", thresholds)
    monkeypatch.setattr(curation_tools, "load_catalog", lambda path: [])
    monkeypatch.setattr(curation_tools.evidence, "select_spots", lambda *a, **k: [{"spot": 0}])
    monkeypatch.setattr(review_module, "run_review", run_review)

    async def main():
        load = asyncio.create_task(server.mcp.call_tool("library_load",
                                                        {"file_path": str(new_msp)}))
        assert await asyncio.to_thread(load_entered.wait, 5)
        review = asyncio.create_task(server.mcp.call_tool("curation_review", {}))
        waiting = asyncio.ensure_future(asyncio.to_thread(review_waiting.wait, 5))
        await asyncio.wait({review, waiting}, return_when=asyncio.FIRST_COMPLETED)
        load_go.set()
        return await asyncio.gather(load, review)

    loaded, _ = asyncio.run(main())
    assert json.loads(_text(loaded))["status"] == "success"
    assert used == [("ok", ["Betaine"])]


def test_library_loaded_before_dataset_load_moves_into_the_new_handle(tmp_path, monkeypatch):  # R19
    from metabolomix.library import tools as library_tools

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    loaded = json.loads(library_tools.library_load(file_path=str(_msp(tmp_path / "lab" / "pos.msp"))))
    assert loaded["status"] == "success", loaded
    store = session_state.session.library.store
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    assert session_state.session.library.store is store
    assert handles.registry.peek(_current()).library.store is store
    assert session_state.session._library.store is None   # 写したのではなく移した
    records, _ = session_state.session.library.store.search(name="Choline")
    assert len(records) == 1


def test_refused_dataset_assignment_keeps_the_library_in_the_session(tmp_path, monkeypatch):  # R19
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.library import tools as library_tools

    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    _open(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))  # 上限を埋める（直近にならない）
    library_tools.library_load(file_path=str(_msp(tmp_path / "lab" / "pos.msp")))
    store = session_state.session.library.store
    assert store is not None
    with pytest.raises(DomainError):
        session_state.session.dataset = DatasetState()
    assert _current() is None
    assert session_state.session.library.store is store


def test_closing_the_current_handle_does_not_fall_back_to_an_older_library(tmp_path,
                                                                            monkeypatch):  # R19
    from metabolomix.library import tools as library_tools

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    library_tools.library_load(file_path=str(_msp(tmp_path / "a" / "old.msp")))
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    library_tools.library_load(file_path=str(_msp(tmp_path / "b" / "new.msp", name="Betaine")))
    handle_tools.dataset_handle_close(_current())
    assert session_state.session.library.store is None
    out = json.loads(library_tools.library_match_feature(precursor_mz=104.107))
    assert out["error"]["code"] == "missing_state"


def test_dataset_load_handle_resolves_alignment_files_when_named(tmp_path):  # Minor 1
    """`dataset_load` の要約に出るハンドル ID を明示しても、省略と同じ `.arf2` / `.dcl` を使う。"""
    from metabolomix.dcl import tools as dcl_tools
    from tests.dcl_fixture import build_dcl_bytes

    stem = "AlignmentResult_2026_02_01_00_00_00"
    folder = _folder(tmp_path / "gui")
    (folder / f"{stem}.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": 760.5851, "rt": 12.0, "spectrum": [(184.0733, 100.0)]},
        {"precursor_mz": 500.1, "rt": 12.0, "spectrum": [(100.0, 10.0)]},
    ]))
    mztab_tools.dataset_load(directory=str(folder))
    hid = _current()
    handle = handles.registry.peek(hid)
    assert handle.arf2_path == str(folder / f"{stem}.arf2")
    assert handle.dcl_path == str(folder / f"{stem}.dcl")
    out = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=0, dataset_handle=hid))
    assert out["status"] == "ok", out
    assert out["precursor_mz"] == 760.5851


def test_call_pinned_to_a_handle_closed_before_the_lock_does_nothing(tmp_path,
                                                                      monkeypatch):  # Minor 2
    """dispatch が直近として掴んだハンドルが、ワーカーがロックを取る前に閉じられた →
    `missing_state`。直近になった別のハンドルにもセッションのライブラリにも触らない。"""
    from metabolomix.core import handle_dispatch
    from metabolomix.library import tools as library_tools

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    first = handles.registry.peek(_current())
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "b" / "AlignResult-1.mzTab")))
    second = _current()
    handles.registry.close(first.handle_id)

    status = json.loads(handle_dispatch._call_pinned(mztab_tools.dataset_status, {}, first))
    assert status["error"]["code"] == "missing_state"
    msp = _msp(tmp_path / "lab" / "pos.msp")
    loaded = json.loads(handle_dispatch._call_pinned(library_tools.library_load,
                                                     {"file_path": str(msp)}, first))
    assert loaded["error"]["code"] == "missing_state"
    assert session_state.session._library.store is None
    assert handles.registry.peek(second).library.store is None
    assert _current() == second
    with session_state.pinned_handle(first.handle_id):
        assert session_state.session.dataset is None
        assert session_state.session.library.store is None
        assert session_state.session.library is not session_state.session._library


def test_dataset_load_refuses_at_the_limit_before_loading(monkeypatch, tmp_path):  # Minor 4
    from metabolomix.dataset import opening

    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    calls = []
    real = opening.open_dataset

    def spy(**kwargs):
        calls.append(kwargs)
        return real(**kwargs)

    monkeypatch.setattr(opening, "open_dataset", spy)
    path = str(write_mztab(tmp_path / "b" / "AlignResult-1.mzTab"))
    assert json.loads(mztab_tools.dataset_load(mztab_path=path))["error"]["code"] \
        == "DATASET_HANDLE_LIMIT"
    out = json.loads(handle_tools.dataset_handle_open(mztab_path=path, library="none"))
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert calls == []
