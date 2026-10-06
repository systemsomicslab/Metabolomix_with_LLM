# tests/test_process_control.py
"""Windowsのプロセス所有権・identity・OS排他の検証。

モックした偽プロセスではなく、実際の無害なPython子プロセス・孫プロセスを起動して
「本当に終了したか」「本当に排他できているか」を確かめる。Job Objectもbyte-range
lockもOSの状態であり、Python側の記録を読んでも検証にならないため。
"""
import contextlib
import ctypes
import ctypes.wintypes as wintypes
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from lipidmix.core.atomic_io import DomainError
from lipidmix.core.process_control import (
    detach_breakaway_mode,
    file_lock,
    launch_detached,
    process_identity,
    same_process,
    start_owned_process,
)
import lipidmix.core.process_control as pc

_REPO_ROOT = str(Path(__file__).resolve().parents[1])

#: 検証中はずっと生きていてほしいだけの無害な子。テストは必ず明示的に終了させる。
_SLEEP = [sys.executable, "-c", "import time; time.sleep(120)"]

windows_only = pytest.mark.skipif(os.name != "nt",
                                  reason="Windows Job Objectの検証")


def _wait_for_file(path: Path, timeout: float = 20.0) -> str:
    """子プロセスが書くファイルが空でなくなるまで待って中身を返す。"""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            text = ""
        if text.strip():
            return text
        time.sleep(0.02)
    raise AssertionError(f"子プロセスが {path} を書かなかった")


def _wait_until_gone(identity: dict, timeout: float = 10.0) -> bool:
    """identityのプロセスが消えるまで待つ。消えたらTrue。"""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not same_process(identity):
            return True
        time.sleep(0.02)
    return not same_process(identity)


def _kill_pid_if_same(identity: dict) -> None:
    """テスト後始末。pid再利用を踏まないようidentity一致時だけ終了させる。"""
    if not identity or not same_process(identity):
        return
    try:
        os.kill(int(identity["pid"]), signal.SIGTERM)
    except OSError:
        pass


def _handle_count() -> int:
    """現在のプロセスが開いているカーネルハンドル数。"""
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.GetProcessHandleCount.argtypes = [wintypes.HANDLE,
                                               ctypes.POINTER(wintypes.DWORD)]
    kernel32.GetProcessHandleCount.restype = wintypes.BOOL
    kernel32.GetCurrentProcess.argtypes = []
    kernel32.GetCurrentProcess.restype = wintypes.HANDLE
    count = wintypes.DWORD()
    assert kernel32.GetProcessHandleCount(kernel32.GetCurrentProcess(),
                                          ctypes.byref(count))
    return int(count.value)


def _write_script(tmp_path: Path, name: str, source: str) -> Path:
    """子プロセス用スクリプトを書く。

    `-c` に改行入りソースを渡さないのは、コマンドライン経由で崩れる余地を
    残さないため。
    """
    path = tmp_path / name
    path.write_text(source, encoding="utf-8")
    return path


# ---------- 代表テスト（Task 2 brief） ----------

@windows_only
def test_termination_reaps_owned_process(tmp_path):
    proc = start_owned_process([sys.executable, "-c", "import time; time.sleep(30)"],
                               cwd=tmp_path, log_path=tmp_path / "child.log")
    identity = proc.identity
    try:
        assert same_process(identity)
        proc.terminate_tree()
        assert proc.wait(timeout=5) != 0
        assert not same_process(identity)
    finally:
        proc.close()


# ---------- identity: pid再利用の見分け ----------

@windows_only
def test_same_pid_with_different_creation_time_is_not_same_process(tmp_path):
    """pidが同じでも生成時刻が違えばFalse。pid再利用で別人を殺さないための要。"""
    proc = start_owned_process(_SLEEP, cwd=tmp_path, log_path=tmp_path / "c.log")
    try:
        real = proc.identity
        assert same_process(real)
        recycled = {"pid": real["pid"], "creation_time": real["creation_time"] + 1}
        assert not same_process(recycled)
    finally:
        proc.terminate_tree()
        proc.close()


@windows_only
def test_process_identity_reports_creation_time_and_pid(tmp_path):
    proc = start_owned_process(_SLEEP, cwd=tmp_path, log_path=tmp_path / "c.log")
    try:
        found = process_identity(proc.pid)
        assert found == proc.identity
        assert found["pid"] == proc.pid
        assert isinstance(found["creation_time"], int)
        assert found["creation_time"] > 0
    finally:
        proc.terminate_tree()
        proc.close()


@windows_only
def test_process_identity_is_none_for_exited_process(tmp_path):
    """終了直後はハンドルが残りOpenProcessが成功しうる。生存判定に混ぜない。"""
    proc = start_owned_process([sys.executable, "-c", "pass"],
                               cwd=tmp_path, log_path=tmp_path / "c.log")
    identity = proc.identity
    try:
        assert proc.wait(timeout=20) == 0
        # OwnedProcessがまだハンドルを保持している状態で問い合わせる。
        assert process_identity(proc.pid) is None
        assert not same_process(identity)
    finally:
        proc.close()


@windows_only
@pytest.mark.parametrize("expected", [
    None, {}, {"pid": 1234}, {"creation_time": 100},
    {"pid": None, "creation_time": 100}, {"pid": 0, "creation_time": 100},
    {"pid": -1, "creation_time": 100}, "not-a-dict",
])
def test_same_process_rejects_malformed_identity(expected):
    assert same_process(expected) is False


# ---------- OwnedProcess: poll / wait ----------

@windows_only
def test_poll_returns_none_while_running_then_exit_code(tmp_path):
    proc = start_owned_process([sys.executable, "-c", "import time; time.sleep(0.5)"],
                               cwd=tmp_path, log_path=tmp_path / "c.log")
    try:
        assert proc.poll() is None
        assert proc.wait(timeout=20) == 0
        assert proc.poll() == 0
    finally:
        proc.close()


@windows_only
def test_wait_timeout_raises_domain_error(tmp_path):
    proc = start_owned_process(_SLEEP, cwd=tmp_path, log_path=tmp_path / "c.log")
    try:
        with pytest.raises(DomainError) as excinfo:
            proc.wait(timeout=0.2)
        assert excinfo.value.code == "PROCESS_WAIT_TIMEOUT"
    finally:
        proc.terminate_tree()
        proc.close()


@windows_only
def test_operations_after_close_raise(tmp_path):
    proc = start_owned_process(_SLEEP, cwd=tmp_path, log_path=tmp_path / "c.log")
    identity = proc.identity
    proc.terminate_tree()
    proc.close()
    proc.close()  # 冪等
    with pytest.raises(DomainError) as excinfo:
        proc.poll()
    assert excinfo.value.code == "PROCESS_HANDLE_CLOSED"
    assert not same_process(identity)


@windows_only
def test_child_stdout_goes_to_log_path(tmp_path):
    log_path = tmp_path / "nested" / "child.log"
    proc = start_owned_process(
        [sys.executable, "-c", "print('hello-from-child')"],
        cwd=tmp_path, log_path=log_path)
    try:
        assert proc.wait(timeout=20) == 0
    finally:
        proc.close()
    assert "hello-from-child" in log_path.read_text(encoding="utf-8")


@windows_only
@pytest.mark.parametrize("command, cwd_exists, code", [
    ([], True, "COMMAND_EMPTY"),
    (["python"], False, "CWD_NOT_FOUND"),
])
def test_launch_argument_validation(tmp_path, command, cwd_exists, code):
    cwd = tmp_path if cwd_exists else tmp_path / "missing"
    with pytest.raises(DomainError) as excinfo:
        start_owned_process(command, cwd=cwd, log_path=tmp_path / "c.log")
    assert excinfo.value.code == code


# ---------- 孫プロセスまで終了させる ----------

_SPAWNER_SRC = """\
import subprocess
import sys
import time
from pathlib import Path

grandchild = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
Path(sys.argv[1]).write_text(str(grandchild.pid), encoding="utf-8")
time.sleep(300)
"""


@windows_only
def test_terminate_tree_also_reaps_grandchild(tmp_path):
    """孫まで終了させる。親のpidだけをkillすると孫が残り、出力を掴んだままになる。"""
    script = _write_script(tmp_path, "spawner.py", _SPAWNER_SRC)
    pid_file = tmp_path / "grandchild.pid"
    proc = start_owned_process([sys.executable, str(script), str(pid_file)],
                               cwd=tmp_path, log_path=tmp_path / "c.log")
    grandchild_identity = None
    try:
        grandchild_pid = int(_wait_for_file(pid_file).strip())
        grandchild_identity = process_identity(grandchild_pid)
        assert grandchild_identity is not None, "孫プロセスが起動していない"

        proc.terminate_tree()

        assert proc.wait(timeout=5) != 0
        assert not same_process(proc.identity)
        assert not same_process(grandchild_identity), "孫プロセスが残っている"
    finally:
        _kill_pid_if_same(grandchild_identity)
        proc.close()


@windows_only
def test_close_without_terminate_still_reaps_tree(tmp_path):
    """KILL_ON_JOB_CLOSEの確認。close()だけでも孤児を残さない。"""
    script = _write_script(tmp_path, "spawner.py", _SPAWNER_SRC)
    pid_file = tmp_path / "grandchild.pid"
    proc = start_owned_process([sys.executable, str(script), str(pid_file)],
                               cwd=tmp_path, log_path=tmp_path / "c.log")
    grandchild_identity = None
    try:
        grandchild_pid = int(_wait_for_file(pid_file).strip())
        grandchild_identity = process_identity(grandchild_pid)
        assert grandchild_identity is not None
        child_identity = proc.identity

        proc.close()

        assert _wait_until_gone(child_identity), "子がJobハンドル解放で終了していない"
        assert _wait_until_gone(grandchild_identity), "孫が残っている"
    finally:
        _kill_pid_if_same(grandchild_identity)
        proc.close()


# ---------- Assign失敗: resumeせずに停止中の子を始末する ----------

def _spy_create_process(monkeypatch, seen: dict):
    """CreateProcessWの結果から実pidを覗くだけのラッパーを挟む。"""
    real_create = pc._CreateProcessW

    def spy(*args):
        ok = real_create(*args)
        if ok:
            seen["pid"] = int(args[9]._obj.dwProcessId)
        return ok

    monkeypatch.setattr(pc, "_CreateProcessW", spy)


@windows_only
def test_assign_failure_terminates_child_without_resuming(tmp_path, monkeypatch):
    """Assignに失敗したら停止中の子を始末する。resumeすると所有できない子が野に出る。"""
    seen: dict = {}
    resumed: list = []
    _spy_create_process(monkeypatch, seen)
    monkeypatch.setattr(pc, "_AssignProcessToJobObject", lambda job, proc: 0)
    monkeypatch.setattr(pc, "_ResumeThread",
                        lambda handle: (resumed.append(handle), 0)[1])

    with pytest.raises(DomainError) as excinfo:
        start_owned_process(_SLEEP, cwd=tmp_path, log_path=tmp_path / "c.log")

    assert excinfo.value.code == "PROCESS_ASSIGN_FAILED"
    assert resumed == [], "Assign失敗後にResumeThreadを呼んでいる"
    assert "pid" in seen, "子プロセスが起動していない（前提が崩れている）"
    # 停止中の子が本当に消えたことを、pidの解決不能で確かめる。
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline and process_identity(seen["pid"]) is not None:
        time.sleep(0.02)
    assert process_identity(seen["pid"]) is None, "停止中の子が残っている"


@windows_only
def test_resume_failure_terminates_child(tmp_path, monkeypatch):
    seen: dict = {}
    _spy_create_process(monkeypatch, seen)
    monkeypatch.setattr(pc, "_ResumeThread", lambda handle: 0xFFFFFFFF)

    with pytest.raises(DomainError) as excinfo:
        start_owned_process(_SLEEP, cwd=tmp_path, log_path=tmp_path / "c.log")

    assert excinfo.value.code == "PROCESS_RESUME_FAILED"
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline and process_identity(seen["pid"]) is not None:
        time.sleep(0.02)
    assert process_identity(seen["pid"]) is None


# ---------- ハンドル解放 ----------

@windows_only
def test_no_handle_leak_across_start_terminate_close_cycles(tmp_path):
    """成功経路をくり返してもハンドルが積み上がらない。"""
    def cycle(index: int) -> None:
        proc = start_owned_process(_SLEEP, cwd=tmp_path,
                                   log_path=tmp_path / f"c{index}.log")
        try:
            proc.terminate_tree()
        finally:
            proc.close()

    cycle(0)  # 初回はctypes・ログ関連の一時的な確保を含むので基準から外す
    baseline = _handle_count()
    for index in range(1, 7):
        cycle(index)
    assert _handle_count() <= baseline + 2


@windows_only
def test_no_handle_leak_on_assign_failure_path(tmp_path, monkeypatch):
    """失敗経路でもJob・プロセス・スレッド・ログの各ハンドルを解放する。"""
    monkeypatch.setattr(pc, "_AssignProcessToJobObject", lambda job, proc: 0)

    def cycle(index: int) -> None:
        with pytest.raises(DomainError):
            start_owned_process(_SLEEP, cwd=tmp_path,
                                log_path=tmp_path / f"f{index}.log")

    cycle(0)
    baseline = _handle_count()
    for index in range(1, 6):
        cycle(index)
    assert _handle_count() <= baseline + 2


# ---------- OSロック ----------

_LOCKER_SRC = """\
import sys
import time
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from lipidmix.core.atomic_io import DomainError
from lipidmix.core.process_control import file_lock

lock_path = Path(sys.argv[2])
result_path = Path(sys.argv[3])
timeout = float(sys.argv[4])
hold_s = float(sys.argv[5])
try:
    with file_lock(lock_path, timeout=timeout, poll_interval=0.02):
        result_path.write_text("acquired", encoding="utf-8")
        time.sleep(hold_s)
except DomainError as exc:
    result_path.write_text(exc.code, encoding="utf-8")
"""


def _run_locker(tmp_path: Path, lock_path: Path, result_name: str,
                timeout: float, hold_s: float = 0.0) -> str:
    script = _write_script(tmp_path, "locker.py", _LOCKER_SRC)
    result_path = tmp_path / result_name
    subprocess.run([sys.executable, str(script), _REPO_ROOT, str(lock_path),
                    str(result_path), str(timeout), str(hold_s)],
                   cwd=str(tmp_path), stdin=subprocess.DEVNULL,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   timeout=60, check=True)
    return result_path.read_text(encoding="utf-8").strip()


def test_file_lock_excludes_a_concurrent_process(tmp_path):
    """別プロセスは、こちらが握っている間ロックを取れない。"""
    lock_path = tmp_path / "index.lock"
    with file_lock(lock_path):
        assert _run_locker(tmp_path, lock_path, "blocked.txt",
                           timeout=0.5) == "LOCK_TIMEOUT"
    # 解放後は同じ別プロセスが取れる。
    assert _run_locker(tmp_path, lock_path, "free.txt", timeout=10) == "acquired"


def test_lock_file_left_on_disk_does_not_mean_locked(tmp_path):
    """ロックの正体はOSの状態。ファイルが残っていること＝ロック中ではない。"""
    lock_path = tmp_path / "index.lock"
    with file_lock(lock_path):
        pass
    assert lock_path.exists()
    assert _run_locker(tmp_path, lock_path, "after.txt", timeout=10) == "acquired"


def test_deleting_owner_record_does_not_release_the_lock(tmp_path):
    """owner記録を消してもロックは解けない。記録の有無で解放を判断しない。"""
    lock_path = tmp_path / "index.lock"
    owner_record = tmp_path / "index.owner"
    with file_lock(lock_path):
        owner_record.write_text("worker-1", encoding="utf-8")
        owner_record.unlink()
        assert _run_locker(tmp_path, lock_path, "still.txt",
                           timeout=0.5) == "LOCK_TIMEOUT"


def test_file_lock_releases_on_exception(tmp_path):
    lock_path = tmp_path / "index.lock"
    with pytest.raises(RuntimeError):
        with file_lock(lock_path):
            raise RuntimeError("boom")
    assert _run_locker(tmp_path, lock_path, "after-error.txt",
                       timeout=10) == "acquired"


def test_file_lock_is_exclusive_within_the_same_process(tmp_path):
    """同一プロセスの別ハンドルからも取れない（OSロックであることの確認）。"""
    lock_path = tmp_path / "index.lock"
    with file_lock(lock_path):
        with pytest.raises(DomainError) as excinfo:
            with file_lock(lock_path, timeout=0.3, poll_interval=0.02):
                pass
        assert excinfo.value.code == "LOCK_TIMEOUT"


def test_file_lock_creates_missing_parent_directory(tmp_path):
    lock_path = tmp_path / "nested" / "deeper" / "index.lock"
    with file_lock(lock_path) as held:
        assert held == lock_path
        assert lock_path.exists()


# ---------- 切り離し起動 ----------

@pytest.mark.parametrize("in_job, flags, expected", [
    (False, 0, "not_in_job"),
    (False, pc.JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE, "not_in_job"),
    (True, 0, "inherited_job"),
    (True, pc.JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE, "unsupported"),
    (True, pc.JOB_OBJECT_LIMIT_BREAKAWAY_OK, "breakaway_ok"),
    (True, pc.JOB_OBJECT_LIMIT_SILENT_BREAKAWAY_OK, "silent_breakaway"),
    (True, pc.JOB_OBJECT_LIMIT_BREAKAWAY_OK
     | pc.JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE, "breakaway_ok"),
    (True, pc.JOB_OBJECT_LIMIT_SILENT_BREAKAWAY_OK
     | pc.JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE, "silent_breakaway"),
])
def test_breakaway_mode_decision_table(in_job, flags, expected):
    assert pc._breakaway_mode(in_job, flags) == expected


# --- テスト自前のJob Object操作（検証対象のコードを流用しない） ---
#
# Windowsには「親が死んだら子も死ぬ」という仕組みが無い。孤児になった子は
# 既定で生き続けるので、親プロセスを終わらせて生存を見ても切り離しの検証には
# ならない。実際に道連れにするのはJob Objectの KILL_ON_JOB_CLOSE なので、
# テスト側でそのJobを作り、その中でlaunch_detachedを走らせる。

_JOB_OBJECT_EXTENDED_LIMIT_INFORMATION_CLASS = 9
_TEST_JOB_LIMIT_BREAKAWAY_OK = 0x00000800
_TEST_JOB_LIMIT_KILL_ON_JOB_CLOSE = 0x00002000
_PROCESS_TERMINATE = 0x0001
_PROCESS_SET_QUOTA = 0x0100
_PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
_SYNCHRONIZE = 0x00100000


class _TEST_JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [
        ("PerProcessUserTimeLimit", ctypes.c_int64),
        ("PerJobUserTimeLimit", ctypes.c_int64),
        ("LimitFlags", wintypes.DWORD),
        ("MinimumWorkingSetSize", ctypes.c_size_t),
        ("MaximumWorkingSetSize", ctypes.c_size_t),
        ("ActiveProcessLimit", wintypes.DWORD),
        ("Affinity", ctypes.POINTER(ctypes.c_ulong)),
        ("PriorityClass", wintypes.DWORD),
        ("SchedulingClass", wintypes.DWORD),
    ]


class _TEST_IO_COUNTERS(ctypes.Structure):
    _fields_ = [(name, ctypes.c_uint64) for name in (
        "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
        "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class _TEST_JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [
        ("BasicLimitInformation", _TEST_JOBOBJECT_BASIC_LIMIT_INFORMATION),
        ("IoInfo", _TEST_IO_COUNTERS),
        ("ProcessMemoryLimit", ctypes.c_size_t),
        ("JobMemoryLimit", ctypes.c_size_t),
        ("PeakProcessMemoryUsed", ctypes.c_size_t),
        ("PeakJobMemoryUsed", ctypes.c_size_t),
    ]


def _job_kernel32():
    """Job Object操作用のkernel32。全関数にargtypes/restypeを宣言する。

    宣言を省くとHANDLEが既定のint32へ切り詰められ、別のJobやプロセスを
    指したまま「検証できたつもり」になる。
    """
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    k.CreateJobObjectW.restype = wintypes.HANDLE
    k.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                          ctypes.c_void_p, wintypes.DWORD]
    k.SetInformationJobObject.restype = wintypes.BOOL
    k.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    k.AssignProcessToJobObject.restype = wintypes.BOOL
    k.IsProcessInJob.argtypes = [wintypes.HANDLE, wintypes.HANDLE,
                                 ctypes.POINTER(wintypes.BOOL)]
    k.IsProcessInJob.restype = wintypes.BOOL
    k.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    k.OpenProcess.restype = wintypes.HANDLE
    k.CloseHandle.argtypes = [wintypes.HANDLE]
    k.CloseHandle.restype = wintypes.BOOL
    return k


@contextlib.contextmanager
def _test_job(limit_flags: int):
    """指定した制限フラグのJobを作り、`(job, close)`を渡す。

    `close`は明示的にJobハンドルを閉じる（＝KILL_ON_JOB_CLOSEを発動させる）
    ためのもの。二重解放しないよう冪等にしてある。テストが途中で失敗しても
    finallyで必ず閉じるので、Jobハンドルもその中のプロセスも残らない。
    """
    k = _job_kernel32()
    job = k.CreateJobObjectW(None, None)
    assert job, f"CreateJobObjectWに失敗した: {ctypes.get_last_error()}"
    state = {"closed": False}

    def close():
        if not state["closed"]:
            state["closed"] = True
            k.CloseHandle(job)

    try:
        limits = _TEST_JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
        limits.BasicLimitInformation.LimitFlags = limit_flags
        assert k.SetInformationJobObject(
            job, _JOB_OBJECT_EXTENDED_LIMIT_INFORMATION_CLASS,
            ctypes.byref(limits), ctypes.sizeof(limits)), \
            f"SetInformationJobObjectに失敗した: {ctypes.get_last_error()}"
        yield k, job, close
    finally:
        close()


def _is_in_job(k, process_handle, job) -> bool:
    flag = wintypes.BOOL()
    assert k.IsProcessInJob(process_handle, job, ctypes.byref(flag)), \
        f"IsProcessInJobに失敗した: {ctypes.get_last_error()}"
    return bool(flag.value)


def _stays_alive(identity: dict, duration: float = 2.0) -> bool:
    """identityのプロセスがduration秒のあいだ生き続けたらTrue。"""
    return not _wait_until_gone(identity, timeout=duration)


_DETACHER_SRC = """\
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from lipidmix.core.atomic_io import DomainError
from lipidmix.core.process_control import launch_detached

out_path = Path(sys.argv[2])
cwd = Path(sys.argv[3])
log_path = Path(sys.argv[4])
worker_script = sys.argv[5]
# argv[6] は任意のゲート。呼び出し側がこのプロセスをJobへ割り当て終えるまで
# 待たせるために使う。待たずに起動すると、まだJobの外にいる状態でworkerを
# 作ってしまい「切り離せた」ように見えるだけの検証になる。
gate = Path(sys.argv[6]) if len(sys.argv) > 6 else None
if gate is not None:
    deadline = time.monotonic() + 60.0
    while not gate.exists():
        if time.monotonic() > deadline:
            sys.exit(3)
        time.sleep(0.02)
try:
    info = launch_detached([sys.executable, worker_script], cwd=cwd,
                           log_path=log_path)
    out_path.write_text(json.dumps(info), encoding="utf-8")
except DomainError as exc:
    out_path.write_text(json.dumps({"error": exc.code}), encoding="utf-8")
"""

_WORKER_SRC = "import time\ntime.sleep(300)\n"


@windows_only
def test_detached_worker_survives_the_job_that_launched_it(tmp_path):
    """切り離したワーカーは、親プロセスの終了でもJobの終了でも死なない。

    親の終了だけを見ても検証にならない。Windowsには親の死を子へ伝える仕組みが
    無く、孤児は放置されるだけなので、切り離しを全部外しても「親が終わった後も
    生きている」は成立してしまう。実際にMCP呼び出しの終了で解析を道連れにするのは
    Job Objectの KILL_ON_JOB_CLOSE なので、ここでは

      BREAKAWAY_OK | KILL_ON_JOB_CLOSE のJobを自前で作る
        → その中でdetacherを走らせ、launch_detachedにworkerを起こさせる
        → detacherの終了を待ってから、Jobハンドルを閉じる

    という順で、ワーカーが本当にJobの外へ出たかを見る。
    `launch_detached`が CREATE_BREAKAWAY_FROM_JOB を付けなければ、ワーカーは
    このJobに残り、Jobを閉じた瞬間に殺される（＝このテストは赤くなる）。
    """
    detacher = _write_script(tmp_path, "detacher.py", _DETACHER_SRC)
    worker_script = _write_script(tmp_path, "worker.py", _WORKER_SRC)
    out_path = tmp_path / "detached.json"
    gate = tmp_path / "assigned.gate"

    identity = None
    try:
        with _test_job(_TEST_JOB_LIMIT_BREAKAWAY_OK
                       | _TEST_JOB_LIMIT_KILL_ON_JOB_CLOSE) as (k, job, close_job):
            # venvのリダイレクタを挟むと、Jobへ入るのはリダイレクタで、detacher
            # 本体はその内側のJob（0x3000）から切り離し判定をしてしまう（GitHub #3）。
            command, env = pc.resolve_python_launch(
                [sys.executable, str(detacher), _REPO_ROOT, str(out_path),
                 str(tmp_path), str(tmp_path / "worker.log"),
                 str(worker_script), str(gate)])
            launcher = subprocess.Popen(
                command, env=env,
                cwd=str(tmp_path), stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            try:
                # AssignProcessToJobObject には SET_QUOTA と TERMINATE が要る。
                launcher_handle = k.OpenProcess(
                    _PROCESS_SET_QUOTA | _PROCESS_TERMINATE
                    | _PROCESS_QUERY_LIMITED_INFORMATION, False, launcher.pid)
                assert launcher_handle, \
                    f"OpenProcessに失敗した: {ctypes.get_last_error()}"
                try:
                    assert k.AssignProcessToJobObject(job, launcher_handle), \
                        f"detacherをJobへ入れられなかった: {ctypes.get_last_error()}"
                    assert _is_in_job(k, launcher_handle, job), \
                        "detacherが自前Jobに入っていない（前提が崩れている）"
                finally:
                    k.CloseHandle(launcher_handle)
                # Jobへの割当が済んでから初めて切り離しを始めさせる。
                gate.write_text("go", encoding="utf-8")
                assert launcher.wait(timeout=60) == 0
            finally:
                if launcher.poll() is None:
                    launcher.kill()
                    launcher.wait(timeout=20)

            info = json.loads(out_path.read_text(encoding="utf-8"))
            assert "error" not in info, f"切り離しに失敗した: {info}"
            identity = info["identity"]
            assert info["pid"] == identity["pid"]
            # このJobの制限では breakaway_ok 以外を選んではならない。
            assert info["breakaway"] == "breakaway_ok"

            worker_handle = k.OpenProcess(
                _PROCESS_QUERY_LIMITED_INFORMATION | _SYNCHRONIZE, False,
                int(info["pid"]))
            assert worker_handle, \
                f"ワーカーのOpenProcessに失敗した: {ctypes.get_last_error()}"
            try:
                assert not _is_in_job(k, worker_handle, job), \
                    "ワーカーが親Jobに残っている（breakawayできていない）"
                # ここでKILL_ON_JOB_CLOSEが発動する。detacherは既に終了済み。
                close_job()
                assert _stays_alive(identity), \
                    "Jobを閉じたらワーカーが道連れになった"
            finally:
                k.CloseHandle(worker_handle)
    finally:
        _kill_pid_if_same(identity)
    assert _wait_until_gone(identity)


@windows_only
def test_detach_from_kill_on_close_job_is_unsupported(tmp_path):
    """親JobがKILL_ON_JOB_CLOSEなら、存続を保証できないので起動しない。"""
    detacher = _write_script(tmp_path, "detacher.py", _DETACHER_SRC)
    worker = _write_script(tmp_path, "worker.py", _WORKER_SRC)
    out_path = tmp_path / "detached.json"

    # start_owned_processのJobはKILL_ON_JOB_CLOSE付き。その中で切り離しを試させる。
    proc = start_owned_process(
        [sys.executable, str(detacher), _REPO_ROOT, str(out_path), str(tmp_path),
         str(tmp_path / "worker.log"), str(worker)],
        cwd=tmp_path, log_path=tmp_path / "detacher.log")
    try:
        assert proc.wait(timeout=60) == 0
    finally:
        proc.close()

    info = json.loads(out_path.read_text(encoding="utf-8"))
    assert info == {"error": "DETACH_UNSUPPORTED"}


@windows_only
def test_detach_breakaway_mode_is_a_known_value():
    assert detach_breakaway_mode() in {
        "not_in_job", "inherited_job", "breakaway_ok", "silent_breakaway",
        "unsupported"}


@windows_only
def test_launch_detached_validates_arguments(tmp_path):
    with pytest.raises(DomainError) as excinfo:
        launch_detached([], cwd=tmp_path, log_path=tmp_path / "c.log")
    assert excinfo.value.code == "COMMAND_EMPTY"


# ---------- venvのリダイレクタを挟まずに本物の解釈系を起動する（GitHub #3） ----------
#
# Windowsのvenvの`Scripts\python.exe`は解釈系そのものではなく、自前のJob Object
# （SILENT_BREAKAWAY_OK | KILL_ON_JOB_CLOSE）を作って本物の`python.exe`を子として
# 起動するリダイレクタである。`[sys.executable, ...]`をそのまま起動すると、
# 返るpid・identity・Job割当・待ち合わせの対象はリダイレクタになり、コードを
# 実行している本物のプロセスではなくなる。


def test_resolve_python_launch_is_unchanged_outside_a_venv(monkeypatch, tmp_path):
    exe = tmp_path / "python.exe"
    exe.write_bytes(b"")
    monkeypatch.setattr(sys, "executable", str(exe))
    monkeypatch.setattr(sys, "_base_executable", str(exe), raising=False)

    command, env = pc.resolve_python_launch([str(exe), "-m", "pkg.worker"])

    assert command == [str(exe), "-m", "pkg.worker"]
    assert env is None


@windows_only
def test_resolve_python_launch_bypasses_the_venv_redirector(monkeypatch, tmp_path):
    """venvのリダイレクタを、venvの見え方を保ったまま本物の解釈系へ置き換える。

    `__PYVENV_LAUNCHER__`はリダイレクタ自身が本物の解釈系へ渡す環境変数で、
    これがあると解釈系は`pyvenv.cfg`を読み、venvの`sys.prefix`・site-packages・
    `sys.executable`をそのまま再現する（CPythonの`multiprocessing`と同じ手順）。
    """
    redirector = tmp_path / "venv" / "Scripts" / "python.exe"
    base = tmp_path / "base" / "python.exe"
    for path in (redirector, base):
        path.parent.mkdir(parents=True)
        path.write_bytes(b"")
    monkeypatch.setattr(sys, "executable", str(redirector))
    monkeypatch.setattr(sys, "_base_executable", str(base), raising=False)
    monkeypatch.setenv("LIPIDMIX_PROBE_INHERITED", "kept")
    before = dict(os.environ)

    # パスの大小文字違いは同じファイルとして扱う（Windowsのパスは大小を区別しない）。
    command, env = pc.resolve_python_launch(
        [str(redirector).upper(), "-m", "pkg.worker", "--flag"])

    assert command == [str(base), "-m", "pkg.worker", "--flag"]
    assert env["__PYVENV_LAUNCHER__"] == str(redirector)
    assert env["LIPIDMIX_PROBE_INHERITED"] == "kept"
    assert dict(os.environ) == before, "呼び出し元の環境を書き換えてはいけない"


@windows_only
def test_environment_block_encodes_characters_outside_the_bmp():
    """BMP外の文字（絵文字など）を含む環境変数があっても環境ブロックを作れる。

    venv経路は親の環境全体をブロックにするので、どれか1つの変数にBMP外の
    文字があるだけで、venvからの起動がすべて失敗してはいけない。wchar_tは
    UTF-16なので、BMP外の1文字はサロゲートペアの2要素になる。
    """
    block = pc._environment_block({"B": "y", "A": "x\U0001F600"})

    decoded = bytes(block).decode("utf-16-le")
    assert decoded.startswith("A=x\U0001F600\0B=y\0\0")


@windows_only
def test_resolve_python_launch_leaves_other_programs_alone(monkeypatch, tmp_path):
    redirector = tmp_path / "venv" / "Scripts" / "python.exe"
    base = tmp_path / "base" / "python.exe"
    for path in (redirector, base):
        path.parent.mkdir(parents=True)
        path.write_bytes(b"")
    monkeypatch.setattr(sys, "executable", str(redirector))
    monkeypatch.setattr(sys, "_base_executable", str(base), raising=False)
    console = str(tmp_path / "MsdialConsoleApp.exe")

    command, env = pc.resolve_python_launch([console, "lcmsdda", "-i", "x"])

    assert command == [console, "lcmsdda", "-i", "x"]
    assert env is None


_VENV_PROBE_SRC = """\
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from lipidmix.core.process_control import launch_detached, start_owned_process

work = Path(sys.argv[2])
child = work / "child.py"
child.write_text(
    "import json, os, sys\\n"
    "from pathlib import Path\\n"
    "Path(sys.argv[1]).write_text(json.dumps("
    "{'pid': os.getpid(), 'prefix': sys.prefix}), encoding='utf-8')\\n",
    encoding="utf-8")


def reported(path):
    deadline = time.monotonic() + 30.0
    while time.monotonic() < deadline:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            text = ""
        if text.strip():
            try:
                return json.loads(text)
            except json.JSONDecodeError:
                pass  # 書きかけを読んだ。書き終わるまで待つ。
        time.sleep(0.02)
    raise SystemExit(f"child did not write {path}")


result = {"executable": sys.executable,
          "base_executable": getattr(sys, "_base_executable", sys.executable),
          "prefix": sys.prefix}
owned = start_owned_process([sys.executable, str(child), str(work / "owned.json")],
                            cwd=work, log_path=work / "owned.log")
try:
    assert owned.wait(timeout=60) == 0
    result["owned"] = {"pid": owned.pid, "child": reported(work / "owned.json")}
finally:
    owned.close()
info = launch_detached([sys.executable, str(child), str(work / "detached.json")],
                       cwd=work, log_path=work / "detached.log")
result["detached"] = {"pid": info["pid"], "child": reported(work / "detached.json")}
(work / "result.json").write_text(json.dumps(result), encoding="utf-8")
"""


@windows_only
def test_launched_pid_is_the_interpreter_that_runs_the_code_even_from_a_venv(tmp_path):
    """venvから起動しても、返るpidはコードを実行している本物の解釈系のもの。

    テスト自身がvenvを作り、その`Scripts\\python.exe`（リダイレクタ）から
    `start_owned_process`/`launch_detached`を呼ばせる。システムの解釈系で
    テストを回していてもvenv経路を必ず通るので、リダイレクタのpidを返す
    退行はここで赤くなる。子の`sys.prefix`がvenvのままであることも見る
    ——リダイレクタを外した代わりにvenvのsite-packagesを失っては意味が無い。
    """
    import venv

    env_dir = tmp_path / "venv"
    venv.create(env_dir, with_pip=False)
    venv_python = env_dir / "Scripts" / "python.exe"
    work = tmp_path / "work"
    work.mkdir()
    probe = _write_script(tmp_path, "venv_probe.py", _VENV_PROBE_SRC)

    completed = subprocess.run([str(venv_python), str(probe), _REPO_ROOT, str(work)],
                               cwd=str(tmp_path), stdin=subprocess.DEVNULL,
                               capture_output=True, text=True, timeout=120)
    assert completed.returncode == 0, completed.stderr

    result = json.loads((work / "result.json").read_text(encoding="utf-8"))
    # 前提: venvのpython.exeはリダイレクタである（でなければ検証にならない）。
    assert os.path.normcase(result["executable"]) != \
        os.path.normcase(result["base_executable"])
    for key in ("owned", "detached"):
        assert result[key]["pid"] == result[key]["child"]["pid"], \
            f"{key}: 返したpidがコードを実行しているプロセスではない（リダイレクタ）"
        assert os.path.normcase(result[key]["child"]["prefix"]) == \
            os.path.normcase(result["prefix"]), f"{key}: 子がvenvを失っている"
