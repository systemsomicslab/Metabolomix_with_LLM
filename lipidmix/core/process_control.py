"""Windowsの解析プロセスの所有権・identity・OS排他。

依存グラフのleaf（stdlibと`lipidmix.core.atomic_io`のみ）。`lipidmix.core.mcp_core` /
`lipidmix.<形式>.tools` / `lipidmix.tools.*` をimportしてはいけない（循環）。

このモジュールが解く問題は三つある。

1. **所有権**: Consoleは孫プロセスを作る。親のpidだけを覚えて後で`kill`しても
   孫が残り、出力ファイルを掴んだままになる。Windows Job Objectへ入れておけば
   `TerminateJobObject`で一族をまとめて終了できる。名前による一括終了
   （`taskkill /IM`）は同名の無関係なプロセスを巻き込むので使わない。
2. **identity**: pidは再利用される。「pidが生きている」だけでは、監視対象が
   別プロセスに入れ替わったことを検出できず、無関係なプロセスを殺しうる。
   `GetProcessTimes`のcreation FILETIMEをpidと組にして同一性を判定する。
3. **排他**: 複数のワーカー・MCP呼び出しが同じ索引を同時に書くのを防ぐ。
   OSのbyte-range lockを使う。「owner記録ファイルを消す＝解除」ではない
   （記録が消えてもロックは残るし、記録が残っていてもロックは解けている）。

起動する子が`sys.executable`のときは`resolve_python_launch`を通す。venvの
`Scripts\\python.exe`は解釈系ではなくリダイレクタで、本物の解釈系を自前のJobの
中で子として走らせる。そのまま起動すると、所有・identity・切り離しの対象が
コードを実行していないリダイレクタになる（GitHub #3）。

`same_process`が「生きている」と答えるための判定に`OpenProcess`の成否だけを
使ってはいけない。実測（Windows 11 / Python 3.14）では、終了済みプロセスでも
そのプロセスハンドルがどこかで開かれている限り`OpenProcess`は成功する。
そのため必ず`WaitForSingleObject(h, 0)`で signaled（＝終了済み）を除外する。
"""
from __future__ import annotations

import contextlib
import os
import subprocess
import sys
import time
from pathlib import Path

from lipidmix.core.atomic_io import DomainError

__all__ = [
    "OwnedProcess",
    "detach_breakaway_mode",
    "file_lock",
    "launch_detached",
    "process_identity",
    "resolve_python_launch",
    "same_process",
    "start_owned_process",
]

# --- CreateProcessW / Job Object のフラグ（winbase.h） ---
CREATE_SUSPENDED = 0x00000004
CREATE_NO_WINDOW = 0x08000000
CREATE_BREAKAWAY_FROM_JOB = 0x01000000
JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x00002000
JOB_OBJECT_LIMIT_BREAKAWAY_OK = 0x00000800
JOB_OBJECT_LIMIT_SILENT_BREAKAWAY_OK = 0x00001000
DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200
CREATE_UNICODE_ENVIRONMENT = 0x00000400

#: venvのリダイレクタが本物の解釈系へ渡す環境変数。解釈系はこれを見て
#: `pyvenv.cfg`を読み、venvの`sys.prefix`・site-packages・`sys.executable`を
#: 再現する（起動直後に自分で消すので、さらに下の子へは漏れない）。
_PYVENV_LAUNCHER_ENV = "__PYVENV_LAUNCHER__"

#: Job Objectの終了コード。TerminateJobObjectへ渡し、所有プロセスの
#: exit_codeとして観測される（0以外であることに意味がある）。
JOB_TERMINATION_EXIT_CODE = 1

_IS_WINDOWS = os.name == "nt"


def _unsupported_platform(operation: str) -> DomainError:
    return DomainError(
        "PLATFORM_UNSUPPORTED",
        f"{operation} はWindowsでのみ利用できます（現在: {sys.platform}）",
        {"operation": operation, "platform": sys.platform},
    )


def _same_path(a: str, b: str) -> bool:
    return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


def resolve_python_launch(command: list[str]) -> tuple[list[str], dict | None]:
    """`[sys.executable, ...]`を、venvのリダイレクタを挟まない起動に直す。

    返り値は`(command, env)`。`env`がNoneなら環境は親のまま継承してよく、
    dictならその環境全体で起動する（`subprocess.Popen(command, env=env)`に
    そのまま渡せる形）。

    なぜ要るか: Windowsのvenvの`Scripts\\python.exe`は解釈系そのものではない。
    自前のJob Object（SILENT_BREAKAWAY_OK | KILL_ON_JOB_CLOSE）を作り、本物の
    `python.exe`を子として起動して終了を待つだけのリダイレクタである。これを
    そのまま起動すると、手に入るpid・プロセスハンドルはリダイレクタのものになり、

    - identityの記録・`same_process`の生死判定・待ち合わせが、コードを実行して
      いないプロセスを見る（リダイレクタが signaled になった時点で、本物の
      workerはまだ終了処理中でありうる）。
    - 自前Jobへの割当・切り離し判定が、リダイレクタの内側のJobに入った本物の
      解釈系に届かない（解釈系からは常にリダイレクタのJob＝0x3000が見える）。

    （GitHub #3。システムの解釈系を直接使う環境では起きない。）

    直し方はCPythonの`multiprocessing`（`popen_spawn_win32`）と同じで、
    `sys._base_executable`（本物の解釈系）を起動し、`__PYVENV_LAUNCHER__`に
    venvの`python.exe`を入れて渡す。これはリダイレクタ自身が本物の解釈系へ
    渡す値そのものなので、子から見たvenv（`sys.prefix`・site-packages・
    `sys.executable`）はリダイレクタ経由と変わらない。PYTHONPATHで
    site-packagesを足す方式は`.pth`の処理や`sys.prefix`が再現されないので採らない。

    書き換えるのは`command[0]`が`sys.executable`そのもので、かつ
    `sys._base_executable`がそれと異なり実在するときだけ。Console本体など
    他のプログラム、明示された別の解釈系、venvでない環境では何もしない。
    """
    command = [str(part) for part in command]
    if not _IS_WINDOWS or not command:
        return command, None
    executable = sys.executable
    base = getattr(sys, "_base_executable", None)
    if (not executable or not base or _same_path(base, executable)
            or not _same_path(command[0], executable)
            or not os.path.isfile(base)):
        return command, None
    env = dict(os.environ)
    env[_PYVENV_LAUNCHER_ENV] = executable
    return [base, *command[1:]], env


if _IS_WINDOWS:  # pragma: no branch - Windows専用の宣言ブロック
    import ctypes
    import ctypes.wintypes as _wt

    _kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)

    # 64bit HANDLE を int32 へ切り詰めないため、全関数にargtypes/restypeを宣言する。
    # 宣言を省くとctypesは戻り値をc_intとして扱い、上位32bitを失ったハンドルで
    # CloseHandleを呼ぶ（＝ハンドルリークと、無関係なハンドルの破壊）。
    _ULONG_PTR = ctypes.c_size_t

    class _SECURITY_ATTRIBUTES(ctypes.Structure):
        _fields_ = [
            ("nLength", _wt.DWORD),
            ("lpSecurityDescriptor", ctypes.c_void_p),
            ("bInheritHandle", _wt.BOOL),
        ]

    class _STARTUPINFOW(ctypes.Structure):
        _fields_ = [
            ("cb", _wt.DWORD),
            ("lpReserved", _wt.LPWSTR),
            ("lpDesktop", _wt.LPWSTR),
            ("lpTitle", _wt.LPWSTR),
            ("dwX", _wt.DWORD),
            ("dwY", _wt.DWORD),
            ("dwXSize", _wt.DWORD),
            ("dwYSize", _wt.DWORD),
            ("dwXCountChars", _wt.DWORD),
            ("dwYCountChars", _wt.DWORD),
            ("dwFillAttribute", _wt.DWORD),
            ("dwFlags", _wt.DWORD),
            ("wShowWindow", _wt.WORD),
            ("cbReserved2", _wt.WORD),
            ("lpReserved2", ctypes.POINTER(ctypes.c_byte)),
            ("hStdInput", _wt.HANDLE),
            ("hStdOutput", _wt.HANDLE),
            ("hStdError", _wt.HANDLE),
        ]

    class _PROCESS_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("hProcess", _wt.HANDLE),
            ("hThread", _wt.HANDLE),
            ("dwProcessId", _wt.DWORD),
            ("dwThreadId", _wt.DWORD),
        ]

    class _IO_COUNTERS(ctypes.Structure):
        _fields_ = [
            ("ReadOperationCount", ctypes.c_ulonglong),
            ("WriteOperationCount", ctypes.c_ulonglong),
            ("OtherOperationCount", ctypes.c_ulonglong),
            ("ReadTransferCount", ctypes.c_ulonglong),
            ("WriteTransferCount", ctypes.c_ulonglong),
            ("OtherTransferCount", ctypes.c_ulonglong),
        ]

    class _JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", _wt.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", _wt.DWORD),
            ("Affinity", _ULONG_PTR),
            ("PriorityClass", _wt.DWORD),
            ("SchedulingClass", _wt.DWORD),
        ]

    class _JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", _JOBOBJECT_BASIC_LIMIT_INFORMATION),
            ("IoInfo", _IO_COUNTERS),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    class _JOBOBJECT_BASIC_ACCOUNTING_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("TotalUserTime", ctypes.c_longlong),
            ("TotalKernelTime", ctypes.c_longlong),
            ("ThisPeriodTotalUserTime", ctypes.c_longlong),
            ("ThisPeriodTotalKernelTime", ctypes.c_longlong),
            ("TotalPageFaultCount", _wt.DWORD),
            ("TotalProcesses", _wt.DWORD),
            ("ActiveProcesses", _wt.DWORD),
            ("TotalTerminatedProcesses", _wt.DWORD),
        ]

    _JobObjectBasicAccountingInformation = 1
    _JobObjectExtendedLimitInformation = 9

    _INFINITE = 0xFFFFFFFF
    _WAIT_OBJECT_0 = 0x00000000
    _WAIT_TIMEOUT = 0x00000102
    _GENERIC_READ = 0x80000000
    _GENERIC_WRITE = 0x40000000
    _FILE_SHARE_ALL = 0x00000007
    _CREATE_ALWAYS = 2
    _OPEN_EXISTING = 3
    _FILE_ATTRIBUTE_NORMAL = 0x00000080
    _INVALID_HANDLE_VALUE = _wt.HANDLE(-1).value
    _STARTF_USESTDHANDLES = 0x00000100
    _PROCESS_QUERY_LIMITED_INFORMATION = 0x00001000
    _SYNCHRONIZE = 0x00100000

    _CreateJobObjectW = _kernel32.CreateJobObjectW
    _CreateJobObjectW.argtypes = [ctypes.POINTER(_SECURITY_ATTRIBUTES), _wt.LPCWSTR]
    _CreateJobObjectW.restype = _wt.HANDLE

    _SetInformationJobObject = _kernel32.SetInformationJobObject
    _SetInformationJobObject.argtypes = [_wt.HANDLE, ctypes.c_int, ctypes.c_void_p,
                                         _wt.DWORD]
    _SetInformationJobObject.restype = _wt.BOOL

    _QueryInformationJobObject = _kernel32.QueryInformationJobObject
    _QueryInformationJobObject.argtypes = [_wt.HANDLE, ctypes.c_int, ctypes.c_void_p,
                                           _wt.DWORD, ctypes.POINTER(_wt.DWORD)]
    _QueryInformationJobObject.restype = _wt.BOOL

    _AssignProcessToJobObject = _kernel32.AssignProcessToJobObject
    _AssignProcessToJobObject.argtypes = [_wt.HANDLE, _wt.HANDLE]
    _AssignProcessToJobObject.restype = _wt.BOOL

    _TerminateJobObject = _kernel32.TerminateJobObject
    _TerminateJobObject.argtypes = [_wt.HANDLE, _wt.UINT]
    _TerminateJobObject.restype = _wt.BOOL

    _IsProcessInJob = _kernel32.IsProcessInJob
    _IsProcessInJob.argtypes = [_wt.HANDLE, _wt.HANDLE, ctypes.POINTER(_wt.BOOL)]
    _IsProcessInJob.restype = _wt.BOOL

    _CreateProcessW = _kernel32.CreateProcessW
    _CreateProcessW.argtypes = [
        _wt.LPCWSTR,                              # lpApplicationName
        _wt.LPWSTR,                               # lpCommandLine（書き換えられる）
        ctypes.POINTER(_SECURITY_ATTRIBUTES),     # lpProcessAttributes
        ctypes.POINTER(_SECURITY_ATTRIBUTES),     # lpThreadAttributes
        _wt.BOOL,                                 # bInheritHandles
        _wt.DWORD,                                # dwCreationFlags
        ctypes.c_void_p,                          # lpEnvironment
        _wt.LPCWSTR,                              # lpCurrentDirectory
        ctypes.POINTER(_STARTUPINFOW),
        ctypes.POINTER(_PROCESS_INFORMATION),
    ]
    _CreateProcessW.restype = _wt.BOOL

    _ResumeThread = _kernel32.ResumeThread
    _ResumeThread.argtypes = [_wt.HANDLE]
    _ResumeThread.restype = _wt.DWORD

    _TerminateProcess = _kernel32.TerminateProcess
    _TerminateProcess.argtypes = [_wt.HANDLE, _wt.UINT]
    _TerminateProcess.restype = _wt.BOOL

    _CloseHandle = _kernel32.CloseHandle
    _CloseHandle.argtypes = [_wt.HANDLE]
    _CloseHandle.restype = _wt.BOOL

    _WaitForSingleObject = _kernel32.WaitForSingleObject
    _WaitForSingleObject.argtypes = [_wt.HANDLE, _wt.DWORD]
    _WaitForSingleObject.restype = _wt.DWORD

    _GetExitCodeProcess = _kernel32.GetExitCodeProcess
    _GetExitCodeProcess.argtypes = [_wt.HANDLE, ctypes.POINTER(_wt.DWORD)]
    _GetExitCodeProcess.restype = _wt.BOOL

    _OpenProcess = _kernel32.OpenProcess
    _OpenProcess.argtypes = [_wt.DWORD, _wt.BOOL, _wt.DWORD]
    _OpenProcess.restype = _wt.HANDLE

    _GetProcessTimes = _kernel32.GetProcessTimes
    _GetProcessTimes.argtypes = [_wt.HANDLE, ctypes.POINTER(_wt.FILETIME),
                                 ctypes.POINTER(_wt.FILETIME),
                                 ctypes.POINTER(_wt.FILETIME),
                                 ctypes.POINTER(_wt.FILETIME)]
    _GetProcessTimes.restype = _wt.BOOL

    _GetCurrentProcess = _kernel32.GetCurrentProcess
    _GetCurrentProcess.argtypes = []
    _GetCurrentProcess.restype = _wt.HANDLE

    _CreateFileW = _kernel32.CreateFileW
    _CreateFileW.argtypes = [_wt.LPCWSTR, _wt.DWORD, _wt.DWORD,
                             ctypes.POINTER(_SECURITY_ATTRIBUTES),
                             _wt.DWORD, _wt.DWORD, _wt.HANDLE]
    _CreateFileW.restype = _wt.HANDLE


def _win_error(code: str, message: str, details: dict | None = None) -> DomainError:
    """直前のGetLastErrorを添えてDomainErrorを組む。"""
    err = ctypes.get_last_error()
    payload = dict(details or {})
    payload["win32_error"] = err
    payload["win32_message"] = ctypes.FormatError(err).strip()
    return DomainError(code, f"{message}（Win32 error {err}: "
                             f"{payload['win32_message']}）", payload)


def _creation_time(handle) -> int:
    """プロセスハンドルの生成FILETIMEを100ns刻みの整数で返す。"""
    creation = _wt.FILETIME()
    exit_ft = _wt.FILETIME()
    kernel_ft = _wt.FILETIME()
    user_ft = _wt.FILETIME()
    if not _GetProcessTimes(handle, ctypes.byref(creation), ctypes.byref(exit_ft),
                            ctypes.byref(kernel_ft), ctypes.byref(user_ft)):
        raise _win_error("PROCESS_TIMES_FAILED", "GetProcessTimesに失敗した")
    return (int(creation.dwHighDateTime) << 32) | int(creation.dwLowDateTime)


def _has_exited(handle) -> bool:
    """ハンドルが signaled（＝プロセス終了済み）かを返す。"""
    return _WaitForSingleObject(handle, 0) == _WAIT_OBJECT_0


def _exit_code(handle) -> int:
    code = _wt.DWORD()
    if not _GetExitCodeProcess(handle, ctypes.byref(code)):
        raise _win_error("PROCESS_EXIT_CODE_FAILED", "GetExitCodeProcessに失敗した")
    return int(code.value)


def _inheritable_security_attributes():
    sa = _SECURITY_ATTRIBUTES()
    sa.nLength = ctypes.sizeof(_SECURITY_ATTRIBUTES)
    sa.lpSecurityDescriptor = None
    sa.bInheritHandle = True
    return sa


def _open_child_stdio(log_path: Path):
    """子へ渡す stdout用ログハンドルと stdin用NULハンドルを開いて返す。

    どちらか一方の作成に失敗したら、既に開いた側も閉じてから送出する
    （ここで漏らすと、失敗を繰り返すたびにハンドルが積み上がる）。
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    sa = _inheritable_security_attributes()
    log_handle = _CreateFileW(str(log_path), _GENERIC_WRITE, _FILE_SHARE_ALL,
                              ctypes.byref(sa), _CREATE_ALWAYS,
                              _FILE_ATTRIBUTE_NORMAL, None)
    if not log_handle or log_handle == _INVALID_HANDLE_VALUE:
        raise _win_error("LOG_OPEN_FAILED", f"ログファイルを開けない: {log_path}",
                         {"log_path": str(log_path)})
    null_handle = _CreateFileW("NUL", _GENERIC_READ, _FILE_SHARE_ALL,
                               ctypes.byref(sa), _OPEN_EXISTING,
                               _FILE_ATTRIBUTE_NORMAL, None)
    if not null_handle or null_handle == _INVALID_HANDLE_VALUE:
        error = _win_error("NUL_OPEN_FAILED", "NULデバイスを開けない")
        _CloseHandle(log_handle)
        raise error
    return log_handle, null_handle


def _startup_info(log_handle, null_handle) -> "_STARTUPINFOW":
    si = _STARTUPINFOW()
    si.cb = ctypes.sizeof(_STARTUPINFOW)
    si.dwFlags = _STARTF_USESTDHANDLES
    si.hStdInput = null_handle
    si.hStdOutput = log_handle
    si.hStdError = log_handle
    return si


def _validate_launch_args(command, cwd: Path) -> tuple[list[str], Path]:
    command = [str(part) for part in command]
    if not command:
        raise DomainError("COMMAND_EMPTY", "起動コマンドが空である")
    cwd = Path(cwd)
    if not cwd.is_dir():
        raise DomainError("CWD_NOT_FOUND", f"作業ディレクトリが存在しない: {cwd}",
                          {"cwd": str(cwd)})
    return command, cwd


def _environment_block(env: dict):
    """CreateProcessW用のUnicode環境ブロック（`名前=値\\0`の列＋終端`\\0`）を作る。

    名前は大文字小文字を無視した順に並べる（Win32の要求）。NULを含む名前・値と、
    `=`を含む名前はブロックを壊すので拒否する。
    """
    entries = []
    for name, value in sorted(env.items(), key=lambda kv: str(kv[0]).upper()):
        name, value = str(name), str(value)
        if not name or "=" in name or "\0" in name or "\0" in value:
            raise DomainError("ENVIRONMENT_INVALID",
                              f"子へ渡せない環境変数がある: {name!r}",
                              {"name": name})
        entries.append(f"{name}={value}\0")
    block = "".join(entries) + "\0"
    buffer = (ctypes.c_wchar * len(block))()
    buffer[:] = block
    return buffer


def _create_process(command: list[str], cwd: Path, log_path: Path,
                    creation_flags: int,
                    env: dict | None = None) -> "_PROCESS_INFORMATION":
    """CreateProcessWを呼び、PROCESS_INFORMATIONを返す。

    stdio用に開いたハンドルは成功・失敗いずれの経路でも必ず閉じる。子は
    CreateProcessWの時点で自分用の複製を受け取っているので、親側の複製を
    残しておく理由がない（残すとログファイルが解放されない）。

    `env`がNoneなら親の環境を継承する。dictならその環境全体で起動する
    （`resolve_python_launch`がvenvを再現するために使う）。
    """
    environment = None
    if env is not None:
        environment = _environment_block(env)
        creation_flags |= CREATE_UNICODE_ENVIRONMENT
    log_handle, null_handle = _open_child_stdio(log_path)
    try:
        si = _startup_info(log_handle, null_handle)
        pi = _PROCESS_INFORMATION()
        cmdline = ctypes.create_unicode_buffer(subprocess.list2cmdline(command))
        created = _CreateProcessW(None, cmdline, None, None, True,
                                  creation_flags, environment, str(cwd),
                                  ctypes.byref(si), ctypes.byref(pi))
        if not created:
            raise _win_error("PROCESS_LAUNCH_FAILED",
                             "CreateProcessWに失敗した",
                             {"command": command, "cwd": str(cwd)})
        return pi
    finally:
        _CloseHandle(log_handle)
        _CloseHandle(null_handle)


def process_identity(pid: int) -> dict | None:
    """生きているプロセスの identity（pidと生成時刻）を返す。無ければNone。

    「生きている」の判定に`OpenProcess`の成否だけを使わない。終了済みでも
    プロセスハンドルが残っている間は`OpenProcess`が成功するため、必ず
    signaled 状態（＝終了済み）を除外する。
    """
    if not _IS_WINDOWS:
        raise _unsupported_platform("process_identity")
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return None
    if pid <= 0:
        return None
    handle = _OpenProcess(_PROCESS_QUERY_LIMITED_INFORMATION | _SYNCHRONIZE,
                          False, pid)
    if not handle:
        return None
    try:
        if _has_exited(handle):
            return None
        return {"pid": pid, "creation_time": _creation_time(handle)}
    finally:
        _CloseHandle(handle)


def same_process(expected: dict) -> bool:
    """`expected` の identity と同じプロセスが今も動いているかを返す。

    pidが同じでも生成時刻が違えば False（pid再利用の見分け）。これが
    Task 16の「報告が途絶えたワーカーを殺してよいか」の判断根拠になる。
    """
    if not isinstance(expected, dict):
        return False
    pid = expected.get("pid")
    expected_creation = expected.get("creation_time")
    if pid is None or expected_creation is None:
        return False
    current = process_identity(pid)
    if current is None:
        return False
    return current["creation_time"] == expected_creation


class OwnedProcess:
    """Job Objectで一族ごと所有する子プロセス。

    `close()`でJobハンドルを閉じると、KILL_ON_JOB_CLOSEにより取り残された
    子孫も終了する。孤児を作らないための最後の砦なので、この制限は外さない。
    """

    def __init__(self, job_handle, process_handle, pid: int, identity: dict,
                 log_path: Path):
        self._job = job_handle
        self._process = process_handle
        self._pid = int(pid)
        self._identity = dict(identity)
        self._log_path = Path(log_path)
        self._closed = False

    @property
    def pid(self) -> int:
        return self._pid

    @property
    def identity(self) -> dict:
        """pidと生成時刻の組（呼び出し側の書き換えを避けるためコピーを返す）。"""
        return dict(self._identity)

    @property
    def log_path(self) -> Path:
        return self._log_path

    def _require_open(self) -> None:
        if self._closed:
            raise DomainError("PROCESS_HANDLE_CLOSED",
                              "close()済みのOwnedProcessは操作できない",
                              {"pid": self._pid})

    def poll(self) -> int | None:
        """終了していれば終了コード、動作中ならNone。"""
        self._require_open()
        status = _WaitForSingleObject(self._process, 0)
        if status == _WAIT_TIMEOUT:
            return None
        if status == _WAIT_OBJECT_0:
            return _exit_code(self._process)
        raise _win_error("PROCESS_WAIT_FAILED", "WaitForSingleObjectに失敗した",
                         {"pid": self._pid})

    def wait(self, timeout: float | None = None) -> int:
        """終了を待って終了コードを返す。timeout超過はDomainError。"""
        self._require_open()
        millis = _INFINITE if timeout is None else max(0, int(timeout * 1000))
        status = _WaitForSingleObject(self._process, millis)
        if status == _WAIT_OBJECT_0:
            return _exit_code(self._process)
        if status == _WAIT_TIMEOUT:
            raise DomainError("PROCESS_WAIT_TIMEOUT",
                              f"プロセスがtimeout={timeout}s以内に終了しなかった",
                              {"pid": self._pid, "timeout_s": timeout})
        raise _win_error("PROCESS_WAIT_FAILED", "WaitForSingleObjectに失敗した",
                         {"pid": self._pid})

    def _active_processes(self) -> int:
        info = _JOBOBJECT_BASIC_ACCOUNTING_INFORMATION()
        returned = _wt.DWORD()
        if not _QueryInformationJobObject(self._job,
                                          _JobObjectBasicAccountingInformation,
                                          ctypes.byref(info), ctypes.sizeof(info),
                                          ctypes.byref(returned)):
            raise _win_error("JOB_QUERY_FAILED",
                             "QueryInformationJobObjectに失敗した",
                             {"pid": self._pid})
        return int(info.ActiveProcesses)

    def terminate_tree(self, timeout: float = 15.0) -> None:
        """Job内の全プロセス（孫を含む）を終了させ、全て消えるまで待つ。

        名前による一括終了（`taskkill /IM`）は使わない。同名の無関係な
        プロセスを巻き込むうえ、所有していないものまで殺してしまう。

        待ち条件が二つあるのは、片方だけでは不足だから。本機での実測では、
        `TerminateJobObject`の直後にJobの`ActiveProcesses`が0になっても、
        肝心のプロセスオブジェクトはまだ signaled でない瞬間がある
        （そこで戻ると、`same_process`が終了済みのはずのプロセスを「生きている」
        と答える）。逆にプロセスハンドルの待ちだけでは、Job内の孫の終了を
        取りこぼす。両方が揃うまで待つ。
        """
        self._require_open()
        if not _TerminateJobObject(self._job, JOB_TERMINATION_EXIT_CODE):
            raise _win_error("JOB_TERMINATE_FAILED", "TerminateJobObjectに失敗した",
                             {"pid": self._pid})
        deadline = time.monotonic() + timeout
        while True:
            if _has_exited(self._process) and self._active_processes() == 0:
                return
            if time.monotonic() >= deadline:
                raise DomainError(
                    "JOB_TERMINATE_TIMEOUT",
                    f"Job内のプロセスがtimeout={timeout}s以内に終了しなかった",
                    {"pid": self._pid, "timeout_s": timeout})
            time.sleep(0.02)

    def close(self) -> None:
        """保持ハンドルを解放する（冪等）。

        Jobハンドルを最後に閉じる。KILL_ON_JOB_CLOSEが効くのはこの時点なので、
        `terminate_tree`を呼ばずに`close`しても孤児は残らない。
        """
        if self._closed:
            return
        self._closed = True
        _CloseHandle(self._process)
        _CloseHandle(self._job)


def start_owned_process(command: list[str], *, cwd: Path,
                        log_path: Path) -> OwnedProcess:
    """子プロセスを専用Job Objectの中で起動する。

    起動順は CreateJobObjectW -> SetInformationJobObject -> CreateProcessW
    （CREATE_SUSPENDED | CREATE_NO_WINDOW）-> AssignProcessToJobObject ->
    ResumeThread。**この順序を崩してはいけない**。停止状態で作らずに起動すると、
    Jobへ割り当てる前に子が孫を作りうる。その孫はJobの外に生まれるので、
    `TerminateJobObject`で回収できない取り残しになる。

    Assignに失敗したら停止中の子を`TerminateProcess`で始末する。resumeして
    しまうと、所有できていないプロセスを野に放つことになる。
    """
    if not _IS_WINDOWS:
        raise _unsupported_platform("start_owned_process")
    command, cwd = _validate_launch_args(command, cwd)
    command, env = resolve_python_launch(command)

    job = _CreateJobObjectW(None, None)
    if not job:
        raise _win_error("JOB_CREATE_FAILED", "CreateJobObjectWに失敗した")
    try:
        limits = _JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
        limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not _SetInformationJobObject(job, _JobObjectExtendedLimitInformation,
                                        ctypes.byref(limits),
                                        ctypes.sizeof(limits)):
            raise _win_error("JOB_LIMIT_FAILED",
                             "SetInformationJobObjectに失敗した")

        pi = _create_process(command, cwd, log_path,
                             CREATE_SUSPENDED | CREATE_NO_WINDOW, env)
        try:
            if not _AssignProcessToJobObject(job, pi.hProcess):
                error = _win_error("PROCESS_ASSIGN_FAILED",
                                   "AssignProcessToJobObjectに失敗した",
                                   {"pid": int(pi.dwProcessId)})
                _TerminateProcess(pi.hProcess, JOB_TERMINATION_EXIT_CODE)
                raise error
            if _ResumeThread(pi.hThread) == 0xFFFFFFFF:
                error = _win_error("PROCESS_RESUME_FAILED", "ResumeThreadに失敗した",
                                   {"pid": int(pi.dwProcessId)})
                _TerminateProcess(pi.hProcess, JOB_TERMINATION_EXIT_CODE)
                raise error
            identity = {"pid": int(pi.dwProcessId),
                        "creation_time": _creation_time(pi.hProcess)}
        except BaseException:
            _CloseHandle(pi.hProcess)
            raise
        finally:
            _CloseHandle(pi.hThread)
    except BaseException:
        _CloseHandle(job)
        raise

    return OwnedProcess(job, pi.hProcess, int(pi.dwProcessId), identity,
                        Path(log_path))


def _breakaway_mode(in_job: bool, limit_flags: int) -> str:
    """親Jobの制限から、切り離し起動が成立するかを判定する（純関数）。

    返す値:
      - `not_in_job`      Jobの外にいる。そのまま切り離せる。
      - `silent_breakaway` 子は自動的にJobから外れる。追加フラグは不要。
      - `breakaway_ok`     CREATE_BREAKAWAY_FROM_JOB を付ければ外れる。
      - `inherited_job`    外れられないが、親JobにKILL_ON_JOB_CLOSEが無いので
                           親の終了で道連れにはならない（存続を保証できる）。
      - `unsupported`      外れられず、かつ親Jobが閉じると殺される。

    `inherited_job` を許すのが要点。Windows 11ではシェル配下のプロセスが
    制限なしのJobに入っていることが普通にあり（本機で実測: in_job=1,
    LimitFlags=0x0）、breakaway権限だけを条件にすると切り離しが常に失敗する。
    危険なのは「外れられない」ことではなく「親Jobが閉じたら殺される」ことなので、
    KILL_ON_JOB_CLOSEの有無で判定する。
    """
    if not in_job:
        return "not_in_job"
    if limit_flags & JOB_OBJECT_LIMIT_SILENT_BREAKAWAY_OK:
        return "silent_breakaway"
    if limit_flags & JOB_OBJECT_LIMIT_BREAKAWAY_OK:
        return "breakaway_ok"
    if limit_flags & JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE:
        return "unsupported"
    return "inherited_job"


def detach_breakaway_mode() -> str:
    """現在のプロセスから見た切り離し可否を返す（`_breakaway_mode`の値）。"""
    if not _IS_WINDOWS:
        raise _unsupported_platform("detach_breakaway_mode")
    in_job = _wt.BOOL()
    if not _IsProcessInJob(_GetCurrentProcess(), None, ctypes.byref(in_job)):
        raise _win_error("JOB_QUERY_FAILED", "IsProcessInJobに失敗した")
    if not in_job.value:
        return "not_in_job"
    info = _JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
    returned = _wt.DWORD()
    # Jobハンドルに NULL を渡すと、現在のプロセスが属するJobを問い合わせる。
    if not _QueryInformationJobObject(None, _JobObjectExtendedLimitInformation,
                                      ctypes.byref(info), ctypes.sizeof(info),
                                      ctypes.byref(returned)):
        raise _win_error("JOB_QUERY_FAILED",
                         "QueryInformationJobObjectに失敗した")
    return _breakaway_mode(True, int(info.BasicLimitInformation.LimitFlags))


def launch_detached(command: list[str], *, cwd: Path, log_path: Path) -> dict:
    """親から切り離してワーカーを起動し、pidとidentityを返す（待たない）。

    切り離しの目的は「MCPの1呼び出しが消えても解析が続く」こと。親Jobから
    抜けられず、かつそのJobが閉じたら殺される環境では目的を果たせないので、
    起動せずに`DETACH_UNSUPPORTED`を送出する。黙って起動すると、親の終了と
    同時に消えるワーカーを「起動できた」と記録してしまう。
    """
    if not _IS_WINDOWS:
        raise _unsupported_platform("launch_detached")
    command, cwd = _validate_launch_args(command, cwd)
    command, env = resolve_python_launch(command)

    mode = detach_breakaway_mode()
    if mode == "unsupported":
        raise DomainError(
            "DETACH_UNSUPPORTED",
            "親Job ObjectがKILL_ON_JOB_CLOSE付きでbreakawayも許可していないため、"
            "親の終了後もワーカーが生き残ることを保証できない",
            {"command": command, "breakaway": mode})

    flags = DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
    if mode == "breakaway_ok":
        flags |= CREATE_BREAKAWAY_FROM_JOB

    pi = _create_process(command, cwd, log_path, flags, env)
    try:
        identity = {"pid": int(pi.dwProcessId),
                    "creation_time": _creation_time(pi.hProcess)}
    finally:
        _CloseHandle(pi.hThread)
        _CloseHandle(pi.hProcess)
    return {"pid": int(pi.dwProcessId), "identity": identity,
            "breakaway": mode, "log_path": str(Path(log_path))}


def _try_lock(fd: int) -> bool:
    """1byteの排他ロックを非ブロッキングで試みる。取れたらTrue。"""
    if _IS_WINDOWS:
        import msvcrt
        os.lseek(fd, 0, os.SEEK_SET)
        try:
            msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)
        except OSError:
            return False
        return True
    import fcntl
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return False
    return True


def _unlock(fd: int) -> None:
    if _IS_WINDOWS:
        import msvcrt
        os.lseek(fd, 0, os.SEEK_SET)
        try:
            msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
        except OSError:
            pass
        return
    import fcntl
    try:
        fcntl.flock(fd, fcntl.LOCK_UN)
    except OSError:
        pass


@contextlib.contextmanager
def file_lock(path: Path, *, timeout: float = 30.0, poll_interval: float = 0.05):
    """OSのbyte-range lockでプロセス間排他を取るコンテキストマネージャ。

    Windowsは`msvcrt.locking`で先頭1byteを、POSIXは`flock`で確保する。
    ロックの正体はOSが握る状態であって、ロックファイルの存在や中身では**ない**。
    ファイルが残っていてもロックは解けているし、記録を消してもロックは解けない。
    「owner記録を消したから解除された」という判断をしてはいけない。

    解放はwithを抜けた時だけ（例外で抜けた場合も含む）。timeout超過は
    `LOCK_TIMEOUT`。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_RDWR | os.O_CREAT
    flags |= getattr(os, "O_BINARY", 0)
    fd = os.open(path, flags, 0o644)
    acquired = False
    try:
        deadline = time.monotonic() + max(0.0, timeout)
        while True:
            if _try_lock(fd):
                acquired = True
                break
            if time.monotonic() >= deadline:
                raise DomainError(
                    "LOCK_TIMEOUT",
                    f"ロックをtimeout={timeout}s以内に取得できなかった: {path}",
                    {"path": str(path), "timeout_s": timeout})
            time.sleep(poll_interval)
        yield path
    finally:
        if acquired:
            _unlock(fd)
        os.close(fd)
