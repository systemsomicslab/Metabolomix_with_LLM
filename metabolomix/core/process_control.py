"""Windowsの解析プロセスの所有権・identity・OS排他。

依存グラフのleaf（stdlibと`metabolomix.core.atomic_io`のみ）。`metabolomix.core.mcp_core` /
`metabolomix.<形式>.tools` / `metabolomix.tools.*` をimportしてはいけない（循環）。

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

`launch_detached`はワーカーを直接起こさず、すぐ終わる中継のプロセス
（`metabolomix.core.detach_relay`）に起こさせる。ワーカーの親が残らないので、
起こした側（MCPサーバ）がプロセスの木ごと消されても（`taskkill /T`）ワーカーは残る。
"""
from __future__ import annotations

import contextlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from metabolomix.core.atomic_io import DomainError

__all__ = [
    "OwnedProcess",
    "commit_memory_status",
    "current_job_state",
    "current_process_in_job",
    "detach_breakaway_mode",
    "file_lock",
    "launch_detached",
    "process_identity",
    "resolve_python_launch",
    "same_process",
    "start_owned_process",
    "wait_for_lock",
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

#: `launch_detached`が起こす中継のモジュール（`python -m`で起こす）。
_RELAY_MODULE = "metabolomix.core.detach_relay"

#: 中継の作業ディレクトリ。`metabolomix`パッケージを含むディレクトリで、`-m`はここから
#: このチェックアウトのコードを読む（pipelineのworkerを`_REPO_ROOT`から起こすのと同じ）。
_PACKAGE_ROOT = Path(__file__).resolve().parents[2]

#: 中継がワーカーを起こし、結果を書いて終わるまでを待つ上限（秒）。普段は1秒もかからない。
_RELAY_TIMEOUT_S = 20.0


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

    class _STARTUPINFOEXW(ctypes.Structure):
        _fields_ = [
            ("StartupInfo", _STARTUPINFOW),
            ("lpAttributeList", ctypes.c_void_p),
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
    _EXTENDED_STARTUPINFO_PRESENT = 0x00080000
    _PROC_THREAD_ATTRIBUTE_HANDLE_LIST = 0x00020002

    _InitializeProcThreadAttributeList = _kernel32.InitializeProcThreadAttributeList
    _InitializeProcThreadAttributeList.argtypes = [ctypes.c_void_p, _wt.DWORD, _wt.DWORD,
                                                   ctypes.POINTER(ctypes.c_size_t)]
    _InitializeProcThreadAttributeList.restype = _wt.BOOL

    _UpdateProcThreadAttribute = _kernel32.UpdateProcThreadAttribute
    _UpdateProcThreadAttribute.argtypes = [ctypes.c_void_p, _wt.DWORD, ctypes.c_size_t,
                                           ctypes.c_void_p, ctypes.c_size_t,
                                           ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]
    _UpdateProcThreadAttribute.restype = _wt.BOOL

    _DeleteProcThreadAttributeList = _kernel32.DeleteProcThreadAttributeList
    _DeleteProcThreadAttributeList.argtypes = [ctypes.c_void_p]
    _DeleteProcThreadAttributeList.restype = None

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

    class _MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [
            ("dwLength", _wt.DWORD),
            ("dwMemoryLoad", _wt.DWORD),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    _GlobalMemoryStatusEx = _kernel32.GlobalMemoryStatusEx
    _GlobalMemoryStatusEx.argtypes = [ctypes.POINTER(_MEMORYSTATUSEX)]
    _GlobalMemoryStatusEx.restype = _wt.BOOL


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
    # `create_unicode_buffer`はUTF-16へ符号化して確保する（BMP外の1文字は
    # サロゲートペアの2要素）。`(c_wchar * len(block))`に1文字ずつ詰めると、
    # 絵文字などを含む変数が1つあるだけでTypeErrorになり起動できない。
    # 埋め込みのNULはそのまま保たれる。
    return ctypes.create_unicode_buffer(block)


def _create_process(command: list[str], cwd: Path, log_path: Path,
                    creation_flags: int,
                    env: dict | None = None) -> "_PROCESS_INFORMATION":
    """CreateProcessWを呼び、PROCESS_INFORMATIONを返す。

    stdio用に開いたハンドルは成功・失敗いずれの経路でも必ず閉じる。子は
    CreateProcessWの時点で自分用の複製を受け取っているので、親側の複製を
    残しておく理由がない（残すとログファイルが解放されない）。

    `env`がNoneなら親の環境を継承する。dictならその環境全体で起動する
    （`resolve_python_launch`がvenvを再現するために使う）。

    **子へ渡すハンドルはログとNULの2つだけ**（`PROC_THREAD_ATTRIBUTE_HANDLE_LIST`）。
    `bInheritHandles=TRUE`だけだと、親の継承可能なハンドルが全部渡る。MCPサーバの
    標準入出力はアプリとつながるパイプで継承可能なので、それが中継・worker・Consoleへ
    渡り、サーバが終わってもアプリの側にパイプの終わりが届かなくなる。
    """
    environment = None
    if env is not None:
        environment = _environment_block(env)
        creation_flags |= CREATE_UNICODE_ENVIRONMENT
    log_handle, null_handle = _open_child_stdio(log_path)
    try:
        siex = _STARTUPINFOEXW()
        siex.StartupInfo = _startup_info(log_handle, null_handle)
        siex.StartupInfo.cb = ctypes.sizeof(_STARTUPINFOEXW)
        inherited = (_wt.HANDLE * 2)(log_handle, null_handle)
        pi = _PROCESS_INFORMATION()
        cmdline = ctypes.create_unicode_buffer(subprocess.list2cmdline(command))
        with _handle_list_attribute(inherited) as attribute_list:
            siex.lpAttributeList = attribute_list
            created = _CreateProcessW(None, cmdline, None, None, True,
                                      creation_flags | _EXTENDED_STARTUPINFO_PRESENT,
                                      environment, str(cwd),
                                      ctypes.byref(siex.StartupInfo), ctypes.byref(pi))
            if not created:
                raise _win_error("PROCESS_LAUNCH_FAILED",
                                 "CreateProcessWに失敗した",
                                 {"command": command, "cwd": str(cwd)})
        return pi
    finally:
        _CloseHandle(log_handle)
        _CloseHandle(null_handle)


@contextlib.contextmanager
def _handle_list_attribute(handles):
    """`handles`だけを子へ渡す属性リストを作り、そのアドレスを渡す（抜けるときに消す）。

    `handles`（HANDLEの配列）は継承可能でなければならず、CreateProcessWが終わるまで
    呼び出し側が保持する。
    """
    size = ctypes.c_size_t(0)
    # 1回目は必要な大きさを聞くだけ（ERROR_INSUFFICIENT_BUFFERで失敗するのが正常）。
    _InitializeProcThreadAttributeList(None, 1, 0, ctypes.byref(size))
    buffer = ctypes.create_string_buffer(size.value)
    if not _InitializeProcThreadAttributeList(buffer, 1, 0, ctypes.byref(size)):
        raise _win_error("PROCESS_ATTRIBUTE_FAILED",
                         "InitializeProcThreadAttributeListに失敗した")
    try:
        if not _UpdateProcThreadAttribute(buffer, 0, _PROC_THREAD_ATTRIBUTE_HANDLE_LIST,
                                          ctypes.byref(handles), ctypes.sizeof(handles),
                                          None, None):
            raise _win_error("PROCESS_ATTRIBUTE_FAILED",
                             "UpdateProcThreadAttributeに失敗した")
        yield ctypes.addressof(buffer)
    finally:
        _DeleteProcThreadAttributeList(buffer)


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


def _in_any_job(process_handle) -> bool:
    """`process_handle` のプロセスが、何段目かを問わずどれかのJobに入っているか。

    `IsProcessInJob` にJobハンドル NULL を渡すと「いずれかのJobの下で動いているか」を答える
    （入れ子のJobの外側まで含む）。どのJobか・そのJobの制限は分からない。
    """
    in_job = _wt.BOOL()
    if not _IsProcessInJob(process_handle, None, ctypes.byref(in_job)):
        raise _win_error("JOB_QUERY_FAILED", "IsProcessInJobに失敗した")
    return bool(in_job.value)


def current_process_in_job() -> bool | None:
    """現在のプロセスがどれかのJobに入っているか（何段目でも）。Windows以外・問い合わせ失敗はNone。

    pipelineのworkerが起動直後に自分を確かめ、runの記録に残すのに使う
    （`launch_detached`の`worker_in_job`と同じ事実を、worker自身から見たもの）。
    """
    if not _IS_WINDOWS:
        return None
    try:
        return _in_any_job(_GetCurrentProcess())
    except DomainError:
        return None


def current_job_state() -> dict:
    """現在のプロセスのJobの状態: `{"in_job": bool | None, "nearest_job": str | None}`。

    `in_job`は`current_process_in_job`（何段目でも）。`nearest_job`はJobに入っているときの
    直近のJobの`_breakaway_mode`の値（`detach_breakaway_mode`）で、`inherited_job`なら直近の
    Jobは閉じても止めない（KILL_ON_JOB_CLOSEが無い）と分かる。それより外側の段は見えない。
    Jobの外・Windows以外・問い合わせ失敗の`nearest_job`はNone。pipelineのworkerが起動直後に
    自分を確かめるのに使う。
    """
    in_job = current_process_in_job()
    nearest_job = None
    if in_job:
        try:
            nearest_job = detach_breakaway_mode()
        except DomainError:
            nearest_job = None
    return {"in_job": in_job, "nearest_job": nearest_job}


def commit_memory_status() -> dict | None:
    """機械全体のcommit（物理メモリ＋ページファイル）の空きと上限をバイトで返す。

    `GlobalMemoryStatusEx` の `ullAvailPageFile`（まだ確保されていないcommitの量）と
    `ullTotalPageFile`（今のcommitの上限）。空きの物理メモリ（`ullAvailPhys`）だけを見ても、
    commitの上限で落ちる（Windowsのエラー1450）ことは予測できない——物理メモリが空いていても
    確保済みのcommitが上限に近ければ、大きな確保は失敗する。Windows以外・問い合わせ失敗はNone。
    """
    if not _IS_WINDOWS:
        return None
    status = _MEMORYSTATUSEX()
    status.dwLength = ctypes.sizeof(_MEMORYSTATUSEX)
    if not _GlobalMemoryStatusEx(ctypes.byref(status)):
        return None
    return {"available_commit_bytes": int(status.ullAvailPageFile),
            "commit_limit_bytes": int(status.ullTotalPageFile),
            "available_physical_bytes": int(status.ullAvailPhys)}


def _detach_flags(mode: str, command: list[str]) -> int:
    """切り離し起動のCreateProcessWフラグを返す（`mode`は`detach_breakaway_mode`の値）。

    `launch_detached`の2段（呼び出し側→中継、中継→ワーカー）で同じ規則を使う。
    breakawayのフラグは`breakaway_ok`のときだけ付ける（`silent_breakaway`は付けなくても
    Jobから外れる）。`unsupported`なら起こさずに`DETACH_UNSUPPORTED`を送出する。
    """
    if mode == "unsupported":
        raise DomainError(
            "DETACH_UNSUPPORTED",
            "親Job ObjectがKILL_ON_JOB_CLOSE付きでbreakawayも許可していないため、"
            "親の終了後もワーカーが生き残ることを保証できない",
            {"command": command, "breakaway": mode})
    flags = DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
    if mode == "breakaway_ok":
        flags |= CREATE_BREAKAWAY_FROM_JOB
    return flags


def _spawn_detached(command: list[str], cwd: Path, log_path: Path,
                    env: dict | None) -> dict:
    """現在のプロセスから1段だけ切り離して起こし、pid・identity・breakaway・worker_in_jobを返す（待たない）。

    中継（`metabolomix.core.detach_relay`）がワーカーを起こすときの本体。`command`は
    `resolve_python_launch`を通した後のもの。identityは子のプロセスハンドルを開いている
    間に取る（ハンドルを持つ間はpidが再利用されないので、pidと生成時刻が同じ子を指す）。

    `breakaway`はこのプロセス（中継）の直近のJobから見た値で、外側のJobは見えない
    （`QueryInformationJobObject(NULL)`は直近のJobしか答えない）。そこで起こした子を
    `IsProcessInJob(子, NULL)`で確かめ、何段目かを問わずJobに残っているかを`worker_in_job`に
    返す（問い合わせに失敗したらNone）。
    """
    mode = detach_breakaway_mode()
    flags = _detach_flags(mode, command)
    pi = _create_process(command, cwd, log_path, flags, env)
    try:
        identity = {"pid": int(pi.dwProcessId),
                    "creation_time": _creation_time(pi.hProcess)}
        try:
            worker_in_job = _in_any_job(pi.hProcess)
        except DomainError:
            worker_in_job = None
    finally:
        _CloseHandle(pi.hThread)
        _CloseHandle(pi.hProcess)
    return {"pid": int(pi.dwProcessId), "identity": identity, "breakaway": mode,
            "worker_in_job": worker_in_job}


def launch_detached(command: list[str], *, cwd: Path, log_path: Path) -> dict:
    """親から切り離してワーカーを起動し、pidとidentityを返す（ワーカーは待たない）。

    切り離しの目的は「MCPの1呼び出しが消えても解析が続く」こと。親Jobから
    抜けられず、かつそのJobが閉じたら殺される環境では目的を果たせないので、
    起動せずに`DETACH_UNSUPPORTED`を送出する。黙って起動すると、親の終了と
    同時に消えるワーカーを「起動できた」と記録してしまう。

    **ワーカーは中継のプロセスから起こす**。デスクトップアプリは呼び出しを打ち切ると
    MCPサーバを作り直す。試し走りの観察から、そのときサーバのプロセスの木ごと消して
    いると見ている（実測はしていない）。木ごとの終了は親子の記録をたどるので、サーバが
    直接起こしたワーカーは切り離しのフラグを付けていても巻き込まれる。
    そこで、このプロセスは中継（`python -m metabolomix.core.detach_relay`）を起こし、
    中継がワーカーを起こしてそのpidとidentityを書き、すぐ終わる。ワーカーの親（中継）は
    もう居ないので、木をたどってもワーカーに届かない。2段とも同じフラグの規則
    （`_detach_flags`。breakawayは許されるときだけ）で起こす。

    戻り値: `{"pid", "identity", "breakaway", "relay_breakaway", "worker_in_job", "log_path"}`。
    pidとidentityは**ワーカーのもの**（中継のものではない）。`breakaway`はこのプロセスの直近の
    Jobから見た値、`relay_breakaway`は中継の直近のJobから見た値（ワーカーを起こした段）。
    `worker_in_job`は起こしたワーカーが何段目かを問わずJobに残っているか（中継が
    `IsProcessInJob(ワーカー, NULL)`で確かめる。分からなければNone）。
    ワーカーの起動の誤り（`PROCESS_LAUNCH_FAILED`など）は中継から同じcodeで返す。

    限界: 直近のJobの制限しか読めないので、`DETACH_UNSUPPORTED`で断れるのは直近（この
    プロセスと中継のそれぞれ1段）がKILL_ON_JOB_CLOSE付きでbreakawayを許さないときだけ。
    それより外側のJobや、breakawayを許さないJob（`inherited_job`）には、中継もワーカーも
    残る。木ごとの終了には残るが、そのJobが閉じる（KILL_ON_JOB_CLOSE）・終了させられる
    （`TerminateJobObject`）と巻き込まれる。`worker_in_job`が真ならこの限界の中にいる。
    """
    if not _IS_WINDOWS:
        raise _unsupported_platform("launch_detached")
    command, cwd = _validate_launch_args(command, cwd)
    command, env = resolve_python_launch(command)
    log_path = Path(log_path)

    mode = detach_breakaway_mode()
    flags = _detach_flags(mode, command)
    launched = _launch_via_relay(command, cwd, log_path, env, flags)
    worker_in_job = launched.get("worker_in_job")
    return {"pid": launched["pid"], "identity": launched["identity"],
            "breakaway": mode, "relay_breakaway": launched.get("breakaway"),
            "worker_in_job": worker_in_job if isinstance(worker_in_job, bool) else None,
            "log_path": str(log_path)}


def _launch_via_relay(command: list[str], cwd: Path, log_path: Path,
                      env: dict | None, flags: int) -> dict:
    """中継を`flags`で起こし、終わるのを待って、中継が書いたワーカーの起動結果を返す。

    中継へ渡すのはワーカーの起動コマンド・cwd・ログ・venvの印（`__PYVENV_LAUNCHER__`の値）と
    結果を書く先だけで、環境変数そのものは渡さない。中継はこのプロセスの環境を受け継ぎ、
    解釈系は起動時に`__PYVENV_LAUNCHER__`を自分の環境から消すので、中継がそれを足し直す
    （`env`は`resolve_python_launch`が作る「このプロセスの環境＋その印」に限られる）。

    中継のstdout/stderrは`<ログ名>.relay<拡張子>`へ結ぶ。どの経路でも空なら消し（中継を
    起こせなかった・止めたときも）、中身があって起動結果が無ければ残して
    `DETACH_RELAY_FAILED`の`details`に末尾を載せる。
    上限（`_RELAY_TIMEOUT_S`）までに終わらなければ中継を止めて`DETACH_RELAY_TIMEOUT`
    （ワーカーは起きているかもしれない。pipelineならworker自身がidentityを記録する）。
    """
    token = os.urandom(6).hex()
    result_path = log_path.with_name(f"{log_path.stem}.relay-{token}.json")
    relay_log = log_path.with_name(f"{log_path.stem}.relay{log_path.suffix}")
    request = {"command": command, "cwd": str(cwd), "log_path": str(log_path),
               "launcher": (env or {}).get(_PYVENV_LAUNCHER_ENV),
               "result_path": str(result_path)}
    relay_command, relay_env = resolve_python_launch(
        [sys.executable, "-m", _RELAY_MODULE, json.dumps(request)])
    details = {"command": command, "relay_log": str(relay_log)}
    try:
        pi = _create_process(relay_command, _PACKAGE_ROOT, relay_log, flags, relay_env)
        _CloseHandle(pi.hThread)
        try:
            details["relay_pid"] = int(pi.dwProcessId)
            status = _WaitForSingleObject(pi.hProcess, int(_RELAY_TIMEOUT_S * 1000))
            if status != _WAIT_OBJECT_0:
                error = (DomainError("DETACH_RELAY_TIMEOUT",
                                     f"中継が{_RELAY_TIMEOUT_S:g}秒以内に終わらなかった"
                                     "（ワーカーは起きているかもしれない）", details)
                         if status == _WAIT_TIMEOUT else
                         _win_error("PROCESS_WAIT_FAILED", "中継の終わりを待てなかった",
                                    details))
                _TerminateProcess(pi.hProcess, JOB_TERMINATION_EXIT_CODE)
                raise error
            exit_code = _exit_code(pi.hProcess)
        finally:
            _CloseHandle(pi.hProcess)
        result = _read_relay_result(result_path)
    finally:
        with contextlib.suppress(OSError):
            result_path.unlink(missing_ok=True)
        _discard_if_empty(relay_log)

    error = result.get("error") if isinstance(result, dict) else None
    if isinstance(error, dict):
        raise DomainError(str(error.get("code") or "DETACH_RELAY_FAILED"),
                          str(error.get("message") or "中継がワーカーを起こせなかった"),
                          dict(error.get("details") or {}))
    if exit_code != 0 or not _is_launch_result(result):
        details.update(exit_code=exit_code, relay_log_tail=_tail(relay_log))
        raise DomainError("DETACH_RELAY_FAILED",
                          f"中継がワーカーの起動結果を返さなかった（終了コード {exit_code}）",
                          details)
    return result


def _discard_if_empty(path: Path) -> None:
    """空のファイルを消す（無い・消せないなら何もしない）。"""
    with contextlib.suppress(OSError):
        if path.stat().st_size == 0:
            path.unlink()


def _read_relay_result(path: Path):
    """中継が書いた結果のJSONを読む。無い・壊れているならNone。"""
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _is_launch_result(result) -> bool:
    """中継の結果が、ワーカーのpidとidentity（同じpidと生成時刻）の組か。"""
    if not isinstance(result, dict):
        return False
    pid = result.get("pid")
    identity = result.get("identity")
    return (isinstance(pid, int) and not isinstance(pid, bool) and pid > 0
            and isinstance(identity, dict) and identity.get("pid") == pid
            and isinstance(identity.get("creation_time"), int))


def _tail(path: Path, limit: int = 2000) -> str:
    """ログの末尾（読めなければ空）。"""
    try:
        return path.read_text(encoding="utf-8", errors="replace")[-limit:]
    except OSError:
        return ""


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


@contextlib.contextmanager
def wait_for_lock(path: Path, *, poll_interval: float = 1.0, on_wait=None,
                  should_abort=None):
    """`file_lock`と同じOSの排他ロックを、取れるまで（上限なしで）待って取る。

    長く待つことが前提の排他（機械全体でMS-DIAL Consoleを同時に1本にする等）に使う。
    待ち始めた（最初の試みで取れなかった）ときに`on_wait()`を1回だけ呼ぶ——待っていることを
    外へ知らせる場所。待つ間は`poll_interval`秒ごとに`should_abort()`を見て、真なら
    `LOCK_WAIT_ABORTED`で抜ける（取消の要求など）。ロックはwithを抜けた時と、持っている
    プロセスが終わった時（OSが外す）に解ける。yieldする値は待った秒数。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_RDWR | os.O_CREAT
    flags |= getattr(os, "O_BINARY", 0)
    fd = os.open(path, flags, 0o644)
    acquired = False
    started = time.monotonic()
    try:
        notified = False
        while not _try_lock(fd):
            if not notified:
                notified = True
                if on_wait is not None:
                    on_wait()
            if should_abort is not None and should_abort():
                raise DomainError("LOCK_WAIT_ABORTED",
                                  f"ロックを待つのを取りやめた: {path}", {"path": str(path)})
            time.sleep(max(0.0, poll_interval))
        acquired = True
        yield time.monotonic() - started
    finally:
        if acquired:
            _unlock(fd)
        os.close(fd)
