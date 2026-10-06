"""MS-DIAL Console の実行体解決とコマンド組み立て。

exe_path は引数で注入できる（テスト用 fake の差し込みに使う）。
既定は `lipidmix.core.user_config` の `msdial.exe`（環境変数 MSDIAL_EXE →
`lipidmix.local.toml` の `[msdial] exe`）から取得し、決められなければ
`MsdialExeNotFoundError`（EnvironmentError）を上げる。

**Console を実際に起動して見張るのは `lipidmix.console.execution.supervise`**。
同期呼出しも切り離しワーカー（`lipidmix.console.worker`）も pipeline の上流工程も
そこを通り、成功・非ゼロ・timeout・取消・起動不能のどの終了経路でも終了証跡と
成果物を残す。

かつてここにあった `run_msdial` / `run_msdial_detached` は削除した。前者は非ゼロ
終了で例外を送出し、呼び出し側（`console_run` の同期分岐）がそこで即 return して
生成物の収集を丸ごと飛ばしていた——**失敗した実行が何の証拠も残さない**という、
`supervise` が解いた欠陥そのもの。後者は誰も監視しないまま Console を放流し、
終了コードも終了理由も回収できなかった。どちらも代わりを増やさずに消す。

`get_exe_path` と `build_msdial_cmd` は `supervise` と計画層の共有部分で、
Console のコマンドラインが 2 通りに分裂しないための唯一の組み立て場所。
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

from lipidmix.core import user_config

_EXE_WHAT = "MS-DIAL Console の実行体（MSDIALCUI.exe / MsdialConsoleApp.exe）"


class MsdialExeNotFoundError(EnvironmentError):
    """Console の実行体を決められない。

    `code` は封筒にそのまま載せる。未設定は `MSDIAL_EXE_NOT_FOUND`、設定ファイルが
    読めなければ `CONFIG_INVALID`。後者も同じ型で送るのは、呼び出し側の既存の
    `except EnvironmentError` をすべてそのまま生かすため。
    """

    def __init__(self, message: str, *, code: str = "MSDIAL_EXE_NOT_FOUND",
                 details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.details = details if details is not None else {}


def get_exe_path() -> str:
    """Console の実行体のパス（環境変数 MSDIAL_EXE → 設定ファイルの `[msdial] exe`）。"""
    try:
        setting = user_config.get_setting("msdial.exe")
    except user_config.ConfigInvalidError as exc:
        raise MsdialExeNotFoundError(exc.message, code=exc.code, details=exc.details()) from exc
    if setting is None:
        raise MsdialExeNotFoundError(
            "MS-DIAL Console の実行体が設定されていません。"
            + user_config.missing_hint("msdial.exe", _EXE_WHAT),
            details=user_config.describe_missing("msdial.exe"))
    if setting.source == "config_file" and not Path(setting.value).is_file():
        # 設定ファイルの値は利用者が書き写したパスなので、typo はここで出どころ付きで返す。
        # 環境変数の値は PATH 解決される素の名前でありうるので存在確認しない。
        raise MsdialExeNotFoundError(
            f"{user_config.setting_label(setting)} が指すファイルがありません: {setting.value}",
            details=user_config.describe_missing("msdial.exe", setting))
    return setting.value


def is_console_exe(
    exe_path: str,
    timeout_s: int = 15,
    *,
    raise_on_os_error: bool = False,
) -> bool:
    """MSDIAL_EXE が Console 実行体かを --help の出力で判定する。

    Console は旧ビルドも master ビルドも --help（旧は引数エラー時の usage）に
    サブコマンド名 `lcms` を含む。GUI の MSDIAL.exe はコンソール出力を持たず、
    ウィンドウを開いたまま返らないため、ここで実行前に弾く。既定では
    起動不能も False にするが、実行経路は `raise_on_os_error=True` で
    起動不能を既存の NOT_FOUND 経路へ渡せる。
    """
    try:
        completed = subprocess.run(
            [exe_path, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            timeout=timeout_s,
            text=True,
            errors="replace",
        )
    except subprocess.TimeoutExpired:
        return False
    except OSError:
        if raise_on_os_error:
            raise
        return False
    return "lcms" in (completed.stdout or "")


def build_msdial_cmd(
    exe: str,
    dataset_root: Path,
    msdial_out_dir: Path,
    method_file: Path,
    save_project: bool = False,
) -> list[str]:
    """Console のコマンドラインを組む。

    待って実行する経路・切り離して実行する経路・`execution.supervise` の
    いずれでも**同じ引数**でなければならないので、組み立てはここ 1 か所に置く。
    """
    cmd = [exe, "lcms", "-i", str(dataset_root), "-o", str(msdial_out_dir),
           "-m", str(method_file)]
    if save_project:
        cmd.append("-p")
    return cmd


def is_process_running(pid: int) -> bool:
    """pid が生きているかを返す（判定できないときは False ＝安全側）。

    **旧 `.detached-state.json` を持つジョブを判定するためだけに残している
    互換 wrapper**。新しい経路は `lipidmix.core.process_control.same_process` を
    使う——pid だけでは pid 再利用を見分けられず、別のプロセスを「同じ実行が
    まだ生きている」と誤読するため。旧 state には creation_time が無いので、
    そこだけはこの弱い判定しか使えない。
    """
    if not pid or pid < 0:
        return False
    if os.name == "nt":
        import subprocess as _sp
        try:
            out = _sp.run(["tasklist", "/FI", f"PID eq {int(pid)}", "/NH"],
                          stdout=_sp.PIPE, stderr=_sp.DEVNULL, stdin=_sp.DEVNULL,
                          timeout=15, text=True, errors="replace")
        except (OSError, _sp.TimeoutExpired):
            return False
        return str(int(pid)) in (out.stdout or "")
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # 自分のものではないが存在はする
    except OSError:
        return False
    return True


#: Console 実行体の名前。GUI（MSDIAL.exe）は候補にしない —— 引数を解釈せず
#: ウィンドウを開いたまま返らないので、設定されると実行時に固まる。
_CONSOLE_EXE_NAMES = ("msdialcui.exe", "msdialconsoleapp.exe")

#: 候補探索の既定の起点。MS-DIAL は zip を展開しただけで使うことが多く、
#: ホーム直下かデスクトップに置かれているのが実情。
def _default_search_roots() -> list[Path]:
    home = Path.home()
    roots = [home, home / "Desktop", home / "Downloads",
             Path("C:/Program Files"), Path("C:/Program Files (x86)")]
    return [r for r in roots if r.is_dir()]


def msdial_exe_candidates(search_roots=None, max_depth: int = 3) -> list[str]:
    """MS-DIAL Console 実行体の候補を返す（設定はしない）。

    MSDIAL_EXE を設定できるのは人間だけなので、サーバにできるのは「どこに
    ありそうか」を示すところまで。深さを制限するのは、ホーム配下を無制限に
    歩くと計画が数十秒止まるため。
    """
    roots = [Path(r) for r in (search_roots if search_roots is not None
                               else _default_search_roots())]
    found: list[str] = []
    for root in roots:
        if not root.is_dir():
            continue
        base_depth = len(root.parts)
        for dirpath, dirnames, filenames in os.walk(root):
            depth = len(Path(dirpath).parts) - base_depth
            if depth >= max_depth:
                dirnames[:] = []
            for name in filenames:
                if name.lower() in _CONSOLE_EXE_NAMES:
                    found.append(str(Path(dirpath) / name))
    return found
