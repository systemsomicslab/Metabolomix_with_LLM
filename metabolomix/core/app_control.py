"""Windows の Application Control が未署名アセンブリを阻むかの事前判定。

stdlib だけに依存する leaf。`metabolomix.core.mcp_core` / `metabolomix.<形式>.tools` /
`metabolomix.tools.*` を import してはいけない。

## なぜ要るか

MS-DIAL Console にプロジェクト保存（`-p`）を指示すると、**全検体の解析と
アライメントを終えた後**に `MsdialIntegrate.Parser.MsdialIntegrateSerializer` の
初期化で `MsdialLcImMsApi` を読む。Windows 11 の Smart App Control が有効だと、
ローカルビルドの未署名アセンブリのここが弾かれる:

    FileLoadException: ... アプリケーションからポリシーによってこのファイルが
    ブロックされました。(HRESULT からの例外:0x800711C7)

`0x800711C7` の下位 16bit は 4551 =「An Application Control policy has blocked
this file」。実測では 41 検体で **13 分半を費やしてから**落ちた。最悪の失敗の
仕方なので、Console を起動する前に判定する。

## 判定規則

**2 条件が揃ったときだけ「塞がれる」と断定する。**

1. Smart App Control が Enforce
2. プロジェクト保存に要るアセンブリが**未署名**

ポリシーの有無だけで止めると、署名済みの公式配布版を使う正当な構成まで
塞いでしまう。Smart App Control が弾くのはまさに未署名バイナリなので、
署名の有無を判定材料に入れることで誤検出が消える。

判定は「ポリシー状態」と「PE のバイト列」の 2 入力に対する純関数にしてあり、
実機の状態に依存せず試験できる。
"""
from __future__ import annotations

import struct
from pathlib import Path

__all__ = [
    "SMART_APP_CONTROL_ENFORCING",
    "has_authenticode_signature",
    "project_save_blocked",
    "smart_app_control_state",
]

#: `VerifiedAndReputablePolicyState` の Enforce。0=Off / 1=Enforce / 2=評価。
SMART_APP_CONTROL_ENFORCING = 1

#: 「呼び出し側が policy_state を渡していない」ことの印。`None` を既定値にすると
#: 「実機から読む」と「読めなかった（＝不明）」が同じ値になり、不明を Enforce と
#: 取り違える余地が残る。
_READ_FROM_MACHINE = object()

_POLICY_KEY = r"SYSTEM\CurrentControlSet\Control\CI\Policy"
_POLICY_VALUE = "VerifiedAndReputablePolicyState"

#: PE のデータディレクトリのうち Certificate Table（Authenticode の埋め込み先）。
_CERTIFICATE_DIRECTORY_INDEX = 4
_PE32_PLUS_MAGIC = 0x20B


def smart_app_control_state() -> int | None:
    """Smart App Control の状態を読む。読めなければ None。

    Windows 以外・レジストリ値が無い（未構成）・読めない、はすべて None。
    **None を Enforce として扱わない**——判定不能を「塞がれている」側へ倒すと、
    関係のない環境で正当な実行を止めてしまう。
    """
    try:
        import winreg
    except ImportError:
        return None
    try:
        with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, _POLICY_KEY) as key:
            value, _kind = winreg.QueryValueEx(key, _POLICY_VALUE)
    except OSError:
        return None
    return value if isinstance(value, int) else None


def has_authenticode_signature(path) -> bool | None:
    """PE に Authenticode 署名が埋め込まれているかを返す。判定不能なら None。

    データディレクトリの Certificate Table（index 4）の size が非ゼロかを見る。
    カタログ署名（ファイル外）はここに現れないので False になりうるが、
    Smart App Control が問題にするローカルビルドの成果物は埋め込み署名も
    カタログ署名も持たないため、この判定で足りる。

    **読めないものを「署名あり」とも「未署名」とも言わない**（None）。
    """
    path = Path(path)
    try:
        data = path.read_bytes()
    except OSError:
        return None
    try:
        if data[:2] != b"MZ":
            return None
        pe_off = struct.unpack_from("<I", data, 0x3C)[0]
        if data[pe_off:pe_off + 4] != b"PE\x00\x00":
            return None
        magic = struct.unpack_from("<H", data, pe_off + 0x18)[0]
        # optional header の後ろにデータディレクトリが並ぶ。PE32 と PE32+ で
        # 先頭までの距離が違う（取り違えると別の値を署名だと読む）。
        directories = pe_off + 0x18 + (0x70 if magic == _PE32_PLUS_MAGIC else 0x60)
        offset = directories + _CERTIFICATE_DIRECTORY_INDEX * 8
        _rva, size = struct.unpack_from("<II", data, offset)
    except (struct.error, IndexError):
        return None
    return size > 0


def project_save_blocked(assembly_path, *, policy_state=_READ_FROM_MACHINE) -> bool:
    """プロジェクト保存が Application Control に塞がれると断定できるかを返す。

    `policy_state` を**省略**すると実機から読む。試験では注入する。明示的な
    `None` は「状態が分からない」の意味で、塞がれているとは断定しない。

    **断定できるときだけ True。** ポリシーが Enforce でない、署名がある、
    判定材料が無い（アセンブリを読めない）はすべて False——「塞がれるかも
    しれない」で止めると、正当な構成まで動かせなくなる。
    """
    state = (smart_app_control_state() if policy_state is _READ_FROM_MACHINE
             else policy_state)
    if state != SMART_APP_CONTROL_ENFORCING:
        return False
    return has_authenticode_signature(assembly_path) is False
