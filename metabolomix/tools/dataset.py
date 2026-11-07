"""データセット入口ツール: list_data_files。

データセットを開く入口 `load_dataset` は `metabolomix.tools.handle_tools`。
DATA_DIR の差し替えは正準の mcp_core.DATA_DIR に対して行う。
"""
from pathlib import Path

from metabolomix.core import mcp_core
from metabolomix.core import path_resolvers
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp

__all__ = ["list_data_files"]

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def list_data_files(
    extension: str | None = None,
    directory: str | None = None,
    all_files: bool = False,
) -> str:
    """データフォルダにある解析対象ファイルを拡張子別に一覧します。

    - directory: 探索するフォルダ。省略時は既定のデータディレクトリ
      (環境変数 LIPIDMIX_DATA_DIR または <project>/data)。
    - extension: 指定するとその拡張子だけに絞る（例 '.pai2', '.arf2', '.EIC.aef'）。
    - all_files: 既定では解析できる拡張子
      (.arf / .arf2 / .pai2 / .dcl / .EIC.aef / .mddata / .mdproject) だけを返す。
      測定生データ（MS-DIAL が読む .abf / .ibf / .cdf / .mzml / .wiff / .raw /
      .d / .wiff2 / .qgd / .lcd / .lrp / .imzml）まで見たいときだけ True にする。
      Agilent・Bruker の .d と Waters の .raw は**フォルダ**が 1 検体なので、
      末尾に `/` を付けて示す。

    返すのは絶対パスで、`pai2_parser(file_path=...)` 等へそのまま渡せる。
    フォルダ全体をまず解析するなら `load_dataset(directory=...)` が入口。
    """
    target_dir = Path(directory).expanduser() if directory else mcp_core.DATA_DIR
    if not target_dir.exists():
        return f"データディレクトリが存在しません: {target_dir}"
    if not target_dir.is_dir():
        return f"指定されたパスはディレクトリではありません: {target_dir}"

    paths = path_resolvers.list_data_files(
        extension=extension, directory=directory, all_files=all_files)
    if not paths:
        scope = extension or ("全ファイル" if all_files else "解析対象の拡張子")
        return (
            f"条件に一致するファイルがありません（ディレクトリ: {target_dir}, 対象: {scope}）。"
            "all_files=True で解析対象外のファイルも確認できます。"
        )

    # 拡張子ごとにまとめる。同じフォルダのパスが延々と並ぶより、何が何件あるかが
    # 先に見えるほうが次の一手（load_dataset か個別パーサか）を選びやすい。
    groups: dict[str, list[str]] = {}
    for path in sorted(paths):
        suffix = next(
            (s for s in path_resolvers.ANALYSABLE_EXTENSIONS if path.endswith(s)),
            Path(path).suffix or "(拡張子なし)",
        )
        groups.setdefault(suffix, []).append(path)

    lines = [f"データディレクトリ: {target_dir}", f"該当ファイル: {len(paths)} 件"]
    for suffix, members in sorted(groups.items()):
        lines.append(f"\n## {suffix}（{len(members)} 件）")
        # フォルダ形式の計測データ（.d / Waters の .raw）は末尾の `/` で示す。
        # 見分けが付かないと「開けないファイル」として扱われる。
        lines.extend(m + "/" if Path(m).is_dir() else m for m in members)
    return "\n".join(lines)
