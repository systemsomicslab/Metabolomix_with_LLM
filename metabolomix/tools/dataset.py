"""データセット入口ツール: list_data_files, load_dataset。

load_dataset は arf2 概観 → arf PCA を一括実行するオーケストレータ。arf/arf2 ツールを
呼ぶため metabolomix.arf.tools / metabolomix.arf2.tools に依存する（依存が最も深く、facade では最後に読み込む）。
DATA_DIR の差し替えは正準の mcp_core.DATA_DIR に対して行う。
"""
from pathlib import Path

from metabolomix.core import mcp_core
from metabolomix.core import path_resolvers
from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.path_resolvers import (
    resolve_arf_file_path,
    resolve_arf2_file_path,
    _describe_batch_selection,
)
from metabolomix.arf.tools import arf_parser
from metabolomix.arf2.tools import arf2_parser

__all__ = ["list_data_files", "load_dataset"]

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
    フォルダ全体をまず解析するなら `load_dataset(directory)` が入口。
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


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def load_dataset(directory: str | None = None) -> str:
    """データフォルダを指定して、最初の標準解析（arf2 概観 → arf 詳細）を一括実行します。

    MS-DIAL出力フォルダを解析する際の **入口** です。フォルダのパスを渡すと:
    1. `.arf2`（データセット全体のカタログ＝概観）を要約し、
    2. サンプル別強度を持つ `.arf` で PCA を実行します。フォルダに DriftSpots.arf と
       PeakProperties.arf が併存する場合は、解析に使う **PeakProperties.arf を自動選択** します。
    以降の `arf_list_classes` / `arf_parser` 等はこのセッション状態をそのまま利用できます。

    複数日付（複数回のMS-DIAL処理＝複数バッチ）のファイルが混在していても解析は
    止まりません。ファイル名の `AlignmentResult_<timestamp>` を見て **最新バッチを
    自動選択** し、選択結果を出力に明示します（旧バッチはスキップ）。

    - directory: MS-DIAL出力フォルダのパス。省略時は既定のデータディレクトリ
      (環境変数 LIPIDMIX_DATA_DIR または <project>/data) を使用します。
      明示した場合は以降のツールの既定探索先もこのフォルダに更新されます。
    """
    if directory:
        target_dir = Path(directory).expanduser()
        if not target_dir.exists():
            return f"データディレクトリが存在しません: {target_dir}"
        if not target_dir.is_dir():
            return f"指定されたパスはディレクトリではありません: {target_dir}"
        mcp_core.DATA_DIR = target_dir  # 以降のツールの既定探索先を更新（正準は mcp_core 側）

    arf2_path = resolve_arf2_file_path()
    arf_path = resolve_arf_file_path()

    blocks: list[str] = [
        f"## 📂 データセット読み込み: {mcp_core.DATA_DIR}\n"
        "標準の初期解析として **arf2（全体概観）→ arf（PeakProperties, サンプル別PCA）** を実行します。\n"
        "この出力（群構造・脂質クラス・極性など）は、解釈に進む前の『実験目的の推測とユーザー確認』"
        "（GATEWAY手順1）の材料になります。"
    ]

    # 意味論ダイジェストは入口の先頭で1回だけ前置する。ここで発火させておくと、
    # 後段で呼ぶ arf2_parser / arf_parser 内の同ガードは caveat_emitted により
    # no-op になり、ダイジェストが arf2 ブロック内へ埋没するのを防げる。
    blocks[0] = session_state.session.maybe_prepend_caveat(blocks[0])

    # 配布先のクローンが origin/main より遅れているときだけ 1 行を前置する。
    # 入口ツールは Markdown を返すので、鍵ではなくブロックとして差す。
    from metabolomix.core import version
    update = version.update_status()
    if update:
        blocks.append(f"🔄 **{update['message']}**")

    batch_note = _describe_batch_selection(mcp_core.DATA_DIR)
    if batch_note:
        blocks.append(batch_note)

    if arf2_path:
        blocks.append(arf2_parser(file_path=arf2_path))
    else:
        blocks.append("⚠️ .arf2 ファイルが見つかりませんでした（全体概観をスキップ）。")

    if arf_path:
        blocks.append(arf_parser(file_path=arf_path))
    else:
        blocks.append(
            "⚠️ 解析対象の .arf（PeakProperties.arf 等）が見つかりませんでした。"
        )

    return "\n\n".join(blocks)
