"""サーバ自身の保守ツール（更新の適用）。

`server_update` 1 本だけ。実体は `metabolomix/core/version.py` の `apply_update`
で、ここはそれを MCP 公開面に薄く載せるだけ。deps: mcp_core / version /
serialization。tools_* / server は import しない。

**なぜ prompt ではなくツールなのか**: Claude Desktop はローカル stdio サーバの
prompt を添付できない回帰を抱えている（anthropics/claude-code#82045・2026-07-28
報告、2026-09-24 時点で OPEN）。tools は影響を受けないので、`/update-...` の
スラッシュコマンドではなくツールとして出し、導線は更新通知の本文に置く
（`version._compute_update_status` がツール名を書く）。
"""
from mcp.types import ToolAnnotations

from metabolomix.core import version
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload

__all__ = ["server_update"]


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True,
    openWorldHint=True),
    structured_output=False)
def server_update() -> str:
    """このMCPサーバ自身を origin/main へ更新する（git 早送り＋必要なら依存の再導入）。

    更新通知（`load_dataset` の冒頭行、`dataset_status` / `console_status` の
    `update_available`）が出たときに呼ぶ。プロセスは触らないので、**適用後は
    利用者が Claude Desktop を再起動する必要がある**——戻り値の `message` を
    そのまま利用者へ伝えること。

    手元の未コミット変更がある・main 以外にいる・早送りできない・オフラインの
    いずれでも何もせずに `status="refused"` と `reason` を返す（例外は投げない）。
    """
    result = dict(version.apply_update())
    result["running_version"] = version.server_version()
    return json_payload(result)
