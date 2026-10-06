#!/bin/sh
# エージェント 1 体につき 1 つの隔離作業ツリーを作る。
#
#   使い方: sh scripts/new-worktree.sh <branch>
#
# なぜ要るか: 全エージェントが同じ作業ツリーを共有すると、片方の未コミット変更が
# もう片方のテスト結果に混ざり、HEAD も断りなく動く。実測で、同じスイートの結果が
# 一度の調査中に 774 → 779 → 782 → 805 とぶれ、失敗件数も 5 → 1 → 0 と動いた。
#
# 追跡外のものは worktree に複製されない（data/ analyses/ .mcp.json docs/HISTRY.md
# docs/task.md）。テストはこれらに依存しない規約なので新品の worktree でも全数通るが、
# MCP サーバとしての実行には効くので、ここで .mcp.json を生成して埋める。
set -e

branch="$1"
if [ -z "$branch" ]; then
    echo "使い方: sh scripts/new-worktree.sh <branch>" >&2
    exit 1
fi

# Windows 形式（C:/...）へ。JSON に書くので MSYS 形式では使えない。
win() { cygpath -m "$1"; }

main=$(cd "$(dirname "$(git rev-parse --git-common-dir)")" && pwd)
slug=$(printf '%s' "$branch" | tr '/' '-')
path="$main/.worktrees/$slug"

if [ -e "$path" ]; then
    echo "既に存在する: $path" >&2
    exit 1
fi

if git -C "$main" show-ref --verify --quiet "refs/heads/$branch"; then
    git -C "$main" worktree add "$path" "$branch"
else
    git -C "$main" worktree add "$path" -b "$branch"
fi

main_w=$(win "$main")
wt_w=$(win "$path")
# 外部資産の場所の設定（lipidmix.local.toml、追跡外）は worktree に来ない。
# main ツリーにあれば LIPIDMIX_CONFIG でそれを指す（設定を分裂させない）。
config_env=""
if [ -f "$main/lipidmix.local.toml" ]; then
    config_env="\"LIPIDMIX_CONFIG\": \"$main_w/lipidmix.local.toml\","
fi

# server.py はこの worktree のものを指す。main を指すと、worktree のコードを
# 編集しながら main の実装を試すことになり、最も気づきにくい形で嘘をつく。
# 一方で蓄積状態は main ツリーへ向ける。worktree ごとに分裂させない。
cat > "$path/.mcp.json" <<JSON
{
  "mcpServers": {
    "ms-data-parser": {
      "type": "stdio",
      "command": "C:/Python314/python.exe",
      "args": [
        "$wt_w/server.py"
      ],
      "env": {
        $config_env
        "LIPIDMIX_DATA_DIR": "$main_w/data",
        "LIPIDMIX_ANALYSES_DIR": "$main_w/analyses",
        "LIPIDMIX_KNOWLEDGE_DIR": "$main_w/knowledge",
        "LIPIDMIX_REPORTS_DIR": "$main_w/reports"
      }
    }
  }
}
JSON

echo "worktree : $path"
echo "branch   : $branch"
echo ".mcp.json: 生成済み（server.py=この worktree / 蓄積状態=main ツリー）"
echo
echo "記録は main ツリー側へ書く（worktree には存在しない。task.md は並列中なので末尾への追記だけ）:"
echo "  $main_w/docs/HISTRY.md"
echo "  $main_w/docs/task.md"
echo
echo "片付け: git worktree remove .worktrees/$slug"
