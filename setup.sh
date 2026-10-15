#!/bin/sh
# ms-data-parser の Python 環境（.venv）と MCP 設定を 1 回で揃える（mac / Linux / Git Bash）。
#
#   sh setup.sh          利用者
#   sh setup.sh --dev    開発者（pytest も入る）
#   sh setup.sh --no-uv  uv があっても標準の venv + pip で作る（検証用）
#
# uv があれば uv.lock どおりに入れる。無ければ Python 3.13 以上で venv を作り、
# requirements.lock.txt（uv.lock の書き出し）を pip で入れる。最後に
# scripts/setup_finish.py が import の自己確認と .mcp.json / lipidmix.local.toml の生成を行う。
# R（パスウェイ投影）は含まない: 使うときだけ Rscript r/setup.R。
set -e
cd "$(dirname "$0")"

dev=0
use_uv=1
for arg in "$@"; do
    case "$arg" in
        --dev) dev=1 ;;
        --no-uv) use_uv=0 ;;
        *) echo "[setup] 不明な引数: $arg" >&2; exit 2 ;;
    esac
done

fail() { echo "[setup] 失敗: $1" >&2; exit 1; }

venv_python() {
    if [ -x .venv/Scripts/python.exe ]; then echo .venv/Scripts/python.exe
    else echo .venv/bin/python; fi
}

if [ "$use_uv" = 1 ] && command -v uv >/dev/null 2>&1; then
    echo "[setup] uv で .venv を作ります（uv.lock どおり）。"
    if [ "$dev" = 1 ]; then group="--group dev"; else group="--no-dev"; fi
    # shellcheck disable=SC2086
    uv sync --locked $group || fail "uv sync が失敗しました。上の出力を確認してください。"
else
    python=""
    for candidate in python3.14 python3.13 python3 python; do
        if command -v "$candidate" >/dev/null 2>&1 &&
           "$candidate" -c "import sys; sys.exit(0 if sys.version_info >= (3, 13) else 1)" 2>/dev/null; then
            python="$candidate"; break
        fi
    done
    [ -n "$python" ] || fail "Python 3.13 以上も uv も見つかりません。どちらかを入れてから再実行してください。 uv: https://docs.astral.sh/uv/getting-started/installation/  Python: https://www.python.org/downloads/"
    echo "[setup] uv が無いので標準の venv で作ります（$python）。"
    "$python" -m venv .venv || fail "python -m venv が失敗しました。"
    if [ "$dev" = 1 ]; then lock=requirements-dev.lock.txt; else lock=requirements.lock.txt; fi
    "$(venv_python)" -m pip install --disable-pip-version-check -r "$lock" ||
        fail "pip install -r $lock が失敗しました。上の出力を確認してください。"
fi

"$(venv_python)" scripts/setup_finish.py || fail "セットアップの確認で問題が見つかりました（上の出力）。"
