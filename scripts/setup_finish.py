"""setup.ps1 / setup.sh の共通の後半。`.venv` の Python で走る。

    .venv の python scripts/setup_finish.py [--root <リポジトリ>]

やること（spec 2026-10-15 §3）:
1. pyproject.toml の直接依存がすべて import できるかを確かめる（できなければ 1 で終わる）。
2. lipidmix.local.toml が無ければ雛形から複写する。
3. .mcp.json が無ければ .venv の Python を指して書く。あれば**触らず**、違いだけを言う。
4. Claude Desktop 用の貼り付け断片と、次の手順を表示する。

既存の設定ファイルを書き換えないのは、利用者が env（参照ライブラリのパス等）を
書き込んでいる可能性があり、上書きすると復元できないため。
"""
from __future__ import annotations

import argparse
import importlib
import json
import re
import shutil
import sys
import tomllib
from pathlib import Path

SERVER_NAME = "ms-data-parser"

#: 配布名（正規化済み）と import 名が違うものだけを載せる。
IMPORT_NAMES = {
    "scikit-learn": "sklearn",
    "pymupdf": "fitz",
}


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def direct_dependencies(root: Path) -> list[str]:
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    return [_norm(re.split(r"[<>=!~;\[ ]", spec, maxsplit=1)[0])
            for spec in project["project"]["dependencies"]]


def missing_imports(distributions: list[str], importer=importlib.import_module) -> list[str]:
    missing = []
    for dist in distributions:
        try:
            importer(IMPORT_NAMES.get(dist, dist.replace("-", "_")))
        except ImportError:
            missing.append(dist)
    return missing


def venv_python(root: Path) -> Path:
    venv = root.resolve() / ".venv"
    if sys.platform == "win32":
        return venv / "Scripts" / "python.exe"
    return venv / "bin" / "python"


def server_entry(root: Path) -> dict:
    return {
        "type": "stdio",
        "command": venv_python(root).as_posix(),
        "args": ["-u", (root.resolve() / "server.py").as_posix()],
        "env": {"PYTHONIOENCODING": "utf-8"},
    }


def ensure_local_config(root: Path) -> str:
    target = root / "lipidmix.local.toml"
    if target.exists():
        return "kept"
    example = root / "lipidmix.example.toml"
    if not example.is_file():
        return "no_example"
    shutil.copyfile(example, target)
    return "created"


def ensure_mcp_json(root: Path) -> tuple[str, str]:
    target = root / ".mcp.json"
    wanted = server_entry(root)
    if not target.exists():
        target.write_text(json.dumps({"mcpServers": {SERVER_NAME: wanted}}, indent=2) + "\n",
                          encoding="utf-8")
        return "created", f"{target} を作りました。"
    try:
        current = json.loads(target.read_text(encoding="utf-8"))["mcpServers"][SERVER_NAME]
        current.get("command")  # dict でない項（文字列など）はここで AttributeError にする
    except (ValueError, KeyError, TypeError, AttributeError):
        return "differs", (f"{target} は既にありますが {SERVER_NAME} の項がありません（触っていません）。"
                           f"次の項を mcpServers に足してください:\n"
                           f"{json.dumps({SERVER_NAME: wanted}, indent=2)}")
    if current.get("command") == wanted["command"]:
        return "kept", f"{target} は既に .venv の Python を指しています。"
    return "differs", (f"{target} は既にあります（触っていません）。{SERVER_NAME} の command を\n"
                       f"  現在: {current.get('command')}\n"
                       f"  推奨: {wanted['command']}\n"
                       f"に手で直してください（env などほかの項はそのままで構いません）。")


def desktop_snippet(root: Path) -> str:
    return json.dumps({"mcpServers": {SERVER_NAME: server_entry(root)}}, indent=2)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    root = parser.parse_args(argv).root.resolve()

    missing = missing_imports(direct_dependencies(root))
    if missing:
        print("[setup] import できない依存があります: " + ", ".join(missing))
        print("[setup] setup をもう一度実行してください。直らなければ出力を添えて報告してください。")
        return 1
    print("[setup] 依存の import を確認しました。")

    config = ensure_local_config(root)
    if config == "created":
        print(f"[setup] {root / 'lipidmix.local.toml'} を作りました。MS-DIAL や参照ライブラリを"
              "使うときは、ここに場所を書いてください（任意）。")
    elif config == "kept":
        print("[setup] lipidmix.local.toml は既にあるので触っていません。")

    _, message = ensure_mcp_json(root)
    print(f"[setup] {message}")

    print("\n[setup] Claude Desktop で使うときは、claude_desktop_config.json の mcpServers に"
          "次を足してください\n（同じ名前のサーバを .mcp.json と二重に登録しないでください）:")
    print(desktop_snippet(root))
    print("\n[setup] 完了。MCP クライアントを再起動してください。")
    print("[setup] パスウェイ投影（pathway_*）を使うときだけ、R を入れて `Rscript r/setup.R` を 1 度実行します。")
    print("[setup] 更新は会話で server_update を呼ぶか、git pull の後に setup をもう一度実行します。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
