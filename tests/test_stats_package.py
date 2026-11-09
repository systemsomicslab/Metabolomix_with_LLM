"""統計エンジン・パッケージの公開面と隔離規則。"""

import ast
import importlib
from pathlib import Path

import metabolomix.analysis.stats as stats

STATS_DIR = Path(stats.__file__).parent
FORBIDDEN_EXACT = {"metabolomix.core.session_state", "metabolomix.core.mcp_core"}
FORBIDDEN_PREFIX = "metabolomix.tools"


def test_all_names_are_importable_submodules():
    for name in stats.__all__:
        module = importlib.import_module(f"metabolomix.analysis.stats.{name}")
        assert module.__name__ == f"metabolomix.analysis.stats.{name}"


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            found.add(node.module)
    return found


def test_no_stats_module_imports_session_mcp_or_tools():
    offenders = []
    for path in sorted(STATS_DIR.glob("*.py")):
        for module in _imported_modules(path):
            if module in FORBIDDEN_EXACT or module == FORBIDDEN_PREFIX or module.startswith(FORBIDDEN_PREFIX + "."):
                offenders.append(f"{path.name}: {module}")
    assert offenders == []
