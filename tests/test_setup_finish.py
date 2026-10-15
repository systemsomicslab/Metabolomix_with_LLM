"""setup.ps1 / setup.sh の共通の後半（scripts/setup_finish.py）を固定する。

縛るのは「利用者の手元を壊さない」こと: 既存の .mcp.json と lipidmix.local.toml は
バイト単位で触らず、違いがあれば言うだけにする（spec 2026-10-15 §3）。
"""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("setup_finish", ROOT / "scripts" / "setup_finish.py")
setup_finish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(setup_finish)


@pytest.fixture
def repo(tmp_path):
    (tmp_path / "pyproject.toml").write_text(
        '[project]\nname = "x"\ndependencies = ["mcp<2", "scikit-learn", "PyMuPDF"]\n',
        encoding="utf-8")
    (tmp_path / "lipidmix.example.toml").write_text("[msdial]\n", encoding="utf-8")
    (tmp_path / "server.py").write_text("", encoding="utf-8")
    return tmp_path.resolve()


def test_direct_dependencies_strip_specifiers(repo):
    assert setup_finish.direct_dependencies(repo) == ["mcp", "scikit-learn", "pymupdf"]


def test_missing_imports_maps_distribution_to_module_names():
    seen = []

    def importer(name):
        seen.append(name)
        if name == "fitz":
            raise ImportError(name)

    missing = setup_finish.missing_imports(["scikit-learn", "pymupdf", "mcp"], importer)
    assert seen == ["sklearn", "fitz", "mcp"]
    assert missing == ["pymupdf"]


def test_venv_python_is_absolute_and_inside_dot_venv(repo):
    python = setup_finish.venv_python(repo)
    assert python.is_absolute()
    assert python.parent.parent == repo / ".venv"


def test_creates_local_config_from_example(repo):
    assert setup_finish.ensure_local_config(repo) == "created"
    assert (repo / "lipidmix.local.toml").read_text(encoding="utf-8") == "[msdial]\n"


def test_keeps_existing_local_config_byte_for_byte(repo):
    (repo / "lipidmix.local.toml").write_bytes(b"[library]\nmsp_positive = 'x'\n")
    assert setup_finish.ensure_local_config(repo) == "kept"
    assert (repo / "lipidmix.local.toml").read_bytes() == b"[library]\nmsp_positive = 'x'\n"


def test_creates_mcp_json_pointing_at_venv_python(repo):
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "created"
    entry = json.loads((repo / ".mcp.json").read_text(encoding="utf-8"))["mcpServers"]["ms-data-parser"]
    assert entry["command"] == setup_finish.venv_python(repo).as_posix()
    assert entry["args"] == ["-u", (repo / "server.py").as_posix()]
    assert entry["env"] == {"PYTHONIOENCODING": "utf-8"}


def test_keeps_matching_mcp_json(repo):
    setup_finish.ensure_mcp_json(repo)
    before = (repo / ".mcp.json").read_bytes()
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "kept"
    assert (repo / ".mcp.json").read_bytes() == before


def test_reports_but_never_rewrites_a_differing_mcp_json(repo):
    original = json.dumps({"mcpServers": {"ms-data-parser": {
        "command": "C:/Python314/python.exe", "args": ["server.py"],
        "env": {"MSDIAL_MSP_POS": "x"}}}}).encode()
    (repo / ".mcp.json").write_bytes(original)
    state, message = setup_finish.ensure_mcp_json(repo)
    assert state == "differs"
    assert "C:/Python314/python.exe" in message
    assert setup_finish.venv_python(repo).as_posix() in message
    assert (repo / ".mcp.json").read_bytes() == original


def test_mcp_json_without_our_server_is_reported_not_rewritten(repo):
    original = b'{"mcpServers": {"other": {"command": "node"}}}'
    (repo / ".mcp.json").write_bytes(original)
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "differs"
    assert (repo / ".mcp.json").read_bytes() == original


def test_mcp_json_with_non_dict_server_entry_is_reported_not_rewritten(repo):
    original = b'{"mcpServers": {"ms-data-parser": "oops"}}'
    (repo / ".mcp.json").write_bytes(original)
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "differs"
    assert (repo / ".mcp.json").read_bytes() == original


def test_desktop_snippet_is_valid_json_with_absolute_paths(repo):
    snippet = json.loads(setup_finish.desktop_snippet(repo))
    entry = snippet["mcpServers"]["ms-data-parser"]
    assert Path(entry["command"]).is_absolute()
    assert Path(entry["args"][-1]).is_absolute()


def test_main_fails_when_a_dependency_is_missing(repo, monkeypatch, capsys):
    monkeypatch.setattr(setup_finish, "missing_imports", lambda names, importer=None: ["rdkit"])
    assert setup_finish.main(["--root", str(repo)]) == 1
    assert "rdkit" in capsys.readouterr().out
    assert not (repo / ".mcp.json").exists()


def test_main_succeeds_and_writes_both_files(repo, monkeypatch):
    monkeypatch.setattr(setup_finish, "missing_imports", lambda names, importer=None: [])
    assert setup_finish.main(["--root", str(repo)]) == 0
    assert (repo / ".mcp.json").is_file()
    assert (repo / "lipidmix.local.toml").is_file()
