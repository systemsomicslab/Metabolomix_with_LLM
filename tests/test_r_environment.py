"""R の仮想環境（renv）の定義が追跡されていることを縛る。R が無くても走る。"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R_DIR = ROOT / "r"


def test_renv_files_are_present():
    for rel in ("renv.lock", ".Rprofile", "setup.R", "renv/activate.R", "renv/settings.json"):
        assert (R_DIR / rel).is_file(), rel


def test_mslipidmapper_is_pinned_from_github():
    lock = json.loads((R_DIR / "renv.lock").read_text(encoding="utf-8"))
    pkg = lock["Packages"]["MSLipidMapper"]
    assert pkg["Source"] == "GitHub"
    assert pkg["RemoteUsername"] == "systemsomicslab"
    assert pkg["RemoteRepo"] == "MSLipidMapper"
    assert pkg["RemoteSha"].startswith("097a88b")


def test_measurement_dependencies_are_locked():
    lock = json.loads((R_DIR / "renv.lock").read_text(encoding="utf-8"))
    for name in ("SummarizedExperiment", "S4Vectors", "jsonlite", "yaml", "renv"):
        assert name in lock["Packages"], name


def test_bridge_points_at_this_environment():
    from metabolomix.pathway import bridge
    assert bridge.RENV_PROJECT == R_DIR
