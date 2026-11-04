"""massbank-context が写して読む見本のバンドル v3 が、今の書き出しと一致し続けること。"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from tests.evidence_bundle_fixture import build_fixture_bundle

FIXTURE = Path(__file__).parent / "fixtures" / "evidence_bundle_v3"
VOLATILE = ("exported_at", "source_analysis_id", "provenance")


def _tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_fixture_matches_current_exporter(tmp_path):
    build_fixture_bundle(tmp_path / "b")
    for name in ("features.tsv", "annotations.tsv", "spectrum_links.tsv",
                 "samples.tsv", "intensities.tsv"):
        assert (tmp_path / "b" / name).read_bytes() == (FIXTURE / name).read_bytes(), name
    fresh = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    kept = json.loads((FIXTURE / "manifest.json").read_text(encoding="utf-8"))
    assert {k: v for k, v in fresh.items() if k not in VOLATILE} == \
           {k: v for k, v in kept.items() if k not in VOLATILE}


def test_fixture_covers_the_cases_massbank_context_reads():
    rows = _tsv(FIXTURE / "intensities.tsv")
    assert {r["detected"] for r in rows} == {"true", "false"}
    assert any(r["value"] == "" for r in rows)
    assert "" in {r["group"] for r in _tsv(FIXTURE / "samples.tsv")}
    manifest = json.loads((FIXTURE / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["contract_version"] == 3 and manifest["comparison"]["polarity"] == "positive"


def test_fixture_has_no_local_paths():
    text = "".join(p.read_text(encoding="utf-8") for p in FIXTURE.iterdir())
    assert ":\\" not in text and ":/" not in text.replace("file:///", "")
