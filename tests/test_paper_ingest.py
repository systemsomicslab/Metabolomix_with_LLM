"""paper_ingest のロジック検証。ネットワークは _http_get_json をモックする。"""

import tempfile
import unittest
from pathlib import Path

from metabolomix.corpus import knowledge_store as ks
from metabolomix.corpus import paper_ingest


def _epmc_result(**over):
    base = {
        "source": "MED",
        "title": "Plasmalogen decline under oxidative stress",
        "abstractText": "We show ether phospholipids decline in macrophages.",
        "doi": "10.1000/Abc",
        "pmid": "111",
        "pubYear": "2024",
        "journalInfo": {"journal": {"title": "J Lipid Res"}, "yearOfPublication": "2024"},
        "pubTypeList": {"pubType": ["research-article"]},
    }
    base.update(over)
    return base


class SanitizeTests(unittest.TestCase):
    def test_removes_control_chars_and_collapses(self):
        self.assertEqual(paper_ingest.sanitize_text("a\x00b\t c\n\nd"), "a b c d")

    def test_truncates(self):
        out = paper_ingest.sanitize_text("x" * 5000, maxlen=100)
        self.assertTrue(out.endswith("[…truncated]"))
        self.assertLessEqual(len(out), 100 + len(" […truncated]"))


class SearchEuropePmcTests(unittest.TestCase):
    def setUp(self):
        self._orig = paper_ingest._http_get_json

    def tearDown(self):
        paper_ingest._http_get_json = self._orig

    def _patch(self, results):
        paper_ingest._http_get_json = lambda url, params=None: {"resultList": {"result": results}}

    def test_filters_preprint_missing_abstract_missing_journal(self):
        self._patch([
            _epmc_result(),  # ok
            _epmc_result(source="PPR", doi="10.1/ppr"),  # preprint -> excluded
            _epmc_result(abstractText="", doi="10.1/noabs"),  # no abstract -> excluded
            _epmc_result(journalInfo={"journal": {}}, doi="10.1/nojourn"),  # no journal -> excluded
        ])
        out = paper_ingest.search_europepmc("plasmalogen", max_results=10)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["doi"], "10.1000/abc")  # lowercased
        self.assertFalse(out[0]["europepmc_retracted"])

    def test_flags_retracted_by_pubtype(self):
        self._patch([_epmc_result(pubTypeList={"pubType": ["Retracted Publication"]})])
        out = paper_ingest.search_europepmc("x")
        self.assertTrue(out[0]["europepmc_retracted"])

    def test_network_error_returns_error_dict(self):
        def boom(url, params=None):
            raise RuntimeError("down")
        paper_ingest._http_get_json = boom
        out = paper_ingest.search_europepmc("x")
        self.assertIn("error", out[0])


class RetractionTests(unittest.TestCase):
    def setUp(self):
        self._orig = paper_ingest._http_get_json

    def tearDown(self):
        paper_ingest._http_get_json = self._orig

    def test_excludes_epmc_retracted(self):
        cands = [{"title": "t", "doi": None, "europepmc_retracted": True}]
        self.assertEqual(paper_ingest.check_retraction(cands), [])

    def test_excludes_crossref_retracted(self):
        paper_ingest._http_get_json = lambda url, params=None: {
            "message": {"updated-by": [{"type": "retraction"}]}
        }
        cands = [{"title": "t", "doi": "10.1/x", "europepmc_retracted": False}]
        self.assertEqual(paper_ingest.check_retraction(cands), [])

    def test_keeps_non_retracted(self):
        paper_ingest._http_get_json = lambda url, params=None: {"message": {}}
        cands = [{"title": "t", "doi": "10.1/x", "europepmc_retracted": False}]
        self.assertEqual(len(paper_ingest.check_retraction(cands)), 1)


class DedupTests(unittest.TestCase):
    def test_drops_existing_doi_and_title(self):
        existing = {"doi:10.1/x", "title:" + paper_ingest._normalize_title("Already Here")}
        cands = [
            {"title": "new one", "doi": "10.1/x"},          # doi dup
            {"title": "Already Here", "doi": "10.1/y"},     # title dup
            {"title": "fresh", "doi": "10.1/z"},            # kept
        ]
        out = paper_ingest.deduplicate(cands, existing)
        self.assertEqual([c["doi"] for c in out], ["10.1/z"])

    def test_drops_intra_batch_duplicate(self):
        cands = [{"title": "same", "doi": "10.1/a"}, {"title": "same2", "doi": "10.1/a"}]
        out = paper_ingest.deduplicate(cands, set())
        self.assertEqual(len(out), 1)


class StageNoteTests(unittest.TestCase):
    def test_requires_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                paper_ingest.stage_note(
                    tmp, title="t", abstract="a", source="", found_for="A/Q1", query="q"
                )

    def test_writes_quarantined_note(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = paper_ingest.stage_note(
                tmp,
                title="Plasmalogen\x00 study",
                abstract="ether lipids decline",
                source="Doe 2024, DOI:10.1/x",
                found_for="2026-06-13-x/Q2",
                query="plasmalogen oxidation",
                relevance_score=0.8,
                doi="10.1/X",
            )
            self.assertTrue(path.is_file())
            self.assertEqual(path.parent.name, "_inbox")
            meta, body = ks.parse_frontmatter(path.read_text(encoding="utf-8"))
            self.assertEqual(meta["status"], "pending")
            self.assertEqual(meta["claim_strength"], "speculative")
            self.assertEqual(meta["type"], "knowledge")
            self.assertEqual(meta["doi"], "10.1/x")  # lowercased
            self.assertEqual(meta["found_for"], "2026-06-13-x/Q2")
            self.assertIn("Plasmalogen study", meta["title"])  # control char removed
            self.assertIn("ether lipids decline", body)

    def test_existing_identifiers_collects_doi_and_title(self):
        with tempfile.TemporaryDirectory() as tmp:
            ks.write_note(
                tmp, "n1",
                {"type": "knowledge", "title": "Known Paper", "doi": "10.1/Known", "source": "S"},
                "body",
            )
            ids = paper_ingest.existing_identifiers(tmp)
            self.assertIn("doi:10.1/known", ids)
            self.assertIn("title:" + paper_ingest._normalize_title("Known Paper"), ids)


if __name__ == "__main__":
    unittest.main()
