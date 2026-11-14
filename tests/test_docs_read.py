"""`docs_read`: output-format の文書をツールで読む（MCP のリソースを読めないクライアント向け）。"""
import json

import pytest

import server  # noqa: F401  ツールの登録
from metabolomix.core import mcp_core, session_state
from metabolomix.tools import resources
from metabolomix.tools.resources import docs_read


@pytest.fixture(autouse=True)
def _fresh_session():
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


def test_no_topic_returns_core_same_as_resource():
    expected = resources.output_format_reference()
    session_state.session.output_format_seen = False
    assert docs_read() == expected
    assert docs_read(topic=None) == expected


def test_core_topic_same_as_resource():
    assert docs_read("core") == resources.output_format_section("core")


@pytest.mark.parametrize("topic", list(mcp_core.OUTPUT_FORMAT_SECTIONS))
def test_every_topic_same_as_resource(topic):
    assert docs_read(topic) == resources.output_format_section(topic)


def test_flag_set_for_none_and_core_only():
    session_state.session.output_format_seen = False
    docs_read("dataset")
    assert session_state.session.output_format_seen is False
    docs_read()
    assert session_state.session.output_format_seen is True
    session_state.session.output_format_seen = False
    docs_read("core")
    assert session_state.session.output_format_seen is True


def test_unknown_topic_is_error_listing_valid_topics():
    payload = json.loads(docs_read("nope"))
    err = payload["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    assert err["details"]["field"] == "topic"
    for topic in mcp_core.OUTPUT_FORMAT_SECTIONS:
        assert topic in err["message"]
    assert err["details"]["valid"] == list(mcp_core.OUTPUT_FORMAT_SECTIONS)


def test_unknown_topic_does_not_set_flag():
    session_state.session.output_format_seen = False
    docs_read("../etc")
    assert session_state.session.output_format_seen is False


def test_outline_lists_level2_headings_with_sizes():
    text = resources.output_format_section("reanalysis")
    out = docs_read("reanalysis", outline=True)
    rows = out.split("heading\tchars\n", 1)[1].strip().splitlines()
    level2 = [line[3:].strip() for line in text.splitlines() if line.startswith("## ")]
    assert [row.split("\t")[0] for row in rows] == level2
    assert all(int(row.split("\t")[1]) > 0 for row in rows)


def test_section_returns_one_section_including_subsections():
    text = resources.output_format_section("reanalysis")
    first = next(line[3:].strip() for line in text.splitlines() if line.startswith("## "))
    body = docs_read("reanalysis", section=first)
    assert body.startswith(f"## {first}")
    assert body in text
    assert len(body) < len(text)
    assert "\n## " not in body


def test_section_matches_case_insensitive_substring_and_reports_misses():
    text = resources.output_format_section("dataset")
    first = next(line[3:].strip() for line in text.splitlines() if line.startswith("## "))
    titles = [line[3:].strip() for line in text.splitlines() if line.startswith("## ")]
    needle = first.upper()
    hits = [t for t in titles if needle.lower() in t.lower()]
    picked = docs_read("dataset", section=needle)
    if len(hits) == 1:
        assert picked.startswith(f"## {first}")
    else:
        assert json.loads(picked)["error"]["details"]["headings"] == hits
    miss = json.loads(docs_read("dataset", section="no such heading xyz"))
    assert miss["error"]["code"] == "DATASET_BAD_REQUEST"
    assert miss["error"]["details"]["field"] == "section"
    assert first in miss["error"]["details"]["headings"]


def test_section_without_topic_reads_the_core():
    core = resources.output_format_reference()
    first = next(line[3:].strip() for line in core.splitlines() if line.startswith("## "))
    assert docs_read(section=first).startswith(f"## {first}")


def test_section_does_not_set_the_seen_flag():
    session_state.session.output_format_seen = False
    docs_read(outline=True)
    assert session_state.session.output_format_seen is False
