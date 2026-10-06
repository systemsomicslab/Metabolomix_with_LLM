"""評価で遭遇したURI・兄弟ファイル・保存先の回帰テスト。"""
from pathlib import Path

import pytest

from metabolomix.arf import tools as arf_tools
from metabolomix.console.validation import map_assays
from metabolomix.core import mcp_core, session_state
from metabolomix.tools.reports import write_report, read_report


@pytest.mark.parametrize('uri,source', [
    ('file://C:/data/a%20b.wiff', 'C:/data/a b.wiff'),
    ('file:///C:/data/a%20b.wiff', 'C:/data/a b.wiff'),
    ('file://nas/share/a%20b.wiff', '//nas/share/a b.wiff'),
    ('file://localhost/C:/data/a.wiff', 'C:/data/a.wiff'),
    ('file:///data/a%20b.wiff', '/data/a b.wiff'),
])
def test_file_uri_preserves_drive_and_unc_host(uri, source):
    parsed = {'metadata': {'ms_run[1]-location': uri, 'assay[1]-ms_run_ref': 'ms_run[1]'}}
    assert map_assays(parsed, [source]) == {'abundance_assay[1]': source}


def test_uri_does_not_match_another_directory_with_same_basename():
    parsed = {'metadata': {'ms_run[1]-location': 'file://C:/other/a.wiff',
                           'assay[1]-ms_run_ref': 'ms_run[1]'}}
    assert map_assays(parsed, ['C:/expected/a.wiff']) != {'abundance_assay[1]': 'C:/expected/a.wiff'}


@pytest.mark.parametrize('stem', ['AlignResult-2026981258', 'AlignmentResult_2026_09_09_18_03_05'])
def test_sibling_arf2_uses_exact_batch_stem(tmp_path, monkeypatch, stem):
    monkeypatch.setattr(session_state, 'session', session_state.AnalysisSession())
    session_state.session.arf.current_file_path = str(tmp_path / (stem + '_PeakProperties.arf'))
    expected = tmp_path / (stem + '.arf2')
    expected.touch()
    (tmp_path / 'AlignmentResult_2026_09_10_20_00_00.arf2').touch()
    assert arf_tools._sibling_arf2_path() == expected
    expected.unlink()
    assert arf_tools._sibling_arf2_path() is None


def test_report_override_wins_over_mutable_input_directory(tmp_path, monkeypatch):
    input_dir = tmp_path / 'input' / 'neg'
    output_dir = tmp_path / 'trial-reports'
    input_dir.mkdir(parents=True)
    monkeypatch.setattr(mcp_core, 'DATA_DIR', input_dir)
    monkeypatch.setenv('LIPIDMIX_REPORTS_DIR', str(output_dir))
    write_report('trial-1', 'synthetic', '保存内容')
    assert (output_dir / 'trial-1.md').is_file()
    assert not (input_dir / 'reports').exists()
    assert '保存内容' in read_report('trial-1')


def test_report_without_override_uses_input_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(mcp_core, 'DATA_DIR', tmp_path)
    monkeypatch.delenv('LIPIDMIX_REPORTS_DIR', raising=False)
    write_report('local-1', 'synthetic', '保存内容')
    assert (tmp_path / 'reports/local-1.md').is_file()
