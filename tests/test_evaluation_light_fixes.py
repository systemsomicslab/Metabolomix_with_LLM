"""評価で遭遇したURI・兄弟ファイル・保存先の回帰テスト。"""
from pathlib import Path

import pytest

from metabolomix.console.validation import map_assays
from metabolomix.core import session_state
from metabolomix.tools.reports import write_report, read_report
from tests.config_fixtures import use_config


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


def _open_analysis_dir(monkeypatch, directory):
    """解析フォルダが `directory` のデータセットを直近にする（記録の置き場所の基準）。"""
    from metabolomix.mztab.dataset_state import DatasetState

    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    ds = DatasetState()
    ds.analysis_dir = str(directory)
    session_state.session.dataset = ds


def test_report_override_wins_over_analysis_directory(tmp_path, monkeypatch):
    input_dir = tmp_path / 'input' / 'neg'
    output_dir = tmp_path / 'trial-reports'
    input_dir.mkdir(parents=True)
    _open_analysis_dir(monkeypatch, input_dir)
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": output_dir})
    write_report('trial-1', 'synthetic', '保存内容')
    assert (output_dir / 'trial-1.md').is_file()
    assert not (input_dir / 'reports').exists()
    assert '保存内容' in read_report('trial-1')


def test_report_without_override_uses_analysis_directory(tmp_path, monkeypatch):
    _open_analysis_dir(monkeypatch, tmp_path)
    write_report('local-1', 'synthetic', '保存内容')
    assert (tmp_path / 'reports/local-1.md').is_file()
