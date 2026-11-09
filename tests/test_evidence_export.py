"""全特徴量を残す新しい根拠バンドル。"""
import csv
import json
import pytest


def test_writer_preserves_unknown_and_quoted_name(tmp_path):
    from metabolomix.analysis.evidence_export import write_evidence_bundle
    comparison = dict(group_a='A', group_b='B', n_a=3, n_b=3,
                      log2fc_sign='positive means group_b is higher', method='Welch',
                      correction='BH', q_threshold=.05, log2fc_threshold=1, preprocess={})
    features = [dict(feature_id='f0', spot_id='0', log2fc=1, p_value=.01, q_value=.02,
                     mean_a=10, mean_b=20, significant=True)]
    result = write_evidence_bundle(comparison=comparison, features=features, annotations=[],
                                   spectrum_links=[], spectra_files=[], output_dir=tmp_path/'out')
    manifest = json.loads((tmp_path/'out'/'manifest.json').read_text())
    assert manifest['population']['n_features_total'] == 1
    assert result['contract_version'] == 2
    assert len(list(csv.DictReader((tmp_path/'out'/'features.tsv').open(), delimiter='\t'))) == 1


def _engine_dataset(*, excluded_features=None):
    """8 検体 × 6 特徴（`make_dataset`）に、統計エンジンの 2 群比較（前半 A・後半 B）を登録する。"""
    from tests.engine_chain import compare_halves
    from tests.pipeline_fixtures import make_dataset
    ds = make_dataset()
    ds.excluded_features = dict(excluded_features or {})
    return ds, compare_halves(ds)


def test_dataset_export_keeps_all_candidates_and_unknown(tmp_path):
    from metabolomix.analysis.evidence_export import export_dataset_evidence
    ds, record = _engine_dataset()
    ds.feature_metadata['F0'] = {'inchikey': 'AAAAAAAAAAAAAA-BBBBBBBBBB-C', 'name': 'primary'}
    ds.feature_metadata['F1'] = {}
    ds.feature_candidates = {'F0': [
        {'rank': 1, 'chemical_name': 'a', 'inchikey': 'AAAAAAAAAAAAAA-BBBBBBBBBB-C'},
        {'rank': 2, 'chemical_name': 'b', 'inchikey': 'CCCCCCCCCCCCCC-DDDDDDDDDD-E'}]}
    exported = export_dataset_evidence(ds, record, tmp_path/'out')
    assert exported['n_features'] == 6        # 注釈の無い F1 も母集団に残る
    with (tmp_path/'out'/'annotations.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='	'))
    f0 = [r for r in rows if r['feature_id'] == 'F0']
    assert len(f0) == 3                       # SME の候補 2 件 ＋ 解決済みの注釈
    assert f0[0]['inchikey'] != f0[1]['inchikey']
    assert not [r for r in rows if r['feature_id'] == 'F1']


def test_writer_does_not_replace_existing_output(tmp_path):
    from metabolomix.analysis.evidence_export import write_evidence_bundle
    (tmp_path/'out').mkdir()
    with pytest.raises(ValueError):
        write_evidence_bundle(comparison={}, features=[], annotations=[], spectrum_links=[],
                               spectra_files=[], output_dir=tmp_path/'out')


def test_resolved_primary_key_survives_raw_sme_candidates(tmp_path):
    from metabolomix.analysis.evidence_export import export_dataset_evidence
    key = 'AAAAAAAAAAAAAA-BBBBBBBBBB-C'
    ds, record = _engine_dataset(excluded_features={'F5': 'bad peak'})
    ds.feature_metadata['F0'] = dict(name='resolved', inchikey=key, inchikey_source='msp',
                                     database_identifier='MSP:123')
    ds.feature_candidates = {'F0': [dict(rank=1, chemical_name='raw', database_identifier='MSP:123')]}
    ds.feature_annotations = {'F0': dict(name='SML-label')}
    export_dataset_evidence(ds, record, tmp_path/'out')
    with (tmp_path/'out'/'annotations.tsv').open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='	') if r['feature_id'] == 'F0']
    assert any(r['inchikey'] == key for r in rows)
    assert any(r['name'] == 'SML-label' for r in rows)
    manifest = json.loads((tmp_path/'out'/'manifest.json').read_text())
    assert manifest['population']['n_features_total'] == 5
    excluded = {e['stage']: e for e in manifest['population']['exclusions']}
    assert excluded['excluded_features']['removed'] == 1
    assert excluded['excluded_features']['names'] == ['F5']
