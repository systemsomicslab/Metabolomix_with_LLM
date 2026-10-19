"""全特徴量を残す新しい根拠バンドル。"""
import csv
import json
from types import SimpleNamespace
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


def test_dataset_export_keeps_all_candidates_and_unknown(tmp_path):
    from metabolomix.analysis.evidence_export import export_dataset_evidence
    ds = SimpleNamespace(exploratory_only=False, feature_metadata={
        'f0': {'inchikey':'AAAAAAAAAAAAAA-BBBBBBBBBB-C','name':'primary'}, 'f1':{}},
        feature_candidates={'f0':[
            {'rank':1, 'chemical_name':'a', 'inchikey':'AAAAAAAAAAAAAA-BBBBBBBBBB-C'},
            {'rank':2, 'chemical_name':'b', 'inchikey':'CCCCCCCCCCCCCC-DDDDDDDDDD-E'}]},
        feature_annotations={}, source_files=[])
    result = dict(kind='two_group', a='A', b='B', n_a=3, n_b=3, q_threshold=.05,
                  log2fc_threshold=1, contract_version=1, log2fc_sign='positive means group_b is higher',
                  results=[dict(feature=f, log2fc=1, p=.01, q=.02, mean_a=10, mean_b=20) for f in ('f0','f1')])
    exported = export_dataset_evidence(ds, result, tmp_path/'out')
    assert exported['n_features'] == 2
    with (tmp_path/'out'/'annotations.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert len(rows) == 3
    assert rows[0]['inchikey'] != rows[1]['inchikey']


def test_writer_does_not_replace_existing_output(tmp_path):
    from metabolomix.analysis.evidence_export import write_evidence_bundle
    (tmp_path/'out').mkdir()
    with pytest.raises(ValueError):
        write_evidence_bundle(comparison={}, features=[], annotations=[], spectrum_links=[],
                               spectra_files=[], output_dir=tmp_path/'out')


def test_resolved_primary_key_survives_raw_sme_candidates(tmp_path):
    from metabolomix.analysis.evidence_export import export_dataset_evidence
    key='AAAAAAAAAAAAAA-BBBBBBBBBB-C'
    ds=SimpleNamespace(exploratory_only=False,feature_metadata={'f0':dict(name='resolved',inchikey=key,
        inchikey_source='msp',database_identifier='MSP:123')},feature_candidates={'f0':[
        dict(rank=1,chemical_name='raw',database_identifier='MSP:123')]},
        feature_annotations={'f0':dict(name='SML-label')},source_files=[],preprocess_id='prep',dataset_id=None,results={'prep':{
            'kind':'preprocess','provenance':{'effective_parameters':{'normalize':'none'}},
            'report':{'features_before':3,'features_after':1,'steps':{'blank_filter':{'removed':2}}}}})
    result=dict(kind='two_group',a='A',b='B',n_a=3,n_b=3,q_threshold=.05,log2fc_threshold=1,
        contract_version=1,log2fc_sign='positive means group_b is higher',
        provenance={'parent_ids':['prep'],'dataset_id':None,'result_id':'diff'},results=[dict(feature='f0',
        log2fc=1,p=.01,q=.02,mean_a=10,mean_b=20)])
    export_dataset_evidence(ds,result,tmp_path/'out')
    with (tmp_path/'out'/'annotations.tsv').open() as stream:
        rows=list(csv.DictReader(stream,delimiter='\t'))
    assert any(r['inchikey']==key for r in rows)
    assert any(r['name']=='SML-label' for r in rows)
    manifest=json.loads((tmp_path/'out'/'manifest.json').read_text())
    assert manifest['population']['n_features_before'] == 3
    assert manifest['population']['exclusions'][0]['removed'] == 2
