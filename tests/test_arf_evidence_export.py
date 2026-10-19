"""ARFの元候補と手動訂正を消さずに出力する。"""
import csv


def test_arf_keeps_unknown_and_rejected(tmp_path):
    from metabolomix.arf.evidence_export import export_arf_evidence
    result = dict(kind='two_group', a='A', b='B', n_a=3, n_b=3, q_threshold=.05,
                  log2fc_threshold=1, contract_version=1, log2fc_sign='positive means group_b is higher',
                  results=[dict(feature=f'Spot_{i}_height', log2fc=1, p=.01, q=.02,
                                mean_a=10, mean_b=20) for i in (0,1)])
    catalog = {0:{'Name':'a', 'InChIKey':'AAAAAAAAAAAAAA-BBBBBBBBBB-C'}, 1:{'Name':'Unknown'}}
    flags = dict(wrong={0}, suspect=set(), assign={}, redundant=set(), digest='fixture')
    out = export_arf_evidence(result, catalog, tmp_path/'out', source_analysis_id='fixture',
                              preprocess={}, flags=flags)
    assert out['n_features'] == 2
    with (tmp_path/'out'/'annotations.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert len(rows) == 1
    assert rows[0]['status'] == 'rejected'


def test_alignment_links_use_index_and_validate_mz_rt(tmp_path, monkeypatch):
    from metabolomix.arf.evidence_export import alignment_links
    from metabolomix.dcl import reader
    path = tmp_path/'AlignmentResult_test.dcl'
    path.write_bytes(b'fixture')
    monkeypatch.setattr(reader,'deserialize_dcl',lambda *a,**k:[
        dict(dcl_index=0,precursor_mz=200,rt=1,n_msms_peaks=1,msms_spectrum=[[100,20]]),
        dict(dcl_index=1,precursor_mz=600,rt=9,n_msms_peaks=1,msms_spectrum=[[100,20]])])
    rows = [dict(feature=f'Spot_{i}_height') for i in range(2)]
    catalog = {i:dict(MassCenter=200,RT=1) for i in range(2)}
    links, statuses = alignment_links(rows,catalog,path,'source')
    assert len(links) == 1
    assert links[0]['spectrum_id'] == 'AlignmentResult_test:0'
    assert statuses['Spot_1_height']['status'] == 'precursor_mismatch'
