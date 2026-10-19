"""ARF比較結果に兄弟ARF2の候補を付け、未同定行も保持する。"""
from metabolomix.arf.identity_join import spot_id_of
from metabolomix.analysis.evidence_export import (
    comparison_of,feature_row,annotation_row,curate_annotations,write_evidence_bundle,exclusion_steps)
from metabolomix.analysis.export_contract import is_significant


def alignment_links(rows,catalog,dcl_path,source_analysis_id):
    """同一語幹のalignment DCLをIDで接続し、m/zとRTは検証だけに使う。"""
    from metabolomix.dcl.reader import deserialize_dcl
    from metabolomix.msdial.alignment_spot import precursor_agrees
    entries = deserialize_dcl(str(dcl_path))
    links, statuses = [], {}
    for row in rows:
        fid = str(row['feature'])
        spot = spot_id_of(fid)
        raw = catalog.get(spot,{})
        if spot is None or not 0 <= spot < len(entries):
            statuses[fid] = dict(status='out_of_range')
            continue
        entry = entries[spot]
        if raw.get('MassCenter') is None or raw.get('RT') is None:
            status = 'metadata_missing'
        elif not precursor_agrees(entry['precursor_mz'],raw['MassCenter'],entry['rt'],raw['RT']):
            status = 'precursor_mismatch'
        elif not entry['msms_spectrum']:
            status = 'no_msms'
        else:
            status = 'ok'
            links.append(dict(feature_id=fid,spectrum_file=dcl_path.name,
                spectrum_id=f'{dcl_path.stem}:{spot}',source_analysis_id=source_analysis_id,
                mapping_evidence='same-alignment-stem; MasterAlignmentID == DCL index; precursor/RT verified'))
        statuses[fid] = dict(status=status,dcl_index=spot)
    return links,statuses


def export_arf_evidence(result, catalog, output_dir, *, source_analysis_id,
                        preprocess, flags=None, provenance=None, dcl_path=None, preprocessing_report=None):
    comparison = comparison_of(result,preprocess)
    features, annotations, mapping = [], [], {}
    for row in result['results']:
        fid = str(row['feature'])
        spot = spot_id_of(fid)
        if spot is None:
            raise ValueError('ARF特徴量のMasterAlignmentIDを確定できません')
        mapping[fid] = spot
        raw = catalog.get(spot,{})
        meta = dict(mz=raw.get('MassCenter'),rt=raw.get('RT'))
        feature = feature_row(row,meta,spot)
        feature['significant'] = is_significant(q=row.get('q'),log2fc=row.get('log2fc'),
            q_threshold=result['q_threshold'],log2fc_threshold=result['log2fc_threshold'])
        features.append(feature)
        name = raw.get('Name') or ''
        if name.lower() != 'unknown' and (name or raw.get('InChIKey')) or raw.get('InChIKey'):
            annotations.append(annotation_row(fid,'msdial_arf2',dict(name=name,
                inchikey=raw.get('InChIKey'),smiles=raw.get('SMILES'),ontology=raw.get('Ontology'))))
    annotations = curate_annotations(annotations,flags,mapping)
    links,statuses = alignment_links(result['results'],catalog,dcl_path,source_analysis_id) if dcl_path else ([],{})
    return write_evidence_bundle(comparison=comparison,features=features,annotations=annotations,
                                  spectrum_links=links,spectra_files=[dcl_path] if links else [],output_dir=output_dir,
                                  source_analysis_id=source_analysis_id,
                                  provenance=dict(provenance or {},spectrum_links_status='verified_alignment' if dcl_path else 'not_provided',
                                                  msms_link_statuses=statuses),preprocessing_report=preprocessing_report,
                                  exclusions=exclusion_steps(preprocessing_report))
