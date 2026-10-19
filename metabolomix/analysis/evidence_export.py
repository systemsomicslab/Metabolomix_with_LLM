"""全特徴量と出典別候補を別表で保存するv2契約。旧v1は変更しない。"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from metabolomix.analysis import export_contract
from metabolomix.analysis.result_state import assert_current
from metabolomix.core.atomic_io import DomainError

FEATURE_COLUMNS = ['feature_id','spot_id','mz','rt','log2fc','p_value','q_value','mean_a','mean_b','significant']
ANNOTATION_COLUMNS = ['annotation_id','feature_id','source','name','inchikey','rank',
                      'provenance_json','curation_json','status','structure_ids_json','scores_json']
LINK_COLUMNS = ['feature_id','spectrum_file','spectrum_id','source_analysis_id','mapping_evidence']
KEY_RE = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')


def clean(value):
    """欠測とnumpyのscalarを標準JSONへ写す。"""
    if isinstance(value, dict):
        return {str(k): clean(v) for k,v in value.items()}
    if isinstance(value, (list,tuple,set)):
        return [clean(v) for v in value]
    if hasattr(value, 'item'):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def packed(value):
    return json.dumps(clean(value), ensure_ascii=False, sort_keys=True, separators=(',',':'), allow_nan=False)


def write_evidence_bundle(*, comparison, features, annotations, spectrum_links,
                           spectra_files, output_dir, source_analysis_id=None,
                           exclusions=None, provenance=None, preprocessing_report=None):
    """同じ親の一時ディレクトリを完成させてから公開する。既存出力は置換しない。"""
    out = Path(output_dir).resolve()
    if out.exists():
        raise ValueError('出力先が既に存在します。新しいディレクトリを指定してください')
    out.parent.mkdir(parents=True, exist_ok=True)
    if not features or len({str(f['feature_id']) for f in features}) != len(features):
        raise ValueError('特徴量は空でなく、一意IDが必要です')
    source_id = source_analysis_id or hashlib.sha256(packed([comparison,features]).encode()).hexdigest()
    with tempfile.TemporaryDirectory(prefix='.'+out.name+'-', dir=out.parent) as temp:
        staging = Path(temp) / 'bundle'
        staging.mkdir()
        files = []
        for role, rows, columns in [('features',features,FEATURE_COLUMNS),
                                     ('annotations',annotations,ANNOTATION_COLUMNS),
                                     ('spectrum_links',spectrum_links,LINK_COLUMNS)]:
            path = staging / (role+'.tsv')
            with path.open('w', encoding='utf-8', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                for raw in rows:
                    row = {k: clean(v) for k,v in raw.items()}
                    for k in columns:
                        if k.endswith('_json'):
                            row[k] = packed(row.get(k[:-5], {})) if k not in row else row[k]
                        elif isinstance(row.get(k), bool):
                            row[k] = str(row[k]).lower()
                    writer.writerow(row)
            files.append(dict(role=role,path=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        for source in spectra_files:
            source = Path(source)
            if (staging/source.name).exists():
                raise ValueError('スペクトルファイル名が衝突しています')
            dest = staging/source.name
            shutil.copyfile(source,dest)
            files.append(dict(role='spectra',path=dest.name,sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
        manifest = dict(contract_name='massbank-context-evidence',contract_version=2,
                        source_analysis_id=source_id,exported_at=datetime.now(timezone.utc).isoformat(),
                        exporter=dict(name='ms-data-parser',version='evidence-v2'),
                        comparison=clean(comparison),population=dict(n_features_total=len(features),
                            population_complete=True,definition='前処理後の全2群比較対象',exclusions=exclusions or [],
                            n_features_before=(preprocessing_report or {}).get('features_before'),
                            preprocessing_report=clean(preprocessing_report),
                            derivation_status='recorded' if preprocessing_report else 'not_available'),
                        provenance=clean(provenance or {}),files=files)
        (staging/'manifest.json').write_text(packed(manifest)+'\n',encoding='utf-8')
        os.rename(staging, out)
    return dict(status='success', output_path=str(out/'manifest.json'),contract_version=2,
                n_features=len(features),n_candidates=len(annotations),
                msms_note='明示的なID対応が提供されたMS/MSだけを含む。空linksは未提供。')


def comparison_of(result, preprocess):
    if result.get('kind') != 'two_group' or result.get('log2fc_sign') != export_contract.LOG2FC_SIGN:
        raise ValueError('現行契約と向きが一致した2群比較が必要です')
    if result.get('contract_version') != export_contract.CONTRACT_VERSION:
        raise ValueError('差次的結果の契約版が不一致です')
    return dict(group_a=result['a'], group_b=result['b'], n_a=result['n_a'], n_b=result['n_b'],
                log2fc_sign=export_contract.LOG2FC_SIGN, method='Welch',correction='BH',
                q_threshold=result['q_threshold'],log2fc_threshold=result['log2fc_threshold'],
                preprocess=clean(preprocess),log_transform=result.get('log_transform'))


def feature_row(row, meta, spot_id=None):
    return dict(feature_id=str(row['feature']),spot_id=str(spot_id if spot_id is not None else row['feature']),
                mz=meta.get('mz'),rt=meta.get('rt'),log2fc=row.get('log2fc'),p_value=row.get('p'),
                q_value=row.get('q'),mean_a=row.get('mean_a'),mean_b=row.get('mean_b'))


def annotation_row(fid, source, raw, index=0):
    key = str(raw.get('inchikey') or '').strip()
    db = str(raw.get('database_identifier') or '')
    if not key and KEY_RE.fullmatch(db):
        key = db
    if key and not KEY_RE.fullmatch(key):
        raise ValueError(f'候補のInChIKeyが不正です: {fid}')
    name = str(raw.get('chemical_name') or raw.get('name') or '').strip()
    return dict(annotation_id=f'{fid}:{source}:{index}',feature_id=str(fid),source=source,
                name=name,inchikey=key,rank=raw.get('rank'),status='active',provenance=clean(raw),curation={},
                structure_ids={k:raw[k] for k in ('inchi','smiles','database_identifier') if raw.get(k)},
                scores={k:raw[k] for k in ('confidence_measure','confidence_value') if raw.get(k) is not None})


def curate_annotations(annotations, flags, fid_to_spot):
    if not flags:
        return annotations
    out = []
    for row in annotations:
        row = dict(row)
        spot = fid_to_spot.get(row['feature_id'])
        actions = [k for k in ('wrong','suspect','redundant','assign') if spot in flags.get(k, ())]
        row['curation'] = dict(actions=actions, flags_digest=flags.get('digest'))
        if 'wrong' in actions or 'redundant' in actions:
            row['status'] = 'rejected'
        elif 'assign' in actions:
            row['status'] = 'superseded'
        out.append(row)
    for fid,spot in fid_to_spot.items():
        assigned = flags.get('assign',{}).get(spot)
        if assigned and spot not in flags.get('wrong', ()) and spot not in flags.get('redundant', ()):
            row = annotation_row(fid,'curation',assigned)
            row['curation'] = dict(actions=['assign'],flags_digest=flags.get('digest'))
            out.append(row)
    return out


def export_dataset_evidence(ds, result, output_dir, *, curation=None):
    if getattr(ds,'exploratory_only',False):
        raise DomainError('EXPLORATORY_ONLY_DATASET','探索専用データセットは比較の根拠として出力できません',{})
    assert_current(ds,result)
    provenance = result.get('provenance') or {}
    preprocess, report = {}, None
    for parent_id in provenance.get('parent_ids', []):
        parent = getattr(ds, 'results', {}).get(parent_id, {})
        if parent.get('kind') == 'preprocess' or parent.get('provenance', {}).get('kind') == 'preprocess':
            preprocess = parent.get('provenance', {}).get('effective_parameters', {})
            report = parent.get('report') or {k:v for k,v in parent.items() if k != 'provenance'}
            break
    comparison = comparison_of(result,preprocess)
    features, annotations, mapping = [], [], {}
    for row in result['results']:
        fid = str(row['feature'])
        meta = ds.feature_metadata.get(row['feature'], {})
        feature = feature_row(row,meta)
        feature['significant'] = export_contract.is_significant(q=row.get('q'),log2fc=row.get('log2fc'),
            q_threshold=result['q_threshold'],log2fc_threshold=result['log2fc_threshold'])
        features.append(feature)
        candidates = getattr(ds,'feature_candidates',{}).get(row['feature'],[]) or []
        if candidates:
            annotations.extend(annotation_row(fid,'msdial_sme',candidate,i) for i,candidate in enumerate(candidates))
            if meta.get('inchikey') or meta.get('name'):
                annotations.append(annotation_row(fid,'msdial_resolved_annotation',meta))
            annotation = getattr(ds,'feature_annotations',{}).get(row['feature'],{})
            if annotation:
                annotations.append(annotation_row(fid,'msdial_sml',annotation))
        else:
            annotation = getattr(ds,'feature_annotations',{}).get(row['feature'],{}) or meta
            if annotation.get('name') or annotation.get('inchikey'):
                annotations.append(annotation_row(fid,'msdial_annotation',annotation))
        try:
            mapping[fid] = int(fid)
        except ValueError:
            mapping[fid] = None
    flags = curation.get('flag_set') if curation and curation.get('state') == 'applied' else None
    annotations = curate_annotations(annotations,flags,mapping)
    source_id = provenance.get('result_id') or hashlib.sha256(packed(result).encode()).hexdigest()
    return write_evidence_bundle(comparison=comparison,features=features,annotations=annotations,
                                  spectrum_links=[],spectra_files=[],output_dir=output_dir,
                                  source_analysis_id=source_id,provenance=dict(result=provenance,
                                      source_files=getattr(ds,'source_files',[]),curation_state=(curation or {}).get('state'),
                                      spectrum_links_status='not_provided'),preprocessing_report=report,
                                  exclusions=exclusion_steps(report))


def exclusion_steps(report):
    """元の段階別件数を写し、重複しうる件数を足し直さない。"""
    return [dict(stage=key,**clean(value)) for key,value in (report or {}).get('steps',{}).items()
            if isinstance(value,dict) and ('removed' in value or 'failed' in value)]
