"""全特徴量と出典別候補を別表で保存する根拠バンドル。旧v1は変更しない。

dataset 経路は v3（v2 に試料表・試料ごとの値と検出・付加体・極性を足したもの）。
試料ごとの値の取り出しが無い ARF 経路は v2 のまま。
"""
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
import numpy as np

from metabolomix.analysis import alignment_qc, export_contract, feature_search
from metabolomix.analysis.dataset_analysis import (PreconditionError, build_dataset_pp_inputs,
                                                   run_dataset_preprocess)
from metabolomix.analysis.result_state import assert_current
from metabolomix.core.atomic_io import DomainError

FEATURE_COLUMNS = ['feature_id','spot_id','mz','rt','log2fc','p_value','q_value','mean_a','mean_b','significant']
ANNOTATION_COLUMNS = ['annotation_id','feature_id','source','name','inchikey','rank',
                      'provenance_json','curation_json','status','structure_ids_json','scores_json']
LINK_COLUMNS = ['feature_id','spectrum_file','spectrum_id','source_analysis_id','mapping_evidence']
BUNDLE_CONTRACT_VERSION = 3
FEATURE_COLUMNS_V3 = [*FEATURE_COLUMNS, 'adduct']
SAMPLE_COLUMNS = ['sample_id','group','role','injection_order']
INTENSITY_COLUMNS = ['feature_id','sample_id','value','detected']
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
                           exclusions=None, provenance=None, preprocessing_report=None,
                           samples=None, intensities=None, intensities_stage=None,
                           detection_unavailable_reason=None):
    """同じ親の一時ディレクトリを完成させてから公開する。既存出力は置換しない。

    `samples` と `intensities` を両方渡せば v3、どちらも無ければ v2（ARF 経路）。
    """
    if (samples is None) != (intensities is None):
        raise ValueError('samples と intensities は両方渡すか、両方省いてください')
    v3 = samples is not None
    version = BUNDLE_CONTRACT_VERSION if v3 else 2
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
        tables = [('features',features,FEATURE_COLUMNS_V3 if v3 else FEATURE_COLUMNS),
                  ('annotations',annotations,ANNOTATION_COLUMNS),
                  ('spectrum_links',spectrum_links,LINK_COLUMNS)]
        if v3:
            tables += [('samples',samples,SAMPLE_COLUMNS),
                       ('intensities',map(_intensity_cell,intensities),INTENSITY_COLUMNS)]
        for role, rows, columns in tables:
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
        manifest = dict(contract_name='massbank-context-evidence',contract_version=version,
                        source_analysis_id=source_id,exported_at=datetime.now(timezone.utc).isoformat(),
                        exporter=dict(name='ms-data-parser',version=f'evidence-v{version}'),
                        comparison=clean(comparison),population=dict(n_features_total=len(features),
                            population_complete=True,definition='前処理後の全2群比較対象',exclusions=exclusions or [],
                            n_features_before=(preprocessing_report or {}).get('features_before'),
                            preprocessing_report=clean(preprocessing_report),
                            derivation_status='recorded' if preprocessing_report else 'not_available'),
                        provenance=clean(provenance or {}),files=files)
        if v3:
            manifest['intensities_stage'] = intensities_stage
            manifest['population']['detection_unavailable_reason'] = detection_unavailable_reason
        (staging/'manifest.json').write_text(packed(manifest)+'\n',encoding='utf-8')
        os.rename(staging, out)
    summary = dict(status='success', output_path=str(out/'manifest.json'),contract_version=version,
                   n_features=len(features),n_candidates=len(annotations),
                   msms_note='明示的なID対応が提供されたMS/MSだけを含む。空linksは未提供。')
    if v3:
        summary.update(n_samples=len(samples),n_intensity_rows=len(intensities),
                       intensities_stage=intensities_stage)
    return summary


def _intensity_cell(row):
    """`intensities.tsv` の 1 行。NaN・None の値は空になる。"""
    return {**row,'value':export_contract.format_number(row['value'],'.10g')}


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


def _injection_order(value):
    """注入順。整数として読めなければ None（`samples.tsv` では空）。"""
    if isinstance(value,bool):
        return None
    try:
        number = float(value)
    except (TypeError,ValueError):
        return None
    return int(number) if number.is_integer() else None


def _samples_table(ds, result):
    """`samples.tsv` の行。全試料を書き、比較に入らない試料の `group` は空。"""
    _, _, _, roles, sample_meta = build_dataset_pp_inputs(ds)
    group_a, group_b = set(result.get('samples_a') or []), set(result.get('samples_b') or [])
    return [dict(sample_id=name,group='a' if name in group_a else 'b' if name in group_b else None,
                 role=roles.get(name),
                 injection_order=_injection_order((sample_meta.get(name) or {}).get('run_order')))
            for name in ds.sample_names]


def _intensity_recipe(ds, result):
    """値を当て直すレシピ。None は親の前処理の行列（`ds.pp_matrix`）を使えるという意味。

    どちらでもなければ止める——別の前処理の値を `post_imputation` と名乗って書かない。
    """
    recipe = feature_search._result_recipe(ds,result)
    parents = (result.get('provenance') or {}).get('parent_ids') or []
    if recipe is None and (getattr(ds,'preprocess_id',None) not in parents
                           or getattr(ds,'pp_matrix',None) is None):
        raise DomainError('EVIDENCE_INTENSITIES_UNAVAILABLE',
                          'この結果を計算した前処理のレシピが分からず、今の前処理の行列もこの結果の親では'
                          'ないので、試料ごとの値を書けません。dataset_differential で計算し直してください。',
                          {'result_id':(result.get('provenance') or {}).get('result_id'),
                           'parent_ids':parents,'preprocess_id':getattr(ds,'preprocess_id',None)})
    return recipe


def _detection_source(ds, samples):
    """検出の表と、使えないときの理由。使えなければ `detected` は全行空（false で埋めない）。"""
    mask = getattr(ds,'detected_mask',None)
    if mask is None:
        return None, (getattr(ds,'feature_qc',None) or {}).get('reason') or 'not_imported'
    mask = np.asarray(mask,dtype=bool)
    if mask.shape != (len(ds.feature_ids),len(ds.sample_names)):
        return None, 'mask_shape_mismatch'
    if not set(samples) <= set(ds.sample_names):
        return None, 'sample_axis_mismatch'
    return mask, None


def _intensity_rows(ds, result):
    """`intensities.tsv` の行（全特徴 × 比較の試料）と `intensities_stage`・検出が無い理由。

    レシピが分かれば補完だけを外して当て直す（補完は前処理の最後の段なので、ほかの段は
    比較に入った値と同じ）。`run_dataset_preprocess` は `ds` を書き換えない。
    値の行列に無い試料の `value` は空のまま。
    """
    recipe = _intensity_recipe(ds,result)
    if recipe is None:
        matrix, names, features = ds.pp_matrix, ds.pp_sample_names, ds.pp_feature_names
        stage = 'post_imputation'
    else:
        matrix, names, features, *_ = run_dataset_preprocess(ds,{**recipe,'impute':'none'})
        stage = 'pre_imputation'
    matrix = np.asarray(matrix,dtype=float)
    wanted = [str(row['feature']) for row in result['results']]
    column = {str(fid):j for j,fid in enumerate(features)}
    if set(wanted) != set(column):
        raise DomainError('EVIDENCE_INTENSITIES_MISMATCH',
                          '前処理を当て直した行列の特徴の集合が、結果の行と一致しません。',
                          {'only_in_result':sorted(set(wanted)-set(column))[:20],
                           'only_in_matrix':sorted(set(column)-set(wanted))[:20]})
    samples = [*(result.get('samples_a') or []),*(result.get('samples_b') or [])]
    if not samples:
        raise DomainError('EVIDENCE_INTENSITIES_UNAVAILABLE',
                          'この結果は群の試料（samples_a / samples_b）を持たないので、試料ごとの値を'
                          '書けません。dataset_differential で計算し直してください。',
                          {'result_id':(result.get('provenance') or {}).get('result_id')})
    position = {name:i for i,name in enumerate(names)}
    mask, reason = _detection_source(ds,samples)
    mask_row = {str(fid):i for i,fid in enumerate(ds.feature_ids)}
    mask_column = {name:i for i,name in enumerate(ds.sample_names)}
    rows = []
    for fid in wanted:
        for name in samples:
            row = position.get(name)
            rows.append(dict(feature_id=fid,sample_id=name,
                             value=None if row is None else matrix[row,column[fid]],
                             detected=None if mask is None
                             else bool(mask[mask_row[fid],mask_column[name]])))
    return rows, stage, reason


def export_dataset_evidence(ds, result, output_dir, *, curation=None, polarity=None):
    if getattr(ds,'exploratory_only',False):
        raise DomainError('EXPLORATORY_ONLY_DATASET','探索専用データセットは比較の根拠として出力できません',{})
    # 値の出所が無い結果は、古いかどうかより先に止める（出所の無さは計算し直す以外に直らない）。
    _intensity_recipe(ds,result)
    assert_current(ds,result)
    provenance = result.get('provenance') or {}
    preprocess, report = {}, None
    for parent_id in provenance.get('parent_ids', []):
        parent = getattr(ds, 'results', {}).get(parent_id, {})
        if parent.get('kind') == 'preprocess' or parent.get('provenance', {}).get('kind') == 'preprocess':
            preprocess = parent.get('provenance', {}).get('effective_parameters', {})
            report = parent.get('report') or {k:v for k,v in parent.items() if k != 'provenance'}
            break
    # 引き継いだ結果の親は ds.results に無いので、値の段はレシピから書く。
    comparison = comparison_of(result,preprocess or feature_search._result_recipe(ds,result) or {})
    comparison['polarity'] = polarity
    try:
        samples = _samples_table(ds,result)
        intensities, stage, detection_reason = _intensity_rows(ds,result)
    except PreconditionError as exc:
        raise DomainError('EVIDENCE_INTENSITIES_UNAVAILABLE',
                          f'試料ごとの値を取り出せません: {exc.message}',exc.details) from exc
    adducts = alignment_qc._feature_adducts(ds)
    features, annotations, mapping = [], [], {}
    for row in result['results']:
        fid = str(row['feature'])
        meta = ds.feature_metadata.get(row['feature'], {})
        feature = feature_row(row,meta)
        feature['significant'] = export_contract.is_significant(q=row.get('q'),log2fc=row.get('log2fc'),
            q_threshold=result['q_threshold'],log2fc_threshold=result['log2fc_threshold'])
        feature['adduct'] = adducts.get(fid)
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
                                  exclusions=exclusion_steps(report),samples=samples,
                                  intensities=intensities,intensities_stage=stage,
                                  detection_unavailable_reason=detection_reason)


def exclusion_steps(report):
    """元の段階別件数を写し、重複しうる件数を足し直さない。"""
    return [dict(stage=key,**clean(value)) for key,value in (report or {}).get('steps',{}).items()
            if isinstance(value,dict) and ('removed' in value or 'failed' in value)]
