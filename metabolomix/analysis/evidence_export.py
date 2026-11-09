"""全特徴量と出典別候補を別表で保存する根拠バンドル。

dataset 経路の v3（v2 に試料表・試料ごとの値と検出・付加体・極性を足したもの）。
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

from metabolomix.analysis import alignment_qc, export_contract
from metabolomix.analysis.dataset_analysis import PreconditionError, build_dataset_pp_inputs
from metabolomix.analysis.result_state import assert_current
from metabolomix.analysis.stats_results import as_two_group, get_matrix, sample_ids_by_name
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
    """manifest の `comparison`。`result` は `as_two_group` の出力で、実際の検定名と `statistics`
    （検定・CI 水準・判定規則・BH の母集団・単位）を持つ。`allow_confounded=true` で完全交絡を越えた比較なら
    `unadjusted_confounded: true` も足す。`n_a` / `n_b` は注入の数で、`samples.tsv` の群 a / b の行数と一致する。"""
    statistics = result['statistics']
    comparison = dict(group_a=result['a'], group_b=result['b'], n_a=result['n_a'], n_b=result['n_b'],
                      log2fc_sign=export_contract.LOG2FC_SIGN,
                      method=statistics['test'],correction='BH',
                      q_threshold=result['q_threshold'],log2fc_threshold=result['log2fc_threshold'],
                      preprocess=clean(preprocess),log_transform=result.get('log_transform'))
    # spec §13: 書き出しのメタ行と同じ定義を manifest にも置く（受け手は未知の欄を許す）
    comparison.update(contrast_id=result.get('contrast_id'),statistics=clean(statistics),
                      effect_size_definition=clean(result.get('effect_size_definition')),
                      mean_definition=export_contract.MEAN_DEFINITION,
                      matrix_id=result.get('matrix_id'),rule_hash=result.get('rule_hash'))
    # allow_confounded=true で完全交絡を越えた比較（as_two_group が来歴に足した印）。
    # 来歴そのもの（manifest の provenance.result.comparison）にも同じ形で載る。
    if ((result.get('provenance') or {}).get('comparison') or {}).get('unadjusted_confounded'):
        comparison['unadjusted_confounded'] = True
    return comparison


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


DATASET_LINK_EVIDENCE = 'same-alignment-stem; SMF_ID == AlignmentID == DCL index; precursor/RT verified'


def link_alignment_spectra(items, catalog, entries, dcl_path, source_analysis_id, mapping_evidence):
    """(特徴 ID, スポット番号) を同じアラインメントの `.dcl` の同じ番号に結ぶ（両経路で共有）。

    結ぶのは番号だけで、`.arf2` のスポットの m/z・RT と `.dcl` の precursor・RT が合うかは
    検証にだけ使う（m/z の近さでは結ばない）。番号が None・範囲外は `out_of_range`。
    戻り値は (links, 特徴ごとの状態)。
    """
    from metabolomix.msdial.alignment_spot import precursor_agrees
    links, statuses = [], {}
    for fid, spot in items:
        if spot is None or not 0 <= spot < len(entries):
            statuses[fid] = dict(status='out_of_range')
            continue
        raw = catalog.get(spot,{})
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
                              mapping_evidence=mapping_evidence))
        statuses[fid] = dict(status=status,dcl_index=spot)
    return links, statuses


def _spot_of(fid):
    """SMF_ID をスポット番号として読む。整数でなければ None。"""
    try:
        return int(fid)
    except (TypeError,ValueError):
        return None


def _id_space_check(ds, items, catalog):
    """SMF_ID と `.arf2` の AlignmentID が同じ番号の空間かを、mzTab と `.arf2` の m/z・RT で確かめる。"""
    from metabolomix.msdial.alignment_spot import precursor_agrees
    checked, mismatched = 0, []
    for fid, spot in items:
        raw = catalog.get(spot) if spot is not None else None
        meta = ds.feature_metadata.get(fid) or {}
        if not raw or raw.get('MassCenter') is None or meta.get('mz') is None:
            continue
        checked += 1
        if not precursor_agrees(meta['mz'],raw['MassCenter'],meta.get('rt'),raw.get('RT')):
            mismatched.append(dict(feature_id=fid,smf_mz=meta.get('mz'),smf_rt=meta.get('rt'),
                                   arf2_mz=raw.get('MassCenter'),arf2_rt=raw.get('RT')))
    return dict(n_checked=checked,n_mismatch=len(mismatched),examples=mismatched[:20])


def dataset_spectrum_links(ds, rows, arf2_path, dcl_path, source_analysis_id):
    """dataset 経路の spectrum_links。戻り値は (links, spectra_files, provenance に足す dict)。

    ARF 経路の `alignment_links` と同じ規則（SMF_ID == AlignmentID == `.dcl` の番号）。
    結べないときは空の links に `spectrum_links_status` と `spectrum_links_reason` を添える:
    - `not_provided`: 同じアラインメントの `.arf2` か `.dcl` が見つからない
    - `id_space_mismatch`: SMF_ID の m/z・RT が同じ番号の `.arf2` のスポットと合わない特徴がある
      （番号の空間が違う。1 件でもあれば全体を結ばない）
    - `id_space_unverified`: m/z を比べられる特徴が 1 件も無く、番号の空間を確かめられない
    結べたときは `verified_alignment`。特徴ごとの状態（`msms_link_statuses`）は ARF 経路と同じ形。
    """
    files = dict(arf2=Path(arf2_path).name if arf2_path else None,
                 dcl=Path(dcl_path).name if dcl_path else None)
    extra = dict(alignment_files=files,spectrum_links_reason=None,msms_link_statuses={},
                 msms_link_status_counts={},id_space_check=None)
    lacking = [name for name,path in (('arf2',arf2_path),('dcl',dcl_path))
               if path is None or not Path(path).is_file()]
    if lacking:
        return [], [], dict(extra,spectrum_links_status='not_provided',
                            spectrum_links_reason='同じアラインメントの '+' と '.join('.'+n for n in lacking)
                            +' が見つからないので、MS/MS を結んでいない（links は未提供で、MS/MS が無いという意味ではない）')
    from metabolomix.arf2.reader import load_catalog
    from metabolomix.dcl.reader import deserialize_dcl
    dcl_path = Path(dcl_path)
    catalog = {spot.get('AlignmentID'): spot for spot in load_catalog(Path(arf2_path))}
    items = [(str(row['feature']),_spot_of(row['feature'])) for row in rows]
    check = _id_space_check(ds,items,catalog)
    extra['id_space_check'] = check
    if check['n_mismatch']:
        return [], [], dict(extra,spectrum_links_status='id_space_mismatch',
                            spectrum_links_reason=f"{check['n_mismatch']} / {check['n_checked']} 特徴で SMF_ID の m/z・RT が"
                            '同じ番号の .arf2 のスポット（AlignmentID）と合わない。番号の空間が違うので MS/MS を結んでいない')
    if not check['n_checked']:
        return [], [], dict(extra,spectrum_links_status='id_space_unverified',
                            spectrum_links_reason='m/z を .arf2 と比べられる特徴が無く、SMF_ID と AlignmentID が'
                            '同じ番号の空間かを確かめられないので、MS/MS を結んでいない')
    links, statuses = link_alignment_spectra(
        items,catalog,deserialize_dcl(str(dcl_path)),dcl_path,source_analysis_id,DATASET_LINK_EVIDENCE)
    for fid,spot in items:
        if spot is None:
            statuses[fid] = dict(status='spot_id_not_integer')
    counts = {}
    for value in statuses.values():
        counts[value['status']] = counts.get(value['status'],0) + 1
    return links, [dcl_path] if links else [], dict(extra,spectrum_links_status='verified_alignment',
                                                    msms_link_statuses=statuses,msms_link_status_counts=counts)


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


def _engine_intensity_rows(ds, record, result):
    """統計エンジンの結果の `intensities.tsv` の行と段・検出が無い理由。

    値は結果が読んだ行列（`ds.stat_matrices[matrix_id]`）の `sample_values`——正規化の後・log2 の前の
    線形値で、欠測は NaN（エンジンは検定に補完しないので段は常に `pre_imputation`）。試料は `samples_a` / `samples_b` の表示名を
    エンジンの試料 ID へ結んで引く。行列が登録に無い（除外・試料の表の変更で古くなった）ときは
    `STATS_MATRIX_MISSING`（ツールが `missing_state` にする）。
    """
    entry = get_matrix(ds,record['matrix_id'])
    if entry is None or entry['prepared'].sample_values is None:
        raise DomainError('STATS_MATRIX_MISSING',
                          f"結果が読んだ行列 {record['matrix_id']} が登録に無いので、試料ごとの値を書けません。",
                          {'matrix_id':record['matrix_id']})
    prepared = entry['prepared']
    column = {str(fid):j for j,fid in enumerate(prepared.feature_ids)}
    wanted = [str(row['feature']) for row in result['results']]
    if set(wanted) != set(column):
        raise DomainError('EVIDENCE_INTENSITIES_MISMATCH',
                          '結果が読んだ行列の特徴の集合が、結果の行と一致しません。',
                          {'only_in_result':sorted(set(wanted)-set(column))[:20],
                           'only_in_matrix':sorted(set(column)-set(wanted))[:20]})
    samples = [*result['samples_a'],*result['samples_b']]
    to_id = sample_ids_by_name(ds)
    position = {sample_id:i for i,sample_id in enumerate(prepared.sample_ids)}
    mask, reason = _detection_source(ds,samples)
    mask_row = {str(fid):i for i,fid in enumerate(ds.feature_ids)}
    mask_column = {name:i for i,name in enumerate(ds.sample_names)}
    values = np.asarray(prepared.sample_values,dtype=float)
    rows = []
    for fid in wanted:
        for name in samples:
            row = position.get(to_id.get(name,name))
            rows.append(dict(feature_id=fid,sample_id=name,
                             value=None if row is None else values[row,column[fid]],
                             detected=None if mask is None
                             else bool(mask[mask_row[fid],mask_column[name]])))
    return rows, 'pre_imputation', reason


def _engine_exclusions(entry):
    """エンジンの行列の登録にある除外（`dataset_exclude`）。名前は 20 件で切り `n_more` を添える。"""
    if entry is None:
        return []
    out = []
    for stage,key in (('excluded_samples','samples'),('excluded_features','features')):
        names = list(entry['exclusions'].get(key) or [])
        out.append(dict(stage=stage,removed=len(names),names=names[:20],n_more=max(0,len(names)-20)))
    return out


def export_dataset_evidence(ds, result, output_dir, *, curation=None, polarity=None,
                            arf2_path=None, dcl_path=None, contrast_id=None):
    """dataset 経路の根拠バンドル v3。`arf2_path` / `dcl_path` は同じアラインメントの `.arf2` / `.dcl`
    （spectrum_links を結ぶのに使う。`dataset_spectrum_links`）。

    `result` は統計エンジンの `compare` の記録（`as_two_group` で 1 つの対比を読む。対比が複数なら
    `contrast_id` が要る）。`comparison.method` は実際の検定名で `comparison.statistics` を足し、
    `intensities.tsv` は行列の `sample_values`（`pre_imputation`）。"""
    if getattr(ds,'exploratory_only',False):
        raise DomainError('EXPLORATORY_ONLY_DATASET','探索専用データセットは比較の根拠として出力できません',{})
    record = result
    assert_current(ds,record)
    result = as_two_group(record,contrast_id)
    provenance = result.get('provenance') or {}
    entry = get_matrix(ds,record['matrix_id'])
    comparison = comparison_of(result,dict(entry['config']) if entry else {})
    comparison['polarity'] = polarity
    try:
        samples = _samples_table(ds,result)
        intensities, stage, detection_reason = _engine_intensity_rows(ds,record,result)
    except PreconditionError as exc:
        raise DomainError('EVIDENCE_INTENSITIES_UNAVAILABLE',
                          f'試料ごとの値を取り出せません: {exc.message}',exc.details) from exc
    adducts = alignment_qc._feature_adducts(ds)
    features, annotations, mapping = [], [], {}
    for row in result['results']:
        fid = str(row['feature'])
        meta = ds.feature_metadata.get(row['feature'], {})
        feature = feature_row(row,meta)
        feature['significant'] = export_contract.row_significant(row,   # 書き出しの significant 列と同じ（R16）
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
    links, spectra, link_provenance = dataset_spectrum_links(ds,result['results'],arf2_path,dcl_path,source_id)
    return write_evidence_bundle(comparison=comparison,features=features,annotations=annotations,
                                  spectrum_links=links,spectra_files=spectra,output_dir=output_dir,
                                  source_analysis_id=source_id,provenance=dict(result=provenance,
                                      source_files=getattr(ds,'source_files',[]),curation_state=(curation or {}).get('state'),
                                      **link_provenance),
                                  exclusions=_engine_exclusions(entry),samples=samples,
                                  intensities=intensities,intensities_stage=stage,
                                  detection_unavailable_reason=detection_reason)
