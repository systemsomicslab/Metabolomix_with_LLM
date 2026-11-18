"""DatasetState からの統計エンジンの入口と、検出・gap-fill・欠測の 3 状態（spec §5・§6）。

このモジュールは DatasetState を**読むだけ**で書き換えない。MCP・セッションには依存しない。
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, replace

import numpy as np

from metabolomix.core.atomic_io import DomainError

#: セルの状態。DETECTED = 実検出、GAP_FILLED = 値はあるが検出は偽（gap-fill）、
#: MISSING = 値が無い（0 以下・非有限）、PRESENT_UNKNOWN = 値はあるが検出の別が分からない。
DETECTED, GAP_FILLED, MISSING, PRESENT_UNKNOWN = 1, 2, 0, 3

_CORE_KEYS = frozenset({"sample_id", "role", "group", "batch", "injection_order",
                        "include", "subject_id", "biological_sample_id", "source_file"})
#: 試料の対応表（sample_manifest）の行が持つ、因子ではない帳簿の欄。来歴（provenance）・食い違い
#: （conflicts）・QC のプール名（qc_pool）は注入ごとに違いうるので、属性にすると同じ個体の反復注入が
#: 別の個体に割れる。
_BOOKKEEPING_KEYS = frozenset({"provenance", "conflicts", "qc_pool"})


@dataclass(frozen=True)
class SampleInfo:
    """試料（注入）1 本の情報。`attributes` は因子・共変量に使う追加の欄。"""
    sample_id: str
    role: str
    group: str | None
    batch: str | None
    injection_order: float | None
    subject_id: str | None
    include: bool
    attributes: dict[str, str] = field(default_factory=dict)
    #: 注入順の出所（`sheet` / `mzTab` / `mddata` / `raw_acquisition_time` / `filename`。値が無ければ None）。
    injection_order_source: str | None = None
    #: 注入順が MS-DIAL 由来で、データセットのその値がファイル名の順と一致する（MS-DIAL の既定かもしれない）。
    injection_order_may_be_file_order: bool = False


@dataclass(frozen=True)
class QuantState:
    """定量値が内部標準などで換算済みか（B-3 の `DatasetState.quant_state` から読む）。"""
    normalized: bool
    method: str | None
    unit: str | None


@dataclass
class StudyData:
    """n_samples × n_features の生の値と、試料の表・検出の別をまとめた入口。"""
    values: np.ndarray
    detected: np.ndarray | None
    sample_ids: list[str]
    feature_ids: list[str]
    samples: list[SampleInfo]
    quant: QuantState | None = None
    dataset_hash: str | None = None
    feature_classes: dict[str, str] = field(default_factory=dict)
    annotated_features: frozenset[str] = frozenset()
    #: 特徴 ID → 注釈の名前（`feature_metadata[fid]["name"]`。名前の無い特徴は持たない）。脂質アッセイの組成が
    #: 脂質と読めない注釈を脂質クラスから外し、MS-DIAL の脂質クラスへ移すのに読む（`composition`）。
    feature_names: dict[str, str] = field(default_factory=dict)
    #: 特徴 ID → 注釈の組成式（`feature_metadata[fid]["formula"]`。無い特徴は持たない）。脂質アッセイの組成が、名前で
    #: 読めない研究室 MSP の脂質を ClassyFire のクラスと組成式から脂質クラスへ救うのに読む（`composition`）。
    feature_formulas: dict[str, str] = field(default_factory=dict)
    #: 割れのまとめの記録（組ごとの特徴 ID）。除外との区別（matrix_id）のために持つ。
    feature_merges: tuple = ()
    #: バッチが不明な試料の、試料名の 8 桁からの候補（試料 ID → `{"value", "source": "filename_date"}`）。
    #: 統計には当てない（前処理が注意として返すだけ）。
    batch_candidates: dict[str, dict] = field(default_factory=dict)
    #: 陽陰を統合したデータセット（`dataset_merge_polarities`）だけ: 特徴 ID → 極性
    #: （`feature_metadata[fid]["origin"]["polarity"]`）。空なら統合ではなく、前処理は全特徴を 1 つの塊として扱う。
    feature_polarity: dict[str, str] = field(default_factory=dict)
    #: 統合したデータセットだけ: 極性 → 試料ごと（`samples` と同じ順）のその極性の測定
    #: （`run_of` の形。None はその試料の測定が分からない）。ここに無い極性は `samples` の値（統合した
    #: データセットの試料の表＝陽イオンの測定）を使う（`polarity_samples`）。
    polarity_runs: dict[str, list] = field(default_factory=dict)


#: 統合したデータセットの試料の表（`samples`）がどの極性の測定か（`quant.merge` は陽イオンの側を土台にする）。
SAMPLE_TABLE_POLARITY = "positive"
#: 試料の測定（注入）ごとに極性で違いうる欄（`polarity_runs` の形）。
RUN_FIELDS = ("batch", "injection_order", "injection_order_source", "injection_order_may_be_file_order")


def run_of(sample: SampleInfo) -> dict:
    """試料の測定の欄（バッチ・注入順とその出所）。統合が陰イオンの測定を残す形。"""
    return {"batch": sample.batch, "injection_order": sample.injection_order,
            "injection_order_source": sample.injection_order_source,
            "injection_order_may_be_file_order": bool(sample.injection_order_may_be_file_order)}


def order_may_be_file_order(sample: SampleInfo) -> bool:
    """注入 1 本の注入順が実際の測定順とは限らない: MS-DIAL の既定（ファイルの順）かもしれないか、試料の表が宣言
    `# injection_order_source = filename` でファイル名の順と書いた値（drift=auto が補正を止める条件・概観の
    ドリフトの兆候の但し書きが同じ規則で読む）。"""
    return bool(sample.injection_order_may_be_file_order) or sample.injection_order_source == "filename"


def polarity_samples(study: StudyData, polarity: str) -> list[SampleInfo]:
    """その極性の特徴の測定としての試料の表（`samples` と同じ順・同じ試料）。

    役割・群・個体・含めるか・属性は `samples` のまま（統合は陽陰で同じ試料を照合している）で、
    バッチと注入順だけをその極性の測定（`polarity_runs`）に替える。`polarity_runs` に無い極性は `samples`。
    """
    runs = study.polarity_runs.get(polarity)
    if runs is None:
        return list(study.samples)
    out = []
    for sample, run in zip(study.samples, runs):
        run = run or {}
        out.append(replace(sample, batch=_str_or_none(run.get("batch")),
                           injection_order=_float_or_none(run.get("injection_order")),
                           injection_order_source=run.get("injection_order_source"),
                           injection_order_may_be_file_order=bool(
                               run.get("injection_order_may_be_file_order"))))
    return out


def merged_feature_polarity(ds) -> dict[str, str]:
    """陽陰を統合したデータセットなら 特徴 ID → 極性（すべての特徴が統合の出自 `origin` を持つときだけ）。
    でなければ空（`quant.select.is_merged` と同じ見分け方）。"""
    metadata = getattr(ds, "feature_metadata", None) or {}
    ids = list(getattr(ds, "feature_ids", None) or [])
    out = {}
    for fid in ids:
        origin = (metadata.get(fid) or {}).get("origin")
        if not isinstance(origin, dict) or not isinstance(origin.get("polarity"), str):
            return {}
        out[str(fid)] = origin["polarity"]
    return out


def _polarity_runs(ds, polarities, display: list[str]) -> dict[str, list]:
    """統合したデータセットの、試料の表の極性（陽イオン）以外の極性の測定（`display` の順）。

    統合が `quant_state["polarity_runs"]` に残した値を試料名で引く。残っていない極性（その版より前の
    統合）は測定が分からない（None。注入順が無いのでドリフト補正はその極性に当たらない）。
    """
    recorded = (getattr(ds, "quant_state", None) or {}).get("polarity_runs") or {}
    out = {}
    for polarity in dict.fromkeys(polarities):
        if polarity == SAMPLE_TABLE_POLARITY:
            continue
        table = recorded.get(polarity) if isinstance(recorded.get(polarity), dict) else {}
        out[polarity] = [table.get(name) for name in display]
    return out


def sample_runs(ds) -> dict[str, dict]:
    """試料名（`ds.sample_names`）→ 統計エンジンがその試料の測定として読む値（`run_of`）。

    `from_dataset` と同じ試料の表から読む（試料の表の行 → .mddata → …。バッチは 2 種類以上あるときだけ）。
    統合（`quant.merge`）が陰イオンの測定を残すのに使う。
    """
    names = _engine_sample_ids(ds)
    samples, _candidates = _samples_table(ds, names)
    display = list(getattr(ds, "sample_names", None) or names)
    return {str(display[i]): run_of(s) for i, s in enumerate(samples)}


def quant_state_of(ds) -> QuantState | None:
    """`ds.quant_state`（lipidmix.quant-state.v1 の dict）から読む。無い・None なら None。"""
    raw = getattr(ds, "quant_state", None)
    if not isinstance(raw, dict):
        return None
    return QuantState(normalized=bool(raw.get("normalized", False)),
                      method=raw.get("method"), unit=raw.get("unit"))


def cell_states(values, detected) -> np.ndarray:
    """セルごとの状態（int8）。有限の正値だけが「値あり」で、それ以外は MISSING。"""
    v = np.asarray(values, dtype=float)
    present = np.isfinite(v) & (v > 0)
    if detected is None:
        return np.where(present, PRESENT_UNKNOWN, MISSING).astype(np.int8)
    d = np.asarray(detected, dtype=bool)
    out = np.full(v.shape, MISSING, dtype=np.int8)
    out[present & d] = DETECTED
    out[present & ~d] = GAP_FILLED
    return out


# 反復注入らしい試料名（plan §15）。茎が数字で終わる 2 段の番号、または inj/rep/tech 接尾辞。
_REPLICATE_A = re.compile(r"^(.*\d)[_\-.](?:inj|rep|tech|r)?\d{1,2}$")
_REPLICATE_B = re.compile(r"(?i)[_\-.](inj|rep|tech)\d+$")


def replicate_name_warnings(sample_ids: list[str]) -> list[dict]:
    """同じ茎の名前が 2 つ以上ある反復注入らしい試料を注意として返す（まとめはしない）。"""
    stems: dict[str, list[str]] = {}
    for name in sample_ids:
        m = _REPLICATE_A.match(name)
        if m:
            stem = m.group(1)
        else:
            m = _REPLICATE_B.search(name)
            if not m:
                continue
            stem = name[:m.start()]
        stems.setdefault(stem, []).append(name)
    return [{"code": "possible_technical_replicates", "stem": stem, "samples": names}
            for stem, names in stems.items() if len(names) >= 2]


def _str_or_none(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _float_or_none(value) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _include(row: dict) -> bool:
    """`include` の読み方は analysis_spec._load と同じ: 欄が無ければ真、あれば Python の真偽値。

    （文字列 "false" や "0" も真になる。sample_metadata_rows の include は取り込み時に bool へ
    正規化済みという前提で、ここで独自に解釈し直して _load と食い違わせない。）
    """
    return bool(row.get("include", True))


def _role(value) -> str | None:
    text = _str_or_none(value)
    return text.lower() if text else None


def _order_and_source(row: dict, ms: dict, sm: dict, mz_order=None) -> tuple[float | None, str | None]:
    """注入順と出所の名前。試料の表の行 → .mddata（msdial_samples）→ sample_meta → mzTab-M の MTD の順。

    試料の表が「宣言した作り方（raw の取得時刻など）では分からない」と印を付けた空欄（`sheet_order_unknown`）は
    不明のまま（None）にする: 表の順位と MS-DIAL の番号は付け方が違い、混ぜるとドリフト補正が順を取り違える。
    """
    from metabolomix.analysis.sample_manifest import (
        injection_order_source_label, order_provenance_label, sheet_order_unknown)
    value = _float_or_none(row.get("injection_order"))
    if value is not None:
        return value, order_provenance_label((row.get("provenance") or {}).get("injection_order"))
    if sheet_order_unknown(row):
        return None, None
    value = _float_or_none(ms.get("injection_order"))
    if value is not None:
        return value, "mddata"
    value = _float_or_none(sm.get("run_order"))
    if value is not None:
        return value, injection_order_source_label(sm.get("run_order_source"))
    value = _float_or_none(mz_order)
    if value is not None:
        return value, "mzTab"
    return None, None


def _row_batch(row: dict) -> str | None:
    """試料の表の行のバッチ。試料名の 8 桁から埋めた値（古い記録の `filename_date`）は採らない。"""
    prov = (row.get("provenance") or {}).get("batch") or {}
    if prov.get("source") == "filename_date":
        return None
    return _str_or_none(row.get("batch"))


def _samples_table(ds, names: list[str]) -> tuple[list[SampleInfo], dict[str, dict]]:
    """試料の表と、バッチの候補。sample_metadata_rows → msdial_samples → roles / sample_meta の順に欄ごと埋める。

    バッチは、試料の表の行（利用者のシート・解析フォルダの行）がそのまま勝ち、無ければ msdial_samples /
    sample_meta の宣言を「2 種類以上あるとき」だけ採る（定数のラベル、MS-DIAL の既定 "1" は情報を持たない。
    `sample_manifest.batch_labels_informative`）。それも無ければ不明（None）のまま。試料名の 8 桁
    （`filename_date_batch`）は候補として返すだけで、ドリフト補正・反復のまとめ・交絡の判定には当てない（A7）。

    注入順は値と出所（`sheet` / `mzTab` / `mddata` / …）を持つ。MS-DIAL 由来の注入順がファイル名の順と
    一致する試料には `injection_order_may_be_file_order` を立てる（MS-DIAL はファイルの順で振るのが既定。A6）。
    試料の表の行が無い（mzTab-M だけから読んだ）データセットは、mzTab-M の MTD（MS:4000088 のバッチ・
    MS:4000089 の注入順）も最後に読む（`dataset_status(view="samples")` の既定の行と同じ値）。
    """
    from metabolomix.analysis.sample_manifest import (
        MSDIAL_ORDER_SOURCES, _mztab_batch_and_order, batch_labels_informative, file_order_keys,
        filename_date_batch, matches_file_order)
    rows = list(getattr(ds, "sample_metadata_rows", None) or [])
    msdial = list(getattr(ds, "msdial_samples", None) or [])
    roles = getattr(ds, "roles", None) or {}
    meta = getattr(ds, "sample_meta", None) or {}
    by_name = {s.get("name"): s for s in msdial}
    display = list(getattr(ds, "sample_names", None) or [])
    plan = []
    for i, name in enumerate(names):
        # 安定 ID を先に、次に表示名（sample_names[i]）で引く。ID ≠ 表示名でも取りこぼさない。
        keys = [name] + ([display[i]] if i < len(display) and display[i] != name else [])
        row = rows[i] if i < len(rows) else {}
        ms = (next((by_name[k] for k in keys if k in by_name), None)
              or (msdial[i] if len(msdial) == len(names) else {}))
        sm = next((meta[k] for k in keys if meta.get(k)), {})
        mz_batch, mz_order = _mztab_batch_and_order(ds, i) if not rows else (None, None)
        role = (_role(row.get("role")) or _role(ms.get("role"))
                or next((r for k in keys if (r := _role(roles.get(k)))), None)
                or _role(sm.get("role")))
        if role is None:
            # どこにも役割が無い試料は表示名（無ければ ID）から推定する。MS-DIAL の Class ID
            # （msdial_samples の group）があれば sample_factors と同じく class_ids として渡す。
            label = keys[-1]
            class_id = _str_or_none(ms.get("group"))
            from metabolomix.analysis.preprocessing import detect_sample_roles
            role = _role(detect_sample_roles(
                [label], {label: class_id} if class_id else None).get(label)) or "sample"
        group = _str_or_none(row.get("group")) or _str_or_none(ms.get("group"))
        declared = (_str_or_none(ms.get("batch"))
                    or (_str_or_none(sm.get("batch")) if sm.get("batch_source") != "filename_date" else None)
                    or _str_or_none(mz_batch))
        order, order_source = _order_and_source(row, ms, sm, mz_order)
        # B-2 の dataset_set_sample_metadata は subject_id を行と sample_meta の両方に置く。
        subject = (_str_or_none(row.get("subject_id")) or _str_or_none(sm.get("subject_id"))
                   or _str_or_none(row.get("biological_sample_id")))
        # 属性は文字列にできる単一の値だけ（dict / list などの入れ物は因子にならないので落とす）
        attrs = {k: s for k, v in row.items()
                 if k not in _CORE_KEYS and k not in _BOOKKEEPING_KEYS
                 and not isinstance(v, (dict, list, tuple, set, frozenset))
                 and (s := _str_or_none(v)) is not None}
        plan.append((name, keys, _row_batch(row), declared, role, group, order, order_source,
                     subject, _include(row), attrs))
    use_declared = batch_labels_informative([p[3] for p in plan])
    # MS-DIAL 由来の注入順がファイルの順と一致するか（鍵は ds.sample_names の順に並ぶ）
    file_keys = file_order_keys(ds) if len(display) == len(names) else list(names)
    msdial_idx = [i for i, p in enumerate(plan) if p[6] is not None and p[7] in MSDIAL_ORDER_SOURCES]
    file_order = matches_file_order([plan[i][6] for i in msdial_idx],
                                    [file_keys[i] for i in msdial_idx])
    out, candidates = [], {}
    for i, (name, keys, row_batch, declared, role, group, order, order_source, subject, include,
            attrs) in enumerate(plan):
        batch = row_batch or (declared if use_declared else None)
        if batch is None:
            candidate = next((d for k in reversed(keys) if (d := filename_date_batch(k))), None)
            if candidate is not None:
                candidates[name] = {"value": candidate, "source": "filename_date"}
        out.append(SampleInfo(sample_id=name, role=role, group=group, batch=batch,
                              injection_order=order, subject_id=subject,
                              include=include, attributes=attrs,
                              injection_order_source=order_source,
                              injection_order_may_be_file_order=file_order and i in msdial_idx))
    return out, candidates


def _fingerprint(ds) -> str | None:
    try:
        from metabolomix.analysis.result_state import dataset_fingerprint
        return dataset_fingerprint(ds)
    except Exception:
        return None


def _resolve_sample_indices(display: list[str], ids: list[str], requested, what: str) -> dict[str, int]:
    """試料名（表示名）または StudyData の試料 ID → 行の位置。名前が先で、無ければ ID
    （`dataset_exclusions._resolve_samples` と同じ規則）。見つからないものは STUDY_DATA_INVALID。"""
    by_name = {n: i for i, n in reversed(list(enumerate(display)))}
    by_id = {n: i for i, n in reversed(list(enumerate(ids)))}
    out: dict[str, int] = {}
    unknown = []
    for item in requested:
        pos = by_name.get(item, by_id.get(item))
        if pos is None:
            unknown.append(item)
        else:
            out[item] = pos
    if unknown:
        raise DomainError("STUDY_DATA_INVALID", f"{what} に存在しない試料がある: {unknown[:20]}",
                          {"field": what, "unknown_samples": unknown})
    return out


def _engine_sample_ids(ds) -> list[str]:
    """エンジンの試料 ID の並び（試料の表の行があれば安定 ID、無ければ表示名）。"""
    rows = getattr(ds, "sample_metadata_rows", None)
    ids = list(getattr(ds, "sample_ids", None) or []) if rows else []
    return ids or list(ds.sample_names)


def _study_samples(ds, names: list[str], exclude_samples,
                   group_overrides) -> tuple[list[SampleInfo], dict[str, dict]]:
    """試料の表に除外（`include=False`）と群の上書きを当てたものと、バッチの候補（`from_dataset` の試料側）。"""
    samples, batch_candidates = _samples_table(ds, names)
    display = list(getattr(ds, "sample_names", None) or names)
    if exclude_samples:
        for pos in _resolve_sample_indices(display, names, list(exclude_samples),
                                           "exclude_samples").values():
            samples[pos] = replace(samples[pos], include=False)
    if group_overrides:
        for item, group in group_overrides.items():
            if not isinstance(group, str) or not group.strip():
                raise DomainError("STUDY_DATA_INVALID", "group_overrides の群は空でない文字列",
                                  {"field": "group_overrides", "sample": item, "value": repr(group)})
        resolved = _resolve_sample_indices(display, names, list(group_overrides),
                                           "group_overrides")
        for item, pos in resolved.items():
            samples[pos] = replace(samples[pos], group=group_overrides[item].strip())
    return samples, batch_candidates


def sample_table(ds, *, exclude_samples=(), group_overrides: dict[str, str] | None = None) -> list[SampleInfo]:
    """`from_dataset(ds, exclude_samples=..., group_overrides=...).samples` と同じ試料の表（値の行列は読まない）。

    試料の表だけが要る呼び手（群の上書きが試料ごとにどの群を当てるかの比べ合わせなど）向け。"""
    return _study_samples(ds, _engine_sample_ids(ds), exclude_samples, group_overrides)[0]


def from_dataset(ds, *, exclude_samples=(), exclude_features=(),
                 group_overrides: dict[str, str] | None = None) -> StudyData:
    """DatasetState から StudyData を作る。ds は書き換えない（配列も写す）。

    `exclude_samples`: 試料名か試料 ID。該当の試料は表に残したまま `include=False`（前処理の選択段が外す）。
    `exclude_features`: 特徴 ID。列ごと落とす。
    `group_overrides`: 試料名か試料 ID → 群。試料の表の群を置き換える。
    名前の解決は `dataset_exclusions._resolve_samples` と同じ（名前が先、無ければ ID）。
    存在しないものは STUDY_DATA_INVALID（何も作らない）。
    """
    values = np.array(np.asarray(ds.feature_matrix, dtype=float).T, copy=True)
    mask = getattr(ds, "detected_mask", None)
    if mask is not None and np.shape(mask) != values.T.shape:
        raise DomainError("STUDY_DATA_INVALID", "detected_mask の形が feature_matrix と合わない",
                          {"feature_matrix_shape": list(values.T.shape),
                           "detected_mask_shape": list(np.shape(mask))})
    detected = None if mask is None else np.array(np.asarray(mask, dtype=bool).T, copy=True)
    names = _engine_sample_ids(ds)
    classes, feature_names, feature_formulas = {}, {}, {}
    for fid, meta in (getattr(ds, "feature_metadata", None) or {}).items():
        ontology = (meta or {}).get("ontology")
        if ontology:
            classes[fid] = str(ontology)
        name = (meta or {}).get("name")
        if isinstance(name, str) and name.strip():
            feature_names[fid] = name
        formula = (meta or {}).get("formula")
        if isinstance(formula, str) and formula.strip():
            feature_formulas[fid] = formula.strip()
    annotated = frozenset(fid for fid, cands in (getattr(ds, "feature_candidates", None) or {}).items()
                          if cands)
    feature_ids = list(ds.feature_ids)
    samples, batch_candidates = _study_samples(ds, names, exclude_samples, group_overrides)
    display = list(getattr(ds, "sample_names", None) or names)
    feature_polarity = merged_feature_polarity(ds)
    if exclude_features:
        known = set(feature_ids)
        unknown = [f for f in exclude_features if f not in known]
        if unknown:
            raise DomainError("STUDY_DATA_INVALID",
                              f"exclude_features に存在しない特徴 ID がある: {unknown[:20]}",
                              {"field": "exclude_features", "unknown_features": unknown})
        drop = set(exclude_features)
        keep = [j for j, f in enumerate(feature_ids) if f not in drop]
        if not keep:
            raise DomainError("STUDY_DATA_INVALID", "exclude_features で特徴が 1 つも残らない",
                              {"field": "exclude_features"})
        values = values[:, keep]
        detected = None if detected is None else detected[:, keep]
        feature_ids = [feature_ids[j] for j in keep]
        classes = {f: c for f, c in classes.items() if f not in drop}
        feature_names = {f: n for f, n in feature_names.items() if f not in drop}
        feature_formulas = {f: v for f, v in feature_formulas.items() if f not in drop}
        annotated = frozenset(f for f in annotated if f not in drop)
        feature_polarity = {f: p for f, p in feature_polarity.items() if f not in drop}
    # 残った特徴の極性の測定だけ（陽イオンは試料の表そのもの）
    polarity_runs = (_polarity_runs(ds, feature_polarity.values(), [str(n) for n in display])
                     if feature_polarity else {})
    return StudyData(values=values, detected=detected, sample_ids=names,
                     feature_ids=feature_ids, samples=samples,
                     quant=quant_state_of(ds), dataset_hash=_fingerprint(ds),
                     feature_classes=classes, annotated_features=annotated,
                     feature_names=feature_names, feature_formulas=feature_formulas,
                     batch_candidates=batch_candidates, feature_polarity=feature_polarity,
                     polarity_runs=polarity_runs)


def merge_feature_columns(study: StudyData, groups: list[list[str]]) -> StudyData:
    """割れた特徴の組ごとに列を 1 本へまとめた StudyData の写し（`study` は変えない）。

    試料ごとに組の中の**大きい方**の値（NaN・0 以下は「無い」とみなし、有る値の最大。全部無ければ NaN）、
    `detected` は組の OR、残す ID は組の先頭（まとめた列は先頭の位置に置く）。行列に無い ID は飛ばし、
    残りが 1 つ以下ならまとめない。組どうしが特徴を共有すれば STUDY_DATA_INVALID。
    割れのまとめの規則は `alignment_qc.merge_matrix_columns`（`feature_search.merged_dataset` が使う）と
    同じ形で、「0 以下は無い」だけを足している。
    """
    ids = [str(f) for f in study.feature_ids]
    position = {fid: i for i, fid in enumerate(ids)}
    seen: set[str] = set()
    plans = []
    for group in groups:
        present = list(dict.fromkeys(str(f) for f in group if str(f) in position))
        if seen.intersection(present):
            raise DomainError("STUDY_DATA_INVALID", "割れの組が特徴を共有している",
                              {"shared": sorted(seen.intersection(present))})
        seen.update(present)
        if len(present) >= 2:
            plans.append(present)
    values = np.array(np.asarray(study.values, dtype=float), copy=True)
    det = None if study.detected is None else np.array(np.asarray(study.detected, bool), copy=True)
    drop: set[int] = set()
    classes = dict(study.feature_classes)
    names = dict(study.feature_names)
    formulas = dict(study.feature_formulas)
    annotated = set(study.annotated_features)
    merges = tuple(study.feature_merges) + tuple(tuple(p) for p in plans)
    for present in plans:
        cols = [position[f] for f in present]
        block = values[:, cols]
        block = np.where(np.isfinite(block) & (block > 0), block, np.nan)
        with np.errstate(all="ignore"):
            values[:, cols[0]] = np.fmax.reduce(block, axis=1)
        if det is not None:
            det[:, cols[0]] = det[:, cols].any(axis=1)
        drop.update(cols[1:])
        head = present[0]
        if head not in classes:
            donor = next((f for f in present[1:] if f in classes), None)
            if donor is not None:
                # クラスを譲った特徴の名前・組成式も一緒に譲る（脂質名の判定がクラスと同じ注釈を読む）
                classes[head] = classes[donor]
                for own in (names, formulas):
                    if donor in own:
                        own[head] = own[donor]
                    else:
                        own.pop(head, None)
        if any(f in annotated for f in present):
            annotated.add(head)
        for f in present[1:]:
            classes.pop(f, None)
            names.pop(f, None)
            formulas.pop(f, None)
            annotated.discard(f)
    keep = [i for i in range(len(ids)) if i not in drop]
    kept_ids = [study.feature_ids[i] for i in keep]
    # 統合したデータセットの極性は、まとめた列の先頭の特徴のもの
    polarity = {f: study.feature_polarity[f] for f in kept_ids if f in study.feature_polarity}
    return replace(study, values=values[:, keep],
                   detected=None if det is None else det[:, keep],
                   feature_ids=kept_ids,
                   feature_classes=classes, annotated_features=frozenset(annotated),
                   feature_names=names, feature_formulas=formulas, feature_merges=merges,
                   feature_polarity=polarity)
