"""サンプル単位の因子トークン（サンプル名 ∪ Class ID）による選択・群分け（純ロジック層）。

MS-DIAL の Class ID は「ユーザーが MS-DIAL 上で入力した1文字列」でしかなく、実験
デザインの全因子を含むとは限らない（時点・複製・測定日はサンプル名にしか無いことが
ある）。そのため本モジュールは Class ID とサンプル名のトークンを統合した空間で
spec を解決し、metabolomix.msdial.classes の Class ID 専用ロジックを一般化する。

依存は metabolomix.msdial.tags（normalize_sample_name）と preprocessing（detect_sample_roles）
のみの leaf。metabolomix.msdial.classes / tools_* / session_state / server は import しない
（metabolomix.msdial.classes → sample_factors の向きに依存させるため。逆向きは循環になる）。
"""
from __future__ import annotations

from dataclasses import dataclass

from metabolomix.analysis import preprocessing
from metabolomix.msdial.tags import normalize_sample_name


@dataclass(frozen=True)
class SampleFacet:
    """1サンプルの選択・群分けに必要な情報一式。

    tokens は casefold 済みの統合トークン集合で、spec 照合の唯一の入力。
    file_id は EIC の file_ids にそのまま渡せる MS-DIAL AnalysisFileId。
    """

    name: str
    file_id: int | None
    class_id: str | None
    role: str
    tokens: frozenset[str]


def split_tokens(value) -> frozenset[str]:
    """`_` 区切りの因子トークン集合（casefold）。空要素は落とす。"""
    if value is None:
        return frozenset()
    return frozenset(token for token in str(value).casefold().split("_") if token)


def sample_tokens(sample_name, class_id=None, extra=None) -> frozenset[str]:
    """サンプル名 ∪ Class ID（∪ 追加ラベル）の統合トークン集合を返す。

    サンプル名は normalize_sample_name で既知の測定ファイル拡張子と末尾12桁の処理
    タイムスタンプを落としてから分割する。処理タイムスタンプは MS-DIAL の再処理
    ごとに変わる識別子で実験因子ではないため、トークン語彙に混ぜない。
    extra は Class ID 以外の群ラベル（session_state の sample_meta["group"] 等）用。
    """
    tokens = set(split_tokens(normalize_sample_name(sample_name, strip_processing_timestamp=True)))
    tokens |= split_tokens(class_id)
    tokens |= split_tokens(extra)
    return frozenset(tokens)


def build_sample_facets(sample_names, class_index=None, sample_meta=None) -> dict[str, SampleFacet]:
    """サンプル名リストから {サンプル名: SampleFacet} を作る（入力順を保持）。

    - class_index: metabolomix.msdial.classes.discover_arf_class_index の戻り。あれば file_id /
      class_id を名前一致で解決する。**None でも成立**し、その場合はサンプル名の
      トークンだけで選択できる（.mddata が無いフォルダでも因子指定が効く）。
    - sample_meta: session_state.session.arf.sample_meta 相当。role と group ラベルの
      供給源。arf_differential は class_index を持たず sample_meta["group"] だけを
      持つ経路があるため、group もトークン源として合流させる。
    - role: sample_meta に明示があればそれを優先し、無ければ
      preprocessing.detect_sample_roles（名前と Class ID のトークン照合）で決める。
    """
    lookup = _class_lookup(class_index)
    meta = sample_meta or {}
    ordered = list(sample_names)

    records: dict[str, dict] = {}
    class_ids: dict[str, str] = {}
    for name in ordered:
        record = lookup.get(normalize_sample_name(name, strip_processing_timestamp=False)) or {}
        records[name] = record
        if record.get("class_id"):
            class_ids[name] = record["class_id"]

    detected = preprocessing.detect_sample_roles(ordered, class_ids)

    facets: dict[str, SampleFacet] = {}
    for name in ordered:
        record = records[name]
        entry = meta.get(name) or {}
        facets[name] = SampleFacet(
            name=name,
            file_id=record.get("file_id"),
            class_id=record.get("class_id"),
            role=entry.get("role") or detected.get(name, "sample"),
            tokens=sample_tokens(name, record.get("class_id"), entry.get("group")),
        )
    return facets


def expand_sample_specs(specs, facets, *, include_roles=("sample",)):
    """各 spec を該当サンプル名のリストへ展開する。

    spec は `_` 区切りのトークン列で、**その全トークンを含む**サンプルに一致する
    （要素内 AND・順不同）。リスト内の複数 spec は互いに独立（呼び出し側で OR 合算
    する）。完全なサンプル名や完全な Class ID を渡しても同じ規則で解決される。

    include_roles: 既定 ("sample",) で QC/blank を落とす。統合トークン空間により
    `class_ids=["cerebellum"]` が `20240311_QC_Cerebellum_...` を名前経由で掴むように
    なったため、群平均・PCA の汚染を既定で防ぐ。落とした分は excluded に残して
    呼び出し側が開示できるようにする。None を渡すと role による絞り込みをしない。

    戻り値 (matches, excluded)。matches={spec: [サンプル名, ...]}（facets の順）。
    role 適用後に1件も残らない spec があれば ValueError（利用可能トークンを添える）。
    """
    allowed = None if include_roles is None else {str(role).casefold() for role in include_roles}
    matches: dict[str, list[str]] = {}
    excluded: dict[str, list[str]] = {}
    for spec in specs:
        wanted = split_tokens(str(spec).strip())
        if not wanted:
            continue
        hits = [facet.name for facet in facets.values() if wanted <= facet.tokens]
        if allowed is None:
            kept, dropped = hits, []
        else:
            kept = [n for n in hits if facets[n].role.casefold() in allowed]
            dropped = [n for n in hits if facets[n].role.casefold() not in allowed]
        if not kept:
            raise ValueError(_no_match_message(spec, facets, dropped))
        matches[spec] = kept
        excluded[spec] = dropped
    return matches, excluded


def _no_match_message(spec, facets, dropped) -> str:
    """一致ゼロの spec に対する説明文。role で全滅した場合はその旨を明示する。

    'matched' の語を含めるのは、既存テスト（expand_class_specs 由来）が
    assertRaisesRegex(ValueError, "matched") で拾っているため。
    """
    if dropped:
        roles = sorted({facets[n].role for n in dropped})
        return (
            f"No sample matched spec '{spec}' after role filtering: "
            f"{len(dropped)} 件が role={'/'.join(roles)} のため除外されました。"
            f"含めるには include_roles=[\"sample\", \"{roles[0]}\"] を指定してください。"
        )
    available = sorted({token for facet in facets.values() for token in facet.tokens})
    return (
        f"No sample matched spec '{spec}'. "
        f"Available tokens: {', '.join(available) if available else 'なし'}"
    )


def assign_factor_groups(facets, group_factors=None, group_levels=None, sep="|"):
    """因子軸ごとに値 spec を解決し、直積ラベルを {サンプル名: ラベル} で返す。

    - group_factors: 因子軸のリスト。各軸は値 spec のリスト（例
      [["control", "ILG"], ["0h", "6h"]] → "control|0h" 等）。値 spec は多トークン可
      （"G_uralensis"）。位置インデックスを使わないのは、G_uralensis のように値が2
      トークンに割れると以降の位置がずれて壊れるため。
    - group_levels: 1軸のときの糖衣で group_factors=[group_levels] と等価。1軸では
      連結が起きないため、従来の出力文字列と完全に一致する（後方互換）。
      group_factors が指定されていればそちらが優先。
    - 両方未指定なら各サンプルの Class ID（無ければ None）をそのままラベルにする。
    - 軸内で2つ以上の値に一致したら ValueError（値が相互排他でない＝指定ミス）。
    - 軸内でどの値にも一致しなければその軸は "other"（除外はしない）。
    - role が sample でないサンプルは直積ではなく role 名（"qc"/"blank"）をラベルに
      する。除外しないのは QC の凝集が前処理品質の判断材料になるため。
    """
    if group_factors:
        axes = [[str(v).strip() for v in axis if str(v).strip()] for axis in group_factors]
        axes = [axis for axis in axes if axis]
    elif group_levels:
        axes = [[str(v).strip() for v in group_levels if str(v).strip()]]
        axes = [axis for axis in axes if axis]
    else:
        axes = []

    groups: dict[str, str | None] = {}
    for name, facet in facets.items():
        if not axes:
            groups[name] = facet.class_id
            continue
        if facet.role != "sample":
            groups[name] = facet.role
            continue
        parts = []
        for axis in axes:
            hits = [value for value in axis if split_tokens(value) <= facet.tokens]
            if len(hits) > 1:
                raise ValueError(
                    f"Sample '{name}' matches multiple values in one factor axis "
                    f"({', '.join(hits)}); values within an axis must be mutually exclusive."
                )
            parts.append(hits[0] if hits else "other")
        groups[name] = sep.join(parts)
    return groups


def token_vocabulary(facets) -> dict:
    """トークン語彙を集計する（何で絞れるかを発見するための一覧）。

    戻り値:
      {"tokens": {token: {"samples": 出現サンプル数,
                          "roles": {role: 件数},
                          "positions": [サンプル名を `_` 分割したときの出現位置, ...]}},
       "by_position": {"0": [token, ...], ...}}

    positions / by_position は因子の並びを推測する手掛かりだが、G_uralensis のように
    値が2トークンに割れると以降がずれるため、**フィルタ指定には使わない**（指定は
    常に値トークンで行う）。Class ID にしか無いトークンは positions が空になる。
    """
    counts: dict[str, dict] = {}
    by_position: dict[int, set[str]] = {}

    for facet in facets.values():
        for token in facet.tokens:
            entry = counts.setdefault(token, {"samples": 0, "roles": {}, "positions": set()})
            entry["samples"] += 1
            entry["roles"][facet.role] = entry["roles"].get(facet.role, 0) + 1
        normalized = normalize_sample_name(facet.name, strip_processing_timestamp=True)
        for position, token in enumerate(normalized.split("_")):
            if not token:
                continue
            by_position.setdefault(position, set()).add(token)
            counts.setdefault(
                token, {"samples": 0, "roles": {}, "positions": set()},
            )["positions"].add(position)

    return {
        "tokens": {
            token: {
                "samples": entry["samples"],
                "roles": dict(sorted(entry["roles"].items())),
                "positions": sorted(entry["positions"]),
            }
            for token, entry in sorted(counts.items())
        },
        "by_position": {str(position): sorted(tokens) for position, tokens in sorted(by_position.items())},
    }


def arf_sample_names(features) -> list[str]:
    """ARF スポット列の AlignedPeakProperties 行に現れるサンプル名を出現順で返す。

    行は MS-DIAL の生 list で index 1 が FileName（AlignmentChromPeakFeature スキーマ）。
    arf_reader を経由せず素の index 参照で済ませ、leaf の依存を増やさない。
    """
    names: list[str] = []
    seen: set[str] = set()
    for spot in features or []:
        for row in (spot or {}).get("AlignedPeakProperties") or []:
            if not isinstance(row, list) or len(row) < 2:
                continue
            value = row[1]
            if isinstance(value, bytes):
                value = value.decode("utf-8", errors="ignore")
            if isinstance(value, str) and value and value not in seen:
                seen.add(value)
                names.append(value)
    return names


def _class_lookup(class_index) -> dict[str, dict]:
    """class_index の records を正規化サンプル名で引ける辞書にする。

    metabolomix.msdial.classes.resolve_sample_class を使わないのは、metabolomix.msdial.classes が本モジュール
    を import する側であり、逆向きの import が循環になるため（session_state.
    _build_sample_meta も同じ理由で同じ引き方をしている）。
    """
    if not class_index:
        return {}
    lookup: dict[str, dict] = {}
    for record in class_index.get("records", []):
        key = normalize_sample_name(record.get("file_name"), strip_processing_timestamp=False)
        if key:
            lookup[key] = record
    return lookup
