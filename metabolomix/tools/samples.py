"""サンプル検索: 因子トークンからサンプル・FileID・実ファイルパスを引く。

ARF に限らず pai2 / dcl / EIC が使える形で条件検索するための入口。ARF のロード前
（＝どのファイルを開くか決める前）でも呼べるよう、ロード済み ARF が無ければ
ディレクトリの .mddata、それも無ければ実ファイル名からサンプル集合を組み立てる。

deps: mcp_core / session_state / path_resolvers / sample_factors / metabolomix.msdial.classes /
metabolomix.msdial.tags。tools_* / server は import しない（循環回避）。
"""
from pathlib import Path

from metabolomix.core import mcp_core
from metabolomix.core import path_resolvers
from metabolomix.msdial import sample_factors
from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.msdial.classes import discover_arf_class_index
from metabolomix.msdial.tags import normalize_sample_name

__all__ = ["sample_search"]

# 既定は「1測定ファイルにつき1個」の per-sample ファイルだけ。.arf / .arf2 /
# .EIC.aef はアラインメント単位でサンプルに1対1対応しないため既定に含めない。
DEFAULT_EXTENSIONS = (".pai2", ".dcl")


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def sample_search(
    specs: list[str] | None = None,
    directory: str | None = None,
    extensions: list[str] | None = None,
    include_roles: list[str] | None = None,
) -> str:
    """因子トークンでサンプルを検索し、FileID と実ファイルパスを返す。

    MS-DIAL の Class ID には実験デザインの全因子が入らない（時点・複製・測定日は
    サンプル名にしか無いことがある）。本ツールは Class ID とサンプル名を統合した
    トークン空間で検索し、pai2 / dcl / EIC が直接使える形で結果を返す。

    - specs: `_` 区切りの因子トークン指定のリスト（要素内 AND・順不同、要素間 OR）。
      例 `["ILG_6h"]`, `["ILG_6h","control_6h"]`。**省略すると全サンプルと
      `token_vocabulary`（何で絞れるかの語彙一覧）を返す**。
    - directory: 探索先。省略時は既定のデータディレクトリ。ロード済み ARF があれば
      そのサンプル名を優先して使う。
    - extensions: 実パスを引く拡張子（既定 `[".pai2", ".dcl"]`）。複数バッチが
      混在していても最新バッチを自動選択する。
    - include_roles: 結果に含める role。**省略時は全 role を返す**（`role` フィールドで
      判断できるため検索側では絞らない）。`["sample"]` で QC/blank を落とせる。

    返り値の各サンプルは `name` / `file_id` / `class_id` / `role` / `tokens` /
    `files{拡張子: 実パス}`。`file_id` は `eic_plot_chromatograms(file_ids=[...])`、
    `files[".pai2"]` は `pai2_parser(file_path=...)` にそのまま渡せる。
    """
    try:
        facets, source, target_dir = _collect_facets(directory)
    except (FileNotFoundError, NotADirectoryError) as exc:
        return _error(str(exc))
    if not facets:
        return _error(
            f"サンプルを特定できませんでした（探索先: {target_dir}）。"
            "MS-DIAL の出力フォルダ（.mddata か .pai2 を含む）を directory に指定するか、"
            "先に arf_parser で ARF を読み込んでください。")

    exts = [str(e) for e in extensions] if extensions else list(DEFAULT_EXTENSIONS)
    roles = tuple(include_roles) if include_roles else None
    vocabulary = sample_factors.token_vocabulary(facets)

    cleaned = [str(s) for s in specs or [] if str(s).strip()]
    if not cleaned:
        kept, dropped = _apply_role_filter(facets, roles)
        return json_payload({
            "status": "success",
            "source": source,
            "directory": str(target_dir) if target_dir else None,
            "specs": [],
            "total_samples": len(facets),
            "matched": len(kept),
            "excluded_by_role": sorted(dropped),
            "samples": [_describe(f, target_dir, exts) for f in kept.values()],
            "token_vocabulary": sample_factors.token_vocabulary(kept),
        })

    try:
        matches, excluded = sample_factors.expand_sample_specs(
            cleaned, facets, include_roles=roles)
    except ValueError as exc:
        return json_payload({
            "status": "error",
            "message": str(exc),
            "token_vocabulary": vocabulary,
        })

    selected = {n for names in matches.values() for n in names}
    ordered = [f for f in facets.values() if f.name in selected]
    return json_payload({
        "status": "success",
        "source": source,
        "directory": str(target_dir) if target_dir else None,
        "specs": cleaned,
        "resolved": matches,
        "total_samples": len(facets),
        "matched": len(ordered),
        "excluded_by_role": sorted({n for names in excluded.values() for n in names}),
        "samples": [_describe(f, target_dir, exts) for f in ordered],
    })


def _error(message: str) -> str:
    return json_payload({"status": "error", "message": message})


def _apply_role_filter(facets, roles):
    """specs 省略時（語彙モード）にも include_roles を適用する。

    specs 指定時は sample_factors.expand_sample_specs が role フィルタを担うが、
    specs 省略の分岐はそこを通らないため facets 自体にここで同じ意味論のフィルタを
    掛ける。roles が None（省略時）なら従来どおり全 role を素通しする。
    戻り値は (残った {サンプル名: SampleFacet}, 除外されたサンプル名の集合)。
    """
    if roles is None:
        return dict(facets), set()
    allowed = {str(role).casefold() for role in roles}
    kept = {name: f for name, f in facets.items() if f.role.casefold() in allowed}
    dropped = {name for name in facets if name not in kept}
    return kept, dropped


def _collect_facets(directory):
    """(facets, source, target_dir) を返す。source は facets の出所の説明。

    ロード済み ARF > ディレクトリの .mddata > ディレクトリの per-sample ファイル名、
    の順に試す。ARF ロード前でも呼べることが本ツールの要件。
    """
    session = session_state.session
    if directory is None and _has_loaded_arf(session):
        # 検索対象は**常にデータセット全体**（features）。filtered_features を優先すると
        # 直前の arf_parser(class_ids=...) の絞り込みで検索空間が黙って痩せ、「どの
        # サンプルが存在するか」を尋ねる本ツールが取りこぼしを起こす。
        names = sample_factors.arf_sample_names(session.arf.features)
        if names:
            facets = sample_factors.build_sample_facets(names, session.arf.class_index)
            return facets, "loaded_arf", _loaded_arf_dir(session)

    target_dir = Path(directory).expanduser() if directory else mcp_core.DATA_DIR
    if not target_dir.exists():
        raise FileNotFoundError(f"データディレクトリが存在しません: {target_dir}")
    if not target_dir.is_dir():
        raise NotADirectoryError(f"指定されたパスはディレクトリではありません: {target_dir}")

    class_index = discover_arf_class_index(target_dir)
    if class_index and class_index.get("records"):
        names = [str(r["file_name"]) for r in class_index["records"] if r.get("file_name")]
        return (sample_factors.build_sample_facets(names, class_index),
                "mddata", target_dir)

    names = _names_from_files(target_dir)
    return sample_factors.build_sample_facets(names, None), "filenames", target_dir


def _has_loaded_arf(session) -> bool:
    return (
        session.arf.features is not None
        and str(session.arf.current_file_path or "").lower().endswith(".arf")
    )


def _loaded_arf_dir(session):
    path = session.arf.current_file_path
    parent = Path(str(path)).parent
    return parent if parent.is_dir() else None


def _names_from_files(target_dir: Path) -> list[str]:
    """per-sample ファイル名から、処理タイムスタンプを除いたサンプル名を復元する。

    .mddata が無いフォルダ向けのフォールバック。表示名は「最新バッチのファイル名から
    タイムスタンプを取り除いた形」にそろえる（MS-DIAL の AnalysisFileName と同じ形）。
    """
    candidates = [str(p) for p in sorted(target_dir.iterdir())
                  if p.is_file() and p.suffix.casefold() in {".pai2", ".dcl"}]
    names: list[str] = []
    seen: set[str] = set()
    for path in path_resolvers._select_latest_batch(candidates):
        stem = Path(path).stem
        key = normalize_sample_name(stem, strip_processing_timestamp=True)
        if not key or key in seen:
            continue
        seen.add(key)
        # 元の大文字小文字を保つため、正規化キーではなくタイムスタンプ除去後の実表記を使う。
        names.append(stem[:-13] if len(stem) > 13 and stem[-13] == "_" and stem[-12:].isdigit() else stem)
    return names


def _describe(facet, target_dir, extensions) -> dict:
    return {
        "name": facet.name,
        "file_id": facet.file_id,
        "class_id": facet.class_id,
        "role": facet.role,
        "tokens": sorted(facet.tokens),
        "files": _resolve_files(facet.name, target_dir, extensions),
    }


def _resolve_files(sample_name, target_dir, extensions) -> dict:
    """サンプル名に対応する実ファイルを拡張子ごとに1つ解決する（最新バッチ優先）。"""
    if target_dir is None or not Path(target_dir).is_dir():
        return {}
    key = normalize_sample_name(sample_name, strip_processing_timestamp=True)
    resolved: dict[str, str] = {}
    for ext in extensions:
        candidates = []
        for p in sorted(Path(target_dir).iterdir()):
            if not p.is_file() or not p.name.casefold().endswith(ext.casefold()):
                continue
            # normalize_sample_name の既知拡張子リストに無い拡張子（.dcl 等）でも
            # 正しく揃うよう、比較対象の拡張子をここで明示的に切り落としてから渡す。
            stem = p.name[:-len(ext)] if ext else p.name
            if normalize_sample_name(stem, strip_processing_timestamp=True) == key:
                candidates.append(str(p))
        if not candidates:
            continue
        latest = path_resolvers._pick_latest(path_resolvers._select_latest_batch(candidates))
        if latest:
            resolved[ext] = latest
    return resolved
