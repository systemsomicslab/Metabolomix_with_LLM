"""データセット入口ツール: list_data_files。

データセットを開く入口 `load_dataset` は `metabolomix.tools.handle_tools`。
フォルダは必ず `directory` で渡す（既定の探索先は無い）。
"""
import csv
from pathlib import Path

from metabolomix.core import path_resolvers
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp

__all__ = ["list_data_files"]

#: 群と取得順の重なりで一覧に載せる群の組の上限（戻り値を肥大させない）。
_PAIR_LIMIT = 30


def _sheet_groups(directory: Path) -> tuple[dict[str, str], str | None, int]:
    """raw のフォルダの実験情報シート（`sample-manifest.tsv`）の `source_file` → 群（生体試料の行だけ）。

    pipeline の受付と違って厳密には検証しない（群を読むだけ）。無い・読めないときは空と理由。群
    `unassigned`（`pipeline_sample_sheet` が因子の値が空・試料の表に行の無い試料に付ける。群ではない）の行は
    判定に入れず、その数を 3 つ目に返す。
    """
    from metabolomix.pipeline.inputs import DEFAULT_MANIFEST_NAME
    from metabolomix.rawmeta.sample_table import UNASSIGNED

    path = directory / DEFAULT_MANIFEST_NAME
    if not path.is_file():
        return {}, None, 0
    try:
        lines = [line for line in path.read_text(encoding="utf-8-sig").splitlines()
                 if not line.startswith("#")]
    except (OSError, UnicodeDecodeError) as exc:
        return {}, f"{DEFAULT_MANIFEST_NAME} を読めない（{type(exc).__name__}）", 0
    reader = csv.DictReader(lines, delimiter="\t")
    if not {"source_file", "group"} <= set(reader.fieldnames or []):
        return {}, f"{DEFAULT_MANIFEST_NAME} に source_file・group の列が無い", 0
    out, unassigned = {}, 0
    for row in reader:
        role = (row.get("role") or "").strip().lower()
        group = (row.get("group") or "").strip()
        if group == UNASSIGNED:
            unassigned += 1
        elif group and role in ("", "sample"):
            out[Path((row.get("source_file") or "").strip()).name] = group
    return out, None, unassigned


def _acquisition_lines(directory: Path, paths: list[str]) -> list[str]:
    """raw ごとの取得時刻（メタデータだけ。スペクトルは解さない）とその順、間隔、群との重なり。"""
    from metabolomix.console.job_manager import is_raw_input
    from metabolomix.rawmeta.acquisition import GAP_MINUTES, acquisition_order, gaps_over, group_overlap

    raws = [Path(p) for p in paths if is_raw_input(Path(p))]
    if not raws:
        return []
    rows = acquisition_order(raws)
    n_read = sum(r["rank"] is not None for r in rows)
    lines = ["\n## 取得時刻（raw のメタデータだけを読む。スペクトルは解さない）",
             f"読めた {n_read} / {len(rows)} 件。順は取得時刻の順（同じ時刻はファイル名）。"
             f"{GAP_MINUTES} 分を超える間隔: {gaps_over(rows)} 箇所",
             "file\tacquired_utc\torder\tsource"]
    for r in rows:
        source = r["source"] if r["rank"] is not None else f"unreadable（{r.get('error') or '理由不明'}）"
        lines.append("\t".join([r["file"], r["start_utc"] or "",
                                "" if r["rank"] is None else str(r["rank"]), source]))
    groups, problem, n_unassigned = _sheet_groups(directory)
    if problem:
        lines.append(f"群と取得順の重なり: 判定しない（{problem}）")
    elif groups:
        judged = group_overlap(rows, groups)
        lines.append(f"\n## 群と取得順の重なり（sample-manifest.tsv の群）: {judged['overlap']}")
        lines.append("complete = 取得順の範囲が重ならない（群と測定順が完全に交絡）/ independent = 重なり、"
                     "どの取得日（UTC）にも両群がいる / partial = それ以外 / unknown = 群のある raw に読めない時刻がある")
        if n_unassigned:
            lines.append(f"群 unassigned（因子の値が空・試料の表に行が無い）の {n_unassigned} 行は判定に入れない")
        lines.append("group_a\tgroup_b\toverlap")
        shown = judged["comparisons"][:_PAIR_LIMIT]
        lines.extend(f"{c['groups'][0]}\t{c['groups'][1]}\t{c['overlap']}" for c in shown)
        if len(judged["comparisons"]) > len(shown):
            lines.append(f"（ほか {len(judged['comparisons']) - len(shown)} 組）")
        ranges = "、".join(f"{g}: {lo}–{hi}" for g, (lo, hi) in judged["group_rank_ranges"].items())
        if ranges:
            lines.append(f"群ごとの取得順の範囲: {ranges}")
    return lines


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def list_data_files(
    directory: str,
    extension: str | None = None,
    all_files: bool = False,
) -> str:
    """データフォルダにある解析対象ファイルを拡張子別に一覧します。

    - directory: 探索するフォルダ（必須）。
    - extension: 指定するとその拡張子だけに絞る（例 '.pai2', '.arf2', '.EIC.aef'）。
    - all_files: 既定では解析できる拡張子
      (.arf / .arf2 / .pai2 / .dcl / .EIC.aef / .mddata / .mdproject) だけを返す。
      測定生データ（MS-DIAL が読む .abf / .ibf / .cdf / .mzml / .wiff / .raw /
      .d / .wiff2 / .qgd / .lcd / .lrp / .imzml）まで見たいときだけ True にする。
      Agilent・Bruker の .d と Waters の .raw は**フォルダ**が 1 検体なので、
      末尾に `/` を付けて示す。

    一覧に測定生データが入るとき（all_files=True か、生データの拡張子を extension に渡したとき）は、
    raw ごとの取得時刻（メタデータだけを読む: .wiff の OLE・Thermo .raw のファイルヘッダ・Bruker .d の
    tdf/sqlite・mzML の startTimeStamp。スペクトルは解さない）と、そこからの順・40 分を超える間隔の数を
    足す。フォルダに実験情報シート `sample-manifest.tsv` があれば、その群と取得順の重なり（群の組ごとに
    complete / independent / partial / unknown。群 `unassigned` の行は入れない）も返す。

    返すのは絶対パスで、`dcl_find_msms(file_path=...)` 等へそのまま渡せる。
    フォルダ全体をまず解析するなら `load_dataset(directory=...)` が入口。
    """
    target_dir = Path(directory).expanduser()
    if not target_dir.exists():
        return f"データディレクトリが存在しません: {target_dir}"
    if not target_dir.is_dir():
        return f"指定されたパスはディレクトリではありません: {target_dir}"

    paths = path_resolvers.list_data_files(
        directory=directory, extension=extension, all_files=all_files)
    if not paths:
        scope = extension or ("全ファイル" if all_files else "解析対象の拡張子")
        return (
            f"条件に一致するファイルがありません（ディレクトリ: {target_dir}, 対象: {scope}）。"
            "all_files=True で解析対象外のファイルも確認できます。"
        )

    # 拡張子ごとにまとめる。同じフォルダのパスが延々と並ぶより、何が何件あるかが
    # 先に見えるほうが次の一手（load_dataset か個別パーサか）を選びやすい。
    groups: dict[str, list[str]] = {}
    for path in sorted(paths):
        suffix = next(
            (s for s in path_resolvers.ANALYSABLE_EXTENSIONS if path.endswith(s)),
            Path(path).suffix or "(拡張子なし)",
        )
        groups.setdefault(suffix, []).append(path)

    lines = [f"データディレクトリ: {target_dir}", f"該当ファイル: {len(paths)} 件"]
    for suffix, members in sorted(groups.items()):
        lines.append(f"\n## {suffix}（{len(members)} 件）")
        # フォルダ形式の計測データ（.d / Waters の .raw）は末尾の `/` で示す。
        # 見分けが付かないと「開けないファイル」として扱われる。
        lines.extend(m + "/" if Path(m).is_dir() else m for m in members)
    lines.extend(_acquisition_lines(target_dir, sorted(paths)))
    return "\n".join(lines)
