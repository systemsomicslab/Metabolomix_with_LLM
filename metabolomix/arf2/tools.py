"""ARF2（MS-DIAL 全体カタログ）ツール群: arf2_parser, arf2_annotate_identities。

deps: mcp_core / session_state / path_resolvers / tool_helpers / arf2_reader /
lipid_identity。tools_* / server は import しない。
"""
from pathlib import Path

from metabolomix.msdial import lipid_identity
from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.core.path_resolvers import resolve_arf2_file_path
from metabolomix.core.tool_helpers import _identity_tables

__all__ = ["arf2_parser", "arf2_annotate_identities"]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf2_parser(file_path: str | None = None) -> str:
    """
    .arf2 ファイル（MS-DIALの全体カタログ）を解析し、データセットの全体像（メタデータ）を要約して返します。
    このファイルにはサンプル個別の強度データは含まれていないため、PCA等の多変量解析は実行できません。
    データ全体の品質や、アノテーション状況の概観を把握するために使用します。
    """
    file_path = resolve_arf2_file_path(file_path)
    if not file_path:
        return "データディレクトリに .arf2 ファイルが見つかりませんでした。"

    from metabolomix.arf2.reader import deserialize, generate_text_summary, summarize_arf2_data

    try:
        # ARF2データの読み込み
        with open(file_path, 'rb') as f:
            deserialized_data = deserialize(f)

        if not deserialized_data:
            return ".arf2 ファイルのパースに失敗したか、データが空です。"

        # 要約テキストの生成
        text_summary = generate_text_summary(deserialized_data)

        # 将来の検索やフィルタリング用に、カタログデータを ARF2 専用スロットへ保持する。
        # ARF（1スポット×1サンプル）とは粒度が違い、サンプル別強度を持たないため、
        # ARF の解析基盤（features / 前処理行列）を置き換えてはいけない。
        session_state.session.arf2.load(file_path, deserialized_data)

        output_text = (
            f"### 📂 ARF2 カタログデータのパース完了: {Path(file_path).name}\n"
            f"このファイルはデータセット全体の要約（平均値等）のみを含んでおり、サンプル別データを持たないためPCAは実行できません。\n\n"
            f"{text_summary}\n\n"
            f"※ 個別のサンプル比較やPCAを行いたい場合は、詳細データを持つ `.arf` (PeakProperties.arf など) を対象に `arf_parser` を使用してください。"
        )

        return session_state.session.maybe_prepend_caveat(output_text, topic="arf2")

    except Exception as e:
        return f"[ERROR] ARF2解析に失敗しました: {str(e)}"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf2_annotate_identities(file_path: str | None = None, max_rows: int = 50) -> str:
    """指定/自動解決の ARF2 スポット注釈を GOSLIN 正規化・RefMet/LIPID MAPS ID・
    MSI レベルで一括標準化し、TSV 表で返す（オフライン）。

    **返すのはファイル先頭から max_rows 件だけ**で、強度順でも MSI 順でもない。
    総スポット数と残り件数はヘッダ行に出るので、カタログ全体を見たと思わないこと。
    絞り込みには `arf_parser(annotation_keyword=...)` や `arf2_parser` の概観を使う。

    ARF2 には MS/MS 取得フラグ・精密質量誤差が無いため、MSI は保守的にクラス上限で
    評価する（`has_msms=False`, バンドは UNKNOWN）。より確度の高い MSI 評価は個別ピークの
    `verify_peak_annotation`（精密質量・アダクト整合を含むドシエ）を参照。
    """
    path = resolve_arf2_file_path(file_path)
    if not path:
        return json_payload({"status": "error", "message": ".arf2 が見つかりません。"})
    from metabolomix.arf2.reader import format_spots_as_table, load_catalog
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError
    spots = load_catalog(path)
    tables = _identity_tables()
    flag_note = None
    try:
        flag_set = curation_apply.flags_for_arf2(path)
    except FlagFileError as exc:
        # 一覧そのものは返す。curation_flag は空にし、壊れた記録を 1 行で名指しする。
        flag_set = {"wrong": set(), "suspect": set(), "assign": {}, "redundant": set()}
        flag_note = f"# curation_flag は空欄: {exc}"
    rows = []
    for spot in spots[:max_rows]:
        raw_name = spot.get("Name") or ""
        name = "" if raw_name.strip().lower() == "unknown" else raw_name
        feat = {"name": name, "ontology": spot.get("Ontology") or "",
                "has_msms": False}
        # ARF2 に MS/MS 取得フラグは無いため has_msms=False（MSI は保守的にクラス上限）
        block = lipid_identity.build_identity_block(
            feat, tables, mass_error_band="UNKNOWN", adduct_band="UNKNOWN")
        spot_id = spot.get("MasterAlignmentID")
        if spot_id in flag_set["wrong"]:
            curation_flag = "wrong"
        elif spot_id in flag_set["suspect"]:
            curation_flag = "suspect"
        elif spot_id in flag_set.get("assign", {}):
            curation_flag = f"assign:{flag_set['assign'][spot_id].get('name') or ''}"
        elif spot_id in flag_set.get("redundant", ()):
            curation_flag = "redundant"
        else:
            curation_flag = ""
        rows.append({"MasterAlignmentID": spot_id,
                     "name": feat["name"], "normalized": block["goslin"]["normalized"],
                     "refmet": block["reference"]["refmet_name"],
                     "lipid_maps_category": block["reference"]["lipid_maps_category"],
                     "msi_level": block["msi"]["level"],
                     "curation_flag": curation_flag})
    if not rows:
        return "ARF2 にスポットがありません。"
    # 列名を1回だけ出す TSV。同じ6キーを行数ぶん繰り返す JSON に対し、実データ
    # 50 行で 9,145 字 → 1,387 字（-85%）になる。
    header = (
        f"# ARF2 同定注釈（{Path(path).name}）\n"
        f"# 総スポット {len(spots)} 件のうち先頭 {len(rows)} 件（ファイル順・"
        f"強度順ではない）。残り {max(0, len(spots) - len(rows))} 件は未表示。\n"
        f"# MSI はクラス上限の保守的推定（ARF2 に MS/MS 取得フラグと質量誤差が無いため）。"
        + (f"\n{flag_note}" if flag_note else "")
    )
    columns = ["MasterAlignmentID", "name", "normalized", "refmet",
              "lipid_maps_category", "msi_level", "curation_flag"]
    return f"{header}\n{format_spots_as_table(rows, columns=columns)}"
