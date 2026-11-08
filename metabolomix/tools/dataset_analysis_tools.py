"""DatasetState に対する解析 MCP ツール。

load_dataset → dataset_preprocess → dataset_compare
→ dataset_export_differential の順に実行する（PCA の `dataset_pca` は `stats_tools`）。解析フォルダも mzTab-M も
`load_dataset` が同じ DatasetState にするので、同じ入力からは同じ数字が出る。dataset_feature_search は読み取り専用で、特徴の検索・目的の物質の
判定に差次的結果の統計を添えて返す（判定の純ロジックは analysis/feature_search.py）。

戻り値には要約だけを載せ、全量（差次的の results/volcano）は
session.dataset に保持する（CLAUDE.md の戻り値肥大禁止）。
"""
from __future__ import annotations

import math
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.analysis import export_contract, stats_results
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.tools._holder import resolve_holder
# 極性の決め方は mztab/polarity.py が唯一の実装（ハンドルの読み込みと共有）。
from metabolomix.mztab.polarity import dataset_polarity as _dataset_polarity

__all__ = [
    "dataset_export_differential",
    "dataset_export_evidence_bundle",
    "dataset_set_sample_metadata",
    "dataset_feature_search",
]

# 欠けている状態 → それを作れるツール。missing_state の required_tools になる。
# **自分自身は入れない**（クライアントが同じ呼び出しを繰り返すループになる）。
_RECOVERY_TOOLS = {
    "dataset": ["load_dataset"],
    "dataset_compare_result": ["dataset_compare"],
}


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True),
    structured_output=False)
def dataset_export_differential(output_path: str,
                                result_id: str | None = None,
                                apply_curation: bool = True,
                                dataset_handle: str | None = None,
                                contrast_id: str | None = None) -> str:
    """指定した差次的結果を InChIKey 付きの 1 ファイルへ書き出す。

    先に dataset_preprocess → dataset_compare を実行しておくこと。
    出力は差次的エクスポートの契約（`analysis/export_contract.py`。15 列 + contract_version
    メタ行）に従うので、下流のパスウェイ解析にそのまま渡せる。

    apply_curation（既定 true）: curation_submit で wrong を付けたスポットを
    同定なしとして扱い、出力から外す。フラグがあればメタ行 # curation = ... で
    適用状況を宣言する。フラグが無ければ出力は変わらない。成功 payload の
    `curation` は `{state, wrong_excluded, suspect, orphaned}`。assign（curation_submit の候補付け）は同定を置き換え（name_source /
    inchikey_source = curation）、redundant は除外する。どちらかがあれば `curation` に
    `assigned` と `redundant_excluded` が足される。
    フラグ記録に読めない行があれば書き出さずに `CURATION_FLAGS_INVALID`。

    result_id: 書き出す結果を名指しする（省略時は直近の差次的結果）。
        **前処理をやり直した後の古い結果は書き出しません**（`STALE_ANALYSIS_RESULT`）。
        古い数字に現在の前処理条件のラベルが付いた TSV は、どちらも正しく見えて
        ずれが分からなくなります。

    InChIKey は DatasetState.feature_metadata から取る（mzTab-M の
    database_identifier / InChI / SMILES 由来。.arf2 との結合は不要）。
    濃縮解析の背景を保つため、有意な行だけでなく InChIKey が付いた全行を出す。
    ontology と msi_level は mzTab-M に対応物が無いため空欄で、その旨をメタ行に
    書く（空欄を「該当なし」と読み違えさせない）。メタ行には result_id と
    preprocess_id も入るので、あとから「どの計算の出力か」を辿れる。

    contrast_id: `dataset_compare` の結果が対比を複数持つときに書き出す対比の ID（`T_vs_C` など）。
        対比が 1 つなら省略できる。複数で省くか未知の ID なら `CONTRAST_REQUIRED`、交互作用の対比は
        2 群の形にならないので `CONTRAST_NOT_TWO_GROUP`（どちらも details に候補）。
        統計エンジンの結果には、15 列の後ろ（`msi_note` の後）に検定・判定規則・単位などのメタ行が
        付く（列と契約の版は同じ。`mean_a` / `mean_b` は前処理後の線形値の算術平均）。

    dataset_handle: `load_dataset` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に load_dataset を実行してください。")

    last = ds.last_differential
    if result_id is not None:
        try:
            result_id, contrast_id = stats_results.split_result_ref(result_id, contrast_id)
        except DomainError as exc:
            return mztab_error(exc.code, exc.message, exc.details or None)
        last = (ds.results or {}).get(result_id)
    if not last or last.get("kind") not in ("two_group", "compare"):
        return _missing(
            "dataset_compare_result",
            "先に dataset_preprocess → dataset_compare を実行してください。"
            if result_id is None else
            f"result_id={result_id} の 2 群比較の結果がありません。")

    from metabolomix.analysis.dataset_export import export_dataset_result
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError, orphaned_warning

    curation = None
    flag_set = None
    lookup = curation_apply.arf2_lookup_kwargs(ds)
    arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p, **lookup) for p in ds.source_files)), None)
    if arf2 is not None:
        try:
            flag_set = curation_apply.flags_for_arf2(arf2)
        except FlagFileError as exc:
            return mztab_error("CURATION_FLAGS_INVALID", str(exc), exc.details())
        if flag_set["n"]:
            # SMF_ID と MasterAlignmentID が同じ空間だと言えるのは、隣接 .arf との
            # 対応が検証済み（feature_axis = "arf_spot_index == mztab_smf_id"）のときだけ。
            mapped = (getattr(ds, "feature_qc", None) or {}).get("source") == "arf"
            state = ("not_applied" if not apply_curation
                     else "applied" if mapped else "unmapped")
            curation = {"state": state, "flag_set": flag_set}

    try:
        info = export_dataset_result(ds, last, Path(output_path), curation=curation,
                                     contrast_id=contrast_id)
    except DomainError as exc:
        if exc.code in ("STALE_ANALYSIS_RESULT", "ANALYSIS_RESULT_NOT_FOUND"):
            # 再計算すれば直る。どのツールを呼べばよいかを機械可読に返す。
            return _missing("dataset_compare_result", exc.message)
        return mztab_error(exc.code, exc.message, exc.details or None)

    curation_stats = info.pop("curation_stats", None)
    info["curation"] = curation_apply.payload_summary(info["curation"], flag_set, curation_stats)
    orphaned = (flag_set or {}).get("orphaned", 0)
    return json_payload({
        "status": "success",
        **info,
        **({"warnings": [orphaned_warning(orphaned)]} if orphaned else {}),
        "log2fc_sign": "log2fc は正なら group_b が高い（上昇）。",
        "note": ("n_unannotated は InChIKey が付かず書き出さなかった行数です。"
                 "「変化が無かった」ではなく「調べていない」行です。"),
    })


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True),
    structured_output=False)
def dataset_set_sample_metadata(manifest_path: str, dataset_handle: str | None = None) -> str:
    """実験情報シート(sample-manifest.v1)を読み込み、検証してから DatasetState へ
    一括適用する（spec §7.3）。

    上流のConsole実行をやり直さずに、群・バッチ・注入順・qc_pool などの誤記入を
    このツール単体で訂正できる。列の意味・拒否条件は sample-manifest.v1 の仕様
    （`docs/superpowers/specs/2026-09-05-raw-folder-pipeline-integrity-and-metadata-design.md`
    §7.2）を参照。

    検証は全件そろってから一括反映する（**一部だけ適用される状態は作らない**）。
    シートに不備があれば `session.dataset` を一切変更せずエラーを返すので、
    そのまま再送してよい。

    任意の列 `subject_id`（見出しの末尾に足す。反復測定の被験者）も読み、`sample_meta` に持つ
    （統計にはまだ使わない）。列が無いシートは今までどおり。

    前処理入力（role/batch/injection_order/qc_pool/include/sample_id/source_file）が
    変われば前処理済み行列ごと・PCA・差次的解析まで無効化する。`group` だけの変更なら
    差次的解析だけを無効化し、PCA は生かしたまま返す（`changed_fields` で確認できる）。

    比較の前提検証は 2 つの比較で同じ向きに守る: pipeline は `resolve_comparison`/`run_comparison`
    （内部関数。v1 の Welch t）、対話の `dataset_compare` は統計エンジン（既定は moderated t）が、
    群の明示を求め・`include=false` と非 sample の役割を外し・完全交絡を `allow_confounded` なしには止める。

    dataset_handle: `load_dataset` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に load_dataset を実行してください。")

    from metabolomix.analysis.dataset_service import apply_sample_manifest

    try:
        summary = apply_sample_manifest(ds, manifest_path)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    except OSError as exc:
        return mztab_error(
            "SAMPLE_MANIFEST_NOT_FOUND", str(exc), {"manifest_path": manifest_path})

    return json_payload({
        "status": "success",
        **summary,
        "next": ("changed_fields に前処理入力が含まれる場合は dataset_preprocess から"
                 "やり直してください。group のみの変更なら dataset_compare を"
                 "再実行するだけで構いません。"),
    })


# ---------- 内部ヘルパ ----------

def _missing(state: str, message: str) -> str:
    return missing_state(state, _RECOVERY_TOOLS[state], message)


# ---------- dataset_feature_search（spec 2026-10-19） ----------

_POLARITIES = ("positive", "negative")

# 丸めの桁（spec §4.2）。m/z 4 桁・ppm 2 桁・RT 3 桁・log2FC 4 桁。p・q・平均は有効 4 桁。
_CANDIDATE_DIGITS = {"theoretical_mz": 4, "mz": 4, "ppm": 2, "rt": 3, "delta_rt": 3}
_WINDOW_DIGITS = {"mz": 4, "rt": 3}
_STATS_DIGITS = {"log2fc": 4}
_STATS_SIGNIFICANT = ("p", "q", "mean_a", "mean_b")


def _round(value, digits: int):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    value = float(value)
    if not math.isfinite(value):
        return value
    rounded = round(value, digits)
    return 0.0 if rounded == 0 else rounded          # "-0.0" を出さない


def _significant(value, digits: int = 4):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    value = float(value)
    if not math.isfinite(value):
        return value
    return float(f"{value:.{digits}g}")


def _cell(value):
    """TSV の 1 セル。数値は丸め済みのものを repr のまま出す（format_spots_as_table の
    固定小数 4 桁に任せると p = 1e-6 が 0 になる）。bool は JSON と同じ小文字。"""
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _table(rows: list[dict], columns, *, digits: dict | None = None,
           significant: tuple = ()) -> str:
    """行を TSV（列名 1 回）にする。行が無ければ空文字。"""
    from metabolomix.arf2.reader import format_spots_as_table

    digits = digits or {}
    cells = []
    for row in rows:
        out = {}
        for column in columns:
            value = row.get(column)
            if column in digits:
                value = _round(value, digits[column])
            elif column in significant:
                value = _significant(value)
            out[column] = _cell(value)
        cells.append(out)
    return format_spots_as_table(cells, columns=list(columns)) if cells else ""


def _curation_note(ds) -> str | None:
    """キュレーションのフラグがあれば「反映していない」注記。読めなくても検索は止めない。

    フラグの見つけ方は dataset_export_differential と同じ（兄弟の `.arf2`）。
    """
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError

    lookup = curation_apply.arf2_lookup_kwargs(ds)
    arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p, **lookup)
                              for p in (getattr(ds, "source_files", None) or {}))), None)
    if arf2 is None:
        return None
    try:
        n_flags = curation_apply.flags_for_arf2(arf2)["n"]
    except (FlagFileError, OSError) as exc:
        return (f"キュレーションのフラグ記録を読めなかった（{exc}）。この検索はキュレーションを"
                "反映していない（SME / SML の生の同定だけを見ている）。")
    if not n_flags:
        return None
    return (f"キュレーションのフラグが {n_flags} 件あるが、この検索は反映していない"
            "（SME / SML の生の同定だけを見ている）。dataset_export_differential の同定とは"
            "食い違いうる。")


def _target_summary(match: dict) -> dict:
    """目標 1 件の判定（候補の行は TSV の側に 1 回だけ載せるので外す）。"""
    keys = ("target_id", "status", "reason", "selected_feature_id", "selected_adduct",
            "selected_evidence", "adduct_tier_used", "counts", "truncated")
    out = {key: match[key] for key in keys}
    if "merged_from" in match:
        out["merged_from"] = match["merged_from"]           # merge_split=true のときだけ
    out["theoretical_mz"] = {adduct: _round(mz, 4)
                             for adduct, mz in match["theoretical_mz"].items()}
    return out


def _rule_targets(normalized: list[dict]) -> list[dict]:
    """戻り値に写す正規化後の目標（理論 m/z だけ 4 桁に丸める。hash は丸める前から作る）。"""
    shown = []
    for target in normalized:
        item = dict(target)
        item["theoretical_mz"] = {adduct: _round(mz, 4)
                                  for adduct, mz in target["theoretical_mz"].items()}
        if "ppm_offset" in target:
            item["ppm_offset"] = {adduct: _round(value, 2)
                                  for adduct, value in target["ppm_offset"].items()}
        shown.append(item)
    return shown


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_feature_search(
    feature_ids: list[str] | None = None,
    min_mz: float | None = None,
    max_mz: float | None = None,
    min_rt: float | None = None,
    max_rt: float | None = None,
    targets: list[dict] | None = None,
    result_ids: list[str] | None = None,
    max_results: int = 20,
    with_group_labels: bool = False,
    merge_split: bool = True,
    alignment_qc_path: str | None = None,
    ppm_offset: dict | None = None,
    sort: str = "mz",
    dataset_handle: str | None = None,
    contrast_id: str | None = None,
) -> str:
    """DatasetState の特徴を探し、目的の物質に当たる特徴を規則で判定する（読み取り専用）。

    2 つの使い方があり、`targets` の有無で切り替わる（同時指定は `DATASET_BAD_REQUEST`）。

    **窓検索**（`targets` なし）: `feature_ids`（SMF_ID の完全一致）・`min_mz`/`max_mz`・
    `min_rt`/`max_rt`（分）の AND。境界を含む。何も指定しなければ `FEATURE_SEARCH_EMPTY_QUERY`
    （全件を返さない。`sort="intensity"` なら条件なしで上位 `max_results` 件）。`feature_ids=[]` は
    「どの ID でもない」条件で 0 件。`dataset_exclude` した特徴は出さない。
    並びは `sort="mz"`（既定）で m/z 昇順、`sort="intensity"` で強度の降順。強度は検出セル
    （gap-fill でない実測。除外した試料の列を除く）の値の**最大**で、ピークの高さ（データセットの
    測定値）の最大——EIC のクロマトグラムの最大ではない。検出セルが無い特徴は最後。検出の印が無いデータセットは全セルの最大で並べ、`notes` に書く。
    戻り値 `features` は TSV（`feature_id mz rt evidence name inchikey adduct n_sme ontology
    n_detected normalized msi_level`。`evidence` は `sme`/`sml`/`none`、`n_detected` は検出セルの数
    （除外した試料を除く。検出状態が無ければ空）、`normalized`（GOSLIN で標準化した名前）と
    `msi_level` は解析フォルダのデータセットだけが持ち、mzTab-M では空。意味は
    `lipidmix://docs/output-format/identity`）。注釈の無い特徴にも統計を付ける。

    **目標照合**（`targets` あり、1 回 50 件まで）: 各要素は `target_id`・`inchikey`（必須。
    先頭 14 文字で照合）と、`formula` か（`expected_mz` と `adduct`）。InChIKey が一意に決まらない
    脂質の種は `inchikey` の代わりに `lipid_species`（`PC 34:4`・`PC 18:3_16:1`。GOSLIN の種の
    レベルに直し、SME 候補の名前の種と照らす。分子種は確かめない。opt-in）。任意で `name` `ppm`
    （既定 10）`rt_min` `rt_tol_min`（既定 0.3）`adduct_tiers` `polarity`。付加体の層を上から、
    証拠の水準 `sme_match` → `sml_match` → `unannotated` の順に見て、qualified がちょうど 1 つの
    ときだけ `resolved`。2 つ以上なら `needs_input`（`multiple_qualified_candidates`）、無ければ
    `needs_input`（`no_qualified_candidate`）。**強度・ppm・RT の近さでは選ばない。**
    `unannotated` は `rt_min` があるときだけ qualified。SME（MS/MS 同定）が別物質なら落とす。
    極性はジョブの宣言 → mzTab-M の付加体の多数決の順に決め、`rule.polarity_source` に書く。
    目標の極性は `polarity`、省けば明示した付加体（`adduct`/`adduct_tiers`）の符号で、
    データセットと違う目標は `not_applicable`（陽陰が混ざる付加体は `FEATURE_TARGET_INVALID`）。
    戻り値 `targets`（目標ごとの判定と段ごとの件数）・`candidates`（TSV。目標ごとに
    `max_results` 件まで）・`rule`（正規化した規則と `rule_hash`。result_id・表示名 `name`・
    目標の並び順を含まないので、同じ目標・同じ極性なら別の session でも同じ hash）。hash は表示用の丸め（理論 m/z を
    4 桁に）より前の値から作るので、返された `rule.targets` から計算し直すことはできない
    ——照合は `rule_hash` の値どうしで行う。

    result_ids: 統計を付ける 2 群比較（`dataset_compare` の結果の対比、v1 の 2 群の結果）。省略時は
        今も有効な 2 群比較の全部（古いものは使わず `notes` に挙げる）。`dataset_compare` の結果は対比ごとに
        1 つの比較として数え、ID は対比が 1 つなら `result_id`、複数なら `result_id#contrast_id`
        （交互作用の対比は 2 群の形にならないので使わない）。名指しした ID が無い・古いときは
        `missing_state`。差次的結果が 1 件も無ければ統計なしで返す。
    contrast_id: `dataset_compare` の結果から使う対比を 1 つに絞る（その対比が無い結果は飛ばす）。
    max_results: 窓検索の行数・目標ごとの候補数の上限（既定 20）。打ち切りは `notes` に書く。
    sort: 窓検索の並び（`mz` 既定 / `intensity`）。`targets` とは併用できない（`DATASET_BAD_REQUEST`。
        目標照合は規則の順で並べる）。

    `stats` は TSV（`feature_id result_id status log2fc p q mean_a mean_b detected_a detected_b`）。
    `status` は `tested` / `not_in_matrix`（前処理で落ちた）/ `not_testable`（p が無い）で、
    後の 2 つは数値欄が空（「差が無かった」ではない）。`detected_*` は群ごとの実検出例数
    （gap-fill を数えない。検出状態が無ければ空）。`results` に BH の母集団
    （`bh_population`）と caveat。キュレーションのフラグは反映しない（あれば `notes`）。
    欄の意味は `lipidmix://docs/output-format/mztab`。

    with_group_labels: True（opt-in。再解析で使う）なら `stats` の末尾に `comparison`
        （pipeline の比較なら comparison_id、そうでなければ `<group_b>_vs_<group_a>`）・
        `group_a`・`group_b`・`log2fc_test_over_ref`（`log2fc` と同じ値を、比較名と同じ向きの
        名前で）を足す。既定 False では列は変わらない。

    merge_split: True（既定）なら、1 つの物質が 2 つに割れたアライメントスポット
        （`dataset_alignment_qc` の `merge_groups`。互いに最も近い mergeable な 2 スポット）を
        試料ごとの大きい方の値で 1 列にまとめた行列で、差次的結果と同じ前処理のレシピ・同じ検定
        （log 変換の Welch 等）・**まとめた後の全特徴**での BH をやり直した値を `stats` に出す。
        まとめる組がデータセットに 1 つでもあれば**すべての行**をやり直した値にする（1 つの表の q は
        1 つの BH の母集団から）。組が無ければやり直さず従来の値。目標照合でも組の 2 スポットを
        1 つの候補（残す ID。`targets[].merged_from`・`candidates` の末尾の `merged_from`）と数える。
        `stats` の末尾に `merged_from`（まとめた 2 つの ID を `;` で。先頭が残す ID＝受け持ちの多い側、同数なら小さい ID で、消える
        側の ID を頼んでも残す ID の行に解く）、`stats_unmerged` にまとめる前の従来の `stats`
        （同じ列・頼んだ ID のまま）、`merge_split` に出所・まとめた組・やり直した結果と BH の
        母集団。False なら従来どおり（`stats_unmerged` も `merged_from` 列も無い）。
    alignment_qc_path: `dataset_alignment_qc` が書いた JSON（今のデータセットで作ったものだけ）。
        省略時はその場で割れを探す（既定の閾値）。`merge_split=False` とは併用できない。

    ppm_offset: 目標照合の質量の窓の中心のずれ（任意）。`{"default": 0.0, "bands": [{"max_mz": 200,
        "offset": -8.3}]}`。中心は「実測 − 理論」の ppm の符号つきの値で、帯は理論 m/z が `max_mz`
        未満の最初の帯（無ければ `default`）。幅（目標の `ppm`）は変えず、候補の `ppm` 列も実測 − 理論の
        まま。`rule.ppm_offset` と各目標の `ppm_offset`（付加体ごと）に出て `rule_hash` に入る。
        窓検索（`targets` なし）では使えない（`DATASET_BAD_REQUEST`）。

    dataset_handle: `load_dataset` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に load_dataset を実行してください。")
    if ppm_offset is not None and targets is None:
        return mztab_error("DATASET_BAD_REQUEST",
                           "ppm_offset は目標照合（targets）のときだけ使える（効かない設定を渡さない）。",
                           {"field": "ppm_offset"})
    if not isinstance(merge_split, bool):
        return mztab_error("DATASET_BAD_REQUEST", f"merge_split は true / false: {merge_split!r}",
                           {"field": "merge_split"})
    if alignment_qc_path is not None and not merge_split:
        return mztab_error("DATASET_BAD_REQUEST",
                           "alignment_qc_path は merge_split=true のときだけ使える"
                           "（効かない設定を渡さない）。", {"field": "alignment_qc_path"})

    if sort not in ("mz", "intensity"):
        return mztab_error("DATASET_BAD_REQUEST", f"sort は mz か intensity: {sort!r}",
                           {"field": "sort"})
    if targets is not None and sort != "mz":
        return mztab_error("DATASET_BAD_REQUEST",
                           "sort は窓検索（targets なし）のときだけ使える（目標照合は規則の順で"
                           "並ぶ。効かない設定を渡さない）。", {"field": "sort"})

    window = {"feature_ids": feature_ids, "min_mz": min_mz, "max_mz": max_mz,
              "min_rt": min_rt, "max_rt": max_rt}
    given_window = [name for name, value in window.items() if value is not None]
    if targets is not None and given_window:
        return mztab_error(
            "DATASET_BAD_REQUEST",
            "targets（目標照合）と窓の引数（feature_ids・m/z・RT）は同時に指定できません"
            "（どちらの意味か曖昧になる）。どちらか一方で呼び直してください。",
            {"window_arguments": given_window})
    if isinstance(max_results, bool) or not isinstance(max_results, int) or max_results < 1:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"max_results は 1 以上の整数で指定してください: {max_results!r}",
                           {"field": "max_results"})

    from metabolomix.analysis import feature_search
    from metabolomix.analysis.mass_accuracy import normalize_ppm_offset
    from metabolomix.analysis.result_state import dataset_fingerprint

    notes: list[str] = []
    payload: dict = {}
    # データの指紋（行列全体の hash）は 1 回だけ計算し、割れの照合・引き継いだ結果の照合・
    # stats_source のすべてで使う。
    dataset_hash = dataset_fingerprint(ds)
    plan = None
    try:
        if targets is None:
            rows, truncated = feature_search.search_window(
                ds, **window, max_results=max_results, sort=sort)
            if feature_ids is not None and not feature_ids:
                notes.append("feature_ids=[] は「どの ID でもない」条件として扱った（0 件）。")
            for fid in (feature_ids or []):
                if fid in ds.excluded_features:
                    notes.append(f"{fid} は除外済み（理由: {ds.excluded_features[fid]}）なので返さない"
                                 "（戻すなら dataset_exclude(mode=\"remove\")）。")
            if sort == "intensity" and ds.detected_mask is None:
                notes.append("このデータセットには検出の印（detected_mask）が無く gap-fill のセルを"
                             "見分けられないので、強度順は全セル（除外した試料を除く）の最大で並べた。")
            if truncated:
                notes.append(f"max_results={max_results} で打ち切った（ほか {truncated} 件）。"
                             "窓を狭めるか max_results を増やしてください。")
            payload["features"] = _table(rows, feature_search.WINDOW_COLUMNS,
                                         digits=_WINDOW_DIGITS)
            stat_features = [row["feature_id"] for row in rows]
        else:
            if isinstance(targets, list) and not targets:
                return mztab_error("FEATURE_TARGET_INVALID", "targets が空です（1 件以上渡す）。",
                                   {"target_id": None, "field": "targets"})
            polarity, polarity_source, polarity_notes = _dataset_polarity(ds)
            notes.extend(polarity_notes)
            offset_spec = normalize_ppm_offset(ppm_offset, code="DATASET_BAD_REQUEST")
            normalized = feature_search.normalize_targets(targets, dataset_polarity=polarity,
                                                          ppm_offset=offset_spec)
            # 特徴の m/z 配列は 1 回だけ作り、全目標で使い回す
            mz_array = feature_search._feature_mz_array(ds)
            merged_map = None
            if merge_split:
                # 割れたスポットの組は判定でも 1 つの候補（残す ID）として数える
                plan = feature_search.merge_plan(ds, alignment_qc_path=alignment_qc_path,
                                                 dataset_hash=dataset_hash)
                merged_map = feature_search.merge_map(plan["groups"])
            matches = [feature_search.match_target(ds, target, max_results=max_results,
                                                   mz_array=mz_array, merged_from=merged_map)
                       for target in normalized]
            cut = [f"{m['target_id']}: {m['truncated']} 件" for m in matches if m["truncated"]]
            if cut:
                notes.append(f"候補を目標ごとに max_results={max_results} 件で打ち切った"
                             f"（{', '.join(cut)}）。")
            candidates = [row for m in matches for row in m["candidates"]]
            payload["targets"] = [_target_summary(m) for m in matches]
            if merge_split:
                shown = [{**row, "merged_from": ";".join(row["merged_from"] or []) or None}
                         for row in candidates]
                payload["candidates"] = _table(
                    shown, feature_search.CANDIDATE_COLUMNS + ("merged_from",),
                    digits=_CANDIDATE_DIGITS)
            else:
                payload["candidates"] = _table(candidates, feature_search.CANDIDATE_COLUMNS,
                                               digits=_CANDIDATE_DIGITS)
            stat_features = list(dict.fromkeys(row["feature_id"] for row in candidates))
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    try:
        chosen = feature_search.select_results(
            ds, None if result_ids is None else list(dict.fromkeys(result_ids)),
            dataset_hash=dataset_hash, contrast_id=contrast_id)
    except DomainError as exc:                    # 対比の指定の誤り（CONTRAST_* / DATASET_BAD_REQUEST）
        return mztab_error(exc.code, exc.message, exc.details or None)
    if result_ids is not None:
        if chosen["missing"] or chosen["stale"]:
            # 名指しの結果を黙って落とすと、呼び出し側は統計が付いたつもりで読む
            return _missing(
                "dataset_compare_result",
                "result_ids に指定した差次的結果が使えません"
                f"（無い・2 群比較でない: {chosen['missing']} / 古い: {chosen['stale']}）。"
                "dataset_compare で計算し直すか、result_ids を省略して今も有効な結果だけを"
                "使ってください。")
    else:
        if chosen["stale"]:
            notes.append("古い差次的結果（前処理・群の訂正前）は統計に使っていない: "
                         + ", ".join(chosen["stale"]))
    results = chosen["results"]
    if not results:
        notes.append("差次的結果（2 群比較）が 1 件も無いので統計を付けていない"
                     "（dataset_compare を実行すると付く）。")

    merge_block: dict = {"enabled": merge_split}
    unmerged = None
    if not merge_split:
        stats = feature_search.feature_stats(ds, stat_features, results,
                                             with_labels=with_group_labels)
    else:
        try:
            if plan is None and results and stat_features:
                # 窓検索で統計を付けないなら割れを探さない（行列全体を掃くので）
                plan = feature_search.merge_plan(ds, alignment_qc_path=alignment_qc_path,
                                                 dataset_hash=dataset_hash)
            plan = plan or {"source": None, "groups": [], "alignment_qc_sha256": None}
            merged = feature_search.merged_feature_stats(ds, stat_features, results,
                                                         groups=plan["groups"],
                                                         with_labels=with_group_labels)
        except DomainError as exc:
            return mztab_error(exc.code, exc.message, exc.details or None)
        stats, unmerged = merged["stats"], merged["stats_unmerged"]
        for row in stats:
            row["merged_from"] = ";".join(row["merged_from"]) if row["merged_from"] else None
        merge_block.update({
            "source": plan["source"],
            "alignment_qc_sha256": plan["alignment_qc_sha256"],
            "n_groups": len(plan["groups"]),
            "groups": list(merged["touched"].values()),
            "renamed": merged["renamed"],
            "recomputed": list(merged["results"]),
            "bh_population": {
                summary["result_id"]: summary["bh_population"]
                for summary in feature_search.result_summaries(merged["view"],
                                                               list(merged["results"].values()))},
        })
        if merged["renamed"]:
            notes.append("割れたスポットをまとめ、消える側の ID を残す ID の行に出した: "
                         + ", ".join(f"{a} → {b}" for a, b in merged["renamed"].items()))
    if any(row["detected_a"] is None for row in stats):
        if getattr(ds, "detected_mask", None) is None:
            notes.append("検出状態（gap-fill の区別）を取り込めていないので detected_a / "
                         "detected_b は空（非ゼロ値は実測と補間を区別しない）。")
        else:
            notes.append("結果の試料と検出状態の列の対応を保証できず、detected_a / detected_b を"
                         "数えていない結果がある。")
    columns = (feature_search.STATS_COLUMNS_LABELED if with_group_labels
               else feature_search.STATS_COLUMNS)
    digits = {**_STATS_DIGITS, "log2fc_test_over_ref": 4}
    if merge_split:
        payload["stats"] = _table(stats, columns + ("merged_from",), digits=digits,
                                  significant=_STATS_SIGNIFICANT)
        payload["stats_unmerged"] = _table(unmerged, columns, digits=digits,
                                           significant=_STATS_SIGNIFICANT)
    else:
        payload["stats"] = _table(stats, columns, digits=digits, significant=_STATS_SIGNIFICANT)
    payload["merge_split"] = merge_block
    payload["results"] = feature_search.result_summaries(ds, results)
    if targets is not None:
        from metabolomix.analysis.analysis_spec import merge_groups_hash

        # merge_split=true は判定でまとめた割れの組の hash を規則の隣に置き、rule_hash にも入れる
        # （false は v7 と同じ規則・同じ hash）
        groups_hash = merge_groups_hash(plan["groups"]) if merge_split else None
        payload["rule"] = {"targets": _rule_targets(normalized), "polarity": polarity,
                           "polarity_source": polarity_source,
                           "merge_split": merge_split, "merge_groups_hash": groups_hash,
                           "rule_hash": feature_search.rule_hash(
                               normalized, merge_groups_hash=groups_hash)}
        if offset_spec is not None:
            # 渡したときだけ出す（省けば従来と同じ戻り値）
            payload["rule"]["ppm_offset"] = offset_spec
    payload["stats_source"] = {
        "result_ids": [(r.get("provenance") or {}).get("result_id") for r in results],
        "dataset_hash": dataset_hash,
    }
    curation = _curation_note(ds)
    if curation:
        notes.append(curation)
    payload["notes"] = notes
    return json_payload(payload)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
def dataset_export_evidence_bundle(output_dir: str, result_id: str, apply_curation: bool = True,
                                   dataset_handle: str | None = None,
                                   contrast_id: str | None = None) -> str:
    """指定した 2 群結果を根拠バンドル v3（samples・intensities・adduct・comparison.polarity）へ出力する。

    `dataset_compare` の結果（`result_id`）は `contrast_id` で対比を 1 つ選ぶ（対比が 1 つなら省略可。
    複数で省くか未知の ID なら `CONTRAST_REQUIRED`、交互作用は `CONTRAST_NOT_TWO_GROUP`）。
    `comparison.method` は実際の検定名（`moderated_t` など。v1 の結果は `Welch`）、`comparison.statistics`
    に検定・CI 水準・判定規則・BH の母集団・単位。`intensities.tsv` は結果が読んだ行列の線形値
    （正規化の後・log2 の前、補完なし＝`pre_imputation`）。

    spectrum_links は同じアラインメントの `.dcl` へ ID（SMF_ID == AlignmentID == `.dcl` の番号）で
    結び、m/z・RT は検証にだけ使う。`.arf2` / `.dcl` はハンドルならそのハンドルのもの、対話の
    セッションなら `load_dataset` と同じ探し方で決める。結べないときは links を空にし、
    manifest の `provenance.spectrum_links_status` と `spectrum_links_reason` に理由を書く。

    dataset_handle: `load_dataset` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
    """
    from metabolomix.analysis.evidence_export import export_dataset_evidence
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    try:
        result_id, contrast_id = stats_results.split_result_ref(result_id, contrast_id)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    result = (getattr(ds,'results',None) or {}).get(result_id) if ds is not None else None
    if result is None:
        return _missing('dataset_compare_result','指定した2群結果がありません。先にdataset_compareを実行してください。')
    try:
        curation = None
        lookup = curation_apply.arf2_lookup_kwargs(ds)
        arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p, **lookup) for p in ds.source_files)), None)
        if arf2 is not None:
            flags = curation_apply.flags_for_arf2(arf2)
            mapped = (getattr(ds,'feature_qc',None) or {}).get('source') == 'arf'
            state = 'not_applied' if not apply_curation else 'applied' if mapped else 'unmapped'
            curation = dict(state=state,flag_set=flags)
        polarity = holder.polarity if dataset_handle else _dataset_polarity(ds)[0]
        if dataset_handle:
            arf2_path, dcl_path = holder.arf2_path, holder.dcl_path
        else:
            from metabolomix.tools.handle_tools import _alignment_files
            arf2_path, dcl_path = _alignment_files(ds, None)
        return json_payload(export_dataset_evidence(ds,result,Path(output_dir),curation=curation,
                                                    polarity=polarity,arf2_path=arf2_path,
                                                    dcl_path=dcl_path,contrast_id=contrast_id))
    except DomainError as exc:
        if exc.code == 'STATS_MATRIX_MISSING':
            return _missing('dataset_compare_result',exc.message+' dataset_compareを実行し直してください。')
        return mztab_error(exc.code,exc.message,exc.details or None)
    except (ValueError,OSError,FlagFileError) as exc:
        return mztab_error('EVIDENCE_EXPORT_FAILED',str(exc))
