"""定量の MCP ツール（Task 11: 内部標準の確認と正規化）。

実体は `metabolomix.quant.normalize`。ここは引数を解いて呼ぶ薄い層。server.py への
取り込みと登録数の更新は Task 14。
"""
from __future__ import annotations

from mcp.types import ToolAnnotations

from metabolomix.core import session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state
from metabolomix.quant import merge, normalize
from metabolomix.quant import select as select_mod
from metabolomix.quant.inputs import quant_error, reference_input
from metabolomix.quant.table import resolve_quant_polarity
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_quant_preview", "dataset_quant_apply", "dataset_merge_polarities",
           "dataset_select_features"]


def _context(dataset_handle, polarity):
    """(holder, ds, polarity, 親ハンドル ID) か、返すべき封筒（str）。"""
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = getattr(holder, "dataset", None)
    if ds is None:
        return missing_state(
            "dataset", ["dataset_load", "dataset_handle_open"],
            "データセットがありません。dataset_load を実行するか、dataset_handle_open で開いた"
            "ハンドルを渡してください。")
    try:
        resolved = resolve_quant_polarity(ds, holder, polarity)
    except DomainError as exc:
        return quant_error(exc.code, exc.message, exc.details or None)
    parent = dataset_handle or getattr(session_state.session, "current_handle_id", None)
    return holder, ds, resolved, parent


def _reference(reference_path, method, concentrations, solvent_volume, cells_millions, unit):
    try:
        return reference_input(reference_path, method=method, concentrations=concentrations,
                               solvent_volume=solvent_volume, cells_millions=cells_millions,
                               unit=unit)
    except DomainError as exc:
        return quant_error(exc.code, exc.message, exc.details or None)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_quant_preview(reference_path: str, bindings: dict | None = None,
                          method: str = "cell_coefficient", concentrations: dict | None = None,
                          solvent_volume: float | None = None,
                          cells_millions: float | None = None, unit: str | None = None,
                          polarity: str | None = None,
                          dataset_handle: str | None = None) -> str:
    """内部標準の参照表で、どのクラスをどの内部標準・係数で割るかを確かめる（値は作らない）。

    reference_path: 内部標準の参照表 CSV（毎回指定。設定の既定は無い）。戻り値は表の sha256・
    単位・ファイル名の末尾だけで、パスと表の中身は返さない。classes はデータセットに現れた
    クラスだけ。内部標準の候補が複数・無いものは unresolved に候補を並べ status は needs_input。
    bindings は {対象クラス: 内部標準の特徴 ID} で候補を選ぶ。濃度から係数を作るなら
    concentrations（{内部標準名: 濃度}）と solvent_volume・cells_millions（unit 任意）。
    polarity は極性が決まらないデータセットのとき必須。規則プラグイン（[quant] rules）が要る。
    dataset_handle: ハンドルのデータセット。省略時は直近のデータセット。
    """
    ctx = _context(dataset_handle, polarity)
    if isinstance(ctx, str):
        return ctx
    _holder, ds, resolved, _parent = ctx
    ref = _reference(reference_path, method, concentrations, solvent_volume, cells_millions, unit)
    if isinstance(ref, str):
        return ref
    return normalize.preview_tool(ds, ref, bindings, resolved)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_quant_apply(reference_path: str | None = None, bindings: dict | None = None,
                        method: str = "cell_coefficient", concentrations: dict | None = None,
                        solvent_volume: float | None = None,
                        cells_millions: float | None = None, unit: str | None = None,
                        polarity: str | None = None, skip: bool = False,
                        dataset_handle: str | None = None) -> str:
    """内部標準で正規化し、新しいハンドルの派生データセットを作る（元は変えない）。

    引数は dataset_quant_preview と同じ。内部標準が未解決（needs_input）なら何もせず
    preview と同じ封筒を返す。成功なら解析フォルダの reports/quant/ に記録を書き、新しい
    dataset_handle・status_counts・単位・record_id を返す（gap-fill のセルも値を計算し、
    detected_mask で印を引き継ぐ。内部標準が 0・負・欠測の試料は欠測にして理由を数える）。
    正規化済みのデータセットには掛けない。skip=True は値を変えず「正規化を飛ばした」記録と
    派生を作る（reference_path は skip=True のときだけ省ける。読まない）。
    dataset_handle: ハンドルのデータセット。省略時は直近のデータセット。
    """
    ctx = _context(dataset_handle, polarity)
    if isinstance(ctx, str):
        return ctx
    holder, ds, resolved, parent = ctx
    ref = None
    if not skip:
        ref = _reference(reference_path, method, concentrations, solvent_volume,
                         cells_millions, unit)
        if isinstance(ref, str):
            return ref
    return normalize.apply_tool(ds, ref, bindings, resolved, holder=holder,
                                parent_handle_id=parent, skip_step=skip)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_merge_polarities(positive_handle: str, negative_handle: str,
                             mapping_csv: str | None = None,
                             record_folder: str | None = None) -> str:
    """陽イオンと陰イオンのデータセットを 1 つに統合し、新しいハンドルを作る（元は変えない）。

    positive_handle / negative_handle: dataset_handle_open などで開いたハンドル。どちらも
    dataset_quant_apply を通した（skip=True で飛ばした記録でもよい）もので、normalized と
    単位が同じでなければ統合しない。特徴 ID は positive:<id> / negative:<id>、行は POS の特徴
    のあとに NEG の特徴、試料の列は共通 ID（POS の順）。試料は名前（極性の印を除いた名前）で
    照合し、食い違えば needs_input（only_in に片側にしか無い試料）で止まる。mapping_csv
    （列: mode, source_sample, sample_id, group, sample_type）で対応づけて再実行する。
    normalized・method・unit のどれかが違っても統合しない。mapping_csv の mode 列は
    positive / negative（pos / neg も可）。group / sample_type 列は規則が検証するだけで
    適用しない（試料の実験情報は POS のものを使う）。
    記録は 2 つの解析フォルダの共通の親の reports/quant/。共通の親が決められなければ
    needs_input で record_folder（既存のフォルダ）を求める。規則プラグイン（[quant] rules）が要る。
    """
    pos = resolve_holder(positive_handle)
    if isinstance(pos, str):
        return pos
    neg = resolve_holder(negative_handle)
    if isinstance(neg, str):
        return neg
    return merge.merge_tool(pos, neg, mapping_csv, record_folder)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_select_features(skip_adduct_filter: bool = False, apply_curation: bool = True,
                            polarity: str | None = None,
                            dataset_handle: str | None = None) -> str:
    """付加体の採用表と除外規則で特徴を選別し、残る行だけの新しいハンドルを作る（元は変えない）。

    順序: キュレーション（apply_curation=True のとき。wrong / redundant は除外、assign は名前・
    クラスの差し替え）→ 内部標準（正規化した記録の内部標準。正規化していなければ内部標準は
    除かない）→ 名前・コメント・MS/MS の除外 → 付加体の採用表（lipidmix.local.toml の
    [quant] adducts）。採用表が未設定なら needs_input（設定するか skip_adduct_filter=True で
    付加体の選別だけ飛ばす）。統合前のデータセットにも使える（mode はそのデータセットの極性。
    決まらなければ polarity が必須。統合後は各特徴の出自の極性）。
    全判断は解析フォルダの reports/quant/<record_id>.audit.json に書き、戻り値は理由ごとの
    件数と例 5 件・record_id・採用表の sha256 とファイル名だけ。規則プラグイン（[quant] rules）が要る。
    dataset_handle: ハンドルのデータセット。省略時は直近のデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    parent = dataset_handle or getattr(session_state.session, "current_handle_id", None)
    return select_mod.select_tool(holder, parent_handle_id=parent,
                                  skip_adduct_filter=skip_adduct_filter,
                                  apply_curation=apply_curation, polarity=polarity)
