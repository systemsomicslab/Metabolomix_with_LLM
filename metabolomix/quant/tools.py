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


def _parent_id(dataset_handle) -> str | None:
    """新しいハンドルの親＝この呼び出しが実際に読んだハンドル。

    `dataset_handle` を省いた呼び出しは `handle_dispatch` が直近のハンドルに固定して回す
    （`session_state.pinned_handle`）。固定の間に別の呼び出しが直近を切り替えうるので、
    `current_handle_id` ではなく固定先を親にする。固定が無いとき（直近のハンドルが無い・
    同期で呼ばれた）だけ直近を読む。
    """
    if dataset_handle:
        return dataset_handle
    return session_state.pinned_handle_id() or getattr(session_state.session,
                                                       "current_handle_id", None)


def _context(dataset_handle, polarity):
    """(holder, ds, polarity, 親ハンドル ID) か、返すべき封筒（str）。"""
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = getattr(holder, "dataset", None)
    if ds is None:
        return missing_state(
            "dataset", ["load_dataset"],
            "データセットがありません。load_dataset を実行するか、開いている"
            "ハンドルを渡してください。")
    try:
        # 統合したデータセットは極性が決まらないので、極性を尋ねる前に統合済みと答える。
        normalize.check_not_merged(ds)
        resolved = resolve_quant_polarity(ds, holder, polarity)
    except DomainError as exc:
        return quant_error(exc.code, exc.message, exc.details or None)
    return holder, ds, resolved, _parent_id(dataset_handle)


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

    reference_path: 内部標準の参照表 CSV（毎回指定。設定の既定は無い）。戻り値は表の sha256・単位・ファイル名の
    末尾だけで、パスと表の中身は返さない。classes はデータセットに現れたクラスだけ。内部標準の候補が複数・無い
    ものは unresolved に候補を並べ、status は needs_input。
    bindings: {内部標準の名前: 特徴 ID}（unresolved の standard_name が鍵）で候補を選ぶ。
    濃度から係数を作るなら concentrations（{内部標準名: 濃度}）と solvent_volume・cells_millions（unit は任意）。
    polarity: 極性が決まらないデータセットでは必須。
    MS-DIAL の解析フォルダ（.arf / .arf2）から読んだデータセットだけ（ほかは QUANT_SOURCE_UNSUPPORTED、陽陰を統合
    したものは QUANT_ALREADY_MERGED）。規則プラグイン（[quant] rules）が要る。
    dataset_handle: 省けば直近のデータセット。
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

    引数は dataset_quant_preview と同じ。内部標準が未解決（needs_input）なら何もせず preview と同じ封筒を返す。
    成功なら解析フォルダの reports/quant/ に記録を書き、新しい dataset_handle・status_counts・単位・record_id を返す
    （gap-fill のセルも値を計算し、detected_mask で印を引き継ぐ。内部標準が 0・負・欠測の試料は欠測にして理由を
    数える）。新しいハンドルの極性は解いた polarity。親が探索専用なら派生も探索専用（戻り値の exploratory）。
    skip: True は値を変えず「正規化を飛ばした」記録と派生を作る（reference_path は skip=True のときだけ省け、読まない）。
    掛けないデータセット: 正規化済み（QUANT_ALREADY_NORMALIZED）・選別済み（QUANT_ALREADY_SELECTED）・陽陰を統合
    したもの（QUANT_ALREADY_MERGED。内部標準は極性ごとなので、正規化は統合の前に極性ごとに行う）。統合し直すと、
    統合したハンドルに当てていた除外とシートは引き継がれない（極性ごとのハンドルの除外は引き継ぐ）ので、統合し直して
    できたハンドルに当て直す。
    dataset_handle: 省けば直近のデータセット。
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

    positive_handle / negative_handle: load_dataset などで開いたハンドル。正規化の状態（normalized・method・unit）が
        違えば統合しない（QUANT_MERGE_STATE_MISMATCH）。dataset_quant_apply を通していない側は未正規化として扱う
        （skip=True を挟まなくてよい）。
    mapping_csv: 試料の対応（列: mode, source_sample, sample_id, group, sample_type。mode は positive / negative、
        pos / neg も可。group / sample_type は規則が検証するだけで適用しない）。試料は名前（極性の印を除いた名前）で
        照合し、食い違えば needs_input（only_in に片側にしか無い試料）で止まるので、これで対応づけて呼び直す。
    record_folder: 記録の置き場（既存のフォルダ）。既定は 2 つの解析フォルダの共通の親の reports/quant/ で、
        決められなければ needs_input で求める。

    特徴 ID は positive:<id> / negative:<id>（行は POS の特徴のあとに NEG）、試料の列は共通 ID（POS の順）。試料の
    表は陽イオンのもの（陽イオンに当てたシートも引き継ぐ）で、陰イオンは試料ごとのバッチと注入順を残す
    （negative_runs に陽イオンとの違い）。極性ごとの注入順・バッチは統合の前にそれぞれのハンドルへ
    dataset_set_sample_metadata で入れる（統合後のハンドルへのシートは陽イオンの測定だけを決める）。
    極性ごとのハンドルの除外（dataset_exclude）は統合後の ID に写して引き継ぐ。前の統合したハンドルに当てたシートと
    除外は、統合し直すと引き継がれない。キュレーションのフラグは写さず、選別が極性ごとの .arf2 からその時点のものを
    読む。どちらかが探索専用なら統合も探索専用。assay_kind は両側から決め、lipid と metabolite で食い違えば unknown
    （直すなら誤っている側の極性のハンドルへ record_objective で記録し直して統合し直す）。規則プラグイン
    （[quant] rules）が要る。2 つのハンドルのロックを待つのは合わせて最長 25 秒（待ちきれなければ何もせず
    HANDLE_BUSY）。戻り値の欄（normalize・negative_runs・exclusions・assay_kind など）と needs_input の答え方は
    `docs_read(topic="quant")` の §4・§5。
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

    順序: キュレーション（apply_curation=True のとき。選別の時点の flags.jsonl を読み直し、wrong / redundant は除外、
    assign は名前・クラスの差し替え）→ 内部標準（正規化した記録の内部標準。正規化していなければ除かない）→ 名前・
    コメント・MS/MS の除外 → 付加体の採用表（lipidmix.local.toml の [quant] adducts）。採用表が未設定なら
    needs_input（設定するか skip_adduct_filter=True で付加体の選別だけ飛ばす）。統合前のデータセットにも使える
    （mode はそのデータセットの極性。決まらなければ polarity が必須。統合後は各特徴の出自の極性）。
    全判断は解析フォルダの reports/quant/<record_id>.audit.json に書き、戻り値は理由ごとの件数と例 5 件・record_id・
    採用表の sha256 とファイル名だけ。規則プラグイン（[quant] rules）が要る。親が探索専用なら派生も探索専用
    （戻り値の exploratory）。
    dataset_handle: 省けば直近のデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    return select_mod.select_tool(holder, parent_handle_id=_parent_id(dataset_handle),
                                  skip_adduct_filter=skip_adduct_filter,
                                  apply_curation=apply_curation, polarity=polarity)
