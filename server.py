"""ms-data-parser MCP サーバのファサード。

実体は metabolomix/ 配下へ機能別に分割されている（依存の浅い順）:
    metabolomix.core     … FastMCP インスタンス・設定・セッション状態・パス解決・共通ヘルパ
    metabolomix.msdial   … MS-DIAL 固有のサイドカー解析（tags.xml / .mddata）・同定・検証
    metabolomix.analysis … 入力形式に依存しない数値処理（前処理・差次的解析）
    metabolomix.plots    … レンダラ中立の payload 組み立て（volcano / eic）
    metabolomix.arf / arf2 / pai2 / dcl / eic
                      … 形式ごとのバイナリパーサ（reader）と MCP ツール（tools）
    metabolomix.corpus   … 蓄積ノートの純ロジック（knowledge_store / paper_ingest）
    metabolomix.tools    … 形式に紐づかない MCP 公開層（入口・サンプル検索・目的・レポート・リソース）

呼び出し連鎖の詳細は docs/workflow/ を参照。

このファイルは薄い層に徹する:
  (1) 各モジュールを import して mcp にツール/リソースを登録する
  (2) テスト・外部が参照する公開面（server.<tool> / server.<helper> / server.mcp /
      server.arf_reader / server.AnalysisSession / server.KNOWLEDGE_DIR 等）を再エクスポートする
  (3) __main__ で mcp.run() する

`.mcp.json` / `.vscode/mcp.json` がこのファイルを絶対パスで指しているため、
ルートから動かしてはならない。

可変状態（DATA_DIR / KNOWLEDGE_DIR / ANALYSES_DIR / session）の**正準**は
metabolomix.core.mcp_core / metabolomix.core.session_state 側にある。ここでの再エクスポートは
読み取り用の束縛にすぎないため、差し替えは必ず正準モジュール
（mcp_core.DATA_DIR / mcp_core.KNOWLEDGE_DIR / mcp_core.ANALYSES_DIR /
session_state.session）に対して行うこと。
"""
import os

# patch.object(server.arf_reader, ...) が共有 module 経由で効くよう、module を公開する。
from metabolomix.arf import reader as arf_reader

# --- 設定・状態・ヘルパの再エクスポート（読み取り用） ---
from metabolomix.core.mcp_core import (
    mcp,
    BASE_DIR,
    OUTPUT_FORMAT_DOC,
    MCP_INSTRUCTIONS,
    KNOWLEDGE_DIR,
    PLAYBOOK_DIR,
    ANALYSES_DIR,
    _build_report_meta,
    _dir_is_writable,
    _first_writable_dir,
)
from metabolomix.core.session_state import AnalysisSession, _build_sample_meta
from metabolomix.core.path_resolvers import (
    resolve_arf_file_path,
    resolve_arf2_file_path,
    resolve_dcl_file_path,
    resolve_eicaef_file_path,
    resolve_pai2_file_path,
    _select_latest_batch,
    _describe_batch_selection,
    _filter_arf_spots,
)
from metabolomix.core.tool_helpers import (
    _identity_tables,
    _build_verification_dossier,
    _pca_scatter_arrays,
    _remember_arf_pca_plot,
    _pp_build_matrix,
    _pp_has_preprocessed,
)

# --- ツール/リソースの登録＋公開面の再エクスポート ---
# import 副作用で @mcp.tool / @mcp.resource が mcp に登録される。star import は各
# モジュールの __all__（＝そのモジュールのツール名）だけを取り込む。
# metabolomix.tools.dataset は arf/arf2 に依存するため最後に読み込む。
from metabolomix.tools import resources as tools_resources  # リソース 4 + テンプレート 3
from metabolomix.tools.objective import *
from metabolomix.tools.reports import *
from metabolomix.pai2.tools import *
from metabolomix.dcl.tools import *
from metabolomix.library.tools import *  # library_load, library_match_feature, library_plot_mirror
from metabolomix.arf.tools import *
from metabolomix.tools.samples import *  # sample_search
from metabolomix.arf2.tools import *
from metabolomix.eic.tools import *
from metabolomix.tools.dataset import *  # list_data_files, load_dataset
from metabolomix.tools.mztab_tools import *  # dataset_load, dataset_status
from metabolomix.tools.console_tools import *  # console_plan, console_run, console_status, job_list
from metabolomix.tools.dataset_analysis_tools import *  # dataset_preprocess, dataset_pca, dataset_differential, dataset_export_differential, dataset_set_sample_metadata
from metabolomix.tools.pipeline_tools import *  # pipeline_plan, pipeline_run, pipeline_status, pipeline_resume, pipeline_cancel
from metabolomix.tools.maintenance import *  # server_update
from metabolomix.tools.curation_tools import *  # curation_review, curation_suggest, curation_submit, curation_flags, curation_view_data


if __name__ == "__main__":
    # origin/main からの遅れを裏で 1 回だけ確認する（通知のみ・自動 pull はしない）。
    # import 副作用にしないのは、テストが server を import しただけで
    # ネットワークを叩き始めないようにするため。
    from metabolomix.core.version import start_update_check
    start_update_check()

    # 既定は stdio（ローカル開発: Claude がサブプロセスとして起動）。
    # 常駐させて HTTP で待受けるときは LIPIDMIX_TRANSPORT=streamable-http を設定する。
    transport = os.environ.get("LIPIDMIX_TRANSPORT", "stdio")
    mcp.run(transport=transport)
