"""全ツールが MCP 標準 annotations を宣言していることの回帰テスト。

汎用 MCP クライアントは annotations だけを見て承認要否とリプレイ安全性を決める。
annotations の無いツールはクライアント側で UNKNOWN 扱いになり、毎回承認待ちで
止まる。新しいツールを足したらここにも足すこと。

readOnlyHint の解釈: 「サーバの外に副作用が無い」＝ファイルとネットワークを
変更しない。サーバ自身の解析セッション状態の更新は副作用に数えない
（セッション状態の依存は missing_state エンベロープで伝えるため、
annotations で二重に表現しない）。
"""
import asyncio
import unittest

import server

READ_ONLY = {"readOnlyHint": True}
EXTERNAL = {"readOnlyHint": True, "openWorldHint": True}
LOCAL_WRITE = {"readOnlyHint": False, "destructiveHint": False, "idempotentHint": True}
# ファイルを書き、かつ冪等でない（同じ引数の再実行で結果が積み増される）。
# ingest_stage は _inbox に新規ノートを増やし、log_search は探索ログに行を追記し、
# update_objective(add_subquestions=...) は既存 Q の最大+1 で採番して小問を追記する。
LOCAL_WRITE_APPEND = {"readOnlyHint": False, "destructiveHint": False, "idempotentHint": False}
# ネットワークへ出たうえでローカルを書き換える。二度目は「最新です」で止まる → 冪等。
REMOTE_WRITE = {"readOnlyHint": False, "destructiveHint": False,
                "idempotentHint": True, "openWorldHint": True}
DESTRUCTIVE = {"readOnlyHint": False, "destructiveHint": True}
# ファイルを消すが、二度目は何も残っていないので同じ状態に落ち着く → 冪等。
LOCAL_DELETE = {"readOnlyHint": False, "destructiveHint": True, "idempotentHint": True}

EXPECTED_ANNOTATIONS = {
    # --- ARF（多サンプル解析）。セッション状態は更新するが外部副作用は無い ---
    "arf_list_tags": READ_ONLY,
    "arf_list_classes": READ_ONLY,
    "arf_list_sample_roles": READ_ONLY,
    "arf_exclude": READ_ONLY,
    "arf_preprocess": READ_ONLY,
    "arf_pca_preprocessed": READ_ONLY,
    "arf_parser": READ_ONLY,
    "arf_differential": READ_ONLY,
    "arf_plot_volcano": READ_ONLY,
    "arf2_parser": READ_ONLY,
    "arf2_annotate_identities": READ_ONLY,
    # --- データセット入口 ---
    "list_data_files": READ_ONLY,
    "load_dataset": READ_ONLY,
    "dataset_load": LOCAL_WRITE,
    "dataset_status": READ_ONLY,
    # --- DatasetState 解析層 ---
    # readOnlyHint はファイル/ネットワークへの副作用の有無。セッション状態の
    # 更新は数えない（このファイル冒頭の規約）。ARF の同一操作も READ_ONLY。
    "dataset_preprocess": READ_ONLY,
    "dataset_pca": READ_ONLY,
    "dataset_differential": READ_ONLY,
    # ファイルを書き、同じ引数なら同じ内容で上書きする → 冪等
    "dataset_export_differential": LOCAL_WRITE,
    # シートを読み検証してからセッション状態を一括反映する。外部副作用はセッション
    # 状態の更新のみ（ファイル書き込みは無い）が、同じシートの再送は同じ状態に
    # 落ち着く → destructive ではなく idempotent。
    "dataset_set_sample_metadata": LOCAL_WRITE,
    # v2統計は名指しした解析行列を読むだけ（行列も session も書き換えない。
    # 結果全量の保持は session.dataset.results への追加で、外界へは書かない）。
    "dataset_statistic": READ_ONLY,
    # 解析行列を1本作って session へ置くだけ（外界へは書かない）。
    "dataset_build_matrix": READ_ONLY,
    # --- Console 実行層（console_plan/run はジョブ状態・ファイルを生成する） ---
    "console_plan": LOCAL_WRITE_APPEND,
    "console_run": LOCAL_WRITE_APPEND,
    "console_status": READ_ONLY,
    # 入力フォルダとメソッドファイルを作る。同じ引数なら同じ結果 → 冪等
    "console_prepare_input": LOCAL_WRITE,
    "console_method_template": LOCAL_WRITE,
    # 読み取り専用だが、同じ入力なら結果は変わらない（idempotent）ことも明示している。
    "console_method_candidates": {"readOnlyHint": True, "destructiveHint": False, "idempotentHint": True},
    # ジョブが記録した生成物を消す。破壊的だが、二度目は何も残っていない → 冪等
    "console_cleanup": LOCAL_DELETE,
    "job_list": READ_ONLY,
    # --- pipeline系（生データフォルダ起点、Task18）。statusだけ読取専用。
    # 他4件は同一request_idの再送・二重cancelが同じ結果に落ち着くためidempotent。
    "pipeline_plan": LOCAL_WRITE,
    "pipeline_run": LOCAL_WRITE,
    "pipeline_status": READ_ONLY,
    "pipeline_resume": LOCAL_WRITE,
    "pipeline_cancel": LOCAL_WRITE,
    # --- DCL / PAI2 / EIC ---
    "dcl_parser": READ_ONLY,
    "dcl_find_msms": READ_ONLY,
    # --- 参照ライブラリ（MS/MS スペクトル照合） ---
    # library_load は SQLite キャッシュを書くので唯一 False。同じ元ファイルなら
    # 同じ store になる（sha256 キー）→ idempotent。
    "library_load": LOCAL_WRITE,
    "library_match_feature": READ_ONLY,
    "library_plot_mirror": READ_ONLY,
    "pai2_parser": READ_ONLY,
    "pai2_inspect_peak": READ_ONLY,
    "verify_peak_annotation": READ_ONLY,
    "eic_parser": READ_ONLY,
    "eic_plot_chromatograms": READ_ONLY,
    "eic_plot_compounds": READ_ONLY,
    "eic_rank_by_max_intensity": READ_ONLY,
    "eic_search_by_mz_range": READ_ONLY,
    "eic_search_by_rt_range": READ_ONLY,
    # --- 読み取り系 ---
    "knowledge_coverage": READ_ONLY,
    "read_report": READ_ONLY,
    "list_reports": READ_ONLY,
    "sample_search": READ_ONLY,
    "ingest_review_queue": READ_ONLY,
    # --- 外部ネットワーク ---
    "paper_search": EXTERNAL,
    # --- ローカル書き出し（同じ引数なら同じパスへ上書き＝idempotent） ---
    "record_objective": LOCAL_WRITE,
    "write_report": LOCAL_WRITE,
    "save_pca_figure": LOCAL_WRITE,
    "save_volcano_figure": LOCAL_WRITE,
    "save_eic_figure": LOCAL_WRITE,
    "arf_export_differential": LOCAL_WRITE,
    # --- knowledge の変更 ---
    "ingest_stage": LOCAL_WRITE_APPEND,
    "ingest_promote": DESTRUCTIVE,
    "ingest_reject": DESTRUCTIVE,
    # --- 追記書き込み（読み取り系ではない: 既存の最大値+1で採番して追記するため冪等でない） ---
    # update_objective は add_subquestions が既存 Q の最大+1 で採番して追記する。
    # log_search / ingest_stage と同じ「追記で結果が積み増される」理由でこの3件が
    # LOCAL_WRITE_APPEND を共有する。
    "update_objective": LOCAL_WRITE_APPEND,
    "log_search": LOCAL_WRITE_APPEND,
    # --- サーバ自身の更新（origin へ出て、作業ツリーを早送りする） ---
    "server_update": REMOTE_WRITE,
    # --- キュレーション。review は呼ぶたびに新しい review-<id>.json/.html を作る（console_plan /
    # console_run と同じ「毎回新しい成果物を生む」性質なので idempotent ではない）。
    # submit は flags.jsonl に追記する。 ---
    "curation_review": LOCAL_WRITE_APPEND,
    "curation_submit": LOCAL_WRITE_APPEND,
    "curation_flags": READ_ONLY,
    "curation_view_data": READ_ONLY,
}


def actual_annotations():
    tools = asyncio.run(server.mcp.list_tools())
    out = {}
    for tool in tools:
        ann = tool.annotations
        out[tool.name] = (
            None
            if ann is None
            else ann.model_dump(mode="json", by_alias=True, exclude_none=True)
        )
    return out


class ToolAnnotationTests(unittest.TestCase):
    def test_expectation_covers_every_registered_tool(self):
        self.assertEqual(set(actual_annotations()), set(EXPECTED_ANNOTATIONS))

    def test_every_tool_declares_annotations(self):
        missing = [name for name, ann in actual_annotations().items() if not ann]
        self.assertEqual(missing, [], f"annotations 未宣言のツール: {missing}")

    def test_annotations_match_the_expected_table(self):
        actual = actual_annotations()
        for name, expected in EXPECTED_ANNOTATIONS.items():
            self.assertEqual(actual[name], expected, name)

    def test_external_network_tools_also_declare_open_world(self):
        """openWorldHint が無いと汎用クライアントがネットワークゲートを迂回する。"""
        actual = actual_annotations()
        for name in ("paper_search", "server_update"):
            with self.subTest(tool=name):
                self.assertIs(actual[name]["openWorldHint"], True)


if __name__ == "__main__":
    unittest.main()
