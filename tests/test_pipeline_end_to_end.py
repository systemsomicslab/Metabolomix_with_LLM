"""fake Consoleによる一気通貫の受入検証（spec §13 D01〜D08・E01〜E03）。

ここでの偽物は**workerが起動するConsoleのコマンドだけ**。engine・loading・
metadata・preprocessing・differential・export・renderはすべて本物を通す。
"""
from __future__ import annotations

import json
import threading
import warnings
from pathlib import Path

import pytest

from metabolomix.analysis.export_contract import CONTRACT_VERSION, EXPORT_COLUMNS, LOG2FC_SIGN
from metabolomix.pipeline import recovery, store
from tests.pipeline_fixtures import (
    DEFAULT_COMPARISON, PipelineHarness, WAIT_TIMEOUT_S, read_contract_tsv,
)

#: 逆向きの比較（reference=treated, test=control）。同じデータで log2FC の符号が
#: 反転することを見て、「向き」のassertが定数を写しただけでないことを保証する。
_REVERSED_COMPARISON = {
    "comparison_id": "control_vs_treated",
    "reference_group": "treated",
    "test_group": "control",
}

#: `tests.pipeline_fixtures._FEATURES` が書く3特徴のInChIKey（全件が背景に残る）。
_EXPECTED_INCHIKEYS = {
    "IPCSVZSSVZVIGE-UHFFFAOYSA-N",
    "XKMRRTOUMJRJIA-UHFFFAOYSA-N",
    "DGGXCMYPQAOAJC-UHFFFAOYSA-N",
}


@pytest.fixture
def pipeline_harness(tmp_path, monkeypatch):
    harness = PipelineHarness(tmp_path, monkeypatch)
    try:
        yield harness
    finally:
        harness.close()


def test_resume_reuses_console_and_completes_exports(pipeline_harness):
    """spec D02/D03: 群未指定なら探索まで進んでneeds_input、シートと比較を足して
    再開したらConsoleを**再実行せず**比較・出力まで完了する。

    **成果物の存在だけでは足りない**（最終レビュー指摘1）。resumeで足した
    比較stageが`report`より後ろへ追記されると、レポートは比較を1件も
    実行していない時点で書かれ、TSVとvolcanoはディスク上に在るのに
    レポートは「missing」「(n/a)」と告げ、runはcompletedを名乗る。
    ここではレポート**本文**が実際の比較結果を載せていることまで見る。
    """
    run = pipeline_harness.start(target="differential", comparisons=[])
    waiting = pipeline_harness.wait(run, expected="needs_input")
    assert waiting["request"]["effective_target"] == "differential"
    assert pipeline_harness.console_start_count(run) == 1
    pipeline_harness.resume_with_groups(run)
    completed = pipeline_harness.wait(run, expected="completed")
    assert pipeline_harness.console_start_count(run) == 1
    tsv = pipeline_harness.output_path(completed, "tsv:treated_vs_control")
    assert Path(tsv).is_file()
    assert json.loads(Path(run).read_text(encoding="utf-8"))["status"] == "completed"

    # runが自分で「不完全なまま書いた」と記録していない。
    codes = {w["code"] for w in completed.get("warnings") or []}
    assert "REPORT_INCOMPLETE_AT_WRITE_TIME" not in codes, completed["warnings"]

    report_text = Path(
        pipeline_harness.output_path(completed, "quality_report")).read_text(encoding="utf-8")
    # 比較行がachieved（missingではない）で、向きもそのまま載っている。
    assert "| treated_vs_control | control | treated | achieved |" in report_text, report_text
    # InChIKey被覆が実数（(n/a)ではない）——TSVを実際に読めた証拠。
    assert "| treated_vs_control | 3 | 3 | 0 |" in report_text, report_text
    # サンプル来歴が実サンプルを名指ししている（最終レビュー指摘2）。
    provenance = report_text.split(
        "## Role / Group / Batch / Order Provenance", 1)[1].split("\n## ", 1)[0]
    assert "(no sample manifest recorded)" not in provenance, provenance
    sample_rows = [line for line in provenance.splitlines() if line.startswith("| S")]
    assert len(sample_rows) == 8, provenance
    assert provenance.count("| control |") == 4, provenance
    assert provenance.count("| treated |") == 4, provenance
    # レポートは比較stageのあとに書かれている（工程順そのものの確認）。
    stages = completed["stages"]
    assert stages["report"]["updated_at"] >= stages["export:treated_vs_control"]["updated_at"]


# briefの受入表は5シナリオすべての期待値をpartialとしているが、実装とspec §9.1の
# 状態機械はここで一致して`failed`を出す——「running --> partial: 有効な出力を残して
# 回復不能な失敗 / running --> failed: 有効な解析出力なし」であり、上流が完了検証を
# 通らなかった時点で`record["results"]`は空（prepare_inputは成果物refを持たない）。
# Consoleが残した中間ファイルは`record["results"]`ではなくConsole jobのartifactsで、
# 「有効な解析出力」ではない。**briefとspecの食い違いはreportへ明記**し、ここでは
# 「ファイルがあるから成功」と読み替えない側（spec §9.1）を採る。
# partialが本当に出る経路は`test_zero_inchikey_*`（有効な探索結果を残したまま
# 必須TSVだけ作れない）で別途検証する。
@pytest.mark.parametrize("scenario, expected, termination, exit_code, error_codes, retained", [
    ("success", "completed", "exited", 0, [], None),
    ("nonzero", "failed", "exited", 1, [], "intermediate.pai2"),
    ("invalid", "failed", "exited", 0, ["MZTAB_STRUCTURE_INVALID"], None),
    ("missing_sample", "failed", "exited", 0, ["SAMPLE_MAPPING_MISSING"], None),
    ("hang", "failed", "timeout", None, [], "intermediate.pai2"),
])
def test_execution_scenarios(pipeline_harness, scenario, expected, termination,
                             exit_code, error_codes, retained):
    """spec A01〜A04: 異常終了・不正mzTab・assay欠落・timeoutの終端状態と証跡。"""
    run = pipeline_harness.start_scenario(scenario, timeout_s=3)
    record = pipeline_harness.wait(run, expected=expected)
    assert record["status"] == expected

    receipt = pipeline_harness.console_receipt(run)
    assert receipt["termination"] == termination
    if exit_code is None:
        assert receipt["exit_code"] != 0, "timeoutを0で補完してはいけない"
    else:
        assert receipt["exit_code"] == exit_code
    for code in error_codes:
        assert code in receipt["validation"]["errors"], receipt["validation"]

    if expected == "completed":
        assert record["stages"]["upstream"]["status"] == "succeeded"
        assert record["upstream"]["verification"]["status"] == "completed"
        return

    # 「終了検証が失敗し、成果物は保持された」——ファイルの有無を成功と読み替えない。
    assert record["stages"]["upstream"]["status"] == "failed"
    assert record["stages"]["upstream"]["error"]["code"] == "MSDIAL_EXECUTION_FAILED"
    assert record["upstream"]["verification"]["status"] != "completed"
    # 下流は自動進行しない（spec A03）。
    assert record["stages"]["load_dataset"]["status"] == "pending"
    assert record["stages"]["report"]["status"] == "pending"
    # 収集は飛ばさない。中間ファイルはjobのartifactsとして残っている。
    assert receipt["collection"]["status"] == "succeeded"
    if retained:
        artifacts = pipeline_harness.console_job(run).artifacts
        assert [a for a in artifacts if a.path.endswith(retained)], \
            f"{retained} が保持されていません: {[a.path for a in artifacts]}"
    # `failed`（`partial`ではない）である根拠: 有効な解析出力が1件も無い。
    assert record["results"] == []


# ---------- 最終レビュー指摘3: rerun_upstreamは前回attemptの証跡を消さない ----------

def test_rerun_upstream_creates_a_new_attempt_and_keeps_the_old_evidence(pipeline_harness):
    """spec §9.1/§9.3・D09: Consoleの再試行は**新しいjob/attempt**を作る。

    timeoutで打ち切ったConsole実行の終了証跡（`execution-result.json`）・
    Console自身のログ（`msdial.log`は`_CREATE_ALWAYS`＝切り詰めで開かれる）・
    再収集に要る実行前スナップショット（`worker.json`の`recovery.befores`）は、
    再実行後も1バイトも変わっていないこと。run_dirがattemptに依らない固定
    パスだと、これらは上書き・切り詰め・置換され「timeoutした事実」を
    後から証明できなくなる。
    """
    from metabolomix.pipeline import service

    run = pipeline_harness.start_scenario("hang", timeout_s=3)
    pipeline_harness.wait(run, expected="failed")

    first_dir = pipeline_harness.console_run_dir(run)
    first_receipt = json.loads((first_dir / "execution-result.json").read_text(encoding="utf-8"))
    assert first_receipt["termination"] == "timeout"
    first_log = (first_dir / "msdial.log").read_bytes()
    first_state = json.loads((first_dir / "worker.json").read_text(encoding="utf-8"))
    assert first_state["recovery"]["befores"], first_state

    pipeline_harness.set_launch_options(run, scenario="success")
    service.resume_pipeline(pipeline_harness.pipeline_root(run), rerun_upstream=True)
    completed = pipeline_harness.wait(run, expected="completed")

    # 明示的な再実行なのでConsoleは2回起動した（自動再試行ではない）。
    assert pipeline_harness.console_start_count(run) == 2
    second_dir = pipeline_harness.console_run_dir(run)
    assert second_dir != first_dir, "再実行が前回と同じrun_dirを使っています"
    assert Path(completed["upstream"]["console_job_path"]).parent == second_dir

    # 前回attemptの証跡はそのまま残っている。
    assert json.loads(
        (first_dir / "execution-result.json").read_text(encoding="utf-8")) == first_receipt
    assert (first_dir / "msdial.log").read_bytes() == first_log
    assert json.loads(
        (first_dir / "worker.json").read_text(encoding="utf-8")) == first_state
    # 新しいattemptは自分の証跡を別に持つ。
    assert json.loads(
        (second_dir / "execution-result.json").read_text(encoding="utf-8"))["termination"] \
        == "exited"


# ---------- E02: 出力TSVをパースして照合する ----------

def test_exported_tsv_matches_the_contract_in_both_directions(pipeline_harness):
    """spec E02: 全背景行・有限値/空欄・15列・方向・result_idの整合。

    向きは定数の写しでは確かめられないので、**同じデータの逆向き比較**を同時に
    走らせ、log2FCが符号だけ反転し絶対値が一致することまで見る
    （`export_contract.LOG2FC_SIGN` = 正ならtest群が高い）。
    """
    pipeline_harness.write_manifest()
    run = pipeline_harness.start(
        target="differential", comparisons=[DEFAULT_COMPARISON, _REVERSED_COMPARISON],
        extra_request={"sample_manifest": "sample-manifest.tsv"})
    record = pipeline_harness.wait(run, expected="completed")

    forward_meta, forward_rows, fieldnames = read_contract_tsv(
        pipeline_harness.output_path(record, "tsv:treated_vs_control"))
    reverse_meta, reverse_rows, reverse_fields = read_contract_tsv(
        pipeline_harness.output_path(record, "tsv:control_vs_treated"))

    # 15列・列順そのもの（別リポジトリ massbank-context との契約）。
    assert fieldnames == EXPORT_COLUMNS
    assert reverse_fields == EXPORT_COLUMNS
    assert forward_meta["contract_version"] == str(CONTRACT_VERSION)
    assert forward_meta["log2fc_sign"] == LOG2FC_SIGN

    # 全背景行が残る（有意なものだけに絞らない）。
    assert len(forward_rows) == 3
    assert {row["inchikey"] for row in forward_rows} == _EXPECTED_INCHIKEYS
    assert forward_meta["n_features_total"] == "3"
    assert forward_meta["n_with_inchikey"] == "3"
    assert forward_meta["n_unannotated"] == "0"

    # 群とその件数。`group_b`（=test_group）が高い方向が正。
    assert (forward_meta["group_a"], forward_meta["n_a"]) == ("control", "4")
    assert (forward_meta["group_b"], forward_meta["n_b"]) == ("treated", "4")
    assert (reverse_meta["group_a"], reverse_meta["n_a"]) == ("treated", "4")
    assert (reverse_meta["group_b"], reverse_meta["n_b"]) == ("control", "4")

    # 合成強度はassay番号とともに単調増加する＝treatedが必ず高い。
    forward_by_key = {row["inchikey"]: float(row["log2fc"]) for row in forward_rows}
    reverse_by_key = {row["inchikey"]: float(row["log2fc"]) for row in reverse_rows}
    assert all(value > 0 for value in forward_by_key.values()), forward_by_key
    assert all(value < 0 for value in reverse_by_key.values()), reverse_by_key
    for key, value in forward_by_key.items():
        assert value == pytest.approx(-reverse_by_key[key], rel=1e-9)

    # 空欄は空欄のまま（「該当なし」ではなく「この経路では取得していない」）。
    for row in forward_rows:
        assert row["ontology"] == ""
        assert row["msi_level"] == ""
        assert row["inchikey_source"] == "database_identifier"
        # 名前の実体は SMF 行ではなく SME 行。旧 "mztab_smf" は事実と違った
        # （spec 2026-09-17-mztab-sml-annotation-design §7.3）。
        assert row["name_source"] == "mztab_sme"
        assert row["significant"] in ("true", "false")
        assert float(row["mz"]) > 0 and float(row["rt"]) > 0
        assert float(row["mean_b"]) > float(row["mean_a"])

    # 2つの比較のTSVが互いに取り違えられていない（result_idも数値も別物）。
    assert forward_meta["result_id"] != reverse_meta["result_id"]
    assert forward_meta["preprocess_id"] == reverse_meta["preprocess_id"]


# ---------- E03: PCA/volcano PNGを描画・確認 ----------

def test_saved_figures_come_from_the_named_results(pipeline_harness):
    """spec E03: 図が「指定データ・群・出所」と一致し、実PNGとして残ること。"""
    from metabolomix.plots.result_output import save_result_figure
    from metabolomix.plots.volcano import render_volcano_plot

    pipeline_harness.write_manifest()
    run = pipeline_harness.start(
        target="differential", comparisons=[DEFAULT_COMPARISON],
        extra_request={"sample_manifest": "sample-manifest.tsv"})
    record = pipeline_harness.wait(run, expected="completed")
    pipeline_root = Path(record["identity"]["pipeline_root"])

    volcano_ref = pipeline_harness.output_ref(record, "volcano:treated_vs_control")
    pca_figure_ref = pipeline_harness.output_ref(record, "pca_figure")
    pca_ref = pipeline_harness.output_ref(record, "pca")

    # 出所: 図refは「どの結果から描いたか」を親IDで名指ししている。その親IDが
    # record["results"]の中で解決できるかは
    # `test_persisted_result_graph_is_joinable_by_result_id`（既知の欠陥）が見る。
    assert pca_figure_ref["parent_ids"] and volcano_ref["parent_ids"]
    assert pca_figure_ref["result_id"].startswith("res_pca_figure_")
    assert volcano_ref["result_id"] == "res_volcano_treated_vs_control"
    assert pca_ref["kind"] == "pca" and volcano_ref["kind"] == "volcano"

    # 実PNGとして残り、記録済みhashと一致する（改竄・破損していない）。
    for ref in (volcano_ref, pca_figure_ref):
        data = (pipeline_root / ref["relative_path"]).read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n"
        assert len(data) > 1000
    assert store.verify_result_refs(pipeline_root, [volcano_ref, pca_figure_ref]) == []

    # 群ラベルと点: 保存された差次的結果そのものを描き直し、図が名指しする群と
    # 実際に描かれる点数を見る（3特徴＝3点）。
    result = pipeline_harness.result_data(record, "differential:treated_vs_control")
    assert (result["a"], result["b"]) == ("control", "treated")
    figure = render_volcano_plot(result)
    try:
        assert figure.axes[0].get_title() == "Volcano (control vs treated)"
        assert figure.axes[0].get_xlabel() == "log2 fold change"
        drawn = sum(len(collection.get_offsets()) for collection in figure.axes[0].collections)
        assert drawn == 3, "volcanoの点が描かれていません"
    finally:
        import matplotlib.pyplot as plt
        plt.close(figure)

    # 描き直しても警告が1件も出ない（豆腐＝`Glyph missing`を含む）。
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        redrawn = save_result_figure(None, result, pipeline_harness.tmp_path / "redraw.png",
                                     kind="volcano")
    assert caught == [], [str(w.message) for w in caught]
    assert redrawn.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def _scatter_ink(png_path) -> int:
    """PNGの中で「散布図のマーカー色」を持つ画素数を数える。

    matplotlibの既定の点の色はC0（青系 #1f77b4）で、軸枠・目盛・文字・注記は
    黒/灰/赤系しか使わない。したがって**青が赤より明確に強い画素**の数は、
    実際に点が描かれたかどうかだけを見る指標になる（枠だけの空図なら0）。
    その「空図なら0」自体もテスト内で実測して確かめる（下記の参照図）。
    """
    import matplotlib.image as mpimg

    image = mpimg.imread(str(png_path))
    return int(((image[..., 2] - image[..., 0]) > 0.1).sum())


def _blank_reference_png(path) -> int:
    """点を一つも描かないPCA図（軸・ラベル・題だけ）の散布図インク量。"""
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(figsize=(8, 6))
    try:
        axes.set_xlabel("PC1")
        axes.set_ylabel("PC2")
        axes.set_title("PCA")
        figure.savefig(str(path), format="png", dpi=120, bbox_inches="tight")
    finally:
        plt.close(figure)
    return _scatter_ink(path)


def test_saved_pca_figure_plots_the_run_samples(pipeline_harness):
    """spec E03: PCA図はこのrunの検体を実際に描いていること。

    「PNGが在る・hashが合う」だけでは空の散布図を見逃す（Task19が見つけた欠陥は
    まさにそれで、`results/pca.png`は枠だけの図なのにrunはcompletedになっていた）。
    ここでは (1) 永続化されたPCA結果自身が検体の座標を持つこと、(2) 保存された
    PNGに**点の画素が実在する**こと、(3) その結果から本番の保存経路で描き直しても
    点が描かれ、警告（豆腐＝`Glyph missing`）が1件も出ないことまで見る。
    """
    from metabolomix.core.tool_helpers import _pca_scatter_arrays, dataset_pca_plot
    from metabolomix.plots.result_output import save_result_figure

    run = pipeline_harness.start(target="exploratory")
    record = pipeline_harness.wait(run, expected="completed")
    pipeline_root = Path(record["identity"]["pipeline_root"])

    # (1) 保存されたPCA結果が、この8検体の座標を持っている。
    pca_result = pipeline_harness.result_data(record, "pca")
    xs, ys, labels, x_label, _y_label, _title = _pca_scatter_arrays(
        dataset_pca_plot(pca_result))
    assert len(xs) == 8, f"PCA結果に検体の座標がありません: {sorted(pca_result)}"
    assert len(ys) == 8
    assert labels == [f"S{i}" for i in range(1, 9)]
    assert x_label.startswith("PC1 ("), x_label   # 寄与率つきの軸ラベル

    # (2) 実際に保存されたPNGに点が描かれている（空の枠ではない）。
    figure_ref = pipeline_harness.output_ref(record, "pca_figure")
    assert store.verify_result_refs(pipeline_root, [figure_ref]) == []
    blank_ink = _blank_reference_png(pipeline_harness.tmp_path / "blank-pca.png")
    assert blank_ink == 0, "指標そのものが壊れています（空図でも点の色を数えている）"
    saved_ink = _scatter_ink(pipeline_root / figure_ref["relative_path"])
    assert saved_ink > 100, f"results/pca.pngが空の散布図です（点の画素数={saved_ink}）"

    # (3) 永続化された結果から本番の保存経路で描き直しても同じ図になる。
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        redrawn = save_result_figure(None, pca_result,
                                     pipeline_harness.tmp_path / "redraw-pca.png", kind="pca")
    assert caught == [], [str(w.message) for w in caught]
    assert redrawn.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert _scatter_ink(redrawn) > 100


# ---------- E01: InChIKey 0件 ----------

def test_zero_inchikey_stays_partial_with_export_background_empty(pipeline_harness):
    """spec E01: 必須TSVを作れないならpartial。有効な探索・比較結果は保持する。"""
    pipeline_harness.write_manifest()
    run = pipeline_harness.start(
        target="differential", comparisons=[DEFAULT_COMPARISON],
        extra_request={"sample_manifest": "sample-manifest.tsv"},
        no_inchikey=True)
    record = pipeline_harness.wait(run, expected="partial")

    codes = {w["code"] for w in record["warnings"]}
    assert "EXPORT_BACKGROUND_EMPTY" in codes
    # 完了へ格下げしない。TSVだけが欠け、他は達成済み。
    assert record["status"] == "partial"
    names = {ref["output_name"] for ref in record["results"]}
    assert "tsv:treated_vs_control" not in names
    assert {"preprocess", "pca", "pca_figure", "differential:treated_vs_control",
            "volcano:treated_vs_control", "quality_report"} <= names
    # InChIKey無しの行を内部の差次的結果から消していない（3特徴すべて残る）。
    result = pipeline_harness.result_data(record, "differential:treated_vs_control")
    assert len(result["results"]) == 3
    # stage自体は成功（成果物の欠落は最終判定＝evaluate_targetが捉える）。
    assert record["stages"]["export:treated_vs_control"]["status"] == "succeeded"


# ---------- GUI project未達 ----------

def test_missing_gui_project_warns_and_blocks_completed(pipeline_harness):
    """save_project=trueなのに.mdprojectが無ければ、必須出力欠落でcompletedにしない。"""
    run = pipeline_harness.start(target="exploratory", save_project=True)
    record = pipeline_harness.wait(run, expected="partial")

    codes = {w["code"] for w in record["warnings"]}
    assert "GUI_PROJECT_UNAVAILABLE" in codes
    names = {ref["output_name"] for ref in record["results"]}
    assert "gui_project" not in names
    assert {"preprocess", "pca", "pca_figure"} <= names
    # 上流自体は検証済み完了（GUI projectだけが欠けている）。
    assert record["stages"]["upstream"]["status"] == "succeeded"
    assert record["stages"]["validate_outputs"]["status"] == "succeeded"


# ---------- D08: ディレクトリ形式raw ----------

def test_directory_style_raw_runs_end_to_end_without_touching_the_source(pipeline_harness):
    """spec D08: 1測定単位がフォルダ（Agilent `.d` 等）でも規定どおり準備して進む。"""
    source = pipeline_harness.source("dirsource", extension=".d", kind="dir")

    def snapshot() -> dict:
        # `runs/`（pipelineの書込先）は元データではないので除く。
        return {str(path.relative_to(source["root"])): path.stat().st_mtime_ns
                for path in sorted(Path(source["root"]).rglob("*"))
                if "runs" not in path.relative_to(source["root"]).parts}

    before = snapshot()

    run = pipeline_harness.start(target="exploratory", source_name="dirsource")
    record = pipeline_harness.wait(run, expected="completed")

    # 1フォルダ＝1測定単位として配置され、mzTabのassayも8件になる。
    entries = [e for e in record["inputs"]["entries"] if e["role"] == "primary"]
    assert [e["kind"] for e in entries] == ["dir"] * 8
    staged = Path(record["identity"]["pipeline_root"]) / "input"
    assert sorted(p.name for p in staged.iterdir()) == [f"S{i}.d" for i in range(8)]
    assert all((staged / f"S{i}.d" / "AcqData.bin").is_file() for i in range(8))
    # 元データは読むだけ（移動・上書き・削除をしない）。
    assert snapshot() == before


# ---------- D04: 同一要求の同時送信 ----------

def test_identical_request_resend_reuses_the_run_and_starts_console_once(pipeline_harness):
    """spec D04: 同一request再送でConsoleは1回だけ、runも増えない。"""
    first = pipeline_harness.start(target="exploratory", request_id="req-1")
    record = pipeline_harness.wait(first, expected="completed")

    second = pipeline_harness.start(target="exploratory", request_id="req-1")
    third = pipeline_harness.start(target="exploratory")  # request_id無しの再送

    assert second == first and third == first
    assert pipeline_harness.console_start_count(first) == 1
    # 再利用したrunの上に2つ目のworkerを起こさない（＝completedのまま）。
    assert pipeline_harness.record(first)["status"] == "completed"
    assert pipeline_harness.record(first)["identity"]["pipeline_id"] == \
        record["identity"]["pipeline_id"]
    runs_parent = Path(record["identity"]["pipeline_root"]).parent
    assert len(list(runs_parent.iterdir())) == 1


def test_two_concurrent_starts_of_the_same_request_launch_console_once(pipeline_harness):
    """spec D04 同時送信: 同一requestの`start_pipeline`が本当に同時に来ても、
    runは1つ・Consoleは1回だけ。

    直前のテストは完了**後**の再送——順次実行なので、D04の「同時」を検証して
    いない（brief追加ケース「同一request同時送信」、spec §13 D04「同時resume」の
    対）。`find_or_create_run`はsource_root単位の`file_lock`で受付そのものを
    直列化するので、同じ`pipeline_path`が返ること自体は`tests/test_pipeline_
    store.py::test_two_processes_accept_concurrently_return_same_run`が既に
    実プロセスで検証済み。ここで検証したいのはその**先**——`start_pipeline`が
    「読み出したstatusがplannedだから起動する」と判断する箇所は`file_lock`の
    **外**にあり、2つの呼び出しがどちらも新規作成直後の`planned`を読めば
    どちらも起動しうる（`metabolomix/pipeline/service.py`の`start_pipeline`docstring
    が名指す「レビュー指摘1」の対象コードそのもの）。実際に二重起動を止めて
    いるのはworker側のowner lock（`engine.run_engine`のworker.lock）である
    ことを、その場しのぎでなく実際の競合で確かめる。

    `store.find_or_create_run`の直後に`threading.Barrier(2)`を挟み、2つの
    `start_pipeline`呼び出しが**どちらも受付を終えてから**起動判定へ進む
    瞬間を強制する。バリアは2者そろわない限り解放されない（各回`timeout`
    超過で`BrokenBarrierError`）——つまりこのテストが例外なく通ること自体が
    「2つの呼び出しが同じ瞬間に起動判定の入口に立っていた」ことの証明になる。
    """
    barrier = threading.Barrier(2, timeout=WAIT_TIMEOUT_S)
    real_find_or_create_run = store.find_or_create_run

    def _synced_find_or_create_run(*args, **kwargs):
        result = real_find_or_create_run(*args, **kwargs)
        # 両方の呼び出しがここへ揃うまで、どちらも起動判定（load_run→status
        # 確認→launch）へ進めない。揃わなければBrokenBarrierErrorで例外送出。
        barrier.wait()
        return result

    pipeline_harness.monkeypatch.setattr(
        store, "find_or_create_run", _synced_find_or_create_run)

    # `PipelineHarness.source()`は`self._sources`への無防備なcheck-then-act。
    # 入力identityはSTATフィンガープリント（サイズ/mtime、内容ではない）なので、
    # 2つのスレッドが初回呼出しとして同時に`start()`（内部で`source("source")`を
    # 呼ぶ）へ入ると、生成競合で2つの実行ディレクトリを割り当てて本テストの
    # 前提（同一入力での同時送信）を偽陽性に壊しうる。事前にメインスレッドから
    # 1回呼び、実質的な生成をスレッド生成前に済ませておく（隣接する同時resume
    # テストが`start()`を先に順次呼んで同じ効果を得ているのと同じ発想）。
    pipeline_harness.source("source")

    results: list = []
    errors: list = []

    def _call():
        try:
            results.append(
                pipeline_harness.start(target="exploratory", request_id="req-simul"))
        except Exception as exc:  # pragma: no cover - 失敗時の診断用
            errors.append(exc)

    threads = [threading.Thread(target=_call) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=WAIT_TIMEOUT_S)

    assert not errors, f"同時起動でエラー（バリアが揃わなかった可能性）: {errors}"
    assert len(results) == 2
    first, second = results
    assert first == second, "2つのstart_pipeline呼び出しが別のrunを指しました"

    # 起動判定そのものが競合したこと（＝どちらも起動を試みたこと）を先に見る。
    workers = pipeline_harness.worker_processes(first)
    assert len(workers) == 2, \
        "2つのstart_pipeline呼び出しが2つのworker起動を試みていません（レースが起きていない）"

    record = pipeline_harness.wait(first, expected="completed")
    assert pipeline_harness.console_start_count(first) == 1
    runs_parent = Path(record["identity"]["pipeline_root"]).parent
    assert len(list(runs_parent.iterdir())) == 1

    # owner lockが片方を`PIPELINE_ALREADY_RUNNING`で拒否するか、遅れて着いた
    # 側がstageを何もせず合流するか（先着が既に全stageを終えていた場合）は
    # 実プロセスの起動・import待ち時間の実測値次第でどちらもありうる——
    # どちらであっても異常終了ではない（returncodeは0か1のいずれかに限る）。
    # ここで縛るべき不変条件は「Console起動が合計1回」であって、2つのworker
    # 間でどちらが実際に停止させられたかではない。
    returncodes = {w.returncode for w in workers}
    assert returncodes <= {0, 1}, f"想定外の終了コード: {returncodes}"


def test_two_concurrent_resumes_of_the_same_run_launch_console_once(pipeline_harness):
    """spec D04 同時resume: 中断済みrunへの`resume_pipeline`が同時に来ても、
    Consoleは合計1回のまま増えない。

    `recovery.prepare_resume`は`state_revision`の楽観排他で**書込みだけ**を
    直列化する（衝突した側は最新を読み直して自分のiterationをやり直す）。
    しかし「statusがplannedになったから起動する」という判断は
    `service.resume_pipeline`側が各呼び出しで独立に行っており、`start_pipeline`
    と全く同じ形の隙間になる——つまりこの2つは**同じ保護機構（worker側の
    owner lock）に帰着する**、独立な仕組みではない。上のテストと対にして
    ここでも実際に確かめる。

    合成: `target=differential`・比較未指定で`needs_input`にした後、シートを
    置いて2つの`resume_pipeline`を同時に投げる（`request_id`は与えない——
    冪等性キャッシュに頼らず、`prepare_resume`のCASそのものを競合させる）。
    `recovery.prepare_resume`の直後にバリアを挟み、上のテストと同じ論法で
    「両方が受付処理を終えてから起動判定に入る」瞬間を強制する。
    """
    cid = DEFAULT_COMPARISON["comparison_id"]
    run = pipeline_harness.start(target="differential", comparisons=[])
    pipeline_harness.wait(run, expected="needs_input")
    assert pipeline_harness.console_start_count(run) == 1
    manifest = pipeline_harness.write_manifest()

    barrier = threading.Barrier(2, timeout=WAIT_TIMEOUT_S)
    real_prepare_resume = recovery.prepare_resume

    def _synced_prepare_resume(*args, **kwargs):
        result = real_prepare_resume(*args, **kwargs)
        barrier.wait()
        return result

    pipeline_harness.monkeypatch.setattr(
        recovery, "prepare_resume", _synced_prepare_resume)

    pipeline_root = pipeline_harness.pipeline_root(run)
    updates = {"sample_manifest": manifest.name, "comparisons": [DEFAULT_COMPARISON]}
    results: list = []
    errors: list = []

    def _call():
        from metabolomix.pipeline import service
        try:
            results.append(service.resume_pipeline(pipeline_root, updates=updates))
        except Exception as exc:  # pragma: no cover - 失敗時の診断用
            errors.append(exc)

    threads = [threading.Thread(target=_call) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=WAIT_TIMEOUT_S)

    assert not errors, f"同時resumeでエラー（バリアが揃わなかった可能性）: {errors}"
    assert len(results) == 2

    workers = pipeline_harness.worker_processes(run)
    assert len(workers) == 3, "初回の起動1つ＋同時resumeが起こす2つのはず"
    resume_workers = workers[1:]

    record = pipeline_harness.wait(run, expected="completed")
    # upstreamはこのpassでskipされる（stage_inputs_unchanged）ので、resumeの
    # 競合があってもConsole起動は初回の1回のまま増えない。
    assert pipeline_harness.console_start_count(run) == 1
    tsv = pipeline_harness.output_path(record, f"tsv:{cid}")
    assert Path(tsv).is_file()

    # 上のstart版と同じ理由で、どちらのworkerが実際に拒否されたかは実測時間
    # 次第——縛るのは「異常終了していない」ことと「Console起動が増えない」こと。
    returncodes = {w.returncode for w in resume_workers}
    assert returncodes <= {0, 1}, f"想定外の終了コード: {returncodes}"


# ---------- D05: 入力・メソッドの変化 ----------

def test_changed_raw_is_detected_on_resume(pipeline_harness):
    """spec D05: 元rawが変わっていたら、無条件に再利用せずINPUT_CHANGEDで止める。"""
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import service

    run = pipeline_harness.start(target="exploratory")
    pipeline_harness.wait(run, expected="completed")

    source_root = Path(pipeline_harness.source("source")["root"])
    (source_root / "S3.wiff").write_text("tampered-raw-content", encoding="ascii")

    with pytest.raises(DomainError) as excinfo:
        service.resume_pipeline(pipeline_harness.pipeline_root(run),
                                updates={"target": "exploratory"})
    assert excinfo.value.code == "INPUT_CHANGED"
    assert "S3.wiff" in excinfo.value.details["changed"]


def test_changed_effective_method_is_detected_on_resume(pipeline_harness):
    """spec D05: Consoleが実際に読む実効メソッドの改変も再開時に検出する。"""
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import service

    run = pipeline_harness.start(target="exploratory")
    record = pipeline_harness.wait(run, expected="completed")

    effective = (Path(record["identity"]["pipeline_root"])
                 / record["inputs"]["method"]["effective_relative_path"])
    effective.write_text(effective.read_text(encoding="ascii") + "Extra: 1\n",
                         encoding="ascii")

    with pytest.raises(DomainError) as excinfo:
        service.resume_pipeline(pipeline_harness.pipeline_root(run),
                                updates={"target": "exploratory"})
    assert excinfo.value.code == "STAGED_INPUT_MISMATCH"


def test_persisted_result_graph_is_joinable_by_result_id(pipeline_harness):
    """spec §6.1/§9.1: `parent_ids`は`record["results"]`の`result_id`で解決できること。

    `result_state.make_provenance` は呼ばれるたびに新しいUUIDを発行するので、
    「同じ数値だから同じ結果」にはならない。下流（別リポジトリ massbank-context）は
    TSVメタ行のresult_idでpipelineの記録と突き合わせるため、ここが切れていると
    「どの解析から出たTSVか」を機械的に辿れない。
    """
    pipeline_harness.write_manifest()
    run = pipeline_harness.start(
        target="differential", comparisons=[DEFAULT_COMPARISON],
        extra_request={"sample_manifest": "sample-manifest.tsv"})
    record = pipeline_harness.wait(run, expected="completed")

    known_ids = {ref["result_id"] for ref in record["results"]}
    dangling = {
        ref["output_name"]: [pid for pid in ref.get("parent_ids") or [] if pid not in known_ids]
        for ref in record["results"]
    }
    dangling = {name: ids for name, ids in dangling.items() if ids}

    differential_ref = pipeline_harness.output_ref(record, "differential:treated_vs_control")
    meta, _rows, _fields = read_contract_tsv(
        pipeline_harness.output_path(record, "tsv:treated_vs_control"))

    assert dangling == {}, f"record['results']で解決できない親ID: {dangling}"
    assert meta["result_id"] == differential_ref["result_id"]
