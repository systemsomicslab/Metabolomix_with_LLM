"""MCP コア: FastMCP インスタンス・共通設定・状態ディレクトリ・レポート先解決。

このモジュールは依存グラフの **leaf**（stdlib / FastMCP / metabolomix.core.user_config のみ）。
metabolomix.<形式>.tools や metabolomix.tools.* を import してはならない（循環回避の絶対ルール）。

`BASE_DIR` はリポジトリルート（このファイルの 2 階層上）。docs/ はこれを起点に解決する。解析の記録（objective とレポート）はデータセットの解析フォルダ
（`DatasetState.analysis_dir`）配下の reports/ に置き、リポジトリ側には書かない。
データの探索先（旧 `DATA_DIR`）は無い——ツールは直近のデータセットのファイルを使う。
"""
import os
from datetime import date as _date
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from metabolomix.core import user_config
from metabolomix.core.atomic_io import DomainError

# このファイルは <root>/metabolomix/core/ にある。docs/
# はリポジトリルート基準で解決するため 2 階層上る。
BASE_DIR = Path(__file__).resolve().parents[2]

# output-format はトピック別に分割してある。一枚岩（約700行）を毎回 pull させると
# 解釈に不要な節まで文脈を食うため、共通核（core）＋パーサ/ツール別トピックに割り、
# 各ツールが自分のトピックだけを指す。節番号は分割前の通し番号を保持している。
OUTPUT_FORMAT_DIR = BASE_DIR / "docs" / "output_format"
OUTPUT_FORMAT_SECTIONS: dict[str, str] = {
    "core": "共通オントロジー・脂質名文法・LLM必須注意（全ツール共通。最初に読む）",
    "dataset": "MS-DIAL 解析フォルダ（.arf/.arf2）から読んだデータセット。特徴×試料の行・gap-fill と detected_mask・.arf2 の代表と同定の列・試料の役割・前処理と差次的解析",
    "dcl": ".dcl（デコンボリューション済み MS/MS）と dcl_* ツール",
    "eic": "dataset_plot の EIC の図（kind=eic / eic_compounds）。.EIC.aef の peak_top と描画契約",
    "identity": "同定信頼度・名前正規化・MSI レベル",
    "mztab": "mzTab-M から読んだデータセット。SME/SML の同定、特徴表、差次的エクスポートの同定列",
    "library": "参照ライブラリ（.dbs/.msp）照合。スコアの意味、-1/0 の区別、意図的に写した瑕疵",
    "curation": "アラインメントのキュレーション（curation_review 等）。判定・理由コードの意味",
    "pathway": "パスウェイ投影（pathway_project / pathway_plot、MSLipidMapper 連携）",
    "reanalysis": "再解析のデータセットのハンドル（load_dataset / dataset_list / dataset_close）・統計と行列の spec・MS/MS の合否・証拠の束",
    "quant": "定量（dataset_quant_* / dataset_merge_polarities / dataset_select_features）。gap-fill の印・単位・needs_input・記録",
}
# 共通核。`lipidmix://docs/output-format` が返す本体。
OUTPUT_FORMAT_DOC = OUTPUT_FORMAT_DIR / "core.md"


def output_format_section_path(topic: str) -> Path:
    """トピック名から本文ファイルのパスを返す。未知/不正なトピックは ValueError。

    topic は MCP リソースの URI 変数（＝呼び出し側由来）なので、宣言済みの名前と
    完全一致するものだけを通す。パス組み立ての唯一の口にして traversal を封じる。
    """
    if topic not in OUTPUT_FORMAT_SECTIONS:
        valid = ", ".join(OUTPUT_FORMAT_SECTIONS)
        raise ValueError(
            f"未知の output-format トピックです: {topic!r}。有効なトピック: {valid}"
        )
    return OUTPUT_FORMAT_DIR / f"{topic}.md"


def _dir_is_writable(directory: Path) -> bool:
    """ディレクトリを作成し、プローブファイルの書き込み/削除で書き込み可否を判定する。"""
    try:
        directory.mkdir(parents=True, exist_ok=True)
        probe = directory / ".write_probe"
        probe.write_text("", encoding="utf-8")
        probe.unlink()
        return True
    except OSError:
        return False


def _first_writable_dir(candidates: list[Path]) -> Path:
    """候補を順に試し、最初に書き込み可能なディレクトリを返す。無ければ OSError。"""
    for cand in candidates:
        if _dir_is_writable(cand):
            return cand
    raise OSError(
        "レポートの書き込み先がありません: "
        + ", ".join(str(c) for c in candidates)
        + f"（{user_config.CONFIG_FILENAME} の [paths] reports_dir に書き込み可能なフォルダを"
        "書いてください）"
    )


def _build_report_meta(analysis_id: str, dataset: str, status: str) -> dict:
    """レポートの frontmatter メタを組み立てる。"""
    return {
        "type": "report",
        "analysis_id": analysis_id,
        "dataset": dataset,
        "date": _date.today().isoformat(),
        "status": status,
    }


#: サーバの指示（クライアントがシステムプロンプトに載せる）。Claude Code は 4,096 字で切るので、
#: 4,000 字以内に保つ（tests/test_instructions_text.py が縛る）。詳細は各ツールの説明と output-format へ委ねる。
MCP_INSTRUCTIONS = """
Parses and analyzes MS-DIAL output (lipidomics and general metabolomics).

MEANINGS — `lipidmix://docs/output-format` (core) and `.../output-format/{topic}` define every field; topics: `dataset` `dcl` `eic` `identity` `mztab` `library` `curation` `pathway` `reanalysis` `quant`. `docs_read(topic=..., section=...)` returns the same text. Do not infer a field from its name. A row is a feature (alignment spot) or one sample's value of it; gap-filled values are not detected peaks; EIC `peak_top` is a coordinate, not an intensity.

ENTRY POINT — two equal entries meet at an MS-DIAL analysis folder. Start with `list_data_files(directory)` (`all_files=True` adds raw files).
1. MS-DIAL output present (.arf/.arf2/.dcl, mzTab-M, a completed job or pipeline run), even beside raw files: `load_dataset(directory=...)` (or `mztab_path=` / `job_path=` / `pipeline_path=`). Each load opens a handle that becomes current; mixed processing batches are fine (the latest is used). On `status: "loading"`, poll `dataset_load_status`. `dataset_overview` and `dataset_status` (groups, roles) show what was loaded.
2. Only raw files (.d and Waters .raw are folders): `pipeline_run` runs MS-DIAL Console (tens of minutes per polarity; writes a run folder); poll `pipeline_status`; at `completed`/`partial`, `load_dataset(pipeline_path=..., polarity=...)`. On `needs_input`, fix the sample sheet via `pipeline_resume`; `rerun_upstream=true` re-runs Console. To reprocess with other Console settings (peak height, filters), call `pipeline_run` again with `request.method_overrides` (a new run). Helpers: `pipeline_plan`, `pipeline_sample_sheet` (from ISA/mwTab), `console_method_candidates` / `console_method_template`.

ASSAY KIND — the format does not say lipid or metabolite. AFTER `load_dataset`, decide it from the evidence (annotations, m/z range, what you were given) and store it with `record_objective(assay_kind=...)` (kept per handle). Lipid-name grammar applies only to `lipid`. Send no filename-derived term to external services.

FIRST LOOK — `dataset_overview` returns facts with IDs, quality (`poor`) before differences; its group direction is provisional until `dataset_compare(reference_group=..., test_group=...)`.

RECORDS — `record_objective`, `write_report`, `save_figure` and `pathway_*` write under `[paths] reports_dir` if set, else `reports/` in the loaded dataset's analysis folder; with no dataset they return `missing_state`. Statistics tools log each call under `reports/stats/` (`record`); cite it in `write_report`.

MS/MS EVIDENCE — `has_msms` only says an MS/MS acquisition exists; spectra are in `.dcl` (`dcl_*`). Score them with `library_match_feature` / `library_msms_qualify`; `dataset_claim_evidence` bundles a claim's numbers with a `rule_hash`. MSI level 2 needs a matched spectrum. `not_found` from `dcl_find_msms` means no MS/MS was acquired, not that fragments are absent.

STATISTICS — `dataset_preprocess` returns a `matrix_id`; `dataset_compare` / `dataset_pca` use the latest one unless given (`matrix_source` says which). Verdicts changed / unchanged / inconclusive rest on q and the CI: "not significant" is not "unchanged". Results carry their warnings (confounding, small n, normalization, drift not corrected). A `needs_input` reply names an ambiguous choice: answer it through a tool in `retry_tools`. A long call may return a `job_id` (`job_status`; meanwhile its handle answers `HANDLE_BUSY`); a long table returns its first rows with `truncated: true` (pass an output path for all).

PARALLEL REANALYSIS — calls without `dataset_handle` act on the current dataset, which every `load_dataset` switches. When agents share this server, pass `dataset_handle`; `dataset_list` / `dataset_close` manage handles. Spec-driven tools (`dataset_analysis_run`, `dataset_claim_evidence`, `library_msms_qualify`) return numbers and a `rule_hash`, never verdicts. `dataset_alignment_qc` finds one compound split over two spots.
""".strip()

mcp = FastMCP(
    "ms-data-parser",
    instructions=MCP_INSTRUCTIONS,
    # HTTPトランスポート時のみ使用。stdioでは無視される。
    host=os.environ.get("LIPIDMIX_HOST", "127.0.0.1"),
    port=int(os.environ.get("LIPIDMIX_PORT", "8000")),
)

def _report_dir_candidates(holder=None) -> list[Path]:
    """記録の書き込み先の候補（先に試すものから）。

    設定ファイルの `[paths] reports_dir`（あれば先。環境変数では指さない。spec 2026-11-06
    §7.4）→ `holder`（`.dataset` を持つ、ハンドルか対話のセッション）のデータセットの
    `analysis_dir` 配下 `reports/`。どちらも無い（設定が無く、データセットも無い）ときは
    `NO_REPORT_DIR` の DomainError（ツールの境界は `mcp_errors.config_invalid_as_envelope` が
    `missing_state("dataset", ["load_dataset"])` にする）。設定ファイルが読めなければ
    `CONFIG_INVALID` の DomainError（黙って解析フォルダへ書き先を変えない）。

    解析の記録（objective とレポート）は生データ（解析フォルダ）の側に置く契約。
    リポジトリ側（`<project>/reports` や `<project>/analyses`）へ退避する経路は
    持たない——コードの置き場に解析記録が溜まるバグの温床だった。`holder` は呼び出し側が
    渡す（このモジュールは leaf なので、ハンドルの登録簿も session も import しない）。
    """
    try:
        setting = user_config.get_setting("paths.reports_dir")
    except user_config.ConfigInvalidError as exc:
        raise DomainError(exc.code, exc.message, exc.details()) from exc
    candidates: list[Path] = []
    if setting is not None:
        candidates.append(Path(setting.value))
    analysis_dir = getattr(getattr(holder, "dataset", None), "analysis_dir", None)
    if analysis_dir:
        candidates.append(Path(analysis_dir) / "reports")
    if not candidates:
        raise DomainError(
            NO_REPORT_DIR,
            "記録の置き場所が決まりません（設定ファイルの [paths] reports_dir が無く、"
            "データセットも開いていません）。先に load_dataset を実行してください。")
    return candidates


#: 記録の置き場所が決まらないときの `DomainError` のコード（`missing_state` の封筒になる）。
NO_REPORT_DIR = "NO_REPORT_DIR"


def _resolve_report_dir(holder=None) -> Path:
    """候補順に書き込み可能なレポートディレクトリを返す。"""
    return _first_writable_dir(_report_dir_candidates(holder))
