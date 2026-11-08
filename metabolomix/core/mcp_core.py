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


MCP_INSTRUCTIONS = """
This server parses and analyzes MS-DIAL outputs — lipidomics and general
(hydrophilic) metabolomics both.

Before interpreting any output from `load_dataset` / `dataset_*` / `curation_*` /
`library_*` / `dcl_*`, you MUST read the MCP resource `lipidmix://docs/output-format`
(the shared core: row granularity, lipid-name grammar, mandatory caveats) and use
it as the authoritative definition. Do not infer a field's meaning from its name
alone. In particular, distinguish a feature (alignment spot) from one sample's
value of it, gap-filled values from detected peaks, and EIC `peak_top`
coordinates from intensity.

The per-topic field definitions are split into topic resources — read the one
matching the output you are about to interpret, not the whole document:
`lipidmix://docs/output-format/{topic}` where topic is one of `dataset` (a
dataset loaded from a MS-DIAL analysis folder: rows are feature x sample,
gap-fill and `detected_mask`, the `.arf2` representative and its identity
columns, sample roles, preprocessing/QC and differential), `dcl` (.dcl
deconvoluted MS/MS), `eic` (EIC figures from `dataset_plot`), `identity`
(annotation confidence), `mztab` (a dataset loaded from mzTab-M: SME vs SML
identification, the feature table, the differential export's identity columns),
`library` (reference-library MS/MS matching: what each score means, `-1` vs `0`,
the deliberately ported upstream quirks), `curation` (alignment curation:
`curation_review`'s verdict, reason-code table, and mandatory caveats — e.g.
`UNKNOWN` is not a mismatch), `pathway` (MSLipidMapper pathway projection: node
table and adjacent opposite pairs — adjacency is a drawing template, not
curated reactions), `reanalysis` (dataset handles for parallel reanalysis
agents, the spec-driven matrix/statistics tools, MS/MS pass/fail rules and the
claim evidence bundle), or `quant` (internal-standard normalization, polarity
merge and feature selection: `dataset_quant_*`, `dataset_merge_polarities`,
`dataset_select_features` - the gap-fill mark in `detected_mask`, units,
`needs_input` envelopes, records). `load_dataset` carries a short semantics digest in its
`caveat` field until you have fetched the core; fetch the topic before
interpreting the matching tool's output.

ENTRY POINT — first check whether MS-DIAL has already been run on the user's
folder, then go downstream if it has, or ask before running it if it has not.
Look before picking a tool: `list_data_files(directory)` (the directory is
required; nothing is searched on its own) lists MS-DIAL output in the folder, and
`all_files=True` also shows raw instrument files.

1. MS-DIAL HAS BEEN RUN (the folder contains MS-DIAL output: .arf/.arf2/
   .dcl/.EIC.aef, an mzTab-M file, a completed analysis-job.json, or a
   completed pipeline run). Go downstream — this is the main path, and it
   applies even when raw files sit in the same folder. Call
   `load_dataset(directory=...)` (or `mztab_path=` / `job_path=` /
   `pipeline_path=`). It opens the dataset as a new handle that EVERY time
   becomes the current one (auto-selecting PeakProperties.arf over
   DriftSpots.arf) and returns the counts, the group and role breakdown and the
   library record — no overview or PCA is run for you. Tools called without
   `dataset_handle` act on that current dataset and use its files; a tool that
   needs a file the dataset does not have says so instead of searching a data
   folder. Look at `dataset_status` (`view=summary|samples|tags`) for the groups,
   roles and tags before proposing a comparison. Re-run MS-DIAL on such a folder
   only if the user explicitly asks for it.
2. NO TRACE OF MS-DIAL (only unprocessed instrument files). Any format
   MS-DIAL itself reads qualifies: .abf .ibf .cdf .mzml .wiff .raw .d .wiff2
   .qgd .lcd .lrp .imzml. For Agilent/Bruker .d and Waters .raw, one
   measurement is a FOLDER, not a file — a directory listing showing only
   sub-folders with those extensions IS raw data, not an empty folder. Running
   MS-DIAL Console is an auxiliary path that takes tens of minutes and writes
   outputs next to the raw data, so ASK THE USER FOR PERMISSION before
   launching it, and launch nothing until they agree. Once they agree, call
   `pipeline_run(dataset_root)`. If they want to see the resolved settings
   (method file, LBM, polarity) first, call `pipeline_plan` before asking.
   `pipeline_run` starts an independent worker that carries the run through
   upstream execution, metadata resolution, preprocessing, PCA, and (if
   comparisons are given or confirmed) differential analysis and reporting,
   returning a compact dispatch receipt (`pipeline_id`/`pipeline_path`) — not
   the results themselves. Poll with `pipeline_status(pipeline_path)` and,
   once it reports `needs_input`, resolve ONLY the missing items with
   `pipeline_resume` (never re-ask about items already confirmed).

MIXED-DATE FOLDERS ARE FINE — a folder may contain files from several MS-DIAL
processing runs (multiple dates/batches). This is NOT a blocker and must not be
treated as unanalyzable: the file selection auto-selects the LATEST batch (by
the `AlignmentResult_<timestamp>` embedded in the filenames). `load_dataset`
reports which batch it selected and which `.arf` it did not choose. Proceed with
analysis; only ask the user if they explicitly want an older batch.

ASSAY KIND — settle whether this is `lipid` (lipidomics) or `metabolite`
(general/hydrophilic metabolomics) before interpreting any compound name. THE
FILE FORMAT DOES NOT TELL YOU: MS-DIAL writes both into the same .arf/.mzTab.
Until it is settled, apply neither lipid shorthand grammar nor lipid-class
knowledge. Settle it AFTER `load_dataset`, because a record needs a loaded
dataset to know where it goes: from the groups, polarity and compound classes the
load shows, infer the likely objective and the biological context (organism/cell
line/treatment) and present 1-2 candidates WITH the evidence behind each — never
a single confident statement. Only the user's confirmation settles it. Then
record it with `record_objective(assay_kind=...)` (create it with analysis_id,
dataset, polarity, groups, comparison and sub_questions; later calls change only
the fields they pass, e.g. `add_subquestions`). That call also switches the
semantics digest this server sends: the rules for the settled kind arrive once
in its reply. For `metabolite`, GOSLIN/LIPID MAPS normalization and the lipid
checks do not apply; identification rests on the candidate set (adducts,
isomers, rank). Never send filename-derived terms to external services before
the user confirms them.

RECORDS — `record_objective`, `write_report`, `save_figure`, `pathway_project` and
`pathway_plot` write, and `read_report` reads, under `[paths] reports_dir` from the
settings file when it is set, otherwise under `reports/` in the analysis folder of the
loaded dataset (never in the repository), so a dataset must be loaded first; without
one they answer `missing_state`.

MS/MS EVIDENCE — the real spectra live in `.dcl`; the acquisition flags that
loaded rows carry (has_msms) only record that an MS/MS acquisition reference
exists. Read the spectrum itself with `dcl_spot_spectrum`, `dcl_diagnostic_ions`
or `dcl_find_msms` (`dcl_*`; they read the dataset's alignment `.dcl`, and a
per-sample `.dcl` needs `file_path`), score it against a reference library with
`library_msms_qualify` / `library_match_feature`, and take the pass/fail numbers
for a claim from `dataset_claim_evidence` (they return numbers and a `rule_hash`,
never verdicts). When you claim identification confidence (MSI Level 2), rest it
on a spectrum seen in the `.dcl` — never on the acquisition flag alone. A
`not_found` from `dcl_find_msms` means no MS/MS was acquired for that precursor,
NOT that the expected fragments are absent.

DIFFERENTIAL ANALYSIS — `dataset_preprocess` (optional; normalization / QC
filtering / drift correction, and it returns a `matrix_id`) builds the matrix that
`dataset_compare` and `dataset_pca` read; without a `matrix_id` they use the
default preprocessing, so fold changes are not dominated by per-sample loading
differences. `dataset_exclude` removes samples or features (reversibly) before it.
Read each feature's verdict — changed / unchanged / inconclusive — which rests on
q and the confidence interval, never on |log2fc| alone. When a tool answers
`needs_input` (not an error), ask the user, then call a tool in `retry_tools`
with the answer in `argument` (for `dataset_analysis_run` the answer goes into
the spec `field`). Always surface the tool's caveats — group/batch confounding,
small n (few replicates), and normalization status — as first-class findings,
never bury them.

PARALLEL REANALYSIS — when several agents share this server (a reanalysis team),
do not rely on the current dataset: `load_dataset` always makes the handle it
opens the current one, which every call without `dataset_handle` uses, and
`library_load` without `dataset_handle` loads into that current dataset, so the
agents would switch it under each other. Pass the `dataset_handle` that
`load_dataset` returned to `dataset_*`, `library_*`, `curation_*` and the
reanalysis tools (`dataset_list` / `dataset_close` manage the open handles);
calls on different handles run in parallel and never touch each other. Decision
numbers come from spec-driven tools (`dataset_analysis_run`,
`dataset_claim_evidence`, `library_msms_qualify`), which return a `rule_hash` and
give identical JSON for an identical spec — see
`lipidmix://docs/output-format/reanalysis`. They return numbers, never verdicts.
`dataset_alignment_qc` finds one compound split over two alignment spots (e.g. the
first injections drifted in RT); run it before reading results.
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
