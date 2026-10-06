# tests/test_console_completion.py
"""Console完了ゲート（metabolomix.console.validation）の検証。

exit_code=0・termination=exited だけでは「完了」ではない。主mzTab-Mが構造的に
妥当で、定量行列に有限値があり、予定した準備済みraw全件がassayへ1対1で
対応していることまで確認して初めてcompletedとする（spec §5.2）。
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from tests.pipeline_fixtures import execution_record, write_mztab


# ---------- completion_status ----------

def test_intermediate_output_is_not_completed():
    from metabolomix.console.validation import completion_status
    receipt = execution_record(exit_code=1)
    validation = {"ok": False, "errors": ["PRIMARY_MZTAB_MISSING"]}
    assert completion_status(receipt, validation, has_artifacts=True) == "partial"


def test_zero_exit_with_invalid_mztab_is_not_completed():
    from metabolomix.console.validation import completion_status
    assert completion_status(execution_record(), {"ok": False}, True) != "completed"


def test_zero_exit_with_valid_mztab_is_completed():
    from metabolomix.console.validation import completion_status
    validation = {"ok": True, "errors": [], "warnings": []}
    assert completion_status(execution_record(), validation, True) == "completed"


def test_failed_execution_without_artifacts_is_failed():
    from metabolomix.console.validation import completion_status
    receipt = execution_record(exit_code=1)
    validation = {"ok": False, "errors": ["PRIMARY_MZTAB_MISSING"]}
    assert completion_status(receipt, validation, has_artifacts=False) == "failed"


def test_non_exited_termination_with_valid_validation_is_not_completed():
    """terminationがexitedでなければ、validationがokでもcompletedにしない。"""
    from metabolomix.console.validation import completion_status
    receipt = execution_record(exit_code=None, termination="timeout")
    validation = {"ok": True, "errors": [], "warnings": []}
    assert completion_status(receipt, validation, has_artifacts=True) == "partial"


def test_bool_exit_code_is_not_treated_as_int():
    """execution.validate_execution_recordと同じ厳密さ: type(rc) is int のみ許容する。"""
    from metabolomix.console.validation import completion_status
    receipt = execution_record(exit_code=True)  # bool は int のサブクラスだが弾く
    validation = {"ok": True, "errors": [], "warnings": []}
    assert completion_status(receipt, validation, has_artifacts=True) == "partial"


# ---------- write_mztab fixtureが実物のreader/validatorに通ることの検証 ----------

def test_write_mztab_passes_real_parser_and_validator(tmp_path):
    from metabolomix.mztab.reader import parse_mztab, extract_abundance_matrix
    from metabolomix.mztab.validator import validate_mztab

    raws = [tmp_path / "raw" / "S1.raw", tmp_path / "raw" / "S2.raw"]
    for r in raws:
        r.parent.mkdir(parents=True, exist_ok=True)
        r.write_bytes(b"\x00")

    mztab_path = write_mztab(tmp_path / "AlignResult-1.mzTab", raws)
    parsed = parse_mztab(mztab_path)
    result = validate_mztab(parsed)

    assert result["ok"] is True, result["errors"]
    assert result["errors"] == []

    matrix, sample_names, feature_ids = extract_abundance_matrix(parsed)
    assert len(feature_ids) >= 3
    assert len(sample_names) == 2
    assert matrix.size > 0
    import numpy as np
    assert np.isfinite(matrix).any()


def test_write_mztab_without_inchikey_writes_null_identifier(tmp_path):
    from metabolomix.mztab.reader import parse_mztab
    raws = [tmp_path / "S1.raw"]
    raws[0].write_bytes(b"\x00")
    mztab_path = write_mztab(tmp_path / "a.mzTab", raws, with_inchikey=False)
    parsed = parse_mztab(mztab_path)
    rows = parsed["sections"]["SMF"]["rows"]
    assert all(r.get("database_identifier") is None for r in rows)


# ---------- map_assays ----------

def _job(tmp_path, *, polarity="negative", measure="peak_height", mztab_rel="AlignResult-1.mzTab"):
    from metabolomix.handoff.schema import AnalysisJob, MztabEntry, SCHEMA_VERSION
    run_dir = tmp_path / "runs" / "job_test"
    run_dir.mkdir(parents=True, exist_ok=True)
    return AnalysisJob(
        schema=SCHEMA_VERSION, job_id="job_test", status="running",
        created_at="2026-09-05T00:00:00+09:00", updated_at="2026-09-05T00:00:00+09:00",
        dataset_root=str(tmp_path), input_count=2,
        software_name="MS-DIAL", software_version="5.5", execution_mode="console",
        method_file="m.txt", omics="lipidomics", polarity=polarity, measure=measure,
        run_dir=str(run_dir),
        primary_mztab_files=[MztabEntry(path=mztab_rel, polarity=polarity, measure=measure,
                                        sha256="", validation={}, root="run_dir")],
    )


def test_map_assays_matches_via_decoded_uri_not_display_name(tmp_path):
    """assay[N]の表示名（S1等）が何であってもURIの実パスだけで結合する。"""
    from metabolomix.mztab.reader import parse_mztab
    from metabolomix.console.validation import map_assays

    raws = [tmp_path / "raw" / "Sample A.raw", tmp_path / "raw" / "Sample B.raw"]
    for r in raws:
        r.parent.mkdir(parents=True, exist_ok=True)
        r.write_bytes(b"\x00")
    mztab_path = write_mztab(tmp_path / "a.mzTab", raws)
    parsed = parse_mztab(mztab_path)

    staged = [str(r) for r in raws]
    result = map_assays(parsed, staged)
    assert result == {
        "abundance_assay[1]": staged[0],
        "abundance_assay[2]": staged[1],
    }


def test_map_assays_decodes_url_characters_and_unicode_whitespace(tmp_path):
    """空白・日本語を含むファイル名でも file URI のURLエンコードを正しくdecodeする。"""
    from metabolomix.mztab.reader import parse_mztab
    from metabolomix.console.validation import map_assays

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    raw = raw_dir / "日本語 サンプル 1.raw"
    raw.write_bytes(b"\x00")
    mztab_path = write_mztab(tmp_path / "a.mzTab", [raw])
    parsed = parse_mztab(mztab_path)

    staged = [str(raw)]
    result = map_assays(parsed, staged)
    assert result == {"abundance_assay[1]": staged[0]}


def test_map_assays_normalizes_windows_case_for_matching(tmp_path):
    """Windowsのパス照合は大文字小文字を正規化する。"""
    from metabolomix.mztab.reader import parse_mztab
    from metabolomix.console.validation import map_assays

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    mztab_path = write_mztab(tmp_path / "a.mzTab", [raw])
    parsed = parse_mztab(mztab_path)

    staged = [str(raw).upper()]
    result = map_assays(parsed, staged)
    assert result == {"abundance_assay[1]": staged[0]}


def test_map_assays_leaves_unmatched_location_as_is_when_not_staged(tmp_path):
    """staged_sourcesに無いassayは、decode済みの実パスのまま残す（黙って捨てない）。"""
    from metabolomix.mztab.reader import parse_mztab
    from metabolomix.console.validation import map_assays

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    mztab_path = write_mztab(tmp_path / "a.mzTab", [raw])
    parsed = parse_mztab(mztab_path)

    result = map_assays(parsed, staged_sources=[])
    assert result == {"abundance_assay[1]": str(raw.resolve())}


# ---------- validate_outputs ----------

def test_validate_outputs_ok_when_everything_matches(tmp_path):
    from metabolomix.console.validation import validate_outputs

    raws = [tmp_path / "raw" / "S1.raw", tmp_path / "raw" / "S2.raw"]
    for r in raws:
        r.parent.mkdir(parents=True, exist_ok=True)
        r.write_bytes(b"\x00")
    job = _job(tmp_path)
    write_mztab(Path(job.run_dir) / "AlignResult-1.mzTab", raws)

    staged = [str(r) for r in raws]
    result = validate_outputs(job, execution_record(), staged)

    assert result["ok"] is True, result["errors"]
    assert result["errors"] == []
    assert result["primary_path"] == str((Path(job.run_dir) / "AlignResult-1.mzTab").resolve())
    assert result["sample_map"] == {
        "abundance_assay[1]": staged[0],
        "abundance_assay[2]": staged[1],
    }


def test_validate_outputs_fails_on_missing_prepared_raw(tmp_path):
    """準備済みrawの一部がassayから欠落しているとok=false。"""
    from metabolomix.console.validation import validate_outputs

    raws = [tmp_path / "raw" / "S1.raw"]
    raws[0].parent.mkdir(parents=True, exist_ok=True)
    raws[0].write_bytes(b"\x00")
    job = _job(tmp_path)
    write_mztab(Path(job.run_dir) / "AlignResult-1.mzTab", raws)

    staged = [str(raws[0]), str(tmp_path / "raw" / "S2_missing.raw")]
    result = validate_outputs(job, execution_record(), staged)

    assert result["ok"] is False
    assert "SAMPLE_MAPPING_MISSING" in result["errors"]


def test_validate_outputs_fails_on_unexpected_extra_sample(tmp_path):
    """mzTabのassayがstaged_sourcesに無いrawを指しているとok=false。"""
    from metabolomix.console.validation import validate_outputs

    raws = [tmp_path / "raw" / "S1.raw", tmp_path / "raw" / "Extra.raw"]
    for r in raws:
        r.parent.mkdir(parents=True, exist_ok=True)
        r.write_bytes(b"\x00")
    job = _job(tmp_path)
    write_mztab(Path(job.run_dir) / "AlignResult-1.mzTab", raws)

    staged = [str(raws[0])]  # Extra.raw は準備済み一覧にない
    result = validate_outputs(job, execution_record(), staged)

    assert result["ok"] is False
    assert "SAMPLE_MAPPING_EXTRA" in result["errors"]


def test_validate_outputs_fails_when_one_raw_has_multiple_assays(tmp_path):
    """1rawから複数assayが出るのは初期版の契約外として停止する。"""
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    # 2つのassayが同じrawを指すmzTabを手で組む（write_mztabは1raw=1assay前提のため）
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uri = raw.resolve().as_uri()
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tmzTab-mode\tComplete
        MTD\tmzTab-type\tQuantification
        MTD\tms_run[1]-location\t{uri}
        MTD\tms_run[2]-location\t{uri}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        MTD\tassay[2]-ms_run_ref\tms_run[2]
        SFH\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\tabundance_assay[1]\tabundance_assay[2]
        SMF\t1\tSML:1\tIPCSVZSSVZVIGE-UHFFFAOYSA-N\tPC 36:2\tnull\tnull\t100.0\t110.0
        SMF\t2\tSML:2\tXKMRRTOUMJRJIA-UHFFFAOYSA-N\tPE 36:2\tnull\tnull\t200.0\t210.0
        SMF\t3\tSML:3\tDGGXCMYPQAOAJC-UHFFFAOYSA-N\tTG 52:3\tnull\tnull\t300.0\t310.0
    """)
    mztab_path.write_text(content, encoding="utf-8")

    staged = [str(raw)]
    result = validate_outputs(job, execution_record(), staged)

    assert result["ok"] is False
    assert "SAMPLE_MAPPING_DUPLICATE" in result["errors"]


def test_validate_outputs_fails_on_ambiguous_primary_candidates(tmp_path):
    """候補が複数あればok=false（正準選択の一意性）。"""
    from metabolomix.handoff.schema import MztabEntry
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    write_mztab(Path(job.run_dir) / "a.mzTab", [raw])
    write_mztab(Path(job.run_dir) / "b.mzTab", [raw])
    job.primary_mztab_files = [
        MztabEntry(path="a.mzTab", polarity="negative", measure="peak_height", sha256=""),
        MztabEntry(path="b.mzTab", polarity="negative", measure="peak_height", sha256=""),
    ]

    result = validate_outputs(job, execution_record(), [str(raw)])
    assert result["ok"] is False
    assert "AMBIGUOUS_PRIMARY_MZTAB" in result["errors"]


def test_validate_outputs_fails_on_empty_numeric_matrix(tmp_path):
    """特徴量0件（数値行列が空）はok=false。"""
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uri = raw.resolve().as_uri()
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tms_run[1]-location\t{uri}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        SFH\tSMF_ID\tabundance_assay[1]
    """)
    mztab_path.write_text(content, encoding="utf-8")

    result = validate_outputs(job, execution_record(), [str(raw)])
    assert result["ok"] is False
    assert "EMPTY_ABUNDANCE_MATRIX" in result["errors"]


def test_validate_outputs_structure_errors_are_stable_codes_not_raw_prose(tmp_path):
    """validate_mztabの日本語エラー文をerrorsへ生で混ぜない。

    Task 4のsupervise()やTask 7のload_dataset_stateはerrorsを文字列比較で分岐する
    契約なので、errorsはSCREAMING_SNAKE_CASEの安定コードのみを持つ必要がある
    （PRIMARY_MZTAB_MISSING等と同じ形）。validate_mztabの人間可読な日本語文は
    捨てずに`structure_errors`で読める。
    """
    import re
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uri = raw.resolve().as_uri()
    # SMF行が1件も無い＝validate_mztabが構造エラー（日本語の人間可読文）を返す
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tms_run[1]-location\t{uri}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        SFH\tSMF_ID\tabundance_assay[1]
    """)
    mztab_path.write_text(content, encoding="utf-8")

    result = validate_outputs(job, execution_record(), [str(raw)])

    assert result["ok"] is False
    assert "MZTAB_STRUCTURE_INVALID" in result["errors"]
    for code in result["errors"]:
        assert re.match(r"^[A-Z][A-Z0-9_]*$", code), f"raw prose leaked into errors: {code!r}"
    # 日本語の人間可読文は破棄されず structure_errors から読める
    assert any("SMF" in e for e in result["structure_errors"])


def test_validate_outputs_fails_on_all_nan_matrix(tmp_path):
    """特徴量はあるが全欠損（有限値なし）はok=false。"""
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uri = raw.resolve().as_uri()
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tms_run[1]-location\t{uri}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        SFH\tSMF_ID\tabundance_assay[1]
        SMF\t1\tnull
        SMF\t2\tnull
        SMF\t3\tnull
    """)
    mztab_path.write_text(content, encoding="utf-8")

    result = validate_outputs(job, execution_record(), [str(raw)])
    assert result["ok"] is False
    assert "NO_FINITE_ABUNDANCE_VALUES" in result["errors"]


def test_validate_outputs_fails_on_measure_mismatch(tmp_path):
    """ファイル名がAreaなのにjobがpeak_heightを宣言しているのは致命的不一致。"""
    from metabolomix.console.validation import validate_outputs

    raws = [tmp_path / "raw" / "S1.raw"]
    raws[0].parent.mkdir(parents=True, exist_ok=True)
    raws[0].write_bytes(b"\x00")
    job = _job(tmp_path, measure="peak_height", mztab_rel="Area_AlignResult-1.mzTab")
    write_mztab(Path(job.run_dir) / "Area_AlignResult-1.mzTab", raws)

    result = validate_outputs(job, execution_record(), [str(raws[0])])
    assert result["ok"] is False
    assert "MEASURE_MISMATCH" in result["errors"]


def test_validate_outputs_fails_on_polarity_majority_mismatch(tmp_path):
    """宣言はnegativeだがアダクト多数決はpositive、は致命的不一致。"""
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path, polarity="negative")
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uri = raw.resolve().as_uri()
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tms_run[1]-location\t{uri}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        SFH\tSMF_ID\tabundance_assay[1]
        SMF\t1\t100.0
        SMF\t2\t200.0
        SMF\t3\t300.0
        SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tadduct_ions
        SML\t1\t1\tnull\t[M+H]1+
        SML\t2\t2\tnull\t[M+NH4]1+
        SML\t3\t3\tnull\t[M+Na]1+
    """)
    mztab_path.write_text(content, encoding="utf-8")

    result = validate_outputs(job, execution_record(), [str(raw)])
    assert result["ok"] is False
    assert "POLARITY_MISMATCH" in result["errors"]


def test_validate_outputs_fails_when_primary_missing(tmp_path):
    from metabolomix.console.validation import validate_outputs
    job = _job(tmp_path)
    job.primary_mztab_files = []
    result = validate_outputs(job, execution_record(), [])
    assert result["ok"] is False
    assert "PRIMARY_MZTAB_MISSING" in result["errors"]
    assert result["primary_path"] is None


# ---------- 制御ファイルの除外（回帰） ----------

def test_control_files_are_not_counted_as_artifacts(tmp_path):
    from metabolomix.console.output_collector import collect_artifacts, snapshot

    run_dir = tmp_path / "runs" / "job_test"
    run_dir.mkdir(parents=True)
    befores = {"run_dir": snapshot(run_dir)}

    (run_dir / "execution-result.json").write_text("{}", encoding="utf-8")
    (run_dir / "worker.json").write_text("{}", encoding="utf-8")
    (run_dir / "AlignResult-1.mzTab").write_text(
        "MTD\tmzTab-version\t2.0.0-M\n", encoding="utf-8")

    mztab_entries, other_artifacts = collect_artifacts({"run_dir": run_dir}, befores)

    all_paths = {a.path for a in other_artifacts} | {"AlignResult-1.mzTab" if mztab_entries else ""}
    assert "execution-result.json" not in all_paths
    assert "worker.json" not in all_paths
    assert len(mztab_entries) == 1


def test_reserved_names_are_limited_to_files_something_actually_writes(tmp_path):
    """除外集合に「誰も書かない予約名」を残さない（最終レビュー保留2）。

    `control.json` はどのコードも書かない予約名だった——将来 MS-DIAL が
    その名前で出力を書けば、収集が黙って捨てて「生成物ゼロ」に化ける。
    実在する運用ファイル名（`metabolomix.console.execution` の定数）だけを
    除外集合に置く。
    """
    from metabolomix.console import execution
    from metabolomix.console.output_collector import _OPERATIONAL_FILES, collect_artifacts, snapshot

    assert "control.json" not in _OPERATIONAL_FILES
    assert execution.RECEIPT_FILENAME in _OPERATIONAL_FILES
    assert execution.SUPERVISION_FILENAME in _OPERATIONAL_FILES
    assert execution.CONSOLE_LOG_FILENAME in _OPERATIONAL_FILES

    # 予約を外した名前は、実際に成果物として拾われる。
    run_dir = tmp_path / "runs" / "job_reserved"
    run_dir.mkdir(parents=True)
    befores = {"run_dir": snapshot(run_dir)}
    (run_dir / "control.json").write_text("{}", encoding="utf-8")

    _mztab_entries, other_artifacts = collect_artifacts({"run_dir": run_dir}, befores)

    assert "control.json" in {a.path for a in other_artifacts}


# ---------- 定量列とassay対応の突き合わせ（レビュー指摘5） ----------

def test_validate_outputs_fails_when_a_declared_assay_has_no_abundance_column(tmp_path):
    """MTDがassay[2]を宣言しているのにSMFへabundance_assay[2]列が無い。

    その検体の定量値はどこにも無いのに、残りの列だけで行列は非空・有限になる
    ——MTDだけで対応表を作ると、定量列が欠けたmzTabがcompletedとして通る。
    """
    from metabolomix.console.validation import validate_outputs

    raws = [tmp_path / "raw" / "S1.raw", tmp_path / "raw" / "S2.raw"]
    for r in raws:
        r.parent.mkdir(parents=True, exist_ok=True)
        r.write_bytes(b"\x00")
    job = _job(tmp_path)
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uris = [r.resolve().as_uri() for r in raws]
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tmzTab-mode\tComplete
        MTD\tmzTab-type\tQuantification
        MTD\tms_run[1]-location\t{uris[0]}
        MTD\tms_run[2]-location\t{uris[1]}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        MTD\tassay[2]-ms_run_ref\tms_run[2]
        SFH\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\tabundance_assay[1]
        SMF\t1\tSML:1\tIPCSVZSSVZVIGE-UHFFFAOYSA-N\tPC 36:2\tnull\tnull\t100.0
        SMF\t2\tSML:2\tXKMRRTOUMJRJIA-UHFFFAOYSA-N\tPE 36:2\tnull\tnull\t200.0
        SMF\t3\tSML:3\tDGGXCMYPQAOAJC-UHFFFAOYSA-N\tTG 52:3\tnull\tnull\t300.0
    """)
    mztab_path.write_text(content, encoding="utf-8")

    result = validate_outputs(job, execution_record(), [str(r) for r in raws])

    assert result["ok"] is False
    assert "ABUNDANCE_COLUMN_MISSING" in result["errors"]


def test_validate_outputs_fails_on_an_abundance_column_without_an_assay(tmp_path):
    """どのassay対応にもぶら下がらない定量列は、どのrawの値か言えない。"""
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    mztab_path = Path(job.run_dir) / "AlignResult-1.mzTab"
    uri = raw.resolve().as_uri()
    content = textwrap.dedent(f"""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tmzTab-mode\tComplete
        MTD\tmzTab-type\tQuantification
        MTD\tms_run[1]-location\t{uri}
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        SFH\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\tabundance_assay[1]\tabundance_assay[2]
        SMF\t1\tSML:1\tIPCSVZSSVZVIGE-UHFFFAOYSA-N\tPC 36:2\tnull\tnull\t100.0\t110.0
        SMF\t2\tSML:2\tXKMRRTOUMJRJIA-UHFFFAOYSA-N\tPE 36:2\tnull\tnull\t200.0\t210.0
        SMF\t3\tSML:3\tDGGXCMYPQAOAJC-UHFFFAOYSA-N\tTG 52:3\tnull\tnull\t300.0\t310.0
    """)
    mztab_path.write_text(content, encoding="utf-8")

    result = validate_outputs(job, execution_record(), [str(raw)])

    assert result["ok"] is False
    assert "ABUNDANCE_COLUMN_UNMAPPED" in result["errors"]


def test_validate_outputs_accepts_a_one_to_one_assay_and_column_mapping(tmp_path):
    """正常な出力（assayと定量列が1対1）はそのまま通る。"""
    from metabolomix.console.validation import validate_outputs

    raws = [tmp_path / "raw" / "S1.raw", tmp_path / "raw" / "S2.raw"]
    for r in raws:
        r.parent.mkdir(parents=True, exist_ok=True)
        r.write_bytes(b"\x00")
    job = _job(tmp_path)
    write_mztab(Path(job.run_dir) / "AlignResult-1.mzTab", raws)

    result = validate_outputs(job, execution_record(), [str(r) for r in raws])

    assert result["ok"] is True, result["errors"]


# ---------- entryとjobの宣言の食い違いを完了ゲートで捕まえる（レビュー指摘P2） ----------

def test_validate_outputs_fails_when_the_entry_measure_differs_from_the_job(tmp_path):
    """ファイル名側を採用したentryがjobの宣言と食い違うなら、完了にしない。

    下流の`select_primary_entry`は**jobの宣言**で候補を選ぶので、この食い違いを
    通すと「completedなのに読めない」出力ができる。
    """
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path)
    write_mztab(Path(job.run_dir) / "AlignResult-1.mzTab", [raw])
    # collect_artifactsがファイル名から作ったentryが、jobの宣言と食い違った状態。
    job.primary_mztab_files[0].measure = "peak_area_above_zero"
    job.primary_mztab_files[0].validation = {
        "conflicts": {"measure": {"filename": "peak_area_above_zero",
                                  "job_declared": "peak_height"}}}

    result = validate_outputs(job, execution_record(), [str(raw)])

    assert result["ok"] is False
    assert "MEASURE_MISMATCH" in result["errors"]


def test_validate_outputs_fails_when_the_entry_polarity_differs_from_the_job(tmp_path):
    from metabolomix.console.validation import validate_outputs

    raw = tmp_path / "raw" / "S1.raw"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(b"\x00")
    job = _job(tmp_path, polarity="negative")
    write_mztab(Path(job.run_dir) / "AlignResult-1.mzTab", [raw])
    job.primary_mztab_files[0].polarity = "positive"

    result = validate_outputs(job, execution_record(), [str(raw)])

    assert result["ok"] is False
    assert "POLARITY_MISMATCH" in result["errors"]
