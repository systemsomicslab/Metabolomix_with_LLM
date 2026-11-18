"""Console の `-i` に渡す試料の表（`metabolomix.console.input_csv`）と、pipeline がそれを使うこと。

フォルダを `-i` に渡すと Console は各ファイルに別の群を振り、`N percent detected in one group` が効かない。
実験情報シートがあれば群・種類・注入順の CSV を解析フォルダに書いて `-i` に渡し、無ければフォルダのまま。
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from metabolomix.analysis.sample_manifest import parse_manifest
from metabolomix.console import input_csv
from metabolomix.console.execution import _default_command
from metabolomix.console.runner import build_msdial_cmd
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline import recovery, service
from metabolomix.pipeline.store import load_run
from tests.config_fixtures import use_config
from tests.pipeline_fixtures import make_source, write_sample_manifest

_FAKE_CONSOLE = Path(__file__).resolve().parent / "fixtures" / "fake_console.py"

_V2_HEADER = ("sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude\t"
              "biological_sample_id")


def _manifest(path: Path, rows: list[tuple]) -> Path:
    """v2 のシート。行は (source_file, role, group, batch, injection_order, include)。"""
    lines = ["# schema = sample-manifest.v2", _V2_HEADER]
    for index, (source, role, group, batch, order, include) in enumerate(rows):
        lines.append("\t".join([f"id{index}", source, role, group, batch, order, "", include, ""]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _raw_folder(root: Path, names: list[str]) -> list[Path]:
    root.mkdir(parents=True, exist_ok=True)
    for name in names:
        (root / name).write_text("raw", encoding="ascii")
    return [(root / name).resolve() for name in names]


def _console_read(path: Path) -> list[dict]:
    """実 Console（`ReadCsvContents`）と同じ読み方: ASCII・`#` の行を飛ばす・`,` で割るだけ。"""
    lines = [line for line in path.read_bytes().decode("ascii").split("\r\n")
             if line and not line.startswith("#")]
    header = lines[0].split(",")
    return [dict(zip(header, line.split(","))) for line in lines[1:]]


# ---------- 行の組み立て ----------

def test_rows_carry_group_type_order_and_batch_from_the_manifest(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["A1.mzML", "A2.mzML", "B1.mzML", "QC1.mzML",
                                            "blank.mzML", "std.mzML"])
    manifest = _manifest(tmp_path / "sheet.tsv", [
        ("A1.mzML", "sample", "ctrl", "run2", "3", "true"),
        ("A2.mzML", "sample", "ctrl", "run2", "5", "false"),     # include=false でも Console には渡す
        ("B1.mzML", "", "treated", "run1", "1", ""),              # 役割の空は Sample
        ("QC1.mzML", "qc", "", "run1", "2", "true"),             # 群の空は役割でまとめる
        ("blank.mzML", "blank", "", "run2", "4", "true"),
        ("std.mzML", "standard", "", "run2", "6", "true"),
    ])
    rows = parse_manifest(manifest, source_root=tmp_path / "input",
                          expected_sources=[p.name for p in raws])

    out, summary = input_csv.build_rows(raws, rows, acquisition="DDA")

    assert [r["file_path"] for r in out] == [str(p) for p in raws]      # 行の順＝file ID の順
    assert [r["file_name"] for r in out] == ["A1", "A2", "B1", "QC1", "blank", "std"]
    assert [r["class_id"] for r in out] == ["ctrl", "ctrl", "treated", "qc", "blank", "standard"]
    assert [r["file_type"] for r in out] == ["Sample", "Sample", "Sample", "QC", "Blank", "Standard"]
    assert [r["analytical_order"] for r in out] == ["3", "5", "1", "2", "4", "6"]
    # バッチは注入の順に最初に出た順で番号を振る（run1 が先に注入された）
    assert [r["batch_order"] for r in out] == ["2", "2", "1", "1", "2", "2"]
    assert {r["acquisition_type"] for r in out} == {"DDA"}
    assert {r["factor"] for r in out} == {"1"}
    assert summary == {"n_files": 6, "classes": {"ctrl": 2, "treated": 1, "qc": 1, "blank": 1,
                                                 "standard": 1},
                       "file_types": {"Sample": 3, "QC": 1, "Blank": 1, "Standard": 1},
                       "acquisition_type": "DDA", "analytical_order_source": "manifest",
                       "batch_order": {"run1": 1, "run2": 2}}


def test_a_missing_injection_order_falls_back_to_the_file_order_for_every_row(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["S0.mzML", "S1.mzML", "S2.mzML"])
    manifest = _manifest(tmp_path / "sheet.tsv", [
        ("S0.mzML", "sample", "g", "", "7", ""),
        ("S1.mzML", "sample", "g", "", "", ""),
        ("S2.mzML", "sample", "g", "", "9", ""),
    ])
    rows = parse_manifest(manifest, source_root=tmp_path / "input",
                          expected_sources=[p.name for p in raws])

    out, summary = input_csv.build_rows(raws, rows, acquisition="SWATH")

    assert [r["analytical_order"] for r in out] == ["1", "2", "3"]   # 番号の付け方を混ぜない
    assert [r["batch_order"] for r in out] == ["1", "1", "1"]
    assert summary["analytical_order_source"] == "file_order"
    assert "batch_order" not in summary
    assert {r["acquisition_type"] for r in out} == {"SWATH"}


def test_class_ids_that_console_cannot_read_are_escaped_without_merging_groups():
    assert input_csv.class_id("ctrl", "sample") == "ctrl"
    assert input_csv.class_id(None, "qc") == "qc"
    assert input_csv.class_id(None, None) == "sample"
    japanese = input_csv.class_id("対照", "sample")
    comma = input_csv.class_id("a,b", "sample")
    assert japanese.isascii() and "," not in comma and comma.isascii()
    # 素のまま書く名前は `\` を含まず、書き換えた名前は必ず含むので、別の群と同じ名前にならない
    literal = input_csv.class_id("\\u5bfe\\u7167", "sample")
    assert len({japanese, literal, comma, input_csv.class_id("a\\x2cb", "sample")}) == 4


def test_acquisition_type_follows_the_method_like_folder_input():
    assert input_csv.acquisition_type({}) == "DDA"
    assert input_csv.acquisition_type({"acquisition type": " swath "}) == "SWATH"
    assert input_csv.acquisition_type({"acquisition type": "AIF"}) == "AIF"
    assert input_csv.acquisition_type({"acquisition type": "None"}) == "DDA"   # Console が受けない値


def test_paths_console_cannot_read_from_a_csv_are_named():
    assert input_csv.unsafe_names([Path("C:/run/input/S0.mzML")]) == []
    assert input_csv.unsafe_names([Path("C:/run/input/試料.mzML"),
                                   Path("C:/run/input/a,b.mzML")]) == ["試料.mzML", "a,b.mzML"]


def test_written_csv_reads_back_the_way_console_reads_it(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["S0.wiff", "S1.wiff"])
    rows = [{"source_file": str(raws[0]), "role": "sample", "group": "対照", "batch": None,
             "injection_order": 2},
            {"source_file": str(raws[1]), "role": "qc", "group": None, "batch": None,
             "injection_order": 1}]
    out, _summary = input_csv.build_rows(raws, rows, acquisition="DDA")
    path = tmp_path / "input" / input_csv.CSV_FILENAME

    digest = input_csv.write_csv(path, out)

    data = path.read_bytes()
    assert digest == hashlib.sha256(data).hexdigest()
    assert data.startswith(b"file_path,file_name,file_type,class_id,acquisition_type,batch_order,"
                           b"analytical_order,factor\r\n")
    assert data.endswith(b"\r\n") and b"\r\n\r\n" not in data      # 空の行は Console を落とす
    read = _console_read(path)
    assert [r["file_path"] for r in read] == [str(p) for p in raws]
    assert [r["file_type"] for r in read] == ["Sample", "QC"]
    assert read[0]["class_id"] == input_csv.class_id("対照", "sample")
    assert all(value != "" for row in read for value in row.values())   # 空の欄は既定に落ちない


# ---------- コマンドライン ----------

def test_command_passes_the_csv_to_i_and_keeps_the_rest_of_the_layout(tmp_path):
    csv_path = tmp_path / "raw" / input_csv.CSV_FILENAME
    cmd = build_msdial_cmd("fake.exe", tmp_path / "raw", tmp_path / "out", tmp_path / "params.txt",
                           save_project=True, input_csv=csv_path)
    assert cmd == ["fake.exe", "lcms", "-i", str(csv_path), "-o", str(tmp_path / "out"),
                   "-m", str(tmp_path / "params.txt"), "-p"]
    folder = build_msdial_cmd("fake.exe", tmp_path / "raw", tmp_path / "out", tmp_path / "params.txt")
    assert folder[2:4] == ["-i", str(tmp_path / "raw")]


def test_supervision_inputs_choose_the_i_argument(tmp_path):
    job = type("Job", (), {"method_file": str(tmp_path / "m.txt"), "save_project": False})()
    csv_path = tmp_path / "input" / input_csv.CSV_FILENAME
    with_csv = _default_command(job, tmp_path / "input", tmp_path / "out",
                                {"exe_path": "fake.exe",
                                 "console_input": {"form": "csv", "path": str(csv_path)}})
    assert with_csv[2:4] == ["-i", str(csv_path)]
    folder = _default_command(job, tmp_path / "input", tmp_path / "out",
                              {"exe_path": "fake.exe",
                               "console_input": {"form": "folder", "reason": "no_manifest"}})
    assert folder[2:4] == ["-i", str(tmp_path / "input")]
    older = _default_command(job, tmp_path / "input", tmp_path / "out", {"exe_path": "fake.exe"})
    assert older == folder


# ---------- pipeline が `-i` を決める ----------

def test_no_manifest_keeps_the_folder_input(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["S0.mzML"])
    (tmp_path / "raw").mkdir()
    result = service._console_input({}, str(tmp_path / "raw"), tmp_path / "input", raws,
                                    tmp_path / "m.txt")
    assert result == {"form": "folder", "reason": "no_manifest"}
    assert not (tmp_path / "input" / input_csv.CSV_FILENAME).exists()


def test_an_unreadable_manifest_falls_back_to_the_folder_with_a_reason(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["S0.mzML"])
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    (raw_root / "sample-manifest.tsv").write_text("not a manifest\n", encoding="utf-8")
    result = service._console_input({}, str(raw_root), tmp_path / "input", raws, tmp_path / "m.txt")
    assert result == {"form": "folder", "reason": "manifest_unreadable",
                      "error_code": "SAMPLE_MANIFEST_INVALID"}


def test_a_path_the_csv_cannot_carry_falls_back_to_the_folder(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["S0.mzML", "試料.mzML"])
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    _manifest(raw_root / "sample-manifest.tsv", [("S0.mzML", "sample", "a", "", "1", ""),
                                                 ("試料.mzML", "sample", "b", "", "2", "")])
    result = service._console_input({}, str(raw_root), tmp_path / "input", raws, tmp_path / "m.txt")
    assert result == {"form": "folder", "reason": "path_not_csv_safe", "files": ["試料.mzML"]}


def test_manifest_writes_the_csv_into_the_analysis_folder(tmp_path):
    raws = _raw_folder(tmp_path / "input", ["S0.mzML", "S1.mzML"])
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    sheet = _manifest(raw_root / "sheet.tsv", [("S0.mzML", "sample", "a", "", "2", ""),
                                               ("S1.mzML", "qc", "", "", "1", "")])
    method = tmp_path / "m.txt"
    method.write_text("Acquisition type: AIF\n", encoding="ascii")

    result = service._console_input({"sample_manifest": "sheet.tsv"}, str(raw_root),
                                    tmp_path / "input", raws, method)

    csv_path = tmp_path / "input" / input_csv.CSV_FILENAME
    assert result["form"] == "csv"
    assert result["path"] == str(csv_path)
    assert result["sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    assert result["manifest_sha256"] == hashlib.sha256(sheet.read_bytes()).hexdigest()
    assert result["classes"] == {"a": 1, "qc": 1}
    assert result["acquisition_type"] == "AIF"
    assert [r["acquisition_type"] for r in _console_read(csv_path)] == ["AIF", "AIF"]


# ---------- 実 worker を通す ----------

def _use_fake_console_script(monkeypatch) -> list:
    """偽 Console へ向け、`build_msdial_cmd` が受けた `input_csv` を記録する。"""
    seen: list = []

    def build(exe, dataset_root, msdial_out_dir, method_file, save_project=False, *, input_csv=None):
        seen.append(input_csv)
        source = input_csv if input_csv is not None else dataset_root
        return [sys.executable, str(_FAKE_CONSOLE), "--scenario", "success",
                "-i", str(source), "-o", str(msdial_out_dir), "-m", str(method_file)]

    monkeypatch.setattr("metabolomix.console.runner.build_msdial_cmd", build)
    return seen


def _run(tmp_path, monkeypatch, *, manifest: bool):
    from metabolomix.pipeline import worker

    source = make_source(tmp_path / "source")
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["exe"]})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    if manifest:
        write_sample_manifest(source["root"] / "sample-manifest.tsv",
                              source_names=[f"S{i}.wiff" for i in range(8)],
                              groups=["control"] * 3 + ["treated"] * 3 + ["", ""],
                              roles=["sample"] * 6 + ["qc", "blank"])
    plan = service.plan_pipeline(source["root"], {"save_project": False})
    pipeline_path = Path(plan["pipeline_path"])
    seen = _use_fake_console_script(monkeypatch)
    result = worker.run_worker(pipeline_path)
    assert result["status"] == "completed", result.get("warnings")
    return source, pipeline_path, seen


def test_worker_passes_the_manifest_groups_to_console_and_records_the_csv(tmp_path, monkeypatch):
    _source, pipeline_path, seen = _run(tmp_path, monkeypatch, manifest=True)

    csv_path = pipeline_path / "negative" / "input" / input_csv.CSV_FILENAME
    assert seen == [csv_path]
    worker_state = json.loads(
        next((pipeline_path / "negative" / "console").glob("attempt-*/worker.json"))
        .read_text(encoding="utf-8"))
    assert worker_state["launch"]["command"][-5] == str(csv_path)        # 実際に起こした `-i`
    read = _console_read(csv_path)
    assert [r["class_id"] for r in read] == ["control"] * 3 + ["treated"] * 3 + ["qc", "blank"]
    assert [r["file_type"] for r in read] == ["Sample"] * 6 + ["QC", "Blank"]
    assert [r["analytical_order"] for r in read] == [str(i) for i in range(1, 9)]

    recorded = load_run(pipeline_path)["units"]["negative"]["upstream"]["console_input"]
    assert recorded["form"] == "csv"
    assert recorded["path"] == "negative/input/console-input.csv"           # pipeline_root 相対
    assert recorded["sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    assert recorded["classes"] == {"control": 3, "treated": 3, "qc": 1, "blank": 1}
    receipt = json.loads((pipeline_path / "results" / "receipt.json").read_text(encoding="utf-8"))
    unit = receipt["units"]["negative"]
    assert unit["console"]["input"] == recorded
    assert unit["method"]["console_input"] == "csv"
    status = recovery.read_status(pipeline_path)
    assert status["receipt"]["method_params"]["negative"]["console_input"] == "csv"


def test_worker_without_a_manifest_keeps_the_folder_and_says_so(tmp_path, monkeypatch):
    _source, pipeline_path, seen = _run(tmp_path, monkeypatch, manifest=False)

    assert seen == [None]
    assert not (pipeline_path / "negative" / "input" / input_csv.CSV_FILENAME).exists()
    recorded = load_run(pipeline_path)["units"]["negative"]["upstream"]["console_input"]
    assert recorded == {"form": "folder", "reason": "no_manifest"}
    receipt = json.loads((pipeline_path / "results" / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["units"]["negative"]["console"]["input"] == recorded
    assert receipt["units"]["negative"]["method"]["console_input"] == "folder"


def test_receipt_warns_when_the_manifest_changed_after_console_used_it(tmp_path, monkeypatch):
    source, pipeline_path, _seen = _run(tmp_path, monkeypatch, manifest=True)
    record = load_run(pipeline_path)
    assert receipt_mod.build_receipt(record, pipeline_path)["warnings"] == []

    write_sample_manifest(source["root"] / "sample-manifest.tsv",
                          source_names=[f"S{i}.wiff" for i in range(8)])   # 群を書き直した

    codes = [w["code"] for w in receipt_mod.build_receipt(record, pipeline_path)["warnings"]]
    assert codes == ["CONSOLE_INPUT_MANIFEST_CHANGED"]


def test_runs_recorded_before_the_csv_input_have_no_console_input(tmp_path):
    assert receipt_mod.unit_method_params(tmp_path, {"method": {}}, {"console_job_path": "x"}) is None
    method = tmp_path / "inputs" / "effective-method.txt"
    method.parent.mkdir()
    method.write_text("Ion mode: negative\n", encoding="ascii")
    inputs = {"method": {"effective_relative_path": "inputs/effective-method.txt"}}
    assert "console_input" not in receipt_mod.unit_method_params(tmp_path, inputs,
                                                                 {"console_job_path": "x"})
    assert "console_input" not in receipt_mod.unit_method_params(tmp_path, inputs)
