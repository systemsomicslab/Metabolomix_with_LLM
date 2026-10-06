"""v2（LC–MSメタボロミクス）のE2E用fixture（plan Task 14）。

**合成だが、モックの成功JSONは返さない。** mzTab-M と `.arf` を実ファイルとして
書き、実物の `parse_mztab` / `build_dataset_state` / `arf_reader.deserialize` に
通す。工程は `metabolomix.pipeline.engine.run_engine` を本物の v2 handler で回し、
`pipeline-run.json` へ保存された成果物を読み戻して検証する。

**どこを差し替えるか。** 上流4工程（`prepare_inputs` / `execute_console` /
`validate_outputs` / `load_dataset`）だけを、Console を起こさない test double へ
替える。ここは v1 と共有する経路で、実子プロセス Console を含む検証は
`tests/test_pipeline_process_lifecycle.py`（v1）と Task 15（実 Console 接続）の
担当。差し替えた `execute_console` は呼ばれた回数を数える——「binding を直して
再開しても Console を起こし直さない」という spec §6.2 の不変条件を、回数で
直接見るため。

**試料構成**（spec §9・plan Task 14）: 3生物群 × 3注入 = 9、QC 3、blank 1、
標準品 2 の計 15 注入。feature は対象2件＋内部標準1件＋無関係1件。
"""
from __future__ import annotations

import json
from pathlib import Path

import lz4.block
import msgpack
import numpy as np

from metabolomix.analysis.sample_manifest import FIELDS_V2
from metabolomix.console.profile_schema import validate_profile
from metabolomix.mztab.dataset_state import build_dataset_state
from metabolomix.mztab.reader import parse_mztab
from metabolomix.pipeline import engine, recovery, request as request_mod, store
from metabolomix.pipeline.service import build_handlers

__all__ = ["MetabolomicsHarness"]

#: 対象2件・内部標準1件・無関係1件。(feature_id, name, inchikey, mz, rt_min)
FEATURES = (
    ("1", "GABA", "BTCSSZJGUNDROE-UHFFFAOYSA-N", 104.0706, 1.2),
    ("2", "GABA-d6", "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N", 110.1100, 1.2),
    ("3", "Proline", "ONIBWKKTOPOVIA-BYPYZUCNSA-N", 116.0706, 2.4),
    ("4", None, None, 500.2000, 5.0),
)

#: 15注入の役割・群・生物試料ID。`sample-manifest.v2` の9列で書く。
def _injections():
    rows = []
    for group_index, group in enumerate(("g1", "g2", "g3")):
        for replicate in range(3):
            index = group_index * 3 + replicate
            rows.append({"sample_id": f"S{index}", "role": "sample", "group": group,
                         "qc_pool": "", "biological_sample_id": f"bio{index}"})
    for index in range(3):
        rows.append({"sample_id": f"QC{index}", "role": "qc", "group": "",
                     "qc_pool": "pool1", "biological_sample_id": ""})
    rows.append({"sample_id": "BLANK0", "role": "blank", "group": "",
                 "qc_pool": "", "biological_sample_id": ""})
    for index in range(2):
        rows.append({"sample_id": f"STD{index}", "role": "standard", "group": "",
                     "qc_pool": "", "biological_sample_id": ""})
    return rows


INJECTIONS = _injections()

#: 対象 feature の合成強度（3群×3注入）。plan Task 14 の値。
TARGET_VALUES = [1., 2., 3., 2., 3., 4., 4., 5., 6.]


def _abundance(feature_id: str, row_index: int, *, qc_fail: bool) -> float:
    """注入 × feature の合成強度。"""
    role = INJECTIONS[row_index]["role"]
    if feature_id == "1":
        if role == "sample":
            return TARGET_VALUES[row_index] * 100.0
        if role == "qc":
            # QC fail 用は [1,10,100]（RSD が跳ね上がる）。
            return ([1., 10., 100.] if qc_fail else [300., 303., 297.])[
                row_index - 9]
        if role == "blank":
            return 1.0
        return 300.0
    if feature_id == "2":                       # 内部標準（全注入で一定）
        return 200.0
    if feature_id == "3":
        return 50.0 + row_index
    return 10.0 + row_index                     # 無関係 feature


# ---------- mzTab ----------

def write_mztab_v2(path: Path, sources: list[Path], *, qc_fail: bool = False,
                   annotated: bool = True, ambiguous: bool = False) -> Path:
    """実物の `parse_mztab` / `build_dataset_state` を通る mzTab-M を書く。

    `ambiguous=True` は、対象と同じ m/z・RT・InChIKey を持つ feature をもう1件
    足す——binding が2候補になり `needs_input` になる入力。
    `annotated=False` は SME を全件 null にする（全 feature 未同定）。
    """
    path = Path(path)
    features = list(FEATURES)
    if ambiguous:
        features.append(("5", "GABA", "BTCSSZJGUNDROE-UHFFFAOYSA-N", 104.0707, 1.21))

    lines = [
        "MTD\tmzTab-version\t2.0.0-M",
        "MTD\tmzTab-mode\tComplete",
        "MTD\tmzTab-type\tQuantification",
        "MTD\tsoftware[1]\t[MS, MS:1003082, MS-DIAL, Msdial console 5.5.241113]",
        "MTD\tquantification_method\t[MS, MS:1001829, label-free raw feature quantitation, ]",
        "MTD\tsmall_molecule-quantification_unit\t[PRIDE, PRIDE:0000429, Abundance, ]",
    ]
    for i, src in enumerate(sources, start=1):
        lines.append(f"MTD\tms_run[{i}]-location\t{Path(src).resolve().as_uri()}")
        lines.append(f"MTD\tassay[{i}]-ms_run_ref\tms_run[{i}]")
        lines.append(f"MTD\tassay[{i}]\t{INJECTIONS[i - 1]['sample_id']}")

    n = len(sources)
    abundance_cols = "\t".join(f"abundance_assay[{i}]" for i in range(1, n + 1))
    lines.append(
        "SFH\tSMF_ID\tSME_ID_REFS\tSML_ID_REFS\tdatabase_identifier\tchemical_name\t"
        f"smiles\tinchi\texp_mass_to_charge\tretention_time_in_seconds\t{abundance_cols}")
    for feature_id, name, inchikey, mz, rt in features:
        values = "\t".join(f"{_abundance(feature_id, j, qc_fail=qc_fail):.1f}"
                           for j in range(n))
        lines.append(
            f"SMF\t{feature_id}\t{feature_id}\tSML:{feature_id}\tnull\t{name or 'null'}"
            f"\tnull\tnull\t{mz}\t{rt * 60.0}\t{values}")

    lines.append("SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tadduct_ions")
    for feature_id, *_ in features:
        lines.append(f"SML\t{feature_id}\t{feature_id}\tnull\tnull")

    lines.append("SEH\tSME_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_name\t"
                 "smiles\tinchi\tadduct_ion\tcharge\trank\t"
                 "best_id_confidence_measure\tbest_id_confidence_value")
    for feature_id, name, inchikey, *_ in features:
        identifier = inchikey if (annotated and inchikey) else "null"
        chemical_name = name if (annotated and name) else "null"
        lines.append(
            f"SME\t{feature_id}\tSMF:{feature_id}\t{identifier}\t{chemical_name}"
            f"\tnull\tnull\t[M+H]1+\t1\t1\tMS-DIAL total score\t0.95")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------- .arf（注入ごとの証拠） ----------

def _peak_row(file_id: int, sample_name: str, *, mz: float, rt: float,
              height: float, gap_filled: bool = False) -> list:
    row = [None] * 38
    row[0] = file_id
    row[1] = sample_name
    row[2] = -2 if gap_filled else file_id + 100
    row[3] = file_id + 1000
    row[10] = {}
    row[15] = [[1, [rt]], [3, [mz]]]
    row[18] = height
    row[20] = height * 2
    row[21] = height * 1.5
    row[22] = mz
    row[23] = 0
    row[24] = ""
    row[37] = [1.0, 20.0]
    return row


def write_arf(path: Path, sources: list[Path], *, qc_fail: bool = False,
              ambiguous: bool = False, standard_missing: bool = False) -> Path:
    """mzTab と同じ feature 順・同じ注入順の `.arf` を書く（注入ごとの実測値）。

    `standard_missing=True` は内部標準を標準品注入で gap-fill にする——分母が
    使えず ratio が locked になる case。
    """
    features = list(FEATURES)
    if ambiguous:
        features.append(("5", "GABA", "BTCSSZJGUNDROE-UHFFFAOYSA-N", 104.0707, 1.21))

    groups = []
    for feature_id, _name, _inchikey, mz, rt in features:
        rows = []
        for index, source in enumerate(sources):
            gap = (standard_missing and feature_id == "2"
                   and INJECTIONS[index]["role"] == "standard")
            rows.append(_peak_row(
                index, Path(source).stem, mz=mz + index * 1e-5,
                rt=rt + index * 1e-3,
                height=_abundance(feature_id, index, qc_fail=qc_fail),
                gap_filled=gap))
        groups.append(rows)

    inner = b"".join(msgpack.packb(item, use_bin_type=True)
                     for item in [[0, 0, *groups]])
    payload = (msgpack.packb(len(inner), use_bin_type=True)
               + lz4.block.compress(inner, store_size=False))
    path = Path(path)
    path.write_bytes(msgpack.packb(msgpack.ExtType(99, payload), use_bin_type=True))
    return path


# ---------- manifest / profile ----------

def write_manifest_v2(path: Path, sources: list[Path]) -> Path:
    lines = ["# schema = sample-manifest.v2", "\t".join(FIELDS_V2)]
    for index, source in enumerate(sources):
        row = INJECTIONS[index]
        lines.append("\t".join([
            row["sample_id"], Path(source).name, row["role"], row["group"], "B1",
            str(index + 1), row["qc_pool"], "true", row["biological_sample_id"]]))
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _null_evidence(paths) -> dict:
    return {p: {"value": None, "reason": "合成fixtureのため未記載",
                "tier": "proposed", "source_uri": None, "source_hash": None,
                "location": None} for p in paths}


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_profile(path: Path, *, statistics: list[dict] | None = None) -> Path:
    """profile と、それが指す method / 依存 / 実行体を**実ファイルで**書く。

    hash は書いたファイルから実測する——`resolve_profile_inputs` は宣言 hash と
    実測 hash の一致を要求するので、固定の偽 hash では profile 解決そのものが
    通らず、E2E が「解決を通っていない」ことに気付けない。
    """
    root = Path(path).parent
    (root / "method").mkdir(parents=True, exist_ok=True)
    (root / "lib").mkdir(parents=True, exist_ok=True)
    (root / "tools").mkdir(parents=True, exist_ok=True)

    method_path = root / "method" / "param.txt"
    method_path.write_text(
        "\n".join([
            "Ion mode: Positive",
            "Target omics: Metabolomics",
            "Acquisition type: DDA",
            # methodの参照キーはmethodファイル自身の親基準で解決される。
            "Msp file path: ../lib/l.msp",
        ]) + "\n",
        encoding="utf-8", newline="\n")
    library_path = root / "lib" / "l.msp"
    library_path.write_text("NAME: synthetic\n", encoding="utf-8", newline="\n")
    exe_path = root / "tools" / "MSDIALCUI.exe"
    exe_path.write_text("fake console executable placeholder", encoding="ascii")

    profile = {
        "schema": "lcms-profile.v1", "profile_id": "synthetic-metabolomics",
        "revision": 1, "omics": "metabolomics",
        "acquisition": {
            "separation": "lc", "acquisition_type": "dda", "polarity": "positive",
            "instrument": None,
            "lc": {"column": None, "mobile_phase_a": None, "mobile_phase_b": None,
                   "flow_rate_ul_min": None, "column_temperature_c": None,
                   "gradient_profile": None},
            "ms_range": {"ms1_low_mz": 50.0, "ms1_high_mz": 1500.0,
                         "ms2_low_mz": None, "ms2_high_mz": None},
            "sample_matrix": None, "scope": "synthetic single-condition batch"},
        "software": {"msdial_version": "5.x", "executable_path": "tools/MSDIALCUI.exe",
                     "executable_sha256": _sha256(exe_path),
                     "adapter_version": "1.0.0"},
        "processing": {
            "method_path": "method/param.txt",
            "method_sha256": _sha256(method_path),
            "measure": "peak_height",
            "dependencies": [{"dependency_id": "msp-lib-1", "kind": "msp",
                              "method_key": "Msp file path", "path": "lib/l.msp",
                              "sha256": _sha256(library_path), "required": True}],
            "effective_settings": {}},
        "analysis_recipe": {
            "statistics": list(statistics or []),
            "internal_standards": [{"target_id": "gaba",
                                    "standard_target_id": "gaba_d6"}]},
        "feature_targets": {
            "gaba": {
                "compound_identifiers": {"name": "GABA",
                                         "inchikey": FEATURES[0][2]},
                "adduct": "[M+H]+", "charge": 1, "expected_mz": FEATURES[0][3],
                "mz_tolerance_ppm": 20.0, "expected_rt_min": FEATURES[0][4],
                "rt_tolerance_min": 0.2,
                "required_evidence": [{"kind": "mass_rt"}]},
            "gaba_d6": {
                "compound_identifiers": {"name": "GABA-d6",
                                         "inchikey": FEATURES[1][2]},
                "adduct": "[M+H]+", "charge": 1, "expected_mz": FEATURES[1][3],
                "mz_tolerance_ppm": 20.0, "expected_rt_min": FEATURES[1][4],
                "rt_tolerance_min": 0.2,
                "required_evidence": [{"kind": "mass_rt"}]}},
        "matrix_recipes": {
            "default": {"base": "peak_height", "normalize": "none",
                        "drift_correct": False, "filter": None, "impute": "none"},
            "ratio": {"base": "internal_standard_ratio", "normalize": "none",
                      "drift_correct": False, "filter": None, "impute": "none"}},
        "qc_policy": {
            "pooled_qc_rsd": {
                "scope": "all_features", "target_ids": None, "required": True,
                "threshold": {"operator": "<=", "value": 30.0},
                "evidence_requirement": False, "minimum_pass_fraction": 0.5},
            "detection_rate": {
                "scope": "all_features", "target_ids": None, "required": False,
                "threshold": {"operator": ">=", "value": 0.5},
                "evidence_requirement": False, "minimum_pass_fraction": 0.5}},
        "evidence": _null_evidence((
            "acquisition.instrument", "acquisition.sample_matrix",
            "acquisition.lc.column", "acquisition.lc.mobile_phase_a",
            "acquisition.lc.mobile_phase_b", "acquisition.lc.flow_rate_ul_min",
            "acquisition.lc.column_temperature_c",
            "acquisition.lc.gradient_profile",
            "acquisition.ms_range.ms2_low_mz",
            "acquisition.ms_range.ms2_high_mz")),
        "validation": {"status": "draft", "scope": "synthetic single-condition batch",
                       "certificate_path": None, "certificate_sha256": None},
    }
    validate_profile(profile)
    path = Path(path)
    path.write_text(json.dumps(profile), encoding="utf-8")
    return path


DEFAULT_STATISTICS = [
    {"statistic_id": "anova", "kind": "anova_tukey", "matrix_recipe_id": "default",
     "transform": "none", "feature_scope": {"mode": "all_eligible"},
     "groups": ["g1", "g2", "g3"], "alpha": 0.05},
    {"statistic_id": "welch", "kind": "welch", "matrix_recipe_id": "ratio",
     "transform": "none", "feature_scope": {"mode": "all_eligible"},
     "reference_group": "g1", "test_group": "g3",
     "q_threshold": 0.05, "log2fc_threshold": 1.0},
]


# ---------- harness ----------

class MetabolomicsHarness:
    """v2 の工程を `run_engine` で回し、保存された成果物を読み戻す harness。

    `launch_count` は fake Console（`execute_console`）が呼ばれた回数。binding の
    訂正で再開しても増えないことを、この数で直接見る。
    """

    def __init__(self, tmp_path):
        self.tmp_path = Path(tmp_path)
        self.source_root = self.tmp_path / "source"
        self.source_root.mkdir(parents=True, exist_ok=True)
        self.launch_count = 0
        self.pipeline_path: Path | None = None
        self._dataset = None
        self._sources: list[Path] = []
        self._arf_path: Path | None = None

    # ---------- 入力の用意 ----------

    def _write_inputs(self, *, qc_fail: bool, annotated: bool, ambiguous: bool,
                      standard_missing: bool) -> None:
        self._sources = []
        for row in INJECTIONS:
            path = self.source_root / f"{row['sample_id']}.wiff"
            path.write_bytes(b"synthetic raw")
            self._sources.append(path)

        self.mztab_path = write_mztab_v2(
            self.source_root / "Height_synthetic.mzTab", self._sources,
            qc_fail=qc_fail, annotated=annotated, ambiguous=ambiguous)
        self._arf_path = write_arf(
            self.source_root / "AlignmentResult_PeakProperties.arf", self._sources,
            qc_fail=qc_fail, ambiguous=ambiguous,
            standard_missing=standard_missing)
        self.manifest_path = write_manifest_v2(
            self.source_root / "sample-manifest.tsv", self._sources)
        self.profile_path = write_profile(
            self.source_root / "profile.json", statistics=DEFAULT_STATISTICS)

    def _write_execution_record(self) -> Path:
        """Console の終了証跡。`recovery.prepare_resume` がこれを読む。

        再開時に「上流が本当に終わっていたか」を判定する材料で、無いと
        `EXECUTION_UNRESOLVED` になる。偽 Console でも**実ファイルとして**書く
        ——record の中だけで完結させると、再開経路の実際の読み取りを
        試験できない。
        """
        from tests.pipeline_fixtures import execution_record

        path = self.mztab_path.parent / "execution-result.json"
        path.write_text(json.dumps(execution_record(
            execution_id="exec-synthetic", exit_code=0,
            outcome="completed")), encoding="utf-8")
        return path

    def _load_dataset(self):
        """実物の `parse_mztab` / `build_dataset_state` で DatasetState を作る。"""
        parsed = parse_mztab(self.mztab_path)
        ds = build_dataset_state(parsed, self.mztab_path.name, self.mztab_path)
        ds.assay_sources = {f"abundance_assay[{i + 1}]": str(path)
                            for i, path in enumerate(self._sources)}
        ds.artifact_paths = {"peak_matrix_source": [str(self._arf_path)]}
        ds.source_verification = "verified"
        return ds

    # ---------- 差し替えた上流4工程 ----------

    def _upstream_handlers(self) -> dict:
        harness = self

        def execute_console(context):
            # 「Console を起こした」回数。binding 訂正の再開で増えないことを見る。
            harness.launch_count += 1
            harness._write_execution_record()
            return {"status": "succeeded", "result_refs": [], "warnings": [],
                    "error": None,
                    "record_updates": {"upstream": {
                        "console_job_path": str(harness.mztab_path),
                        "execution_id": "exec-synthetic",
                        "verification": {"status": "completed"}}}}

        def validate_outputs(context):
            return {"status": "succeeded", "result_refs": [], "warnings": [],
                    "error": None}

        def load_dataset(context):
            ds = harness._load_dataset()
            harness._dataset = ds
            context["runtime"]["dataset"] = ds
            return {"status": "succeeded", "result_refs": [], "warnings": [],
                    "error": None}

        return {"execute_console": execute_console,
                "validate_outputs": validate_outputs,
                "load_dataset": load_dataset}

    def _handlers(self) -> dict:
        handlers = dict(build_handlers())
        handlers.update(self._upstream_handlers())
        return handlers

    # ---------- 実行 ----------

    def _request(self, *, statistics=None, standard_assays=None,
                 feature_bindings=None) -> dict:
        explicit = {
            "schema": "pipeline-request.v2",
            "profile_file": str(self.profile_path),
            "execution_purpose": "validation",
            "sample_manifest": str(self.manifest_path),
            "statistics": statistics if statistics is not None else DEFAULT_STATISTICS,
            "standard_assays": standard_assays or {},
        }
        if feature_bindings:
            explicit["feature_bindings"] = feature_bindings
        profile = json.loads(self.profile_path.read_text(encoding="utf-8"))
        return request_mod.resolve_request(self.source_root, explicit,
                                           profile=validate_profile(profile))

    def run(self, *, binding_mode: str = "unique", annotated: bool = True,
            qc_fail: bool = False, standard_missing: bool = False,
            statistics=None) -> dict:
        """入力を書き、run を作り、engine を1回分回す。"""
        self._write_inputs(qc_fail=qc_fail, annotated=annotated,
                           ambiguous=(binding_mode == "ambiguous"),
                           standard_missing=standard_missing)
        request = self._request(statistics=statistics)
        inputs = {"fingerprint": "synthetic-inputs", "source_root": str(self.source_root)}
        self.pipeline_path = store.create_run(self.source_root, request, inputs)
        return self._run_engine()

    def _run_engine(self) -> dict:
        result = engine.run_engine(self.pipeline_path, self._handlers())
        result = dict(result)
        result["allowed_binding_update"] = self._binding_update()
        return result

    def resume(self, updates: dict) -> dict:
        """要求を更新して再開する（Console は起こし直さない）。"""
        recovery.prepare_resume(self.pipeline_path, updates=updates)
        return self._run_engine()

    # ---------- 観測 ----------

    def record(self) -> dict:
        return store.load_run(self.pipeline_path)

    def result_data(self, output_name: str):
        """保存済み成果物を**実ファイルから**読み戻す（最新の1件）。"""
        record = self.record()
        matches = [r for r in (record.get("results") or [])
                   if r.get("output_name") == output_name]
        if not matches:
            return None
        path = Path(record["identity"]["pipeline_root"]) / matches[-1]["relative_path"]
        if path.suffix == ".json":
            return json.loads(path.read_text(encoding="utf-8")).get("data")
        return path.read_text(encoding="utf-8")

    def output_names(self) -> list[str]:
        return [r.get("output_name") for r in (self.record().get("results") or [])]

    def stage_status(self, stage_id: str) -> str:
        return self.record()["stages"][stage_id]["status"]

    def _binding_update(self) -> dict:
        """needs_input のとき、利用者が渡せる binding 訂正（qualified候補の1件）。"""
        bindings = self.result_data("feature_bindings")
        if bindings is None:
            return {}
        update = {}
        for target_id, binding in (bindings.get("bindings") or {}).items():
            if binding.get("status") == "resolved":
                continue
            qualified = [c["feature_id"] for c in binding.get("candidates") or []
                         if c.get("qualified")]
            if qualified:
                update[target_id] = {
                    "feature_id": sorted(qualified)[0],
                    "reason": "合成fixture: 最初のqualified候補を採用"}
        if not update:
            return {}
        # 要求側の payload は「どの dataset に対する選択か」を自分で宣言する。
        return {"dataset_hash": bindings.get("dataset_hash"),
                "selections": update}
