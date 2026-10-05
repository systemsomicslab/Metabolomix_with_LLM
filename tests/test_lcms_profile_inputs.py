"""`lipidmix.console.profiles` / `profile_adapter` の統合テスト（spec §5, §5.1, §6.1）。

`lipidmix.console.profile_schema` の構造検証は `tests/test_lcms_profile_schema.py`。
ここは実ファイル（method・依存・実行体・raw）を使う統合層: `load_profile` →
`resolve_profile_inputs` → `snapshot_profile` の一気通貫、raw構成ファイルの
hash（§6.1）、method_key単位の原本→実効値の差分記録（§5.1）を検証する。
`method_file.METABOLOMICS_REFERENCE_KEYS`はmetabolomics/adapter経路
（`resolve_profile_inputs`）だけが使う拡張キー集合で、lipidomics v1経路
（`pipeline.inputs.inspect_inputs`）の既定挙動（`REFERENCE_KEYS`＝LBMのみ）は
本タスクで変えていない——それを固定するテストもここに含む。fixtureはこの
テスト自身が全て`tmp_path`配下に組み立てる。
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from lipidmix.console import method_file as method_file_mod
from lipidmix.console import profile_adapter
from lipidmix.console import profiles
from lipidmix.core.atomic_io import DomainError


# ---------- fixtureヘルパ ----------

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content.encode("ascii"))
    return path


def _write_fake_exe(base: Path) -> Path:
    """実行体+同梱DLL/設定を並べる（実行環境manifestの列挙対象）。"""
    exe_dir = base / "msdial_app"
    exe = _write(exe_dir / "MSDIALCUI.exe", "fake-msdial-console-exe")
    _write(exe_dir / "MsdialCore.dll", "fake-msdialcore-dll")
    _write(exe_dir / "MSDIALCUI.exe.config", "<configuration/>")
    _write(exe_dir / "MSDIALCUI.pdb", "not-a-companion")  # .pdbは対象外（無視される）
    return exe


def _msp_dependency(msp_path: Path, *, required: bool = True) -> dict:
    return {
        "dependency_id": "msp-lib-1", "kind": "msp", "method_key": method_file_mod.MSP_KEY,
        "path": str(msp_path), "sha256": _sha256(msp_path), "required": required,
    }


def _text_dependency(text_path: Path, *, required: bool = True) -> dict:
    return {
        "dependency_id": "text-db-1", "kind": "text_identification",
        "method_key": method_file_mod.TEXT_DB_KEY,
        "path": str(text_path), "sha256": _sha256(text_path), "required": required,
    }


def _build_profile(tmp_path, *, method_lines: str, dependencies: list[dict],
                    exe: Path, polarity: str = "positive", msdial_version: str = "5.x") -> dict:
    method = _write(tmp_path / "method" / "params.txt", method_lines)
    return {
        "schema": "lcms-profile.v1",
        "profile_id": "kanzo-lcms-metabolomics",
        "revision": 1,
        "omics": "metabolomics",
        "acquisition": {
            "separation": "lc",
            "acquisition_type": "dda",
            "polarity": polarity,
            "instrument": None,
            "lc": {
                "column": None, "mobile_phase_a": None, "mobile_phase_b": None,
                "flow_rate_ul_min": None, "column_temperature_c": None, "gradient_profile": None,
            },
            "ms_range": {
                "ms1_low_mz": 50.0, "ms1_high_mz": 1500.0,
                "ms2_low_mz": None, "ms2_high_mz": None,
            },
            "sample_matrix": None,
            "scope": "single LC method, DDA acquisition, peak height only",
        },
        "software": {
            "msdial_version": msdial_version,
            "executable_path": str(exe),
            "executable_sha256": _sha256(exe),
            "adapter_version": "1.0.0",
        },
        "processing": {
            "method_path": str(method),
            "method_sha256": _sha256(method),
            "measure": "peak_height",
            "dependencies": dependencies,
            "effective_settings": {},
        },
        "analysis_recipe": {"statistics": [], "internal_standards": []},
        "feature_targets": {},
        "matrix_recipes": {
            "default": {"base": "peak_height", "normalize": "none",
                        "drift_correct": False, "filter": None, "impute": "none"},
        },
        "qc_policy": {},
        "evidence": {
            "acquisition.instrument": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.sample_matrix": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.lc.column": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.lc.mobile_phase_a": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.lc.mobile_phase_b": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.lc.flow_rate_ul_min": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.lc.column_temperature_c": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.lc.gradient_profile": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.ms_range.ms2_low_mz": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
            "acquisition.ms_range.ms2_high_mz": {
                "value": None, "reason": "テスト用fixtureのため未記載。", "tier": "proposed",
                "source_uri": None, "source_hash": None, "location": None,
            },
        },
        "validation": {
            "status": "draft",
            "scope": "single LC method, DDA acquisition, peak height only",
            "certificate_path": None, "certificate_sha256": None,
        },
    }


def _nested_root(tmp_path: Path) -> Path:
    """profile/依存/実行体を tmp_path 直下ではなく1段ネストして置く（他テストとの
    兄弟フォルダ混在事故を避ける既存の流儀 — tests/test_console_plan_method.py参照）。

    直下に生データ（`.wiff`）を1件置く——`resolve_profile_inputs`が
    `pipeline.inputs`の既存raw形式選択を再利用するため、raw構成ファイルが
    ちょうど1形式・1件無いと形式選択自体が失敗する（`DATASET_SELECTION_REQUIRED`
    等）。method/依存/実行体はすべて`method/` `library/` `msdial_app/`という
    サブフォルダの下に置くため（`_build_profile`/`_write_fake_exe`参照）、
    直下形式選択（非再帰）とは衝突しない。
    """
    root = tmp_path / "fixture"
    root.mkdir(exist_ok=True)
    (root / "sample1.wiff").write_bytes(b"raw-instrument-data-v1")
    return root


# ---------- brief記載のRED: hash_files ----------

def test_hash_detects_stat_preserving_change(tmp_path):
    """statを保ったまま内容だけ変えても検出する（brief記載のコード例そのもの）。"""
    p = tmp_path / "raw.wiff"
    p.write_bytes(b"aaaa")
    stat = p.stat()
    old = profiles.hash_files([p])
    p.write_bytes(b"bbbb")
    os.utime(p, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert profiles.hash_files([p]) != old


def test_hash_files_matches_plain_sha256(tmp_path):
    p = tmp_path / "a.txt"
    p.write_bytes(b"hello world")
    assert profiles.hash_files([p])[str(p)] == hashlib.sha256(b"hello world").hexdigest()


def test_hash_files_detects_mutation_during_single_call(tmp_path, monkeypatch):
    """brief「二回のstatが異なるhash計算はINPUT_CHANGEDとして拒否する」の、
    単一呼び出し内での実際の検出経路（読み込みの前後でstatを比較する分岐）を
    直接確認する（controller裁定 fix round 2 minor 1）。

    2回の`hash_files`呼び出しを比べるだけの既存テスト（
    `test_hash_detects_stat_preserving_change`）は、内容が変わればhashも
    当然変わるので、この前後stat比較の分岐が無くても通ってしまう——ここでは
    `Path.stat`をモックし、`hash_files`内の1回の呼び出しの最中に
    （読み込みの前後で）statが変化したように見せかけ、その分岐自体が
    `INPUT_CHANGED`を送出することを確認する。
    """
    p = tmp_path / "raw.wiff"
    p.write_bytes(b"aaaa")
    real_stat = Path.stat
    calls = {"n": 0}

    def fake_stat(self, *args, **kwargs):
        calls["n"] += 1
        real = real_stat(self, *args, **kwargs)
        if calls["n"] == 1:
            return real
        # 2回目（読み込み後のstat）だけ、実際には起きていないサイズ変化を装う。
        return SimpleNamespace(st_size=real.st_size + 1, st_mtime_ns=real.st_mtime_ns)

    monkeypatch.setattr(Path, "stat", fake_stat)
    with pytest.raises(DomainError, match="INPUT_CHANGED"):
        profiles.hash_files([p])
    assert calls["n"] == 2  # 実際に前後2回statしていることの確認


# ---------- adapter_capabilities ----------

def test_adapter_capabilities_returns_msp_lbm_text_rt_keys():
    capabilities = profile_adapter.adapter_capabilities("msdial5")
    assert capabilities["adapter_id"] == "msdial5"
    assert capabilities["dependency_keys"] == {
        "msp": "Msp file path",
        "lbm": "Lbm file path",
        "text_identification": "Text DB file path",
        "rt_reference": "Compounds library file path for RT correction",
    }
    assert "wiff" in capabilities["raw_formats"]
    # 注入ごとの証拠は`.arf`のAlignedPeakProperties（スポット内に注入行を持つ）
    # から取る。`.pai2`はアライメント特徴への対応を持たないため、名前やm/z近傍の
    # 再結合が要り、spec §9.1が禁じる「名前だけの結合」になる。
    assert capabilities["evidence_reader"] == "arf"
    assert capabilities["evidence_rt_unit"] == "minute"


def test_adapter_capabilities_unknown_id_rejected():
    with pytest.raises(DomainError, match="PROFILE_ADAPTER_UNSUPPORTED"):
        profile_adapter.adapter_capabilities("msdial99")


def test_adapter_capabilities_returns_independent_copies():
    """呼び出し側がdictを書き換えてもレジストリ自体は汚染されない。"""
    first = profile_adapter.adapter_capabilities("msdial5")
    first["dependency_keys"]["msp"] = "tampered"
    second = profile_adapter.adapter_capabilities("msdial5")
    assert second["dependency_keys"]["msp"] == "Msp file path"


# ---------- load_profile ----------

def test_profile_dependencies_and_raw_inventory_use_separate_roots(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                             dependencies=[], exe=exe)
    raw_root = tmp_path / "batch"
    raw_root.mkdir()
    (raw_root / "batch.wiff").write_bytes(b"different-batch")
    plan = profiles.resolve_profile_inputs(profile, root, raw_root=raw_root)
    assert [entry["relative_path"] for entry in plan["raw_files"]["files"]] == ["batch.wiff"]
    assert plan["source_root"] == str(raw_root)
    assert Path(plan["method"]["source_path"]).is_relative_to(root)


def test_routine_requires_observed_certificate_evidence(tmp_path):
    root = _nested_root(tmp_path)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                             dependencies=[], exe=_write_fake_exe(root))
    certificate = {}
    profile["validation"] = {"status": "validated", "scope": "test",
        "certificate_path": "certificate.json", "certificate_sha256": profiles.canonical_hash(certificate)}
    (root / "certificate.json").write_text(json.dumps(certificate), encoding="utf-8")
    path = root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        profiles.load_profile(path, "routine")
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        profiles.certificate_routine_overrides(profile, path)


@pytest.mark.parametrize("failure", [None, "failed_criterion", "changed_output", "changed_profile"])
def test_routine_validates_independent_evidence(tmp_path, failure):
    root = _nested_root(tmp_path)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                             dependencies=[], exe=_write_fake_exe(root))
    observed = {}
    certificate = {"schema": "lcms-profile-validation.v1",
        "profile_content_sha256": profiles.profile_content_hash(profile),
        "criteria": [{"criterion_id": "check", "status": "fail" if failure == "failed_criterion" else "pass",
                      "required": True, "detail": None}],
        "performed_by": "fixture", "performed_at": "2026-09-16T00:00:00Z", "scope": "test"}
    for cert_key, observed_key in [("dependency_hashes", "dependencies"),
            ("fixed_input_hashes", "fixed_inputs"), ("reference_file_hashes", "reference_files"),
            ("validation_output_hashes", "validation_outputs")]:
        evidence = root / observed_key
        evidence.write_bytes(observed_key.encode())
        digest = profiles.hash_files([evidence])[str(evidence)]
        certificate[cert_key] = {"evidence": digest}
        if failure == "changed_output" and observed_key == "validation_outputs":
            evidence.write_bytes(b"changed")
        observed[observed_key] = {"evidence": profiles.hash_files([evidence])[str(evidence)]}
    profile["validation"] = {"status": "validated", "scope": "test", "certificate_path": "certificate.json",
        "certificate_sha256": profiles.canonical_hash(certificate)}
    if failure == "changed_profile":
        profile["revision"] += 1
    (root / "certificate.json").write_text(json.dumps(certificate), encoding="utf-8")
    path = root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    if failure:
        with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
            profiles.load_profile(path, "routine", observed_hashes=observed)
    else:
        assert profiles.load_profile(path, "routine", observed_hashes=observed)["validation"]["status"] == "validated"


def test_load_profile_rejects_unknown_purpose(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe)
    path = root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    with pytest.raises(DomainError, match="PROFILE_PURPOSE_INVALID"):
        profiles.load_profile(path, "bogus")


def test_load_profile_routine_requires_validated_status(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe)
    path = root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    with pytest.raises(DomainError, match="PROFILE_NOT_VALIDATED"):
        profiles.load_profile(path, "routine")
    # purpose="validation"はdraftのままでよい。
    loaded = profiles.load_profile(path, "validation")
    assert loaded["validation"]["status"] == "draft"


def test_load_profile_does_not_absolutize_paths(tmp_path):
    """`profile_content_hash`の環境非依存性を守るため、パスは宣言どおり返す。"""
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe)
    path = root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    loaded = profiles.load_profile(path, "validation")
    assert loaded["processing"]["method_path"] == profile["processing"]["method_path"]


# ---------- resolve_profile_inputs: 正常系 ----------

def test_resolve_profile_inputs_metabolomics_without_lbm(tmp_path):
    """LBMなしMSP/TXT: metabolomicsはlbmを依存に含めなくてよい。"""
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    text_db = _write(root / "library" / "text1.txt", "text-db-content")
    dependencies = [_msp_dependency(msp), _text_dependency(text_db)]
    profile = _build_profile(
        root, method_lines="Ion mode: Positive\nTarget omics: Metabolomics\nAcquisition type: DDA\n",
        dependencies=dependencies, exe=exe)

    plan = profiles.resolve_profile_inputs(profile, root)

    assert plan["adapter"]["adapter_id"] == "msdial5"
    kinds_present = {d["kind"]: d["present"] for d in plan["dependencies"]}
    assert kinds_present == {"msp": True, "text_identification": True}
    assert plan["execution_environment"]["companions"]["MsdialCore.dll"]
    assert "MSDIALCUI.pdb" not in plan["execution_environment"]["companions"]
    assert plan["polarity"] == "positive"
    # raw指紋（spec §6.1「上流指紋にraw全構成ファイルSHA-256、採用形式...を含める」）。
    assert plan["raw_files"]["selected_format"] == "wiff"
    raw_entry = next(e for e in plan["raw_files"]["files"] if e["relative_path"] == "sample1.wiff")
    assert raw_entry["sha256"] == _sha256(root / "sample1.wiff")


def test_resolve_profile_inputs_optional_dependency_missing_is_not_fatal(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dep = _msp_dependency(msp, required=False)
    dep["path"] = str(root / "library" / "missing.msp")  # 実在しない、かつrequired=False
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[dep], exe=exe)

    plan = profiles.resolve_profile_inputs(profile, root)
    assert plan["dependencies"][0]["present"] is False


# ---------- RED: 必須依存欠落 ----------

def test_resolve_profile_inputs_required_dependency_missing_is_incomplete(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dep = _msp_dependency(msp, required=True)
    dep["path"] = str(root / "library" / "missing.msp")
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[dep], exe=exe)

    with pytest.raises(DomainError, match="PROFILE_INCOMPLETE"):
        profiles.resolve_profile_inputs(profile, root)


# ---------- RED: キー名を推測してMSPをLBM欄へ入れない ----------

def test_resolve_profile_inputs_rejects_msp_in_lbm_slot(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dep = _msp_dependency(msp)
    dep["method_key"] = method_file_mod.LBM_KEY  # kind=mspなのにLBM欄を騙る
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[dep], exe=exe)

    with pytest.raises(DomainError, match="PROFILE_METHOD_CONFLICT"):
        profiles.resolve_profile_inputs(profile, root)


# ---------- RED: 依存の method_key 衝突（同じ枠を2件が取り合う） ----------

def test_resolve_profile_inputs_rejects_two_dependencies_same_method_key(tmp_path):
    """methodファイルの1つのキーは1行しか書けない。2つのdependencyが同じ
    method_key（ここではkindも同じ`msp`）を宣言したら、後勝ちで片方を黙って
    捨てるのではなく、一意に決まらないこと自体をエラーにする。

    schema側（`validate_profile`）は`dependency_id`の重複だけを拒否し、`kind`や
    `method_key`の重複は禁止していない（`lipidmix/console/profile_schema.py`の
    `_validate_dependency`参照）——ここでの衝突検出はTask 2（`resolve_profile_inputs`）
    が新たに持つべき責務であり、Task 1のスキーマ検証を重複させるものではない。
    """
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp_a = _write(root / "library" / "lib_a.msp", "msp-content-a")
    msp_b = _write(root / "library" / "lib_b.msp", "msp-content-b")
    dep_a = _msp_dependency(msp_a)
    dep_a["dependency_id"] = "msp-lib-a"
    dep_b = _msp_dependency(msp_b)
    dep_b["dependency_id"] = "msp-lib-b"
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[dep_a, dep_b], exe=exe)

    with pytest.raises(DomainError, match="PROFILE_METHOD_CONFLICT"):
        profiles.resolve_profile_inputs(profile, root)


def test_resolve_profile_inputs_allows_two_dependencies_different_kind(tmp_path):
    """`kind`の重複そのものは禁止しない（spec §5.1・schemaのどちらも`kind`の
    重複を禁じていない）。衝突判定の基準は`kind`ではなく`method_key`——
    2件が実際に異なるmethodファイルの行（method_key）を占めるなら共存できる。

    msdial5アダプタは`kind`→`method_key`が1対1（`profile_adapter.py`の
    `dependency_keys`）なので、「`kind`が同じで`method_key`が異なる」組み合わせは
    このアダプタでは構成できない——衝突判定を`kind`単位ではなく`method_key`単位に
    した設計そのものが、この1対1の下では「`kind`が異なれば必ず`method_key`も
    異なり、同じ`kind`なら必ず`method_key`も同じ」という関係を素直に反映する。
    """
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    text_db = _write(root / "library" / "text1.txt", "text-db-content")
    dependencies = [_msp_dependency(msp), _text_dependency(text_db)]
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=dependencies, exe=exe)

    plan = profiles.resolve_profile_inputs(profile, root)
    assert {d["kind"] for d in plan["dependencies"]} == {"msp", "text_identification"}


# ---------- RED: DDA/極性矛盾 ----------

def test_resolve_profile_inputs_rejects_polarity_conflict(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(
        root, method_lines="Ion mode: Negative\n", dependencies=[], exe=exe, polarity="positive")

    with pytest.raises(DomainError, match="PROFILE_METHOD_CONFLICT"):
        profiles.resolve_profile_inputs(profile, root)


def test_resolve_profile_inputs_rejects_acquisition_type_conflict(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(
        root, method_lines="Ion mode: Positive\nAcquisition type: SWATH\n",
        dependencies=[], exe=exe, polarity="positive")

    with pytest.raises(DomainError, match="PROFILE_METHOD_CONFLICT"):
        profiles.resolve_profile_inputs(profile, root)


def test_resolve_profile_inputs_rejects_omics_conflict(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(
        root, method_lines="Ion mode: Positive\nTarget omics: Lipidomics\n",
        dependencies=[], exe=exe, polarity="positive")

    with pytest.raises(DomainError, match="PROFILE_METHOD_CONFLICT"):
        profiles.resolve_profile_inputs(profile, root)


# ---------- RED: 未知adapter ----------

def test_resolve_profile_inputs_unknown_adapter_version(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe, msdial_version="99.0")

    with pytest.raises(DomainError, match="PROFILE_ADAPTER_UNSUPPORTED"):
        profiles.resolve_profile_inputs(profile, root)


# ---------- RED: 同サイズ/mtime改変 (dependency/method/exe hash mismatch) ----------

def test_resolve_profile_inputs_rejects_tampered_method_content(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe)
    # profileが記録した後にメソッドファイルの中身だけ差し替える（サイズ違いでも可）。
    method_path = Path(profile["processing"]["method_path"])
    stat = method_path.stat()
    method_path.write_bytes(b"Ion mode: Negative\n")
    os.utime(method_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    with pytest.raises(DomainError, match="INPUT_CHANGED"):
        profiles.resolve_profile_inputs(profile, root)


# ---------- RED (fix round 2, finding 1 / A06): raw構成ファイルのSHA-256 ----------

def test_resolve_profile_inputs_raw_file_hash_detects_stat_preserving_tamper(tmp_path):
    """A06: 生データがsize/mtimeを保ったまま改変されても検出できること。

    method/依存/実行体だけでなく、実際の計測ファイル（.wiff等）自身の内容が
    改変されたことを、statではなくhash（`profiles.hash_files`経由）で検出する。
    `resolve_profile_inputs`を改変前後で2回呼び、`plan["raw_files"]`内の
    同じ相対パスのsha256が変わっていることを確認する——`test_hash_detects_
    stat_preserving_change`と同じ性質を、実際の統合経路（`resolve_profile_inputs`
    が返すraw指紋）を通して確認する。
    """
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe)

    raw_path = root / "sample1.wiff"
    stat = raw_path.stat()
    plan_before = profiles.resolve_profile_inputs(profile, root)
    before_entry = next(e for e in plan_before["raw_files"]["files"]
                        if e["relative_path"] == "sample1.wiff")

    # 同サイズの別内容（A06はサイズ・mtimeを保ったままの改変を想定する）。
    tampered = b"raw-instrument-data-v2"
    assert len(tampered) == stat.st_size
    raw_path.write_bytes(tampered)
    os.utime(raw_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert raw_path.stat().st_size == stat.st_size  # サイズも保ったまま（A06どおり）

    plan_after = profiles.resolve_profile_inputs(profile, root)
    after_entry = next(e for e in plan_after["raw_files"]["files"]
                       if e["relative_path"] == "sample1.wiff")

    assert before_entry["sha256"] != after_entry["sha256"]


def test_resolve_profile_inputs_raw_files_include_sidecars(tmp_path):
    """既存のsidecar列挙（`.wiff`+`.wiff.scan`等）をそのまま再利用しているかを確認する。"""
    root = _nested_root(tmp_path)
    # .wiff.scanは`_companions_of`のCOMPANION_RULESに実在する随伴ファイル。
    (root / "sample1.wiff.scan").write_bytes(b"scan-sidecar")
    exe = _write_fake_exe(root)
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=[], exe=exe)

    plan = profiles.resolve_profile_inputs(profile, root)
    relative_paths = {e["relative_path"] for e in plan["raw_files"]["files"]}
    assert "sample1.wiff.scan" in relative_paths


# ---------- snapshot_profile ----------

def test_snapshot_profile_writes_effective_method_with_absolute_dependency_paths(tmp_path):
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dependencies = [_msp_dependency(msp)]
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=dependencies, exe=exe)
    plan = profiles.resolve_profile_inputs(profile, root)

    run_dir = tmp_path / "run1"
    snapshot = profiles.snapshot_profile(plan, run_dir)

    effective = run_dir / snapshot["effective_method_relative_path"]
    assert effective.is_file()
    text = effective.read_text(encoding="ascii")
    assert f"Msp file path: {msp.resolve()}" in text or f"Msp file path: {msp}" in text


# ---------- RED (fix round 2, finding 2): 差分(原本宣言値→実効値)を記録する ----------

def test_snapshot_profile_records_method_key_diff_with_stale_original(tmp_path):
    """spec §5.1「原本hash、書換え後hash、差分を記録する」。method原本が古い/
    別の値を宣言していた場合、snapshotはその宣言値（before）と実効値（after）を
    method_keyごとに残す——`snapshot_profile`が計算するoverridesを書き込みにしか
    使わず捨てていた、という指摘（controller裁定 fix round 2 finding 2）への対応。
    """
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    # method原本と同じディレクトリに実在する「古い」msp（内部整合性チェック
    # ＝Finding 3のresolve_method_references(METABOLOMICS_REFERENCE_KEYS)を
    # 通すため、宣言値は解決可能である必要がある。ここでの主眼はその値が
    # profileのdependencyとは別物であることを示す差分の記録）。
    _write(root / "method" / "stale.msp", "stale-msp-content")
    dependencies = [_msp_dependency(msp)]
    profile = _build_profile(
        root,
        method_lines="Ion mode: Positive\nMsp file path: stale.msp\n",
        dependencies=dependencies, exe=exe)
    plan = profiles.resolve_profile_inputs(profile, root)

    snapshot = profiles.snapshot_profile(plan, tmp_path / "run1")

    diff_by_key = {d["method_key"]: d for d in snapshot["method_overrides"]}
    assert diff_by_key["Msp file path"]["original_value"] == "stale.msp"
    assert diff_by_key["Msp file path"]["new_value"] == str(msp)


def test_snapshot_profile_records_method_key_diff_when_key_absent_originally(tmp_path):
    """method原本にキー自体が無かった場合はoriginal_valueがNone（空文字と区別する）。"""
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dependencies = [_msp_dependency(msp)]
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=dependencies, exe=exe)
    plan = profiles.resolve_profile_inputs(profile, root)

    snapshot = profiles.snapshot_profile(plan, tmp_path / "run1")

    diff_by_key = {d["method_key"]: d for d in snapshot["method_overrides"]}
    assert diff_by_key["Msp file path"]["original_value"] is None
    assert diff_by_key["Msp file path"]["new_value"] == str(msp)


def test_snapshot_profile_does_not_modify_originals(tmp_path):
    """原本不変: method原本・依存原本は snapshot_profile 後もバイト同一。"""
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dependencies = [_msp_dependency(msp)]
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=dependencies, exe=exe)
    method_path = Path(profile["processing"]["method_path"])
    before_method = method_path.read_bytes()
    before_msp = msp.read_bytes()

    plan = profiles.resolve_profile_inputs(profile, root)
    profiles.snapshot_profile(plan, tmp_path / "run1")

    assert method_path.read_bytes() == before_method
    assert msp.read_bytes() == before_msp


def test_snapshot_profile_plan_identity_hash_is_stable_across_output_roots(tmp_path):
    """別出力先の計画hash同一: run_dirが変わってもplan_identity_hashは変わらない。

    実効コピーの中身（依存への絶対参照）はrun_dir自身の場所には依存しないため
    バイト内容は同じになりうる——ここで守るべき不変条件は「run_dirが変わっても
    計画の同一性(plan_identity_hash)は変わらない」ことそのもの（spec §6.1）で
    あり、実行証跡ハッシュがrun_dirごとに違うことまでは要求しない。run_dir・
    実効コピーの絶対path自体はrun_dirごとに異なることを別途確認する。
    """
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    msp = _write(root / "library" / "lib1.msp", "msp-content")
    dependencies = [_msp_dependency(msp)]
    profile = _build_profile(root, method_lines="Ion mode: Positive\n",
                              dependencies=dependencies, exe=exe)
    plan = profiles.resolve_profile_inputs(profile, root)

    snapshot_a = profiles.snapshot_profile(plan, tmp_path / "run_a")
    snapshot_b = profiles.snapshot_profile(plan, tmp_path / "run_b")

    assert snapshot_a["plan_identity_hash"] == snapshot_b["plan_identity_hash"]
    assert snapshot_a["plan_identity_hash"] == profiles.canonical_hash(plan)
    assert snapshot_a["run_dir"] != snapshot_b["run_dir"]
    effective_a = Path(snapshot_a["run_dir"]) / snapshot_a["effective_method_relative_path"]
    effective_b = Path(snapshot_b["run_dir"]) / snapshot_b["effective_method_relative_path"]
    assert effective_a != effective_b


# ---------- Finding 3 (fix round 2): REFERENCE_KEYSの拡張はmetabolomics側限定 ----------

def test_inspect_inputs_ignores_unresolvable_msp_file_path_like_before(tmp_path, monkeypatch):
    """v1 lipidomics経路の既存挙動を固定するテスト（A01「既存lipidomics実行が
    変わらない」）。`method_file.REFERENCE_KEYS`は`Lbm file path`だけの既定値の
    ままなので、lipidomicsのメソッドファイルが（本タスクと無関係な理由で）
    解決不能な`Msp file path`を宣言していても、`inspect_inputs`は今までどおり
    成功する——`METHOD_REFERENCE_UNRESOLVED`にはならない。

    2026-10-05 から、相対宣言は実在を問わず原本基準の絶対パスとして`overrides`
    に載る（LC-MS の Console は相対値を自分の cwd＝run_dir 基準で読むので、
    相対のまま写すと原本の隣にあっても届かない）。実在しない宣言では Console
    はどちらにしても黙って飛ばすので、実行結果は変わらない。
    """
    from lipidmix.core import session_state
    session_state.session = session_state.AnalysisSession()
    from lipidmix.pipeline.inputs import inspect_inputs

    # このテストが見たいのは「exe と同じフォルダの .lbm2 が使われる」経路。
    # `resolve_lbm` の優先順位は 明示引数 → メソッド宣言 → ビルド生成物 →
    # `MSDIAL_LBM` → exe と同じフォルダ なので、環境変数が設定されていると
    # 下の `lib.lbm2` へ到達する前に打ち切られる。実機の `MSDIAL_LBM` が
    # 実在しないファイルを指していると `LBM_NOT_FOUND` で落ちるため、
    # 周囲の環境から切り離す（同じ手当てが test_console_plan_method.py ほかにある）。
    monkeypatch.delenv("MSDIAL_LBM", raising=False)

    dataset_root = tmp_path / "dataset"
    dataset_root.mkdir()
    (dataset_root / "a.wiff").touch()

    exe_dir = tmp_path / "msdial_app"
    exe_dir.mkdir()
    exe = exe_dir / "MSDIALCUI.exe"
    exe.touch()
    lbm = exe_dir / "lib.lbm2"
    lbm.touch()
    monkeypatch.setattr("lipidmix.console.runner.is_console_exe", lambda *a, **k: True)

    # `lib1.msp`は実在しない（解決不能な宣言）。REFERENCE_KEYS拡張前と同じく
    # 無視されるはず。
    method = dataset_root / "params.txt"
    method.write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\nMsp file path: lib1.msp\n",
        encoding="ascii")

    request = {"method_file": str(method), "polarity": "negative"}
    plan = inspect_inputs(dataset_root, request, exe_path=exe)

    # 例外を投げずに計画が返ること自体が「止めない」ことの証拠
    # （inspect_inputsの戻り値に"status"キーは無い——それはMCPツール層
    # (`console_plan`)が付けるもので、ここでは戻り値自体で確認する）。
    assert Path(plan["method"]["overrides"][method_file_mod.MSP_KEY]) == \
        (dataset_root / "lib1.msp").resolve()


def test_resolve_profile_inputs_rejects_unresolvable_declared_reference_in_method_file(tmp_path):
    """metabolomics/adapter経路限定でREFERENCE_KEYSの拡張集合
    （`METABOLOMICS_REFERENCE_KEYS`）を使う。method原本が（profileのdependency
    宣言とは別に）`Text DB file path`を解決不能な値で宣言していれば、
    `resolve_profile_inputs`はそれを見逃さず`METHOD_REFERENCE_UNRESOLVED`にする
    ——profile自身のdependency検証（`_resolve_dependency`）が通っていても、
    method原本自体の内部整合性は別途確認する。
    """
    root = _nested_root(tmp_path)
    exe = _write_fake_exe(root)
    profile = _build_profile(
        root,
        method_lines=(
            "Ion mode: Positive\n"
            "Text DB file path: does_not_exist.txt\n"
        ),
        dependencies=[], exe=exe)

    with pytest.raises(DomainError, match="METHOD_REFERENCE_UNRESOLVED"):
        profiles.resolve_profile_inputs(profile, root)
