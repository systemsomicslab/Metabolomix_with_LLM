"""Profile resolution must agree across file requests and worker processes."""
import os
import json

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import request as request_mod, service
from tests.test_metabolomics_engine import _profile


def test_request_file_validation_purpose_is_used_when_loading_profile(tmp_path):
    profile_path = _profile(tmp_path)
    (tmp_path / "analysis-request.json").write_text(json.dumps({
        "schema": "pipeline-request.v2", "profile_file": str(profile_path),
        "execution_purpose": "validation",
    }), encoding="utf-8")
    arguments = service._profile_arguments(tmp_path, None)
    resolved = request_mod.resolve_request(tmp_path, None, **arguments)
    assert resolved["execution_purpose"] == "validation"
    assert arguments["profile"]["validation"]["status"] == "draft"


def test_relative_profile_is_resolved_against_dataset_root(tmp_path, monkeypatch):
    _profile(tmp_path)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    explicit = {"schema": "pipeline-request.v2", "profile_file": "profile.json",
                "execution_purpose": "validation"}
    arguments = service._profile_arguments(tmp_path, explicit)
    resolved = request_mod.resolve_request(tmp_path, explicit, **arguments)
    assert resolved["profile_file"] == str(tmp_path / "profile.json")
    assert resolved["value_sources"]["profile_file"] == "explicit"
    assert explicit["profile_file"] == "profile.json"


@pytest.mark.parametrize("value", [None, "", 3, []])
def test_invalid_explicit_profile_does_not_fall_back_to_request_file(tmp_path, value):
    profile_path = _profile(tmp_path)
    (tmp_path / "analysis-request.json").write_text(json.dumps({
        "schema": "pipeline-request.v2", "profile_file": str(profile_path),
        "execution_purpose": "validation",
    }), encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        service._profile_arguments(tmp_path, {"schema": "pipeline-request.v2",
                                              "profile_file": value})
    assert exc.value.code == "PIPELINE_REQUEST_INVALID"


def test_explicit_routine_overrides_file_validation(tmp_path):
    profile_path = _profile(tmp_path)
    (tmp_path / "analysis-request.json").write_text(json.dumps({
        "schema": "pipeline-request.v2", "profile_file": str(profile_path),
        "execution_purpose": "validation",
    }), encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        service._profile_arguments(tmp_path, {"execution_purpose": "routine"})
    assert exc.value.code == "PROFILE_NOT_VALIDATED"


# 「公開入口が v2 要求を v1 の推定経路へ落とさない」ことの検証は
# `tests/test_pipeline_inputs_v2.py` へ移した。以前はここで
# `PIPELINE_V2_UPSTREAM_UNAVAILABLE` による停止そのものを縛っていたが、停止は
# 「v1 へ落ちない」ための手段であって目的ではない。経路が繋がった今は、
# `inspect_inputs` と環境設定の実行体が v2 要求では一度も呼ばれないことを
# 直接縛っている（同じ主張を2箇所に置かない）。


def test_snapshot_uses_profile_directory_and_persists_raw_hashes(tmp_path):
    from metabolomix.pipeline import metabolomics_handlers
    from tests.test_lcms_profile_inputs import _build_profile, _write_fake_exe, _write, _msp_dependency

    profile_root = tmp_path / "configuration"
    profile_root.mkdir()
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    _write(raw_root / "sample.mzML", "synthetic raw")
    exe = _write_fake_exe(profile_root)
    library = _write(profile_root / "lib.msp", "synthetic library")
    profile = _build_profile(profile_root, method_lines="Ion mode: Positive\n",
                             dependencies=[_msp_dependency(library)], exe=exe)
    profile["processing"]["method_path"] = "method/params.txt"
    path = profile_root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    context = {"request": {"profile_file": str(path), "execution_purpose": "validation"},
               "identity": {"source_root": str(raw_root)}, "runtime": {},
               "pipeline_root": tmp_path / "run", "pipeline_id": "boundary"}
    outcome = metabolomics_handlers.snapshot_profile_outcome(context, {"status": "succeeded"})
    plan = context["runtime"]["profile_plan"]
    assert plan["raw_files"]["files"][0]["relative_path"] == "sample.mzML"
    ref = next(r for r in outcome["result_refs"] if r["output_name"] == "execution_manifest")
    # Read the actual persisted artifact, not a mocked persistence callback.
    artifact = json.loads((context["pipeline_root"] / ref["relative_path"]).read_text(encoding="utf-8"))
    assert artifact["data"]["raw"]["files"][0]["sha256"]


# ---------- 再解決で raw の内容が変わっていたら止める（2026-09-17 Stage B の #3）----------

def _profile_context(tmp_path):
    """profile / raw / exe を実ファイルで用意し、snapshot_profile_outcome 用の
    context を返す。"""
    from tests.test_lcms_profile_inputs import _build_profile, _write_fake_exe, _write, _msp_dependency

    profile_root = tmp_path / "configuration"
    profile_root.mkdir()
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    _write(raw_root / "sample.mzML", "synthetic raw")
    exe = _write_fake_exe(profile_root)
    library = _write(profile_root / "lib.msp", "synthetic library")
    profile = _build_profile(profile_root, method_lines="Ion mode: Positive\n",
                             dependencies=[_msp_dependency(library)], exe=exe)
    profile["processing"]["method_path"] = "method/params.txt"
    path = profile_root / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    return {"request": {"profile_file": str(path), "execution_purpose": "validation"},
            "identity": {"source_root": str(raw_root)}, "runtime": {},
            "pipeline_root": tmp_path / "run", "pipeline_id": "rawcheck"}, raw_root


def test_raw_content_change_between_runs_stops_with_input_changed(tmp_path):
    """再解決した raw hash が、前回固定した profile snapshot と食い違えば止める。

    v2 は既に全 raw を hash して snapshot へ保存している。突き合わせないと、
    入力が差し替わっても気付かないまま解析が進む。stat（size/mtime）だけを
    見る v1 の `verify_inputs` は仕様上 raw の内容 hash を検証しない契約なので、
    ここで v2 が自分の持つ hash を使う。
    """
    from metabolomix.pipeline import metabolomics_handlers
    from tests.test_lcms_profile_inputs import _write

    context, raw_root = _profile_context(tmp_path)
    first = metabolomics_handlers.snapshot_profile_outcome(context, {"status": "succeeded"})

    # 同じ長さのまま内容だけ差し替える（size と mtime では気付けない改変）。
    target = raw_root / "sample.mzML"
    stat = target.stat()
    _write(target, "synthetic RAW")          # 同じ 13 文字
    os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert target.stat().st_size == stat.st_size
    assert target.stat().st_mtime_ns == stat.st_mtime_ns

    context["runtime"] = {}
    context["results"] = first["result_refs"]
    with pytest.raises(DomainError) as exc:
        metabolomics_handlers.snapshot_profile_outcome(context, {"status": "succeeded"})
    assert exc.value.code == "INPUT_CHANGED"
    assert "sample.mzML" in json.dumps(exc.value.details, ensure_ascii=False)


def test_an_unchanged_rerun_passes_the_raw_check(tmp_path):
    """内容が同じ再実行は素通りする（照合そのものが誤検出しないこと）。"""
    from metabolomix.pipeline import metabolomics_handlers

    context, _raw_root = _profile_context(tmp_path)
    first = metabolomics_handlers.snapshot_profile_outcome(context, {"status": "succeeded"})

    context["runtime"] = {}
    context["results"] = first["result_refs"]
    again = metabolomics_handlers.snapshot_profile_outcome(context, {"status": "succeeded"})
    assert again["status"] == "succeeded"
