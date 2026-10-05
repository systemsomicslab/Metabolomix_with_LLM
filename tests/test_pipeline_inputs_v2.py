"""v2 の入力配置計画は profile だけを情報源にする。

v1 の `inspect_inputs` はフォルダから method を推定し LBM を必須にし、実行体を
環境設定から取る。v2 でそこへ落ちると、profile が固定したのとは別の条件で
Console が走る——`PIPELINE_V2_UPSTREAM_UNAVAILABLE` で公開経路を塞いでいたのは
これを避けるためだった。ここで縛るのは「落ちないこと」そのもの。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from lipidmix.console.profile_schema import validate_profile
from lipidmix.core.atomic_io import DomainError
from lipidmix.pipeline import inputs as inputs_mod
from tests.metabolomics_fixtures import write_profile


def _setup(tmp_path, monkeypatch) -> tuple[Path, dict, dict]:
    # fixture の実行体は placeholder なので、Console 判定（実際に --help を
    # 起動する）は既存テストと同じ流儀で差し替える。
    monkeypatch.setattr("lipidmix.console.runner.is_console_exe", lambda *a, **k: True)
    source_root = tmp_path / "source"
    source_root.mkdir()
    for name in ("A", "B"):
        (source_root / f"{name}.wiff").write_bytes(b"synthetic raw")
    profile_path = write_profile(source_root / "profile.json")
    profile = validate_profile(json.loads(profile_path.read_text(encoding="utf-8")))
    # draft profile なので purpose は validation。省略すると既定の routine になり、
    # 証明書が要る経路（R2 未接続）へ入ってしまう。
    request = {"schema": "pipeline-request.v2", "profile_file": str(profile_path),
               "execution_purpose": "validation"}
    return source_root, request, profile


def test_plan_has_the_same_shape_as_the_v1_plan(tmp_path, monkeypatch):
    source_root, request, profile = _setup(tmp_path, monkeypatch)
    plan = inputs_mod.plan_from_profile(source_root, request, profile)
    assert set(plan) == {"source_root", "selected_format", "entries", "raw_stat",
                         "companions", "method", "lbm", "exe", "polarity", "unverified"}
    assert plan["selected_format"] == "wiff"
    assert len(plan["raw_stat"]) == 2


def test_method_and_executable_come_from_the_profile(tmp_path, monkeypatch):
    source_root, request, profile = _setup(tmp_path, monkeypatch)
    plan = inputs_mod.plan_from_profile(source_root, request, profile)
    assert plan["method"]["sha256"] == profile["processing"]["method_sha256"]
    assert plan["exe"]["sha256"] == profile["software"]["executable_sha256"]
    assert plan["polarity"]["source"] == "profile"
    assert plan["polarity"]["value"] == profile["acquisition"]["polarity"]


def test_every_declared_dependency_becomes_a_method_override(tmp_path, monkeypatch):
    source_root, request, profile = _setup(tmp_path, monkeypatch)
    plan = inputs_mod.plan_from_profile(source_root, request, profile)
    for dependency in profile["processing"]["dependencies"]:
        assert dependency["method_key"] in plan["method"]["overrides"]


def test_a_changed_dependency_stops_before_any_plan_is_returned(tmp_path, monkeypatch):
    source_root, request, profile = _setup(tmp_path, monkeypatch)
    target = source_root / profile["processing"]["dependencies"][0]["path"]
    target.write_text("tampered", encoding="utf-8")
    with pytest.raises(DomainError) as excinfo:
        inputs_mod.plan_from_profile(source_root, request, profile)
    assert excinfo.value.code == "INPUT_CHANGED"


# ---------- 受付が v1 経路へ落ちないこと ----------

def test_v2_request_does_not_go_through_the_v1_inspection(tmp_path, monkeypatch):
    from lipidmix.pipeline import service

    source_root, request, _profile = _setup(tmp_path, monkeypatch)

    def _fail(*args, **kwargs):
        raise AssertionError("v2 要求が v1 の inspect_inputs へ落ちた")

    monkeypatch.setattr(inputs_mod, "inspect_inputs", _fail)
    result = service.plan_pipeline(source_root, request)
    assert result["status"] in {"planned", "needs_input"}, result


def test_v2_request_does_not_touch_the_environment_executable(tmp_path, monkeypatch):
    from lipidmix.pipeline import service

    source_root, request, _profile = _setup(tmp_path, monkeypatch)

    def _fail():
        raise AssertionError("v2 要求が環境設定の実行体を参照した")

    monkeypatch.setattr(service, "_resolve_exe_path", _fail)
    service.plan_pipeline(source_root, request)


def test_a_v2_request_without_a_profile_stops_instead_of_falling_back(tmp_path, monkeypatch):
    from lipidmix.pipeline import service

    source_root, request, _profile = _setup(tmp_path, monkeypatch)
    Path(request["profile_file"]).unlink()

    def _fail(*args, **kwargs):
        raise AssertionError("profile が無い v2 要求が v1 の inspect_inputs へ落ちた")

    monkeypatch.setattr(inputs_mod, "inspect_inputs", _fail)
    with pytest.raises(DomainError) as excinfo:
        service.plan_pipeline(source_root, request)
    assert excinfo.value.code != "PIPELINE_V2_UPSTREAM_UNAVAILABLE"
    assert excinfo.value.code in {"PIPELINE_REQUEST_INVALID", "PROFILE_NOT_FOUND"}


# ---------- Console 起動（job v3） ----------

def _upstream_context(tmp_path, monkeypatch):
    """受付 → prepare_inputs まで通し、`execute_console` の context を返す。

    engine を回さず handler を直接呼ぶのは、ここで見たいのが「job に何を書くか」
    だけだから。prepare_inputs が固定した成果物 ref を record へ載せて渡す
    （engine が commit_stage_outcome でやることを、この試験の範囲だけ手で行う）。
    """
    from lipidmix.pipeline import engine as engine_mod
    from lipidmix.pipeline import request as request_mod
    from lipidmix.pipeline import service, store

    source_root, request, _profile = _setup(tmp_path, monkeypatch)
    receipt = service.plan_pipeline(source_root, request)
    pipeline_root = Path(receipt["pipeline_path"])
    resolved = request_mod.resolve_request(
        source_root, request, **service._profile_arguments(source_root, request))
    record = store.load_run(pipeline_root)

    prepared = service._handle_prepare_inputs_v2(engine_mod.make_context(
        record, record["stages"]["prepare_inputs"], {}, resolved))
    assert prepared["status"] == "succeeded", prepared
    record["results"] = list(prepared["result_refs"])
    context = engine_mod.make_context(
        record, record["stages"]["execute_console"], {}, resolved)
    return context, pipeline_root


def test_v2_upstream_writes_a_v3_job_with_the_profile_snapshot(tmp_path, monkeypatch):
    from lipidmix.console import execution as console_execution
    from lipidmix.pipeline import service

    context, pipeline_root = _upstream_context(tmp_path, monkeypatch)
    monkeypatch.setattr(console_execution, "supervise",
                        lambda job_path, cancel_path=None: {
                            "execution_id": "exec-test", "termination": "exited",
                            "exit_code": 0})

    service._handle_upstream(context)

    job_path = next(pipeline_root.rglob("analysis-job.json"))
    data = json.loads(job_path.read_text(encoding="utf-8"))
    assert data["schema"] == "analysis-job.v3"
    # dataclass の field は profile_snapshot、wire 上のキーは profile。
    assert data["profile"]["profile_id"] == "synthetic-metabolomics"
    assert data["profile"]["dependencies"], "依存manifestがsnapshotに無い"
    assert data["project"]["omics"] == "metabolomics"
    assert data["project"]["measure"] == "peak_height"


def test_upstream_reads_the_snapshot_from_persisted_results_not_from_runtime(
        tmp_path, monkeypatch):
    """再開で prepare_inputs が skip されると runtime は空になる。

    そのとき snapshot を runtime から取っていると、再開後に書く job だけが
    snapshot を失い、v3 のはずの job が黙って v2 になる。成果物が無いなら
    「無い」と言って止まること（前の工程を先に通す）を縛る。
    """
    from lipidmix.console import execution as console_execution
    from lipidmix.pipeline import service

    context, _pipeline_root = _upstream_context(tmp_path, monkeypatch)
    monkeypatch.setattr(console_execution, "supervise",
                        lambda job_path, cancel_path=None: {
                            "execution_id": "exec-test", "termination": "exited",
                            "exit_code": 0})
    context["results"] = []          # 成果物が record に無い状態
    context["runtime"] = {}          # runtime も空（再開直後と同じ）

    outcome = service._handle_upstream(context)
    assert outcome["status"] == "needs_input"
    assert outcome["error"]["code"] == "PROFILE_SNAPSHOT_MISSING"


# ---------- profile 外の相対宣言も絶対パスで渡す ----------

def _declare_extra_relative_text_db(source_root: Path, profile: dict) -> Path:
    """method 原本に profile が依存として宣言していない相対キーを足す。

    LC-MS の Console は宣言パスを解決せず cwd 基準で読むので、profile の依存
    だけを上書きすると、この行は相対のまま Console へ届いて黙って無視される。
    """
    import hashlib
    method = source_root / profile["processing"]["method_path"]
    (source_root / "lib" / "t.txt").write_text("synthetic\n", encoding="utf-8")
    method.write_text(method.read_text(encoding="utf-8")
                      + "Text DB file path: ../lib/t.txt\n",
                      encoding="utf-8", newline="\n")
    profile["processing"]["method_sha256"] = hashlib.sha256(method.read_bytes()).hexdigest()
    return (source_root / "lib" / "t.txt").resolve()


def test_a_relative_method_declaration_outside_the_profile_becomes_absolute(tmp_path, monkeypatch):
    source_root, request, profile = _setup(tmp_path, monkeypatch)
    expected = _declare_extra_relative_text_db(source_root, profile)
    plan = inputs_mod.plan_from_profile(source_root, request, profile)
    assert Path(plan["method"]["overrides"]["Text DB file path"]) == expected


def test_the_profile_snapshot_writes_the_same_absolute_declaration(tmp_path, monkeypatch):
    from lipidmix.console import method_file as method_file_mod
    from lipidmix.console import profiles as profiles_mod
    source_root, _request, profile = _setup(tmp_path, monkeypatch)
    expected = _declare_extra_relative_text_db(source_root, profile)
    plan = profiles_mod.resolve_profile_inputs(profile, source_root, raw_root=source_root)
    snapshot = profiles_mod.snapshot_profile(plan, tmp_path / "run")
    effective = tmp_path / "run" / snapshot["effective_method_relative_path"]
    keys = method_file_mod.read_method_keys(effective)
    assert Path(keys["text db file path"]) == expected
    assert {"method_key": "Text DB file path", "original_value": "../lib/t.txt",
            "new_value": str(expected)} in snapshot["method_overrides"]
