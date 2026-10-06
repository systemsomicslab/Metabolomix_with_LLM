# tests/test_execution_record.py
import json

import pytest

from tests.pipeline_fixtures import execution_record
from metabolomix.core.atomic_io import DomainError, atomic_write_json, canonical_hash
from metabolomix.console.execution import validate_execution_record


# ---------- validate_execution_record: 正常系 ----------

def test_valid_record_round_trips_as_dict():
    result = validate_execution_record(execution_record())
    assert result["execution_id"] == "exec-1"
    assert result["exit_code"] == 0


def test_null_exit_code_allowed_when_not_exited():
    """回収不能なexit_codeはnullが正しい。0を補ってはいけない（spec 5.1）。"""
    result = validate_execution_record(
        execution_record(termination="timeout", exit_code=None)
    )
    assert result["exit_code"] is None


@pytest.mark.parametrize("termination", ["worker_lost", "launch_failed"])
def test_unknown_pid_allowed_for_worker_lost_and_launch_failed(termination):
    """起動失敗・監視喪失では、起動情報（hash等）とは別にpidの不明を許可する。"""
    result = validate_execution_record(
        execution_record(termination=termination, exit_code=None,
                         pid=None, process_identity=None)
    )
    assert result["pid"] is None


# ---------- validate_execution_record: 異常系 ----------

def test_exit_code_bool_is_not_an_integer_receipt():
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(exit_code=True))


def test_exited_without_exit_code_is_rejected():
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(exit_code=None))


def test_missing_exit_code_key_raises_domain_error_not_key_error():
    """exit_codeキー自体が無い場合もKeyErrorではなくDomainErrorであるべき。"""
    record = execution_record()
    del record["exit_code"]
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(record)


@pytest.mark.parametrize("field", ["started_at", "ended_at"])
def test_tz_aware_non_utc_offset_is_rejected(field):
    """tzを持っていてもUTC以外のオフセットは拒否する（spec 5.1: UTC日時）。"""
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(**{field: "2026-09-05T09:00:00+09:00"}))


def test_wrong_schema_is_rejected():
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(schema="console-execution.v0"))


def test_unknown_termination_is_rejected():
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(termination="crashed"))


@pytest.mark.parametrize("timeout_s", [0, -1, 1.5, True])
def test_non_positive_or_non_int_timeout_is_rejected(timeout_s):
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(timeout_s=timeout_s))


def test_ended_before_started_is_rejected():
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(
            started_at="2026-09-05T00:00:10Z", ended_at="2026-09-05T00:00:00Z",
        ))


@pytest.mark.parametrize("field", ["execution_id", "job_id"])
def test_empty_id_is_rejected(field):
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(**{field: ""}))


@pytest.mark.parametrize("field", ["command_sha256", "method_sha256", "exe_sha256"])
def test_hash_not_64_hex_chars_is_rejected(field):
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(**{field: "not-a-hash"}))


def test_pid_and_process_identity_mismatch_is_rejected():
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(
            pid=1234, process_identity={"pid": 9999, "creation_time": 100}
        ))


def test_exited_with_null_pid_is_rejected():
    """起動情報のpid省略が許されるのはworker_lost/launch_failedだけ。"""
    with pytest.raises(DomainError, match="EXECUTION_RECORD_INVALID"):
        validate_execution_record(execution_record(pid=None, process_identity=None))


# ---------- canonical_hash ----------

def test_canonical_hash_is_deterministic_regardless_of_key_order():
    a = canonical_hash({"x": 1, "y": 2})
    b = canonical_hash({"y": 2, "x": 1})
    assert a == b
    assert len(a) == 64


def test_canonical_hash_rejects_nan():
    with pytest.raises(ValueError):
        canonical_hash({"x": float("nan")})


# ---------- atomic_write_json ----------

def test_atomic_write_json_creates_readable_file(tmp_path):
    path = tmp_path / "state.json"
    atomic_write_json(path, {"status": "completed"})
    assert json.loads(path.read_text(encoding="utf-8")) == {"status": "completed"}


def test_atomic_write_json_rejects_nan_and_leaves_no_file(tmp_path):
    path = tmp_path / "state.json"
    with pytest.raises(ValueError):
        atomic_write_json(path, {"x": float("nan")})
    assert not path.exists()


def test_failed_replace_keeps_old_json(tmp_path, monkeypatch):
    path = tmp_path / "state.json"
    path.write_text('{"status":"running"}', encoding="utf-8")

    def fail_replace(*args):
        raise OSError("locked")

    monkeypatch.setattr("metabolomix.core.atomic_io.os.replace", fail_replace)
    with pytest.raises(OSError):
        atomic_write_json(path, {"status": "completed"})
    assert json.loads(path.read_text(encoding="utf-8"))["status"] == "running"


def test_failed_replace_leaves_no_stray_temp_file(tmp_path, monkeypatch):
    path = tmp_path / "state.json"
    path.write_text('{"status":"running"}', encoding="utf-8")

    def fail_replace(*args):
        raise OSError("locked")

    monkeypatch.setattr("metabolomix.core.atomic_io.os.replace", fail_replace)
    with pytest.raises(OSError):
        atomic_write_json(path, {"status": "completed"})
    leftovers = [p for p in tmp_path.iterdir() if p.name != "state.json"]
    assert leftovers == []
