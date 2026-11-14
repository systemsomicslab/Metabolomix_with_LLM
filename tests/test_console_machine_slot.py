"""機械全体で MS-DIAL Console を 1 本にする排他と、起こす前の commit の空き（B8-T66）。"""
from __future__ import annotations

import threading

import pytest

from metabolomix.console import machine_slot
from metabolomix.core import process_control
from tests.config_fixtures import use_config

_GB = 1024 ** 3


def _commit(monkeypatch, available_gb: float, limit_gb: float = 64) -> None:
    monkeypatch.setattr(process_control, "commit_memory_status", lambda: {
        "available_commit_bytes": int(available_gb * _GB),
        "commit_limit_bytes": int(limit_gb * _GB),
        "available_physical_bytes": int(available_gb * _GB)})


def test_lock_lives_next_to_the_config_file(tmp_path, monkeypatch):
    config = use_config(monkeypatch, tmp_path, {})
    assert machine_slot.console_lock_path() == config.parent / "data" / ".msdial-console.lock"


def test_low_commit_is_a_warning_with_the_default_threshold(monkeypatch):
    _commit(monkeypatch, 3.0)
    measured, warning = machine_slot.check_commit_memory()
    assert measured == {"available_commit_gb": 3.0, "commit_limit_gb": 64.0,
                        "min_free_commit_gb": machine_slot.DEFAULT_MIN_FREE_COMMIT_GB}
    assert warning["code"] == "LOW_COMMIT_MEMORY"
    assert warning["details"] == measured
    assert "[msdial] min_free_commit_gb" in warning["message"]


def test_enough_commit_is_measured_without_a_warning(monkeypatch):
    _commit(monkeypatch, 40.0)
    measured, warning = machine_slot.check_commit_memory()
    assert measured["available_commit_gb"] == 40.0
    assert warning is None


def test_threshold_comes_from_the_config(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"msdial.min_free_commit_gb": 32})
    _commit(monkeypatch, 20.0)
    measured, warning = machine_slot.check_commit_memory()
    assert measured["min_free_commit_gb"] == 32
    assert warning is not None


def test_zero_threshold_skips_the_measurement(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"msdial.min_free_commit_gb": 0})
    _commit(monkeypatch, 0.5)
    assert machine_slot.check_commit_memory() == (None, None)


def test_unmeasurable_commit_is_silent(monkeypatch):
    monkeypatch.setattr(process_control, "commit_memory_status", lambda: None)
    assert machine_slot.check_commit_memory() == (None, None)


def test_slot_records_its_owner_and_a_second_taker_waits(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {})
    inside = threading.Event()
    release = threading.Event()

    def first():
        with machine_slot.console_slot(owner={"pipeline_id": "p1", "unit": "negative"}):
            inside.set()
            release.wait(10)

    thread = threading.Thread(target=first)
    thread.start()
    try:
        assert inside.wait(10)
        assert machine_slot.read_owner()["pipeline_id"] == "p1"
        notices = []

        def on_wait():
            notices.append(machine_slot.read_owner())
            threading.Timer(0.2, release.set).start()

        with machine_slot.console_slot(owner={"pipeline_id": "p2", "unit": "positive"},
                                       on_wait=on_wait, poll_interval=0.05) as waited:
            assert waited > 0
            assert machine_slot.read_owner()["pipeline_id"] == "p2"
        assert [n["pipeline_id"] for n in notices] == ["p1"]
    finally:
        release.set()
        thread.join(10)


def test_slot_wait_can_be_abandoned(tmp_path, monkeypatch):
    from metabolomix.core.atomic_io import DomainError

    use_config(monkeypatch, tmp_path, {})
    with process_control.file_lock(machine_slot.console_lock_path()):
        with pytest.raises(DomainError) as caught:
            with machine_slot.console_slot(owner={}, should_abort=lambda: True,
                                           poll_interval=0.01):
                pass
    assert caught.value.code == "LOCK_WAIT_ABORTED"
