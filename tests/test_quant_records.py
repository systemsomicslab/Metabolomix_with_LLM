"""定量の工程の記録（metabolomix.quant.records）。fixture は tmp に合成する。"""
from __future__ import annotations

import errno
import json
import sys

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.quant import records


def _write(folder, step="normalize", **kw):
    args = dict(inputs={"x/a.mzTab": "aaa"}, parent_ids=[], settings={"rules_version": "1"},
                summary={"n": 1})
    args.update(kw)
    return records.write_record(folder, step, **args)


def test_write_then_read_roundtrip(tmp_path):
    r1 = _write(tmp_path)
    r2 = _write(tmp_path, "merge", parent_ids=[r1["record_id"]])
    got, warnings = records.read_records(tmp_path)
    assert warnings == []
    assert [r["record_id"] for r in got] == [r1["record_id"], r2["record_id"]]
    assert r1["schema"] == records.RECORD_SCHEMA
    assert r1["record_id"].startswith("q_")
    assert (tmp_path / "reports" / "quant" / f"{r1['record_id']}.json").is_file()
    for key in ("step", "created_at", "inputs", "parent_ids", "settings", "summary", "skipped"):
        assert key in r1


def test_never_overwrites(tmp_path, monkeypatch):
    monkeypatch.setattr(records, "_new_record_id", lambda: "q_fixed")
    _write(tmp_path)
    with pytest.raises(DomainError):
        _write(tmp_path)
    got, _ = records.read_records(tmp_path)
    assert len(got) == 1


def test_stale_when_source_sha_changes(tmp_path):
    r = _write(tmp_path, inputs={"x/a.mzTab": "aaa"})
    assert records.is_stale(r, {"x/a.mzTab": "aaa"}, set()) is None
    assert records.is_stale(r, {"x/a.mzTab": "bbb"}, set()) == "input_changed:a.mzTab"
    assert records.is_stale(r, {}, set()) == "input_changed:a.mzTab"


def test_stale_propagates_to_children(tmp_path):
    parent = _write(tmp_path)
    child = _write(tmp_path, "merge", parent_ids=[parent["record_id"]])
    pid = parent["record_id"]
    assert records.is_stale(child, {"x/a.mzTab": "aaa"}, {pid}) is None
    assert records.is_stale(child, {"x/a.mzTab": "aaa"}, set(), stale_ids={pid}) == f"parent_stale:{pid}"
    assert records.is_stale(child, {"x/a.mzTab": "aaa"}, set()) == f"parent_missing:{pid}"


def test_merged_folder_common_parent(tmp_path):
    (tmp_path / "POS").mkdir()
    (tmp_path / "NEG").mkdir()
    assert records.merged_record_folder(tmp_path / "POS", tmp_path / "NEG") == tmp_path.resolve()


def test_merged_folder_refuses_nested(tmp_path):
    (tmp_path / "POS").mkdir()
    for a, b in ((tmp_path, tmp_path / "POS"), (tmp_path / "POS", tmp_path),
                 (tmp_path, tmp_path)):
        with pytest.raises(DomainError) as e:
            records.merged_record_folder(a, b)
        assert e.value.code == "QUANT_RECORD_LOCATION_REQUIRED"


@pytest.mark.skipif(sys.platform != "win32", reason="ドライブがあるのは Windows だけ")
def test_merged_folder_refuses_cross_drive():
    with pytest.raises(DomainError) as e:
        records.merged_record_folder(__import__("pathlib").Path("C:/a/POS"),
                                     __import__("pathlib").Path("D:/a/NEG"))
    assert e.value.code == "QUANT_RECORD_LOCATION_REQUIRED"


def test_corrupt_record_is_listed_as_warning(tmp_path):
    ok = _write(tmp_path)
    d = records.record_dir(tmp_path)
    (d / "broken.json").write_text("{not json", encoding="utf-8")
    (d / "other.json").write_text(json.dumps({"schema": "x.v9"}), encoding="utf-8")
    got, warnings = records.read_records(tmp_path)
    assert [r["record_id"] for r in got] == [ok["record_id"]]
    assert any("broken.json" in w for w in warnings)
    assert any("other.json" in w for w in warnings)


def test_read_records_missing_folder(tmp_path):
    assert records.read_records(tmp_path) == ([], [])


@pytest.mark.parametrize("err", [errno.EPERM, getattr(errno, "ENOTSUP", errno.EOPNOTSUPP)])
def test_link_unsupported_falls_back_to_exclusive_create(tmp_path, monkeypatch, err):
    """exFAT / FAT の外付けドライブは os.link を持たない。上書きしない排他作成に落とす。"""
    def no_link(src, dst):
        raise OSError(err, "hard links not supported")

    monkeypatch.setattr(records.os, "link", no_link)
    r = _write(tmp_path)
    directory = records.record_dir(tmp_path)
    path = directory / f"{r['record_id']}.json"
    assert json.loads(path.read_text(encoding="utf-8"))["record_id"] == r["record_id"]
    assert [p.name for p in directory.iterdir()] == [path.name]  # 一時ファイルを残さない
    with pytest.raises(DomainError) as exc:
        records._create_exclusive(path, {"x": 1})
    assert exc.value.code == "QUANT_RECORD_EXISTS"
    assert json.loads(path.read_text(encoding="utf-8"))["record_id"] == r["record_id"]
