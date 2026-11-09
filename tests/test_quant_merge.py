"""POS と NEG の統合（merge / dataset_merge_polarities）。人工の解析フォルダと偽の規則で。"""
from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pytest

from metabolomix.core import handles, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.quant import derive as qd
from metabolomix.quant import merge as qm
from metabolomix.quant import normalize as qn
from metabolomix.core import user_config
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.quant_fake_rules import write_fake_rules

POLS = ("positive", "negative")


def _point(tmp_path, monkeypatch):
    folder = write_fake_rules(tmp_path / "rules")
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text(f"[quant]\nrules = '{folder}'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))


def _samples(names, group="A"):
    return [{"file_id": i, "name": n, "class_id": group} for i, n in enumerate(names)]


def _load(folder, names, spots, group="A"):
    write_msdial_folder(folder, spots=spots, samples=_samples(names, group))
    return load_msdial_folder(select_alignment(folder))


def _spot(name, heights, ontology="PC", mz=500.0):
    return {"name": name, "ontology": ontology, "heights": heights, "mz": mz}


def _state(ds, unit="nmol", record_id=None):
    ds.quant_state = qd.quant_state(
        normalized=unit is not None, method="internal_standard" if unit else None, unit=unit,
        record_id=record_id or f"q_{id(ds)}", steps=["normalize"] if unit else [],
        skipped=[] if unit else ["normalize"], internal_standard_ids=["0"] if unit else [])
    return ds


def _pair(tmp_path, pos_names=("S1_pos", "S2_pos", "S3_pos"),
          neg_names=("S1_neg", "S2_neg", "S3_neg"),
          root="run", unit="nmol"):
    pos = _load(tmp_path / root / "pos", list(pos_names),
                [_spot("IS-A", [10.0, 20.0, 50.0]), _spot("PC 34:1", [1.0, 2.0, 5.0], mz=600.0)])
    neg = _load(tmp_path / root / "neg", list(neg_names),
                [_spot("IS-B", [30.0, 40.0, 60.0]), _spot("PE 36:2", [3.0, 4.0, 6.0], "PE", 700.0)])
    pos.per_sample_mz = {"1": {n: 600.0 + i for i, n in enumerate(pos_names)}}
    neg.per_sample_mz = {"1": {n: 700.0 + i for i, n in enumerate(neg_names)}}
    pos.feature_tags = {"1": ["pos-tag"]}
    neg.feature_tags = {"1": ["neg-tag"]}
    pos.curation_flags = {"flags": {"1": "keep"}}
    return _state(pos, unit, "q_pos"), _state(neg, unit, "q_neg")


def _merge(pos, neg, **kw):
    kw.setdefault("polarities", POLS)
    return qm.merge(pos, neg, **kw)


def test_merge_by_polarity_token(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    merged, info = _merge(pos, neg, mapping_rows=None)
    assert info["rule"] == "polarity_token_removed"
    assert merged.sample_names == ["S1", "S2", "S3"]
    assert merged.feature_matrix.shape == (4, 3)
    assert pos.sample_names == ["S1_pos", "S2_pos", "S3_pos"]  # 元は不変
    assert pos.quant_state["steps"] == ["normalize"]


def test_neg_columns_follow_pos_order(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path, neg_names=("S2_neg", "S1_neg", "S3_neg"))
    merged, _ = _merge(pos, neg, mapping_rows=None)
    assert merged.sample_names == ["S1", "S2", "S3"]
    row = merged.feature_ids.index("negative:0")
    assert list(merged.feature_matrix[row]) == [40.0, 30.0, 60.0]  # S1 は NEG の 2 列目
    assert list(merged.feature_matrix[merged.feature_ids.index("positive:0")]) == [10.0, 20.0, 50.0]
    assert merged.detected_mask.shape == merged.feature_matrix.shape


def test_merge_stops_on_mismatched_samples(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path, neg_names=("S1_neg", "S2_neg", "S9_neg"))
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_MERGE_SAMPLE_MISMATCH"
    only = exc.value.details["only_in"]
    # 除いた名前の集合も揃わないので規則は名前を変えない（exact_name）。元の名前がそのまま出る
    assert exc.value.details["rule"] == "exact_name"
    assert only["positive"] == ["S1_pos", "S2_pos", "S3_pos"]
    assert only["negative"] == ["S1_neg", "S2_neg", "S9_neg"]
    assert exc.value.details["only_in_ids"] == only
    # 混在（トークンが付く側と付かない側）でも黙って部分集合にしない
    pos, neg = _pair(tmp_path, neg_names=("S1_neg", "S2_neg", "S3"), root="run2")
    merged, _ = _merge(pos, neg, mapping_rows=None)  # 付かない名前は正確一致で S2
    assert merged.sample_names == ["S1", "S2", "S3"]
    pos, neg = _pair(tmp_path, neg_names=("S1_neg", "S2_neg", "x3"), root="run3")
    with pytest.raises(DomainError):
        _merge(pos, neg, mapping_rows=None)


def _rows(pairs):
    rows = []
    for mode, src, sid in pairs:
        rows.append({"mode": mode, "source_sample": src, "sample_id": sid, "group": "A",
                     "sample_type": "sample"})
    return rows


def test_merge_with_mapping_csv(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path, neg_names=("zz", "yy", "ww"))
    rows = _rows([("positive", "S1_pos", "X1"), ("positive", "S2_pos", "X2"),
                  ("positive", "S3_pos", "X3"), ("negative", "zz", "X2"),
                  ("negative", "yy", "X1"), ("negative", "ww", "X3")])
    merged, info = _merge(pos, neg, mapping_rows=rows, mapping_input=("map.csv", "ab" * 32))
    assert info["rule"] == "mapping_csv" and merged.sample_names == ["X1", "X2", "X3"]
    row = merged.feature_ids.index("negative:0")
    assert list(merged.feature_matrix[row]) == [40.0, 30.0, 60.0]  # X1 は yy（NEG の 2 列目）
    rec = json.loads((tmp_path / "run" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["inputs"]["mapping:map.csv"] == "ab" * 32


def test_read_mapping(tmp_path):
    path = tmp_path / "map.csv"
    path.write_text("mode,source_sample,sample_id,group,sample_type\npos,a,S1,g,x\n",
                    encoding="utf-8")
    rows, sha = qm.read_mapping(path)
    assert rows == [{"mode": "positive", "source_sample": "a", "sample_id": "S1", "group": "g",
                     "sample_type": "x"}] and len(sha) == 64
    path.write_text("mode,source_sample\npos,a\n", encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        qm.read_mapping(path)
    assert exc.value.code == "QUANT_BAD_REQUEST"
    with pytest.raises(DomainError) as exc:
        qm.read_mapping(tmp_path / "nope.csv")
    assert exc.value.code == "QUANT_MAPPING_NOT_FOUND"


def test_merge_refuses_same_polarity_and_mixed_units(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    with pytest.raises(DomainError) as exc:
        qm.merge(pos, neg, mapping_rows=None, polarities=("positive", "positive"))
    assert exc.value.code == "QUANT_MERGE_POLARITY_INVALID"
    with pytest.raises(DomainError) as exc:
        qm.merge(pos, neg, mapping_rows=None, polarities=("negative", "positive"))
    assert exc.value.code == "QUANT_MERGE_POLARITY_INVALID"
    _state(neg, "pmol")  # 単位違い
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_MERGE_STATE_MISMATCH"
    _state(neg, None)  # 片方だけ正規化
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_MERGE_STATE_MISMATCH"
    assert not (tmp_path / "run" / "reports").exists()  # 何も書かない


def test_mapping_mismatch_names_source_samples(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path, neg_names=("zz", "yy", "ww"))
    rows = _rows([("positive", "S1_pos", "X1"), ("positive", "S2_pos", "X2"),
                  ("positive", "S3_pos", "X3"), ("negative", "zz", "X1"),
                  ("negative", "yy", "X2"), ("negative", "ww", "X9")])
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=rows)
    d = exc.value.details
    assert d["only_in"] == {"positive": ["S3_pos"], "negative": ["ww"]}
    assert d["only_in_ids"] == {"positive": ["X3"], "negative": ["X9"]}


def test_read_mapping_accepts_mode_aliases(tmp_path):
    path = tmp_path / "map.csv"
    path.write_text("mode,source_sample,sample_id,group,sample_type\n"
                    "POS,a,S1,g,x\nNeg,b,S1,g,x\nPositive,c,S2,g,x\nnegative,d,S2,g,x\n",
                    encoding="utf-8")
    rows, _ = qm.read_mapping(path)
    assert [r["mode"] for r in rows] == ["positive", "negative", "positive", "negative"]


def test_lock_timeout_is_busy(tmp_path, monkeypatch):
    import threading
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    monkeypatch.setattr(qm, "LOCK_TIMEOUT", 0.2)
    held, release = threading.Event(), threading.Event()
    blocker = handles.registry.get(max(hp, hn))

    def hold():
        with blocker.lock:
            held.set()
            release.wait(5)
    t = threading.Thread(target=hold)
    t.start()
    held.wait(5)
    try:
        out = json.loads(_tool()(hp, hn))
    finally:
        release.set()
        t.join()
    assert out["error"]["code"] == "QUANT_HANDLE_BUSY"
    ok = json.loads(_tool()(hp, hn))  # 取ったロックは戻っている
    assert ok["status"] == "merged"


def test_merge_by_exact_name(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path, pos_names=("a", "b", "c"), neg_names=("b", "c", "a"))
    merged, info = _merge(pos, neg, mapping_rows=None)
    assert info["rule"] == "exact_name" and merged.sample_names == ["a", "b", "c"]
    row = merged.feature_ids.index("negative:0")
    assert list(merged.feature_matrix[row]) == [60.0, 30.0, 40.0]  # a は NEG の 3 列目


def test_merge_requires_quant_state(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    neg.quant_state = None
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_NORMALIZE_STEP_REQUIRED"
    assert "dataset_quant_apply" in exc.value.message and "skip" in exc.value.message


def test_merged_ids_and_origin(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    merged, _ = _merge(pos, neg, mapping_rows=None, pos_handle_id="dh_p", neg_handle_id="dh_n")
    assert merged.feature_ids == ["positive:0", "positive:1", "negative:0", "negative:1"]
    origin = merged.feature_metadata["negative:1"]["origin"]
    assert origin["polarity"] == "negative" and origin["feature_id"] == "1"
    assert origin["handle_id"] == "dh_n"
    assert origin["source_files"] and all("/" not in f and "\\" not in f
                                          for f in origin["source_files"])
    assert merged.feature_metadata["positive:1"]["name"] == "PC 34:1"
    # 副次データは両側の ID に付け替わる（derive の絞り込みで落ちない）
    assert merged.feature_tags == {"positive:1": ["pos-tag"], "negative:1": ["neg-tag"]}
    assert set(merged.feature_candidates) <= set(merged.feature_ids)
    assert merged.per_sample_mz["positive:1"] == {"S1": 600.0, "S2": 601.0, "S3": 602.0}
    assert merged.per_sample_mz["negative:1"] == {"S1": 700.0, "S2": 701.0, "S3": 702.0}
    assert set(pos.feature_tags) == {"1"}  # 元は不変
    assert merged.msdial_samples[0]["name"] == "S1"
    assert merged.roles.keys() == {"S1", "S2", "S3"} and merged.sample_meta.keys() == {"S1", "S2", "S3"}


def test_merged_quant_state_and_record(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    merged, info = _merge(pos, neg, mapping_rows=None)
    st = merged.quant_state
    assert st["normalized"] is True and st["unit"] == "nmol" and st["method"] == "internal_standard"
    assert st["steps"] == ["normalize", "merge"] and st["record_id"] == info["record_id"]
    assert st["internal_standard_ids"] == ["positive:0", "negative:0"]
    assert "curation_by_polarity" not in st  # 読み込み時の写しは持たない（選別がその場で読む）
    arf2s = st["arf2_by_polarity"]
    assert arf2s["positive"] == next(p for p in pos.source_files if p.endswith(".arf2"))
    assert arf2s["negative"] == next(p for p in neg.source_files if p.endswith(".arf2"))
    assert merged.curation_flags is None
    rec = json.loads((tmp_path / "run" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["step"] == "merge" and rec["parent_ids"] == ["q_pos", "q_neg"]
    assert any(p.endswith(".arf2") for p in rec["inputs"]) and len(rec["inputs"]) >= 4
    assert rec["settings"]["rules"]["digest"]


def test_merge_of_skipped_normalization(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    pos.quant_state = neg.quant_state = None
    p2, _ = qn.skip(pos, "positive")
    n2, _ = qn.skip(neg, "negative")
    merged, _ = _merge(p2, n2, mapping_rows=None)
    assert merged.quant_state["normalized"] is False
    assert merged.quant_state["skipped"] == ["normalize"]
    assert merged.quant_state["steps"] == ["merge"]


def test_record_in_common_parent(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    _merged, info = _merge(pos, neg, mapping_rows=None)
    assert (tmp_path / "run" / "reports" / "quant" / f"{info['record_id']}.json").is_file()
    assert not (tmp_path / "run" / "pos" / "reports").exists()


def test_nested_folders_need_record_folder(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos = _load(tmp_path / "run", ["S1_pos", "S2_pos", "S3_pos"], [_spot("PC 34:1", [1.0, 2.0, 3.0])])
    neg = _load(tmp_path / "run" / "neg", ["S1_neg", "S2_neg", "S3_neg"], [_spot("PE 36:2", [3.0, 4.0, 5.0], "PE")])
    _state(pos), _state(neg)
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_RECORD_LOCATION_REQUIRED"
    where = tmp_path / "records"
    where.mkdir()
    _merged, info = _merge(pos, neg, mapping_rows=None, folder=where)
    assert (where / "reports" / "quant" / f"{info['record_id']}.json").is_file()


def test_group_or_role_disagreement_stops(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    neg.msdial_samples[0]["group"] = "other"
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_MERGE_SAMPLE_CONFLICT"
    assert exc.value.details["conflicts"][0]["sample_id"] == "S1"
    pos, neg = _pair(tmp_path, root="run2")
    neg.roles[neg.sample_names[1]] = "qc"
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_MERGE_SAMPLE_CONFLICT"


def test_rules_unavailable(tmp_path, monkeypatch):
    pos, neg = _pair(tmp_path)
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(tmp_path / "absent.toml"))
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "quant_rules_unavailable"


# ---- derive の小さな変更 ----

def test_derive_rejects_duplicate_sample_index_and_renames_msdial_samples(tmp_path, monkeypatch):
    pos, _neg = _pair(tmp_path)
    ones = np.ones((len(pos.feature_ids), 2))
    mask = np.ones(ones.shape, dtype=bool)
    with pytest.raises(ValueError):
        qd.derive(pos, matrix=ones, detected=mask, sample_names=["a", "b"], sample_index=[0, 0],
                  state=pos.quant_state)
    new = qd.derive(pos, matrix=ones, detected=mask, sample_names=["a", "b"],
                    sample_index=[1, 0], state=pos.quant_state)
    assert [m["name"] for m in new.msdial_samples] == ["a", "b"]
    assert pos.msdial_samples[0]["name"] == "S1_pos"


# ---- ツール ----

def _open(ds, polarity):
    session_state.session.dataset = ds
    holder = handles.registry.get(session_state.session.current_handle_id)
    holder.polarity = polarity
    return holder.handle_id


def _tool():
    def run(hp, hn, mapping_csv=None, record_folder=None):
        return qm.merge_tool(handles.registry.get(hp), handles.registry.get(hn),
                             mapping_csv, record_folder)
    return run


def test_tool_merges_into_new_handle(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    out = json.loads(_tool()(hp, hn))
    assert out["status"] == "merged" and out["dataset_handle"] not in (hp, hn)
    assert out["rule"] == "polarity_token_removed" and out["n_samples"] == 3
    assert str(tmp_path) not in json.dumps(out)
    merged = handles.registry.get(out["dataset_handle"])
    assert merged.polarity is None  # 陽陰の両方なので極性は決まらない
    assert merged.source["parents"] == [hp, hn]
    assert merged.dataset.feature_metadata["negative:0"]["origin"]["handle_id"] == hn
    assert merged.arf2_path is None
    assert session_state.session.current_handle_id == out["dataset_handle"]
    assert merged.dataset.quant_state["record_id"] == out["record_id"]
    assert handles.registry.get(hp).dataset is pos  # 元のハンドルは残る


def test_tool_mismatch_needs_input_and_mapping(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path, neg_names=("zz", "yy", "ww"))
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    out = json.loads(_tool()(hp, hn))
    assert out["status"] == "needs_input" and out["reason"] == "sample_mismatch"
    assert out["rule"] == "exact_name"
    assert set(out["only_in"]) == {"positive", "negative"}
    assert sorted(out["only_in"]["negative"]) == ["ww", "yy", "zz"]
    assert "mode, source_sample, sample_id, group, sample_type" in out["how_to_answer"]
    assert not (tmp_path / "run" / "reports").exists()
    mapping = tmp_path / "map.csv"
    mapping.write_text("mode,source_sample,sample_id,group,sample_type\n"
                       "pos,S1_pos,X1,A,s\npos,S2_pos,X2,A,s\npos,S3_pos,X3,A,s\n"
                       "neg,zz,X1,A,s\nneg,yy,X2,A,s\nneg,ww,X3,A,s\n",
                       encoding="utf-8")
    ok = json.loads(_tool()(hp, hn, mapping_csv=str(mapping)))
    assert ok["status"] == "merged" and ok["rule"] == "mapping_csv"
    rec = json.loads((tmp_path / "run" / "reports" / "quant"
                      / f"{ok['record_id']}.json").read_text(encoding="utf-8"))
    assert "mapping:map.csv" in rec["inputs"]
    assert str(tmp_path) not in json.dumps(ok)


def test_tool_nested_needs_record_folder(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos = _state(_load(tmp_path / "run", ["S1_pos", "S2_pos", "S3_pos"], [_spot("PC 34:1", [1.0, 2.0, 3.0])]))
    neg = _state(_load(tmp_path / "run" / "neg", ["S1_neg", "S2_neg", "S3_neg"], [_spot("PE 36:2", [3.0, 4.0, 5.0], "PE")]))
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    out = json.loads(_tool()(hp, hn))
    assert out["status"] == "needs_input" and out["reason"] == "record_folder_required"
    assert "record_folder" in out["how_to_answer"]
    bad = json.loads(_tool()(hp, hn, record_folder=str(tmp_path / "missing")))
    assert bad["error"]["code"] == "QUANT_RECORD_FOLDER_INVALID"
    where = tmp_path / "records"
    where.mkdir()
    ok = json.loads(_tool()(hp, hn, record_folder=str(where)))
    assert ok["status"] == "merged"
    assert (where / "reports" / "quant" / f"{ok['record_id']}.json").is_file()


def test_tool_refusals_are_coded(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    hp, hn = _open(pos, "positive"), _open(neg, "positive")  # 同じ極性
    assert json.loads(_tool()(hp, hn))["error"]["code"] == "QUANT_MERGE_POLARITY_INVALID"
    hn2 = _open(_state(_pair(tmp_path, root="r2")[1], None), "negative")
    assert json.loads(_tool()(hp, hn2))["error"]["code"] in (
        "QUANT_MERGE_POLARITY_INVALID", "QUANT_MERGE_STATE_MISMATCH")


def test_tool_normalize_step_required(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    neg.quant_state = None
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    err = json.loads(_tool()(hp, hn))["error"]
    assert err["code"] == "QUANT_NORMALIZE_STEP_REQUIRED" and "dataset_quant_apply" in err["message"]


def test_tool_rules_unavailable_envelope(tmp_path, monkeypatch):
    pos, neg = _pair(tmp_path)
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(tmp_path / "absent.toml"))
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    assert json.loads(_tool()(hp, hn))["error"]["code"] == "quant_rules_unavailable"


def test_merge_tool_registers_in_subprocess():
    code = (
        "from metabolomix.core.mcp_core import mcp\n"
        "import metabolomix.quant.tools as t\n"
        "tm = mcp._tool_manager._tools\n"
        "n = 'dataset_merge_polarities'\n"
        "assert n in t.__all__ and n in tm\n"
        "assert tm[n].annotations is not None and tm[n].output_schema is None\n"
        "assert {'positive_handle','negative_handle','mapping_csv','record_folder'} "
        "<= set(tm[n].parameters['properties'])\n"
        "print('ok')\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and "ok" in out.stdout, out.stderr


def test_merge_refuses_non_msdial_source(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    neg.source_format = "mztab"
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_SOURCE_UNSUPPORTED"
    assert not (tmp_path / "run" / "reports").exists()


def test_merge_requires_normalize_step_not_just_any_state(tmp_path, monkeypatch):
    """quant_state があっても、正規化（か、その skip）を通っていなければ統合しない。"""
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    for ds, rid in ((pos, "q_p"), (neg, "q_n")):
        ds.quant_state = qd.quant_state(normalized=False, method=None, unit=None, record_id=rid,
                                        steps=["select"], skipped=[])
    with pytest.raises(DomainError) as exc:
        _merge(pos, neg, mapping_rows=None)
    assert exc.value.code == "QUANT_NORMALIZE_STEP_REQUIRED"


def test_tool_checks_handle_capacity_before_writing(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    pos, neg = _pair(tmp_path)
    hp, hn = _open(pos, "positive"), _open(neg, "negative")
    cfg = tmp_path / "lipidmix.local.toml"  # _point が書いた設定に上限を足す（B-2 で環境変数は無い）
    cfg.write_text(cfg.read_text(encoding="utf-8") + "\n[server]\nmax_dataset_handles = 2\n",
                   encoding="utf-8")
    user_config._cache.clear()
    out = json.loads(_tool()(hp, hn))
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert not (tmp_path / "run" / "reports").exists()


def test_read_mapping_non_utf8_is_bad_request(tmp_path):
    path = tmp_path / "deep" / "map.csv"
    path.parent.mkdir()
    path.write_bytes("mode,source_sample,sample_id,group,sample_type\npos,試料1,S1,g,x\n"
                     .encode("cp932"))
    with pytest.raises(DomainError) as exc:
        qm.read_mapping(path)
    assert exc.value.code == "QUANT_BAD_REQUEST"
    assert "map.csv" in exc.value.message and str(path.parent) not in exc.value.message
