"""アッセイ種別はハンドル（データセット）ごと（棚卸し C1-M02-M03・A9-M04-G06。dataset_close の戻り値は E8-T31）。

種別は `DatasetState.assay_kind` に住み、`session.assay_kind` は直近のハンドルの値を読む。
`record_objective(assay_kind=...)` が書いたハンドルだけが切り替わり、目的ファイルに残った種別は
同じ解析フォルダを `load_dataset` で開き直すと読み戻す（食い違えば読み戻さない）。
fixture はテスト自身が tmp に作る。
"""
from __future__ import annotations

import json

from metabolomix.core import handles, session_state
from metabolomix.tools import handle_tools, mztab_tools, objective, overview_tools
from tests.test_load_dataset import _folder, _load, fresh_state  # noqa: F401  (fixture)


def test_kind_is_per_handle(tmp_path, fresh_state):
    a = _load(directory=str(_folder(tmp_path / "a")))["dataset_handle"]
    out = objective.record_objective("exp-a", assay_kind="lipid")
    assert "ハンドル " + a in out
    b_out = _load(directory=str(_folder(tmp_path / "b")))
    b = b_out["dataset_handle"]
    # 新しいハンドルは unknown（前のハンドルの種別を持ち越さない）。直近が b なのでセッションも unknown。
    assert b_out["assay_kind"] == "unknown"
    assert session_state.session.assay_kind == "unknown"
    assert handles.registry.peek(a).assay_kind == "lipid"
    # 直近でないハンドルへも書ける（dataset_handle）。直近の種別は動かない。
    objective.record_objective("exp-b", assay_kind="metabolite", dataset_handle=a)
    assert handles.registry.peek(a).assay_kind == "metabolite"
    assert handles.registry.peek(b).assay_kind == "unknown"


def test_dataset_list_and_status_show_the_kind(tmp_path, fresh_state):
    a = _load(directory=str(_folder(tmp_path / "a")))["dataset_handle"]
    objective.record_objective("exp", assay_kind="lipid")
    b = _load(directory=str(_folder(tmp_path / "b")))["dataset_handle"]
    listing = json.loads(handle_tools.dataset_list())["handles"].splitlines()
    header = listing[0].split("\t")
    assert "assay_kind" in header
    column = header.index("assay_kind")
    kinds = {row.split("\t")[0]: row.split("\t")[column] for row in listing[1:]}
    assert kinds == {a: "lipid", b: "unknown"}
    assert json.loads(mztab_tools.dataset_status(dataset_handle=a))["assay_kind"] == "lipid"
    assert json.loads(mztab_tools.dataset_status())["assay_kind"] == "unknown"
    assert json.loads(overview_tools.dataset_overview(dataset_handle=a))["assay_kind"] == "lipid"


def test_reopening_reads_the_kind_back_from_the_objective(tmp_path, fresh_state):
    folder = _folder(tmp_path / "a")
    first = _load(directory=str(folder))["dataset_handle"]
    objective.record_objective("exp", assay_kind="lipid")
    assert (folder / "reports" / "exp.objective.md").exists()
    handle_tools.dataset_close(first)
    # 新しい会話（サーバの再起動の後と同じく、まだダイジェストを出していないセッション）で開き直す。
    session_state.session = session_state.AnalysisSession()

    out = _load(directory=str(folder))
    assert out["assay_kind"] == "lipid"
    assert any("目的ファイル" in note and "lipid" in note for note in out.get("notes") or []), out
    assert session_state.session.assay_kind == "lipid"
    # ダイジェストも読み戻した種別で届く。
    assert session_state._ASSAY_CAVEAT_LINE["lipid"] in (out.get("caveat") or "")


def test_conflicting_objectives_are_not_read_back(tmp_path, fresh_state):
    folder = _folder(tmp_path / "a")
    first = _load(directory=str(folder))["dataset_handle"]
    objective.record_objective("exp-1", assay_kind="lipid")
    objective.record_objective("exp-2", assay_kind="metabolite")
    handle_tools.dataset_close(first)

    out = _load(directory=str(folder))
    assert out["assay_kind"] == "unknown"
    assert any("食い違う" in note for note in out.get("notes") or []), out


def test_objective_elsewhere_counts_only_for_its_own_folder(tmp_path, fresh_state):
    """共有の置き場所（`[paths] reports_dir`）の目的は、frontmatter の analysis_dir が同じフォルダのものだけ。"""
    shared = tmp_path / "shared"
    shared.mkdir()
    a = _folder(tmp_path / "a")
    b = _folder(tmp_path / "b")
    (shared / "other.objective.md").write_bytes((
        "---\ntype: objective\nassay_kind: lipid\nanalysis_dir: " + json.dumps(str(b)) + "\n---\n"
    ).encode("utf-8"))
    ds_a = type("D", (), {"analysis_dir": str(a)})()
    ds_b = type("D", (), {"analysis_dir": str(b)})()
    assert objective.recall_assay_kind(ds_a, [shared, a / "reports"]) == (None, None)
    kind, note = objective.recall_assay_kind(ds_b, [shared, b / "reports"])
    assert kind == "lipid" and "other.objective.md" in note


def test_unknown_objective_does_not_override(tmp_path, fresh_state):
    folder = _folder(tmp_path / "a")
    first = _load(directory=str(folder))["dataset_handle"]
    objective.record_objective("exp", dataset="NEG")    # 種別を渡さない作成は unknown
    handle_tools.dataset_close(first)
    out = _load(directory=str(folder))
    assert out["assay_kind"] == "unknown"
    assert not any("目的ファイル" in note for note in out.get("notes") or [])


def test_session_kind_without_handles_stays_on_the_session(fresh_state):
    session = session_state.session
    assert session.current_handle_id is None
    session.assay_kind = "metabolite"
    assert session.assay_kind == "metabolite"


# ---------- dataset_close の戻り値（E8-T31） ----------


def test_dataset_close_reports_what_was_closed(tmp_path, fresh_state):
    a = _load(directory=str(_folder(tmp_path / "alpha")))["dataset_handle"]
    b = _load(directory=str(_folder(tmp_path / "beta")))["dataset_handle"]
    out = json.loads(handle_tools.dataset_close(a))
    assert out["closed"] == a and out["was_current"] is False
    assert out["source"]["name"] == "alpha"
    assert out["source"]["kind"]
    assert str(tmp_path) not in json.dumps(out, ensure_ascii=False)   # 置き場所は出さない
    assert out["open_handles"] == 1
    out = json.loads(handle_tools.dataset_close(b))
    assert out["was_current"] is True and out["open_handles"] == 0
