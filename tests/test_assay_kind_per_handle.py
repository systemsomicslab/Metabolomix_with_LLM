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


# ---------- 意味論ダイジェストは種別ごとにセッションで 1 回（レビュー指摘） ----------


def test_digest_once_per_kind_across_handles(tmp_path, fresh_state):
    """pos を開く → lipid を記録 → neg を開く → lipid を記録、でダイジェストは 2 回だけ。

    種別はハンドルごとだが、ダイジェストは種別ごとにセッションで 1 回（2 つ目のハンドルを開いても
    unknown 用を、その種別を記録しても lipid 用を出し直さない）。
    """
    def digest(text: str) -> bool:
        return text.startswith(session_state.SEMANTICS_CAVEAT)

    out_a = _load(directory=str(_folder(tmp_path / "pos")))
    assert digest(out_a.get("caveat") or "")
    assert digest(objective.record_objective("exp", assay_kind="lipid"))
    out_b = _load(directory=str(_folder(tmp_path / "neg")))
    assert out_b["assay_kind"] == "unknown"
    assert "caveat" not in out_b
    assert not digest(objective.record_objective("exp", assay_kind="lipid"))
    # まだ出していない種別になったときは、もう 1 回だけ届く。
    third = objective.record_objective("exp-m", assay_kind="metabolite")
    assert digest(third) and session_state._ASSAY_CAVEAT_LINE["metabolite"] in third
    assert not digest(objective.record_objective("exp-m2", assay_kind="metabolite"))


# ---------- データセットを開く前に記録した種別（レビュー指摘） ----------


def test_kind_recorded_before_any_dataset_carries_to_the_first_load(tmp_path, fresh_state,
                                                                    monkeypatch):
    """`[paths] reports_dir` があればデータセットが無くても目的を書ける。そのとき渡した種別はセッションに置き、
    次に直近になるハンドルへ移す（ライブラリと同じ規則。2 本目のハンドルへは移さない）。"""
    from tests.config_fixtures import use_config

    use_config(monkeypatch, tmp_path, {"paths.reports_dir": tmp_path / "shared"})
    out = objective.record_objective("exp", assay_kind="lipid")
    assert "データセットが無いのでセッションに置いた" in out
    assert "セッション のデータセット" not in out
    assert session_state.session.assay_kind == "lipid"

    first = _load(directory=str(_folder(tmp_path / "a")))
    assert first["assay_kind"] == "lipid"
    assert any("引き継いだ" in note for note in first.get("notes") or []), first
    assert session_state.session._assay_kind == "unknown"   # 移した（写していない）
    second = _load(directory=str(_folder(tmp_path / "b")))
    assert second["assay_kind"] == "unknown"
    assert handles.registry.peek(first["dataset_handle"]).assay_kind == "lipid"


def test_carried_kind_that_disagrees_with_the_folder_is_not_applied(tmp_path, fresh_state,
                                                                     monkeypatch):
    """目的ファイルから読み戻した種別と、開く前に記録した種別が食い違えば、どちらも採らない。"""
    folder = _folder(tmp_path / "a")
    first = _load(directory=str(folder))["dataset_handle"]
    objective.record_objective("exp-folder", assay_kind="metabolite")   # 解析フォルダ配下 reports/
    handle_tools.dataset_close(first)
    assert session_state.session.dataset is None      # 直近を閉じた（ハンドルが無い）

    from tests.config_fixtures import use_config
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": tmp_path / "shared"})
    objective.record_objective("exp-session", assay_kind="lipid")       # データセットが無いうちに
    out = _load(directory=str(folder))
    assert out["assay_kind"] == "unknown"
    assert any("食い違う" in note for note in out.get("notes") or []), out
    assert session_state.session._assay_kind == "unknown"


def test_carried_kind_stays_on_the_session_when_the_handle_limit_stops_the_load(tmp_path,
                                                                                fresh_state,
                                                                                monkeypatch):
    """上限で登録できなければ、種別はセッションに残る（ハンドルごと捨てられる）。"""
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.mztab.dataset_state import DatasetState
    from tests.config_fixtures import use_config

    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    session = session_state.session
    session.dataset = DatasetState()          # 1 本目で上限
    session.dataset = None                    # 直近だけ外す（ハンドルは開いたまま）
    session.assay_kind = "lipid"              # 直近が無いのでセッション自身へ
    applied = []
    extra = handles.DatasetHandle(DatasetState(), polarity=None, polarity_source=None,
                                  source={"kind": "session"})
    try:
        session.adopt_handle(extra, take_session_assay_kind=lambda kind, analysis_id: applied.append(
            (kind, analysis_id)))
    except DomainError as exc:
        assert exc.code == "DATASET_HANDLE_LIMIT"
    else:
        raise AssertionError("上限で止まるはず")
    assert session._assay_kind == "lipid"
    assert session.current_handle_id is None


# ---------- pathway_project は渡したハンドルの種別を見る（レビュー指摘。脂質の関門は担当 9 で外し、
# 種別は notes.assay_kind の注意だけに効く。注意の読み先は tests/test_pathway_tools.py が縛る） ----------


def test_pathway_gate_reads_the_passed_handle(fresh_state, monkeypatch):
    from metabolomix.mztab.dataset_state import DatasetState
    from metabolomix.pathway import bridge
    from metabolomix.pathway import tools as pathway_tools
    from tests.test_stats_results import _ds

    class _Reached(Exception):
        pass

    def stop(*_args, **_kwargs):
        raise _Reached

    monkeypatch.setattr(bridge, "run_bridge", stop)
    a = _ds(groups=("C", "P"), per_group=2, reps=1, n_qc=0)
    a.assay_kind = "lipid"
    a.feature_annotations = {"f0": {"name": "PC 34:1", "ontology": "PC"},
                             "f1": {"name": "LPC 18:1", "ontology": "LPC"}}
    session_state.session.dataset = a
    handle_a = session_state.session.current_handle_id
    session_state.session.dataset = DatasetState()        # 直近は unknown のハンドル
    assert session_state.session.assay_kind == "unknown"
    try:
        out = pathway_tools.pathway_project(
            group_a=["C0_r0", "C1_r0"], group_b=["P0_r0", "P1_r0"],
            group_a_label="C", group_b_label="P", dataset_handle=handle_a)
    except _Reached:
        return                                             # 関門を越えてブリッジまで来た
    error = (json.loads(out).get("error") or {}) if out.lstrip().startswith("{") else {}
    assert error.get("state") != "assay_kind", out


# ---------- dataset_close の入力の名前（レビュー指摘） ----------


def test_dataset_close_names_a_pipeline_run_by_its_folder(fresh_state, tmp_path):
    from metabolomix.mztab.dataset_state import DatasetState

    run = tmp_path / "run-20261114"
    handle = handles.DatasetHandle(
        DatasetState(), polarity="positive", polarity_source="test",
        source={"kind": "pipeline_path", "path": str(run / "pipeline-run.json"),
                "polarity": "positive"})
    handles.registry.add(handle)
    out = json.loads(handle_tools.dataset_close(handle.handle_id))
    assert out["source"]["name"] == "run-20261114"
    assert out["source"]["polarity"] == "positive"


# ---------- 共有の置き場所でも目的は解析フォルダごと（レビュー指摘 R14） ----------


def _meta(path) -> dict:
    from metabolomix.tools import notes

    meta, _ = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
    return meta


def _reopen(folder) -> dict:
    """新しい会話（サーバの再起動の後と同じ）でフォルダを開き直す。"""
    for handle in handles.registry.list():
        handle_tools.dataset_close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()
    return _load(directory=str(folder))


def test_same_analysis_id_on_two_folders_keeps_each_kind(tmp_path, fresh_state, monkeypatch):
    """`[paths] reports_dir` を共有しても、同じ analysis_id を 2 つのデータセットで記録した種別が混ざらない。

    pos の目的は共有の置き場所に、neg の目的は（共有の置き場所の同名の目的が pos のものなので）neg の
    解析フォルダ配下 reports/ に作る。開き直すと、それぞれが自分の種別を読み戻す。
    """
    from tests.config_fixtures import use_config

    shared = tmp_path / "shared"
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": shared})
    pos_dir, neg_dir = _folder(tmp_path / "pos"), _folder(tmp_path / "neg")
    pos = _load(directory=str(pos_dir))["dataset_handle"]
    neg = _load(directory=str(neg_dir))["dataset_handle"]
    objective.record_objective("study", assay_kind="lipid", dataset_handle=pos)
    out = objective.record_objective("study", assay_kind="metabolite", dataset_handle=neg)
    assert "objective を記録" in out
    assert "別の解析フォルダ（pos）" in out and str(tmp_path) not in out
    assert _meta(shared / "study.objective.md")["assay_kind"] == "lipid"          # pos のものは動かない
    assert objective._same_dir(_meta(shared / "study.objective.md")["analysis_dir"], pos_dir)
    own = neg_dir / "reports" / "study.objective.md"
    assert _meta(own)["assay_kind"] == "metabolite"
    assert objective._same_dir(_meta(own)["analysis_dir"], neg_dir)
    # 以後の neg の更新も neg の目的へ（pos の目的は書き換えない）。
    updated = objective.record_objective("study", confirmed_objective="neg の確定", dataset_handle=neg)
    assert "objective を更新" in updated
    assert _meta(own)["confirmed_objective"] == "neg の確定"
    assert not _meta(shared / "study.objective.md").get("confirmed_objective")
    assert handles.registry.peek(pos).assay_kind == "lipid"
    assert handles.registry.peek(neg).assay_kind == "metabolite"

    again_pos = _reopen(pos_dir)
    assert again_pos["assay_kind"] == "lipid"
    assert any("lipid" in note and "目的ファイル" in note for note in again_pos.get("notes") or [])
    again_neg = _reopen(neg_dir)
    assert again_neg["assay_kind"] == "metabolite"
    assert any("metabolite" in note and "目的ファイル" in note for note in again_neg.get("notes") or [])


def test_same_kind_on_two_folders_is_read_back_by_both(tmp_path, fresh_state, monkeypatch):
    """1 つの研究の 2 極性に同じ analysis_id と同じ種別を記録する自然な使い方でも、neg が unknown に戻らない。"""
    from tests.config_fixtures import use_config

    use_config(monkeypatch, tmp_path, {"paths.reports_dir": tmp_path / "shared"})
    pos_dir, neg_dir = _folder(tmp_path / "pos"), _folder(tmp_path / "neg")
    pos = _load(directory=str(pos_dir))["dataset_handle"]
    neg = _load(directory=str(neg_dir))["dataset_handle"]
    objective.record_objective("study", assay_kind="lipid", dataset_handle=pos)
    objective.record_objective("study", assay_kind="lipid", dataset_handle=neg)
    assert _reopen(neg_dir)["assay_kind"] == "lipid"
    assert _reopen(pos_dir)["assay_kind"] == "lipid"


def test_legacy_shared_objective_without_a_folder_is_still_updated(tmp_path, fresh_state, monkeypatch):
    """frontmatter に analysis_dir の無い（データセットを開く前に書いた・古い）共有の目的は、今までどおり
    最初に書いたデータセットが引き取る（別のフォルダの目的とはみなさない）。"""
    from tests.config_fixtures import use_config

    shared = tmp_path / "shared"
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": shared})
    shared.mkdir()
    (shared / "study.objective.md").write_bytes(
        b"---\ntype: objective\nanalysis_id: study\nassay_kind: unknown\n---\n")
    pos_dir = _folder(tmp_path / "pos")
    _load(directory=str(pos_dir))
    out = objective.record_objective("study", assay_kind="lipid")
    assert "objective を更新" in out
    assert objective._same_dir(_meta(shared / "study.objective.md")["analysis_dir"], pos_dir)
    assert not (pos_dir / "reports" / "study.objective.md").exists()


# ---------- 開く前に記録した種別は目的ファイルにも残す（レビュー指摘 R15） ----------


def test_kind_recorded_before_any_dataset_is_read_back_on_reopen(tmp_path, fresh_state, monkeypatch):
    from tests.config_fixtures import use_config

    shared = tmp_path / "shared"
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": shared})
    objective.record_objective("exp", assay_kind="lipid")
    assert "analysis_dir" not in _meta(shared / "exp.objective.md")     # まだどのフォルダのものでもない
    folder = _folder(tmp_path / "a")
    first = _load(directory=str(folder))
    assert first["assay_kind"] == "lipid"
    assert any("引き継いだ" in note and "exp.objective.md" in note and "残した" in note
               for note in first.get("notes") or []), first
    meta = _meta(shared / "exp.objective.md")
    assert meta["assay_kind"] == "lipid" and objective._same_dir(meta["analysis_dir"], folder)
    assert session_state.session._assay_kind_objective is None           # 引き継いだら忘れる

    # 同じ会話で閉じて開き直しても、新しい会話で開き直しても読み戻す。
    handle_tools.dataset_close(first["dataset_handle"])
    same = _load(directory=str(folder))
    assert same["assay_kind"] == "lipid"
    assert any("目的ファイル" in note and "読み戻した" in note for note in same.get("notes") or [])
    fresh = _reopen(folder)
    assert fresh["assay_kind"] == "lipid"
    assert any("読み戻した" in note for note in fresh.get("notes") or []), fresh
    # 2 本目のフォルダへは引き継がない（目的は 1 本目のもの）。
    assert _reopen(_folder(tmp_path / "b"))["assay_kind"] == "unknown"


def test_carried_kind_without_its_objective_says_it_is_not_saved(tmp_path, fresh_state):
    """目的を介さずにセッションへ置いた種別（どの目的のものか分からない）は、引き継ぐが残せないと言う。"""
    session_state.session.assay_kind = "lipid"
    out = _load(directory=str(_folder(tmp_path / "a")))
    assert out["assay_kind"] == "lipid"
    assert any("引き継いだ" in note and "record_objective(assay_kind=...)" in note
               for note in out.get("notes") or []), out


def test_carried_kind_is_not_saved_through_a_busy_handle(tmp_path, monkeypatch):
    """ハンドルのロックが塞がっていれば、目的へ書くのは上限まで待ってあきらめ、記録し直すよう言う。"""
    import threading

    from metabolomix.core import handle_dispatch
    from metabolomix.mztab.dataset_state import DatasetState

    monkeypatch.setattr(handle_dispatch, "HANDLE_LOCK_WAIT", 0.05)
    ds = DatasetState()
    ds.analysis_dir = str(tmp_path / "a")
    handle = handles.DatasetHandle(ds, polarity=None, polarity_source=None, source={"kind": "session"})
    carry = handle_tools._SessionKindCarry(handle, None, None)
    carry("lipid", "exp")
    held, release = threading.Event(), threading.Event()

    def hold():
        with handle.lock:
            held.set()
            release.wait(5)

    thread = threading.Thread(target=hold, daemon=True)
    thread.start()
    try:
        assert held.wait(5)
        carry.persist()
    finally:
        release.set()
        thread.join(5)
    assert "塞がっていた" in carry.note and "record_objective(assay_kind=...)" in carry.note
    assert not (tmp_path / "a" / "reports").exists()


def test_remembered_objective_is_dropped_when_the_kind_changes_elsewhere(fresh_state):
    session = session_state.session
    session.assay_kind = "lipid"
    session.remember_session_kind_objective("exp", "lipid")
    assert session._assay_kind_objective == "exp"
    session.assay_kind = "metabolite"            # 目的を介さずに書き換えた
    assert session._assay_kind_objective is None
    session.remember_session_kind_objective("exp", "lipid")   # 種別がもう lipid でなければ覚えない
    assert session._assay_kind_objective is None


# ---------- データセットが無いときも、別の解析フォルダの共有の目的を書き換えない（レビュー指摘 R14 の残り） ----------


def test_record_before_load_leaves_another_folders_shared_objective_alone(tmp_path, fresh_state,
                                                                          monkeypatch):
    """データセットを開く前の記録（解析フォルダの無い呼び出し）も、共有の置き場所にある別の解析フォルダの
    目的を書き換えない。まだどのフォルダのものでもない目的として別の名前で作り、次に開いたデータセットへ
    結び付ける。元のフォルダは開き直しても自分の種別を読み戻す。"""
    from tests.config_fixtures import use_config

    shared = tmp_path / "shared"
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": shared})
    pos_dir, neg_dir = _folder(tmp_path / "pos"), _folder(tmp_path / "neg")
    pos = _load(directory=str(pos_dir))["dataset_handle"]
    objective.record_objective("study", assay_kind="lipid", dataset_handle=pos)
    owned = shared / "study.objective.md"
    assert objective._same_dir(_meta(owned)["analysis_dir"], pos_dir)

    # サーバの再起動の後、まだ何も開かずに同じ analysis_id へ記録する。
    handle_tools.dataset_close(pos)
    session_state.session = session_state.AnalysisSession()
    out = objective.record_objective("study", assay_kind="metabolite", confirmed_objective="次の確定")
    assert "objective を記録" in out, out
    assert "別の解析フォルダ（pos）" in out and str(tmp_path) not in out
    assert _meta(owned)["assay_kind"] == "lipid"                 # pos の目的は動かない
    assert not _meta(owned).get("confirmed_objective")
    unbound = shared / "study-2.objective.md"
    meta = _meta(unbound)
    assert meta["analysis_id"] == "study" and meta["assay_kind"] == "metabolite"
    assert meta["confirmed_objective"] == "次の確定" and "analysis_dir" not in meta
    # 開く前の続きの記録も、その結び付いていない目的へ（pos の目的ではない）。
    more = objective.record_objective("study", add_subquestions=["追記の小問"])
    assert "objective を更新: study-2.objective.md" in more, more
    assert "追記の小問" not in owned.read_text(encoding="utf-8")

    # 次に開いた neg が種別を引き継ぎ、その目的に neg の解析フォルダが付く。
    first = _load(directory=str(neg_dir))
    assert first["assay_kind"] == "metabolite"
    assert any("引き継いだ" in note and "study-2.objective.md" in note and "残した" in note
               for note in first.get("notes") or []), first
    assert objective._same_dir(_meta(unbound)["analysis_dir"], neg_dir)
    assert _meta(owned)["assay_kind"] == "lipid"

    again_pos = _reopen(pos_dir)
    assert again_pos["assay_kind"] == "lipid"
    assert any("study.objective.md" in note for note in again_pos.get("notes") or []), again_pos
    assert _reopen(neg_dir)["assay_kind"] == "metabolite"


def test_numbered_objective_is_not_taken_for_its_file_name(tmp_path, fresh_state, monkeypatch):
    """番号付きの名前（`study-2.objective.md`）の目的は analysis_id が study のもの。analysis_id "study-2" の
    記録はそれを自分の目的とみなさず、上書きもしない。"""
    from tests.config_fixtures import use_config

    shared = tmp_path / "shared"
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": shared})
    pos_dir = _folder(tmp_path / "pos")
    pos = _load(directory=str(pos_dir))["dataset_handle"]
    objective.record_objective("study", assay_kind="lipid", dataset_handle=pos)
    handle_tools.dataset_close(pos)
    session_state.session = session_state.AnalysisSession()
    objective.record_objective("study", confirmed_objective="結び付いていない")
    numbered = shared / "study-2.objective.md"
    assert _meta(numbered)["analysis_id"] == "study"

    # データセットの無い呼び出しも、解析フォルダのある呼び出しも、study-2 を別の目的として扱う。
    out = objective.record_objective("study-2", confirmed_objective="別の研究")
    assert "objective を記録" in out, out
    assert _meta(numbered)["confirmed_objective"] == "結び付いていない"
    assert _meta(shared / "study-2-2.objective.md")["analysis_id"] == "study-2"
    other = _load(directory=str(_folder(tmp_path / "other")))["dataset_handle"]
    out = objective.record_objective("study-2", assay_kind="metabolite", dataset_handle=other)
    assert "objective を更新: study-2-2.objective.md" in out, out
    assert _meta(numbered)["confirmed_objective"] == "結び付いていない"
