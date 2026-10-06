"""console_plan のメソッドファイル解決（LBM 自動解決・自動保存パラメータの探索）。

純ロジックは tests/test_console_method_file.py。ここは MCP ツール層の結線を見る。
"""
from __future__ import annotations

import json as _json
from pathlib import Path


def _fake_exe_with_lbm(tmp_path, *names):
    """MSDIAL_EXE と同じフォルダに脂質ライブラリを置く（GUI と同じ配置）。"""
    d = tmp_path / "msdial_app"
    d.mkdir(exist_ok=True)
    for n in names:
        (d / n).touch()
    return d / "MSDIALCUI.exe"


def _plan_ready(tmp_path, monkeypatch, exe):
    from metabolomix.core import session_state
    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", str(exe))
    monkeypatch.delenv("MSDIAL_LBM", raising=False)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "a.wiff").touch()


def _nested_root(tmp_path):
    """dataset_root を 1 段ネストする。

    tmp_path をそのまま dataset_root にすると、兄弟フォルダ探索が pytest の共有親
    （pytest-<N>/）を走査して他テストの一時フォルダを拾う。
    """
    root = tmp_path / "dataset"
    root.mkdir(exist_ok=True)
    (root / "a.wiff").touch()
    return root


def test_console_plan_rejects_lipidomics_without_lbm(tmp_path, monkeypatch):
    """LBM 空のまま走らせると警告なしで同定 0 件になる（実測 60 サンプル 31 分）。"""
    exe = _fake_exe_with_lbm(tmp_path)  # .lbm2 を置かない
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: \n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "LBM_NOT_FOUND"


def test_console_plan_resolves_lbm_from_exe_directory(tmp_path, monkeypatch):
    """GUI と同じく MSDIAL_EXE と同じフォルダの *.lbm2 を 1 件だけ自動採用する。"""
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: \n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["status"] == "planned"
    assert parsed["lbm"]["source"] == "exe_dir"
    assert Path(parsed["lbm"]["path"]).name == "Msp_lipids.lbm2"


def test_console_plan_writes_effective_method_file_without_touching_source(tmp_path, monkeypatch):
    """解決した LBM は run_dir の写しに書く。ユーザーのファイルは変えない。"""
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    original = "Ion mode: Negative\nLbm file path: \n"
    method.write_text(original, encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    effective = Path(parsed["method_file"])
    assert effective != method
    assert effective.parent == Path(parsed["run_dir"])
    assert "Msp_lipids.lbm2" in effective.read_text(encoding="ascii")
    assert method.read_text(encoding="ascii") == original


def test_console_plan_job_points_at_effective_method_file(tmp_path, monkeypatch):
    """console_run は job.method_file を Console へ渡すので、そこが写しでなければ意味がない。"""
    from metabolomix.console.job_manager import load_job
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: \n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    job = load_job(Path(parsed["job_path"]))
    assert "Msp_lipids.lbm2" in Path(job.method_file).read_text(encoding="ascii")


def test_console_plan_rejects_ambiguous_lbm(tmp_path, monkeypatch):
    """GUI も候補が 1 件でなければ実行を止める。"""
    exe = _fake_exe_with_lbm(tmp_path, "a.lbm2", "b.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "LBM_AMBIGUOUS"


def test_console_plan_lbm_file_argument_overrides(tmp_path, monkeypatch):
    exe = _fake_exe_with_lbm(tmp_path, "installed.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    chosen = tmp_path / "chosen.lbm2"
    chosen.touch()
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      lbm_file=str(chosen)))
    assert Path(parsed["lbm"]["path"]).name == "chosen.lbm2"
    assert parsed["lbm"]["source"] == "argument"


def test_console_plan_metabolomics_does_not_require_lbm(tmp_path, monkeypatch):
    exe = _fake_exe_with_lbm(tmp_path)
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      omics="metabolomics"))
    assert parsed["status"] == "planned"
    assert parsed["lbm"]["source"] == "not_required"


def test_console_plan_rejects_polarity_mismatch(tmp_path, monkeypatch):
    """メソッドの Ion mode と宣言極性がずれたまま走ると、別極性の結果が黙って出る。"""
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="positive", measure="peak_height"))
    assert parsed["error"]["code"] == "METHOD_FILE_POLARITY_MISMATCH"


def test_console_plan_discovers_auto_saved_method_file(tmp_path, monkeypatch):
    """GUI は実行のたびに *_param_<ts>.txt を自動保存する。渡さなくても見つける。"""
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    (root / "Dataset_2026_param_202605151055.txt").write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\nLbm file path: \n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root),
                                      polarity="negative", measure="peak_height"))
    assert parsed["status"] == "planned"
    assert parsed["method_source"]["discovered_from"].endswith(
        "Dataset_2026_param_202605151055.txt")


def test_console_plan_discovery_ignores_other_polarity(tmp_path, monkeypatch):
    """NEG の自動保存パラメータを POS の計画に流用してはいけない。

    候補として提示はするが、採用はしない（METHOD_FILE_CHOICE_REQUIRED）。
    """
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    (root / "Dataset_2026_param_202605151055.txt").write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root),
                                      polarity="positive", measure="peak_height"))
    assert parsed["error"]["code"] == "METHOD_FILE_CHOICE_REQUIRED"
    assert parsed["error"]["details"]["candidates"][0]["usable"] == "needs_polarity_conversion"


def test_console_plan_missing_method_file_reports_no_candidates(tmp_path, monkeypatch):
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "METHOD_FILE_NOT_GIVEN"


def test_console_plan_offers_a_sibling_folder_candidate(tmp_path, monkeypatch):
    """実データの形（POS の隣に GUI 処理済みの NEG がある）。黙って採用せず提示する。"""
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    neg = tmp_path / "NEG"
    neg.mkdir()
    (neg / "Dataset_2026_param_202605151055.txt").write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\n", encoding="ascii")

    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root), polarity="positive"))
    error = parsed["error"]
    assert error["code"] == "METHOD_FILE_CHOICE_REQUIRED"
    entry = error["details"]["candidates"][0]
    assert entry["origin"] == "sibling"
    assert entry["usable"] == "needs_polarity_conversion"
    assert "console_method_template" in error["required_tools"]


def test_console_plan_adopts_a_sibling_of_the_same_polarity(tmp_path, monkeypatch):
    """極性が一致していれば兄弟フォルダのものでも採用してよい（変換が要らない）。"""
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    other = tmp_path / "OTHER"
    other.mkdir()
    (other / "Dataset_2026_param_202605151055.txt").write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\nLbm file path: \n", encoding="ascii")

    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root), polarity="negative"))
    assert parsed["status"] == "planned"
    assert parsed["method_source"]["discovered_from"].endswith(
        "Dataset_2026_param_202605151055.txt")


def test_console_plan_choice_required_truncates_many_candidates(tmp_path, monkeypatch):
    """兄弟フォルダが多いラボ配置でも METHOD_FILE_CHOICE_REQUIRED の戻り値を
    無制限に太らせない。n_candidates 相当のメッセージ中の件数は総数のまま。"""
    from metabolomix.console.method_file import MAX_REPORTED_CANDIDATES
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    n_siblings = MAX_REPORTED_CANDIDATES + 2
    for i in range(n_siblings):
        sib = tmp_path / f"SIB{i}"
        sib.mkdir()
        (sib / f"d{i}_param_1.txt").write_text(
            "Ion mode: Negative\nTarget omics: Lipidomics\n", encoding="ascii")

    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root), polarity="positive"))
    error = parsed["error"]
    assert error["code"] == "METHOD_FILE_CHOICE_REQUIRED"
    assert f"{n_siblings} 件" in error["message"]
    assert len(error["details"]["candidates"]) == MAX_REPORTED_CANDIDATES
    assert error["details"]["truncated"] is True


def test_console_plan_not_given_reports_where_it_looked(tmp_path, monkeypatch):
    exe = _fake_exe_with_lbm(tmp_path, "x.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    root = _nested_root(tmp_path)
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(root), polarity="positive"))
    assert parsed["error"]["code"] == "METHOD_FILE_NOT_GIVEN"
    assert str(root) in parsed["error"]["details"]["searched"]


# ---------- 相対パス宣言は実効メソッドへ絶対パスで書く ----------
#
# LC-MS の Console（`ConfigParser.ReadForLcmsParameter`）は宣言パスを解決せず、
# 相対値は Console プロセスの cwd（＝run_dir）基準になる。見つからなければ
# ライブラリを黙って飛ばし、同定 0 件で完走する。


def test_console_plan_writes_a_relative_lbm_declaration_as_an_absolute_path(tmp_path, monkeypatch):
    from metabolomix.console import method_file as method_file_mod
    from metabolomix.console.job_manager import load_job
    exe = _fake_exe_with_lbm(tmp_path)
    _plan_ready(tmp_path, monkeypatch, exe)
    proj = tmp_path / "proj"
    proj.mkdir()
    lib = proj / "lab.lbm2"
    lib.touch()
    method = proj / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: lab.lbm2\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["lbm"]["source"] == "method_file"
    job = load_job(Path(parsed["job_path"]))
    effective = Path(job.method_file)
    assert effective != method
    declared = method_file_mod.read_method_keys(effective)["lbm file path"]
    assert Path(declared).is_absolute()
    assert Path(declared) == lib.resolve()


def test_console_plan_writes_relative_msp_and_text_db_as_absolute_paths(tmp_path, monkeypatch):
    from metabolomix.console import method_file as method_file_mod
    from metabolomix.console.job_manager import load_job
    exe = _fake_exe_with_lbm(tmp_path)
    _plan_ready(tmp_path, monkeypatch, exe)
    proj = tmp_path / "proj"
    proj.mkdir()
    (proj / "lib.msp").touch()
    (proj / "db.txt").touch()
    method = proj / "params.txt"
    method.write_text("Ion mode: Negative\nMsp file path: lib.msp\nText DB file path: db.txt\n",
                      encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      omics="metabolomics"))
    assert parsed["status"] == "planned"
    keys = method_file_mod.read_method_keys(Path(load_job(Path(parsed["job_path"])).method_file))
    assert Path(keys["msp file path"]) == (proj / "lib.msp").resolve()
    assert Path(keys["text db file path"]) == (proj / "db.txt").resolve()


def test_console_plan_metabolomics_writes_a_declared_lbm_as_an_absolute_path(tmp_path, monkeypatch):
    """Console は Target omics に関わらず LBM を読んで脂質を同定する。"""
    from metabolomix.console import method_file as method_file_mod
    from metabolomix.console.job_manager import load_job
    exe = _fake_exe_with_lbm(tmp_path)
    _plan_ready(tmp_path, monkeypatch, exe)
    proj = tmp_path / "proj"
    proj.mkdir()
    (proj / "lab.lbm2").touch()
    method = proj / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: lab.lbm2\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height",
                                      omics="metabolomics"))
    assert parsed["lbm"]["source"] == "method_file"
    keys = method_file_mod.read_method_keys(Path(load_job(Path(parsed["job_path"])).method_file))
    assert Path(keys["lbm file path"]) == (proj / "lab.lbm2").resolve()


def test_console_method_template_writes_relative_paths_as_absolute(tmp_path, monkeypatch):
    """別フォルダへ書き出すので、相対宣言は原本基準の絶対パスに固定する。"""
    from metabolomix.console import method_file as method_file_mod
    monkeypatch.delenv("MSDIAL_EXE", raising=False)
    proj = tmp_path / "proj"
    proj.mkdir()
    (proj / "lab.lbm2").touch()
    src = proj / "params.txt"
    src.write_text("Ion mode: Negative\nLbm file path: lab.lbm2\nMsp file path: lib.msp\n",
                   encoding="ascii")
    out = tmp_path / "elsewhere" / "pos.txt"
    from metabolomix.tools.console_tools import console_method_template
    parsed = _json.loads(console_method_template(out_path=str(out), polarity="positive",
                                                 based_on=str(src)))
    assert parsed["status"] == "written"
    keys = method_file_mod.read_method_keys(out)
    assert Path(keys["lbm file path"]) == (proj / "lab.lbm2").resolve()
    assert Path(keys["msp file path"]) == (proj / "lib.msp").resolve()


def test_console_plan_rejects_a_non_ascii_absolute_path_before_creating_the_job(tmp_path, monkeypatch):
    """Console はメソッドを `Encoding.ASCII` で読む。`?` に化けた絶対パスでは
    ライブラリが見つからず、黙って同定 0 件になる。"""
    exe = _fake_exe_with_lbm(tmp_path)
    _plan_ready(tmp_path, monkeypatch, exe)
    proj = tmp_path / "解析"
    proj.mkdir()
    (proj / "lab.lbm2").touch()
    method = proj / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: lab.lbm2\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "METHOD_ENCODING_UNSUPPORTED"
    assert not (tmp_path / "runs").exists()


def test_console_method_template_rejects_a_non_ascii_absolute_path(tmp_path, monkeypatch):
    monkeypatch.delenv("MSDIAL_EXE", raising=False)
    proj = tmp_path / "解析"
    proj.mkdir()
    (proj / "lab.lbm2").touch()
    src = proj / "params.txt"
    src.write_text("Ion mode: Negative\nLbm file path: lab.lbm2\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_method_template
    parsed = _json.loads(console_method_template(out_path=str(tmp_path / "pos.txt"),
                                                 polarity="positive", based_on=str(src)))
    assert parsed["error"]["code"] == "METHOD_ENCODING_UNSUPPORTED"
