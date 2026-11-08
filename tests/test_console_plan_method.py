"""console_plan のメソッドファイル解決（LBM 自動解決・自動保存パラメータの探索）。

純ロジックは tests/test_console_method_file.py。ここは MCP ツール層の結線を見る。
"""
from __future__ import annotations

import json as _json
from pathlib import Path

from tests.config_fixtures import use_config


def _fake_exe_with_lbm(tmp_path, *names):
    """Console の実行体（[msdial] exe）と同じフォルダに脂質ライブラリを置く（GUI と同じ配置）。"""
    d = tmp_path / "msdial_app"
    d.mkdir(exist_ok=True)
    for n in names:
        (d / n).touch()
    exe = d / "MSDIALCUI.exe"
    exe.touch()
    return exe


def _plan_ready(tmp_path, monkeypatch, exe):
    from metabolomix.core import session_state
    session_state.session = session_state.AnalysisSession()
    use_config(monkeypatch, tmp_path, {"msdial.exe": exe})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
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
    """GUI と同じく Console の実行体と同じフォルダの *.lbm2 を 1 件だけ自動採用する。"""
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


def test_console_plan_shows_only_the_name_of_a_library_lbm(tmp_path, monkeypatch):
    """`[library] lbm` だけで解いた LBM は、console_plan の戻り値にディレクトリを出さない
    （欄 `lbm.path` は残し、中身をファイル名にする。spec 2026-10-07 §3.3）。実効メソッドには
    絶対パスを書く（Console が読むため）。"""
    exe = tmp_path / "app" / "MSDIALCUI.exe"   # 実行体のフォルダに LBM は無い
    exe.parent.mkdir()
    exe.touch()
    _plan_ready(tmp_path, monkeypatch, exe)
    configured = tmp_path / "secret_share" / "lab_lipids.lbm2"
    configured.parent.mkdir()
    configured.write_bytes(b"fake-lab-lbm")
    use_config(monkeypatch, tmp_path, {"library.lbm": configured})
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan

    raw = console_plan(dataset_root=str(tmp_path), method_file=str(method),
                       polarity="negative", measure="peak_height")
    parsed = _json.loads(raw)
    assert parsed["status"] == "planned"
    assert parsed["lbm"] == {"path": "lab_lipids.lbm2", "source": "config_file"}
    assert "secret_share" not in raw
    effective = Path(parsed["method_file"]).read_text(encoding="ascii")
    assert str(configured) in effective


def test_console_plan_lbm_not_found_hides_the_configured_directory(tmp_path, monkeypatch):
    """`[library] lbm` が無いファイルを指すときの LBM_NOT_FOUND もファイル名だけを示す。"""
    exe = tmp_path / "app" / "MSDIALCUI.exe"
    exe.parent.mkdir()
    exe.touch()
    _plan_ready(tmp_path, monkeypatch, exe)
    missing = tmp_path / "secret_share" / "lab_lipids.lbm2"
    use_config(monkeypatch, tmp_path, {"library.lbm": missing})
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan

    raw = console_plan(dataset_root=str(tmp_path), method_file=str(method),
                       polarity="negative", measure="peak_height")
    parsed = _json.loads(raw)
    assert parsed["error"]["code"] == "LBM_NOT_FOUND"
    assert "lab_lipids.lbm2" in parsed["error"]["message"]
    assert "[library] lbm" in parsed["error"]["message"]
    assert "secret_share" not in raw


def test_console_plan_non_ascii_library_lbm_hides_the_directory(tmp_path, monkeypatch):
    """ASCII 以外を含む `[library] lbm` の置き場所も、METHOD_ENCODING_UNSUPPORTED の
    エラー文に出さない（ファイル名だけ。研究室 MSP と同じ）。"""
    exe = tmp_path / "app" / "MSDIALCUI.exe"
    exe.parent.mkdir()
    exe.touch()
    _plan_ready(tmp_path, monkeypatch, exe)
    configured = tmp_path / "secret_共有" / "lab_lipids.lbm2"
    configured.parent.mkdir()
    configured.write_bytes(b"fake-lab-lbm")
    use_config(monkeypatch, tmp_path, {"library.lbm": configured})
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan

    raw = console_plan(dataset_root=str(tmp_path), method_file=str(method),
                       polarity="negative", measure="peak_height")
    parsed = _json.loads(raw)
    assert parsed["error"]["code"] == "METHOD_ENCODING_UNSUPPORTED"
    assert "lab_lipids.lbm2" in parsed["error"]["message"]
    assert "secret_" not in raw


def test_console_lbm_reads_library_lbm(tmp_path, monkeypatch):
    """LBM の設定は `[library] lbm` 1 つ（決めたこと 4）。解決の順は変わらない:
    引数 → メソッドの宣言 → ビルド生成物 → 設定 → 実行体のフォルダ。設定は 5 番目の段の
    実行体のフォルダより先に、引数と宣言より後に効く。"""
    exe = _fake_exe_with_lbm(tmp_path, "installed.lbm2")  # 実行体のフォルダにも候補がある
    _plan_ready(tmp_path, monkeypatch, exe)
    configured = tmp_path / "lab" / "configured.lbm2"
    configured.parent.mkdir()
    configured.touch()
    use_config(monkeypatch, tmp_path, {"library.lbm": configured})
    from metabolomix.tools.console_tools import console_plan

    def plan(method_text, **kwargs):
        method = tmp_path / "params.txt"
        method.write_text(method_text, encoding="ascii")
        return _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                        polarity="negative", measure="peak_height", **kwargs))

    parsed = plan("Ion mode: Negative\nLbm file path: \n")
    assert parsed["status"] == "planned"
    # [library] lbm は研究室の資産。戻り値は置き場所を出さずファイル名だけ（どの経路でも）。
    assert parsed["lbm"] == {"path": "configured.lbm2", "source": "config_file"}

    declared = tmp_path / "declared.lbm2"
    declared.touch()
    assert plan(f"Ion mode: Negative\nLbm file path: {declared}\n")["lbm"]["source"] == "method_file"
    chosen = tmp_path / "chosen.lbm2"
    chosen.touch()
    assert plan("Ion mode: Negative\n", lbm_file=str(chosen))["lbm"]["source"] == "argument"

    use_config(monkeypatch, tmp_path, {"library.lbm": None})
    assert plan("Ion mode: Negative\n")["lbm"]["source"] == "exe_dir"


def test_console_plan_rejects_an_old_msdial_lbm_key(tmp_path, monkeypatch):
    """`[msdial] lbm` だけが残っていたら黙って実行体のフォルダへ落ちず、CONFIG_INVALID で止まる。"""
    exe = _fake_exe_with_lbm(tmp_path, "installed.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    use_config(monkeypatch, tmp_path, {"msdial.lbm": tmp_path / "old.lbm2"})
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "CONFIG_INVALID"
    assert "[library] lbm" in parsed["error"]["message"]
    assert not (tmp_path / "runs").exists()


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
    proj = tmp_path / "解析"
    proj.mkdir()
    (proj / "lab.lbm2").touch()
    src = proj / "params.txt"
    src.write_text("Ion mode: Negative\nLbm file path: lab.lbm2\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_method_template
    parsed = _json.loads(console_method_template(out_path=str(tmp_path / "pos.txt"),
                                                 polarity="positive", based_on=str(src)))
    assert parsed["error"]["code"] == "METHOD_ENCODING_UNSUPPORTED"


# ---------- library_mode="msp_only"（LBM 不使用・研究室 MSP のみ） ----------

import os

from tests.lab_msp_fixtures import write_lab_msp_config


def _effective_keys(path) -> dict[str, list[str]]:
    keys: dict[str, list[str]] = {}
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            keys.setdefault(key.strip().lower(), []).append(value.strip())
    return keys


def _msp_only_plan(tmp_path, monkeypatch, method_text, **kwargs):
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")  # auto なら拾われる LBM
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text(method_text, encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    kwargs.setdefault("polarity", "negative")
    kwargs.setdefault("library_mode", "msp_only")
    return _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                    measure="peak_height", **kwargs))


def test_console_plan_msp_only_blanks_lbm_and_writes_the_lab_msp(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(
        tmp_path, monkeypatch,
        "Ion mode: Negative\nLbm file path: \nMsp file path: C:/public/other.msp\n"
        "Text DB file path: other.txt\n")
    assert parsed["status"] == "planned"
    assert parsed["library_mode"] == "msp_only"
    assert parsed["lbm"] == {"path": None, "source": "disabled"}
    assert parsed["msp"] == {"file": "lab_neg.msp", "source": "config_file"}
    assert parsed["removed_declarations"] == ["Text DB file path"]
    keys = _effective_keys(parsed["method_file"])
    assert keys["lbm file path"] == [""]
    assert keys["msp file path"] == [os.path.abspath(paths["negative"])]
    assert keys["text db file path"] == [""]
    assert "lab_library" not in _json.dumps(parsed, ensure_ascii=False)  # 置き場所を返さない


def test_console_plan_msp_only_picks_the_positive_library_for_positive(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Positive\n", polarity="positive")
    assert parsed["msp"]["file"] == "lab_pos.msp"
    assert _effective_keys(parsed["method_file"])["msp file path"] == [os.path.abspath(paths["positive"])]


def test_console_plan_msp_only_for_metabolomics(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\nLbm file path: C:/x/y.lbm2\n",
                            omics="metabolomics")
    assert parsed["status"] == "planned"
    assert parsed["removed_declarations"] == ["Lbm file path"]
    assert _effective_keys(parsed["method_file"])["lbm file path"] == [""]


def test_console_plan_msp_only_ignores_the_old_env_var(tmp_path, monkeypatch):
    """spec 2026-11-06 §7.4: 資産の場所の環境変数は読まない（以前は設定ファイルに勝った）。"""
    write_lab_msp_config(tmp_path, monkeypatch)
    env_lib = tmp_path / "env_dir" / "env_neg.msp"
    env_lib.parent.mkdir()
    env_lib.write_text("NAME: e\n", encoding="ascii")
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(env_lib))
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["msp"] == {"file": "lab_neg.msp", "source": "config_file"}


def test_console_plan_msp_only_argument_wins(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    given = tmp_path / "given.msp"
    given.write_text("NAME: g\n", encoding="ascii")
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n", msp_file=str(given))
    assert parsed["msp"] == {"file": "given.msp", "source": "argument"}


def test_console_plan_msp_only_stops_when_not_configured(tmp_path, monkeypatch):
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["error"]["code"] == "MSP_NOT_CONFIGURED"
    assert not (tmp_path / "runs").exists()  # ジョブを作らない


def test_console_plan_msp_only_stops_when_the_file_is_missing(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    paths["negative"].unlink()
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["error"]["code"] == "MSP_NOT_FOUND"
    assert "secret_share" not in _json.dumps(parsed, ensure_ascii=False)
    assert not (tmp_path / "runs").exists()


def test_console_plan_msp_only_rejects_lbm_file(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n", lbm_file="C:/x/y.lbm2")
    assert parsed["error"]["code"] == "LIBRARY_MODE_CONFLICT"


def test_console_plan_auto_rejects_msp_file(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n",
                            library_mode="auto", msp_file="C:/x/y.msp")
    assert parsed["error"]["code"] == "LIBRARY_MODE_CONFLICT"


def test_console_plan_rejects_an_unknown_library_mode(tmp_path, monkeypatch):
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n", library_mode="lbm_only")
    assert parsed["error"]["code"] == "JOB_NOT_PLANNED"


def test_console_plan_msp_only_rejects_a_non_ascii_library_path(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="非公開の置き場")
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["error"]["code"] == "METHOD_ENCODING_UNSUPPORTED"
    assert "非公開の置き場" not in _json.dumps(parsed, ensure_ascii=False)
    assert not (tmp_path / "runs").exists()


def test_console_plan_auto_reports_msp_not_used(tmp_path, monkeypatch):
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: \n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["library_mode"] == "auto"
    assert parsed["lbm"]["source"] == "exe_dir"
    assert parsed["msp"] == {"file": None, "source": "not_used"}
    assert parsed["removed_declarations"] == []


def test_console_method_template_msp_only(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    src = tmp_path / "neg_param.txt"
    src.write_text("Ion mode: Negative\nLbm file path: C:/x/y.lbm2\nText DB file path: db.txt\n",
                   encoding="ascii")
    out = tmp_path / "pos_param.txt"
    from metabolomix.tools.console_tools import console_method_template
    parsed = _json.loads(console_method_template(
        out_path=str(out), polarity="positive", based_on=str(src), library_mode="msp_only"))
    assert parsed["status"] == "written"
    assert parsed["lbm"] == {"path": None, "source": "disabled"}
    assert parsed["msp"] == {"file": "lab_pos.msp", "source": "config_file"}
    assert parsed["removed_declarations"] == ["Lbm file path", "Text DB file path"]
    keys = _effective_keys(out)
    assert keys["lbm file path"] == [""]
    assert keys["msp file path"] == [os.path.abspath(paths["positive"])]
    assert keys["ion mode"] == ["Positive"]


def test_console_method_template_msp_only_next_hint_keeps_the_mode(tmp_path, monkeypatch):
    """最終レビュー #7: 雛形の次の一手が auto の console_plan を勧めない（LBM が自動補完され契約が破れる）。"""
    write_lab_msp_config(tmp_path, monkeypatch)
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    src = tmp_path / "neg_param.txt"
    src.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_method_template
    parsed = _json.loads(console_method_template(
        out_path=str(tmp_path / "pos_param.txt"), polarity="positive", based_on=str(src),
        library_mode="msp_only"))
    assert "library_mode='msp_only'" in parsed["next"]
    auto = _json.loads(console_method_template(
        out_path=str(tmp_path / "pos_auto.txt"), polarity="positive", based_on=str(src)))
    assert "library_mode" not in auto["next"]
