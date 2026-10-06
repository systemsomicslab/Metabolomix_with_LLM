"""入力一覧・メソッド選択・解析専用配置（spec §4）を検証する。

対象は metabolomix.pipeline.inputs の inspect_inputs / select_method /
stage_inputs / verify_inputs。実rawや既存成果物は使わず、
tests.pipeline_fixtures.make_source が作る合成データだけを使う。

fake.exe は実行しない（is_console_exe を差し替える）。fake.lbm2 も
Console には読ませない合成バイト列。
"""
from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from metabolomix.console import method_file as method_file_mod
from metabolomix.core import app_control
from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline.inputs import inspect_inputs, select_method, stage_inputs, verify_inputs
from metabolomix.pipeline.request import resolve_request
from tests.pipeline_fixtures import make_source


def _allow_fake_exe(monkeypatch) -> None:
    """fake.exe を実際には実行させず、Console 実行体として通す。"""
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)


def _method_lbm_exe(root: Path, *, ion_mode: str = "negative") -> dict:
    """S0..S7.wiff を持たない、method/lbm/exeだけの最小フォルダを作る。

    ディレクトリraw・混在形式などmake_sourceの形と違う raw を使うテスト用。
    """
    root.mkdir(parents=True, exist_ok=True)
    lbm = root / "fake.lbm2"
    lbm.write_bytes(b"fake-lbm-library")
    method = root / "lab_param_202609050001.txt"
    method.write_text(
        f"Ion mode: {ion_mode}\nTarget omics: Lipidomics\nLbm file path: fake.lbm2\n",
        encoding="ascii", newline="\n")
    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")
    return {"root": root, "method": method, "lbm": lbm, "exe": exe}


def _snapshot_tree(root: Path) -> dict[str, tuple[int, int]]:
    """フォルダ直下ファイルの (size, mtime_ns) を集めた辞書（改変検出用）。"""
    return {p.name: (p.stat().st_size, p.stat().st_mtime_ns)
            for p in sorted(root.iterdir()) if p.is_file()}


def _make_junction(link_path: Path, target: Path) -> None:
    result = subprocess.run(
        ["cmd", "/c", "mklink", "/J", str(link_path), str(target)],
        capture_output=True)
    assert result.returncode == 0, f"mklink /J に失敗しました: {result.returncode}"


# ---------- brief記載のRED ----------

def test_distinct_methods_are_not_selected_by_mtime():
    candidates = [
        {"path": "own.txt", "sha256": "a" * 64, "mtime": 1, "polarity": "negative"},
        {"path": "sibling.txt", "sha256": "b" * 64, "mtime": 2, "polarity": "negative"},
    ]
    with pytest.raises(DomainError, match="METHOD_FILE_CHOICE_REQUIRED"):
        select_method(candidates)


def test_byte_identical_copies_are_one_method():
    candidates = [{"path": name, "sha256": "a" * 64, "polarity": "negative"}
                  for name in ["own.txt", "copy.txt"]]
    assert select_method(candidates)["sha256"] == "a" * 64


# ---------- 基本の一気通貫 ----------

def test_inspect_and_stage_the_default_source(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])

    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    assert plan["selected_format"] == "wiff"
    assert plan["polarity"] == {"value": "negative", "source": "method_declaration"}
    assert "polarity_from_method_declaration_not_verified_from_raw" in plan["unverified"]
    assert plan["method"]["source_path"] == str(src["method"].resolve())
    # 元メソッドの `Lbm file path: fake.lbm2` は相対参照（原本フォルダ基準で解決可能）。
    # 実効コピーは原本と別ディレクトリに置かれるため、コピー後も相対参照の意味を
    # 変えない（spec §4.3）よう、絶対パスへ書き換える上書きが必要になる
    # （空dictにはならない）。
    assert plan["method"]["overrides"] == {method_file_mod.LBM_KEY: str(src["lbm"].resolve())}

    pipeline_root = tmp_path / "pipeline_run"
    snapshot = stage_inputs(plan, pipeline_root)

    input_dir = pipeline_root / "input"
    for i in range(8):
        assert (input_dir / f"S{i}.wiff").is_file()
    effective = pipeline_root / "inputs" / "effective-method.txt"
    assert effective.is_file()
    assert snapshot["method"]["effective_relative_path"] == "inputs/effective-method.txt"
    assert snapshot["method"]["effective_sha256"]

    verify_inputs(snapshot)  # 直後の再検査は何も検出しない


def test_effective_method_lbm_reference_resolves_from_its_own_directory(tmp_path, monkeypatch):
    """相対 `Lbm file path` を実効コピーへ verbatim コピーすると、コピー先の
    別ディレクトリでは解決できない（spec §4.3 違反）。この後退を検出するのは
    `overrides` dict のアサートだけでは不十分——TDDでこのバグを見逃した原因が
    まさにそれ（overrideの有無しか見ておらず、実効ファイルの中身を読み戻して
    いなかった）。ここでは実際に書き出された `effective-method.txt` を読み戻し、
    `Lbm file path` が**絶対パス**で、実在するファイルを指すことを検証する。

    LC-MS の Console（`ConfigParser.ReadForLcmsParameter`）は宣言パスを解決
    しないので、相対値はメソッドファイル基準ではなく Console プロセスの cwd
    基準になる（メソッド基準で解くのは GC-MS 経路の `ResolvePathFromMethodFile`
    だけ。MsdialWorkbench afd5f95）。相対のまま残すと、どの基準で読んでも
    原本の意図とずれる。
    """
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])

    pipeline_root = tmp_path / "pipeline_run"
    stage_inputs(plan, pipeline_root)

    effective = pipeline_root / "inputs" / "effective-method.txt"
    keys = method_file_mod.read_method_keys(effective)
    declared = keys.get(method_file_mod.LBM_KEY.lower())
    assert declared, "実効メソッドに Lbm file path が書かれていない"

    candidate = Path(declared)
    assert candidate.is_absolute(), f"実効コピーの Lbm file path が相対のまま: {declared!r}"
    assert candidate.is_file(), (
        f"実効コピーの Lbm file path が解決できません: 宣言値={declared!r}")


def test_effective_method_writes_relative_msp_and_text_db_as_absolute_paths(tmp_path, monkeypatch):
    """LBM 以外のライブラリ宣言も、相対のままでは Console の cwd 基準で読まれる。

    既存 lipidomics 経路の厳格な参照検査（`REFERENCE_KEYS`＝LBM のみ）は広げない
    ——実在しない古い宣言でも止めず、原本基準の絶対パスに固定するだけ。
    """
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    (src["root"] / "lib.msp").write_text("NAME: x\n", encoding="ascii")
    method = src["method"]
    method.write_text(method.read_text(encoding="ascii")
                      + "Msp file path: lib.msp\nText DB file path: gone.txt\n",
                      encoding="ascii", newline="\n")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])

    pipeline_root = tmp_path / "pipeline_run"
    stage_inputs(plan, pipeline_root)
    keys = method_file_mod.read_method_keys(pipeline_root / "inputs" / "effective-method.txt")
    assert Path(keys["msp file path"]) == (src["root"] / "lib.msp").resolve()
    assert Path(keys["text db file path"]) == (src["root"] / "gone.txt").resolve()


def test_source_is_byte_and_stat_identical_after_staging(tmp_path, monkeypatch):
    """束縛条件: 元rawフォルダは移動・上書き・削除されない（バイト同一・stat同一）。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    before_bytes = {p.name: p.read_bytes() for p in sorted(src["root"].iterdir()) if p.is_file()}
    before_stat = _snapshot_tree(src["root"])

    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    stage_inputs(plan, tmp_path / "pipeline_run")

    after_bytes = {p.name: p.read_bytes() for p in sorted(src["root"].iterdir()) if p.is_file()}
    after_stat = _snapshot_tree(src["root"])
    assert after_bytes == before_bytes
    assert after_stat == before_stat
    assert sorted(p.name for p in src["root"].iterdir()) == sorted(before_stat)


def test_uses_hardlinks_on_the_same_volume(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"
    stage_inputs(plan, pipeline_root)

    original = src["root"] / "S0.wiff"
    staged = pipeline_root / "input" / "S0.wiff"
    assert original.stat().st_ino == staged.stat().st_ino


def test_falls_back_to_copy_when_link_fails(tmp_path, monkeypatch):
    """同一ボリュームでリンクできない状況（他プロセス・別ボリューム等）ではコピーする。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"

    def _always_fail(_src, _dst):
        raise OSError("同一ボリュームではない、と仮定する")

    stage_inputs(plan, pipeline_root, link_fn=_always_fail)
    original = src["root"] / "S0.wiff"
    staged = pipeline_root / "input" / "S0.wiff"
    assert original.stat().st_ino != staged.stat().st_ino
    assert staged.read_text(encoding="ascii") == original.read_text(encoding="ascii")


# ---------- 出力先のバリエーション ----------

def test_different_output_destination(tmp_path, monkeypatch):
    """出力先を元フォルダと無関係な場所に変えても問題なく配置できる。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])

    other_root = tmp_path / "elsewhere" / "not_related_to_source" / "pipeline_run"
    snapshot = stage_inputs(plan, other_root)
    assert Path(snapshot["pipeline_root"]) == other_root.resolve()
    assert (other_root / "input" / "S0.wiff").is_file()


def test_readonly_source_with_separate_output_root(tmp_path, monkeypatch):
    """元フォルダが読み取り専用でも、別のoutput_rootへ配置できる（書込隔離の確認）。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    for p in src["root"].iterdir():
        if p.is_file():
            os.chmod(p, stat.S_IREAD)
    try:
        request = resolve_request(src["root"])
        plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
        pipeline_root = tmp_path / "separate_output_root" / "pipeline_run"
        snapshot = stage_inputs(plan, pipeline_root)
        assert (pipeline_root / "input" / "S0.wiff").is_file()
        verify_inputs(snapshot)
    finally:
        for p in src["root"].iterdir():
            if p.is_file():
                os.chmod(p, stat.S_IWRITE | stat.S_IREAD)


# ---------- D07: 同内容の複製 と 異内容の複数method ----------

def test_byte_identical_method_duplicates_are_aggregated(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    # 兄弟フォルダに全く同じ内容のメソッドファイルを置く（method_search_dirsが拾う）。
    sibling = tmp_path / "sibling_project"
    sibling.mkdir()
    duplicate = sibling / "lab_param_202609060002.txt"
    duplicate.write_bytes(src["method"].read_bytes())
    # 複製のLbm file pathは相対参照なので、複製先にも同名・同内容のlbm2を置く
    # （参照先ハッシュも一致させ、「同一メソッド」として集約されることを見る）。
    (sibling / "fake.lbm2").write_bytes(src["lbm"].read_bytes())

    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    assert plan["method"]["source_path"] in (str(src["method"].resolve()), str(duplicate.resolve()))


def test_distinct_method_contents_require_choice(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    sibling = tmp_path / "sibling_project"
    sibling.mkdir()
    # Ion modeが違う=内容が違うメソッド。極性を明示しなければどちらも候補に残る。
    _method_lbm_exe(sibling, ion_mode="negative")
    # 上のヘルパーはfake.exe/fake.lbm2も作るが、method選択には影響しない。
    other_method = sibling / "lab_param_202609060003.txt"
    other_method.write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\nMinimum peak height: 999\n",
        encoding="ascii", newline="\n")

    request = resolve_request(src["root"])
    with pytest.raises(DomainError, match="METHOD_FILE_CHOICE_REQUIRED"):
        inspect_inputs(src["root"], request, exe_path=src["exe"])


def test_differing_reference_path_bases_are_not_one_method(tmp_path, monkeypatch):
    """バイト同一のメソッドでも、宣言参照(Lbm file path)の解決先が違えば別扱い。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    sibling = tmp_path / "sibling_project"
    sibling.mkdir()
    duplicate = sibling / "lab_param_202609060004.txt"
    duplicate.write_bytes(src["method"].read_bytes())  # 同じ相対参照文字列
    # sibling側のfake.lbm2は「存在するが中身が違う」——参照文字列は同じでも
    # 解決先の実効内容が異なる。
    (sibling / "fake.lbm2").write_bytes(b"different-lbm-content")

    request = resolve_request(src["root"])
    with pytest.raises(DomainError, match="METHOD_FILE_CHOICE_REQUIRED"):
        inspect_inputs(src["root"], request, exe_path=src["exe"])


# ---------- 随伴ファイル: wiff.scan保持・wiff2除外・PAI2誤混入防止 ----------

def test_wiff_scan_is_kept_and_wiff2_and_pai2_are_excluded(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "raw"
    fixture = _method_lbm_exe(root)
    (root / "S0.wiff").write_text("primary", encoding="ascii")
    (root / "S0.wiff.scan").write_text("scan", encoding="ascii")
    (root / "S0.wiff2").write_text("other-format", encoding="ascii")
    (root / "S0.timeseries.data").write_text("meta", encoding="ascii")
    # 既存の解析出力が同じstemで残っているケース（誤って随伴に化けてはいけない）。
    (root / "S0.pai2").write_bytes(b"\x00pai2-binary")
    (root / "S0_tags.xml").write_text("<tags/>", encoding="ascii")

    request = resolve_request(root, {"keep_extension": "wiff"})
    plan = inspect_inputs(root, request, exe_path=fixture["exe"])
    assert plan["companions"]["S0.wiff"] == ["S0.wiff.scan", "S0.timeseries.data"]

    snapshot = stage_inputs(plan, tmp_path / "pipeline_run")
    staged_names = {p.name for p in (Path(snapshot["pipeline_root"]) / "input").iterdir()}
    assert staged_names == {"S0.wiff", "S0.wiff.scan", "S0.timeseries.data"}
    assert "S0.wiff2" not in staged_names
    assert "S0.pai2" not in staged_names
    assert "S0_tags.xml" not in staged_names


# ---------- D08: ディレクトリraw・配置済み不一致・リンク不可 ----------

def test_directory_style_raw_is_copied_wholesale(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "raw"
    fixture = _method_lbm_exe(root)
    dir_raw = root / "sample_a.d"
    dir_raw.mkdir()
    (dir_raw / "chromatogram.raw").write_bytes(b"agilent-inner-file")
    (dir_raw / "meta.xml").write_text("<meta/>", encoding="ascii")

    request = resolve_request(root, {"keep_extension": "d"})
    plan = inspect_inputs(root, request, exe_path=fixture["exe"])
    assert plan["selected_format"] == "d"
    assert any(e["name"] == "sample_a.d" and e["kind"] == "dir" for e in plan["entries"])

    snapshot = stage_inputs(plan, tmp_path / "pipeline_run")
    staged_dir = Path(snapshot["pipeline_root"]) / "input" / "sample_a.d"
    assert (staged_dir / "chromatogram.raw").read_bytes() == b"agilent-inner-file"
    assert (staged_dir / "meta.xml").is_file()
    # 元ディレクトリは無変更。
    assert (dir_raw / "chromatogram.raw").read_bytes() == b"agilent-inner-file"

    verify_inputs(snapshot)


def test_staged_mismatch_is_detected_on_resume(tmp_path, monkeypatch):
    """配置済みファイルが元入力と食い違えば、再開時に単純skipしない。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"

    def _always_fail(_src, _dst):
        raise OSError("copy2経路を強制するため常に失敗させる")

    stage_inputs(plan, pipeline_root, link_fn=_always_fail)  # 独立コピーとして配置

    staged = pipeline_root / "input" / "S0.wiff"
    staged.write_text("corrupted-after-staging", encoding="ascii")  # 配置済み側だけ改変

    with pytest.raises(DomainError, match="STAGED_INPUT_MISMATCH"):
        stage_inputs(plan, pipeline_root, link_fn=_always_fail)


def test_resume_fills_in_missing_staged_file(tmp_path, monkeypatch):
    """再開時、配置済み側から欠けているファイルだけを補って配置する。"""
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"
    stage_inputs(plan, pipeline_root)

    (pipeline_root / "input" / "S1.wiff").unlink()
    snapshot = stage_inputs(plan, pipeline_root)
    assert (pipeline_root / "input" / "S1.wiff").is_file()
    verify_inputs(snapshot)


def test_unexpected_extra_file_in_destination_is_rejected(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"
    (pipeline_root / "input").mkdir(parents=True)
    (pipeline_root / "input" / "unexpected_leftover.wiff").write_text("x", encoding="ascii")

    with pytest.raises(DomainError, match="STAGED_INPUT_MISMATCH"):
        stage_inputs(plan, pipeline_root)


def test_insufficient_free_space_stops_before_copying(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"

    def _tiny_free_space(_path):
        return SimpleNamespace(total=1, used=0, free=1)

    with pytest.raises(DomainError, match="INSUFFICIENT_FREE_SPACE"):
        stage_inputs(plan, pipeline_root, disk_usage=_tiny_free_space)
    assert not (pipeline_root / "input").exists() or not any((pipeline_root / "input").iterdir())


def test_reparse_point_escaping_source_root_is_rejected(tmp_path, monkeypatch):
    """ディレクトリraw内のジャンクションが元root外を指す場合は拒否する。"""
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "raw"
    fixture = _method_lbm_exe(root)
    outside = tmp_path / "outside_the_source_root"
    outside.mkdir()
    (outside / "secret.txt").write_text("must not be copied in", encoding="ascii")

    trap = root / "trap.d"
    trap.mkdir()
    _make_junction(trap / "escape", outside)

    request = resolve_request(root, {"keep_extension": "d"})
    with pytest.raises(DomainError, match="INPUT_ESCAPES_SOURCE_ROOT"):
        inspect_inputs(root, request, exe_path=fixture["exe"])


# ---------- raw stat変化 / 元method不変 ----------

def test_verify_inputs_detects_changed_raw_stat(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    snapshot = stage_inputs(plan, tmp_path / "pipeline_run")

    (src["root"] / "S0.wiff").write_text("changed-after-staging-with-more-bytes",
                                         encoding="ascii")
    with pytest.raises(DomainError, match="INPUT_CHANGED"):
        verify_inputs(snapshot)


def test_verify_inputs_detects_changed_method_content(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    snapshot = stage_inputs(plan, tmp_path / "pipeline_run")

    src["method"].write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\nLbm file path: fake.lbm2\n"
        "Minimum peak height: 12345\n",
        encoding="ascii", newline="\n")
    with pytest.raises(DomainError, match="INPUT_CHANGED"):
        verify_inputs(snapshot)


def test_original_method_file_is_never_modified(tmp_path, monkeypatch):
    """LBM解決でoverrideが要る場合でも、元メソッドファイルは書き換えない。"""
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "raw"
    root.mkdir()
    for i in range(2):
        (root / f"S{i}.wiff").write_text(f"raw-{i}", encoding="ascii")
    lbm_dir = tmp_path / "lbm_elsewhere"
    lbm_dir.mkdir()
    lbm = lbm_dir / "elsewhere.lbm2"
    lbm.write_bytes(b"lbm-in-a-different-folder")
    method = root / "lab_param_202609070005.txt"
    # Lbm file path を宣言しない=resolve_lbmがlbm_file引数へフォールバックし、
    # 実効メソッドにoverrideが書かれる（元メソッドは変わらないことを検証する）。
    method.write_text("Ion mode: negative\nTarget omics: Lipidomics\n",
                      encoding="ascii", newline="\n")
    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")
    before = method.read_bytes()

    request = resolve_request(root, {"lbm_file": str(lbm)})
    plan = inspect_inputs(root, request, exe_path=exe)
    assert plan["method"]["overrides"]  # 実際にoverrideが発生するケース
    stage_inputs(plan, tmp_path / "pipeline_run")

    assert method.read_bytes() == before


# ---------- METHOD_REFERENCE_UNRESOLVED / METHOD_ENCODING_UNSUPPORTED ----------

def test_unresolvable_declared_reference_stops_before_running(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "raw"
    root.mkdir()
    (root / "S0.wiff").write_text("raw-0", encoding="ascii")
    method = root / "lab_param_202609080006.txt"
    method.write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\n"
        "Lbm file path: this_file_does_not_exist.lbm2\n",
        encoding="ascii", newline="\n")
    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")

    request = resolve_request(root)
    with pytest.raises(DomainError, match="METHOD_REFERENCE_UNRESOLVED"):
        inspect_inputs(root, request, exe_path=exe)


def test_non_ascii_override_value_is_rejected_not_substituted(tmp_path, monkeypatch):
    """実効メソッドへ書けない文字は`?`置換せずMETHOD_ENCODING_UNSUPPORTEDにする。"""
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "raw"
    root.mkdir()
    (root / "S0.wiff").write_text("raw-0", encoding="ascii")
    method = root / "lab_param_202609090007.txt"
    method.write_text("Ion mode: negative\nTarget omics: Lipidomics\n",
                      encoding="ascii", newline="\n")
    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")

    lbm_dir = tmp_path / "ライブラリ"  # 非ASCIIパス
    lbm_dir.mkdir()
    lbm = lbm_dir / "fake.lbm2"
    lbm.write_bytes(b"lbm-under-non-ascii-directory")

    request = resolve_request(root, {"lbm_file": str(lbm)})
    plan = inspect_inputs(root, request, exe_path=exe)
    with pytest.raises(DomainError, match="METHOD_ENCODING_UNSUPPORTED"):
        stage_inputs(plan, tmp_path / "pipeline_run")


def test_relative_lbm_reference_resolving_to_non_ascii_path_is_rejected(tmp_path, monkeypatch):
    """相対 `Lbm file path` を絶対化した結果が非ASCIIパスになる場合も、
    `?`置換ではなくMETHOD_ENCODING_UNSUPPORTEDで止まる（宣言済み参照を絶対化
    する今回の修正が新しく通す書き込み経路であり、既存の非ASCIIガードが
    ここにも及ぶことを確認する）。
    """
    _allow_fake_exe(monkeypatch)
    root = tmp_path / "ライブラリ" / "raw"  # 原本フォルダ自体が非ASCIIパス配下
    root.mkdir(parents=True)
    (root / "S0.wiff").write_text("raw-0", encoding="ascii")
    lbm = root / "fake.lbm2"
    lbm.write_bytes(b"fake-lbm-library")
    method = root / "lab_param_202609100008.txt"
    # 宣言自体はASCII（相対パス）だが、原本ディレクトリが非ASCIIのため
    # 絶対解決した結果は非ASCIIパスになる。
    method.write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\nLbm file path: fake.lbm2\n",
        encoding="ascii", newline="\n")
    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")

    request = resolve_request(root)
    plan = inspect_inputs(root, request, exe_path=exe)
    assert plan["method"]["overrides"]  # 相対宣言なので絶対パスへの上書きが発生する
    with pytest.raises(DomainError, match="METHOD_ENCODING_UNSUPPORTED"):
        stage_inputs(plan, tmp_path / "pipeline_run")


# ---------- 上流が入力フォルダへ書く生成物（2026-09-17 Stage B の #1）----------

def test_restaging_tolerates_the_artifacts_msdial_writes_into_the_input_folder(
        tmp_path, monkeypatch):
    """`rerun_upstream` の再配置が、MS-DIAL 自身の生成物で塞がらないこと。

    MS-DIAL は `-o` だけでなく **`-i` 側**にも `.arf` / `.pai2` / `.dcl` /
    `_tags.xml` を書く（`console/execution.py` が両ルートを snapshot するのは
    そのため）。pipeline では `-i` が staged input そのものなので、Console が
    一度でも走れば必ずこれらが入力フォルダに現れる。これを「想定外のファイル」
    として拒むと、やり直しの唯一の正規手順が使えなくなる。
    """
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"
    stage_inputs(plan, pipeline_root)

    input_dir = pipeline_root / "input"
    for name in ("AlignResult-20260917.arf", "AlignResult-20260917.arf2",
                 "S0_20260917.pai2", "S0_20260917.dcl", "S0_20260917_tags.xml",
                 "AlignResult-20260917.EIC.aef"):
        (input_dir / name).write_text("upstream artifact", encoding="ascii")

    stage_inputs(plan, pipeline_root)   # 例外を投げないこと


def test_restaging_still_refuses_a_file_it_cannot_account_for(tmp_path, monkeypatch):
    """上流生成物と判定できないファイルは従来どおり拒む。

    「入力フォルダに知らないものがある」を一律に許すと、別バッチの生データや
    取り違えたコピーが紛れ込んでも気付けなくなる。
    """
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"])
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    pipeline_root = tmp_path / "pipeline_run"
    stage_inputs(plan, pipeline_root)

    (pipeline_root / "input" / "stray-sample.wiff").write_text("よそのバッチ", encoding="utf-8")

    with pytest.raises(DomainError, match="STAGED_INPUT_MISMATCH"):
        stage_inputs(plan, pipeline_root)


# ---------- プロジェクト保存が Application Control に塞がれる場合（2026-09-17 #2）----------

def _unsigned_dll(path):
    """Certificate Table が空の最小 PE を書く（＝未署名）。"""
    import struct
    pe_off = 0x80
    b = bytearray(0x400)
    b[0:2] = b"MZ"
    struct.pack_into("<I", b, 0x3C, pe_off)
    b[pe_off:pe_off + 4] = b"PE\x00\x00"
    struct.pack_into("<H", b, pe_off + 4, 0x8664)
    struct.pack_into("<H", b, pe_off + 0x18, 0x20B)
    dd = pe_off + 0x18 + 0x70
    struct.pack_into("<II", b, dd + 4 * 8, 0, 0)
    path.write_bytes(bytes(b))
    return path


def test_planning_refuses_project_save_when_app_control_would_block_it(
        tmp_path, monkeypatch):
    """Console を起動する**前**に止める。

    Smart App Control が有効なマシンでは、プロジェクト保存に要る未署名
    アセンブリの読み込みが弾かれる。しかもそれが起きるのは MS-DIAL が全検体の
    解析とアライメントを終えた**後**で、実測では 13 分半を費やしてから落ちた。
    計画時に判定すれば 1 秒で分かる。
    """
    _allow_fake_exe(monkeypatch)
    monkeypatch.setattr(app_control, "smart_app_control_state", lambda: 1)
    src = make_source(tmp_path / "raw")
    _unsigned_dll(Path(src["exe"]).parent / "MsdialLcImMsApi.dll")
    request = resolve_request(src["root"])
    request["save_project"] = True

    with pytest.raises(DomainError) as exc:
        inspect_inputs(src["root"], request, exe_path=src["exe"])
    assert exc.value.code == "PROJECT_SAVE_BLOCKED"
    # 対処を封筒に載せる（クライアントが次の一手を推測しなくて済むように）。
    assert "save_project" in exc.value.message


def test_planning_proceeds_when_the_project_is_not_requested(tmp_path, monkeypatch):
    """`save_project=false` なら塞がれる経路に入らないので止めない。"""
    _allow_fake_exe(monkeypatch)
    monkeypatch.setattr(app_control, "smart_app_control_state", lambda: 1)
    src = make_source(tmp_path / "raw")
    _unsigned_dll(Path(src["exe"]).parent / "MsdialLcImMsApi.dll")
    request = resolve_request(src["root"])
    request["save_project"] = False

    inspect_inputs(src["root"], request, exe_path=src["exe"])   # 例外を投げない


def test_planning_proceeds_when_the_policy_is_not_enforcing(tmp_path, monkeypatch):
    """ポリシーが無効なら未署名でも止めない（ポリシー単独で塞がない）。"""
    _allow_fake_exe(monkeypatch)
    monkeypatch.setattr(app_control, "smart_app_control_state", lambda: 0)
    src = make_source(tmp_path / "raw")
    _unsigned_dll(Path(src["exe"]).parent / "MsdialLcImMsApi.dll")
    request = resolve_request(src["root"])
    request["save_project"] = True

    inspect_inputs(src["root"], request, exe_path=src["exe"])   # 例外を投げない
