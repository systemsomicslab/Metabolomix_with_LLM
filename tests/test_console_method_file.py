"""MS-DIAL Console のメソッドファイル解決（LBM 自動解決・既存パラメータ探索）。

MS-DIAL GUI は脂質ライブラリ（`*.lbm2`）をアプリフォルダから自動で拾い、
ユーザーには選ばせない（`MethodSettingModelFactory.cs` /
`DataBaseSettingModel.TrySetLbmLibrary`）。一方 Console はメソッドファイルの
`Lbm file path:` しか見ない（`CommonProcess.cs`）。そして GUI が実行のたびに
自動保存するパラメータ（`MethodModelBase.AutoParametersSave`）は、GUI が
`param.LbmFilePath` に一度も代入しないため**必ず空**で出る。

ここはその段差を埋める層のテスト。GUI と同じ規則（アプリフォルダにちょうど 1 件）
を再現し、GUI が拒否する状況でだけ拒否する。
"""
from __future__ import annotations

import time
from pathlib import Path

import pytest

from metabolomix.console.method_file import (
    KEY_PARAMS_MAX_CANDIDATES,
    LbmResolution,
    MethodCandidate,
    discover_method_candidates,
    extract_key_params,
    find_lbm_files,
    method_search_dirs,
    past_run_method_files,
    read_method_keys,
    relative_path_overrides,
    resolve_lbm,
    scan_dir_for_method_files,
    write_effective_method_file,
)
from metabolomix.core.user_config import Setting


def _env_lbm(value) -> Setting:
    """環境変数 MSDIAL_LBM から来た値（resolve_lbm は出どころごと受け取る）。"""
    return Setting(key="msdial.lbm", value=str(value), source="env",
                   env_var="MSDIAL_LBM", config_file=None)


# ---------- read_method_keys ----------

def test_read_method_keys_parses_colon_pairs(tmp_path):
    p = tmp_path / "param.txt"
    p.write_text("Ion mode: Positive\nTarget omics: Lipidomics\n", encoding="ascii")
    keys = read_method_keys(p)
    assert keys["ion mode"] == "Positive"
    assert keys["target omics"] == "Lipidomics"


def test_read_method_keys_skips_comments_and_blanks(tmp_path):
    p = tmp_path / "param.txt"
    p.write_text("# header\n\nIon mode: Negative\n", encoding="ascii")
    assert read_method_keys(p) == {"ion mode": "Negative"}


def test_read_method_keys_splits_on_first_delimiter(tmp_path):
    """ConfigParser は最初の ':' か '=' のうち先に来るほうで割る。"""
    p = tmp_path / "param.txt"
    p.write_text("File ID=0: 1\n", encoding="ascii")
    assert read_method_keys(p) == {"file id": "0: 1"}


def test_read_method_keys_keeps_empty_value(tmp_path):
    """`Lbm file path:` が空で存在する、が本件の中心なので落としてはいけない。"""
    p = tmp_path / "param.txt"
    p.write_text("Lbm file path: \nIon mode: Positive\n", encoding="ascii")
    keys = read_method_keys(p)
    assert "lbm file path" in keys
    assert keys["lbm file path"] == ""


def test_read_method_keys_preserves_windows_path_value(tmp_path):
    """値側の ':' で割ってはいけない（C:\\... が壊れる）。"""
    p = tmp_path / "param.txt"
    p.write_text("Lbm file path: C:\\lib\\x.lbm2\n", encoding="ascii")
    assert read_method_keys(p)["lbm file path"] == "C:\\lib\\x.lbm2"


# ---------- find_lbm_files ----------

def test_find_lbm_files_matches_lbm_and_lbm2(tmp_path):
    (tmp_path / "a.lbm").touch()
    (tmp_path / "b.lbm2").touch()
    (tmp_path / "c.txt").touch()
    (tmp_path / "d.lbmx").touch()
    found = {p.name for p in find_lbm_files(tmp_path)}
    assert found == {"a.lbm", "b.lbm2"}


def test_find_lbm_files_is_top_directory_only(tmp_path):
    """GUI は SearchOption.TopDirectoryOnly。"""
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "deep.lbm2").touch()
    assert find_lbm_files(tmp_path) == []


def test_find_lbm_files_missing_directory_is_empty(tmp_path):
    assert find_lbm_files(tmp_path / "nope") == []


# ---------- resolve_lbm ----------

def _exe_dir_with(tmp_path, *names):
    d = tmp_path / "msdial"
    d.mkdir()
    for n in names:
        (d / n).touch()
    return d / "MSDIALCUI.exe"


def test_resolve_lbm_prefers_explicit_method_file_value(tmp_path):
    lib = tmp_path / "explicit.lbm2"
    lib.touch()
    exe = _exe_dir_with(tmp_path, "other.lbm2")
    res = resolve_lbm({"lbm file path": str(lib)}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.source == "method_file"
    assert Path(res.path) == lib
    assert res.error_code is None


def test_resolve_lbm_resolves_relative_value_against_method_file(tmp_path):
    """相対宣言はメソッドファイルの親基準で解く。

    LC-MS の Console（`ConfigParser.ReadForLcmsParameter`）はパスを解決せず、
    相対値はプロセスの cwd 基準になる。メソッド基準で解くのは GC-MS 経路の
    `ResolvePathFromMethodFile` だけ。ここでメソッド基準に解いた結果は、
    呼び出し側が実効メソッドへ**絶対パスで**書き戻して初めて Console に届く。
    """
    lib = tmp_path / "lib.lbm2"
    lib.touch()
    exe = _exe_dir_with(tmp_path, "other.lbm2")
    res = resolve_lbm({"lbm file path": "lib.lbm2"}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert Path(res.path) == lib
    assert res.source == "method_file"


def test_resolve_lbm_errors_when_declared_path_is_missing(tmp_path):
    exe = _exe_dir_with(tmp_path, "other.lbm2")
    res = resolve_lbm({"lbm file path": str(tmp_path / "gone.lbm2")}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.error_code == "LBM_NOT_FOUND"
    assert res.path is None


def test_resolve_lbm_falls_back_to_env(tmp_path):
    lib = tmp_path / "env.lbm2"
    lib.touch()
    exe = _exe_dir_with(tmp_path)
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe),
                      lbm_setting=_env_lbm(str(lib)))
    assert res.source == "env"
    assert Path(res.path) == lib


def test_resolve_lbm_errors_when_env_path_is_missing(tmp_path):
    exe = _exe_dir_with(tmp_path, "one.lbm2")
    res = resolve_lbm({}, tmp_path / "param.txt", omics="lipidomics",
                      exe_path=str(exe), lbm_setting=_env_lbm(str(tmp_path / "gone.lbm2")))
    assert res.error_code == "LBM_NOT_FOUND"


def test_resolve_lbm_takes_the_value_from_the_config_file(tmp_path):
    lib = tmp_path / "conf.lbm2"
    lib.touch()
    exe = _exe_dir_with(tmp_path)
    setting = Setting(key="msdial.lbm", value=str(lib), source="config_file",
                      env_var="MSDIAL_LBM", config_file=str(tmp_path / "lipidmix.local.toml"))
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe), lbm_setting=setting)
    assert res.source == "config_file"
    assert Path(res.path) == lib


def test_resolve_lbm_names_the_config_key_when_its_file_is_missing(tmp_path):
    exe = _exe_dir_with(tmp_path, "one.lbm2")
    setting = Setting(key="msdial.lbm", value=str(tmp_path / "gone.lbm2"), source="config_file",
                      env_var="MSDIAL_LBM", config_file=str(tmp_path / "lipidmix.local.toml"))
    res = resolve_lbm({}, tmp_path / "param.txt", omics="lipidomics",
                      exe_path=str(exe), lbm_setting=setting)
    assert res.error_code == "LBM_NOT_FOUND"
    assert "lipidmix.local.toml の [msdial] lbm" in res.message


def test_resolve_lbm_names_the_env_var_when_its_file_is_missing(tmp_path):
    exe = _exe_dir_with(tmp_path, "one.lbm2")
    res = resolve_lbm({}, tmp_path / "param.txt", omics="lipidomics",
                      exe_path=str(exe), lbm_setting=_env_lbm(tmp_path / "gone.lbm2"))
    assert "環境変数 MSDIAL_LBM" in res.message


def test_resolve_lbm_uses_exe_directory_when_exactly_one(tmp_path):
    """GUI が `Assembly.GetExecutingAssembly().Location` の隣を見るのと同じ。"""
    exe = _exe_dir_with(tmp_path, "only.lbm2")
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.source == "exe_dir"
    assert Path(res.path).name == "only.lbm2"


def test_resolve_lbm_errors_when_exe_directory_has_none(tmp_path):
    exe = _exe_dir_with(tmp_path)
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.error_code == "LBM_NOT_FOUND"


def test_resolve_lbm_errors_when_exe_directory_is_ambiguous(tmp_path):
    """GUI も 1 件でなければ MessageBox で止める（DatasetParameterSettingModel）。"""
    exe = _exe_dir_with(tmp_path, "a.lbm2", "b.lbm2")
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.error_code == "LBM_AMBIGUOUS"
    assert len(res.candidates) == 2


def test_resolve_lbm_not_required_for_metabolomics(tmp_path):
    """宣言が無ければ、metabolomics では GUI 流の自動補完（exe フォルダ等）をしない。"""
    exe = _exe_dir_with(tmp_path, "only.lbm2")
    res = resolve_lbm({}, tmp_path / "param.txt", omics="metabolomics",
                      exe_path=str(exe), lbm_setting=_env_lbm(str(tmp_path / "x.lbm2")))
    assert res.error_code is None
    assert res.source == "not_required"
    assert res.path is None


# ---------- find_method_candidates の behavior 移行先 ----------
# find_method_candidates は極性不一致を「落とす」挙動を持っていたが、それは
# discover_method_candidates が意図的に変えた点そのもの（needs_polarity_conversion
# として提示する）。ここでは find_method_candidates 固有だった他の挙動
# （自動保存パラメータの発見・mtime 順・has_lbm・重複排除・バイナリのスキップ）を
# discover_method_candidates / scan_dir_for_method_files に対して migrate する。

def _param(dir_: Path, name: str, ion_mode: str, *, lbm: str = "", omics: str = "Lipidomics"):
    p = dir_ / name
    p.write_text(
        f"Ion mode: {ion_mode}\nTarget omics: {omics}\nLbm file path: {lbm}\n",
        encoding="ascii")
    return p


def test_discover_finds_auto_saved_param_files(tmp_path):
    """GUI は実行のたびに `<project>_param_<endtimestamp>.txt` を自動保存する。"""
    root = tmp_path / "POS"
    root.mkdir()
    _param(root, "Dataset_2026_05_15_param_202605151055.txt", "Negative")
    (root / "unrelated.txt").write_text("hello\n", encoding="ascii")
    found, _searched = discover_method_candidates(root, polarity=None)
    assert [Path(c.path).name for c in found] == ["Dataset_2026_05_15_param_202605151055.txt"]
    assert found[0].ion_mode == "negative"


def test_discover_sorts_newest_first_within_the_same_usable_group(tmp_path):
    """usable が揃う候補どうしは mtime の新しい順。"""
    root = tmp_path / "POS"
    root.mkdir()
    old = _param(root, "old_param_1.txt", "Positive")
    new = _param(root, "new_param_2.txt", "Positive")
    import os
    os.utime(old, (time.time() - 5000, time.time() - 5000))
    found, _searched = discover_method_candidates(root, polarity="positive")
    assert [Path(c.path).name for c in found] == [new.name, old.name]


def test_scan_dir_for_method_files_reports_missing_lbm(tmp_path):
    _param(tmp_path, "a_param_1.txt", "Positive", lbm="")
    _param(tmp_path, "b_param_2.txt", "Positive", lbm="C:\\x.lbm2")
    by_name = {Path(c.path).name: c for c in scan_dir_for_method_files(tmp_path, "same_dir")}
    assert by_name["a_param_1.txt"].has_lbm is False
    assert by_name["b_param_2.txt"].has_lbm is True


def test_discover_deduplicates_overlapping_dirs(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    _param(root, "a_param_1.txt", "Positive")
    found, _searched = discover_method_candidates(root, polarity="positive", search_dirs=[root])
    assert len(found) == 1


def test_scan_dir_for_method_files_skips_unreadable_binary(tmp_path):
    p = tmp_path / "bad_param_1.txt"
    p.write_bytes(b"\x00\x01\x02binary")
    assert scan_dir_for_method_files(tmp_path, "same_dir") == []


# ---------- write_effective_method_file ----------

def test_write_effective_method_file_replaces_existing_key(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("Ion mode: Positive\nLbm file path: \nSolvent type: CH3COONH4\n",
                   encoding="ascii")
    dest = tmp_path / "out" / "effective.txt"
    write_effective_method_file(src, dest, {"Lbm file path": "C:\\lib\\x.lbm2"})
    lines = dest.read_text(encoding="ascii").splitlines()
    assert lines == ["Ion mode: Positive",
                     "Lbm file path: C:\\lib\\x.lbm2",
                     "Solvent type: CH3COONH4"]


def test_write_effective_method_file_appends_absent_key(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("Ion mode: Positive\n", encoding="ascii")
    dest = tmp_path / "effective.txt"
    write_effective_method_file(src, dest, {"Lbm file path": "C:\\lib\\x.lbm2"})
    assert "Lbm file path: C:\\lib\\x.lbm2" in dest.read_text(encoding="ascii").splitlines()


def test_write_effective_method_file_is_ascii_with_lf(tmp_path):
    """ConfigParser は StreamReader(path, Encoding.ASCII) で 1 行ずつ読む。"""
    src = tmp_path / "src.txt"
    src.write_text("Ion mode: Positive\n", encoding="ascii")
    dest = tmp_path / "effective.txt"
    write_effective_method_file(src, dest, {"Lbm file path": "x.lbm2"})
    raw = dest.read_bytes()
    assert b"\r\n" not in raw
    raw.decode("ascii")  # 非 ASCII が混ざっていないこと


def test_write_effective_method_file_leaves_source_untouched(tmp_path):
    src = tmp_path / "src.txt"
    original = "Ion mode: Positive\nLbm file path: \n"
    src.write_text(original, encoding="ascii")
    write_effective_method_file(src, tmp_path / "effective.txt",
                               {"Lbm file path": "x.lbm2"})
    assert src.read_text(encoding="ascii") == original


# ---------- resolve_lbm: ビルドツリー探索 ----------
#
# この環境の Console 実行体は MsdialWorkbench のビルド生成物
# （tests/MSDIAL5/MsdialCoreTestApp/bin/Debug/net48/MSDIALCUI.exe）で、
# その exe フォルダに .lbm2 は無い。ライブラリは GUI アプリ側のビルド出力
# （src/MSDIAL5/MsdialGuiApp/bin/Debug/**/*.lbm2）にあるため、exe フォルダ
# 探索だけでは原理的に当たらない。インストール版（MSDIAL_LBM が指す）より
# ビルド生成物を優先する。

def _build_tree(tmp_path, *, msdial5_tfms=("", "net48", "net481"),
                msdial4_tfms=("net48",), exe_tfm="net48"):
    """MsdialWorkbench のビルド生成物を模したツリーを作り、Console exe のパスを返す。"""
    repo = tmp_path / "MsdialWorkbench"
    for tfm in msdial5_tfms:
        d = repo / "src/MSDIAL5/MsdialGuiApp/bin/Debug"
        if tfm:
            d = d / tfm
        d.mkdir(parents=True, exist_ok=True)
        (d / "Msp20260116160945_NCDK_conventional_converted_dev.lbm2").touch()
    for tfm in msdial4_tfms:
        d = repo / "src/MSDIAL4/MsDial/bin/Debug" / tfm
        d.mkdir(parents=True, exist_ok=True)
        (d / "Msp20221205132019_conventional_converted.lbm2").touch()
    exe_dir = repo / "tests/MSDIAL5/MsdialCoreTestApp/bin/Debug" / exe_tfm
    exe_dir.mkdir(parents=True, exist_ok=True)
    exe = exe_dir / "MSDIALCUI.exe"
    exe.touch()
    return repo, exe


def test_resolve_lbm_prefers_the_build_tree_over_the_installed_env_path(tmp_path):
    repo, exe = _build_tree(tmp_path)
    installed = tmp_path / "MSDIAL.v5.5-net48" / "Msp20251120132005_NCDK_dev.lbm2"
    installed.parent.mkdir(parents=True)
    installed.touch()
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe),
                      lbm_setting=_env_lbm(str(installed)))
    assert res.source == "build_tree"
    assert res.error_code is None
    assert Path(res.path).parent == repo / "src/MSDIAL5/MsdialGuiApp/bin/Debug/net48"


def test_resolve_lbm_build_tree_picks_the_tfm_matching_the_console_exe(tmp_path):
    """同名コピーが TFM ごとに並ぶので、exe 自身の TFM と揃える。"""
    repo, exe = _build_tree(tmp_path, msdial5_tfms=("", "net472", "net48", "net481"),
                            exe_tfm="net481")
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert Path(res.path).parent == repo / "src/MSDIAL5/MsdialGuiApp/bin/Debug/net481"


def test_resolve_lbm_build_tree_falls_back_to_the_plain_debug_copy(tmp_path):
    """exe の TFM に対応するコピーが無ければ素の Debug/ を使う。"""
    repo, exe = _build_tree(tmp_path, msdial5_tfms=("", "net472"), exe_tfm="net48")
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert Path(res.path).parent == repo / "src/MSDIAL5/MsdialGuiApp/bin/Debug"


def test_resolve_lbm_build_tree_is_not_ambiguous_for_tfm_copies(tmp_path):
    """同名の TFM 別コピーは『候補 4 件』ではない。LBM_AMBIGUOUS にしない。"""
    _repo, exe = _build_tree(tmp_path, msdial5_tfms=("", "net472", "net48", "net481"))
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.error_code is None
    assert res.candidates == ()


def test_resolve_lbm_build_tree_excludes_msdial4(tmp_path):
    """MSDIAL4 の conventional ライブラリは NCDK 無しの別世代。Console(MSDIAL5) に混ぜない。"""
    _repo, exe = _build_tree(tmp_path, msdial5_tfms=(), msdial4_tfms=("net48", "net481"))
    installed = tmp_path / "installed.lbm2"
    installed.touch()
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe),
                      lbm_setting=_env_lbm(str(installed)))
    assert res.source == "env"
    assert Path(res.path) == installed


def test_resolve_lbm_build_tree_is_ambiguous_for_distinct_libraries(tmp_path):
    """TFM 複製ではなく、ファイル名が異なる .lbm2 が同じフォルダに複数ある場合は
    found[0]（アルファベット順）を黙って選ばず LBM_AMBIGUOUS にする。"""
    repo = tmp_path / "MsdialWorkbench"
    debug_net48 = repo / "src/MSDIAL5/MsdialGuiApp/bin/Debug/net48"
    debug_net48.mkdir(parents=True)
    (debug_net48 / "Msp_a.lbm2").touch()
    (debug_net48 / "Msp_b.lbm2").touch()
    exe_dir = repo / "tests/MSDIAL5/MsdialCoreTestApp/bin/Debug/net48"
    exe_dir.mkdir(parents=True)
    exe = exe_dir / "MSDIALCUI.exe"
    exe.touch()
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.error_code == "LBM_AMBIGUOUS"
    assert res.path is None
    assert len(res.candidates) == 2


def test_resolve_lbm_method_file_declaration_still_beats_the_build_tree(tmp_path):
    """メソッドファイルの明示宣言はビルドツリーより強い（MS-DIAL の規則を維持）。"""
    _repo, exe = _build_tree(tmp_path)
    declared = tmp_path / "declared.lbm2"
    declared.touch()
    res = resolve_lbm({"lbm file path": str(declared)}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe))
    assert res.source == "method_file"
    assert Path(res.path) == declared


# ---------- key_params ----------

def test_extract_key_params_keeps_only_the_decision_relevant_keys():
    keys = {"ion mode": "Positive", "target omics": "Lipidomics",
            "minimum peak height": "1000", "smoothing level": "3",
            "ms1 tolerance for centroid": "0.01"}
    out = extract_key_params(keys)
    assert out == {"Ion mode": "Positive", "Target omics": "Lipidomics",
                   "Minimum peak height": "1000",
                   "MS1 tolerance for centroid": "0.01"}


def test_extract_key_params_uses_the_spelling_from_the_real_param_file():
    """実ファイルの綴りは `MS1 mass range begin`。`Mass range begin` ではない。"""
    keys = {"ms1 mass range begin": "0", "ms1 mass range end": "2000",
            "retention time tolerance for alignment": "0.1",
            "ms1 tolerance for alignment": "0.015"}
    assert set(extract_key_params(keys)) == {
        "MS1 mass range begin", "MS1 mass range end",
        "Retention time tolerance for alignment", "MS1 tolerance for alignment"}


def test_extract_key_params_summarises_the_adduct_list():
    """POS の実値は 37 種・約 700 文字。候補 10 件で 7 KB になるので要約する。"""
    adducts = ",".join(["[M+H]+", "[M+NH4]+", "[M+Na]+"] + [f"[M+X{i}]+" for i in range(34)])
    out = extract_key_params({"searched adduct ions": adducts})
    assert out["Searched adduct ions"] == "37 種（先頭: [M+H]+, [M+NH4]+, [M+Na]+）"


def test_extract_key_params_omits_absent_keys():
    assert extract_key_params({"ion mode": "Negative"}) == {"Ion mode": "Negative"}


def test_method_candidate_defaults_keep_existing_callers_working():
    c = MethodCandidate(path="p", ion_mode="positive", omics="lipidomics",
                        has_lbm=False, mtime=0.0)
    assert (c.origin, c.usable, c.key_params) == ("same_dir", "direct", None)


# ---------- scan_dir_for_method_files ----------

def test_scan_dir_for_method_files_labels_the_origin(tmp_path):
    p = tmp_path / "Dataset_2026_05_15_10_12_46_param_202605151055.txt"
    p.write_text("Ion mode: Negative\nTarget omics: Lipidomics\n", encoding="ascii")
    found = scan_dir_for_method_files(tmp_path, "sibling")
    assert [c.origin for c in found] == ["sibling"]
    assert found[0].ion_mode == "negative"


def test_scan_dir_for_method_files_does_not_filter_by_polarity(tmp_path):
    """絞り込みは呼び出し側の仕事。ここで落とすと別極性の候補を提示できない。"""
    for name, ion in (("a_param_1.txt", "Positive"), ("b_param_2.txt", "Negative")):
        (tmp_path / name).write_text(f"Ion mode: {ion}\n", encoding="ascii")
    assert len(scan_dir_for_method_files(tmp_path, "same_dir")) == 2


def test_scan_dir_for_method_files_skips_unreadable_directory(tmp_path):
    assert scan_dir_for_method_files(tmp_path / "nope", "same_dir") == []


# ---------- 探索先の列挙 ----------

def test_method_search_dirs_lists_self_then_siblings(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    (tmp_path / "NEG").mkdir()
    (tmp_path / "OTHER").mkdir()
    pairs = method_search_dirs(root)
    assert pairs[0] == (root, "same_dir")
    assert sorted(str(d.name) for d, o in pairs if o == "sibling") == ["NEG", "OTHER"]


def test_method_search_dirs_excludes_itself_from_the_siblings(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    (tmp_path / "NEG").mkdir()
    assert [d for d, o in method_search_dirs(root) if o == "sibling"] == [tmp_path / "NEG"]


def test_method_search_dirs_appends_explicit_dirs_as_given(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    extra = tmp_path / "elsewhere"
    extra.mkdir()
    pairs = method_search_dirs(root, [str(extra)])
    assert (extra, "given") in pairs


def test_method_search_dirs_does_not_recurse(tmp_path):
    """兄弟までで止める。深く掘ると無関係なパラメータが候補を濁らせる。"""
    root = tmp_path / "POS"
    root.mkdir()
    deep = tmp_path / "NEG" / "inner"
    deep.mkdir(parents=True)
    assert deep not in [d for d, _o in method_search_dirs(root)]


def test_past_run_method_files_reads_the_recorded_method(tmp_path):
    """analysis-job.json の software.method_file が実際に使ったメソッドを指す。"""
    import json
    used = tmp_path / "param_POS_generated.txt"
    used.write_text("Ion mode: Positive\n", encoding="ascii")
    run = tmp_path / "runs" / "job_1"
    run.mkdir(parents=True)
    (run / "analysis-job.json").write_text(
        json.dumps({"software": {"method_file": str(used)}}), encoding="utf-8")
    assert past_run_method_files(tmp_path) == [used]


def test_past_run_method_files_skips_records_pointing_at_a_deleted_file(tmp_path):
    import json
    run = tmp_path / "runs" / "job_1"
    run.mkdir(parents=True)
    (run / "analysis-job.json").write_text(
        json.dumps({"software": {"method_file": str(tmp_path / "gone.txt")}}), encoding="utf-8")
    assert past_run_method_files(tmp_path) == []


def test_past_run_method_files_ignores_a_broken_job_file(tmp_path):
    run = tmp_path / "runs" / "job_1"
    run.mkdir(parents=True)
    (run / "analysis-job.json").write_text("{not json", encoding="utf-8")
    assert past_run_method_files(tmp_path) == []


def test_past_run_method_files_ignores_a_top_level_json_array(tmp_path):
    """valid JSON だが object ではない（[]・"x"・3 等）は例外にせず読み飛ばす。"""
    run = tmp_path / "runs" / "job_1"
    run.mkdir(parents=True)
    (run / "analysis-job.json").write_text("[]", encoding="utf-8")
    assert past_run_method_files(tmp_path) == []


def test_past_run_method_files_ignores_a_non_dict_software_value(tmp_path):
    """`software` が object でない（配列等）ときも .get() で落ちずに読み飛ばす。

    空配列は既存の `or {}` が偶然拾うため、非空の配列で確実に .get() へ届かせる。
    """
    import json
    run = tmp_path / "runs" / "job_1"
    run.mkdir(parents=True)
    (run / "analysis-job.json").write_text(
        json.dumps({"software": ["not", "a", "dict"]}), encoding="utf-8")
    assert past_run_method_files(tmp_path) == []


# ---------- discover_method_candidates ----------

def _write_param(directory, name, ion="Positive", omics="Lipidomics"):
    directory.mkdir(parents=True, exist_ok=True)
    p = directory / name
    p.write_text(f"Ion mode: {ion}\nTarget omics: {omics}\n"
                 "Minimum peak height: 1000\nLbm file path: \n", encoding="ascii")
    return p


def test_discover_marks_the_other_polarity_as_needing_conversion(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    _write_param(tmp_path / "NEG", "d_param_1.txt", ion="Negative")
    found, _searched = discover_method_candidates(root, polarity="positive")
    assert [(c.origin, c.usable) for c in found] == [("sibling", "needs_polarity_conversion")]


def test_discover_puts_directly_usable_candidates_first(tmp_path):
    root = tmp_path / "POS"
    _write_param(root, "own_param_2.txt", ion="Positive")
    _write_param(tmp_path / "NEG", "neg_param_1.txt", ion="Negative")
    found, _searched = discover_method_candidates(root, polarity="positive")
    assert [c.usable for c in found] == ["direct", "needs_polarity_conversion"]


def test_discover_attaches_key_params_for_a_small_candidate_set(tmp_path):
    root = tmp_path / "POS"
    _write_param(root, "own_param_1.txt")
    found, _searched = discover_method_candidates(root, polarity="positive")
    assert found[0].key_params["Minimum peak height"] == "1000"


def test_discover_omits_key_params_beyond_the_cap(tmp_path):
    root = tmp_path / "POS"
    for i in range(KEY_PARAMS_MAX_CANDIDATES + 1):
        _write_param(root, f"p{i}_param_{i}.txt")
    found, _searched = discover_method_candidates(root, polarity="positive")
    assert len(found) == KEY_PARAMS_MAX_CANDIDATES + 1
    assert all(c.key_params is None for c in found)


def test_discover_prefers_past_run_origin_for_a_duplicate_path(tmp_path):
    """過去 run のメソッドがユーザーのフォルダにある元ファイルを指すことがある。"""
    import json
    root = tmp_path / "POS"
    used = _write_param(root, "own_param_1.txt")
    run = root / "runs" / "job_1"
    run.mkdir(parents=True)
    (run / "analysis-job.json").write_text(
        json.dumps({"software": {"method_file": str(used)}}), encoding="utf-8")
    found, _searched = discover_method_candidates(root, polarity="positive")
    assert len(found) == 1
    assert found[0].origin == "past_run"


def test_discover_filters_by_omics(tmp_path):
    root = tmp_path / "POS"
    _write_param(root, "lip_param_1.txt", omics="Lipidomics")
    _write_param(root, "met_param_2.txt", omics="Metabolomics")
    found, _searched = discover_method_candidates(root, polarity="positive", omics="lipidomics")
    assert [Path(c.path).name for c in found] == ["lip_param_1.txt"]


def test_discover_reports_where_it_looked(tmp_path):
    root = tmp_path / "POS"
    root.mkdir()
    (tmp_path / "NEG").mkdir()
    _found, searched = discover_method_candidates(root, polarity="positive")
    assert str(root) in searched
    assert str(tmp_path / "NEG") in searched


# ---------- 相対パス宣言（LC-MS Console は cwd 基準で読む） ----------
#
# MS-DIAL Console の LC-MS 経路（`ConfigParser.ReadForLcmsParameter`）は
# `Lbm file path` / `Msp file path` / `Text DB file path` 等を**解決しない**。
# 相対値はプロセスの cwd 基準になり、見つからなければ `CommonProcess.ParseLibraries`
# がそのライブラリを黙って飛ばす（同定 0 件で完走）。メソッド基準で解くのは
# GC-MS 経路の `ResolvePathFromMethodFile` だけ（MsdialWorkbench afd5f95）。


def test_resolve_lbm_returns_an_absolute_path_for_a_relative_declaration(tmp_path, monkeypatch):
    """メソッドファイル自体が相対で渡されても、Console へ渡す値は絶対パスにする。"""
    (tmp_path / "proj").mkdir()
    lib = tmp_path / "proj" / "lib.lbm2"
    lib.touch()
    monkeypatch.chdir(tmp_path)
    res = resolve_lbm({"lbm file path": "lib.lbm2"}, Path("proj") / "param.txt",
                      omics="lipidomics", exe_path=None)
    assert res.error_code is None
    assert Path(res.path).is_absolute()
    assert Path(res.path) == lib


def test_resolve_lbm_honours_a_declared_lbm_for_metabolomics(tmp_path):
    """Console は Target omics に関わらず `Lbm file path` を読み、脂質同定に使う
    （`LcmsProcess` が LBM annotator を `TargetOmics.Lipidomics` 固定で組む）。"""
    lib = tmp_path / "lib.lbm2"
    lib.touch()
    res = resolve_lbm({"lbm file path": "lib.lbm2"}, tmp_path / "param.txt",
                      omics="metabolomics", exe_path=None)
    assert res.error_code is None
    assert res.source == "method_file"
    assert Path(res.path) == lib


def test_resolve_lbm_rejects_a_missing_declared_lbm_for_metabolomics(tmp_path):
    res = resolve_lbm({"lbm file path": "gone.lbm2"}, tmp_path / "param.txt",
                      omics="metabolomics", exe_path=None)
    assert res.error_code == "LBM_NOT_FOUND"


def test_resolve_lbm_honours_the_argument_for_metabolomics(tmp_path):
    lib = tmp_path / "chosen.lbm2"
    lib.touch()
    res = resolve_lbm({}, tmp_path / "param.txt", omics="metabolomics",
                      exe_path=None, override=str(lib))
    assert res.source == "argument"
    assert Path(res.path) == lib


def test_relative_path_overrides_absolutizes_library_paths_against_the_method(tmp_path):
    method = tmp_path / "proj" / "param.txt"
    keys = {
        "lbm file path": "lib.lbm2",
        "msp file path": r"..\libs\a.msp",
        "text db file path": "t.txt",
        "isotope text db file path": "iso.txt",
        "compounds library file path for target detection": "target.txt",
        "compounds library file path for rt correction": "rt.txt",
        "rt correction peak selection file path": "sel.tsv",
    }
    out = relative_path_overrides(keys, method)
    base = method.parent
    assert {k.lower(): v for k, v in out.items()} == {
        "lbm file path": str((base / "lib.lbm2").resolve()),
        "msp file path": str((tmp_path / "libs" / "a.msp").resolve()),
        "text db file path": str((base / "t.txt").resolve()),
        "isotope text db file path": str((base / "iso.txt").resolve()),
        "compounds library file path for target detection": str((base / "target.txt").resolve()),
        "compounds library file path for rt correction": str((base / "rt.txt").resolve()),
        "rt correction peak selection file path": str((base / "sel.tsv").resolve()),
    }


def test_relative_path_overrides_absolutizes_annotator_settings_paths(tmp_path):
    """設定表のパスは Console 自身がメソッドファイル基準で解く。実効コピーを別フォルダに
    書くとその基準がずれるので、原本基準の絶対パスへ固定する。"""
    method = tmp_path / "param.txt"
    keys = {"msp annotator settings file path": "msp.tsv",
            "text db annotator settings file path": "text.tsv"}
    out = {k.lower(): v for k, v in relative_path_overrides(keys, method).items()}
    assert out == {"msp annotator settings file path": str((tmp_path / "msp.tsv").resolve()),
                   "text db annotator settings file path": str((tmp_path / "text.tsv").resolve())}


def test_relative_path_overrides_leaves_absolute_empty_and_unknown_keys_alone(tmp_path):
    keys = {"msp file path": str(tmp_path / "abs.msp"),
            "lbm file path": "",
            "ion mode": "Positive",
            "some other file": "x.txt"}
    assert relative_path_overrides(keys, tmp_path / "param.txt") == {}


def test_write_effective_method_file_replaces_every_duplicate_line(tmp_path):
    """ConfigParser は全行を順に読み、後の行が勝つ。最初の 1 行だけ差し替えると
    後ろの重複行が上書きを打ち消す。"""
    src = tmp_path / "src.txt"
    src.write_text("Msp file path: a.msp\nIon mode: Positive\nmsp file path: b.msp\n",
                   encoding="ascii")
    dest = tmp_path / "effective.txt"
    write_effective_method_file(src, dest, {"Msp file path": "C:\\lib\\x.msp"})
    lines = dest.read_text(encoding="ascii").splitlines()
    assert lines == ["Msp file path: C:\\lib\\x.msp", "Ion mode: Positive",
                     "Msp file path: C:\\lib\\x.msp"]


# ---------- library_mode="msp_only": MSP の解決と実効メソッドの上書き ----------

import os

from metabolomix.console import method_file as _mf
from metabolomix.core import path_resolvers as _path_resolvers
from metabolomix.core import user_config as _user_config


def _msp_setting(path, *, source="config_file", key="library.msp_negative"):
    return _user_config.Setting(
        key=key, value=str(path), source=source,
        env_var=_user_config.SETTINGS[key],
        config_file=None if source == "env" else "C:/cfg/lipidmix.local.toml")


def test_msp_setting_keys_are_the_single_source_for_library_lookup():
    assert _user_config.MSP_SETTING_KEYS == {
        "positive": "library.msp_positive", "negative": "library.msp_negative"}
    assert _path_resolvers.LIBRARY_SETTING_KEYS is _user_config.MSP_SETTING_KEYS


def test_resolve_msp_prefers_the_argument(tmp_path):
    arg = tmp_path / "given.msp"
    arg.write_text("NAME: a\n", encoding="ascii")
    other = tmp_path / "configured.msp"
    other.write_text("NAME: b\n", encoding="ascii")
    res = _mf.resolve_msp(str(arg), "negative", _msp_setting(other))
    assert res == _mf.MspResolution(path=os.path.abspath(arg), source="argument")


def test_resolve_msp_uses_the_setting_and_reports_where_it_came_from(tmp_path):
    lib = tmp_path / "lab_neg.msp"
    lib.write_text("NAME: a\n", encoding="ascii")
    from_file = _mf.resolve_msp(None, "negative", _msp_setting(lib))
    assert from_file.path == os.path.abspath(lib) and from_file.source == "config_file"
    from_env = _mf.resolve_msp(None, "negative", _msp_setting(lib, source="env"))
    assert from_env.source == "env"


def test_resolve_msp_not_configured_names_the_setting_key_and_env_var():
    res = _mf.resolve_msp(None, "positive", None)
    assert res.error_code == "MSP_NOT_CONFIGURED"
    assert res.path is None and res.source == "not_configured"
    assert "msp_positive" in res.message and "MSDIAL_MSP_POS" in res.message


def test_resolve_msp_not_found_shows_only_the_file_name(tmp_path):
    hidden_dir = tmp_path / "secret_lab_share"
    missing = hidden_dir / "lab_neg.msp"
    res = _mf.resolve_msp(None, "negative", _msp_setting(missing))
    assert res.error_code == "MSP_NOT_FOUND"
    assert "lab_neg.msp" in res.message
    assert "secret_lab_share" not in res.message
    arg = _mf.resolve_msp(str(missing), "negative", None)
    assert arg.error_code == "MSP_NOT_FOUND" and arg.source == "argument"
    assert "secret_lab_share" not in arg.message


def test_msp_only_overrides_blank_lbm_set_msp_and_clear_other_identification_keys():
    method_keys = {
        "ion mode": "Negative",
        "lbm file path": "C:/libs/lipids.lbm2",
        "msp file path": "C:/public/other.msp",
        "text db file path": "db.txt",                    # 小文字の綴りでも拾う
        "msp annotator settings file path": "msp.tsv",
        "isotope text db file path": "iso.txt",           # 同定用ではないので残す
        "compounds library file path for rt correction": "rt.txt",
    }
    overrides, removed = _mf.msp_only_overrides(method_keys, "C:/lab/lab_neg.msp")
    assert overrides[_mf.LBM_KEY] == ""
    assert overrides[_mf.MSP_KEY] == "C:/lab/lab_neg.msp"
    assert overrides[_mf.TEXT_DB_KEY] == ""
    assert overrides["MSP annotator settings file path"] == ""
    assert "Isotope text DB file path" not in overrides
    assert _mf.RT_REFERENCE_KEY not in overrides
    assert removed == [_mf.LBM_KEY, _mf.TEXT_DB_KEY, "MSP annotator settings file path"]


def test_msp_only_overrides_does_not_add_keys_the_method_never_declared():
    overrides, removed = _mf.msp_only_overrides({"ion mode": "Positive", "lbm file path": ""},
                                                "C:/lab/lab_pos.msp")
    assert overrides == {_mf.LBM_KEY: "", _mf.MSP_KEY: "C:/lab/lab_pos.msp"}
    assert removed == []


def test_blank_override_replaces_every_duplicate_line(tmp_path):
    src = tmp_path / "params.txt"
    src.write_text(
        "Ion mode: Negative\n"
        "Lbm file path: C:/libs/a.lbm2\n"
        "MSP FILE PATH: C:/public/one.msp\n"
        "Lbm file path: C:/libs/b.lbm2\n"
        "msp file path: C:/public/two.msp\n",
        encoding="ascii")
    overrides, _ = _mf.msp_only_overrides(_mf.read_method_keys(src), "C:/lab/lab_neg.msp")
    dest = _mf.write_effective_method_file(src, tmp_path / "out" / "effective.txt", overrides)
    lines = dest.read_text(encoding="ascii").splitlines()
    assert lines.count("Lbm file path: ") == 2
    assert lines.count("Msp file path: C:/lab/lab_neg.msp") == 2
    assert not any("one.msp" in line or "two.msp" in line or ".lbm2" in line for line in lines)
