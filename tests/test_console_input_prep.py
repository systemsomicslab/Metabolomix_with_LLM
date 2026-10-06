"""単一フォーマットの入力フォルダを作る（MIXED_RAW_FORMATS の解き方）。

MS-DIAL の `SupportMsRawDataExtension` は `wiff` と `wiff2` を**別フォーマット**として
数えるため、SCIEX が 1 測定につき両方出すのが普通の環境では必ず混在する。混在すると
`AnalysisFilesParser.ReadInput` が `Console.ReadLine()` で Y/N を聞き、stdin を塞いだ
実行では異常終了する（続行できても 60 サンプルが 120 解析ファイルになる）。

エラー封筒は「片方だけ残したフォルダを作れ」と正しく言うが、MCP クライアントには
フォルダを作る手段が無い。ここがその手段。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from metabolomix.console.input_prep import prepare_single_format_input


def _sciex_pair(root: Path, stem: str):
    """SCIEX が 1 測定につき出す一式。"""
    (root / f"{stem}.wiff").write_text("primary", encoding="ascii")
    (root / f"{stem}.wiff.scan").write_text("scan", encoding="ascii")
    (root / f"{stem}.wiff2").write_text("newer", encoding="ascii")
    (root / f"{stem}.timeseries.data").write_text("meta", encoding="ascii")


def test_keeps_only_the_requested_raw_format(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "sample_a")
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    names = {p.name for p in Path(result.out_dir).iterdir()}
    assert "sample_a.wiff" in names
    assert "sample_a.wiff2" not in names


def test_carries_companion_files_of_the_kept_format(tmp_path):
    """.wiff は .wiff.scan が無いと読めない。随伴ファイルを置いていってはいけない。"""
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "sample_a")
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    names = {p.name for p in Path(result.out_dir).iterdir()}
    assert "sample_a.wiff.scan" in names
    assert "sample_a.timeseries.data" in names


def test_result_counts_primary_and_companions(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    _sciex_pair(src, "b")
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    assert result.primary == 2
    assert result.companions == 4


def test_uses_hardlinks_on_the_same_volume(tmp_path):
    """4GB の実体コピーを避ける。生データは読むだけなのでリンクで足りる。"""
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    assert result.mode == "hardlink"
    original = src / "a.wiff"
    linked = Path(result.out_dir) / "a.wiff"
    assert original.stat().st_ino == linked.stat().st_ino


def test_leaves_the_source_folder_untouched(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    before = sorted(p.name for p in src.iterdir())
    prepare_single_format_input(src, "wiff", tmp_path / "out")
    assert sorted(p.name for p in src.iterdir()) == before


def test_output_folder_has_exactly_one_raw_format(tmp_path):
    """作った先が console_plan の raw_input_summary を 1 種類で通ること。"""
    from metabolomix.console.job_manager import raw_input_summary
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    _sciex_pair(src, "b")
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    assert raw_input_summary(Path(result.out_dir)) == {"wiff": 2}


def test_rejects_extension_absent_from_source(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    with pytest.raises(ValueError, match="raw"):
        prepare_single_format_input(src, "raw", tmp_path / "out")


def test_rejects_unknown_extension(tmp_path):
    """MS-DIAL が計測ファイルとして数えない拡張子を主として選んでも意味がない。"""
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    with pytest.raises(ValueError):
        prepare_single_format_input(src, "scan", tmp_path / "out")


def test_is_idempotent(tmp_path):
    """再実行で既存リンクを壊さない（掃除と作り直しを往復するため）。"""
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    out = tmp_path / "out"
    prepare_single_format_input(src, "wiff", out)
    result = prepare_single_format_input(src, "wiff", out)
    assert result.skipped_existing == 3
    assert (out / "a.wiff").read_text(encoding="ascii") == "primary"


def test_refuses_to_write_into_the_source_folder(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    with pytest.raises(ValueError):
        prepare_single_format_input(src, "wiff", src)


def test_same_stem_analysis_artifacts_are_not_pulled_in_as_companions(tmp_path):
    """PAI2/タグ等、同じstemの既存解析出力を随伴ファイルとして誤って拾わない。

    旧実装は「primaryの先頭トークンで始まり、計測拡張子でないもの」を随伴と
    見なしていたため、同じstemのPAI2やタグXMLまで拾ってしまっていた
    （metabolomix/console/input_prep.py の COMPANION_RULES 導入前の欠陥）。
    """
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "sample_a")
    (src / "sample_a.pai2").write_bytes(b"\x00pai2-binary")
    (src / "sample_a_tags.xml").write_text("<tags/>", encoding="ascii")
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    names = {p.name for p in Path(result.out_dir).iterdir()}
    assert "sample_a.pai2" not in names
    assert "sample_a_tags.xml" not in names
    assert result.companions == 2  # .wiff.scan と .timeseries.data だけ


def test_directories_in_source_are_ignored(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _sciex_pair(src, "a")
    (src / "a.d").mkdir()  # Agilent の .d はフォルダ
    result = prepare_single_format_input(src, "wiff", tmp_path / "out")
    assert not (Path(result.out_dir) / "a.d").exists()


# ---------- フォルダ形式の計測データ（Agilent/Bruker の .d、Waters の .raw） ----------
#
# MS-DIAL の `AnalysisFilesParser.ReadFolderContents` は、拡張子が `.raw` / `.d` の
# **ディレクトリ**を 1 検体の計測データとして受ける（`isVendorDirectory`）。
# `.wiff` のようなファイル 1 つが 1 検体という前提だけで書くと、この形式の
# データフォルダでは「.d がありません」と嘘をつく。


def _vendor_dir(root: Path, name: str) -> Path:
    """Agilent の .d を模した、入れ子を持つ計測フォルダ。"""
    d = root / name
    (d / "AcqData").mkdir(parents=True)
    (d / "AcqData" / "MSScan.bin").write_text("scan", encoding="ascii")
    (d / "AcqData" / "Devices.xml").write_text("<Devices/>", encoding="ascii")
    (d / "Contents.xml").write_text("<Contents/>", encoding="ascii")
    return d


def test_keeps_directory_format_raw(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    for stem in ("a", "b", "c"):
        _vendor_dir(src, f"{stem}.d")
    (src / "stray.mzml").write_text("mzml", encoding="ascii")

    result = prepare_single_format_input(src, "d", tmp_path / "out")

    out = Path(result.out_dir)
    assert result.primary == 3
    assert sorted(p.name for p in out.iterdir()) == ["a.d", "b.d", "c.d"]
    assert (out / "a.d" / "AcqData" / "MSScan.bin").read_text(encoding="ascii") == "scan"
    assert not (out / "stray.mzml").exists()


def test_directory_format_raw_inner_files_are_hardlinked(tmp_path):
    """.d は数GB になりうる。中身まで実体コピーするとこのツールの意味が消える。

    ディレクトリ**そのもの**はリンクにしない（実体のフォルダを作る）。
    リンクするのは中のファイルだけ。
    """
    src = tmp_path / "raw"
    src.mkdir()
    _vendor_dir(src, "a.d")

    result = prepare_single_format_input(src, "d", tmp_path / "out")

    out = Path(result.out_dir)
    assert (out / "a.d").is_dir() and not (out / "a.d").is_symlink()
    inner = out / "a.d" / "AcqData" / "MSScan.bin"
    assert inner.stat().st_nlink == 2
    assert inner.stat().st_ino == (src / "a.d" / "AcqData" / "MSScan.bin").stat().st_ino
    assert result.linked == 3 and result.copied == 0
    assert result.mode == "hardlink"


def test_directory_format_raw_leaves_the_source_untouched(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _vendor_dir(src, "a.d")
    before = sorted(str(p.relative_to(src)) for p in src.rglob("*"))

    prepare_single_format_input(src, "d", tmp_path / "out")

    assert sorted(str(p.relative_to(src)) for p in src.rglob("*")) == before


def test_directory_format_raw_is_idempotent(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    _vendor_dir(src, "a.d")

    first = prepare_single_format_input(src, "d", tmp_path / "out")
    second = prepare_single_format_input(src, "d", tmp_path / "out")

    assert first.linked == 3
    assert second.linked == 0 and second.copied == 0
    assert second.skipped_existing == 3
    assert (Path(second.out_dir) / "a.d" / "Contents.xml").exists()


def test_waters_raw_directory_is_kept_as_a_measurement(tmp_path):
    """Waters の `.raw` はフォルダ。Thermo の `.raw` はファイル。拡張子は同じ。"""
    src = tmp_path / "raw"
    src.mkdir()
    d = src / "sample.raw"
    d.mkdir()
    (d / "_FUNC001.DAT").write_text("func", encoding="ascii")

    result = prepare_single_format_input(src, "raw", tmp_path / "out")

    assert result.primary == 1
    assert (Path(result.out_dir) / "sample.raw" / "_FUNC001.DAT").exists()


def test_non_vendor_directory_is_not_a_measurement(tmp_path):
    """`.d` / `.raw` 以外の拡張子はフォルダだと計測データにならない（上流の規則）。

    `DataAccess.IsDataFormatSupported` は `File.Exists` を要求するので、
    `backup.mzml` という名前のフォルダは MS-DIAL の入力にならない。
    """
    src = tmp_path / "raw"
    src.mkdir()
    (src / "backup.mzml").mkdir()
    (src / "real.mzml").write_text("mzml", encoding="ascii")

    result = prepare_single_format_input(src, "mzml", tmp_path / "out")

    assert result.primary == 1
    assert (Path(result.out_dir) / "real.mzml").exists()
    assert not (Path(result.out_dir) / "backup.mzml").exists()
