"""単一フォーマットの入力フォルダを作る（`MIXED_RAW_FORMATS` の解き方）。

MS-DIAL の `SupportMsRawDataExtension` は `wiff` と `wiff2` を**別フォーマット**として
数える。SCIEX は 1 測定につき両方を出すので、生データフォルダはそのままでは必ず
混在し、`AnalysisFilesParser.ReadInput` が `Console.ReadLine()` で Y/N を聞く。
stdin を塞いだ実行では NullReferenceException で終了コード 1 になり、Y と答えて
進めた場合は 60 サンプルが 120 解析ファイルとして扱われる
（docs/HISTRY.md 2026-09-03(5) §4）。

`console_plan` のエラー封筒は「解析に使うほうだけを残したフォルダを作って指定して
ください」と正しく言うが、**MCP クライアントにはフォルダを作る手段が無い**。
ここがその手段。生データは読むだけなのでハードリンクで足りる（実測 4.1GB の
フォルダを実体コピーせずに済む）。

Agilent / Bruker の `.d` と Waters の `.raw` は**フォルダそのものが 1 検体**で、
MS-DIAL もそう読む（`AnalysisFilesParser.ReadFolderContents` の
`isVendorDirectory`）。この形式はフォルダ構造を実体として作り直し、中の
ファイルだけをリンクする —— ディレクトリ自体をリンクやジャンクションに
すると、隔離したように見えて元フォルダと同じものを指し続ける。
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

from metabolomix.console.job_manager import _RAW_EXTENSIONS, is_raw_input


@dataclass(frozen=True)
class PrepareResult:
    """作った入力フォルダの内訳。"""

    out_dir: str
    primary: int
    companions: int
    linked: int
    copied: int
    skipped_existing: int
    mode: str  # hardlink | copy | mixed | none


#: 随伴ファイルの明示規則（形式ごと）。要素は (基準, 追加suffix)。
#: 基準 "name" は primary のファイル名全体に、"stem" は主拡張子を除いた部分に
#: suffix を足した名前を候補にし、その名前がフォルダに実在するときだけ随伴として扱う。
#:
#: 以前は「primary の先頭トークンで始まり、計測拡張子でない」という前方一致で
#: 判定していたが、これは同じ stem を持つ PAI2/DCL/ARF/タグファイル
#: （例: `sample_a.pai2`）まで随伴として誤って拾ってしまう欠陥だった。
#: 明示登録した完全一致の名前だけを随伴とすることでこれを防ぐ。
COMPANION_RULES: dict[str, tuple[tuple[str, str], ...]] = {
    "wiff": (("name", ".scan"), ("stem", ".timeseries.data")),
}


def _companions_of(entries: list[Path], primary: Path, ext: str) -> list[Path]:
    """primary と同じ測定に属する随伴ファイルを、`COMPANION_RULES` の完全一致でだけ返す。

    `.wiff` は `.wiff.scan` が無いと読めず、`.timeseries.data` も同じ測定の一部。
    一方 `.wiff2` は**別フォーマット**なので連れて行ってはいけない —— それを
    置いていくことがこの関数の目的そのものだから（`COMPANION_RULES` に無い名前は
    そもそも候補に上がらない）。
    """
    by_name = {p.name: p for p in entries}
    out = []
    for base, suffix in COMPANION_RULES.get(ext, ()):
        name = (primary.name if base == "name" else primary.stem) + suffix
        found = by_name.get(name)
        if found is not None and found != primary:
            out.append(found)
    return out


def _replicate_dir(source: Path, target: Path) -> tuple[int, int, int]:
    """フォルダ形式の raw を複製し、`(リンク数, コピー数, 既存で飛ばした数)` を返す。

    ディレクトリは実体として作る。リンクするのは中のファイルだけ —— `.d` は
    数 GB になりうるのでファイルまで実体コピーするとこのツールの意味が消えるが、
    ディレクトリそのものをリンクすると元フォルダと同じ実体を指したまま
    「隔離した」と名乗ることになる。

    数え方はファイル形式の raw と揃える（1 検体ではなく 1 ファイルで数える）。
    `copytree` に任せず自分で歩くのは、リンクとコピーの内訳を数えるためと、
    **既にあるフォルダの中まで見て欠けているファイルだけを補う**ため。
    フォルダの存在だけで飛ばすと、前回が途中で止まって不完全な `.d` が
    そのまま完成扱いになる（`copytree` は `dirs_exist_ok` でも欠落の検出を
    しない）。
    """
    linked = copied = skipped = 0
    for dirpath, _dirnames, filenames in os.walk(source):
        rel = Path(dirpath).relative_to(source)
        (target / rel).mkdir(parents=True, exist_ok=True)
        for name in filenames:
            src_file = Path(dirpath) / name
            dst_file = target / rel / name
            if dst_file.exists():
                skipped += 1
                continue
            try:
                os.link(src_file, dst_file)
                linked += 1
            except OSError:
                shutil.copy2(src_file, dst_file)
                copied += 1
    return linked, copied, skipped


def prepare_single_format_input(
    dataset_root, keep_extension: str, out_dir,
) -> PrepareResult:
    """`keep_extension` の計測ファイルと随伴ファイルだけの入力フォルダを作る。

    同一ボリュームならハードリンク、またげないときはコピーへフォールバックする。
    元フォルダは読むだけで一切変更しない。既に同名がある出力先は上書きしない
    （掃除と作り直しを往復するので冪等であることが要る）。
    """
    src = Path(dataset_root).expanduser()
    dest = Path(out_dir).expanduser()
    ext = keep_extension.lower().lstrip(".")

    if ext not in _RAW_EXTENSIONS:
        raise ValueError(
            f"MS-DIAL が計測ファイルとして数えない拡張子です: {keep_extension!r}  "
            f"対象: {', '.join(sorted(_RAW_EXTENSIONS))}")
    if not src.is_dir():
        raise ValueError(f"dataset_root が存在しません: {src}")
    if dest.resolve() == src.resolve():
        raise ValueError(
            "出力先が入力フォルダと同じです。混在を解消できないうえ、"
            "元フォルダを変更することになります。別のフォルダを指定してください。")

    all_entries = sorted(src.iterdir())
    file_entries = [p for p in all_entries if p.is_file()]
    # 随伴の探索はファイルだけを見る。フォルダ形式の raw は自己完結しているので
    # 随伴を持たない（`COMPANION_RULES` にも登録が無い）。
    primaries = [p for p in all_entries
                 if p.suffix.lower().lstrip(".") == ext and is_raw_input(p)]
    if not primaries:
        raise ValueError(
            f"入力フォルダに .{ext} がありません: {src}  "
            "list_data_files(all_files=True) で実際の拡張子を確認してください。")

    to_copy: list[Path] = []
    companions = 0
    for primary in primaries:
        to_copy.append(primary)
        found = _companions_of(file_entries, primary, ext)
        companions += len(found)
        to_copy.extend(found)

    dest.mkdir(parents=True, exist_ok=True)
    linked = copied = skipped = 0
    for source in to_copy:
        target = dest / source.name
        if source.is_dir():
            # 既にあっても中まで歩く（途中で止まった前回の続きを埋める）。
            d_linked, d_copied, d_skipped = _replicate_dir(source, target)
            linked += d_linked
            copied += d_copied
            skipped += d_skipped
            continue
        if target.exists():
            skipped += 1
            continue
        try:
            os.link(source, target)
            linked += 1
        except OSError:
            # ボリュームをまたぐ／リンク非対応のファイルシステム。
            shutil.copy2(source, target)
            copied += 1

    if linked and copied:
        mode = "mixed"
    elif linked:
        mode = "hardlink"
    elif copied:
        mode = "copy"
    else:
        mode = "none"

    return PrepareResult(
        out_dir=str(dest),
        primary=len(primaries),
        companions=companions,
        linked=linked,
        copied=copied,
        skipped_existing=skipped,
        mode=mode,
    )
