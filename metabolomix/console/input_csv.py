"""MS-DIAL Console の `-i` に渡す試料の表（CSV）。

`-i` にフォルダを渡すと、Console（`AnalysisFilesParser.ReadFolderContents`）は各ファイルに別の
Class ID（"0"、"1"、…）を振り、種類はすべて Sample にする。整列の `N percent detected in one group`
は群ごとの本数で閾値を決めるので、1 本ずつの群では効かない。CSV（`ReadCsvContents`）なら群
（`class_id`）・種類（`file_type`）・注入順（`analytical_order`）・バッチ（`batch_order`）を渡せる。

上流の読み方（MsdialWorkbench afd5f9522 の `tests/MSDIAL5/MsdialCoreTestApp/Parser/AnalysisFilesParser.cs`
と `Process/CommonProcess.cs`。現用の Console のビルドにも同じ見出しの文字列がある）:

- ASCII で読み、`#` で始まる行を飛ばし、最初の行を見出しにして `,` で割るだけ（引用符を解さない）。
  空の行は file_path が空になって落ちるので書かない。
- 空の欄は空のまま使う（`file_name` と `class_id` は既定に落ちない）ので、どの列も必ず書く。
- `file_type` は Sample / Standard / QC / Blank、`acquisition_type` は DDA / SWATH / AIF（大小は問わない。
  読めなければ Sample・DDA）。フォルダ入力ではメソッドの `Acquisition type`（無ければ DDA）を全ファイルに
  当てるが、CSV では列の値をそのまま使うので、メソッドの値を書く。
- 行の順が file ID（0 始まり）。メソッドの `Alignment reference file ID` はこの添字を指す。
- `.mddata`（プロジェクト保存）は CSV のフォルダに、アラインメント（`.arf` ほか）は先頭のファイルの
  フォルダに書かれる。フォルダ入力と同じ場所に出すため、CSV は raw と同じフォルダに置く。
- `file_type` が QC / Blank だと、メソッドの `Remove feature based on peak height fold-change` などの
  種類を見る絞り込みが効くようになる（フォルダ入力ではすべて Sample なので効かなかった）。

stdlib と `core.atomic_io` だけに依存する（pipeline の worker が import する）。
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

from metabolomix.core.atomic_io import DomainError, atomic_write_bytes

__all__ = ["CSV_FILENAME", "COLUMNS", "acquisition_type", "build_rows", "class_id",
           "unsafe_names", "write_csv"]

#: `-i` に渡す CSV の名前（raw と同じフォルダ＝解析フォルダに置く）。
CSV_FILENAME = "console-input.csv"

#: `ReadCsvContents` の見出し（この名前で列を探す）。
COLUMNS = ("file_path", "file_name", "file_type", "class_id", "acquisition_type",
           "batch_order", "analytical_order", "factor")

#: 実験情報シートの role → Console の `AnalysisFileType`。`unknown` と空は Sample（フォルダ入力と同じ）。
_FILE_TYPES = {"sample": "Sample", "qc": "QC", "blank": "Blank", "standard": "Standard"}

#: メソッドの `Acquisition type` の値（`ConfigParser` が受ける 3 つ）→ CSV に書く値。
_ACQUISITION_TYPES = {"dda": "DDA", "swath": "SWATH", "aif": "AIF"}
ACQUISITION_TYPE_KEY = "Acquisition type"


def acquisition_type(method_keys: dict[str, str]) -> str:
    """メソッドの `Acquisition type`（`read_method_keys` の小文字キー）。無い・読めない値は DDA。

    Console はフォルダ入力のとき、この値（無ければ DDA）を全ファイルに当てる（`CommonProcess.SetProjectProperty`）。
    """
    value = (method_keys.get(ACQUISITION_TYPE_KEY.lower()) or "").strip().lower()
    return _ACQUISITION_TYPES.get(value, "DDA")


def class_id(group: str | None, role: str | None) -> str:
    """CSV の `class_id`: シートの群。群が空なら役割（`qc` / `blank` / `sample` …）でまとめる。

    ASCII で書けない・`,` や `\\` を含む名前は `unicode_escape` で書く（`,` は `\\x2c`）。素のまま書く名前は
    `\\` を含まず、書き換えた名前は必ず `\\` を含むので、違う名前が同じ群にまとまらない。
    """
    label = group if group else (role or "sample")
    if label.isascii() and label.isprintable() and "," not in label and "\\" not in label:
        return label
    return label.encode("unicode_escape").decode("ascii").replace(",", "\\x2c")


def _name_without_extension(name: str) -> str:
    """.NET の `Path.GetFileNameWithoutExtension` と同じ（最後の `.` から後ろを落とす）。"""
    cut = name.rfind(".")
    return name[:cut] if cut >= 0 else name


def _csv_safe(text: str) -> bool:
    return text.isascii() and not any(ch in text for ch in ",\r\n")


def unsafe_names(raw_paths) -> list[str]:
    """CSV に書けない（ASCII でない・`,` や改行を含む）パスの raw の名前。

    Console は CSV を ASCII で読み `,` で割るだけなので、こうしたパスは CSV では渡せない。
    """
    return [Path(p).name for p in raw_paths
            if not (_csv_safe(str(p)) and _csv_safe(_name_without_extension(Path(p).name)))]


def _norm(path) -> str:
    return os.path.normcase(os.path.normpath(str(path)))


def _batch_numbers(labels: list, order_keys: list) -> tuple[list[int], dict[str, int]]:
    """バッチの名前を、注入の順に最初に出てきた順で 1, 2, … にする。どれも空なら全部 1。"""
    if all(label is None for label in labels):
        return [1] * len(labels), {}
    numbers: dict[str | None, int] = {}
    for index in sorted(range(len(labels)), key=lambda i: order_keys[i]):
        numbers.setdefault(labels[index], len(numbers) + 1)
    return ([numbers[label] for label in labels],
            {("" if label is None else str(label)): n for label, n in numbers.items()})


def build_rows(raw_paths, manifest_rows: list[dict], *, acquisition: str) -> tuple[list[dict], dict]:
    """Console が読む raw の順に、CSV の行と記録用の要約を作る。

    `raw_paths` は Console に渡す raw の絶対パス（`job_manager.list_raw_inputs` の順。これが file ID の順）。
    `manifest_rows` は同じフォルダを基準に解した実験情報シートの行（`parse_manifest`。全 raw に 1 行ずつ）。
    `include` は見ない（フォルダ入力と同じく全 raw を Console に渡す）。

    - `analytical_order`: 全行に注入順があればシートの値、1 行でも欠ければ並びの順（1 始まり。Console の既定と同じ）。
    - `batch_order`: シートにバッチがあれば注入の順に 1, 2, …（空のバッチも 1 つに数える）、無ければ 1。
    """
    by_path = {_norm(row.get("source_file") or ""): row for row in manifest_rows}
    paths = [Path(p) for p in raw_paths]
    missing = [p.name for p in paths if _norm(p) not in by_path]
    if missing:
        raise DomainError("SAMPLE_MANIFEST_INVALID",
                          f"実験情報シートに行の無い raw があります: {', '.join(missing[:5])}",
                          {"missing_sources": missing})
    rows = [by_path[_norm(p)] for p in paths]

    orders = [row.get("injection_order") for row in rows]
    from_manifest = all(isinstance(o, int) and not isinstance(o, bool) for o in orders)
    analytical = orders if from_manifest else list(range(1, len(rows) + 1))
    batches, batch_map = _batch_numbers([row.get("batch") for row in rows],
                                        [(analytical[i], i) for i in range(len(rows))])

    out: list[dict] = []
    for index, (path, row) in enumerate(zip(paths, rows)):
        out.append({
            "file_path": str(path),
            "file_name": _name_without_extension(path.name),
            "file_type": _FILE_TYPES.get(row.get("role") or "sample", "Sample"),
            "class_id": class_id(row.get("group"), row.get("role")),
            "acquisition_type": acquisition,
            "batch_order": str(batches[index]),
            "analytical_order": str(analytical[index]),
            "factor": "1",
        })

    classes: dict[str, int] = {}
    file_types: dict[str, int] = {}
    for row in out:
        classes[row["class_id"]] = classes.get(row["class_id"], 0) + 1
        file_types[row["file_type"]] = file_types.get(row["file_type"], 0) + 1
    summary = {"n_files": len(out), "classes": classes, "file_types": file_types,
               "acquisition_type": acquisition,
               "analytical_order_source": "manifest" if from_manifest else "file_order"}
    if batch_map:
        summary["batch_order"] = batch_map
    return out, summary


def write_csv(path: Path, rows: list[dict]) -> str:
    """CSV を ASCII・CRLF で原子的に書き、sha256 を返す。"""
    lines = [",".join(COLUMNS)] + [",".join(row[column] for column in COLUMNS) for row in rows]
    data = ("\r\n".join(lines) + "\r\n").encode("ascii")
    atomic_write_bytes(Path(path), data)
    return hashlib.sha256(data).hexdigest()
