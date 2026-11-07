"""データファイルのパス解決とバッチ選択（ファイル名の処理タイムスタンプ基準）。

依存は mcp_core・user_config と stdlib のみの下位レイヤ。tools_* / server は import しない。
`list_data_files` はここに純関数として置き、MCP ツールとしての登録は上位
（metabolomix.tools.dataset）が担う。

**データの探索先（旧 `DATA_DIR`）は無い**（B-2b）。フォルダは呼び出し側が渡し、パスを省いた
ツールは直近のデータセットのファイルを使う（`metabolomix.dataset.artifacts`）。バッチ選択の
関数（`_select_latest_batch` ほか）は解析フォルダのアダプタ（`dataset/alignment_files.py`）が使う。
"""
import os
from pathlib import Path


# MS-DIALのアライメント結果ファイル名に埋め込まれる処理タイムスタンプ。
# 例: AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf
#     → 再アライメントすると新しいタイムスタンプのセットが増える（＝旧版/新版の重複）。
import re

_ALIGNMENT_TIMESTAMP_RE = re.compile(r"(\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2})")
# サンプル名等に付く12〜14桁の連続タイムスタンプ（例: _202605151012）も拾う。
_COMPACT_TIMESTAMP_RE = re.compile(r"(\d{12,14})")


def _recency_key(path: str) -> tuple[str, float]:
    """ファイルの「新しさ」の並べ替えキー。

    第一に**ファイル名に埋め込まれた処理タイムスタンプ**（コピーでも保たれる）、
    第二に更新時刻(mtime)。タイムスタンプ無しは空文字となり mtime で比較される。
    """
    name = os.path.basename(path)
    match = _ALIGNMENT_TIMESTAMP_RE.search(name)
    if match:
        timestamp = match.group(1).replace("_", "")
    else:
        compact = _COMPACT_TIMESTAMP_RE.search(name)
        timestamp = compact.group(1) if compact else ""
    try:
        mtime = os.path.getmtime(path)
    except OSError:
        mtime = 0.0
    return (timestamp, mtime)


def _pick_latest(paths: list[str]) -> str | None:
    """同種ファイルが重複（旧版/新版）する場合に最新版のパスを返す。

    時刻を持つ候補があれば、時刻の無い候補より先に選ぶ（時刻の無い候補どうしは mtime）。
    """
    if not paths:
        return None
    return max(_timestamped_first(paths), key=_recency_key)


def _timestamped_first(paths: list[str]) -> list[str]:
    """ファイル名に処理時刻を持つ候補があればそれだけ、無ければ全件。

    `_select_latest_batch` は時刻の無いファイルを温存するので、その後に種別（PeakProperties）の
    優先を当てると、時刻の無い PeakProperties が時刻つきの別バッチに勝ってしまう。
    種別の優先はこの絞り込みのあとに当てる。
    """
    stamped = [p for p in paths if _batch_key(p)]
    return stamped or list(paths)


def _batch_key(path: str) -> str:
    """ファイル名に埋め込まれた処理タイムスタンプ（＝バッチ識別子）を返す。

    MS-DIAL は1回の処理で生成する全ファイルに同一の `AlignmentResult_<timestamp>`
    接頭辞を付ける。そのタイムスタンプ（無ければ12-14桁連番）を正規化して返す。
    どちらも持たないファイルは空文字（＝バッチ不明）となる。mtime は見ない
    （コピーでも保たれるファイル名側の識別子だけでバッチを束ねるため）。
    """
    name = os.path.basename(path)
    match = _ALIGNMENT_TIMESTAMP_RE.search(name)
    if match:
        return match.group(1).replace("_", "")
    compact = _COMPACT_TIMESTAMP_RE.search(name)
    return compact.group(1) if compact else ""


def _select_latest_batch(paths: list[str]) -> list[str]:
    """複数バッチ（処理タイムスタンプ）が混在する場合に最新バッチへ絞る。

    - 埋め込みタイムスタンプを持つファイルがあれば、その最大値に一致する
      ファイル群だけを残す（＝旧バッチを除外）。
    - タイムスタンプを持たないファイルはバッチ判定不能なので除外せず温存する
      （誤って解析対象を失わない安全側）。
    - 全ファイルが無タイムスタンプなら全件そのまま返す（現状互換）。
    """
    if not paths:
        return []
    keyed = [(p, _batch_key(p)) for p in paths]
    timestamps = [k for _, k in keyed if k]
    if not timestamps:
        return list(paths)
    latest = max(timestamps)
    return [p for p, k in keyed if k == latest or not k]


# このサーバのパーサが読める拡張子。MS-DIAL の出力フォルダには測定生データ
# （.wiff / .wiff.scan / .wiff2 / .timeseries.data / .txt）が同居し、実データでは
# 495 ファイル中 7 割以上がそれだった。既定でそこまで列挙すると、入口ツールの
# 戻り値だけで 46,977 字（≒1万数千トークン）を占める。
ANALYSABLE_EXTENSIONS: tuple[str, ...] = (
    ".arf", ".arf2", ".pai2", ".dcl", ".EIC.aef", ".mddata", ".mdproject",
)


def list_data_files(
    directory: str,
    extension: str | None = None,
    all_files: bool = False,
) -> list[str]:
    """フォルダ内のファイルの絶対パス一覧を返す（純関数）。

    - directory: 探索するフォルダ（必須。既定の探索先は無い）。
    - extension: 指定するとその拡張子（例 '.pai2', '.arf2', '.EIC.aef'）だけに絞る。
    - all_files: extension 未指定のとき、解析対象外の拡張子まで含めるか。

    見つからない場合・ディレクトリが無い場合は**空リスト**を返す。以前はここで
    エラー文面を1要素のリストとして返していたため、呼び出し側5箇所が
    「メッセージをパスとして掴まない」よう `os.path.isfile` で防御していた。
    """
    target_dir = Path(directory).expanduser()
    if not target_dir.is_dir():
        return []

    if extension is not None:
        suffixes: tuple[str, ...] = (extension,)
    elif all_files:
        suffixes = ()
    else:
        suffixes = ANALYSABLE_EXTENSIONS

    # フォルダ形式の計測データ（Agilent/Bruker の `.d`・Waters の `.raw`）は
    # フォルダそのものが 1 検体。一律に is_file() で落とすと、その形式だけが
    # 入った生データフォルダが「空」に見え、入口で「生データが無い」と誤読される。
    # 既定（ANALYSABLE_EXTENSIONS）には現れない ——`.d` はこのサーバのパーサが
    # 読める形式ではなく、Console 経路への入力だから。
    include_vendor_dirs = extension is not None or all_files
    # 「どのフォルダが計測データか」の正準は console 層（上流の規則をそこで
    # 写している）。同じ規則をここに書き写すとドリフトするので借りる。
    # module 先頭に置かないのは、このモジュールを mcp_core と stdlib だけに
    # 依存する下位レイヤのまま保つため（冒頭の docstring の約束）。
    from metabolomix.console.job_manager import is_raw_input

    file_paths = []
    for file in sorted(target_dir.iterdir()):
        if not file.is_file():
            if not (include_vendor_dirs and is_raw_input(file)):
                continue
        if suffixes and not str(file).endswith(suffixes):
            continue
        file_paths.append(str(file.absolute()))
    return file_paths


class LibraryPathError(Exception):
    """参照ライブラリの指定先が無い。

    `code` は機械可読（`LIBRARY_NOT_FOUND`）。`message` にはファイル名だけを載せ、ライブラリの
    置き場所（ディレクトリ）は載せない——戻り値は LLM の文脈に入る。
    """

    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details if details is not None else {}


def resolve_library_path(file_path: str | None = None) -> str | None:
    """参照ライブラリの明示パスを確かめる。省略（None）なら None。

    **探索はしない**（旧版はデータディレクトリの `*_Loaded.msp2.dbs` / `*.msp` と極性の
    設定を探した。B-2b で消した）。省いたときの既定はデータセットが覚えたライブラリで、
    `library_load` が決める。無いファイルなら `LIBRARY_NOT_FOUND`（黙って別のライブラリへ
    落ちない）。
    """
    if not file_path:
        return None
    if os.path.isfile(file_path):
        return file_path
    raise LibraryPathError(
        "LIBRARY_NOT_FOUND",
        f"指定された参照ライブラリ（{os.path.basename(file_path)}）がありません。",
    )
