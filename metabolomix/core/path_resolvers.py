"""データファイルのパス解決とバッチ選択（ファイル名の処理タイムスタンプ基準）。

依存は mcp_core（DATA_DIR）・user_config（外部資産の場所）と stdlib のみの下位レイヤ。tools_* / server は import
しない。`list_data_files` はここに純関数として置き、MCP ツールとしての登録は上位
（metabolomix.tools.dataset）が担う — こうすることで resolve_* → list_data_files → DATA_DIR という
参照が下位で閉じ、metabolomix.tools.dataset との循環を避けられる。

DATA_DIR は実行時に差し替わるため `mcp_core.DATA_DIR` を動的参照する。
"""
import os
from pathlib import Path

from metabolomix.core import mcp_core
from metabolomix.core import user_config

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
    """同種ファイルが重複（旧版/新版）する場合に最新版のパスを返す。"""
    if not paths:
        return None
    return max(paths, key=_recency_key)


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


def _describe_batch_selection(directory: Path) -> str | None:
    """フォルダ内に複数バッチが混在する場合、最新バッチを自動選択した旨の注記を返す。

    バッチが1つ（または判別不能）なら None を返し、注記を出さない。
    """
    try:
        names = [f.name for f in directory.iterdir() if f.is_file()]
    except OSError:
        return None
    # 解析対象は .arf/.arf2 のみ。プロジェクト/メタ(.mddata/.mdproject/.msp2)や
    # per-sample(.pai2/.dcl)も AlignmentResult 形式のタイムスタンプを持つため、
    # 拡張子で解析対象に限定しないと告知バッチが実際に解析する .arf バッチとズレる
    # （例: POS で .mddata の 2026_06_17 を告知するが解析 .arf は 2024_06_13）。
    # 正規化キー -> 表示用タイムスタンプ（アンダースコア付きの読みやすい形）
    display: dict[str, str] = {}
    for name in names:
        low = name.lower()
        if not (low.endswith(".arf") or low.endswith(".arf2")):
            continue
        match = _ALIGNMENT_TIMESTAMP_RE.search(name)
        if match:
            display[match.group(1).replace("_", "")] = match.group(1)
    if len(display) <= 1:
        return None
    latest_key = max(display)
    n_old = len(display) - 1
    return (
        f"🗂️ フォルダ内に複数バッチ（{len(display)} 件の処理タイムスタンプ）を検出しました。"
        f"最新バッチ **{display[latest_key]}** を自動選択し、旧バッチ {n_old} 件はスキップします。"
    )


# このサーバのパーサが読める拡張子。MS-DIAL の出力フォルダには測定生データ
# （.wiff / .wiff.scan / .wiff2 / .timeseries.data / .txt）が同居し、実データでは
# 495 ファイル中 7 割以上がそれだった。既定でそこまで列挙すると、入口ツールの
# 戻り値だけで 46,977 字（≒1万数千トークン）を占める。
ANALYSABLE_EXTENSIONS: tuple[str, ...] = (
    ".arf", ".arf2", ".pai2", ".dcl", ".EIC.aef", ".mddata", ".mdproject",
)


def list_data_files(
    extension: str | None = None,
    directory: str | None = None,
    all_files: bool = False,
) -> list[str]:
    """データディレクトリ内のファイルの絶対パス一覧を返す（純関数）。

    - directory: 探索するディレクトリ。省略時は既定のデータディレクトリ
      (環境変数 LIPIDMIX_DATA_DIR または <project>/data)。
    - extension: 指定するとその拡張子（例 '.pai2', '.arf2', '.EIC.aef'）だけに絞る。
    - all_files: extension 未指定のとき、解析対象外の拡張子まで含めるか。

    見つからない場合・ディレクトリが無い場合は**空リスト**を返す。以前はここで
    エラー文面を1要素のリストとして返していたため、呼び出し側5箇所が
    「メッセージをパスとして掴まない」よう `os.path.isfile` で防御していた。
    """
    target_dir = Path(directory).expanduser() if directory else mcp_core.DATA_DIR
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


def _resolve_data_file(
    extension: str,
    file_path: str | None = None,
    prefer_suffix: str | None = None,
) -> str | None:
    """拡張子から解析対象ファイルを1つ選ぶ（全リゾルバの共通実装）。

    明示パスがあればそれを優先し、無ければ (1) 複数バッチ混在なら最新バッチへ絞り、
    (2) `prefer_suffix` に一致するものがあればそちらを優先し、(3) 重複時は最新版を
    採る。この 3 段は形式によらず同じなので、拡張子と優先条件だけを変えて共有する。

    `file_path` がフォルダなら、既定のデータディレクトリの代わりにその中を探す。
    LLM は `load_dataset` に渡したフォルダをそのまま `file_path` に渡しがちで、
    フォルダを返すと読み手が open して Windows では EACCES になる。フォルダの中に
    無いときは None を返し、既定ディレクトリの別データセットへは落ちない。
    """
    search_dir = None
    if file_path and os.path.isdir(file_path):
        search_dir = file_path
    elif file_path and os.path.exists(file_path):
        return file_path

    real_paths = [p for p in list_data_files(extension=extension, directory=search_dir)
                  if os.path.isfile(p)]
    if not real_paths:
        return None
    real_paths = _select_latest_batch(real_paths)  # 複数バッチ混在時は最新バッチへ
    if prefer_suffix:
        preferred = [p for p in real_paths if p.lower().endswith(prefer_suffix)]
        if preferred:
            return _pick_latest(preferred)
    return _pick_latest(real_paths)  # 重複時は最新版


def resolve_arf_file_path(file_path: str | None = None) -> str | None:
    """.arfファイルのパスを解決するヘルパー。

    MS-DIAL出力フォルダには DriftSpots.arf と PeakProperties.arf が併存しうるが、
    PCA等に使うサンプル別強度を持つのは **PeakProperties.arf** の方。両者がある場合は
    PeakProperties.arf を自動選択する（無ければ先頭にフォールバック）。
    """
    return _resolve_data_file(".arf", file_path, prefer_suffix="peakproperties.arf")


def resolve_arf2_file_path(file_path: str | None = None) -> str | None:
    """.arf2ファイルのパスを解決するヘルパー"""
    return _resolve_data_file(".arf2", file_path)


def resolve_eicaef_file_path(file_path: str | None = None) -> str | None:
    """EIC.aefファイルのパスを解決するヘルパー"""
    return _resolve_data_file(".EIC.aef", file_path)


def resolve_pai2_file_path(file_path: str | None = None) -> str | None:
    """.pai2ファイルのパスを解決するヘルパー。

    複数日付（複数バッチ）が混在していても最新バッチの .pai2 を自動選択する。
    """
    return _resolve_data_file(".pai2", file_path)


def resolve_dcl_file_path(file_path: str | None = None) -> str | None:
    """.dcl ファイルのパスを解決するヘルパー（MSDecResult / MS-MS 本体）。

    .dcl は測定ファイル1つにつき1個あるため、指定なしでは最新バッチの先頭を返す。
    特定サンプルの MS/MS が欲しいときは file_path を明示するか、.pai2 と同名の
    兄弟ファイルを引く dcl_reader.find_dcl_for_pai2 を使う。
    """
    return _resolve_data_file(".dcl", file_path)


#: 極性 → 研究室参照ライブラリ（`.msp`）の設定キー。値は環境変数か
#: `lipidmix.local.toml` の `[library]` で指す（`metabolomix.core.user_config`）。
#: ライブラリ本体はリポジトリの外（外部流出禁止の資産）に置き、場所だけを教える。
LIBRARY_SETTING_KEYS = {"positive": "library.msp_positive", "negative": "library.msp_negative"}
#: 同じ値を指す環境変数（正準は user_config.SETTINGS。既存の参照のため名前を残す）。
LIBRARY_ENV_VARS = {mode: user_config.SETTINGS[key] for mode, key in LIBRARY_SETTING_KEYS.items()}


class LibraryPathError(Exception):
    """参照ライブラリを 1 つに決められない（または指定先が無い）。

    `code` は機械可読（`LIBRARY_NOT_FOUND` / `MSP_ENV_NOT_FOUND` /
    `MSP_AMBIGUOUS` / `INVALID_ION_MODE` / `CONFIG_INVALID`）。`message` と `details`
    にはファイル名・設定キー・環境変数名・設定ファイルのパスだけを載せ、ライブラリの
    置き場所（ディレクトリ）は載せない——戻り値は LLM の文脈に入る。
    `MSP_ENV_NOT_FOUND` は値が設定ファイルから来た場合も同じコード（互換のため据え置き）。
    """

    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details if details is not None else {}


def _library_setting(ion_mode: str) -> user_config.Setting | None:
    """極性の設定（環境変数 → 設定ファイル）。設定ファイルが読めなければ CONFIG_INVALID。"""
    try:
        return user_config.get_setting(LIBRARY_SETTING_KEYS[ion_mode])
    except user_config.ConfigInvalidError as exc:
        raise LibraryPathError(exc.code, exc.message, exc.details()) from exc


def _library_from_setting(ion_mode: str) -> str | None:
    """極性の設定が指すパス。未設定なら None、指す先が無ければ例外。"""
    setting = _library_setting(ion_mode)
    if setting is None:
        return None
    if not os.path.isfile(setting.value):
        label = user_config.setting_label(setting)
        raise LibraryPathError(
            "MSP_ENV_NOT_FOUND",
            f"{label} が指すファイル（{os.path.basename(setting.value)}）がありません。"
            f"{label} を確認してください。",
            user_config.describe_missing(LIBRARY_SETTING_KEYS[ion_mode], setting),
        )
    return setting.value


def resolve_library_path(file_path: str | None = None, *, ion_mode: str | None = None) -> str | None:
    """参照ライブラリのパスを解決するヘルパー。

    順に:

    1. `file_path` の明示（無いファイルなら `LIBRARY_NOT_FOUND`。以前は data
       ディレクトリの探索へ黙って落ちて別のライブラリを掴んでいた）。
    2. `ion_mode` を指定したら、その極性の設定（環境変数 `MSDIAL_MSP_POS` /
       `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の `[library] msp_positive` /
       `msp_negative`。`LIBRARY_SETTING_KEYS`）。未設定なら 3 以降へ落ちる。
    3. data ディレクトリの `*_Loaded.msp2.dbs`（その run が実際に使った参照）。
    4. 極性の設定。1 つだけ設定されていればそれ、両方なら `MSP_AMBIGUOUS`
       （`ion_mode` の指定を求める）。出どころ（環境変数か設定ファイルか）は問わない。
    5. data ディレクトリの `*.msp`。複数あれば `MSP_AMBIGUOUS`——pos / neg の
       ように並ぶファイルを更新日時で黙って選ぶと、極性違いで照合しても候補が
       少し減るだけで誤りに気づけない。

    何も見つからなければ None。**`*.msp2` と `*.lbm2` は候補にしない**
    ——前者は ASCII `.msp` を指定したときだけ中身が入るので0バイトのことがあり
    （脂質経路では常に0）、後者は脂質専用の in-silico ライブラリで本機能
    （親水性メタボロミクス）の入口としては出さない。`.dbs` は複数バッチ混在時に
    最新バッチへ絞り、同一種別内の重複は最新版（`_pick_latest`）を採る。
    """
    if ion_mode is not None:
        normalized = ion_mode.strip().lower()
        if normalized not in LIBRARY_ENV_VARS:
            raise LibraryPathError(
                "INVALID_ION_MODE",
                f"ion_mode は {' / '.join(LIBRARY_ENV_VARS)} のどちらかを指定してください"
                f"（受け取った値: {ion_mode!r}）。",
            )
        ion_mode = normalized

    if file_path:
        if os.path.isfile(file_path):
            return file_path
        raise LibraryPathError(
            "LIBRARY_NOT_FOUND",
            f"指定された参照ライブラリ（{os.path.basename(file_path)}）がありません。",
        )

    if ion_mode is not None:
        from_setting = _library_from_setting(ion_mode)
        if from_setting:
            return from_setting

    dbs_paths = [p for p in list_data_files(extension=".msp2.dbs") if os.path.isfile(p)]
    if dbs_paths:
        dbs_paths = _select_latest_batch(dbs_paths)
        return _pick_latest(dbs_paths)

    settings = {mode: _library_setting(mode) for mode in LIBRARY_SETTING_KEYS}
    configured = [mode for mode, setting in settings.items() if setting is not None]
    if len(configured) > 1:
        names = " / ".join(user_config.setting_label(settings[mode]) for mode in configured)
        raise LibraryPathError(
            "MSP_AMBIGUOUS",
            f"参照ライブラリが極性ごとに設定されています（{names}）。"
            f"ion_mode（{' / '.join(configured)}）を指定してください。",
        )
    if configured:
        return _library_from_setting(configured[0])

    msp_paths = sorted(p for p in list_data_files(extension=".msp") if os.path.isfile(p))
    if not msp_paths:
        return None
    if len(msp_paths) > 1:
        names = ", ".join(os.path.basename(p) for p in msp_paths)
        raise LibraryPathError(
            "MSP_AMBIGUOUS",
            f"データディレクトリに .msp が複数あります（{names}）。file_path で 1 つ指定してください。",
        )
    return msp_paths[0]


def _filter_arf_spots(
    features: list[dict],
    min_intensity: float = 0.0,
    annotation_keyword: str | None = None,
    *,
    catalog: list[dict] | None = None,
    spot_ids: list[int] | None = None,
    ontologies: list[str] | None = None,
) -> list[dict]:
    """同一バッチの統合注釈で選択する。元スポットとピーク値は変更しない。"""
    if spot_ids is not None and (not spot_ids or any(type(i) is not int or i < 0 for i in spot_ids)):
        raise ValueError("spot_idsには0以上の整数を1件以上指定してください。")
    if ontologies is not None and (not ontologies or any(not isinstance(v, str) or not v.strip() for v in ontologies)):
        raise ValueError("ontologiesには空でない脂質クラス名を1件以上指定してください。")
    ids = set(spot_ids) if spot_ids is not None else None
    classes = {v.strip().casefold() for v in ontologies} if ontologies is not None else None
    by_id = {}
    for entry in catalog or []:
        sid = entry.get("MasterAlignmentID")
        if sid in by_id:
            raise ValueError(f"ARF2のMasterAlignmentIDが重複しています: {sid}")
        by_id[sid] = entry
    filtered_spots = []
    keyword = annotation_keyword.lower() if annotation_keyword else None
    for spot in features:
        if ids is not None and spot.get("MasterAlignmentID") not in ids:
            continue
        if catalog is not None:
            spot = dict(spot)
            entry = by_id.get(spot.get("MasterAlignmentID"), {})
            arf_name = spot.get("Name") or ""
            arf2_name = entry.get("Name") or ""
            has_name = bool(arf2_name.strip()) and arf2_name.strip().lower() != "unknown"
            spot["Name"] = arf2_name if has_name else arf_name
            spot["Ontology"] = entry.get("Ontology") or spot.get("Ontology")
            spot["arf_name"] = arf_name
            spot["arf2_name"] = arf2_name
            spot["annotation_source"] = "arf2" if has_name else "arf"
            spot["annotation_conflict"] = bool(has_name and arf_name != arf2_name)
        if classes is not None and (spot.get("Ontology") or "").strip().casefold() not in classes:
            continue
        height = spot.get("HeightAverage")
        if height is not None and height < min_intensity:
            continue

        if keyword:
            name = spot.get("Name", "")
            if not name or keyword not in name.lower():
                continue

        filtered_spots.append(spot)
    return filtered_spots
