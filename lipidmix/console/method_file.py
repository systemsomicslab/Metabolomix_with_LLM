"""MS-DIAL Console のメソッドファイルを解決する純ロジック。

**なぜこの層が要るか**（docs/HISTRY.md 2026-09-04(2)）:

MS-DIAL 5 の GUI と Console は、脂質ライブラリ（LBM）の持ち方が違う。

- GUI: `DataBaseModel` のリストで持ち、アプリフォルダ
  （`Assembly.GetExecutingAssembly().Location` の親）の `*.lbm?` を 1 件
  自動で拾う（`MethodSettingModelFactory` / `DataBaseSettingModel.TrySetLbmLibrary`）。
  ユーザーには選ばせない。
- Console: レガシーの単一フィールド `param.LbmFilePath` しか読まない
  （`MsdialCoreTestApp/Process/CommonProcess.cs` の
  `if (ErrorHandler.IsFileExist(param.LbmFilePath))`）。

そして GUI は `param.LbmFilePath` に**一度も代入しない**（MSDIAL5 の src 全体で
このフィールドへの代入はプロパティ自身の setter だけ）。`ParameterToString()` は
そのレガシーフィールドを書き出すので、GUI 由来のパラメータは
**`Lbm file path:` が構造的に必ず空**になる。空のまま Console に渡すと、
警告もエラーも出さずに同定 0 件で完走する。

この層は GUI と同じ規則を再現して段差を埋める。**GUI が拒否する状況でだけ拒否する**
（アプリフォルダの `*.lbm?` がちょうど 1 件でなければ GUI も MessageBox で止める。
`DatasetParameterSettingModel.Prepare`）。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass, replace
from pathlib import Path

from lipidmix.core.atomic_io import DomainError
from lipidmix.core.user_config import Setting, setting_label

# GUI の DataBaseSettingViewModel が使う判定と同じ（`@"\.lbm\d*"`）。
# .NET の `GetFiles(dir, "*.lbm?")` が拾う .lbm / .lbm2 に一致し、.lbmx は拾わない。
_LBM_SUFFIX = re.compile(r"\.lbm\d*$", re.IGNORECASE)

# GUI の AutoParametersSave が付ける名前は `<project>_param_<yyyyMMddHHmm>.txt`。
_PARAM_FILENAME = re.compile(r"_param_\d+\.txt$", re.IGNORECASE)

LBM_KEY = "Lbm file path"

#: metabolomics向け依存キー（spec §5.1）。MS-DIAL 5 Console
#: (`tests/MSDIAL5/MsdialCoreTestApp`) のソースで実際に読まれることを確認済み
#: （証拠: MsdialWorkbench コミット afd5f9522fa2f11990e1ad51f88b22eef00a087c,
#: 2026-09-08 時点 `master`）。
#:
#: - `MSP_KEY`: `ParameterBase.MspFilePath` / `ConfigParser.cs` の
#:   `case "msp file path": param.MspFilePath = value;`。`CommonProcess.cs`が
#:   `LibraryHandler.ReadMsLibrary(param.MspFilePath, ...)` で読み、
#:   `DataBaseSource.Msp` として同定に使う。
#: - `TEXT_DB_KEY`: `ParameterBase.TextDBFilePath` / `ConfigParser.cs` の
#:   `case "text db file path": param.TextDBFilePath = value;`。同じく
#:   `CommonProcess.cs`が`LibraryHandler.ReadMsLibrary(param.TextDBFilePath, ...)`
#:   で読み、`DataBaseSource.Text` として同定に使う（msp/lbmと並ぶテキスト形式の
#:   同定用データベース——`lcms-profile.v1`の`kind="text_identification"`に対応）。
#: - `RT_REFERENCE_KEY`: `ParameterBase.CompoundListForRtCorrectionPath` /
#:   `ConfigParser.cs`の`case "compounds library file path for rt correction"`。
#:   `RetentionTimeCorrectionProcess.cs`が「RT correction anchor library」として
#:   読み、宣言が無い・ファイルが無い場合はConsole自身がそこで停止する
#:   （`kind="rt_reference"`に対応）。
#:
#: 見つかったが対応させなかったキー: `IsotopeTextDBFilePath`
#: （"Isotope text DB file path"）はアイソトープ追跡専用で、
#: `msp/lbm/text_identification/rt_reference`のどれにも該当しないため
#: 未登録（推測でkindへ割り当てない）。`CompoundListInTargetModePath`
#: （"Compounds library file path for target detection"）はターゲットモード
#: 検出専用（同定ではなく検出）で、同定用の`text_identification`とは役割が違うため
#: 別扱いのまま残す。
MSP_KEY = "Msp file path"
#: 表示名は`ParameterBase.ParameterToString()`の`"Text DB file path"`
#: （`ParameterBase.cs:535`）に合わせる。`ConfigParser`の読込は
#: `.lower()`後の比較なので大文字小文字は実害がないが、証拠に合わせておく。
TEXT_DB_KEY = "Text DB file path"
RT_REFERENCE_KEY = "Compounds library file path for RT correction"

#: LC-MS の Console が**ファイルパスとして**読むキー（`ConfigParser.ReadCommonParameter`
#: の `//File paths` 節）。`ReadForLcmsParameter` はこれらを解決しないので、相対値は
#: Console プロセスの cwd 基準になる——メソッド基準で解くのは GC-MS 経路の
#: `ResolveGcmsFilePaths` → `ResolvePathFromMethodFile` だけ（MsdialWorkbench
#: afd5f95）。見つからなければ `CommonProcess.ParseLibraries` が黙って飛ばす
#: （`IsFileExist` 判定。警告もエラーも出ない）。
CONSOLE_PATH_KEYS: tuple[str, ...] = (
    LBM_KEY, MSP_KEY, TEXT_DB_KEY,
    "Isotope text DB file path",
    "Compounds library file path for target detection",
    RT_REFERENCE_KEY,
    "RT correction peak selection file path",
)

#: annotator 設定表（TSV）を指すキー。こちらは Console 自身が**渡されたメソッド
#: ファイルの親**基準で解く（`ConfigParser.ReadMspAnnotatorSettings` /
#: `ReadTextAnnotatorSettings`）。実効コピーを別フォルダに書くと基準がずれるので、
#: 原本基準の絶対パスに固定する。綴りの別名は Console の `case` 列挙どおり。
SETTINGS_PATH_KEYS: tuple[str, ...] = (
    "MSP annotator settings file path",
    "MSP annotation settings file path",
    "MSP search settings file path",
    "Text annotator settings file path",
    "Text library annotator settings file path",
    "Text DB annotator settings file path",
    "Text annotation settings file path",
)


def _absolute(path: Path) -> str:
    """絶対パス文字列にする。シンボリックリンクやドライブ割当ては解かない。"""
    return os.path.abspath(path)


def relative_path_overrides(method_keys: dict[str, str], method_file: Path) -> dict[str, str]:
    """相対で宣言されたパスキーを、原本メソッドの親基準の絶対パスへ書き換える上書きを返す。

    実効メソッドは原本と別フォルダに書かれ、Console はそれを run_dir を cwd に
    して読む。相対値のままだと、`CONSOLE_PATH_KEYS` は cwd 基準、
    `SETTINGS_PATH_KEYS` は実効コピー基準で読まれ、どちらも原本の意図とずれる。

    参照先の実在は問わない（止めるかどうかは呼び出し側の方針。既存 lipidomics
    経路は LBM 以外の古い宣言で止めない）。絶対・空のキーと未知のキーは触らない。
    """
    out: dict[str, str] = {}
    for key in (*CONSOLE_PATH_KEYS, *SETTINGS_PATH_KEYS):
        declared = (method_keys.get(key.lower()) or "").strip()
        if not declared or Path(declared).is_absolute():
            continue
        out[key] = _absolute(method_file.parent / declared)
    return out

ION_MODE_KEY = "Ion mode"
ADDUCT_KEY = "Searched adduct ions"

#: 極性ごとのアダクト標準セット。ラボの実パラメータ（2_lipidome_lcms/NEG と
#: 20240915_spleen/POS）から採った。両者は極性・アダクト・スレッド数以外は
#: 同一だったので、この 2 行の差し替えが極性変換のすべて。
STANDARD_ADDUCTS = {
    "positive": (
        "[M+H]+,[M+NH4]+,[M+Na]+,[M+CH3OH+H]+,[M+K]+,[M+Li]+,[M+ACN+H]+,[M+H-H2O]+,"
        "[M+H-2H2O]+,[M+2Na-H]+,[M+IsoProp+H]+,[M+ACN+Na]+,[M+2K-H]+,[M+DMSO+H]+,"
        "[M+2ACN+H]+,[M+IsoProp+Na+H]+,[M-C6H10O4+H]+,[M-C6H10O5+H]+,[M-C6H8O6+H]+,"
        "[2M+H]+,[2M+NH4]+,[2M+Na]+,[2M+3H2O+2H]+,[2M+K]+,[2M+ACN+H]+,[2M+ACN+Na]+,"
        "[M+2H]2+,[M+H+NH4]2+,[M+H+Na]2+,[M+H+K]2+,[M+ACN+2H]2+,[M+2Na]2+,"
        "[M+2ACN+2H]2+,[M+3ACN+2H]2+,[M+3H]3+,[M+2H+Na]3+,[M+H+2Na]3+,[M+3Na]3+"),
    "negative": (
        "[M-H]-,[M-H2O-H]-,[M+Na-2H]-,[M+Cl]-,[M+K-2H]-,[M+HCOO]-,[M+CH3COO]-,"
        "[M+C2H3N+Na-2H]-,[M+Br]-,[M+TFA-H]-,[M-C6H10O4-H]-,[M-C6H10O5-H]-,"
        "[M-C6H8O6-H]-,[M+CH3COONa-H]-,[2M-H]-,[2M+FA-H]-,[2M+Hac-H]-,[3M-H]-,"
        "[M-2H]2-,[M-3H]3-"),
}


# メソッドファイルは実測 9KB 程度。誤って .mddata（GB 級）を渡されても
# 落ちないよう上限を置く（_looks_like_method_text が先に弾くが、単体でも安全に）。
_MAX_METHOD_BYTES = 1 << 20


@dataclass(frozen=True)
class LbmResolution:
    """LBM の解決結果。`error_code` が非 None なら呼び出し側は停止する。"""

    path: str | None
    source: str  # argument | method_file | build_tree | env | config_file | exe_dir | not_required
    error_code: str | None = None
    message: str | None = None
    candidates: tuple[str, ...] = ()


@dataclass(frozen=True)
class MethodCandidate:
    """GUI が自動保存した既存パラメータファイル 1 件。"""

    path: str
    ion_mode: str | None
    omics: str | None
    has_lbm: bool
    mtime: float
    # どこで見つかったか。UI のグルーピングと、同じパスが二重に出たときの優先に使う。
    origin: str = "same_dir"          # same_dir | sibling | past_run | given
    # 求める極性に対してそのまま使えるか。別極性は console_method_template を経由させる。
    # polarity を指定せずに探索したときは「一致するか」を判定しようがないので None
    # （「一致する」を意味する "direct" と取り違えられないよう、値を出さない）。
    usable: str | None = "direct"     # direct | needs_polarity_conversion | None
    key_params: dict[str, str] | None = None


# 候補どうしの差を読むのに要る少数キー。全キー（実測 287 行 / 11.7 KB）を候補ごとに
# 返すと戻り値が肥大する。綴りは実ファイル（param_POS_generated.txt）準拠。
KEY_PARAM_KEYS: tuple[str, ...] = (
    "Ion mode",
    "Target omics",
    "Minimum peak height",
    "Retention time begin",
    "Retention time end",
    "MS1 mass range begin",
    "MS1 mass range end",
    "MS1 tolerance for centroid",
    "Retention time tolerance for alignment",
    "MS1 tolerance for alignment",
    "Searched adduct ions",
)

# これを超える候補数では key_params を付けない。比較表は絞ってから引き直す。
KEY_PARAMS_MAX_CANDIDATES = 10

# 一覧として実際に返す候補数の上限。KEY_PARAMS_MAX_CANDIDATES とは意味の違う
# 別の定数にしてある — 1 つに束ねると、どちらか片方のつもりの変更がもう片方の
# 挙動（key_params を付けるかどうか）まで静かに動かしてしまう。実測でラボの
# レイアウト（兄弟フォルダ数十 × GUI 自動保存複数）は候補 240 件・69KB になり得る。
# ソート済み（direct 優先・新しい順）の先頭から切るので、有用な候補ほど残る。
MAX_REPORTED_CANDIDATES = 10

_ADDUCT_PREVIEW = 3


def extract_key_params(method_keys: dict[str, str]) -> dict[str, str]:
    """判断に効くキーだけを、実ファイルの綴りで取り出す。

    `Searched adduct ions` は POS の実値が 37 種・約 700 文字あるので要約する。
    極性の違いは先頭 3 種で判別できる。
    """
    out: dict[str, str] = {}
    for key in KEY_PARAM_KEYS:
        value = method_keys.get(key.lower())
        if value is None or value == "":
            continue
        if key == ADDUCT_KEY:
            items = [t.strip() for t in value.split(",") if t.strip()]
            head = ", ".join(items[:_ADDUCT_PREVIEW])
            out[key] = f"{len(items)} 種（先頭: {head}）"
        else:
            out[key] = value
    return out


def read_method_keys(path: Path) -> dict[str, str]:
    """`key: value` を小文字キーの辞書にして返す。

    MS-DIAL の `ConfigParser.ReadForLcmsParameter` と同じ割り方をする:
    `#` 始まりを飛ばし、**最初に現れる `:` か `=`** で 1 回だけ割る。
    値側の `:`（`C:\\...`）で割ってはいけないので `split(sep, 1)` にする。

    値が空の行も残す。`Lbm file path:` が「空で存在する」ことが、
    GUI 由来のパラメータを見分ける手がかりそのものだから。
    """
    try:
        raw = path.read_bytes()[:_MAX_METHOD_BYTES]
    except OSError:
        return {}
    if b"\x00" in raw:
        return {}
    text = raw.decode("ascii", errors="replace")

    keys: dict[str, str] = {}
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        positions = [i for i in (stripped.find(":"), stripped.find("=")) if i > 0]
        if not positions:
            continue
        cut = min(positions)
        key = stripped[:cut].strip().lower()
        if key:
            keys[key] = stripped[cut + 1:].strip()
    return keys


def find_lbm_files(directory: Path) -> list[Path]:
    """directory 直下の `*.lbm` / `*.lbm2` を返す（GUI と同じ TopDirectoryOnly）。"""
    try:
        entries = sorted(directory.iterdir())
    except OSError:
        return []
    return [p for p in entries if p.is_file() and _LBM_SUFFIX.search(p.name)]


# MsdialWorkbench のビルド生成物から .lbm2 を引くための座標。
# Console 実行体は tests/MSDIAL5/MsdialCoreTestApp/bin/Debug/<TFM>/MSDIALCUI.exe に
# 出るが、その exe フォルダに .lbm2 は無い。ライブラリは GUI アプリ側の
# src/MSDIAL5/MsdialGuiApp/bin/Debug/<TFM>/ に出るため、exe フォルダ探索
# （GUI と同じ TopDirectoryOnly）だけでは原理的に当たらない。
# MSDIAL4（src/MSDIAL4/MsDial/...）は見ない。あちらの conventional ライブラリは
# NCDK 無しの別世代で、MSDIAL5 の Console に食わせると同定結果が静かに変わる。
_BUILD_LBM_RELATIVE = ("src", "MSDIAL5", "MsdialGuiApp", "bin", "Debug")
_BUILD_TREE_MAX_ANCESTORS = 10


def _pick_build_lbm(debug_dir: Path, exe_tfm: str) -> tuple[Path | None, tuple[Path, ...]]:
    """`bin/Debug` 配下から 1 本選ぶ。

    同じライブラリが TFM ごとに複製されるので「候補が複数あるから決められない」
    とは扱わない（そう扱うと必ず LBM_AMBIGUOUS で止まる）。exe 自身の TFM に
    揃えるのが最も安全で、無ければ素の `Debug/` 直下、それも無ければ最新の
    mtime を採る。

    ただし、選んだ 1 フォルダの中に**ファイル名の異なる** `.lbm2` が複数あるときは
    TFM 複製ではなく別ライブラリなので、`found[0]`（アルファベット順）を黙って
    選ぶと識別結果を静かに変える。戻り値の第 2 要素にその候補を入れて返し、
    呼び出し側（`resolve_lbm`）で `exe_dir` 探索と同じ LBM_AMBIGUOUS にする。

    Returns: (選んだパス | None, 曖昧だったときの候補 | 空タプル)
    """
    for directory in (debug_dir / exe_tfm, debug_dir):
        found = find_lbm_files(directory)
        if found:
            if len({p.name for p in found}) > 1:
                return None, tuple(found)
            return found[0], ()

    others: list[Path] = []
    try:
        entries = sorted(debug_dir.iterdir())
    except OSError:
        return None, ()
    for entry in entries:
        if entry.is_dir():
            others.extend(find_lbm_files(entry))
    if not others:
        return None, ()
    return max(others, key=lambda p: p.stat().st_mtime), ()


def find_build_tree_lbm(exe_path: str | None) -> tuple[Path | None, tuple[Path, ...]]:
    """Console exe を起点に MsdialWorkbench のビルド生成物内の .lbm2 を返す。

    exe フォルダから上へ辿り、`src/MSDIAL5/MsdialGuiApp/bin/Debug` を持つ階層を
    リポジトリルートと見なす。見つからなければ (None, ())（＝ビルド運用ではない）。
    第 2 要素が非空なら、選んだフォルダに名前の異なる `.lbm2` が複数あり
    1 本に決められないことを示す（`_pick_build_lbm` 参照）。
    """
    if not exe_path:
        return None, ()
    exe_dir = Path(exe_path).expanduser().parent
    for ancestor in [exe_dir, *exe_dir.parents][:_BUILD_TREE_MAX_ANCESTORS + 1]:
        debug_dir = ancestor.joinpath(*_BUILD_LBM_RELATIVE)
        if debug_dir.is_dir():
            return _pick_build_lbm(debug_dir, exe_dir.name)
    return None, ()


def resolve_lbm(
    method_keys: dict[str, str],
    method_file: Path,
    omics: str,
    exe_path: str | None,
    lbm_setting: Setting | None = None,
    override: str | None = None,
) -> LbmResolution:
    """脂質ライブラリのパスを GUI と同じ規則で解決する。

    順に: 明示引数 → メソッドファイルの宣言 → **ビルド生成物** → `lbm_setting`
    （環境変数 MSDIAL_LBM → 設定ファイルの `[msdial] lbm`。呼び出し側が
    `user_config.get_setting("msdial.lbm")` で引いて渡す）→ Console の実行体と同じ
    フォルダ。ビルド生成物を `lbm_setting` より上に置くのは、
    この環境の Console がソースからのビルドで、ライブラリもそのツリー内の
    新しいものを使うため（インストール版より優先する）。
    返すパスは常に絶対パス。呼び出し側はそれを実効メソッドへ書き戻す
    （宣言が相対でも原本のままにしない。`relative_path_overrides` 参照）。

    明示引数とメソッドの宣言は omics を問わず採る——Console は `Target omics`
    に関わらず `Lbm file path` を読み（`CommonProcess.ParseLibraries`）、LBM
    annotator を `TargetOmics.Lipidomics` 固定で組んで脂質を同定する
    （`LcmsProcess`）。宣言されたのに見つからなければ止める（Console は黙って
    飛ばすので）。GUI 流の自動補完（ビルド生成物 → `lbm_setting`（`[msdial] lbm` / `MSDIAL_LBM`）→ exe フォルダ）
    だけを lipidomics に限る。metabolomics で宣言も引数も無ければ
    `not_required`（LBM を足さない）。
    """
    if override:
        candidate = Path(override).expanduser()
        if candidate.is_file():
            return LbmResolution(path=_absolute(candidate), source="argument")
        return LbmResolution(
            path=None, source="argument", error_code="LBM_NOT_FOUND",
            message=f"lbm_file が指すファイルがありません: {override}")

    declared = (method_keys.get(LBM_KEY.lower()) or "").strip()
    if declared:
        # LC-MS の Console（`ConfigParser.ReadForLcmsParameter`）は宣言を解決せず、
        # 相対値は Console プロセスの cwd 基準で読む。メソッド基準で解くのは
        # GC-MS 経路の `ResolvePathFromMethodFile` だけ（MsdialWorkbench afd5f95）。
        # ここではメソッド基準で解き、呼び出し側が絶対パスを実効メソッドへ書く。
        candidate = Path(declared)
        if not candidate.is_absolute():
            candidate = method_file.parent / candidate
        if candidate.is_file():
            return LbmResolution(path=_absolute(candidate), source="method_file")
        return LbmResolution(
            path=None, source="method_file", error_code="LBM_NOT_FOUND",
            message=(f"メソッドファイルが指す脂質ライブラリが見つかりません: {declared}  "
                     f"（{method_file} 基準で解決: {candidate}）"))

    if omics != "lipidomics":
        return LbmResolution(path=None, source="not_required")

    from_build, build_ambiguous = find_build_tree_lbm(exe_path)
    if build_ambiguous:
        return LbmResolution(
            path=None, source="build_tree", error_code="LBM_AMBIGUOUS",
            message=(
                f"ビルド生成物の脂質ライブラリ候補が {len(build_ambiguous)} 件あり、"
                f"どれを使うか決められません: {build_ambiguous[0].parent}  "
                "TFM ごとの複製（同名コピー）ではなく、名前の異なる .lbm2 が"
                "同じフォルダに複数あります。MS-DIAL GUI も 1 件でなければ実行を"
                "止めます。1 件だけ残すか、lbm_file 引数で明示してください。"),
            candidates=tuple(str(p) for p in build_ambiguous))
    if from_build is not None:
        return LbmResolution(path=_absolute(from_build), source="build_tree")

    if lbm_setting is not None:
        candidate = Path(lbm_setting.value)
        if candidate.is_file():
            return LbmResolution(path=_absolute(candidate), source=lbm_setting.source)
        return LbmResolution(
            path=None, source=lbm_setting.source, error_code="LBM_NOT_FOUND",
            message=f"{setting_label(lbm_setting)} が指すファイルがありません: {lbm_setting.value}")

    if not exe_path:
        return LbmResolution(
            path=None, source="exe_dir", error_code="LBM_NOT_FOUND",
            message="脂質ライブラリ（.lbm2）を解決できません。Console の実行体（[msdial] exe / MSDIAL_EXE）が未設定です。")

    exe_dir = Path(exe_path).expanduser().parent
    found = find_lbm_files(exe_dir)
    if len(found) == 1:
        return LbmResolution(path=_absolute(found[0]), source="exe_dir")
    if not found:
        return LbmResolution(
            path=None, source="exe_dir", error_code="LBM_NOT_FOUND",
            message=(
                "脂質ライブラリ（.lbm2）が見つかりません。MS-DIAL GUI は"
                "アプリフォルダの *.lbm2 を自動で使いますが、Console はメソッド"
                "ファイルの `Lbm file path:` しか読まず、GUI が書き出す"
                "パラメータはこの行が必ず空です。空のまま実行すると"
                "**警告なしで同定 0 件**になります。"
                f"探した場所: {exe_dir}  "
                "MS-DIAL のインストールフォルダにある .lbm2 のパスを設定ファイル"
                "（lipidmix.local.toml）の [msdial] lbm か環境変数 MSDIAL_LBM に"
                "設定するか、lbm_file 引数で渡してください。"))
    return LbmResolution(
        path=None, source="exe_dir", error_code="LBM_AMBIGUOUS",
        message=(
            f"脂質ライブラリの候補が {len(found)} 件あり、どれを使うか決められません: {exe_dir}  "
            "MS-DIAL GUI も 1 件でなければ実行を止めます。1 件だけ残すか、"
            "lbm_file 引数で明示してください。"),
        candidates=tuple(str(p) for p in found))


def scan_dir_for_method_files(directory: Path, origin: str) -> list[MethodCandidate]:
    """1 フォルダ直下の `*_param_<ts>.txt` を候補にする。絞り込みはしない。

    絞り込み（極性・omics）を呼び出し側に残すのは、別極性の候補を
    「使えないから消す」のではなく「変換が要る」と提示するため。
    """
    out: list[MethodCandidate] = []
    try:
        entries = sorted(Path(directory).iterdir())
    except OSError:
        return out
    for p in entries:
        if not p.is_file() or not _PARAM_FILENAME.search(p.name):
            continue
        keys = read_method_keys(p)
        if not keys:
            continue  # バイナリ／読めない
        out.append(MethodCandidate(
            path=str(p),
            ion_mode=(keys.get("ion mode") or "").strip().lower() or None,
            omics=(keys.get("target omics") or "").strip().lower() or None,
            has_lbm=bool((keys.get(LBM_KEY.lower()) or "").strip()),
            mtime=p.stat().st_mtime,
            origin=origin,
        ))
    return out


RUNS_SUBDIR_NAME = "runs"


def method_search_dirs(dataset_root: Path, search_dirs=None) -> list[tuple[Path, str]]:
    """メソッドファイルを探すフォルダを優先順に返す。

    `dataset_root` 直下 → 兄弟フォルダ直下 → 明示追加。**再帰しない** — 深く掘ると
    無関係なプロジェクトのパラメータが候補に混ざり、比較表が意味を失う。
    """
    root = Path(dataset_root).expanduser()
    pairs: list[tuple[Path, str]] = [(root, "same_dir")]
    try:
        siblings = sorted(p for p in root.parent.iterdir() if p.is_dir())
    except OSError:
        siblings = []
    for sibling in siblings:
        if sibling.resolve() == root.resolve():
            continue
        pairs.append((sibling, "sibling"))
    for extra in (search_dirs or []):
        pairs.append((Path(extra).expanduser(), "given"))
    return pairs


def past_run_method_files(dataset_root: Path) -> list[Path]:
    """過去 run が実際に使ったメソッドファイルを返す。

    `analysis-job.json` の `software.method_file` から引く。**ファイル名で拾わない** —
    `run_dir/effective-method.txt` は「原本からの書き換え（LBM の補完・相対パスの
    絶対化）が要るとき」だけ書かれるので、グロブでは取りこぼす。
    """
    runs = Path(dataset_root).expanduser() / RUNS_SUBDIR_NAME
    out: list[Path] = []
    try:
        entries = sorted(runs.iterdir())
    except OSError:
        return out
    for run_dir in entries:
        job = run_dir / "analysis-job.json"
        try:
            record = json.loads(job.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        # 有効な JSON でも object とは限らない（`[]` / `"x"` / `3` 等）。そのまま
        # `.get()` すると AttributeError が MCP 境界まで漏れる。
        if not isinstance(record, dict):
            continue
        software = record.get("software")
        if not isinstance(software, dict):
            software = {}
        declared = software.get("method_file")
        if not isinstance(declared, str):
            declared = ""
        declared = declared.strip()
        if not declared:
            continue
        path = Path(declared)
        if path.is_file():
            out.append(path)
    return out


def write_effective_method_file(src: Path, dest: Path, overrides: dict[str, str]) -> Path:
    """`src` を写して `overrides` のキーを差し替えたメソッドファイルを `dest` に書く。

    **元ファイルは触らない**。ユーザーのデータフォルダにある GUI 由来の
    パラメータを書き換えると、次に GUI で開いたときの整合が取れなくなる。
    出力は ASCII / LF（`ConfigParser` は `StreamReader(path, Encoding.ASCII)`）。

    同じキーの行が複数あれば**全部**差し替える。`ConfigParser` は全行を順に
    読んで後の行が勝つので、最初の 1 行だけでは後ろの重複行が上書きを打ち消す。
    """
    by_lower = {key.lower(): key for key in overrides}
    applied: set[str] = set()
    lines: list[str] = []
    for line in src.read_text(encoding="ascii", errors="replace").splitlines():
        stripped = line.strip()
        override_key = None
        if stripped and not stripped.startswith("#"):
            positions = [i for i in (stripped.find(":"), stripped.find("=")) if i > 0]
            if positions:
                override_key = by_lower.get(stripped[:min(positions)].strip().lower())
        if override_key is None:
            lines.append(line)
        else:
            lines.append(f"{override_key}: {overrides[override_key]}")
            applied.add(override_key)

    for key, value in overrides.items():
        if key not in applied:
            lines.append(f"{key}: {value}")

    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    return dest


def discover_method_candidates(
    dataset_root,
    *,
    polarity: str | None = None,
    omics: str | None = "lipidomics",
    search_dirs=None,
) -> tuple[list[MethodCandidate], list[str]]:
    """候補を集めて `usable` と `key_params` を付ける。探した場所も返す。

    極性で**落とさない**。別極性は `needs_polarity_conversion` として提示し、
    console_method_template を経由させる（黙って別解析にしないため）。
    """
    root = Path(dataset_root).expanduser()
    by_path: dict[Path, MethodCandidate] = {}
    searched: list[str] = []

    for directory, origin in method_search_dirs(root, search_dirs):
        searched.append(str(directory))
        for candidate in scan_dir_for_method_files(directory, origin):
            by_path.setdefault(Path(candidate.path).resolve(), candidate)

    runs_dir = root / RUNS_SUBDIR_NAME
    if runs_dir.is_dir():
        searched.append(str(runs_dir))
    for path in past_run_method_files(root):
        keys = read_method_keys(path)
        if not keys:
            continue
        # 出所の情報量が多い past_run を優先して上書きする（同じパスが
        # same_dir としても拾われうる）。
        by_path[path.resolve()] = MethodCandidate(
            path=str(path),
            ion_mode=(keys.get("ion mode") or "").strip().lower() or None,
            omics=(keys.get("target omics") or "").strip().lower() or None,
            has_lbm=bool((keys.get(LBM_KEY.lower()) or "").strip()),
            mtime=path.stat().st_mtime,
            origin="past_run",
        )

    selected = [c for c in by_path.values()
                if omics is None or c.omics == omics]
    annotated: list[MethodCandidate] = []
    attach_params = len(selected) <= KEY_PARAMS_MAX_CANDIDATES
    for candidate in selected:
        usable = None
        if polarity is not None:
            usable = ("direct" if candidate.ion_mode == polarity
                      else "needs_polarity_conversion")
        key_params = (extract_key_params(read_method_keys(Path(candidate.path)))
                      if attach_params else None)
        annotated.append(replace(candidate, usable=usable, key_params=key_params))

    annotated.sort(key=lambda c: (c.usable != "direct", -c.mtime))
    return annotated, searched


# ---------- pipeline専用: 既知の参照キー登録（spec §4.3） ----------

#: pipeline（Task13以降）が原本（method_file）基準で絶対解決する既知の参照キー。
#: ここに無いキーの値はパスと決め付けて解決・コピーしない
#: （spec §4.3「未知キーの値をパスと決め付けてコピーしない」）。
#: LBM自体は既存の resolve_lbm（build_tree/env/exe_dir へのフォールバックを持つ）
#: が別途解決するが、「宣言されているのに解決できない」場合の検出はここが担う。
#:
#: **既存lipidomics v1経路（`pipeline/inputs.py`の`inspect_inputs`）が使う既定値は
#: 変えない**（controller裁定 fix round 2 finding 3, A01「既存lipidomics実行が
#: 変わらない」）。metabolomics向けの3キーを黙って混ぜると、GUIが書き出した
#: lipidomicsのメソッドファイルが偶然（無関係の理由で）`Msp file path`等に
#: 古い/無効な値を持っているだけで、これまで無視されていたものが突然
#: `METHOD_REFERENCE_UNRESOLVED`で落ちるようになる——既存lipidomics実行への
#: 無言の副作用であり、避ける。metabolomics側で広い集合を使いたい呼び出し元は
#: `METABOLOMICS_REFERENCE_KEYS`を明示的に渡す（`lipidmix/console/profiles.py`の
#: `resolve_profile_inputs`参照）。
REFERENCE_KEYS: frozenset[str] = frozenset({LBM_KEY.lower()})

#: metabolomicsプロファイル（spec §5.1）が使う既知参照キーの広い集合。
#: `resolve_method_references`/`method_reference_fingerprint`へ明示的に
#: `reference_keys=METABOLOMICS_REFERENCE_KEYS`として渡したときだけ有効になる
#: ——`REFERENCE_KEYS`（既定値・lipidomics v1が暗黙に使う）とは別の集合に
#: しておくことで、v1呼び出し側のコードを一切変更せずに済む。
METABOLOMICS_REFERENCE_KEYS: frozenset[str] = frozenset({
    LBM_KEY.lower(), MSP_KEY.lower(), TEXT_DB_KEY.lower(), RT_REFERENCE_KEY.lower(),
})


def method_reference_fingerprint(
    method_keys: dict[str, str], method_file: Path, *,
    reference_keys: frozenset[str] | None = None,
) -> dict[str, dict]:
    """既知の参照キーごとに宣言値の解決結果を返す（副作用なし・例外を投げない弱い版）。

    select_method の候補グルーピングに使う。値が空、またはキー自体が
    `reference_keys`（省略時は`REFERENCE_KEYS`——既存lipidomics v1呼び出し側の
    挙動を変えない既定値）に無ければそのキーは結果に含めない。解決できた場合は
    参照先の内容ハッシュ（`sha256`）を持ち、できなければ `resolved=False` と
    宣言値だけを持つ——「同じ相対パス文字列でも解決元ディレクトリが違えば
    実効参照が異なりうる」ことを、この関数の呼び出し元（`method_file`引数に
    候補ごとの実ファイルパスを渡す）が自然に表現する。
    """
    keys = REFERENCE_KEYS if reference_keys is None else reference_keys
    out: dict[str, dict] = {}
    for key in keys:
        declared = (method_keys.get(key) or "").strip()
        if not declared:
            continue
        candidate = Path(declared)
        if not candidate.is_absolute():
            candidate = method_file.parent / candidate
        if candidate.is_file():
            digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
            out[key] = {"declared": declared, "resolved_path": str(candidate.resolve()),
                        "sha256": digest, "resolved": True}
        else:
            out[key] = {"declared": declared, "resolved_path": str(candidate),
                        "sha256": None, "resolved": False}
    return out


def resolve_method_references(
    method_keys: dict[str, str], method_file: Path, *,
    reference_keys: frozenset[str] | None = None,
) -> dict[str, Path]:
    """既知の参照キーを原本(method_file)基準で絶対解決する（pipeline専用・厳格版）。

    未対応（宣言されているのに解決できない）の参照キーがあれば
    `DomainError("METHOD_REFERENCE_UNRESOLVED", ...)` を送出し、元のまま実行しない
    （spec §4.3）。呼び出し側は「最終的に採用したメソッドファイル」に対して
    これを呼ぶ想定——候補選別段階では弱い版（`method_reference_fingerprint`）を使う。
    `reference_keys`は`method_reference_fingerprint`と同じ意味（省略時は
    `REFERENCE_KEYS`＝lipidomics v1の既定挙動）。
    """
    fingerprint = method_reference_fingerprint(method_keys, method_file,
                                               reference_keys=reference_keys)
    resolved: dict[str, Path] = {}
    for key, info in fingerprint.items():
        if not info["resolved"]:
            raise DomainError(
                "METHOD_REFERENCE_UNRESOLVED",
                f"メソッドの参照キー {key!r} が指すファイルが見つかりません: "
                f"{info['declared']}（{method_file} 基準で解決: {info['resolved_path']}）",
                {"key": key, "declared": info["declared"],
                 "resolved_path": info["resolved_path"]},
            )
        resolved[key] = Path(info["resolved_path"])
    return resolved
