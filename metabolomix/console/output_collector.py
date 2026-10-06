"""ジョブ生成物の収集と役割付け。

実行前後のディレクトリスナップショット差分で今回ジョブの生成物を特定する。
更新時刻ではなく「実行前に存在しなかった or サイズ変化したファイル」を使う。

mzTab エントリの polarity / measure は 2 つの独立な証拠から決める:
  1. ファイル名の信号（MS-DIAL GUI の `Height_` / `Area_` prefix、`Neg` / `Pos` トークン）
  2. ジョブが宣言した値（console_plan でユーザーが指定したもの）
ファイル名は**実物の性質**を語り、宣言は**意図**でしかないので、両方あって食い違う
ときはファイル名を採り、食い違い自体を validation に記録する。ファイル名が黙って
いるときだけ宣言値で埋める。どちらも無いときにだけ既定へ落とす。

**信号なしを既定値と混ぜてはいけない**——旧実装は極性トークンを持たない
`Height_AlignmentResult_<timestamp>.mzTab`（MS-DIAL のアライメント出力名の実物）を
無条件に positive としていたため、negative で計画したジョブの全エントリが positive
と記録され、analysis-job.json と console_status が嘘を表示し続けた。
"""
from __future__ import annotations

import os
import re
from pathlib import Path

from metabolomix.handoff.schema import Artifact, MztabEntry, sha256_file

# ジョブ運用のためにランディレクトリへ書かれるファイル。MS-DIAL の生成物ではない。
# msdial.log は supervise が必ず作るため、除外しないと「出力ゼロ」を検出できない。
# analysis-job.json は実行中に status 遷移で書き換わるため、差分に混入する。
# execution-result.json / worker.json は監視ワーカーが書く終了証跡・監視状態
# （`metabolomix.console.execution` の RECEIPT_FILENAME / SUPERVISION_FILENAME）。
# control.json はここに**置かない**——実際に書く側がどこにも居ない予約名は、
# 将来 MS-DIAL がその名前で出力を書いたときに黙って捨てる罠になる
# （worker.json は実在するので残す）。worker.log は切り離しワーカー自身の出力、
# job.lock は job 単位の排他ロックの実体ファイル。pipeline-owner.json は
# Task 14（pipeline/store.py）が Console 起動前に書く所有権 sidecar
# （このジョブをどの pipeline が所有しているかの参照だけを持つ）。
# いずれも MS-DIAL の生成物ではないので、除外しないと「実行の結果できたファイル」
# として誤収集され、生成物ゼロの実行が partial に化ける。
_OPERATIONAL_FILES = frozenset({
    "msdial.log", "analysis-job.json",
    "execution-result.json", "worker.json",
    "worker.log", "job.lock", "pipeline-owner.json",
})

# dataset_root を撮るときに降りないディレクトリ名。job_manager の RUNS_SUBDIR と
# 同じ値だが、収集層はジョブ管理を知らないままにしておくため定数を複製する。
RUNS_SUBDIR = "runs"

# 拡張子パターン → (role, format) のマッピング（長い拡張子を先に評価する）。
# Console の実出力とエクスポートを拡張子だけで判定する。どのルートかは
# Artifact.root が持つため、ここでは拡張子だけを見る。
_ROLE_MAP: list[tuple[str, str, str]] = [
    ("_tags.xml",  "peak_tags",          "tagsxml"),
    (".qa.tsv",    "quality_matrix",     "qatsv"),
    (".msp2.dbs",  "library_cache",      "msp2dbs"),
    (".EIC.aef",   "chromatogram",       "eicaef"),
    (".mzTab",     "primary_mztab",      "mztab"),
    (".mdproject", "gui_project",        "mdproject"),
    (".mdalign",   "alignment_table",    "mdalign"),
    (".mdpeak",    "sample_peak_table",  "mdpeak"),
    (".mdmsp",     "msms_spectra",       "mdmsp"),
    (".mddata",    "project_data",       "mddata"),
    (".msp2",      "library_snapshot",   "msp2"),
    (".arf2",      "spot_catalog",       "arf2"),
    (".arf",       "peak_matrix_source", "arf"),
    (".pai2",      "sample_peaks",       "pai2"),
    (".dcl",       "msms_evidence",      "dcl"),
]

# ファイル名の極性トークン。`Neg_` / `_NEG` / `.pos.` のような**語**として現れた
# ものだけを信号とみなす。部分文字列一致にすると Negev のような無関係な語を
# 極性と誤読する。
_POLARITY_TOKEN_RE = re.compile(r"(?:^|[^a-z])(neg|pos)(?:[^a-z]|$)", re.IGNORECASE)
# 定量種別は MS-DIAL GUI の prefix でのみ判断する（validator._HEIGHT_PREFIX_RE と同規則）。
_HEIGHT_PREFIX_RE = re.compile(r"^Height_", re.IGNORECASE)
_AREA_PREFIX_RE = re.compile(r"^Area_", re.IGNORECASE)
# spec §8.1: normalized value は「実装と実データ検証を追加するまで未対応」。
# peak_height / peak_area_above_zero のどちらでもないので、正準候補から外す。
_NORMALIZED_PREFIX_RE = re.compile(r"^Normalized", re.IGNORECASE)

UNSUPPORTED_MZTAB_ROLE = "unsupported_mztab"

_DEFAULT_POLARITY = "positive"
_DEFAULT_MEASURE = "peak_height"


def snapshot(
    directory: Path,
    exclude_dir_names: frozenset[str] | set[str] = frozenset(),
) -> dict[str, int]:
    """ディレクトリ以下の全ファイルを {相対パス文字列: サイズ} で返す。

    exclude_dir_names に指定した名前のディレクトリは走査しない。
    """
    result: dict[str, int] = {}
    if not directory.exists():
        return result
    for root, dirs, files in os.walk(directory):
        dirs[:] = [d for d in dirs if d not in exclude_dir_names]
        for name in files:
            fp = Path(root) / name
            try:
                result[str(fp.relative_to(directory))] = fp.stat().st_size
            except (OSError, ValueError):
                pass
    return result


def collect_artifacts(
    roots: dict[str, Path],
    befores: dict[str, dict[str, int]],
    *,
    declared_polarity: str | None = None,
    declared_measure: str | None = None,
) -> tuple[list[MztabEntry], list[Artifact]]:
    """複数ルートを before スナップショットと比較し、新規・変化ファイルを収集する。

    declared_polarity / declared_measure は analysis-job.json の project 値
    （console_plan でユーザーが宣言したもの）。ファイル名が黙っている軸を埋め、
    食い違う軸を記録するために使う。省略時は従来どおり推定と既定値だけで決める。

    Returns
    -------
    (mztab_entries, other_artifacts):
        mztab_entries: 正準候補になれる .mzTab を MztabEntry として返す。
        other_artifacts: それ以外。未対応 measure の .mzTab（Normalized*）も
            role=unsupported_mztab としてここに入る——記録は残すが、
            dataset_load が正準として選べないようにするため。
    """
    mztab_entries: list[MztabEntry] = []
    other_artifacts: list[Artifact] = []

    for root_name, root_path in roots.items():
        before = befores.get(root_name, {})
        exclude = {RUNS_SUBDIR} if root_name == "dataset_root" else frozenset()
        after = snapshot(root_path, exclude_dir_names=exclude)
        new_or_changed = {
            rel: size for rel, size in after.items()
            if (rel not in before or before[rel] != size)
            and Path(rel).name not in _OPERATIONAL_FILES
        }
        for rel_str in sorted(new_or_changed):
            fp = root_path / rel_str
            if not fp.is_file():
                continue
            role, fmt = _assign_role(rel_str)
            # role の付かないファイルはハッシュしない。大きな副産物もある。
            checksum = sha256_file(fp) if role != "unknown" else ""

            if fmt == "mztab":
                if _NORMALIZED_PREFIX_RE.match(fp.name):
                    other_artifacts.append(Artifact(
                        path=rel_str, role=UNSUPPORTED_MZTAB_ROLE, format=fmt,
                        sha256=checksum, root=root_name,
                    ))
                    continue
                polarity, measure, validation = _resolve_mztab_meta(
                    fp.name, declared_polarity, declared_measure)
                mztab_entries.append(MztabEntry(
                    path=rel_str, polarity=polarity, measure=measure,
                    sha256=checksum, validation=validation, root=root_name,
                ))
            else:
                other_artifacts.append(Artifact(
                    path=rel_str, role=role, format=fmt,
                    sha256=checksum, root=root_name,
                ))

    return mztab_entries, other_artifacts


def is_upstream_artifact(name: str) -> bool:
    """`name` が MS-DIAL の生成物として役割を割り当てられるかを返す。

    `_ROLE_MAP` が唯一の出所。MS-DIAL は `-o` だけでなく `-i` 側にも
    `.arf` / `.pai2` / `.dcl` / `_tags.xml` 等を書くので、「入力フォルダに
    入力以外がある」を異常と決め付ける側（`pipeline.inputs.stage_inputs`）が
    この判定を借りる。同じ拡張子表を2か所に書くとドリフトする。

    役割が付かないもの（生データ・よそのバッチのコピー等）は False。
    """
    return _assign_role(name)[0] != "unknown"


def _assign_role(rel_str: str) -> tuple[str, str]:
    lower = rel_str.lower()
    for suffix, role, fmt in _ROLE_MAP:
        if lower.endswith(suffix.lower()):
            return role, fmt
    return "unknown", Path(rel_str).suffix.lstrip(".") or "bin"


def _infer_mztab_meta(filename: str) -> tuple[str | None, str | None]:
    """ファイル名から極性と定量種別を推定する。**推定できなければ None を返す。**

    None は「このファイル名は当該軸について何も語っていない」の意味。既定値へ
    落とすのは呼び出し側の判断であり、ここで既定値を返すと「positive と書いてある」
    と「何も書いていない」が区別できなくなる。

    MS-DIAL GUI の命名規則: `Height_AlignmentResult_...mzTab` / `Area_AlignmentResult_...mzTab`。
    Console のアライメント出力名には極性トークンが無いことが実データで確認済み。
    """
    m = _POLARITY_TOKEN_RE.search(filename)
    polarity: str | None = None
    if m:
        polarity = "negative" if m.group(1).lower() == "neg" else "positive"

    measure: str | None = None
    if _NORMALIZED_PREFIX_RE.match(filename):
        measure = None      # spec §8.1 未対応。height/area のどちらでもない
    elif _HEIGHT_PREFIX_RE.match(filename):
        measure = "peak_height"
    elif _AREA_PREFIX_RE.match(filename):
        measure = "peak_area_above_zero"

    return polarity, measure


def _resolve_mztab_meta(
    filename: str,
    declared_polarity: str | None,
    declared_measure: str | None,
) -> tuple[str, str, dict]:
    """ファイル名の推定とジョブの宣言値から (polarity, measure, validation) を決める。"""
    inferred_polarity, inferred_measure = _infer_mztab_meta(filename)

    polarity, polarity_source = _pick(inferred_polarity, declared_polarity, _DEFAULT_POLARITY)
    measure, measure_source = _pick(inferred_measure, declared_measure, _DEFAULT_MEASURE)

    validation: dict = {
        "polarity_source": polarity_source,
        "measure_source": measure_source,
    }
    conflicts: dict = {}
    if inferred_polarity and declared_polarity and inferred_polarity != declared_polarity:
        conflicts["polarity"] = {"filename": inferred_polarity, "job_declared": declared_polarity}
    if inferred_measure and declared_measure and inferred_measure != declared_measure:
        conflicts["measure"] = {"filename": inferred_measure, "job_declared": declared_measure}
    if conflicts:
        validation["conflicts"] = conflicts
    return polarity, measure, validation


def _pick(inferred: str | None, declared: str | None, default: str) -> tuple[str, str]:
    """ファイル名 → ジョブ宣言 → 既定値 の順に採用し、採用元の名前も返す。"""
    if inferred:
        return inferred, "filename"
    if declared:
        return declared, "job_declared"
    return default, "default"


#: mzTab から読む出所情報の上限。MTD は先頭にまとまっているが、SML の
#: アダクト集計は全行を舐める必要があるので、行数ではなくバイト数で抑える。
_PROVENANCE_MAX_BYTES = 8 << 20


def _read_head(path, limit: int = _PROVENANCE_MAX_BYTES) -> str | None:
    try:
        with open(path, "rb") as fh:
            return fh.read(limit).decode("utf-8", errors="replace")
    except OSError:
        return None


def read_software_version(mztab_path) -> str | None:
    """mzTab の `MTD software[1]` から MS-DIAL の版数を読む。

    実物は `[MS, MS:1003082, MS-DIAL, Msdial console 5.5.241113]` の形で、
    4 番目の要素が版数。`analysis-job` の software.version は Console 実行から
    知りようがないので、成果物側の記録を採る。
    """
    text = _read_head(mztab_path, 1 << 16)
    if text is None:
        return None
    for line in text.splitlines():
        if not line.startswith("MTD"):
            continue
        parts = line.split("\t")
        if len(parts) < 3 or not parts[1].strip().startswith("software["):
            continue
        term = parts[2].strip().strip("[]")
        fields = [f.strip() for f in term.split(",")]
        if len(fields) >= 4 and fields[3]:
            return fields[3]
    return None


def read_adduct_polarity(mztab_path) -> dict:
    """SML のアダクト表記から極性を多数決で推定する。

    Console のアライメント出力名には極性トークンが無いため、`polarity_source` は
    `job_declared` にしかならない（仕様）。宣言ミスを検出する手段がそれだと無い
    ので、**中身**から裏取りする。推定なので `polarity_source` は書き換えず、
    別フィールドとして持つ。
    """
    empty = {"adduct_majority": None, "n_positive": 0, "n_negative": 0}
    text = _read_head(mztab_path)
    if text is None:
        return empty

    column = None
    n_pos = n_neg = 0
    for line in text.splitlines():
        if line.startswith("SMH"):
            header = [h.strip() for h in line.split("\t")]
            if "adduct_ions" in header:
                column = header.index("adduct_ions")
            continue
        if column is None or not line.startswith("SML"):
            continue
        parts = line.split("\t")
        if column >= len(parts):
            continue
        value = parts[column].strip()
        if not value or value == "null":
            continue
        # `[M+H]1+` / `[M-H]1-` の末尾が電荷の符号。
        if value.endswith("+"):
            n_pos += 1
        elif value.endswith("-"):
            n_neg += 1

    if not (n_pos or n_neg):
        return empty
    majority = "positive" if n_pos > n_neg else ("negative" if n_neg > n_pos else None)
    return {"adduct_majority": majority, "n_positive": n_pos, "n_negative": n_neg}


# ---------- 実行後の所見（証跡だけでは読めないもの）----------
# ここに置くのは、収集した成果物そのものを読まないと分からないことだけ。
# 同期実行も切り離しワーカーも metabolomix.console.worker.annotate_job 経由で
# 通るので、実行方式によって所見が変わることはない。ワーカーは MCP ツール層を
# import しないため、この層に置く必要がある。

def artifact_abs_path(job, root: str, rel: str) -> Path:
    """analysis-job.v2 の root に応じて生成物の絶対パスを解決する。"""
    base = Path(job.dataset_root) if root == "dataset_root" else Path(job.run_dir)
    return (base / rel).resolve()


def record_mztab_provenance(job, mztab_entries) -> None:
    """成果物の mzTab から、ジョブ側で分からない出所情報を採る。

    - `software.version`: Console 実行からは知りようがない。mzTab の
      `MTD software[1]` に `Msdial console 5.5.241113` が入っている。
    - 極性の裏取り: Console のアライメント出力名には極性トークンが無いので
      `polarity_source` は `job_declared` にしかならない（仕様）。宣言ミスを
      検出する手段が他に無いため、アダクトの多数決で**別フィールドとして**
      裏取りする。`polarity_source` は書き換えない —— 推定を出所として
      記録すると、そちらのほうが嘘になる。
    """
    for entry in mztab_entries:
        path = artifact_abs_path(job, getattr(entry, "root", "run_dir"), entry.path)
        if not job.software_version:
            version = read_software_version(path)
            if version:
                job.software_version = version
        crosscheck = read_adduct_polarity(path)
        majority = crosscheck["adduct_majority"]
        crosscheck["agrees"] = None if majority is None else (majority == entry.polarity)
        entry.validation["polarity_crosscheck"] = crosscheck


def polarity_crosscheck_warnings(mztab_entries) -> list[str]:
    """アダクトから推定した極性が宣言と食い違うエントリを warning にする。"""
    warnings: list[str] = []
    for entry in mztab_entries:
        check = entry.validation.get("polarity_crosscheck") or {}
        if check.get("agrees") is False:
            warnings.append(
                f"{entry.path}: アダクトの多数決は {check['adduct_majority']} ですが、"
                f"ジョブの宣言は {entry.polarity} です"
                f"（陽性 {check['n_positive']} / 陰性 {check['n_negative']} 行）。"
                "Console 出力のファイル名に極性が入らないため宣言をそのまま記録して"
                "いますが、console_plan の polarity かメソッドファイルの Ion mode を"
                "確認してください。")
    return warnings


def meta_conflict_warnings(mztab_entries) -> list[str]:
    """ファイル名と宣言値が食い違ったエントリを warning にする。

    collect_artifacts はファイル名側を採用する（実物の性質を語るのはファイル）。
    採用の事実だけを validation に残して黙っていると、ユーザーは自分が
    console_plan で宣言した値と違うものを解析していることに気付けない。
    """
    warnings: list[str] = []
    for entry in mztab_entries:
        for axis, detail in (entry.validation.get("conflicts") or {}).items():
            warnings.append(
                f"{entry.path}: {axis} がジョブの宣言と食い違います"
                f"（ファイル名={detail['filename']} / 宣言={detail['job_declared']}）。"
                "ファイル名側を採用しました。"
            )
    return warnings


def unsupported_mztab_warnings(artifacts) -> list[str]:
    """未対応 measure の .mzTab（Normalized*）を拾ったことを伝える。

    記録は残すが正準候補にはしない（spec §8.1）。黙って捨てると
    「出力があるのに dataset_load が読めない」という説明不能な状態になる。
    """
    paths = [a.path for a in artifacts if a.role == UNSUPPORTED_MZTAB_ROLE]
    if not paths:
        return []
    return [
        "未対応の定量種別（normalized）の mzTab を検出しました: "
        + ", ".join(paths)
        + "。peak_height / peak_area_above_zero のみ対応するため、"
        "正準候補（primary_mztab_files）には含めていません。"
    ]


def missing_per_sample_output_warnings(artifacts) -> list[str]:
    """サンプル別ファイル（.pai2）が 1 つも出ていないことを伝える。

    MS-DIAL Console は .pai2 を**生データフォルダ側**に書く（-o ではない。
    docs/HISTRY.md 2026-09-03(6) の実走で確認）。collect_artifacts が両ルートを
    見るようになったので、この検査は「本当に出ていない」ときだけ発火する。
    .pai2 が無いと pai2_parser / dcl_find_msms が読むものが無く、MS/MS 根拠の
    経路が丸ごと空になる。アライメント結果だけは出ているので実行は成功扱いの
    まま、「後で MS/MS を辿れない」ことだけ先に知らせる。
    """
    if any(a.format == "pai2" for a in artifacts):
        return []
    return [
        "サンプル別ファイル（.pai2）が 1 つも生成されていません。"
        "MS/MS 根拠（pai2_parser / dcl_find_msms）を辿る経路が使えません。"
    ]
