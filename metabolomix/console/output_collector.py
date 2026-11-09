"""ジョブ生成物の収集と役割付け。

実行前後のディレクトリスナップショット差分で今回ジョブの生成物を特定する。
更新時刻ではなく「実行前に存在しなかった or サイズ変化したファイル」を使う。

mzTab エントリの polarity / measure は 2 つの独立な証拠から決める:
  1. ファイル名の信号（MS-DIAL GUI の `Height_` / `Area_` prefix、`Neg` / `Pos` トークン）
  2. ジョブが宣言した値（pipeline の要求でユーザーが指定したもの）
ファイル名は**実物の性質**を語り、宣言は**意図**でしかないので、両方あって食い違う
ときはファイル名を採り、食い違い自体を validation に記録する。ファイル名が黙って
いるときだけ宣言値で埋める。どちらも無いときにだけ既定へ落とす。

**信号なしを既定値と混ぜてはいけない**——旧実装は極性トークンを持たない
`Height_AlignmentResult_<timestamp>.mzTab`（MS-DIAL のアライメント出力名の実物）を
無条件に positive としていたため、negative で計画したジョブの全エントリが positive
と記録され、analysis-job.json と状態の表示が嘘を述べ続けた。
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
# （worker.json は実在するので残す）。worker.log はワーカー自身の出力、
# job.lock は job 単位の排他ロックの実体ファイル。pipeline-owner.json は
# 旧版の pipeline が書いた所有権 sidecar（既存の run に残っていても収集しない）。
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
    （要求でユーザーが宣言したもの）。ファイル名が黙っている軸を埋め、
    食い違う軸を記録するために使う。省略時は従来どおり推定と既定値だけで決める。

    Returns
    -------
    (mztab_entries, other_artifacts):
        mztab_entries: 正準候補になれる .mzTab を MztabEntry として返す。
        other_artifacts: それ以外。未対応 measure の .mzTab（Normalized*）も
            role=unsupported_mztab としてここに入る——記録は残すが、
            load_dataset が正準として選べないようにするため。
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
# pipeline の worker も同じ `execution.supervise` を通るので、実行方式によって所見が
# 変わることはない。ワーカーは MCP ツール層を import しないため、この層に置く必要がある。

def artifact_abs_path(job, root: str, rel: str) -> Path:
    """analysis-job.v2 の root に応じて生成物の絶対パスを解決する。"""
    base = Path(job.dataset_root) if root == "dataset_root" else Path(job.run_dir)
    return (base / rel).resolve()
