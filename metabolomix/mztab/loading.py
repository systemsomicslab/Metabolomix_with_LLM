"""mzTab-M の読み込み（session に依存しない）。spec §7。

MCP ツールも pipeline ワーカーもここを通る。**`session_state` も `mcp_core` も
import しない**——グローバルな現在状態を触る関数はワーカーから呼べないし、呼べても
複数の解析が同じスロットを奪い合う。読み込みの結果は戻り値の DatasetState で、
それをどこへ置くかは呼び出し側が決める。

エラーは MCP の封筒ではなく `DomainError` で返す。封筒は MCP 層の表現なので、
数値層とワーカーがそれを組み立てるのは筋が違う（ツール層が code を見て封筒へ畳む）。

読んだデータセットには、**どれだけ信用してよいか**を 2 つの軸で書く:

``source_verification``
    ``verified``            終了証跡が検証を通り、ジョブ・実行・ファイル hash が
                            互いに一致した完了実行の出力。
    ``legacy_unverified``   ジョブ経由だが、証跡が無い/一致しない（旧実行など）。
    ``direct_unverified``   .mzTab を直接読んだ。実行の出所を語れない。
``exploratory_only``
    完了していない実行（partial / failed / running）の出力。探索には使えるが、
    2 群比較と TSV 出力は拒否する——欠けた検体を欠測ではなく「その群には無い」と
    読み違えたまま結論が出てしまう。

`completed` という文字列だけで `verified` にはしない。旧実装はそれを唯一の根拠に
していたが、旧経路の completed は「ファイルが増えた」以上の意味を持たなかった。
"""
from __future__ import annotations

import json
from pathlib import Path

from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import evidence as mztab_evidence
from metabolomix.mztab.dataset_state import build_dataset_state
from metabolomix.mztab.reader import parse_mztab
from metabolomix.mztab.validator import validate_mztab

__all__ = ["artifact_abs_path", "load_dataset_state", "select_primary_entry"]

#: 完了していない実行の出力に付ける警告。読み手が最初に見る位置へ差す。
_INCOMPLETE_WARNING = (
    "このジョブは status={status} です（completed ではありません）。"
    "MS-DIAL Console の実行は最後まで到達しておらず、読み込んだ mzTab-M は"
    "中断時点の生成物です。特徴量・サンプルが欠けている可能性があるため、"
    "console_status の error と msdial.log を確認してください。")


def artifact_abs_path(job, root: str, rel: str) -> Path:
    """生成物の相対パスを、記録された出所ルートから絶対パスへ戻す。

    MS-DIAL Console は run_dir と生データフォルダの両方へ生成物を出すため、
    analysis-job（v2 以降）が記録する root に応じて解決する。
    """
    base = Path(job.dataset_root) if root == "dataset_root" else Path(job.run_dir)
    return (base / rel).resolve()


def select_primary_entry(job):
    """ジョブの宣言（polarity + measure）で正準 mzTab エントリを一意に選ぶ。

    一意に決まらないときは選ばずに停止する（spec §7「暗黙の単数選択と
    newest_modified_time は廃止する」）。どれを読むかは解析結果そのものを
    変えるので、LLM にもファイル名の辞書順にも決めさせない。
    """
    candidates = job.primary_mztab_files
    by_measure = [e for e in candidates if e.measure == job.measure]
    if not by_measure:
        raise DomainError(
            "QUANTIFICATION_CONFLICT",
            f"ジョブは measure={job.measure} を宣言していますが、"
            "その定量種別の mzTab-M が生成物にありません。"
            "別種別のファイルを代わりに読むと、宣言と違う数値で解析することになります。"
            "意図的に別種別を読むなら mztab_path で明示してください。",
            {"declared_measure": job.measure, "candidates": _describe(candidates)})

    matched = [e for e in by_measure if e.polarity == job.polarity]
    if not matched:
        raise DomainError(
            "POLARITY_MISMATCH",
            f"ジョブは polarity={job.polarity} を宣言していますが、"
            f"measure={job.measure} の候補にその極性がありません。"
            "極性が違えば検出される脂質クラスが変わるため、代替で読みません。"
            "意図的に別極性を読むなら mztab_path で明示してください。",
            {"declared_polarity": job.polarity, "declared_measure": job.measure,
             "candidates": _describe(candidates)})

    if len(matched) > 1:
        raise DomainError(
            "AMBIGUOUS_PRIMARY_MZTAB",
            f"polarity={job.polarity} / measure={job.measure} の候補が "
            f"{len(matched)} 件あり、一意に決まりません。"
            "mztab_path でどれを読むか明示してください。",
            {"candidates": _describe(matched)})
    return matched[0]


def load_dataset_state(*, mztab_path: Path | None = None,
                       job_path: Path | None = None,
                       allow_incomplete: bool = False):
    """mzTab-M を読み、DatasetState を返す（どこへも保存しない）。

    `mztab_path` と `job_path` はどちらか一方だけ。`allow_incomplete=True` を
    渡したときだけ、完了していないジョブの出力を読める（その場合 ds は
    `exploratory_only=True` になる）。
    """
    if type(allow_incomplete) is not bool:
        # 真偽値以外を真偽として読むと、"false" という文字列が True になる。
        raise DomainError("DATASET_BAD_REQUEST",
                          f"allow_incomplete は bool です: {allow_incomplete!r}",
                          {"allow_incomplete": repr(allow_incomplete)})
    if mztab_path and job_path:
        raise DomainError(
            "DATASET_BAD_REQUEST",
            "mztab_path と job_path を同時に指定できません。どちらか一方だけを"
            "使ってください（mztab_path: .mzTab への直接パス。job_path:"
            " analysis-job.json のパス。console_run 完了後は job_path を推奨）。")
    if not mztab_path and not job_path:
        raise DomainError(
            "DATASET_BAD_REQUEST",
            "mztab_path または job_path のどちらか一方を指定してください"
            "（mztab_path: .mzTab への直接パス。job_path: analysis-job.json の"
            "パス。console_run 完了後は job_path を推奨）。")

    if job_path:
        return _load_from_job(Path(job_path).expanduser(),
                              allow_incomplete=allow_incomplete)
    return _load_from_mztab(Path(mztab_path).expanduser())


# ---------- 直接読み ----------

def _load_from_mztab(path: Path):
    if not path.is_file():
        raise DomainError("MZTAB_NOT_FOUND",
                          f"mzTab-M ファイルが見つかりません: {path}",
                          {"path": str(path)})
    parse_result = _parse_and_validate(path)
    ds = build_dataset_state(parse_result, path.name, path)
    # 検出状態（gap-fill）は mzTab-M に無い。隣接 `.arf` から補えるかを試す。
    # 取り込めなくても解析は続けられるので、例外ではなく warning で伝える。
    mztab_evidence.attach_to_dataset(ds, path)
    ds.source_verification = "direct_unverified"
    ds.exploratory_only = False
    return ds


# ---------- ジョブ経由 ----------

def _load_from_job(job_path: Path, *, allow_incomplete: bool):
    if not job_path.is_file():
        raise DomainError("MZTAB_NOT_FOUND",
                          f"analysis-job.json が見つかりません: {job_path}",
                          {"path": str(job_path)})
    from metabolomix.handoff.schema import AnalysisJob
    try:
        job = AnalysisJob.load(job_path)
    except (ValueError, KeyError) as exc:
        raise DomainError("MZTAB_NOT_FOUND",
                          f"analysis-job.json の読み込みに失敗しました: {exc}",
                          {"path": str(job_path)}) from exc

    if job.status != "completed" and not allow_incomplete:
        # 「読めなかった」ではなく「読ませない」。中断された実行の生成物を
        # 完了品として下流へ流すと、欠けた検体が欠測ではなく実測 0 に見える。
        raise DomainError(
            "INCOMPLETE_ANALYSIS_JOB",
            f"未完了の解析ジョブです（status={job.status}）。"
            "中断時点の生成物なので、既定では読み込みません。"
            "探索目的で読むなら allow_incomplete=True を明示してください"
            "（2 群比較と差次的エクスポートは拒否されます）。",
            {"job_id": job.job_id, "status": job.status, "error": job.error})

    if not job.primary_mztab_files:
        raise DomainError(
            "MZTAB_NOT_FOUND",
            f"ジョブ {job.job_id} に primary_mztab_files がありません。"
            "console_run を先に実行してください。",
            {"job_id": job.job_id, "status": job.status})

    entry = select_primary_entry(job)
    mztab_abs = artifact_abs_path(job, getattr(entry, "root", "run_dir"), entry.path)
    if not mztab_abs.is_file():
        raise DomainError(
            "MZTAB_NOT_FOUND",
            f"mzTab-M ファイルが見つかりません: {mztab_abs}"
            f"（job_id={job.job_id}, entry.path={entry.path}）",
            {"job_id": job.job_id, "path": str(mztab_abs)})

    parse_result = _parse_and_validate(mztab_abs)
    ds = build_dataset_state(parse_result, mztab_abs.name, mztab_abs)

    ds.exploratory_only = job.status != "completed"
    if ds.exploratory_only:
        ds.validation_result.setdefault("warnings", []).insert(
            0, _INCOMPLETE_WARNING.format(status=job.status))

    ds.job_path = str(job_path)
    for art in job.artifacts:
        abs_p = str(artifact_abs_path(job, getattr(art, "root", "run_dir"), art.path))
        ds.artifact_paths.setdefault(art.role, []).append(abs_p)

    verified, reason = _verify_source(job, entry, mztab_abs, ds)
    ds.source_verification = "verified" if verified else "legacy_unverified"
    if not verified:
        ds.validation_result.setdefault("warnings", []).append(
            f"この出力は実行証跡で裏取りできていません（{reason}）。"
            "単体の解析には使えますが、パイプラインは検証済みの出力だけを受け付けます。")

    ds.assay_sources = _assay_sources(job, parse_result)

    # artifact_paths を入れ終えてから呼ぶ。handoff が記録した peak_matrix_source を
    # 候補の先頭に使えるのは、この時点以降だけ。
    mztab_evidence.attach_to_dataset(ds, mztab_abs)
    return ds


def _parse_and_validate(path: Path) -> dict:
    parse_result = parse_mztab(path)
    validation = validate_mztab(parse_result)
    if not validation["ok"]:
        raise DomainError(
            "MZTAB_STRUCTURE_INVALID",
            f"mzTab-M 構造不正 ({path.name}): " + "; ".join(validation["errors"]),
            {"errors": validation["errors"], "warnings": validation["warnings"]})
    return parse_result


def _verify_source(job, entry, mztab_abs: Path, ds) -> tuple[bool, str]:
    """終了証跡とファイル hash で「この出力の出所」を裏取りする。

    `completed` という文字列は根拠にしない。旧経路の completed は「実行後に
    ファイルが増えた」以上の意味を持たず、途中で落ちた実行でも付いた。
    """
    from metabolomix.console.execution import receipt_path, validate_execution_record
    from metabolomix.handoff.schema import sha256_file

    if job.status != "completed":
        return False, f"status={job.status}"
    try:
        record = validate_execution_record(
            json.loads(receipt_path(Path(job.run_dir)).read_text(encoding="utf-8")))
    except (OSError, ValueError, DomainError):
        return False, "終了証跡（execution-result.json）が無いか不正"
    if record["job_id"] != job.job_id:
        return False, "終了証跡が別のジョブのもの"
    if record["termination"] != "exited" or record["exit_code"] != 0:
        return False, (f"termination={record['termination']} "
                       f"exit_code={record['exit_code']}")
    recorded = (entry.sha256 or "").strip()
    if not recorded:
        return False, "ジョブに mzTab の hash が記録されていない"
    try:
        if sha256_file(mztab_abs) != recorded:
            return False, "mzTab の内容がジョブ記録の hash と一致しない"
    except OSError:
        return False, "mzTab の hash を計算できない"
    return _verify_effective_method(job)


def _verify_effective_method(job) -> tuple[bool, str]:
    """analysis-job.v3（検証済みプロファイル実行）の実効メソッドを裏取りする。

    spec §6.1 は上流指紋の確認を「計画時、実行直前、**完了検証時**」の 3 点で
    求める。ここがその完了検証時にあたる——実際に Console へ渡ったメソッド
    （依存の絶対パスで書き換えた実効コピー）が消えている・書き換わっている出力を
    `verified` と名乗らせない。原本ではなく実効コピーを見るのは、走ったのが
    そちらだから（spec §5.1「原本hash、書換え後hash、差分を記録する」）。

    profile snapshot を持たないジョブ（v1 / v2）はこの検査の対象外で、判定は
    従来どおり——既存経路の結果は変わらない。

    依存ファイル本体（LBM 等）はここで hash し直さない。数百 MB になりうるものを
    読み込みのたびに読み直すのは重く、その検証はパイプラインの `validate_outputs`
    工程（実行直前・完了検証）の持ち場である。
    """
    from metabolomix.handoff.schema import sha256_file

    snapshot = getattr(job, "profile_snapshot", None)
    if not snapshot:
        return True, ""
    relative = snapshot.get("effective_method_relative_path")
    recorded = (snapshot.get("effective_method_sha256") or "").strip()
    if not relative or not recorded:
        return False, "profile snapshot に実効メソッドの記録が無い"
    effective = (Path(job.run_dir) / relative).resolve()
    try:
        actual = sha256_file(effective)
    except OSError:
        return False, f"実効メソッドを読めない（{relative}）"
    if actual != recorded:
        return False, f"実効メソッドの内容が profile snapshot の hash と一致しない（{relative}）"
    return True, ""


def _assay_sources(job, parse_result: dict) -> dict:
    """assay 列 → 実行時に予定した raw の対応表。

    `map_assays` は Task 3 の完了ゲートと同じ関数を使う（表示名では結合せず、
    `assay[N]-ms_run_ref` → `ms_run[N]-location` を辿る）。同じ対応関係を 2 通りに
    計算すると、完了判定と解析側で別の答えが出る。

    console 層の import を関数内に置くのは、mzTab 層からの依存を module import 時に
    固定しないため（層をまたぐのはこの 1 か所だけにする）。
    """
    from metabolomix.console.execution import read_supervision_state
    from metabolomix.console.validation import map_assays

    inventory = (read_supervision_state(Path(job.run_dir)).get("inputs", {})
                 .get("raw_inventory") or [])
    return map_assays(parse_result, [str(p) for p in inventory])


def _describe(entries) -> list[dict]:
    return [{"path": e.path, "polarity": e.polarity, "measure": e.measure}
            for e in entries]
