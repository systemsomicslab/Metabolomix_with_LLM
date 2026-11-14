"""特定の差次的結果から契約 TSV を作る（spec §8）。

列定義の正準は `metabolomix/analysis/export_contract.py`（別リポジトリ massbank-context
との契約でもある）。ここが決めるのは**どの結果を、どの来歴で書き出すか**だけで、
列の追加・改名・並べ替えはしない。

書き出す前に、その結果が今の前処理から出たものかを確かめる（`assert_current`）。前処理を
やり直した後の古い結果を書き出すと、TSV の数字と**現在の**前処理条件が並んだファイルができる。
どちらも正しく見えるので、読み手にはずれが分からない。

探索専用のデータセット（中断された実行の出力・検証を通っていない pipeline の単位）は断らずに書き出し、
メタ行 `# exploratory = true` と理由を足す（欠けた検体が「その群には無い」ように見えうることを
ファイルの中に残す）。

メタ行の前処理条件は、その結果が使った行列の登録（`ds.stat_matrices` の前処理の設定）から取る。
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

from metabolomix.analysis import export_contract
from metabolomix.analysis.result_state import assert_current
from metabolomix.analysis.stats_results import as_two_group, get_matrix
from metabolomix.core.atomic_io import DomainError

__all__ = ["exploratory_reason", "export_dataset_result"]


def exploratory_reason(ds) -> str | None:
    """探索専用のデータセットなら理由（1 行。`=` とタブを含めない）、そうでなければ None。"""
    if not getattr(ds, "exploratory_only", False):
        return None
    reason = getattr(ds, "exploratory_reason", None) or "incomplete analysis run (exploratory only)"
    return " ".join(str(reason).replace("=", ":").split())


def export_dataset_result(ds, result: dict, path: Path, curation: dict | None = None,
                          contrast_id: str | None = None,
                          exclude_feature_ids=None) -> dict:
    """差次的結果を InChIKey 付きの契約 TSV として書き出す。

    `result` は統計エンジンの `compare` の記録（`as_two_group` で 1 つの対比を読む。対比が複数なら
    `contrast_id` が要り、無ければ `CONTRAST_REQUIRED`）。2×2 の交互作用の対比も書ける（`mean_a` /
    `mean_b` は空欄、メタ行に定義）。エンジンのメタ行を `msi_note` の後ろへ足す
    （`export_contract.stats_meta_lines`）。15 列・列の順序・契約の版は変えない。

    `exclude_feature_ids` は検定の**後で**書き出しから落とす特徴 ID（外来の分子など）。q は計算し直さない
    （BH の母集団は検定のときのまま）ことと落とした件数をメタ行に書く。結果に無い ID は戻り値の
    `exclude_feature_ids_not_found` に挙げる（誤りにはしない）。

    Returns
    -------
    dict
        `output_path` ・件数・contract_version などの要約。
    """
    assert_current(ds, result)
    record = result
    # 比較でない結果・対比の曖昧さはここで DomainError。書き出しは交互作用の対比も読む。
    result = as_two_group(record, contrast_id, allow_interaction=True)

    q_threshold = result["q_threshold"]
    log2fc_threshold = result["log2fc_threshold"]

    dropped_ids = {str(f) for f in (exclude_feature_ids or [])}
    result_ids = {str(row["feature"]) for row in result["results"]}
    not_found = sorted(dropped_ids - result_ids)
    n_dropped = len(dropped_ids & result_ids)

    rows: list[dict] = []
    n_unannotated = 0
    annotations = getattr(ds, "feature_annotations", {}) or {}
    for row in result["results"]:
        fid = row["feature"]
        if str(fid) in dropped_ids:
            continue                            # 検定の後で落とす（q は計算し直さない）
        meta = ds.feature_metadata.get(fid, {})
        # 行単位で出所を1つに決める。決め手は InChIKey を供給した側——列ごとに
        # 選ぶと、name が SML 由来で inchikey が SME 由来、という食い違った行ができる。
        override = None
        if curation is not None and curation["state"] == "applied":
            from metabolomix.curation.apply import identity_for
            override = identity_for(_as_spot_id(fid), curation["flag_set"])
        inchikey = str(meta.get("inchikey") or "").strip()
        source = meta
        name_source = "mztab_sme"               # 名前の実体は SME 行（旧 "mztab_smf" は誤り）
        if override is not None:
            # curation の assign が同定を置き換える。SME / SML の分岐は通らない。
            inchikey = override["inchikey"].strip()
            source = {"name": override["name"], "inchikey_source": "curation"}
            name_source = "curation"
        elif not inchikey:
            annotation = annotations.get(fid) or {}
            if annotation.get("ambiguous"):
                annotation = {}
            candidate_key = str(annotation.get("inchikey") or "").strip()
            if candidate_key:
                inchikey = candidate_key
                source = annotation
                name_source = "mztab_sml"       # MS1 照合のみ。MS/MS 裏付けではない
        if not inchikey:
            n_unannotated += 1
            continue
        rows.append({
            "spot_id": fid,                     # mzTab-M の SMF_ID（メタ行で id_space を宣言）
            "name": (source.get("name") or "").strip(),
            "name_source": name_source,
            "ontology": override["ontology"] if override is not None else "",  # mzTab-M に対応物なし
            "inchikey": inchikey,
            "inchikey_source": source.get("inchikey_source") or "",
            "msi_level": None,                  # 同上（.arf2 由来の注釈確度が無い）
            "mz": meta.get("mz"),
            "rt": meta.get("rt"),
            "log2fc": row.get("log2fc"),
            "p_value": row.get("p"),
            "q_value": row.get("q"),
            "mean_a": row.get("mean_a"),
            "mean_b": row.get("mean_b"),
            # as_two_group が付けた印（判定が changed のときだけ真。R16）
            "significant": export_contract.row_significant(
                row, q_threshold=q_threshold, log2fc_threshold=log2fc_threshold),
        })

    curation_line = None
    stats = None
    if curation is not None:
        from metabolomix.curation.apply import meta_line
        flag_set = curation["flag_set"]
        if curation["state"] == "applied":
            wrong = {str(s) for s in flag_set["wrong"]}
            redundant = {str(s) for s in flag_set.get("redundant", ())}
            suspect = {str(s) for s in flag_set["suspect"]}
            before = len(rows)
            rows_kept = [r for r in rows if str(r["spot_id"]) not in wrong]
            n_wrong = before - len(rows_kept)
            rows = [r for r in rows_kept if str(r["spot_id"]) not in redundant]
            n_redundant = len(rows_kept) - len(rows)
            n_unannotated += n_wrong + n_redundant
            stats = {"wrong_excluded": n_wrong, "redundant_excluded": n_redundant,
                     "suspect": sum(1 for r in rows if str(r["spot_id"]) in suspect),
                     "assigned": sum(1 for r in rows if r.get("name_source") == "curation")}
        curation_line = meta_line(curation["state"], flag_set, stats)

    n_total = len(result["results"])
    if not rows:
        raise DomainError(
            "NO_ANNOTATED_FEATURES",
            "InChIKey が付いた特徴が 0 件のため書き出しません。"
            "下流のパスウェイ解析に使える背景集合がありません。",
            {"n_features_total": n_total, "n_with_inchikey": 0,
             "n_unannotated": n_unannotated, "n_excluded_after_test": n_dropped})

    reason = exploratory_reason(ds)
    extra = [*_exploratory_lines(reason),
             *_post_test_exclusion_lines(n_dropped, exclude_feature_ids is not None)]
    lines = [*_meta_lines(ds, result, rows, n_total, n_unannotated,
                          curation_line=curation_line,
                          preprocess=_engine_preprocess(ds, record), extra_lines=extra),
             "\t".join(export_contract.EXPORT_COLUMNS)]
    lines += [export_contract.format_row(r) for r in rows]
    out = _atomic_write_text(Path(path), "\n".join(lines) + "\n")

    info = {
        "output_path": str(out),
        "contract_version": export_contract.CONTRACT_VERSION,
        "result_id": result["provenance"]["result_id"],
        "contrast_id": result["contrast_id"],
        "contrast_kind": result["contrast_kind"],
        "group_a": result["a"],
        "group_b": result["b"],
        "n_features_total": n_total,
        "n_with_inchikey": len(rows),
        "n_unannotated": n_unannotated,
        "p_adjust_method": result["statistics"]["p_adjust_method"],
        "curation": curation["state"] if curation else None,
        "curation_stats": stats if curation else None,
    }
    if exclude_feature_ids is not None:
        info["n_excluded_after_test"] = n_dropped
        info["exclude_feature_ids_not_found"] = not_found
    if reason is not None:
        info["exploratory"] = {"exploratory_only": True, "reason": reason}
    return info


# ---------- 内部 ----------

def _as_spot_id(fid):
    """mzTab-M の SMF_ID（文字列）を、curation フラグの spot_id（int の MasterAlignmentID）と
    突き合わせられる形にする。数字でなければそのまま返す。"""
    try:
        return int(fid)
    except (TypeError, ValueError):
        return fid


def _unadjusted_confounded_lines(prov: dict) -> list[str]:
    """spec §7.4(R16): allow_confounded=true で継続した比較は、TSV付随メタにも
    その旨を残す。図（ASCII 但し書き）・レポートと合わせて3つの経路の1つ。

    真偽値の行（機械可読）と、`as_two_group` が来歴へ足した日本語の説明
    （`result["provenance"]["warnings"]`）の両方を書く。未調整でなければ何も
    足さない——常に行があると、本当に注意が要る比較が埋もれる。
    """
    comparison = prov.get("comparison") or {}
    if not comparison.get("unadjusted_confounded"):
        return []
    explanation = next((w for w in prov.get("warnings") or [] if "未調整" in w), None)
    lines = ["# unadjusted_confounded = true"]
    if explanation:
        lines.append(f"# unadjusted_confounded_note = {explanation}")
    return lines


def _exploratory_lines(reason: str | None) -> list[str]:
    """探索専用のデータセットの印（機械可読の真偽と理由）。探索専用でなければ何も足さない。"""
    if reason is None:
        return []
    return [f"# exploratory = true\texploratory_reason = {reason}"]


def _post_test_exclusion_lines(n_dropped: int, requested: bool) -> list[str]:
    """検定の後で落とした特徴の件数と、q を計算し直していないことの宣言（指定が無ければ何も足さない）。"""
    if not requested:
        return []
    return [f"# n_excluded_after_test = {n_dropped}"
            "\texcluded_after_test_note = rows dropped after testing; q_value not recomputed"
            " (the BH population still includes them)"]


def _engine_preprocess(ds, record: dict) -> dict | str:
    """エンジンの結果の前処理条件（行列の登録にある前処理の設定）。無ければ不明と明示する。"""
    entry = get_matrix(ds, record["matrix_id"])
    return dict(entry["config"]) if entry else f"unknown（行列 {record['matrix_id']} が登録に無い）"


def _meta_lines(ds, result, rows, n_total, n_unannotated,
                *, curation_line, preprocess, extra_lines=()) -> list[str]:
    prov = result["provenance"]
    return export_contract.build_meta(
        group_a=result["a"], n_a=result["n_a"],
        group_b=result["b"], n_b=result["n_b"],
        q_threshold=result["q_threshold"],
        log2fc_threshold=result["log2fc_threshold"],
        log_transform=result.get("log_transform"),
        n_features_total=n_total, n_with_inchikey=len(rows),
        n_unannotated=n_unannotated,
        # mzTab-M に .arf2 由来の注釈確度が無いことを、空欄の意味とあわせて宣言する。
        msi_note=("# ontology / msi_level は mzTab-M に対応物が無いため空欄"
                  "（『該当なし』ではなく『この経路では取得していない』）"),
        source_lines=[
            f"# source_mztab = {'; '.join(ds.source_files) or ''}",
            f"# source_job = {ds.job_path or ''}",
            "# id_space = mztab_smf_id",
            f"# result_id = {prov['result_id']}",
            f"# preprocess_id = {'; '.join(prov.get('parent_ids') or [])}",
            f"# source_verification = {getattr(ds, 'source_verification', '')}",
            *_unadjusted_confounded_lines(prov),
            *([curation_line] if curation_line else []),
        ],
        preprocess_line=f"# preprocess = {preprocess}",
        # 契約の版を上げずに足す行は msi_note の後ろだけ（探索専用・検定後の除外はさらにその末尾）
        stats_lines=[*export_contract.stats_meta_lines(result), *extra_lines],
    )


def _atomic_write_text(path: Path, text: str) -> Path:
    """同じ親の一時ファイルへ書いてから置換する。

    途中で落ちた出力を「新しい完成品」として残さない。前の TSV が下流で使われて
    いる最中に、半端な内容へ差し替わるのも防ぐ。
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp",
                                    dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return path
