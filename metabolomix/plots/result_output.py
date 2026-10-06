"""どの結果から図を描くかを決め、確定した図として保存する（spec §8）。

旧実装は、ARF と mzTab-M の両方に結果があると**黙って ARF を採っていた**。
新しく読み込んだ mzTab の解析をしているつもりでも、前のデータの ARF 結果が図に
なる。戻り値の `source=arf` は「選んだ結果」の報告であって「選んでよいか」の
確認ではない。どちらを使うかは図の数字そのものを変えるので、曖昧なら描かずに
止めて、呼び出し側に指定させる。

保存は一時ファイルへ書いてから置換する。描画の途中で落ちたときに、壊れた PNG が
「保存できた図」として残り、そのままレポートへ貼られるのを防ぐ。
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.tool_helpers import _pca_scatter_arrays, dataset_pca_plot
from metabolomix.plots.volcano import render_volcano_plot

__all__ = ["figure_annotations", "save_result_figure", "select_result"]

#: 図に載せる但し書きの文面。付ける条件は「読み手が数字を誤解しうるとき」だけで、
#: 常に付けると本当に注意が要る図が埋もれる。
#:
#: **ASCII で書く**。matplotlib の既定フォント（DejaVu Sans）は日本語の字形を
#: 持たないので、日本語で書くと豆腐（□□□）が並ぶ。読めない但し書きは無いのと
#: 同じで、しかも「何か警告がある」ことだけ伝わって内容が分からない分たちが悪い。
_EXPLORATORY_NOTE = ("NOTE: exploratory only - drawn from an incomplete analysis "
                     "run; samples/features may be missing.")

#: spec §7.4(R16): allow_confounded=true で継続した比較は、図の中にもその旨が要る。
#: 同じくASCII——`run_comparison`が`result["provenance"]["warnings"]`へ残す日本語の
#: 注記は図には使えない（豆腐になる）ので、ここだけ別文面で書く。
_UNADJUSTED_CONFOUNDED_NOTE = (
    "NOTE: allow_confounded=true - group and batch are fully confounded; "
    "this comparison is UNADJUSTED for batch effect.")


def select_result(candidates: list[dict], *, source: str = "auto",
                  result_id: str | None = None) -> dict:
    """図・出力に使う結果を 1 つに決める。決まらなければ描かずに止める。

    候補は `{source, result_id, dataset_id, valid, result}`。`valid=False` は
    古くなった結果（前処理をやり直した後など）で、指名されても採らない。

    `source="auto"` は「どれでもよい」ではなく「候補が 1 つに定まるなら任せる」。
    2 つ以上残ったら `AMBIGUOUS_RESULT_SOURCE` で止める——優先順位を決め打ちすると、
    それが暗黙の既定になって「なぜこの図なのか」が説明できなくなる。
    """
    eligible = [c for c in candidates if c["valid"]
                and (source == "auto" or c["source"] == source)
                and (result_id is None or c["result_id"] == result_id)]
    if len(eligible) > 1:
        raise DomainError(
            "AMBIGUOUS_RESULT_SOURCE",
            "図に使う結果を指定してください（source か result_id）。"
            "複数の解析結果が同時に有効なので、どれを描くかは選べません。",
            {"candidates": [{"source": c["source"], "result_id": c["result_id"],
                             "dataset_id": c["dataset_id"]} for c in eligible]})
    if not eligible:
        raise DomainError(
            "ANALYSIS_RESULT_NOT_FOUND",
            "指定に合う有効な解析結果がありません"
            "（前処理をやり直した後の古い結果は選べません）。",
            {"requested_source": source, "requested_result_id": result_id,
             "available": [{"source": c["source"], "result_id": c["result_id"],
                            "valid": c["valid"]} for c in candidates]})
    return eligible[0]


def figure_annotations(ds, result: dict | None = None) -> list[str]:
    """図へ載せる但し書きを返す（無ければ空）。

    PNG はレポートへ貼られた後は単体で読まれる。「探索専用のデータから描いた」も
    「群とバッチが交絡したまま未調整で解析した」も、戻り値ではなく図の中に要る。

    `result` は任意（PCA 図など、比較の来歴を持たない図からは呼ばれない）。
    """
    notes: list[str] = []
    if ds is not None and getattr(ds, "exploratory_only", False):
        notes.append(_EXPLORATORY_NOTE)
    comparison_prov = ((result or {}).get("provenance") or {}).get("comparison") or {}
    if comparison_prov.get("unadjusted_confounded"):
        notes.append(_UNADJUSTED_CONFOUNDED_NOTE)
    return notes


def save_result_figure(ds, result: dict, path: Path, *, kind: str,
                       title: str | None = None) -> Path:
    """指定された結果から図を描き、確定した PNG として保存する。

    `ds` は但し書きの判断にだけ使う（ARF 経路では None）。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fig = _render(result, kind=kind, title=title)
    try:
        notes = figure_annotations(ds, result)
        if notes:
            fig.text(0.5, 0.005, "\n".join(notes), ha="center", va="bottom",
                     fontsize=8, color="#a33")
        _atomic_savefig(fig, path)
    finally:
        plt.close(fig)
    return path


# ---------- 内部 ----------

def _pca_plot_payload(result: dict) -> dict:
    """PCA 結果を、描画が食う散布図 dict（`points` を持つ形）へ揃える。

    PCA の結果は経路によって**形が違う**。ARF 経路（`session.arf.last_pca_plot`）は
    最初から `points` を持つが、DatasetState 経路
    （`analysis.dataset_service.pca_dataset` → `dataset_analysis.run_dataset_pca`）は
    `scores`（name / role / PC1..PCn）しか持たない。ここで揃えずに `points` だけを
    読むと、DatasetState 経路では点が 0 件のまま**枠だけの図**が保存され、しかも
    ファイルは生成されて hash も通るので「必須出力は達成」と数えられてしまう
    （pipeline の `results/pca.png` が実際にそうなっていた）。

    どちらの形でも点が取れなければ描かずに止める。空の散布図は「点が無いデータ」
    ではなく「入力の形が想定と違う」の症状で、図として出せば必ず誤読される。
    """
    if (result or {}).get("points"):
        return result
    projected = dataset_pca_plot(result) if result else None
    if projected is None:
        raise DomainError(
            "PCA_FIGURE_EMPTY",
            "PCA 図に描ける点がありません（points も PC1/PC2 を持つ scores も"
            "見つかりません）。空の散布図は保存しません。",
            {"keys": sorted(result or {})})
    return projected


def _render(result: dict, *, kind: str, title: str | None):
    if kind == "pca":
        xs, ys, labels, x_label, y_label, plot_title = _pca_scatter_arrays(
            _pca_plot_payload(result))
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(xs, ys, alpha=0.6)
        for x, y, label in zip(xs, ys, labels):
            if label:
                ax.annotate(str(label), (x, y), fontsize=8)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title or plot_title)
        return fig
    if kind == "volcano":
        # 描画は arf_plot_volcano（画像返し）と共有する。色や軸がツール間でずれると
        # 「画面で見た図」と「レポートに貼った図」が別物になるため。
        return render_volcano_plot(result, title=title)
    raise DomainError("UNSUPPORTED_FIGURE_KIND",
                      f"未対応の図の種類です: {kind}", {"kind": kind})


def _atomic_savefig(fig, path: Path) -> None:
    """同じ親ディレクトリの一時ファイルへ書いてから置換する。

    描画の途中で落ちても、既にある図を壊れた PNG で置き換えない。
    """
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp",
                                    dir=str(path.parent))
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        fig.savefig(tmp, format="png", dpi=120, bbox_inches="tight")
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
