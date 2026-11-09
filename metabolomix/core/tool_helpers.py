"""ツール本体が使う小さな純ロジック（identity の対応表・PCA の散布図の整形）。

依存: lipid_identity（下位レイヤ）。tools_* / server / session は import しない。
"""
from metabolomix.msdial import lipid_identity


_IDENTITY_TABLES = None


def _identity_tables():
    """同梱の RefMet/LIPID MAPS 対応表を1回だけ読み込みキャッシュする。"""
    global _IDENTITY_TABLES
    if _IDENTITY_TABLES is None:
        _IDENTITY_TABLES = lipid_identity.load_reference_tables()
    return _IDENTITY_TABLES


def _pca_scatter_arrays(plot: dict):
    """`dataset_pca_plot` の散布図 dict から配列とラベルを取り出す（純ロジック）。"""
    points = plot.get("points", [])
    xs = [float(p["x"]) for p in points]
    ys = [float(p["y"]) for p in points]
    labels = [p.get("label") for p in points]
    return (
        xs, ys, labels,
        plot.get("x_label", "PC1"),
        plot.get("y_label", "PC2"),
        plot.get("title", "PCA"),
    )


def is_stat_pca(record: dict | None) -> bool:
    """統計エンジンの PCA の記録（`stats_results.register_pca` が置く `kind="stat_pca"`）か。"""
    return (record or {}).get("kind") == "stat_pca"


#: 統計エンジンの PCA で、QC の点につける群の名前（凡例に出る）。
QC_LABEL = "QC"


def stat_pca_rows(record: dict) -> list[dict]:
    """`stat_pca` の記録（stats-pca.v1 を `engine` に持つ）を、スコア図が読む行（個体 → QC の順）にする。

    行は `{id, name, role, group, fitted, PC1, PC2, ...}`。`id` はエンジンの個体 ID（QC は QC の試料 ID）、
    `name` はその表示名（記録の `labels`。無ければ ID）、`group` は
    記録の `groups`（色分けの表示だけに使った群。QC は `QC`）、`fitted` は軸の計算に使った行か
    （QC は投影だけなので False）。
    """
    engine = record.get("engine") or {}
    groups = record.get("groups") or {}
    labels = record.get("labels") or {}
    scores = engine.get("scores") or {}
    rows = []
    for names, matrix, role in ((engine.get("unit_ids") or [], scores.get("units") or [], "sample"),
                                (engine.get("qc_ids") or [], scores.get("qc") or [], "qc")):
        for name, values in zip(names, matrix):
            row = {"id": name, "name": labels.get(name, name), "role": role, "fitted": role == "sample",
                   "group": groups.get(name) if role == "sample" else QC_LABEL}
            row.update({f"PC{k + 1}": float(v) for k, v in enumerate(values)})
            rows.append(row)
    return rows


def dataset_pca_plot(last_pca: dict, title: str | None = None) -> dict | None:
    """統計エンジンの PCA の記録（`kind="stat_pca"`）を `_pca_scatter_arrays` が食う散布図 dict に射影する。

    `engine` の個体と QC のスコア（QC は `fitted=False`・群 `QC`）を、描画側が読む 1 つの形へ寄せる。
    PCA の記録でない・PC2 が無い（主成分 1 つ）場合は散布図にできないので None を返す。
    """
    if not is_stat_pca(last_pca):
        return None
    scores = stat_pca_rows(last_pca)
    evr = (last_pca.get("engine") or {}).get("explained_variance_ratio") or []
    points = []
    for row in scores:
        x, y = row.get("PC1"), row.get("PC2")
        if x is None or y is None:
            continue
        points.append({"x": float(x), "y": float(y), "label": row.get("name"),
                       "group": row.get("group") or "(none)", "fitted": bool(row.get("fitted"))})
    if not points:
        return None
    if title is None:
        title = "PCA"
    x_label = f"PC1 ({evr[0] * 100:.2f}%)" if len(evr) > 0 else "PC1"
    y_label = f"PC2 ({evr[1] * 100:.2f}%)" if len(evr) > 1 else "PC2"
    return {"title": title, "x_label": x_label, "y_label": y_label, "points": points}
