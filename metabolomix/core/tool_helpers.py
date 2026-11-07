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


def dataset_pca_plot(last_pca: dict, title: str = "PCA (mzTab-M)") -> dict | None:
    """DatasetState の PCA 結果を `_pca_scatter_arrays` が食う散布図 dict に射影する。

    `ds.last_pca["scores"]`（name / role / PC1 / PC2…）を、描画側が読む 1 つの形へ寄せる。

    PC2 が無い（主成分 1 つ）場合は散布図にできないので None を返す。
    """
    scores = (last_pca or {}).get("scores") or []
    points = []
    for row in scores:
        x, y = row.get("PC1"), row.get("PC2")
        if x is None or y is None:
            continue
        points.append({"x": float(x), "y": float(y), "label": row.get("name")})
    if not points:
        return None
    evr = (last_pca or {}).get("explained_variance_ratio") or []
    x_label = f"PC1 ({evr[0] * 100:.2f}%)" if len(evr) > 0 else "PC1"
    y_label = f"PC2 ({evr[1] * 100:.2f}%)" if len(evr) > 1 else "PC2"
    return {"title": title, "x_label": x_label, "y_label": y_label, "points": points}
