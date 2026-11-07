"""統計エンジン（構成見直し C）。

入力形式にも MCP にもセッションにも依存しない数値処理の層。
`metabolomix.core.session_state`・`metabolomix.core.mcp_core`・`metabolomix.tools.*` を import しない。

入口: `stats.compare.compare`・`stats.preprocess.preprocess`・`stats.pca.pca`・
`stats.composition.composition`・`stats.spec_v2.run_spec_v2`・`stats.data.from_dataset`。
サブモジュールは遅延させず import しない（`pca` モジュールと `pca` 関数の名前衝突を避ける）。
"""

__all__ = [
    "compare",
    "composition",
    "confounding",
    "data",
    "design",
    "limma",
    "pca",
    "preprocess",
    "spec_v2",
]
