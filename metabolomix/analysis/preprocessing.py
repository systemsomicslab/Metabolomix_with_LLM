"""試料名からの役割（sample / qc / blank）の推定（純ロジック層、MCP 非依存）。

ファイル名と MS-DIAL の Class ID のトークンで役割を決める。前処理そのもの（ブランク・正規化・
ドリフト補正・足切り）は統計エンジン（`analysis/stats/preprocess.py`）が持つ。
"""

from __future__ import annotations

DEFAULT_ROLE_TOKENS: dict[str, set[str]] = {
    "qc": {"qc"},
    "blank": {"blank"},
}


def _segments(text: str) -> set[str]:
    return {seg for seg in str(text).casefold().split("_") if seg}


def detect_sample_roles(
    sample_names: list[str],
    class_ids: dict[str, str] | None = None,
    config: dict | None = None,
) -> dict[str, str]:
    """各サンプルを sample/qc/blank に分類する。

    ファイル名と Class ID の `_` 区切りセグメントを大小無視でトークン照合する。
    blank を qc より優先評価する（曖昧語衝突を避けるため決定的順序）。
    """
    cfg = config or {}
    tokens = {
        "qc": {t.casefold() for t in cfg.get("qc_tokens", DEFAULT_ROLE_TOKENS["qc"])},
        "blank": {t.casefold() for t in cfg.get("blank_tokens", DEFAULT_ROLE_TOKENS["blank"])},
    }
    class_ids = class_ids or {}
    roles: dict[str, str] = {}
    for name in sample_names:
        segs = _segments(name) | _segments(class_ids.get(name, ""))
        if segs & tokens["blank"]:
            roles[name] = "blank"
        elif segs & tokens["qc"]:
            roles[name] = "qc"
        else:
            roles[name] = "sample"
    return roles
