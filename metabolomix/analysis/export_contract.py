"""差次的エクスポートの契約（列定義・版・数値書式・メタ行）。

このファイルの EXPORT_COLUMNS と CONTRACT_VERSION は**別リポとの契約**。
massbank-context の load_differential → pathway_activity がこの形を前提に読む。
列の増減・改名・順序変更は CONTRACT_VERSION の引き上げと下流の同時更新なしに
やってはいけない。

ARF 経路（arf_export_differential）と DatasetState 経路
（dataset_export_differential）の両方がここを参照する。同じ契約の実装が
2 箇所にあると、片方だけ直った状態で下流が壊れる。

依存は stdlib のみ。metabolomix.* を import しない leaf。
"""
from __future__ import annotations

import math
from datetime import datetime, timezone

EXPORT_COLUMNS = [
    "spot_id", "name", "name_source", "ontology", "inchikey",
    "inchikey_source", "msi_level", "mz", "rt", "log2fc",
    "p_value", "q_value", "mean_a", "mean_b", "significant",
]
CONTRACT_VERSION = 1
LOG2FC_SIGN = "positive means group_b is higher"


def format_number(value, format_spec: str) -> str:
    """有限の数値だけを書き出し、欠測・NaN・inf は空欄にする。"""
    if value is None:
        return ""
    number = float(value)
    if not math.isfinite(number):
        return ""
    return format(number, format_spec)


def is_significant(*, q, log2fc, q_threshold, log2fc_threshold) -> bool:
    """契約 15 列目 `significant` の唯一の判定。

    dict ではなく**値**をキーワードで受ける。ARF 経路は結果行に `q_value`、
    DatasetState 経路は `q` というキーで同じ量を持っており、キー名を契約側に
    持ち込むと経路ごとの分岐が再発する。呼び出し側が自分のキーから値を取り出す。

    欠測・NaN・inf は有意にしない。下流（massbank-context）はこの真偽値を
    そのまま濃縮対象の選別に使うため、判定不能を偽として畳む。
    """
    if q is None or log2fc is None:
        return False
    q = float(q)
    log2fc = float(log2fc)
    if not (math.isfinite(q) and math.isfinite(log2fc)):
        return False
    return q <= q_threshold and abs(log2fc) >= log2fc_threshold


def build_meta(*, group_a, n_a, group_b, n_b,
               q_threshold, log2fc_threshold, log_transform,
               n_features_total, n_with_inchikey, n_unannotated,
               msi_note, source_lines=(), preprocess_line=None) -> list[str]:
    """メタ行ブロックを契約の順序で組む。

    行の**順序も契約の一部**なので、経路ごとに固有の行（source 系・preprocess）は
    自由な位置に足させず、決められたスロットに差す:

        1. contract_version
        2. exported_at
        3. *source_lines        ← 経路固有（# source_arf / # source_mztab など）
        4. group_a / n_a
        5. group_b / n_b
        6. log2fc_sign
        7. q_threshold / log2fc_threshold / log_transform
        8. preprocess_line      ← 経路固有（省略可）
        9. n_features_total / n_with_inchikey / n_unannotated
       10. msi_note             ← 文面が経路で違う（.arf2 由来 か mzTab-M か）

    この順序は arf_export_differential の現行出力と一致させてある。
    変えると既存ファイルとの差分が出るので、CONTRACT_VERSION を上げずに
    順序を触ってはいけない。
    """
    lines = [
        f"# contract_version = {CONTRACT_VERSION}",
        f"# exported_at = {datetime.now(timezone.utc).isoformat()}",
        *source_lines,
        f"# group_a = {group_a}\tn_a = {n_a}",
        f"# group_b = {group_b}\tn_b = {n_b}",
        f"# log2fc_sign = {LOG2FC_SIGN}",
        f"# q_threshold = {q_threshold}\tlog2fc_threshold = {log2fc_threshold}"
        f"\tlog_transform = {str(bool(log_transform)).lower()}",
    ]
    if preprocess_line is not None:
        lines.append(preprocess_line)
    lines += [
        f"# n_features_total = {n_features_total}"
        f"\tn_with_inchikey = {n_with_inchikey}"
        f"\tn_unannotated = {n_unannotated}",
        msi_note,
    ]
    return lines


def format_row(row: dict) -> str:
    """EXPORT_COLUMNS の順に 1 行を組む。row は 15 列すべてのキーを持つこと。"""
    return "\t".join([
        str(row["spot_id"]),
        row["name"] or "",
        row["name_source"] or "",
        row["ontology"] or "",
        row["inchikey"] or "",
        row["inchikey_source"] or "",
        str(row["msi_level"]) if row["msi_level"] is not None else "",
        format_number(row["mz"], ".4f"),
        format_number(row["rt"], ".4f"),
        format_number(row["log2fc"], ".6f"),
        format_number(row["p_value"], ".6g"),
        format_number(row["q_value"], ".6g"),
        format_number(row["mean_a"], ".6g"),
        format_number(row["mean_b"], ".6g"),
        "true" if row["significant"] else "false",
    ])
