"""差次的エクスポートの契約（列定義・版・数値書式・メタ行）。

このファイルの EXPORT_COLUMNS と CONTRACT_VERSION は**別リポとの契約**。
massbank-context の load_differential → pathway_activity がこの形を前提に読む。
列の増減・改名・順序変更は CONTRACT_VERSION の引き上げと下流の同時更新なしに
やってはいけない。

旧 ARF 経路と DatasetState 経路
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
    """契約 15 列目 `significant` の q・log2fc の規則（唯一の定義）。結果の行から判定するときは
    `row_significant` を使う（統計エンジンの行は判定 `verdict` も見る。R16）。

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


#: 統計エンジンの結果の `significant`（R16）をメタ行 `significant_rule` に書く文面（ASCII）。
SIGNIFICANT_RULE = "q<=q_threshold & |log2fc|>=log2fc_threshold & verdict=changed"


def row_significant(row: dict, *, q_threshold, log2fc_threshold) -> bool:
    """結果の行 1 つの `significant`（書き出しの 15 列目・根拠バンドル・volcano の色・新しい知見の候補）。

    行が `significant` を持てばその値（`stats_results.as_two_group` が記録の閾値で 1 回だけ付けたもの。
    消費側はそれを読み、閾値から計算し直さない）。持たなければ `is_significant` で、`verdict` を持つ行
    （統計エンジン）は加えて `verdict == "changed"` のときだけ有意（R16: 検出が足りず量を判定できない
    `inconclusive` の行は、q が小さく |log2fc| が大きくても有意にしない）。`verdict` の無い v1 の行は
    `is_significant` のまま。
    """
    flag = row.get("significant")
    if flag is not None:
        return bool(flag)
    if not is_significant(q=row.get("q"), log2fc=row.get("log2fc"),
                          q_threshold=q_threshold, log2fc_threshold=log2fc_threshold):
        return False
    verdict = row.get("verdict")
    return verdict is None or verdict == "changed"


def build_meta(*, group_a, n_a, group_b, n_b,
               q_threshold, log2fc_threshold, log_transform,
               n_features_total, n_with_inchikey, n_unannotated,
               msi_note, source_lines=(), preprocess_line=None,
               stats_lines=()) -> list[str]:
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
       11. *stats_lines         ← 統計エンジンの結果のときだけ（`stats_meta_lines`）。**末尾に足す**

    11 は契約の版を上げずに足した行で、`msi_note` の後ろにだけ置く（既存の行の位置は動かさない）。
    下流（massbank-context の `parse_differential`）はタブ区切りの `# key = value` を全部拾い、
    知らないキーは読み飛ばす。既存のキーと同じ名前は使わない（`test_export_contract` が固定）。

    この順序は `dataset_export_differential` の現行出力と一致させてある。
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
        *stats_lines,
    ]
    return lines


#: `statistics.subject_source`（`stats/preprocess.py` `_collapse`）→ メタ行の `unit`（3 値）。
#: `subject_id`＝全試料が被験者 ID を持ち技術反復を 1 個体にまとめた（比較の単位は被験者）、
#: `injection`＝被験者 ID が無い・まとめを切った（単位は注入）、`mixed`＝混在（そのまま書く）。
_UNIT_OF_SUBJECT_SOURCE = {"subject_id": "subject", "injection": "injection", "mixed": "mixed"}

MEAN_DEFINITION = "linear arithmetic mean over units after preprocessing"


def is_engine_result(result: dict) -> bool:
    """統計エンジンの結果（`as_two_group` の出力）か。来歴の `kind` が `compare` のもの。
    `statistics` の有無では見ない（v1 の結果が同名の欄を持つ可能性に依存しない）。"""
    return ((result or {}).get("provenance") or {}).get("kind") == "compare"


def stats_meta_lines(result: dict) -> list[str]:
    """統計エンジンの結果（`as_two_group` の出力）のメタ行。v1 の結果は [] 。

    `msi_note` の後に足す行（契約の版は据え置き）。キーは既存のメタ行と重ならない
    （`effect_size_definition` は pipeline v2 の文字列版と同名なので、エンジンの結果ではこちらだけが書く）。
    `mean_a` / `mean_b` は個体（単位）ごとの前処理後の線形値の算術平均（`mean_definition`）、
    `n_a` / `n_b` は注入の数のまま、`n_units_a` / `n_units_b` が比較に入った単位の数。
    `contrast_id` は `matrix_id` の行の末尾（どの対比を書いたかを残す）。`bh_n` は BH の母集団が空
    （対比の n が無い）なら 0。`significant_rule` は 15 列目 `significant` の規則（R16。`row_significant`）。
    """
    if not is_engine_result(result):
        return []
    stats = result["statistics"]
    rule = stats["verdict_rule"]
    bh = stats["bh_population"]
    unit = unit_of(stats.get("subject_source"))
    bh_n = bh.get("n") or 0
    return [
        f"# test = {stats['test']}\teffect_size_definition = {result['effect_size_definition']['name']}"
        f"\tci_level = {stats['ci_level']}",
        f"# verdict_q_threshold = {rule['q_threshold']}\tequivalence_fc = {rule['equivalence_fc']}"
        f"\tbh_population = {bh['name']}\tbh_n = {bh_n}\tq_interpretation = {bh['q_interpretation']}",
        f"# unit = {unit}\tn_units_a = {stats['n_units_a']}\tn_units_b = {stats['n_units_b']}",
        f"# mean_definition = {MEAN_DEFINITION}",
        f"# matrix_id = {result['matrix_id']}\trule_hash = {result['rule_hash']}"
        f"\tcontrast_id = {result['contrast_id']}",
        f"# significant_rule = {SIGNIFICANT_RULE}",
    ]


def unit_of(subject_source) -> str:
    """`subject_source` → メタ行の `unit`（`subject` / `injection` / `mixed` のどれか）。
    None・未知の値は `injection`（個体にまとめていない扱い。"None" という文字列は書かない）。"""
    return _UNIT_OF_SUBJECT_SOURCE.get(subject_source, "injection")


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
