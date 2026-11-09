"""偽の規則プラグイン（パッケージ `msdial_rules` の面だけを持つ小さな実装）。

本物の規則の写しではない。MCP の配管（変換・記録・封筒・ハンドル）を確かめるための
最小の実装で、`docs/quant-rules-plugin.md` の 4 関数と `API_VERSION` / `RULES_VERSION` だけを
表で決めた単純な規則で満たす。名前（`IS-A` など）もクラス（`PC` `PE`）も人工。
`write_fake_rules(tmp_path)` が tmp に `msdial_rules/{__init__,api}.py` を書き出し、
そのフォルダを返す（設定の `[quant] rules` に向けて使う）。
"""
from __future__ import annotations

from pathlib import Path

_INIT = '"""fake msdial_rules (test double)."""\n'

_API = '''"""fake rules for tests: not the real rules."""
import csv
import io

API_VERSION = {api_version}
RULES_VERSION = "fake-rules-1"
_ANY = "Any others"
_ALIAS = {{"Cer_NS": "Cer-NS", "Cer_AS": "Cer-AS"}}  # class aliases (the real rules resolve aliases)


def _canon(cls):
    return _ALIAS.get(cls, cls)


class FakeRuleError(ValueError):
    def __init__(self, message, code, details=None):
        super().__init__(message)
        self.code = code
        self.details = details or {{}}


def _rows(text):
    first = text.splitlines()[0] if text.strip() else ""
    delim = "\t" if "\t" in first else ","
    return list(csv.DictReader(io.StringIO(text), delimiter=delim))


def _table(reference_text):
    rows = _rows(reference_text)
    need = {{"TargetClass", "StandardName", "Coefficient", "Unit"}}
    if not rows or not need <= set(rows[0]):
        raise FakeRuleError("reference table needs TargetClass,StandardName,Coefficient,Unit",
                            "reference_invalid")
    return {{_canon(r["TargetClass"]): r for r in rows}}


def normalize(features, n_samples, reference_text, method="cell_coefficient", bindings=None):
    bindings = bindings or {{}}
    table = _table(reference_text)
    standards = {{r["StandardName"] for r in table.values()}}
    for f in features:
        if len(f["heights"]) != n_samples:
            raise FakeRuleError("heights length differs from n_samples", "bad_features",
                                {{"feature_id": f["feature_id"]}})
    roles = ["internal_standard" if f["name"] in standards else "target" for f in features]
    classes = sorted({{_canon(f["lipid_class"]) for f, r in zip(features, roles) if r == "target"}})
    assignments, any_others, key_of = {{}}, [], {{}}
    for cls in classes:
        row = table.get(cls)
        if row is None:
            row = table.get(_ANY)
            if row is not None:
                any_others.append(cls)
        if row is None:
            key_of[cls] = None  # class_unmapped
            continue
        key = row["TargetClass"]
        key_of[cls] = key
        if key in assignments:
            continue
        cands = [f for f in features if f["name"] == row["StandardName"]]
        entry = {{"standard_name": row["StandardName"], "feature_id": None,
                 "coefficient": float(row["Coefficient"]), "concentration": None,
                 "standard_adduct": None, "unit": row["Unit"], "status": "ok",
                 "candidates": [{{k: c[k] for k in ("feature_id", "name", "adduct", "mz", "rt")}}
                                for c in cands]}}
        chosen = None
        if len(cands) == 1:
            chosen = cands[0]
        elif len(cands) > 1:
            chosen = next((c for c in cands if c["feature_id"] == bindings.get(row["StandardName"])), None)
            if chosen is None:
                entry["status"] = "standard_ambiguous"
        else:
            entry["status"] = "standard_missing"
        if chosen is not None:
            entry["feature_id"] = chosen["feature_id"]
            entry["standard_adduct"] = chosen["adduct"]
        assignments[key] = entry
    by_id = {{f["feature_id"]: f for f in features}}
    values, statuses, fallback, by_feature = [], [], [], []
    for f, role in zip(features, roles):
        if role != "target":
            values.append([None] * n_samples)
            statuses.append([role] * n_samples)
            fallback.append(False)
            by_feature.append(None)
            continue
        cls = _canon(f["lipid_class"])
        key = key_of[cls]
        by_feature.append(key)
        fallback.append(cls in any_others)
        if key is None:
            values.append([None] * n_samples)
            statuses.append(["class_unmapped"] * n_samples)
            continue
        a = assignments[key]
        if a["status"] != "ok":
            values.append([None] * n_samples)
            statuses.append([a["status"]] * n_samples)
            continue
        std = by_id[a["feature_id"]]
        vrow, srow = [], []
        for s in range(n_samples):
            h, ref = f["heights"][s], std["heights"][s]
            if ref is None or ref <= 0:
                vrow.append(None)
                srow.append("standard_invalid_height")
            elif h is None:
                vrow.append(None)
                srow.append("missing")
            else:
                vrow.append(h / ref * a["coefficient"])
                srow.append("ok" if f["detected"][s] else "gap_filled")
        values.append(vrow)
        statuses.append(srow)
    unit = next(iter(table.values()))["Unit"]
    return {{"unit": unit, "normalization_method": method, "values": values,
            "statuses": statuses, "roles": roles, "assignments": assignments,
            "any_others_classes": any_others, "fallback_by_feature": fallback,
            "target_class_by_feature": by_feature}}


def coefficient_reference_text(reference_text, concentrations, solvent_volume, cells_millions,
                               unit=None):
    rows = _rows(reference_text)
    if not rows or not {{"TargetClass", "StandardName"}} <= set(rows[0]):
        raise FakeRuleError("reference table needs TargetClass,StandardName", "reference_invalid")
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\\n")
    w.writerow(["TargetClass", "StandardName", "Coefficient", "Unit", "StandardAdduct"])
    for r in rows:
        name = r["StandardName"]
        if name not in concentrations:
            raise FakeRuleError(f"no concentration for {{name}}", "concentration_missing",
                                {{"standard_name": name}})
        w.writerow([r["TargetClass"], name, concentrations[name] * solvent_volume / cells_millions,
                    unit or "pmol / 10^6 cells", r.get("StandardAdduct", "")])
    return out.getvalue()


def _strip(name):
    for suffix in ("_pos", "_neg"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def match_samples(names_by_mode, mapping_rows=None):
    """Mirrors msdial_rules: full position-aligned ids, rule choice, mismatched_samples."""
    if mapping_rows is not None:
        by_key = {{(r["mode"], r["source_sample"]): r for r in mapping_rows}}
        expected = {{(m, n) for m, names in names_by_mode.items() for n in names}}
        if len(by_key) != len(mapping_rows) or set(by_key) != expected:
            raise FakeRuleError("mapping must cover exactly the input samples", "mapping_invalid")
        ids = {{m: [by_key[(m, n)]["sample_id"] for n in names]
               for m, names in names_by_mode.items()}}
        rows = {{m: {{n: by_key[(m, n)] for n in names}} for m, names in names_by_mode.items()}}
        rule = "mapping_csv"
    else:
        rows = None
        if len({{frozenset(v) for v in names_by_mode.values()}}) == 1:
            rule, ids = "exact_name", {{m: list(v) for m, v in names_by_mode.items()}}
        else:
            stripped = {{m: [_strip(n) for n in v] for m, v in names_by_mode.items()}}
            if any("" in v or len(set(v)) != len(v) for v in stripped.values()):
                raise FakeRuleError("sample names are ambiguous", "sample_names_ambiguous")
            if len({{frozenset(v) for v in stripped.values()}}) != 1:
                rule, ids = "exact_name", {{m: list(v) for m, v in names_by_mode.items()}}
            else:
                rule, ids = "polarity_token_removed", stripped
    mismatched = {{}}
    for m, v in ids.items():
        others = [set(o) for k, o in ids.items() if k != m]
        missing = [i for i in v if not all(i in o for o in others)]
        if missing:
            mismatched[m] = missing
    return {{"rule": rule, "ids_by_mode": ids, "mismatched": mismatched, "rows_by_mode": rows}}


def _exclusion_reason(name, comment, msms_assigned):
    """Mirrors msdial_rules.selection.exclusion_reason (same order and reason names)."""
    if name.strip() == "Unknown":
        return "name_unknown"
    if "w/o" in name:
        return "name_without_ms2"
    if "RIKEN" in name:
        return "name_riken"
    if "(d" in name:
        return "name_deuterated"
    if "low score" in name.casefold():
        return "name_low_score"
    if comment.strip().startswith("x"):
        return "comment_excluded"
    if msms_assigned.strip().casefold() == "false":
        return "msms_not_assigned"
    return None


def select(features, adduct_text):
    """Mirrors msdial_rules.selection.select_feature: internal standard, then name /
    comment / MS-MS exclusions, then the adduct table (None = filter_disabled)."""
    chosen = None
    if adduct_text is not None:
        chosen = set()
        for r in csv.DictReader(io.StringIO(adduct_text)):
            if r["IsSelected"].strip().lower() == "true":
                chosen.add((_canon(r["Class"].strip()), r["Adduct"].strip(),
                            r["Ion mode"].strip().lower()))
    out = []
    for f in features:
        reason = _exclusion_reason(f["name"], f["comment"], f["msms_assigned"])
        if f["role"] == "internal_standard":
            action, reason = "excluded", "internal_standard"
        elif reason:
            action = "excluded"
        elif chosen is None:
            action, reason = "retained", "filter_disabled"
        elif (_canon(f["ontology"]), f["adduct"], f["mode"]) in chosen:
            action, reason = "retained", "reference_selected"
        else:
            action, reason = "excluded", "reference_not_selected"
        out.append({{"feature_id": f["feature_id"], "action": action, "reason": reason,
                    "original_class": f["ontology"], "canonical_class": _canon(f["ontology"]),
                    "adduct": f["adduct"], "mode": f["mode"],
                    "class_rule_version": RULES_VERSION}})
    return out
'''


#: conftest が `--quant-rules` 付きのときだけ置く本物の規則のフォルダ（既定は None＝偽を書く）。
REAL_RULES_FOLDER: Path | None = None


def write_fake_rules(tmp_path, api_version: int = 1, extra: str = "") -> Path:
    """tmp に `msdial_rules/` を書いてそのフォルダを返す。`extra` は api.py の末尾に足す。

    `--quant-rules` 付きの定量スイートでは、書かずに本物の規則のフォルダを返す。
    """
    if REAL_RULES_FOLDER is not None and api_version == 1 and not extra:
        return REAL_RULES_FOLDER
    folder = Path(tmp_path) / "msdial_rules"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "__init__.py").write_text(_INIT, encoding="utf-8")
    (folder / "api.py").write_text(_API.format(api_version=api_version) + extra, encoding="utf-8")
    return folder
