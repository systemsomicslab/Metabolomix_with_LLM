"""脂質アッセイの注釈の判定（`lipid_identity.lipid_annotation_class` / `lipid_name_parses`）を、実データの mzTab-M の
全注釈で書き出し、2 つの書き出しが同じかを比べる CLI（テストではない。速くする前後で判定が変わらないことを確かめる）。

    python scripts/verify_lipid_judgement.py dump <mzTab-M> <out.json>
    python scripts/verify_lipid_judgement.py compare <before.json> <after.json>

`dump` は SML 行の (chemical_name, opt_global_Ontology, chemical_formula) の組ごとに、組成（`composition`）が見るのと
同じ引数で `lipid_annotation_class` を呼び（名前の無い注釈は組成が判定しないので除く）、名前ごとに `lipid_name_parses` も
呼ぶ。プロセスを起こしたばかりの冷えた状態での時間（名前の判定だけ・注釈の判定全体）と GOSLIN を呼んだ回数も書く。
`compare` は判定の表が 1 つ残らず同じかを見て、違えば例を出して 1 で終わる。
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path


def _sml_rows(path: Path):
    header = None
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("SMH\t"):
                header = line.rstrip("\n").split("\t")
            elif line.startswith("SML\t") and header is not None:
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def dump(mztab: Path, out: Path) -> None:
    from metabolomix.msdial import lipid_identity as li

    triples = sorted({(r.get("chemical_name", ""), r.get("opt_global_Ontology", ""), r.get("chemical_formula", ""))
                      for r in _sml_rows(mztab)})
    names = sorted({t[0] for t in triples if li.has_annotation_name(t[0])})
    calls = {"n": 0}
    original = li._parse_lipid

    def counted(text):
        calls["n"] += 1
        return original(text)

    li._parse_lipid = counted
    try:
        start = time.perf_counter()
        parses = {name: li.lipid_name_parses(name) for name in names}
        t_names = time.perf_counter() - start
        n_names = calls["n"]
        start = time.perf_counter()
        judged = {"\t".join(t): list(li.lipid_annotation_class(*t)) for t in triples if li.has_annotation_name(t[0])}
        t_all = time.perf_counter() - start
    finally:
        li._parse_lipid = original
    out.write_text(json.dumps({
        "n_triples": len(judged), "n_names": len(names),
        "seconds_lipid_name_parses": round(t_names, 2), "seconds_lipid_annotation_class_after": round(t_all, 2),
        "goslin_calls_names": n_names, "goslin_calls_total": calls["n"],
        "lipid_name_parses": parses, "lipid_annotation_class": judged}, ensure_ascii=False, sort_keys=True),
        encoding="utf-8")
    print(f"triples {len(judged)}  names {len(names)}  lipid_name_parses {t_names:.1f} s ({n_names} GOSLIN calls)  "
          f"then lipid_annotation_class {t_all:.1f} s ({calls['n'] - n_names} more GOSLIN calls)")


def compare(before: Path, after: Path) -> int:
    a = json.loads(before.read_text(encoding="utf-8"))
    b = json.loads(after.read_text(encoding="utf-8"))
    bad = 0
    for key in ("lipid_name_parses", "lipid_annotation_class"):
        if a[key].keys() != b[key].keys():
            print(f"{key}: different keys ({len(a[key])} vs {len(b[key])})")
            bad += 1
            continue
        diff = [k for k in a[key] if a[key][k] != b[key][k]]
        print(f"{key}: {len(a[key])} entries, {len(diff)} differ")
        for k in diff[:20]:
            print(f"   {k!r}: {a[key][k]} -> {b[key][k]}")
        bad += len(diff)
    for key in ("seconds_lipid_name_parses", "seconds_lipid_annotation_class_after", "goslin_calls_total"):
        print(f"{key}: {a.get(key)} -> {b.get(key)}")
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "dump":
        dump(Path(sys.argv[2]), Path(sys.argv[3]))
    elif len(sys.argv) == 4 and sys.argv[1] == "compare":
        sys.exit(compare(Path(sys.argv[2]), Path(sys.argv[3])))
    else:
        print(__doc__)
        sys.exit(2)
