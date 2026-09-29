# scripts/verify_curation_suggest.py
"""候補付けの質を実データで測る（spec §11.1）。研究データのパスは引数で渡す（ここに書かない）。

1. 正解ありの検証: 元レビューの判定が ok・IsReferenceMatched・フラグ無しのスポットを正解とみなし、
   代表を隠して ② を回したとき、正解のレコードが 1 位・3 位以内に入る割合（分子種 / 和組成）。
   ※ ② が測るのは「MS-DIAL 自身の確信した選択（verdict ok かつ IsReferenceMatched かつ無フラグ）を
   再検索が再現する頻度」＝一致度であって、正解率ではない。
2. ④ の再現: MS-DIAL の FoundInUpperMsMs リンク（注釈付き同士）のうち、軽い側 X（断片候補）と重い側 Y
   （前駆体）の向きに限り、Δm/z の関係表で名指しできた割合。リンクは相互なので、向きを限らないと
   同じ対が 2 回数えられ、断片になりえない「重い側 X」の向きも混ざる。
非公開ヘルパ（suggest._scoring / relations._explanations）を直接呼ぶ。open_store は照合用キャッシュを
LIPIDMIX_LIBRARY_CACHE_DIR に作る／更新することがある（.arf2 の隣には何も書かない）。
使い方: python scripts/verify_curation_suggest.py --arf2 <AlignmentResult_*.arf2> --library <*_Loaded.msp2.dbs>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lipidmix.arf2.ion_features import load_ion_features            # noqa: E402
from lipidmix.arf2.match_results import load_spot_annotations       # noqa: E402
from lipidmix.arf2.reader import load_catalog                       # noqa: E402
from lipidmix.curation import candidates, evidence, judge, relations, review, suggest, trend  # noqa: E402
from lipidmix.library.store import open_store                       # noqa: E402
from lipidmix.msdial.analysis_params import resolve_analysis_params # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arf2", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--limit", type=int, default=400)
    args = parser.parse_args(argv)
    th = judge.resolve_thresholds(None)
    store = open_store(args.library)
    catalog = load_catalog(args.arf2)
    annotations = load_spot_annotations(args.arf2)
    spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
    base = review.run_review(args.arf2, spots, store=store, ms2_tol=0.025, th=th, file_ids=None,
                             max_traces=1, selection={})
    all_truth = [s for s in base["spots"] if s["verdict"] == "ok" and (s.get("match") or {}).get("is_reference_matched")
                 and not s.get("flag")]
    truth = all_truth[: args.limit]  # カタログ順の先頭 N 件（無作為ではない）
    by_id = {c["MasterAlignmentID"]: c for c in catalog}
    print(f"正解候補（ok ∧ IsReferenceMatched ∧ 無フラグ）{len(all_truth)} 件のうち、カタログ順の先頭 {len(truth)} 件を使う")
    scoring = suggest._scoring(store)
    evs, _ = evidence.collect(args.arf2, [by_id[t["spot_id"]] for t in truth], store=store, ms2_tol=scoring["ms2_tol"], th=th,
                              max_traces=1, keep_measured=True)
    hits = {"species_top1": 0, "species_top3": 0, "sum_top1": 0, "sum_top3": 0, "n": 0, "n_sum": 0}
    for ev in evs:
        rep = annotations[ev["spot_id"]]["representative"]
        out = candidates.build_library_candidates(ev, msdial_matches=[], representative=None, store=store,
                                                  scoring=scoring, trends={"classes": {}}, th=th,
                                                  exclude_current=False, top_n=3)
        names = [c["name"] for c in out["candidates"]]
        sums = [c["sum_name"] for c in out["candidates"]]
        truth_sum = trend.sum_composition(rep["name"])
        hits["n"] += 1
        hits["n_sum"] += bool(truth_sum)  # 和組成が読めない名前（脂質名でない等）は和組成の分母に入れない
        # MS-DIAL の名前は "和組成|分子種" のことがあるので、| で分けた各部分のどれかに一致すれば正解とみなす
        truth_names = {part.strip() for part in rep["name"].split("|")}
        hits["species_top1"] += bool(names[:1]) and names[0] in truth_names
        hits["species_top3"] += any(name in truth_names for name in names)
        hits["sum_top1"] += bool(truth_sum) and sums[:1] == [truth_sum]
        hits["sum_top3"] += bool(truth_sum) and truth_sum in sums
    print("② 一致度の検証（代表を隠した再検索が MS-DIAL の確信した選択を再現するか。正解率ではない）")
    print(f"  和組成が読めず和組成の分母から除いた件数	{hits['n'] - hits['n_sum']}（分子種は {hits['n']} 件、和組成は {hits['n_sum']} 件が分母）")
    for key in ("species_top1", "species_top3", "sum_top1", "sum_top3"):
        denominator = hits["n_sum"] if key.startswith("sum") else hits["n"]
        print(f"  {key}\t{hits[key]}/{denominator}\t{hits[key] / max(denominator, 1):.1%}")

    features = load_ion_features(args.arf2)
    params = resolve_analysis_params(args.arf2, next((s.get("IonMode") for s in catalog), None))
    named = total = 0
    for spot_id, feature in features.items():
        for link in feature["links"]:
            if link["kind"] != "found_in_upper_msms" or link["spot_id"] not in by_id:
                continue
            x, y = by_id[spot_id], by_id[link["spot_id"]]
            if not x["MassCenter"] < y["MassCenter"]:  # 軽い側 X が断片、重い側 Y が前駆体の向きだけ数える
                continue
            if (annotations.get(y["MasterAlignmentID"]) or {}).get("representative") is None:
                continue
            total += 1
            found = relations._explanations(
                {"mz": x["MassCenter"]}, {"mz": y["MassCenter"], "adduct": y.get("AdductType")},
                params["searched_adducts"], 0.010)
            named += bool(found)
    print(f"④ FoundInUpperMsMs（軽い側 X → 重い側の注釈付き Y、対ごとに 1 回）のうち関係表で名指しできた割合\t{named}/{total}\t{named / max(total, 1):.1%}")
    store.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
