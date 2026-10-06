"""実データ突き合わせ（受け入れゲート）— 移植した MS/MS 照合スコアを MS-DIAL 自身の
mzTab-M 出力（`id_confidence_measure[4..8]`）と突き合わせる CLI。

spec: docs/superpowers/sdd/2026-09-19-msms-spectral-matching タスク 7。

やること（brief 準拠）:
  1. `.dbs` から store を作る（`search_params` の `ms2_tolerance` を bin_width に使う）。
  2. mzTab の SME 行を読む（`metabolomix.mztab.reader.parse_mztab`）。
  3. 各 SME 行の `spectra_ref`（`ms_run[N]:ms1scanID=X ms2scanID=Y` の `| ` 区切りリスト）
     が指す測定を `.dcl` から引く。**`ms2scanID` は上流 `MSDecResult.RawSpectrumID`
     と同一**（`AddSpectraRef` が `properties[i].MS2RawSpectrumID` を書き、
     `MSDecResult.RawSpectrumID = targetSpecID = chromPeakFeature.MS2RawSpectrumID`
     — `Ms2Dec.cs` / `DataAccess.GetTargetCEIndexForMS2RawSpectrum`）。
     こちらの `.dcl` reader の `raw_spec_id` フィールドがこれに対応する
     （`scan_id` ではない。`MsdecResultsWriter.SaveScanData` が
     `ScanID` の直後に `RawSpectrumID` を書く順）。
  4. 名前が一致する参照レコードを store から引く。
  5. `match_spectrum` を回し、`[4] Simple` `[5] Weighted` `[6] Reverse`
     `[7] count` `[8] percentage` と突き合わせる。
  6. 相対誤差 1e-4 以内を一致とする。

**`spot.MatchResults` は「代表ファイル」1 個の測定スペクトルだけから作られる**
（`DataObjConverter.SetRepresentativeProperty` → `spot.MatchResults.MergeContainers(
representative.MatchResults)`。`representative` は `GetRepresentativeFileID` が選ぶ、
全ファイル中でトータルスコア最高・次点でピーク高最大の 1 ファイル）。
`spectra_ref` に列挙されるのは「そのファイル自身の**独立**トップ候補がこの SME 行と
同じ候補である」ファイル全部（`AddSpectraRef` 呼び出し側の条件）で、代表ファイルは
論理的にこの中に含まれるはずだが、どれが代表かは mzTab からは分からない。
そこで **列挙された全ファイルの実測スペクトルを順に試し、`[4][5][6]` が許容誤差内で
再現できたものを「代表ファイルの測定」とみなす**（brief 「近似で進めるな」は
spectra_ref の**解決方法自体**の話で、これは解決した候補の中から機械的な基準
（再現するかどうか）で選ぶだけなので該当しない）。

この run は脂質（Lipidomics）なので `[7]` `[8]`（matched peaks の count/percentage）は
上流の `GetLipidomicsMatchedPeaksScores`（クラス固有の診断イオン規則）で計算されており、
こちらが移植した汎用版とは値が合わない見込み（brief 記載の既知の限界）。`[4][5][6]` は
omics 分岐の外で計算されるので一致するはず。
"""
from __future__ import annotations

import argparse
import difflib
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from metabolomix.analysis import spectral_match
from metabolomix.dcl import reader as dcl_reader
from metabolomix.library import dbs as dbs_reader
from metabolomix.library import store as library_store
from metabolomix.mztab import reader as mztab_reader

_MS_RUN_LOCATION_RE = re.compile(r"^ms_run\[(\d+)\]-location$")
_SPECTRA_REF_ENTRY_RE = re.compile(
    r"^ms_run\[(\d+)\]:ms1scanID=(-?\d+)\s+ms2scanID=(-?\d+)$"
)

# id_confidence_measure のインデックス。MachineCategory.LCMS の並び
# （MztabFormatExport.cs ~1229 行目）: [1]=TotalScore [2]=RT similarity
# [3]=m/z similarity [4]=Simple dot product [5]=Weighted dot product
# [6]=Reverse dot product [7]=Matched peaks count [8]=Matched peaks percentage
_DOT_PRODUCT_INDICES = [4, 5, 6]  # -1/0 を同一視してよい 3 種
_MATCHED_PEAKS_INDICES = [7, 8]  # 素のフィールド。同一視しない
_METRIC_NAMES = {
    4: "simple_dot_product",
    5: "weighted_dot_product",
    6: "reverse_dot_product",
    7: "matched_peaks_count",
    8: "matched_peaks_percentage",
}
_REL_TOL = 1e-4


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="移植した match_spectrum を mzTab-M の id_confidence_measure[4..8] と突き合わせる。"
    )
    p.add_argument("--dbs", required=True, help="参照ライブラリ .dbs のパス")
    p.add_argument("--mztab", required=True, help="MS-DIAL が出力した mzTab-M ファイルのパス")
    p.add_argument("--dcl-dir", required=True, help=".dcl（測定 MS/MS）が並ぶディレクトリ")
    p.add_argument("--limit", type=int, default=200, help="比較する SME 行数の上限（既定 200）")
    return p.parse_args()


def _strip_file_uri(location: str) -> str:
    if location.startswith("file://"):
        return location[len("file://"):]
    return location


def build_ms_run_to_wiff_stem(metadata: dict) -> dict[int, str]:
    """`ms_run[N]-location` から N -> 生データの basename（拡張子なし）を作る。"""
    mapping = {}
    for key, value in metadata.items():
        m = _MS_RUN_LOCATION_RE.match(key)
        if not m:
            continue
        n = int(m.group(1))
        path = _strip_file_uri(value)
        mapping[n] = Path(path).stem
    return mapping


def find_dcl_for_stem(dcl_dir: Path, stem: str) -> Path | None:
    """`<stem>_<timestamp>.dcl` を dcl_dir から探す。複数あれば最新（文字列最大）を使う。"""
    candidates = sorted(dcl_dir.glob(f"{stem}_*.dcl"))
    if not candidates:
        return None
    return candidates[-1]


def parse_spectra_ref(spectra_ref: str | None) -> list[tuple[int, int, int]]:
    """`ms_run[N]:ms1scanID=X ms2scanID=Y` の `| ` 区切りリストを
    `[(N, ms1scanID, ms2scanID), ...]` に分解する。`null`/空なら空リスト。"""
    if not spectra_ref or spectra_ref == "null":
        return []
    entries = []
    for token in spectra_ref.split("|"):
        token = token.strip()
        m = _SPECTRA_REF_ENTRY_RE.match(token)
        if not m:
            raise ValueError(f"spectra_ref のエントリを解釈できません: {token!r}")
        entries.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    return entries


class DclCache:
    """`.dcl` ファイルをファイルパス単位で読み直さないためのキャッシュ。
    `raw_spec_id -> [entry, ...]`（重複はまず無いはずだが理論上あり得るのでリストで持つ）。"""

    def __init__(self, dcl_dir: Path, ms_run_to_stem: dict[int, str]):
        self._dcl_dir = dcl_dir
        self._ms_run_to_stem = ms_run_to_stem
        self._by_ms_run: dict[int, dict[int, list[dict]] | None] = {}
        self.missing_ms_runs: set[int] = set()

    def _load(self, ms_run: int) -> dict[int, list[dict]] | None:
        if ms_run in self._by_ms_run:
            return self._by_ms_run[ms_run]
        stem = self._ms_run_to_stem.get(ms_run)
        if stem is None:
            self._by_ms_run[ms_run] = None
            self.missing_ms_runs.add(ms_run)
            return None
        dcl_path = find_dcl_for_stem(self._dcl_dir, stem)
        if dcl_path is None:
            self._by_ms_run[ms_run] = None
            self.missing_ms_runs.add(ms_run)
            return None
        results = dcl_reader.deserialize_dcl(str(dcl_path), include_spectrum=True, top_n_peaks=None)
        index: dict[int, list[dict]] = {}
        for r in results:
            index.setdefault(r["raw_spec_id"], []).append(r)
        self._by_ms_run[ms_run] = index
        return index

    def lookup(self, ms_run: int, ms2_scan_id: int) -> list[dict]:
        index = self._load(ms_run)
        if index is None:
            return []
        return index.get(ms2_scan_id, [])


def resolve_references(store, exp_mz: float, ion_mode: str, chemical_name: str, mz_tol: float):
    """名前 + precursor m/z 窓で参照レコード候補を集める。`([record, ...], reason)` を返す。

    見つからなければ `([], reason)`。

    **同名でも異なるスペクトルを持つ参照レコードが実在する**（例: `ST 24:1;O5` /
    `[M-H]-` に 26 本・8 本・5 本ピークの 3 レコードが同居していた。恐らく CE 違い/
    情報源違いの in-silico 予測スペクトル）。この場合どれが実際に使われたかは
    名前だけでは決まらないので、**「一致」ではなく「候補」を返し**、呼び出し側が
    測定側の候補（`spectra_ref` の各エントリ）と総当たりして `id_confidence_measure`
    を再現できる組を探す（brief の「近似で進めるな」は解決方法自体の話で、これは
    実データという正解と突き合わせて機械的に選ぶだけなので該当しない）。
    """
    candidates = store.candidates(exp_mz, mz_tol=mz_tol, ion_mode=ion_mode)
    if not candidates:
        return [], "no_candidates_in_window"

    # 名前の一致は **`O-`/`P-` の表記ゆらぎを吸収して**見る（`ether_name_variants`）。
    # 素朴な文字列一致だと、MS-DIAL が採点した `P-` のレコードではなく同じ m/z 窓に
    # 居る別レコード（`O-` の名前を持つ、スペクトルも違うもの）を掴む。
    wanted = ether_name_variants(chemical_name)
    exact = [c for c in candidates if c["name"] in wanted]
    if len(exact) == 1:
        return exact, "exact"
    if len(exact) > 1:
        # 同義名で複数当たるのは想定内（`O-` 側と `P-` 側が両方ライブラリに居る）。
        # どちらが採点されたかは名前では決まらないので、呼び出し側が測定側と
        # 総当たりして再現する組を選ぶ。
        return exact, "exact_ambiguous"

    if len(candidates) == 1:
        return candidates, "single_candidate_fallback"

    # 名前が一致しない（脂質の sn 位置精緻化などで書き換わっている可能性）。
    # m/z 窓の候補から名前が最も近いものを選ぶ（brief 「代替を使ってよい」）。
    # ここは正解と突き合わせようがない（本当に名前が引けない）ので 1 件に絞り、
    # 呼び出し側でこの区分だけ headline 統計から外して別集計する。
    best = max(
        candidates,
        key=lambda c: difflib.SequenceMatcher(None, c["name"] or "", chemical_name).ratio(),
    )
    return [best], "fuzzy_fallback"


_ETHER_CHAIN_RE = re.compile("([OP])-([0-9]+):([0-9]+)")


def ether_name_variants(name: str) -> set[str]:
    """エーテル脂質の `O-` / `P-` 表記ゆらぎを吸収した同義名の集合を返す。

    `O-n:m`（プラスマニル＝1-O-alkyl）と `P-n:(m-1)`（プラスマローゲン＝
    1-O-alk-1'-enyl）は**同じ化学種の別表記**で、ライブラリには両方の名前の
    レコードが別々に入っている（スペクトルも別物）。**MS-DIAL は `P-` の
    レコードを採点し、`GetRefinedLipidAnnotationLevel` を通して `O-` 表記で
    mzTab へ書き出す**ので、mzTab の名前で素朴に引くと採点されたのとは
    別のレコードに当たる（HISTRY 2026-09-22(5)。これを落とすと `EtherLPE` の
    5 件が「移植の誤り」に化けて見える）。

    鎖が複数ある名前（`PE O-16:0_18:1`）では `O-`/`P-` が付いた鎖だけ書き換える。
    二重結合数が 0 の `O-n:0` は相方が `P-n:-1` になるので作らない。
    """
    variants = {name}
    if not name:
        return variants

    def _swap(match: re.Match) -> str:
        kind, carbon, double_bonds = match.group(1), match.group(2), int(match.group(3))
        if kind == "O":
            return match.group(0) if double_bonds == 0 else f"P-{carbon}:{double_bonds - 1}"
        return f"O-{carbon}:{double_bonds + 1}"

    swapped = _ETHER_CHAIN_RE.sub(_swap, name)
    if swapped != name:
        variants.add(swapped)
    return variants


def _clamp_sentinel(value: float) -> float:
    return 0.0 if value == -1.0 else value


def _relative_error(ours: float, reported: float) -> float:
    if abs(reported) < 1e-9 and abs(ours) < 1e-9:
        return 0.0
    denom = max(abs(reported), abs(ours), 1e-12)
    return abs(reported - ours) / denom


def compare_row_against_measured(computed: dict, row: dict) -> dict:
    """`computed`（`match_spectrum` の戻り値）と SME 行の `id_confidence_measure[N]` を
    突き合わせ、指標ごとの一致可否・相対誤差を返す。"""
    result = {}
    ours_by_index = {
        4: computed["simple_dot_product"],
        5: computed["weighted_dot_product"],
        6: computed["reverse_dot_product"],
        7: computed["matched_peaks_count"],
        8: computed["matched_peaks_percentage"],
    }
    for idx in _DOT_PRODUCT_INDICES + _MATCHED_PEAKS_INDICES:
        reported_raw = row.get(f"id_confidence_measure[{idx}]")
        if reported_raw is None:
            result[idx] = None
            continue
        reported = float(reported_raw)
        ours = ours_by_index[idx]
        ours_cmp = _clamp_sentinel(ours) if idx in _DOT_PRODUCT_INDICES else ours
        rel_err = _relative_error(ours_cmp, reported)
        result[idx] = {
            "ours_raw": ours,
            "ours_compared": ours_cmp,
            "reported": reported,
            "rel_err": rel_err,
            "match": rel_err <= _REL_TOL,
        }
    return result


def main() -> int:
    args = parse_args()
    dbs_path = Path(args.dbs)
    mztab_path = Path(args.mztab)
    dcl_dir = Path(args.dcl_dir)

    print(f"[INFO] .dbs を読み込み中: {dbs_path}")
    store = library_store.open_store(dbs_path)
    search_params = dbs_reader.read_storage_meta(dbs_path)["search_params"]
    ms2_tol = search_params["ms2_tolerance"]
    ms1_tol = search_params["ms1_tolerance"]
    mass_begin = search_params["mass_range_begin"]
    mass_end = search_params["mass_range_end"]
    relative_amp_cutoff = search_params["relative_amp_cutoff"]
    absolute_amp_cutoff = search_params["absolute_amp_cutoff"]
    print(
        f"[INFO] search_params: ms1_tol={ms1_tol} ms2_tol={ms2_tol} "
        f"mass_begin={mass_begin} mass_end={mass_end} "
        f"relative_amp_cutoff={relative_amp_cutoff} absolute_amp_cutoff={absolute_amp_cutoff}"
    )
    print(f"[INFO] store: {store.summary()}")

    print(f"[INFO] mzTab を読み込み中: {mztab_path}")
    parsed = mztab_reader.parse_mztab(mztab_path)
    ms_run_to_stem = build_ms_run_to_wiff_stem(parsed["metadata"])
    print(f"[INFO] ms_run 数: {len(ms_run_to_stem)}")

    sme_rows = parsed["sections"].get("SME", {}).get("rows", [])
    print(f"[INFO] SME 行数（全体）: {len(sme_rows)} / 比較対象上限: {args.limit}")

    dcl_cache = DclCache(dcl_dir, ms_run_to_stem)

    # 「名前で自信を持って引けた」行（exact / exact_ambiguous）を headline 統計に、
    # 「名前で引けず代替に頼った」行（fuzzy_fallback / single_candidate_fallback）は
    # 参照レコード自体が正しい保証が無いので別集計にする（brief 「代替を使ってよいが
    # 何をしたか書け」に対応）。
    _CONFIDENT_REASONS = {"exact", "exact_ambiguous"}

    rows_compared = 0
    rows_skipped: Counter[str] = Counter()
    resolution_reasons: Counter[str] = Counter()
    metric_match_counts: Counter[int] = Counter()
    metric_total_counts: Counter[int] = Counter()
    fallback_metric_match_counts: Counter[int] = Counter()
    fallback_metric_total_counts: Counter[int] = Counter()
    mismatches: list[dict] = []
    fallback_mismatches: list[dict] = []

    for row in sme_rows[: args.limit]:
        chemical_name = row.get("chemical_name")
        exp_mz_raw = row.get("exp_mass_to_charge")
        charge_raw = row.get("charge")
        spectra_ref_raw = row.get("spectra_ref")

        if chemical_name is None or exp_mz_raw is None:
            rows_skipped["missing_name_or_mz"] += 1
            continue

        try:
            exp_mz = float(exp_mz_raw)
        except ValueError:
            rows_skipped["bad_exp_mz"] += 1
            continue

        ion_mode = "negative" if (charge_raw or "").strip().startswith("-") else "positive"

        try:
            entries = parse_spectra_ref(spectra_ref_raw)
        except ValueError as exc:
            rows_skipped["unparsable_spectra_ref"] += 1
            print(f"[WARN] SME {row.get('SME_ID')}: {exc}")
            continue

        if not entries:
            rows_skipped["no_spectra_ref"] += 1
            continue

        references, reason = resolve_references(store, exp_mz, ion_mode, chemical_name, ms1_tol)
        resolution_reasons[reason] += 1
        if not references:
            rows_skipped[f"reference_not_found:{reason}"] += 1
            continue

        # spectra_ref の各測定 × 参照候補（同名でもスペクトルが違う重複がある）を
        # 総当たりし、[4][5][6] を再現できる組を探す。見つからなければ、[4][5][6]の
        # 誤差合計が最小の組を診断用に残す。
        best_debug = None
        best_error_sum = None
        matched_debug = None
        for ms_run, ms1_scan_id, ms2_scan_id in entries:
            dcl_hits = dcl_cache.lookup(ms_run, ms2_scan_id)
            if not dcl_hits:
                continue
            # 同じ raw_spec_id が複数あれば precursor m/z が近いものを選ぶ。
            dcl_entry = min(dcl_hits, key=lambda r: abs(r["precursor_mz"] - exp_mz))
            mz_delta = abs(dcl_entry["precursor_mz"] - exp_mz)

            for reference in references:
                computed = spectral_match.match_spectrum(
                    dcl_entry["msms_spectrum"],
                    reference["spectrum"],
                    ms2_tol=ms2_tol,
                    mass_begin=mass_begin,
                    mass_end=mass_end,
                    relative_amp_cutoff=relative_amp_cutoff,
                    absolute_amp_cutoff=absolute_amp_cutoff,
                )
                comparison = compare_row_against_measured(computed, row)
                dot_product_ok = all(
                    comparison[idx] is not None and comparison[idx]["match"]
                    for idx in _DOT_PRODUCT_INDICES
                )
                error_sum = sum(
                    comparison[idx]["rel_err"] for idx in _DOT_PRODUCT_INDICES
                    if comparison[idx] is not None
                )
                debug = {
                    "ms_run": ms_run,
                    "ms2_scan_id": ms2_scan_id,
                    "mz_delta": mz_delta,
                    "comparison": comparison,
                    "reference": reference,
                }
                if best_error_sum is None or error_sum < best_error_sum:
                    best_error_sum = error_sum
                    best_debug = debug
                if dot_product_ok:
                    matched_debug = debug
                    break
            if matched_debug is not None:
                break

        if best_debug is None:
            rows_skipped["no_dcl_hit_for_any_spectra_ref_entry"] += 1
            continue

        rows_compared += 1
        confident = reason in _CONFIDENT_REASONS
        used = matched_debug or best_debug
        match_counts = metric_match_counts if confident else fallback_metric_match_counts
        total_counts = metric_total_counts if confident else fallback_metric_total_counts
        target_mismatches = mismatches if confident else fallback_mismatches
        for idx in _DOT_PRODUCT_INDICES + _MATCHED_PEAKS_INDICES:
            entry = used["comparison"].get(idx)
            if entry is None:
                continue
            total_counts[idx] += 1
            if entry["match"]:
                match_counts[idx] += 1
            else:
                target_mismatches.append(
                    {
                        "sme_id": row.get("SME_ID"),
                        "chemical_name": chemical_name,
                        "metric": _METRIC_NAMES[idx],
                        "ours": entry["ours_compared"],
                        "reported": entry["reported"],
                        "rel_err": entry["rel_err"],
                        "compound_class": used["reference"].get("compound_class"),
                        "ms_run": used["ms_run"],
                        "mz_delta": used["mz_delta"],
                        "resolution_reason": reason,
                        "used_best_effort": matched_debug is None,
                    }
                )

    def print_mismatches(title, mism):
        if not mism:
            return
        print()
        print(title)
        dot_mismatches = [m for m in mism if m["metric"] in (
            "simple_dot_product", "weighted_dot_product", "reverse_dot_product")]
        other_mismatches = [m for m in mism if m not in dot_mismatches]
        ordered = sorted(dot_mismatches, key=lambda m: -m["rel_err"]) + \
            sorted(other_mismatches, key=lambda m: -m["rel_err"])
        for m in ordered[:10]:
            print(
                f"  SME={m['sme_id']} name={m['chemical_name']!r} class={m['compound_class']!r} "
                f"metric={m['metric']} ours={m['ours']:.6g} reported={m['reported']:.6g} "
                f"rel_err={m['rel_err']:.3g} ms_run={m['ms_run']} resolved_by={m['resolution_reason']}"
            )

    print()
    print("===== 結果 =====")
    print(f"比較した SME 行数: {rows_compared}")
    print(f"スキップした行の内訳: {dict(rows_skipped)}")
    print(f"参照レコード解決方法の内訳: {dict(resolution_reasons)}")
    if dcl_cache.missing_ms_runs:
        print(f"[WARN] .dcl が見つからなかった ms_run: {sorted(dcl_cache.missing_ms_runs)}")

    print()
    print("[headline] 名前で自信を持って解決できた行（exact / exact_ambiguous）の一致率:")
    print("指標ごとの一致率（相対誤差 1e-4 以内）:")
    for idx in _DOT_PRODUCT_INDICES + _MATCHED_PEAKS_INDICES:
        total = metric_total_counts[idx]
        matched = metric_match_counts[idx]
        rate = (matched / total * 100.0) if total else float("nan")
        print(f"  [{idx}] {_METRIC_NAMES[idx]:<28} {matched}/{total} ({rate:.1f}%)")
    print_mismatches("不一致（上位 10 件、[4][5][6] 優先、headline）:", mismatches)

    fallback_total_rows = sum(1 for r in resolution_reasons
                               if r in ("fuzzy_fallback", "single_candidate_fallback")
                               for _ in range(resolution_reasons[r]))
    if fallback_total_rows:
        print()
        print(
            "[参考] 名前で引けず代替（fuzzy_fallback / single_candidate_fallback）に "
            "頼った行 — 参照レコード自体の正しさを保証できないので headline から除外:"
        )
        for idx in _DOT_PRODUCT_INDICES + _MATCHED_PEAKS_INDICES:
            total = fallback_metric_total_counts[idx]
            matched = fallback_metric_match_counts[idx]
            rate = (matched / total * 100.0) if total else float("nan")
            print(f"  [{idx}] {_METRIC_NAMES[idx]:<28} {matched}/{total} ({rate:.1f}%)")
        print_mismatches("不一致（上位 10 件、[4][5][6] 優先、参考）:", fallback_mismatches)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
