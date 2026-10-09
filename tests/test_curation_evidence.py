import pytest

from metabolomix.arf2.reader import load_catalog
from metabolomix.curation import evidence, judge
from metabolomix.library import store as library_store
from tests.curation_fixtures import write_alignment_set


@pytest.fixture()
def dataset(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg")
    s = library_store.open_store(paths["msp"])
    yield paths, s
    s.close()


def test_sibling_files_share_the_alignment_stem(dataset):
    paths, _ = dataset
    files = evidence.sibling_files(paths["arf2"])
    assert files["dcl"].name == "AlignmentResult_x.dcl"
    assert files["eic"].name == "AlignmentResult_x.EIC.aef"
    assert files["arf"].name == "AlignmentResult_x_PeakProperties.arf"


def test_select_spots_skips_unknown_and_filters(dataset):
    paths, _ = dataset
    catalog = load_catalog(paths["arf2"])
    assert [s["MasterAlignmentID"] for s in evidence.select_spots(catalog, ontology=None, name_contains=None)] == [0, 1]
    assert [s["MasterAlignmentID"] for s in evidence.select_spots(catalog, ontology=None, name_contains="36:2")] == [1]
    assert evidence.select_spots(catalog, ontology=["PE"], name_contains=None) == []


def test_collect_builds_the_evidence_contract(dataset):
    paths, s = dataset
    catalog = load_catalog(paths["arf2"])
    spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
    th = judge.resolve_thresholds(None)
    evs, stats = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025, th=th)
    by_id = {e["spot_id"]: e for e in evs}
    first = by_id[0]
    assert first["reference"]["name"] == "PC 34:1"
    assert first["ppm_basis"] == "reference"
    assert abs(first["ppm"]) < 1.0
    assert first["drt"] == pytest.approx(0.0, abs=1e-6)
    assert first["mirror"]["measured"] and first["mirror"]["reference"]
    assert first["rescore"]["weighted_dot_product"] > 0.5
    assert [smp["detected"] for smp in first["eic"]["samples"]] == [True, True, False]
    assert first["eic"]["samples"][0]["representative"] is True
    assert first["eic_shape"]["band"] == "PASS"
    second = by_id[1]
    assert second["name_prefix"] == "low score"
    assert second["drt"] == pytest.approx(1.5, abs=1e-3)
    assert stats["n_reference_resolved"] == 2


def test_collect_attaches_the_measured_and_theoretical_ms1_isotopes(dataset):
    paths, s = dataset
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    evs, _ = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                              th=judge.resolve_thresholds(None))
    iso = {e["spot_id"]: e for e in evs}[0]["isotopes"]
    assert [p[1] for p in iso["measured"]] == [100.0, 46.0, 12.0]
    assert [p[2] for p in iso["measured"]] == [10000, 4600, 1200]
    assert iso["measured"][0][0] == pytest.approx(760.5851, abs=1e-4)
    assert iso["theoretical"]["basis"] == "formula+adduct"
    assert iso["theoretical"]["relative"][1] == pytest.approx(47.1, abs=0.5)


def test_eic_points_are_trimmed_around_the_peak(dataset):
    paths, s = dataset
    catalog = load_catalog(paths["arf2"])
    spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
    evs, _ = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                              th=judge.resolve_thresholds(None))
    sample = evs[0]["eic"]["samples"][0]
    xs = [p[0] for p in sample["points"]]
    width = sample["right"] - sample["left"]
    assert min(xs) >= sample["left"] - 1.5 * width - 1e-6
    assert max(xs) <= sample["right"] + 1.5 * width + 1e-6


def test_short_eic_trace_is_unchanged_apart_from_rounding_intensity():
    points = [[0.0, 1.4], [0.01, 2.6], [0.02, 3.0]]
    out = evidence._downsample_points(points, 0.0, 0.02)
    assert out == [[0.0, 1], [0.01, 3], [0.02, 3]]


def test_long_eic_trace_is_cut_to_the_cap_and_keeps_the_apex():
    points = [[round(0.01 * i, 3), float(i)] for i in range(200)]
    points[150] = [points[150][0], 9000.0]     # 頂点を端から離れた位置に置く
    out = evidence._downsample_points(points, points[0][0], points[-1][0])
    assert len(out) <= evidence.EIC_MAX_POINTS
    xs = [p[0] for p in out]
    assert points[0][0] in xs and points[-1][0] in xs
    assert points[150][0] in xs
    assert all(isinstance(p[1], int) for p in out)


def test_collect_computes_eic_shape_on_the_full_series_but_downsamples_the_payload(
        tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg", eic_points=61)
    s = library_store.open_store(paths["msp"])
    captured = {}
    original_spot_shape = evidence.spot_shape

    def spy(samples, th):
        captured.setdefault("lengths", []).extend(
            len(sample["chromatogram"]) for sample in samples)
        return original_spot_shape(samples, th)

    monkeypatch.setattr(evidence, "spot_shape", spy)
    try:
        catalog = load_catalog(paths["arf2"])
        spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
        evs, _ = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                  th=judge.resolve_thresholds(None))
    finally:
        s.close()
    assert max(captured["lengths"]) > evidence.EIC_MAX_POINTS   # 形状は間引き前の全点で計算
    first = evs[0]
    assert first["eic_shape"]["band"] == "PASS"
    for sample in first["eic"]["samples"]:
        assert len(sample["points"]) <= evidence.EIC_MAX_POINTS   # payload は間引き後


def test_mirror_payload_is_cut_to_the_cap_keeping_all_matched_peaks():
    measured = [[100.0 + i * 0.01, 1.0 + (i % 7)] for i in range(400)]
    matched_measured_mz = [measured[10][0], measured[250][0], measured[390][0]]
    mirror = {
        "measured": measured,
        "reference": [[100.1, 500.0], [101.0, 10.0]],
        "matched_mz": [100.1],
        "matched_measured_mz": matched_measured_mz,
        "unscored_mz": [], "scored_peak_count": len(measured), "unscored_peak_count": 0,
        "labels": [],
    }
    cut = evidence._cut_mirror_for_payload(mirror)
    assert len(cut["measured"]) <= evidence.MIRROR_MAX_PEAKS
    kept_mz = {p[0] for p in cut["measured"]}
    assert set(matched_measured_mz) <= kept_mz
    # 判定に使うフィールドは間引きの影響を受けない(満スペクトルの値のまま)
    assert cut["scored_peak_count"] == len(measured)
    assert cut["matched_mz"] == [100.1]


def test_missing_reference_is_reported_not_raised(dataset, tmp_path):
    paths, _ = dataset
    other = tmp_path / "other.msp"
    other.write_text("NAME: X\nPRECURSORMZ: 100.0\nIONMODE: Positive\nNum Peaks: 0\n", encoding="utf-8")
    s = library_store.open_store(other)
    try:
        catalog = load_catalog(paths["arf2"])
        spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
        evs, stats = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                      th=judge.resolve_thresholds(None))
        assert all(e["reference"] is None for e in evs)
        assert evs[0]["ppm_basis"] == "formula"
        assert stats["n_reference_resolved"] == 0
    finally:
        s.close()


# ---------- C1: 間引きは全区間に等間隔で散らす ----------

@pytest.mark.parametrize("n", [41, 60, 79, 80, 300, 529])
def test_downsample_spreads_points_across_the_whole_trace(n):
    import math
    points = [[round(0.01 * i, 3), 10.0 + (i % 5)] for i in range(n)]
    apex = (3 * n) // 4                        # 頂点を右半分に置く(旧実装は右半分を落とした)
    points[apex] = [points[apex][0], 1e6]
    left_i, right_i = n // 3, (2 * n) // 3 + 1
    left, right = points[left_i][0] + 0.001, points[right_i][0] - 0.001
    out = evidence._downsample_points(points, left, right)
    assert len(out) <= evidence.EIC_MAX_POINTS
    index_of = {p[0]: i for i, p in enumerate(points)}
    kept = [index_of[p[0]] for p in out]
    assert kept == sorted(kept)
    for anchor in (0, n - 1, apex, left_i, right_i):
        assert anchor in kept, (n, anchor)
    bound = math.ceil(2 * n / evidence.EIC_MAX_POINTS)
    gaps = [b - a for a, b in zip(kept, kept[1:])]
    assert max(gaps) <= bound, (n, max(gaps), bound)


def test_downsample_leaves_traces_at_or_below_the_cap_unchanged_apart_from_rounding():
    points = [[round(0.01 * i, 3), i + 0.4] for i in range(evidence.EIC_MAX_POINTS)]
    out = evidence._downsample_points(points, 0.1, 0.2)
    assert out == [[p[0], round(p[1])] for p in points]


# ---------- M3: file_id の無い行を飛ばす / 未知の file_ids を先に弾く ----------

def test_choose_file_ids_skips_rows_without_a_file_id():
    rows = [{}, {"file_id": 1, "height": 5.0, "is_gap_filled": False}, {"height": 9.0}]
    assert evidence._choose_file_ids(rows, None, None, 12) == [1]


def test_collect_rejects_file_ids_absent_from_the_alignment(dataset):
    paths, s = dataset
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    with pytest.raises(evidence.UnknownFileIdsError) as info:
        evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                         th=judge.resolve_thresholds(None), file_ids=[0, 98, 99])
    assert info.value.missing == [98, 99]


# ---------- M10: 兄弟ファイルが欠けたら missing_files に出る ----------

def test_collect_reports_missing_eic_and_dcl_siblings(dataset):
    paths, s = dataset
    folder = paths["arf2"].parent
    (folder / "AlignmentResult_x.EIC.aef").unlink()
    (folder / "AlignmentResult_x.dcl").unlink()
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    evs, stats = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                  th=judge.resolve_thresholds(None))
    assert stats["missing_files"] == ["dcl", "eic"]
    assert all(e["eic"]["samples"] == [] and e["mirror"] is None for e in evs)


def test_collect_reports_the_mz_difference_in_mda(dataset):
    paths, s = dataset
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    evs, _ = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                              th=judge.resolve_thresholds(None))
    first = {e["spot_id"]: e for e in evs}[0]
    expected = (first["rep_mz"] - first["reference"]["precursor_mz"]) * 1000
    assert first["dmz_mda"] == pytest.approx(expected, abs=0.01)


def test_mz_difference_falls_back_to_the_formula(dataset, tmp_path):
    paths, _ = dataset
    other = tmp_path / "other.msp"
    other.write_text("NAME: X\nPRECURSORMZ: 100.0\nIONMODE: Positive\nNum Peaks: 0\n", encoding="utf-8")
    s = library_store.open_store(other)
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        evs, _ = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                  th=judge.resolve_thresholds(None))
    finally:
        s.close()
    assert evs[0]["ppm_basis"] == "formula"
    assert evs[0]["dmz_mda"] is not None
    assert evs[0]["dmz_mda"] == pytest.approx(evs[0]["ppm"] * evs[0]["rep_mz"] / 1000, abs=0.05)
