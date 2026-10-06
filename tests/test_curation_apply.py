from metabolomix.curation import apply, flags


def _flag(tmp_path, spot_id, flag):
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    if not arf2.exists():
        arf2.write_bytes(b"x")
    flags.FlagStore(flags.curation_dir(arf2)).append(
        [{"spot_id": spot_id, "flag": flag}], alignment=flags.alignment_key(arf2),
        review_id="r", source="user")
    return arf2


def test_no_flags_means_no_meta_line(tmp_path):
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")
    state = apply.flags_for_arf2(arf2)
    assert state["n"] == 0
    assert apply.meta_line("applied", state, {"wrong_excluded": 0, "suspect": 0}) is None


def test_wrong_rows_are_dropped_and_suspect_kept(tmp_path):
    _flag(tmp_path, 1, "wrong")
    arf2 = _flag(tmp_path, 2, "suspect")
    state = apply.flags_for_arf2(arf2)
    rows = [{"spot_id": 1}, {"spot_id": 2}, {"spot_id": 3}]
    kept, stats = apply.filter_rows(rows, state, key=lambda r: r["spot_id"])
    assert [r["spot_id"] for r in kept] == [2, 3]
    assert stats == {"wrong_excluded": 1, "redundant_excluded": 0, "suspect": 1, "assigned": 0}
    line = apply.meta_line("applied", state, stats)
    assert line.startswith("# curation = applied\t")
    assert "curation_wrong_excluded = 1" in line
    assert f"curation_flags_sha256 = {state['digest']}" in line


def test_not_applied_is_declared(tmp_path):
    arf2 = _flag(tmp_path, 1, "wrong")
    line = apply.meta_line("not_applied", apply.flags_for_arf2(arf2), None)
    assert line.startswith("# curation = not_applied\t")


def test_arf2_for_mztab_matches_the_batch_stem(tmp_path):
    (tmp_path / "AlignmentResult_2026_01_01.arf2").write_bytes(b"x")
    (tmp_path / "AlignmentResult_2026_02_02.arf2").write_bytes(b"x")
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2026_01_01_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    assert apply.arf2_for_mztab(mztab).name == "AlignmentResult_2026_01_01.arf2"
    other = tmp_path / "unrelated.mzTab"
    other.write_text("", encoding="utf-8")
    assert apply.arf2_for_mztab(other) is None


# ---------- export_dataset_result への反映（Step 5） ----------

def _prepared_dataset():
    from metabolomix.analysis.dataset_service import compare_dataset, preprocess_dataset
    from tests.pipeline_fixtures import make_dataset
    ds = make_dataset()
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})
    result = compare_dataset(ds, ds.sample_names[:4], ds.sample_names[4:])
    return ds, result


def test_export_dataset_result_with_no_curation_has_no_meta_line(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()

    info = export_dataset_result(ds, result, tmp_path / "out.tsv", curation=None)

    text = (tmp_path / "out.tsv").read_text(encoding="utf-8")
    assert not any(l.startswith("# curation") for l in text.splitlines())
    assert info["curation"] is None


def test_export_dataset_result_applied_curation_drops_the_wrong_row(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    curation = {"state": "applied",
               "flag_set": {"wrong": {"F0"}, "suspect": set(),
                            "digest": "d", "n": 1}}

    info = export_dataset_result(ds, result, tmp_path / "out.tsv", curation=curation)

    text = (tmp_path / "out.tsv").read_text(encoding="utf-8")
    body = [l for l in text.splitlines() if l and not l.startswith("#")][1:]
    assert not any(l.split("\t")[0] == "F0" for l in body)
    meta = [l for l in text.splitlines() if l.startswith("# curation")]
    assert len(meta) == 1
    assert meta[0].startswith("# curation = applied\t")
    assert "curation_wrong_excluded = 1" in meta[0]
    assert info["curation"] == "applied"


def test_export_dataset_result_unmapped_curation_keeps_all_rows(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    curation = {"state": "unmapped",
               "flag_set": {"wrong": {"F0"}, "suspect": set(),
                            "digest": "d", "n": 1}}

    info = export_dataset_result(ds, result, tmp_path / "out.tsv", curation=curation)

    text = (tmp_path / "out.tsv").read_text(encoding="utf-8")
    body = [l for l in text.splitlines() if l and not l.startswith("#")][1:]
    assert any(l.split("\t")[0] == "F0" for l in body)
    meta = [l for l in text.splitlines() if l.startswith("# curation")]
    assert len(meta) == 1
    assert meta[0].startswith("# curation = unmapped\t")
    assert info["curation"] == "unmapped"


# ---------- I7: 以前の版のアラインメントに付いたフラグ(orphaned) ----------

def test_flags_recorded_against_an_earlier_version_are_counted_as_orphaned(tmp_path):
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"new")
    store = flags.FlagStore(flags.curation_dir(arf2))
    old = {"alignment_file": arf2.name, "alignment_sha256": "0" * 64}
    store.append([{"spot_id": 1, "flag": "wrong"}, {"spot_id": 2, "flag": "suspect"},
                  {"spot_id": 3, "flag": "wrong"}], alignment=old, review_id="r", source="user")
    store.append([{"spot_id": 3, "flag": "clear"}], alignment=old, review_id="r", source="user")
    other_file = {"alignment_file": "AlignmentResult_other.arf2", "alignment_sha256": "1" * 64}
    store.append([{"spot_id": 9, "flag": "wrong"}], alignment=other_file, review_id="r", source="user")
    state = apply.flags_for_arf2(arf2)
    assert state["n"] == 0 and state["wrong"] == set()
    assert state["orphaned"] == 2


def test_no_flags_file_skips_hashing_the_alignment(tmp_path, monkeypatch):
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")

    def boom(_):
        raise AssertionError("フラグが無いのに .arf2 を hash した")

    monkeypatch.setattr(apply, "alignment_key", boom)
    state = apply.flags_for_arf2(arf2)
    assert state["n"] == 0 and state["orphaned"] == 0


# ---------- M6: バッチ語幹の後ろに数字が続く別バッチを拾わない ----------

def test_arf2_for_mztab_does_not_match_a_stem_that_is_a_prefix_of_another(tmp_path):
    (tmp_path / "AlignmentResult_2026_01_01_2_3.arf2").write_bytes(b"x")
    (tmp_path / "AlignmentResult_2026_01_01_2_30.arf2").write_bytes(b"x")
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2_30_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    assert apply.arf2_for_mztab(mztab).name == "AlignmentResult_2026_01_01_2_30.arf2"
    short = tmp_path / "Height_AlignmentResult_2026_01_01_2_3_09.mzTab"
    short.write_text("", encoding="utf-8")
    assert apply.arf2_for_mztab(short).name == "AlignmentResult_2026_01_01_2_3.arf2"


# ---------- assign / redundant（候補付けの判断の反映） ----------

from metabolomix.curation import apply as curation_apply  # noqa: E402


def _flag_set(**kw):
    base = {"wrong": set(), "suspect": set(), "assign": {}, "redundant": set(), "digest": "d", "n": 0,
            "orphaned": 0}
    base.update(kw)
    base["n"] = len(base["wrong"]) + len(base["suspect"]) + len(base["assign"]) + len(base["redundant"])
    return base


def test_meta_line_unchanged_without_decisions():
    fs = _flag_set(wrong={1}, suspect={2})
    line = curation_apply.meta_line("applied", fs, {"wrong_excluded": 1, "suspect": 1,
                                                    "redundant_excluded": 0, "assigned": 0})
    assert line == ("# curation = applied\tcuration_flags = 2\tcuration_wrong_excluded = 1\t"
                    "curation_suspect = 1\tcuration_flags_sha256 = d")


def test_meta_line_with_decisions():
    fs = _flag_set(assign={3: {"name": "PC 34:1", "ontology": "PC", "inchikey": "K"}}, redundant={4})
    line = curation_apply.meta_line("applied", fs, {"wrong_excluded": 0, "suspect": 0,
                                                    "redundant_excluded": 1, "assigned": 1})
    assert "curation_assigned = 1\tcuration_redundant_excluded = 1\tcuration_flags_sha256 = d" in line


def test_override_identity_makes_an_unannotated_spot_exportable():
    catalog = {3: {"MasterAlignmentID": 3, "Name": "Unknown", "Ontology": "", "InChIKey": ""}}
    fs = _flag_set(assign={3: {"name": "PE 36:2", "ontology": "PE", "inchikey": "KEY-PE362"}})
    out = curation_apply.override_identity(catalog, fs)
    assert out[3]["Name"] == "PE 36:2" and out[3]["InChIKey"] == "KEY-PE362" and out[3]["_curation"] == "assign"
    assert catalog[3]["Name"] == "Unknown"                       # 元の dict は変えない


def test_filter_rows_excludes_wrong_and_redundant_and_counts_assigned():
    rows = [{"spot_id": i} for i in range(5)]
    fs = _flag_set(wrong={0}, redundant={1}, suspect={2}, assign={3: {"name": "x", "ontology": "", "inchikey": "k"}})
    kept, stats = curation_apply.filter_rows(rows, fs, key=lambda r: r["spot_id"])
    assert [r["spot_id"] for r in kept] == [2, 3, 4]
    assert stats == {"wrong_excluded": 1, "redundant_excluded": 1, "suspect": 1, "assigned": 1}
