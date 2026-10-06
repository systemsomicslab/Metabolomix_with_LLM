# tests/test_dataset_state.py
import math
import textwrap
import pytest
from metabolomix.mztab.reader import parse_mztab
from metabolomix.mztab.dataset_state import DatasetState, build_dataset_state
from metabolomix.core import session_state

# SMF は同定を持たず SME_ID_REFS で SME を指す（mzTab-M 2.0.0-M の実形状）。
# 構造・名称を SMF 行へ直接書いたフィクスチャは形式として存在せず、
# それに合わせた実装は実データで同定を 1 件も拾えなくなる。
_MZTAB_CONTENT = textwrap.dedent("""\n    MTD	mzTab-version	2.0.0-M
    MTD	mzTab-mode	Complete
    MTD	mzTab-type	Quantification
    MTD	ms_run[1]-location	file:///s1.raw
    MTD	assay[1]-ms_run_ref	ms_run[1]
    SFH	SMF_ID	SME_ID_REFS	exp_mass_to_charge	retention_time_in_seconds	abundance_assay[1]
    SMF	1	1	786.6	300.0	12345.6
    SMF	2	null	880.8	360.0	0.0
    SEH	SME_ID	database_identifier	smiles	inchi	chemical_name	rank
    SME	1	IPCSVZSSVZVIGE-UHFFFAOYSA-N	null	null	PC 36:2	1
""")


@pytest.fixture
def mztab_file(tmp_path):
    p = tmp_path / "Height_test.mzTab"
    p.write_text(_MZTAB_CONTENT, encoding="utf-8")
    return p


# SME が構造だけを持ち InChIKey を持たない形。RDKit が無いとここが 0 件になる。
_MZTAB_SMILES_ONLY = _MZTAB_CONTENT.replace(
    "SME	1	IPCSVZSSVZVIGE-UHFFFAOYSA-N	null	null	PC 36:2	1",
    "SME	1	null	CCCCCCCC	null	PC 36:2	1",
)


@pytest.fixture
def mztab_smiles_file(tmp_path):
    p = tmp_path / "Height_smiles.mzTab"
    p.write_text(_MZTAB_SMILES_ONLY, encoding="utf-8")
    return p


def test_build_dataset_state_creates_instance(mztab_file):
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    assert isinstance(ds, DatasetState)


def test_dataset_state_source_format(mztab_file):
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    assert ds.source_format == "mztab"


def test_dataset_state_feature_matrix_shape(mztab_file):
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    assert ds.feature_matrix.shape == (2, 1)
    assert ds.feature_matrix[0, 0] == pytest.approx(12345.6)


def test_dataset_state_inchikey_from_database_identifier(mztab_file):
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    f1 = ds.feature_metadata["1"]
    assert f1["inchikey"] == "IPCSVZSSVZVIGE-UHFFFAOYSA-N"
    assert f1["inchikey_source"] == "database_identifier"


def test_dataset_state_inchikey_none_when_absent(mztab_file):
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    f2 = ds.feature_metadata["2"]
    assert f2["inchikey"] is None
    assert f2["inchikey_source"] == "none"


def test_dataset_state_inchikey_coverage(mztab_file):
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    cov = ds.inchikey_coverage
    assert cov["total_features"] == 2
    assert cov["with_inchikey"] == 1
    assert cov["by_source"]["database_identifier"] == 1
    assert cov["by_source"]["none"] == 1


def test_session_has_dataset_slot():
    sess = session_state.AnalysisSession()
    assert hasattr(sess, "dataset")
    assert sess.dataset is None


def test_session_dataset_does_not_affect_arf_slot():
    sess = session_state.AnalysisSession()
    sess.dataset = "dummy"
    # ARF スロットはそのまま
    assert sess.arf.features is None


def test_dataset_state_has_analysis_fields():
    from metabolomix.mztab.dataset_state import DatasetState
    ds = DatasetState()
    assert ds.pp_matrix is None
    assert ds.pp_sample_names == []
    assert ds.pp_feature_names == []
    assert ds.roles == {}
    assert ds.sample_meta == {}
    assert ds.preprocessing_recipe == {}
    assert ds.last_pca is None
    assert ds.last_differential is None


# ---------- sample_names 解決（assay 表示名） ----------
# 実際の MS-DIAL 出力は MTD assay[N] の素の行に表示名を持つ
# （例: `MTD  assay[1]  20220901_RAW_control_0h_1_NEG`）。
# 以前はこれを拾えず sample_names が abundance_assay[N] のまま残っていた。

_MZTAB_WITH_NAMES = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tmzTab-mode\tComplete
    MTD\tmzTab-type\tQuantification
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tms_run[2]-location\tfile:///s2.raw
    MTD\tassay[1]\tzzz_last_NEG
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    MTD\tassay[2]\taaa_first_NEG
    MTD\tassay[2]-ms_run_ref\tms_run[2]
    SMF\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\tabundance_assay[1]\tabundance_assay[2]
    SMF\t1\tSML:1\tnull\tPC 36:2\tnull\tnull\t100.0\t200.0
    SMF\t2\tSML:2\tnull\tTG 54:3\tnull\tnull\t50.0\t60.0
""")


@pytest.fixture
def mztab_named_file(tmp_path):
    p = tmp_path / "Height_named.mzTab"
    p.write_text(_MZTAB_WITH_NAMES, encoding="utf-8")
    return p


def test_dataset_state_resolves_bare_assay_display_name(mztab_named_file):
    """素の `assay[N]` 行の表示名が sample_names に反映される。"""
    pr = parse_mztab(mztab_named_file)
    ds = build_dataset_state(pr, mztab_named_file.name, str(mztab_named_file))
    assert ds.sample_names == ["zzz_last_NEG", "aaa_first_NEG"]


def test_dataset_state_display_name_reachable_from_assay_metadata(mztab_named_file):
    """表示名は assay_metadata からも参照できる（dataset_status 等の将来利用のため）。"""
    pr = parse_mztab(mztab_named_file)
    ds = build_dataset_state(pr, mztab_named_file.name, str(mztab_named_file))
    assert ds.assay_metadata["assay[1]"]["name"] == "zzz_last_NEG"
    assert ds.assay_metadata["assay[2]"]["name"] == "aaa_first_NEG"
    # ms_run_ref など既存の suffix 情報も引き続き読める
    assert ds.assay_metadata["assay[1]"]["ms_run_ref"] == "ms_run[1]"


def test_dataset_state_sample_names_aligned_with_feature_matrix_columns(mztab_named_file):
    """sample_names[i] は feature_matrix[:, i]（= abundance_assay[i+1] 列）と対応していること。

    assay 番号昇順のまま列は並ぶので、表示名のアルファベット順（aaa < zzz）とは
    一致しない列0=zzz_last_NEG を選んで、順序の取り違えが起きていないことを確認する。
    """
    pr = parse_mztab(mztab_named_file)
    ds = build_dataset_state(pr, mztab_named_file.name, str(mztab_named_file))
    assert ds.sample_names[0] == "zzz_last_NEG"
    assert ds.feature_matrix[:, 0].tolist() == [100.0, 50.0]
    assert ds.sample_names[1] == "aaa_first_NEG"
    assert ds.feature_matrix[:, 1].tolist() == [200.0, 60.0]


def test_dataset_state_sample_names_fallback_without_display_name(mztab_file):
    """表示名が無い assay（既存フィクスチャ全体がこれに該当）は列識別子のままにする。"""
    pr = parse_mztab(mztab_file)
    ds = build_dataset_state(pr, mztab_file.name, str(mztab_file))
    assert ds.sample_names == ["abundance_assay[1]"]


_MZTAB_DUPLICATE_NAMES = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tmzTab-mode\tComplete
    MTD\tmzTab-type\tQuantification
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tms_run[2]-location\tfile:///s2.raw
    MTD\tassay[1]\tdup_sample
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    MTD\tassay[2]\tdup_sample
    MTD\tassay[2]-ms_run_ref\tms_run[2]
    SMF\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\tabundance_assay[1]\tabundance_assay[2]
    SMF\t1\tSML:1\tnull\tPC 36:2\tnull\tnull\t100.0\t200.0
""")


@pytest.fixture
def mztab_duplicate_names_file(tmp_path):
    p = tmp_path / "Height_dup.mzTab"
    p.write_text(_MZTAB_DUPLICATE_NAMES, encoding="utf-8")
    return p


def test_dataset_state_duplicate_display_name_warns_but_does_not_raise(mztab_duplicate_names_file):
    """不正ファイルで表示名が重複しても例外にはせず、warning に記録する。

    重複を検出できないと group_a/group_b の名前解決（sample_names.index 相当）が
    どちらのアッセイを指しているか曖昧になり、群選択が黙って誤る恐れがある。
    """
    pr = parse_mztab(mztab_duplicate_names_file)
    ds = build_dataset_state(pr, mztab_duplicate_names_file.name, str(mztab_duplicate_names_file))
    assert ds.sample_names == ["dup_sample", "dup_sample"]
    assert any("dup_sample" in w for w in ds.validation_result["warnings"])


# ---------- end-to-end: 表示名解決が QC/blank ロール検出を機能させる ----------
# detect_sample_roles はサンプル名のトークンで QC/blank を判定するため、
# sample_names が abundance_assay[N] のままだと QC-RSD フィルタ・ブランク除去が
# 常に不発になる（バグの結論2）。build_dataset_state → run_dataset_preprocess の
# 実経路を通して初めて「直った」と言える。

_MZTAB_WITH_QC_NAME = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tmzTab-mode\tComplete
    MTD\tmzTab-type\tQuantification
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tms_run[2]-location\tfile:///s2.raw
    MTD\tms_run[3]-location\tfile:///s3.raw
    MTD\tassay[1]\tcontrol_1
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    MTD\tassay[2]\ttreat_1
    MTD\tassay[2]-ms_run_ref\tms_run[2]
    MTD\tassay[3]\tQC_1
    MTD\tassay[3]-ms_run_ref\tms_run[3]
    SMF\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\tabundance_assay[1]\tabundance_assay[2]\tabundance_assay[3]
    SMF\t1\tSML:1\tnull\tPC 36:2\tnull\tnull\t100.0\t200.0\t150.0
    SMF\t2\tSML:2\tnull\tTG 54:3\tnull\tnull\t50.0\t60.0\t55.0
""")


@pytest.fixture
def mztab_qc_named_file(tmp_path):
    p = tmp_path / "Height_qc.mzTab"
    p.write_text(_MZTAB_WITH_QC_NAME, encoding="utf-8")
    return p


def test_dataset_preprocess_detects_qc_role_via_resolved_display_name(mztab_qc_named_file):
    """dataset_load 相当（build_dataset_state）→ dataset_preprocess 相当（run_dataset_preprocess）
    を実際に通し、QC_1 という表示名を持つアッセイが role="qc" に分類されることを確認する。
    直る前は sample_names が abundance_assay[3] のままで、detect_sample_roles が
    "qc" トークンを見つけられず全サンプルが "sample" のままだった。
    """
    from metabolomix.analysis.dataset_analysis import run_dataset_preprocess

    pr = parse_mztab(mztab_qc_named_file)
    ds = build_dataset_state(pr, mztab_qc_named_file.name, str(mztab_qc_named_file))
    assert ds.sample_names == ["control_1", "treat_1", "QC_1"]

    _, _, _, roles, _, _ = run_dataset_preprocess(ds, {})

    assert roles["QC_1"] == "qc"
    assert roles["control_1"] == "sample"
    assert roles["treat_1"] == "sample"


# --- 実 mzTab-M 形状（SMF は構造情報を持たず、SME を SME_ID_REFS で参照する） ---
# 旧フィクスチャ (_MZTAB_CONTENT) は database_identifier / chemical_name / smiles を
# SMF 行に直接置いていたが、mzTab-M 2.0.0-M ではこれらは **SME セクション専用**で、
# SMF が持つのは SMF_ID / SME_ID_REFS / exp_mass_to_charge /
# retention_time_in_seconds / abundance_assay[N] である。実データ
# (MS-DIAL Console 出力) もこの形。実装がフィクスチャに合わせて誤った列名を
# 読んでいたため、実ファイルでは全特徴の name/mz/rt/inchikey が None になっていた。
_MZTAB_REAL_SHAPE = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tmzTab-mode\tComplete
    MTD\tmzTab-type\tQuantification
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SFH\tSMF_ID\tSME_ID_REFS\tadduct_ion\texp_mass_to_charge\tcharge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\t1\t[M-H]1-\t101.06208\t-1\t472.701\t12345.6
    SMF\t2\tnull\t[M-H]1-\t157.12462\t-1\t480.0\t2000.0
    SMF\t3\t3|4\t[M-H]1-\t200.00000\t-1\t600.0\t500.0
    SEH\tSME_ID\tdatabase_identifier\tsmiles\tinchi\tchemical_name\trank
    SME\t1\tIPCSVZSSVZVIGE-UHFFFAOYSA-N\tnull\tnull\tPC 36:2\t1
    SME\t3\tsomedb:FA 9:0\tnull\tnull\tFA 9:0\t2
    SME\t4\tFBUKVWPVBMHYJY-UHFFFAOYSA-N\tnull\tnull\tFA 5:0\t1
""")


@pytest.fixture
def mztab_real_shape(tmp_path):
    p = tmp_path / "Height_real_shape.mzTab"
    p.write_text(_MZTAB_REAL_SHAPE, encoding="utf-8")
    return p


def _build(path):
    return build_dataset_state(parse_mztab(path), path.name, str(path))


def test_mz_read_from_smf_exp_mass_to_charge(mztab_real_shape):
    ds = _build(mztab_real_shape)
    assert ds.feature_metadata["1"]["mz"] == pytest.approx(101.06208)


def test_rt_converted_from_seconds_to_minutes(mztab_real_shape):
    """rt は分で持つ。ARF 経路が分で書き、同じエクスポート契約の rt 列を共有するため。"""
    ds = _build(mztab_real_shape)
    assert ds.feature_metadata["1"]["rt"] == pytest.approx(472.701 / 60.0)


def test_inchikey_resolved_through_sme_id_refs(mztab_real_shape):
    ds = _build(mztab_real_shape)
    f1 = ds.feature_metadata["1"]
    assert f1["inchikey"] == "IPCSVZSSVZVIGE-UHFFFAOYSA-N"
    assert f1["inchikey_source"] == "database_identifier"


def test_name_resolved_through_sme_id_refs(mztab_real_shape):
    ds = _build(mztab_real_shape)
    assert ds.feature_metadata["1"]["name"] == "PC 36:2"


def test_best_ranked_sme_wins_when_multiple_refs(mztab_real_shape):
    """SME_ID_REFS が複数あるときは rank が小さい（＝上位）証拠を採る。"""
    ds = _build(mztab_real_shape)
    f3 = ds.feature_metadata["3"]
    assert f3["name"] == "FA 5:0"
    assert f3["inchikey"] == "FBUKVWPVBMHYJY-UHFFFAOYSA-N"


def test_feature_without_sme_ref_has_no_identity(mztab_real_shape):
    ds = _build(mztab_real_shape)
    f2 = ds.feature_metadata["2"]
    assert f2["inchikey"] is None
    assert f2["inchikey_source"] == "none"
    assert f2["mz"] == pytest.approx(157.12462)


def test_inchikey_coverage_counts_sme_resolved_features(mztab_real_shape):
    ds = _build(mztab_real_shape)
    cov = ds.inchikey_coverage
    assert cov["total_features"] == 3
    assert cov["with_inchikey"] == 2
    assert cov["by_source"]["none"] == 1


def test_dataset_warnings_precede_parser_level_warnings(tmp_path):
    """データセット層の warning（表示名の重複など）を先頭に置く。

    パーサ層の良性 warning は実データで数百件になり得る。後ろに append すると、
    件数を絞って表示する dataset_load の要約に載る余地がなくなる。
    「見えない warning」は無いのと同じなので、順序で優先度を表す。
    """
    content = textwrap.dedent("""\
        MTD\tmzTab-version\t2.0.0-M
        MTD\tmzTab-mode\tComplete
        MTD\tmzTab-type\tQuantification
        MTD\tms_run[1]-location\tfile:///s1.raw
        MTD\tassay[1]\tdup_sample
        MTD\tassay[1]-ms_run_ref\tms_run[1]
        MTD\tassay[2]\tdup_sample
        MTD\tassay[2]-ms_run_ref\tms_run[2]
        SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]\tabundance_assay[2]
        SMF\t1\t1\t786.6\t300.0\t100.0\t200.0
        SEH\tSME_ID\tdatabase_identifier\tsmiles\tinchi\tchemical_name\trank
        SME\t1\tnull\tnull\tnull\tPC 36:2\t1\t
    """)
    p = tmp_path / "Height_dup_and_trailing.mzTab"
    p.write_text(content, encoding="utf-8")

    pr = parse_mztab(p)
    ds = build_dataset_state(pr, p.name, str(p))
    warnings = ds.validation_result["warnings"]
    assert "dup_sample" in warnings[0]
    assert any("trailing" in w.lower() for w in warnings)


_MZTAB_WITH_CUSTOM = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tassay[1]\tS_first
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    MTD\tassay[1]-custom[1]\t[MS,MS:4000088,batch label,B1]
    MTD\tassay[1]-custom[2]\t[MS,MS:4000089,injection sequence label,7]
    MTD\tassay[2]\tS_second
    MTD\tassay[2]-ms_run_ref\tms_run[2]
    MTD\tassay[2]-custom[1]\t[MS,MS:4000088,batch label,B2]
    MTD\tassay[2]-custom[2]\t[MS,MS:4000089,injection sequence label,3]
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]\tabundance_assay[2]
    SMF\t1\tnull\t786.6\t300.0\t10.0\t20.0
""")


@pytest.fixture
def mztab_with_custom(tmp_path):
    p = tmp_path / "AlignResult-1.mzTab"
    p.write_text(_MZTAB_WITH_CUSTOM, encoding="utf-8")
    return p


def test_parse_cv_term_extracts_accession_and_value():
    from metabolomix.mztab.dataset_state import _parse_cv_term
    assert _parse_cv_term("[MS,MS:4000089,injection sequence label,7]") == ("MS:4000089", "7")


def test_parse_cv_term_returns_none_for_garbage():
    from metabolomix.mztab.dataset_state import _parse_cv_term
    assert _parse_cv_term("not a cv term") == (None, None)
    assert _parse_cv_term("[MS,MS:4000089]") == (None, None)
    assert _parse_cv_term(None) == (None, None)


def test_build_dataset_state_reads_injection_order(mztab_with_custom):
    ds = _build(mztab_with_custom)
    assert ds.assay_metadata["assay[1]"]["run_order"] == 7
    assert ds.assay_metadata["assay[2]"]["run_order"] == 3


def test_build_dataset_state_reads_batch_label(mztab_with_custom):
    ds = _build(mztab_with_custom)
    assert ds.assay_metadata["assay[1]"]["batch"] == "B1"
    assert ds.assay_metadata["assay[2]"]["batch"] == "B2"


def test_build_dataset_state_records_sample_assay_ids(mztab_with_custom):
    ds = _build(mztab_with_custom)
    assert ds.sample_names == ["S_first", "S_second"]
    assert ds.sample_assay_ids == ["assay[1]", "assay[2]"]


def test_build_dataset_state_without_custom_terms_has_none(mztab_file):
    ds = _build(mztab_file)
    assert ds.assay_metadata["assay[1]"].get("run_order") is None
    assert ds.assay_metadata["assay[1]"].get("batch") is None


def test_rdkit_available_reports_bool():
    from metabolomix.mztab.identity import rdkit_available
    assert isinstance(rdkit_available(), bool)


def test_inchikey_coverage_reports_rdkit_availability(mztab_file):
    ds = _build(mztab_file)
    assert "rdkit_available" in ds.inchikey_coverage


# RDKit 非搭載環境の derive_inchikey と同じ挙動（database_identifier だけ拾う）。
# rdkit_available だけを False にしても実 RDKit が導出してしまうため、
# 導出側も同時に落として実環境を再現する。
def _no_rdkit_derive(database_identifier, inchi, smiles):
    from metabolomix.mztab.identity import _INCHIKEY_RE
    if database_identifier and _INCHIKEY_RE.match(database_identifier.strip()):
        return database_identifier.strip(), "database_identifier"
    return None, "none"


def _without_rdkit(monkeypatch):
    monkeypatch.setattr("metabolomix.mztab.identity.rdkit_available", lambda: False)
    monkeypatch.setattr("metabolomix.mztab.dataset_state.derive_inchikey", _no_rdkit_derive)


def test_dataset_state_warns_when_rdkit_missing(mztab_smiles_file, monkeypatch):
    """SMILES があるのに導出できていない件数を名指しする。"""
    _without_rdkit(monkeypatch)
    ds = _build(mztab_smiles_file)
    assert ds.inchikey_coverage["rdkit_available"] is False
    rdkit_warnings = [w for w in ds.validation_result["warnings"] if "RDKit" in w]
    assert len(rdkit_warnings) == 1
    assert "1 件" in rdkit_warnings[0]


def test_dataset_state_no_rdkit_warning_when_nothing_derivable(mztab_file, monkeypatch):
    """SMILES / InChI が無いなら RDKit があっても結果は変わらない。警告しない。"""
    _without_rdkit(monkeypatch)
    ds = _build(mztab_file)
    assert ds.inchikey_coverage["rdkit_available"] is False
    assert ds.inchikey_coverage["with_inchikey"] == 1
    assert not any("RDKit" in w for w in ds.validation_result.get("warnings", []))


# ---------- SML 由来の注釈（spec 2026-09-17-mztab-sml-annotation-design）----------

# Text DB 運用の実形状: SME セクションが 0 行で、同定は SML にしか出ない。
# MS-DIAL の `ShouldWriteSmeLine` が `IsTextDbBasedRepresentative` を除外するため
# （MztabFormatExport.cs）。ヘッダ行は実形式の `SMH` で書く。
_SML_ONLY_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tmzTab-mode\tComplete
    MTD\tmzTab-type\tQuantification
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_formula\tsmiles\tinchi\tchemical_name\tadduct_ions\treliability\tbest_id_confidence_measure\tbest_id_confidence_value
    SML\t1\t1\tTextDB:GABA\tnull\tnull\tnull\tGABA\t[M+H]1+\tannotated by user-defined text library\t[,, MS-DIAL algorithm matching score, ]\t0.999982
    SML\t2\t2\tnull\tnull\tnull\tnull\tnull\t[M+H]1+\tnull\tnull\tnull
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\tnull\t104.07066\t72.0\t5000.0
    SMF\t2\tnull\t880.8\t360.0\t100.0
""")


def _build_text(tmp_path, content, name="Height_sml.mzTab"):
    """本文から mzTab を書き起こして DatasetState を作る。

    既存の `_build(path)` とは別物（あちらは fixture が作った path を取る）。
    同名にすると既存テストの呼び出しを静かに壊す。
    """
    p = tmp_path / name
    p.write_text(content, encoding="utf-8")
    return build_dataset_state(parse_mztab(p), p.name, p)


def test_an_sml_annotation_is_recorded_for_its_feature(tmp_path):
    ds = _build_text(tmp_path, _SML_ONLY_MZTAB)

    annotation = ds.feature_annotations["1"]
    assert annotation["name"] == "GABA"
    assert annotation["database_identifier"] == "TextDB:GABA"
    assert annotation["adduct"] == "[M+H]1+"
    assert annotation["confidence_value"] == pytest.approx(0.999982)
    assert annotation["ambiguous"] is False


def test_an_sml_row_without_identity_is_not_recorded(tmp_path):
    ds = _build_text(tmp_path, _SML_ONLY_MZTAB)

    # SML 2 は chemical_name も database_identifier も null。
    assert "2" not in ds.feature_annotations


def test_the_annotation_never_holds_an_inchi_key(tmp_path):
    """MS-DIAL は SML の `inchi` を常に null で書く（MztabFormatExport.cs:393）。

    キーを置くと「取得していない」と「無い」の区別を偽るので、持たない。
    """
    ds = _build_text(tmp_path, _SML_ONLY_MZTAB)

    assert "inchi" not in ds.feature_annotations["1"]


def test_the_evidence_slots_stay_untouched_by_an_sml_annotation(tmp_path):
    """SML 注釈は証拠スロットへ一切入らない（spec §5 の不変条件）。"""
    ds = _build_text(tmp_path, _SML_ONLY_MZTAB)

    assert ds.feature_metadata["1"]["name"] is None
    assert ds.feature_metadata["1"]["inchikey"] is None
    assert ds.feature_candidates["1"] == []


_MULTI_REF_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tsmiles\tchemical_name\tadduct_ions
    SML\t1\t1|2\tTextDB:GABA\tnull\tGABA\t[M+H]1+
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\tnull\t104.07\t72.0\t1.0
    SMF\t2\tnull\t126.05\t72.0\t2.0
""")

_AMBIGUOUS_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tsmiles\tchemical_name\tadduct_ions
    SML\t1\t1\tTextDB:GABA\tnull\tGABA\t[M+H]1+
    SML\t2\t1\tTextDB:Alanine\tnull\tAlanine\t[M+H]1+
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\tnull\t104.07\t72.0\t1.0
""")

_UNKNOWN_REF_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tsmiles\tchemical_name\tadduct_ions
    SML\t1\t99\tTextDB:GABA\tnull\tGABA\t[M+H]1+
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\tnull\t104.07\t72.0\t1.0
""")

_SMILES_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tsmiles\tchemical_name\tadduct_ions
    SML\t1\t1\tTextDB:GABA\tNCCCC(=O)O\tGABA\t[M+H]1+
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\tnull\t104.07\t72.0\t1.0
""")


def test_one_sml_pointing_at_two_features_annotates_both(tmp_path):
    """同一分子が複数 adduct で検出された形。両方に同じ注釈が付くのが正しい。"""
    ds = _build_text(tmp_path, _MULTI_REF_MZTAB, "Height_multi.mzTab")

    assert ds.feature_annotations["1"]["name"] == "GABA"
    assert ds.feature_annotations["2"]["name"] == "GABA"


def test_two_sml_rows_on_one_feature_do_not_pick_a_name(tmp_path):
    """どれが正しいか決められないときは選ばない（リポジトリ共通の方針）。"""
    ds = _build_text(tmp_path, _AMBIGUOUS_MZTAB, "Height_ambiguous.mzTab")

    annotation = ds.feature_annotations["1"]
    assert annotation["ambiguous"] is True
    assert annotation["name"] is None
    assert sorted(annotation["sml_ids"]) == ["1", "2"]


def test_an_ambiguous_annotation_is_reported_once(tmp_path):
    """警告は件数ではなく種類で 1 件に集約する。"""
    ds = _build_text(tmp_path, _AMBIGUOUS_MZTAB, "Height_ambiguous2.mzTab")

    hits = [w for w in ds.validation_result.get("warnings", []) if "複数の SML" in w]
    assert len(hits) == 1


def test_an_sml_pointing_at_a_missing_feature_is_reported(tmp_path):
    ds = _build_text(tmp_path, _UNKNOWN_REF_MZTAB, "Height_unknown.mzTab")

    assert ds.feature_annotations == {}
    hits = [w for w in ds.validation_result.get("warnings", []) if "SMF_ID_REFS" in w]
    assert len(hits) == 1


def test_an_inchikey_is_derived_from_the_sml_smiles(tmp_path):
    """MS-DIAL の SML で InChIKey に到達しうるのは smiles 経由だけ。

    `database_identifier` は必ず `<db>:<name>` 形式（MztabFormatExport.cs:407）、
    `inchi` は常に null（同 :393）。
    """
    # `importorskip("rdkit")` ではガードにならない。`rdkit` パッケージ本体は純 Python
    # なので import は通り、native な `rdchem` を引く `rdkit.Chem` だけが落ちる環境が
    # ある（Windows の Application Control は**ファイル単位**で効く）。必要なのは
    # 「パッケージが在るか」ではなく「InChI を導出できるか」なので、本番と同じ判定を使う。
    from metabolomix.mztab.identity import rdkit_available
    if not rdkit_available():
        pytest.skip("RDKit の native DLL を読めない環境では smiles 経路を検証できない")
    ds = _build_text(tmp_path, _SMILES_MZTAB, "Height_smiles.mzTab")

    annotation = ds.feature_annotations["1"]
    assert annotation["inchikey"] == "BTCSSZJGUNDROE-UHFFFAOYSA-N"
    assert annotation["inchikey_source"] == "smiles_derived"


def test_the_coverage_separates_evidence_from_ms1_annotation(tmp_path):
    ds = _build_text(tmp_path, _SML_ONLY_MZTAB, "Height_cov.mzTab")

    by = ds.inchikey_coverage["identified_by"]
    assert by["sme"] == 0
    assert by["sml_only"] == 1
    assert by["none"] == 1


_SUBSCORE_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    MTD\tid_confidence_measure[1]\t[,, MS-DIAL algorithm matching score, ]
    MTD\tid_confidence_measure[2]\t[,, Retention time similarity, ]
    MTD\tid_confidence_measure[3]\t[,, m/z similarity, ]
    MTD\tid_confidence_measure[4]\t[,, Simple dot product, ]
    MTD\tid_confidence_measure[5]\t[,, Weighted dot product, ]
    MTD\tid_confidence_measure[6]\t[,, Reverse dot product, ]
    MTD\tid_confidence_measure[7]\t[,, Matched peaks count, ]
    MTD\tid_confidence_measure[8]\t[,, Matched peaks percentage, ]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_name\tadduct_ions
    SML\t1\t1\tnull\tGABA\t[M+H]1+
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\t1\t104.07\t72.0\t1.0
    SEH\tSME_ID\tevidence_input_id\tdatabase_identifier\tchemical_name\tadduct_ion\texp_mass_to_charge\tbest_id_confidence_measure\tbest_id_confidence_value\tid_confidence_measure[1]\tid_confidence_measure[2]\tid_confidence_measure[3]\tid_confidence_measure[4]\tid_confidence_measure[5]\tid_confidence_measure[6]\tid_confidence_measure[7]\tid_confidence_measure[8]\trank
    SME\t1\t1\tnull\tGABA\t[M+H]1+\t104.07\t[,, MS-DIAL algorithm matching score, ]\t0.91\t0.91\tnull\t0.99\t0.95\t0.93\t0.88\t7\t0.7\t1
""")


def test_the_sub_scores_are_recovered_from_the_sme_row(tmp_path):
    """mzTab は [4..8] に個別スコアを出している。total だけ読むと捨てることになる。"""
    ds = _build_text(tmp_path, _SUBSCORE_MZTAB, "Height_subscores.mzTab")

    measures = ds.feature_annotations["1"]["confidence_measures"]

    assert measures["simple_dot_product"] == pytest.approx(0.95)
    assert measures["weighted_dot_product"] == pytest.approx(0.93)
    assert measures["reverse_dot_product"] == pytest.approx(0.88)
    assert measures["matched_peaks_count"] == pytest.approx(7.0)
    assert measures["matched_peaks_percentage"] == pytest.approx(0.7)
    assert measures["mz_similarity"] == pytest.approx(0.99)
    assert "retention_time_similarity" not in measures       # null の列は入れない


def test_the_existing_total_score_reading_is_unchanged(tmp_path):
    """best_id_confidence_value は下流の契約。サブスコア回収で壊さない。"""
    ds = _build_text(tmp_path, _SUBSCORE_MZTAB, "Height_total.mzTab")
    assert ds.feature_annotations["1"]["confidence_value"] == pytest.approx(0.91)


# manualAssigned 相当（上流 SetIdConfidenceMeasure が付け足す9本目）＋対応表に無い
# 名前の両方を一度に張る fixture。range(1, 9) のような決め打ちへ「簡略化」されると
# ここが真っ先に壊れる（Task 11 レビュー Important 2）。
_NINTH_MEASURE_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t2.0.0-M
    MTD\tms_run[1]-location\tfile:///s1.raw
    MTD\tassay[1]-ms_run_ref\tms_run[1]
    MTD\tid_confidence_measure[1]\t[,, MS-DIAL algorithm matching score, ]
    MTD\tid_confidence_measure[2]\t[,, Retention time similarity, ]
    MTD\tid_confidence_measure[3]\t[,, m/z similarity, ]
    MTD\tid_confidence_measure[4]\t[,, Simple dot product, ]
    MTD\tid_confidence_measure[5]\t[,, Weighted dot product, ]
    MTD\tid_confidence_measure[6]\t[,, Reverse dot product, ]
    MTD\tid_confidence_measure[7]\t[,, Matched peaks count, ]
    MTD\tid_confidence_measure[8]\t[,, Matched peaks percentage, ]
    MTD\tid_confidence_measure[9]\t[MS, MS:1001058, quality estimation by manual validation, ]
    SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_name\tadduct_ions
    SML\t1\t1\tnull\tGABA\t[M+H]1+
    SFH\tSMF_ID\tSME_ID_REFS\texp_mass_to_charge\tretention_time_in_seconds\tabundance_assay[1]
    SMF\t1\t1\t104.07\t72.0\t1.0
    SEH\tSME_ID\tevidence_input_id\tdatabase_identifier\tchemical_name\tadduct_ion\texp_mass_to_charge\tbest_id_confidence_measure\tbest_id_confidence_value\tid_confidence_measure[1]\tid_confidence_measure[2]\tid_confidence_measure[3]\tid_confidence_measure[4]\tid_confidence_measure[5]\tid_confidence_measure[6]\tid_confidence_measure[7]\tid_confidence_measure[8]\tid_confidence_measure[9]\trank
    SME\t1\t1\tnull\tGABA\t[M+H]1+\t104.07\t[,, MS-DIAL algorithm matching score, ]\t0.91\t0.91\tnull\t0.99\t0.95\t0.93\t0.88\t7\t0.7\t0.5\t1
""")


def test_a_manually_assigned_ninth_measure_and_an_unmapped_name_are_not_dropped(tmp_path):
    """manualAssigned は9本目を追加するだけ（SetIdConfidenceMeasure）で、対応表に無い
    名前もある。range(1, 9) のような固定範囲へ「簡略化」すると、この9本目が
    真っ先に静かに消える。位置ではなく宣言を読んでいることをここで縛る。"""
    ds = _build_text(tmp_path, _NINTH_MEASURE_MZTAB, "Height_ninth.mzTab")

    measures = ds.feature_annotations["1"]["confidence_measures"]

    # 対応表に無い宣言名は、そのままスネークケース化されて残る（捨てない）。
    assert measures["quality_estimation_by_manual_validation"] == pytest.approx(0.5)
    # 9本目が増えても 1..8 の対応関係はずれない。
    assert measures["simple_dot_product"] == pytest.approx(0.95)
    assert measures["matched_peaks_percentage"] == pytest.approx(0.7)
    assert "retention_time_similarity" not in measures
