from tests.curation_fixtures import arf2_spot_raw, match_result, write_arf2
from metabolomix.arf2.spot_fields import load_spot_fields, spot_fields


def test_spot_fields_reads_identity_keys(tmp_path):
    raw = arf2_spot_raw(spot_id=7, name="PC 34:1", mz=760.5851, rt=12.0,
                        adduct="[M+H]+", matches=[match_result()])
    raw[25] = "checked"; raw[7] = 3; raw[58] = False
    f = spot_fields(raw)
    assert f["master_id"] == 7 and f["name"] == "PC 34:1"
    assert f["adduct"] == "[M+H]+" and f["ontology"] == "PC"
    assert f["comment"] == "checked" and f["internal_standard_id"] == 3
    assert f["rt"] == 12.0 and f["ion_mode"] == "Positive"
    assert f["representative"]["name"] == "PC 34:1"
    assert f["is_text_db_representative"] is False and f["is_manually_modified"] is False


def test_text_db_and_manual_representative_flags():
    text = arf2_spot_raw(spot_id=1, matches=[match_result({26: 16})])
    assert spot_fields(text)["is_text_db_representative"] is True
    in_list = arf2_spot_raw(spot_id=2, matches=[match_result()])
    in_list[56][2] = [match_result()]          # コンテナ Key2 TextDbBasedMatchResults
    assert spot_fields(in_list)["is_text_db_representative"] is True
    manual = arf2_spot_raw(spot_id=3, matches=[match_result({26: 4 | 64})])
    assert spot_fields(manual)["is_manually_modified"] is True


def test_empty_strings_become_none_and_no_match():
    raw = arf2_spot_raw(spot_id=4, ontology="", formula="", adduct="", name="")
    raw[15] = ""; raw[16] = ""
    f = spot_fields(raw)
    assert f["name"] is None and f["formula"] is None and f["adduct"] is None
    assert f["ontology"] is None and f["smiles"] is None and f["inchikey"] is None
    assert f["representative"] is None and f["n_candidates"] == 0
    assert f["is_text_db_representative"] is False


def test_load_spot_fields_keeps_file_order(tmp_path):
    path = write_arf2(tmp_path / "A.arf2", [arf2_spot_raw(spot_id=i) for i in (0, 1, 2)])
    assert [f["master_id"] for f in load_spot_fields(path)] == [0, 1, 2]


def test_empty_spot_array_is_reported():
    import pytest

    from metabolomix.arf2.match_results import spot_annotation
    from metabolomix.core.atomic_io import DomainError

    for call in (spot_fields, spot_annotation):
        for raw in ([], [None]):
            with pytest.raises(DomainError) as exc:
                call(raw)
            assert exc.value.code == "ARF2_SPOT_MALFORMED"


def test_spot_fields_values():
    """Key ごとの値（docs/schema/AlignmentSpotProperty.md の Key 番号表に照らす）。"""
    raw = arf2_spot_raw(spot_id=9, name="PC 34:1", mz=760.5851, rt=12.5, ontology="PC",
                        adduct="[M+H]+", formula="C42H82NO8P", ion_mode=1,
                        matches=[match_result()])
    raw[1] = 99          # Key1 AlignmentID は master_id ではない
    raw[7] = 4           # Key7 InternalStandardAlignmentID
    raw[15] = "CCSMILES"  # Key15 SMILES
    raw[16] = "INCHIKEY-X"  # Key16 InChIKey
    raw[25] = "note"     # Key25 Comment
    raw[49] = 0.75       # Key49 FillParcentage
    raw[58] = True       # Key58 IsBlankFilteredByPostCurator
    f = spot_fields(raw)
    assert f["master_id"] == 9                 # Key0
    assert f["mz"] == 760.5851                 # Key5 MassCenter
    assert f["rt"] == 12.5                     # Key4 TimesCenter
    assert f["ion_mode"] == "Negative"         # Key11（1=Negative）
    assert f["name"] == "PC 34:1"              # Key12
    assert f["formula"] == "C42H82NO8P"        # Key13 index0
    assert f["ontology"] == "PC"               # Key14
    assert f["smiles"] == "CCSMILES" and f["inchikey"] == "INCHIKEY-X"
    assert f["adduct"] == "[M+H]+"             # Key54 index2
    assert f["comment"] == "note"              # Key25
    assert f["internal_standard_id"] == 4      # Key7
    assert f["fill_percentage"] == 0.75        # Key49
    assert f["blank_filtered"] is True         # Key58
    assert f["n_candidates"] == 1              # Key56 の結果の一覧
