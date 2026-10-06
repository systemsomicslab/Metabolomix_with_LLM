from metabolomix.msdial.analysis_params import (
    DEFAULT_ADDUCTS, find_param_file, read_analysis_params, resolve_analysis_params)

PARAM = """# Project information
Ion mode: Negative
Searched adduct ions: [M-H]-,[M+HCOO]-,[M+CH3COO]-,[2M-H]-
MS1 tolerance for centroid: 0.01
Retention time tolerance for alignment: 0.1
"""


def test_read_params(tmp_path):
    path = tmp_path / "Dataset_2026_09_09_17_28_59_param_202609091800.txt"
    path.write_text(PARAM, encoding="utf-8")
    assert read_analysis_params(path) == {
        "ion_mode": "Negative", "searched_adducts": ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[2M-H]-"],
        "rt_tolerance_alignment": 0.1, "ms1_tolerance": 0.01}


def test_find_param_file_picks_the_newest_name(tmp_path):
    (tmp_path / "Dataset_a_param_202601010000.txt").write_text(PARAM, encoding="utf-8")
    (tmp_path / "Dataset_a_param_202609091800.txt").write_text(PARAM, encoding="utf-8")
    arf2 = tmp_path / "AlignmentResult_x.arf2"
    assert find_param_file(arf2).name == "Dataset_a_param_202609091800.txt"


def test_resolve_falls_back_to_defaults_without_a_param_file(tmp_path):
    resolved = resolve_analysis_params(tmp_path / "AlignmentResult_x.arf2", "Negative")
    assert resolved == {"searched_adducts": DEFAULT_ADDUCTS["Negative"], "rt_window": 0.1,
                        "source": "default", "path": None}


def test_resolve_uses_the_param_file(tmp_path):
    (tmp_path / "Dataset_a_param_1.txt").write_text(PARAM.replace("0.1\n", "0.08\n"), encoding="utf-8")
    resolved = resolve_analysis_params(tmp_path / "AlignmentResult_x.arf2", "Negative")
    assert resolved["source"] == "param_file" and resolved["rt_window"] == 0.08
    assert resolved["searched_adducts"][0] == "[M-H]-"
