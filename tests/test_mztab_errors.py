import json
from metabolomix.core.mcp_errors import mztab_error, missing_state, MISSING_STATE


def test_mztab_error_returns_json_string():
    result = mztab_error("MZTAB_NOT_FOUND", "ファイルが見つかりません")
    parsed = json.loads(result)
    assert parsed["error"]["code"] == "MZTAB_NOT_FOUND"
    assert parsed["error"]["message"] == "ファイルが見つかりません"
    assert parsed["error"].get("details") is None


def test_mztab_error_with_details():
    result = mztab_error("QUANTIFICATION_CONFLICT", "定量値が衝突", {"mtd_unit": "height", "filename_prefix": "Area"})
    parsed = json.loads(result)
    assert parsed["error"]["details"]["mtd_unit"] == "height"


def test_mztab_error_code_is_not_missing_state():
    result = mztab_error("MZTAB_STRUCTURE_INVALID", "構造不正")
    parsed = json.loads(result)
    assert parsed["error"]["code"] != MISSING_STATE


def test_missing_state_still_works_after_change():
    result = missing_state("preprocessed_matrix", ["arf_preprocess"], "テスト")
    parsed = json.loads(result)
    assert parsed["error"]["code"] == MISSING_STATE
