import numpy as np
import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.quant.table import FeatureRow, apply_curation, resolve_quant_polarity, feature_rows
from tests.msdial_folder_fixtures import write_msdial_folder


def _samples(n=3):
    return [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(n)]


def _load(tmp_path, spots):
    write_msdial_folder(tmp_path, spots=spots, samples=_samples())
    return load_msdial_folder(select_alignment(tmp_path))


class _Holder:
    def __init__(self, polarity=None):
        self.polarity = polarity


def test_rows_follow_workflow_conventions(tmp_path):
    ds = _load(tmp_path, [
        {"name": "PC 34:1", "ontology": "PC", "mz": 760.5851, "heights": [100.0, 0.0, 50.0],
         "gap_filled": [False, True, False]},
        {"name": "", "mz": 500.1, "heights": [10.0, 20.0, 30.0]}])
    rows = feature_rows(ds)
    assert [r.feature_id for r in rows] == ds.feature_ids == ["0", "1"]
    assert rows[0].name == "PC 34:1" and rows[0].ontology == "PC"
    assert rows[0].detected == (True, False, True)
    assert rows[0].msms_assigned in ("True", "False")
    assert rows[1].name == "Unknown"
    assert isinstance(rows[0].heights, tuple) and len(rows[0].heights) == 3


def test_none_fields_and_nan_heights(tmp_path):
    ds = _load(tmp_path, [{"name": "PC 34:1"}])
    ds.feature_metadata["0"].update({"ontology": None, "adduct": None, "comment": None})
    ds.feature_matrix[0, 1] = np.nan
    r = feature_rows(ds)[0]
    assert (r.ontology, r.adduct, r.comment) == ("", "", "")
    assert r.heights[1] is None


def test_needs_detection_mask(tmp_path):
    ds = _load(tmp_path, [{"name": "a"}])
    ds.detected_mask = None
    with pytest.raises(DomainError) as exc:
        feature_rows(ds)
    assert exc.value.code == "QUANT_NEEDS_DETECTION"
    ds2 = _load(tmp_path / "b", [{"name": "a"}])
    ds2.feature_matrix = None
    with pytest.raises(DomainError) as exc:
        feature_rows(ds2)
    assert exc.value.code == "QUANT_NEEDS_DETECTION"


def test_polarity_required_when_unknown(tmp_path):
    ds = _load(tmp_path, [{"name": "a"}])
    ds.ion_mode_majority = None
    with pytest.raises(DomainError) as exc:
        resolve_quant_polarity(ds, object(), None)
    assert exc.value.code == "QUANT_POLARITY_REQUIRED"
    assert "positive" in exc.value.message


def test_polarity_from_ion_mode(tmp_path):
    ds = _load(tmp_path, [{"name": "a"}])
    ds.ion_mode_majority = "negative"
    assert resolve_quant_polarity(ds, object(), None) == "negative"
    assert resolve_quant_polarity(ds, _Holder("Positive"), None) == "positive"
    assert resolve_quant_polarity(ds, _Holder("pos"), "NEG") == "negative"
    assert resolve_quant_polarity(ds, _Holder(""), None) == "negative"


def test_polarity_invalid_explicit_values(tmp_path):
    ds = _load(tmp_path, [{"name": "a"}])
    ds.ion_mode_majority = "negative"
    with pytest.raises(DomainError) as exc:
        resolve_quant_polarity(ds, object(), "postive")
    assert exc.value.code == "QUANT_POLARITY_INVALID" and exc.value.details == {"polarity": "postive"}
    with pytest.raises(DomainError) as exc:
        resolve_quant_polarity(ds, _Holder("garbage"), None)
    assert exc.value.code == "QUANT_POLARITY_INVALID"


def _row(fid, name="X", ontology="PC", adduct="[M+H]+"):
    return FeatureRow(fid, name, ontology, adduct, 500.0, 10.0, "", "False", (1.0,), (True,))


def test_curation_wrong_redundant_assign():
    rows = [_row("1"), _row("2"), _row("3"), _row("4"), _row("5")]
    flags = {"wrong": {1}, "suspect": {5}, "redundant": {2},
             "assign": {3: {"name": "PE 36:2", "ontology": "PE", "adduct": "[M+Na]+"}, 1: {"name": "ignored"},
                        4: {"ontology": "PS"}}}
    out, audit = apply_curation(rows, flags)
    assert [r.feature_id for r in out] == ["3", "4", "5"]
    assert (out[1].name, out[1].ontology) == ("Unknown", "PS")  # name blanked -> Unknown
    r3 = out[0]
    assert (r3.name, r3.ontology, r3.adduct) == ("PE 36:2", "PE", "[M+H]+")  # adduct never changed
    assert out[2] == rows[4]  # suspect untouched
    by_id = {a["feature_id"]: a for a in audit}
    assert by_id["1"] == {"feature_id": "1", "action": "excluded", "reason": "curation_wrong"}
    assert by_id["2"]["reason"] == "curation_redundant"
    assert by_id["3"]["action"] == "replaced" and by_id["3"]["reason"] == "curation_assign"
    assert by_id["3"]["from"]["name"] == "X" and by_id["3"]["to"] == {"name": "PE 36:2", "ontology": "PE"}
    assert len(audit) == 4  # wrong wins over assign: no replaced audit for 1


def test_curation_none_flags_is_noop():
    rows = [_row("1")]
    assert apply_curation(rows, None) == (rows, [])
