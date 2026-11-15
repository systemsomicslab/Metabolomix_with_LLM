"""`core/polarity_names.normalize_polarity` の別名と、読めない値を誤りにすること。"""
import pytest

from metabolomix.core.polarity_names import normalize_polarity


@pytest.mark.parametrize("value, expected", [
    ("positive", "positive"), ("POS", "positive"), (" pos ", "positive"), ("+", "positive"),
    ("negative", "negative"), ("Neg", "negative"), ("-", "negative"),
])
def test_aliases(value, expected):
    assert normalize_polarity(value) == expected


@pytest.mark.parametrize("value", [None, "", "both", "positiv", 1, True])
def test_unreadable_values_raise(value):
    with pytest.raises(ValueError, match="ion_mode"):
        normalize_polarity(value, field="ion_mode")
