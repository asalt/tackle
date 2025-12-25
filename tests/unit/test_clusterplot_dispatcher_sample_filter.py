import pandas as pd

from tackle.clusterplot_dispatcher import _series_matches_includes


def test_series_matches_includes_numeric_metadata_with_string_values():
    series = pd.Series([0.0, 1.0, 2.0, 1.0], index=["s1", "s2", "s3", "s4"])

    mask = _series_matches_includes(series, ["0", "2"])

    assert mask.to_dict() == {"s1": True, "s2": False, "s3": True, "s4": False}


def test_series_matches_includes_object_numeric_metadata_with_string_values():
    series = pd.Series([0.0, 1.0, 2.0, ""], index=["s1", "s2", "s3", "s4"], dtype=object)

    mask = _series_matches_includes(series, ["1"])

    assert mask.to_dict() == {"s1": False, "s2": True, "s3": False, "s4": False}

