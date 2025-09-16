import sys
from pathlib import Path

import pandas as pd
import pandas.testing as pdt

sys.path.append(str(Path(__file__).resolve().parents[2]))

from tackle.utils import (
    clean_categorical,
    get_default_color_mapping,
    set_pandas_datatypes,
)


def test_clean_categorical_converts_non_string_categories_to_strings():
    numeric_categories = pd.Series([1, 2, 1], dtype="category")
    label_categories = pd.Series(["low", "high", "low"], dtype="category")
    df = pd.DataFrame(
        {
            "numeric_category": numeric_categories,
            "label_category": label_categories,
            "value": [0.1, 0.2, 0.3],
        }
    )

    cleaned = clean_categorical(df)

    assert list(cleaned.columns) == [
        "numeric_category",
        "label_category",
        "value",
    ]
    assert isinstance(cleaned["numeric_category"].dtype, pd.CategoricalDtype)
    assert isinstance(cleaned["label_category"].dtype, pd.CategoricalDtype)
    assert cleaned["numeric_category"].cat.categories.tolist() == ["1", "2"]
    assert cleaned["label_category"].cat.categories.tolist() == ["high", "low"]
    pdt.assert_series_equal(cleaned["value"], df["value"])


def test_set_pandas_datatypes_handles_mixed_inputs():
    df = pd.DataFrame(
        {
            "numeric": ["1", "2", "3"],
            "boolean": ["TRUE", "FALSE", "NA"],
            "text": ["alpha", "beta", "gamma"],
            "plex": ["10", "20", "30"],
        }
    )

    converted = set_pandas_datatypes(df)

    assert converted["numeric"].tolist() == [1, 2, 3]
    assert pd.api.types.is_integer_dtype(converted["numeric"].dtype)

    assert converted["boolean"].tolist() == ["True", "False", "NA"]
    assert converted["boolean"].dtype == object

    assert converted["text"].tolist() == ["alpha", "beta", "gamma"]
    assert converted["text"].dtype == object

    assert converted["plex"].tolist() == ["10", "20", "30"]
    assert converted["plex"].dtype == object


def test_get_default_color_mapping_returns_none_for_float_series():
    float_series = pd.Series([0.5, 1.2, 0.3], name="float_values")

    color_mapping = get_default_color_mapping(float_series)

    assert color_mapping is None


def test_get_default_color_mapping_returns_colors_for_categorical_series():
    categorical_series = pd.Series(
        ["control", "treated", "control"], dtype="category"
    )

    color_mapping = get_default_color_mapping(categorical_series)

    assert color_mapping is not None
    assert set(color_mapping.keys()) == {"control", "treated"}
    assert all(isinstance(color, str) for color in color_mapping.values())
    assert all(color.startswith("#") for color in color_mapping.values())
