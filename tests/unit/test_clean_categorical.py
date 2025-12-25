import pandas as pd

from tackle.utils import clean_categorical


def test_clean_categorical_multiple_columns():
    df = pd.DataFrame(
        {
            "first": pd.Series([1, 2, 1], dtype="category"),
            "second": pd.Series([3, 4, 3], dtype="category"),
        }
    )

    # Sanity check: categories start as non-strings.
    assert any(not isinstance(cat, str) for cat in df["first"].cat.categories)
    assert any(not isinstance(cat, str) for cat in df["second"].cat.categories)

    cleaned = clean_categorical(df)

    for column in ("first", "second"):
        assert all(isinstance(cat, str) for cat in cleaned[column].cat.categories)
