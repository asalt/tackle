import pandas as pd

from tackle.utils import normalize


def test_normalize_trim_mean_discards_lowest_and_highest_quartiles():
    df = pd.DataFrame(
        {
            "GeneID": ["1", "2", "3", "4"],
            "iBAQ_dstrAdj": [1.0, 2.0, 3.0, 100.0],
        }
    )

    out = normalize(df, trim_mean=True)

    expected_norm = 2.5
    expected = df["iBAQ_dstrAdj"] / expected_norm
    pd.testing.assert_series_equal(out.reset_index(drop=True), expected.reset_index(drop=True))
