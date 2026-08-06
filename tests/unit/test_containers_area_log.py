from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from tackle.containers import Data


def _area_data(*, fill_na_zero: bool) -> Data:
    data = Data.__new__(Data)
    data.df_filtered = pd.DataFrame(
        [[0.0, 10.0, np.nan, 20.0]],
        index=pd.MultiIndex.from_tuples(
            [("g1", "area")],
            names=["GeneID", "Metric"],
        ),
        columns=["zero", "minimum", "na", "higher"],
    )
    data.config = {}
    data.funcats = None
    data.geneid_subset = None
    data.impute_missing_values = False
    data.fill_na_zero = fill_na_zero
    data.norm_info = None
    data.export_all = False
    data.normed = False
    data.batch = None
    data._cache = SimpleNamespace(maybe_save=lambda _data: None)
    return data


@pytest.mark.parametrize("fill_na_zero", [True, False])
def test_area_log_uses_subminimum_floor_and_preserves_explicit_mask(fill_na_zero):
    data = _area_data(fill_na_zero=fill_na_zero)

    data.set_area_dfs()

    expected_minimum = np.log10(10.0 / 9.0)
    expected_higher = np.log10(20.0 / 9.0)
    assert data.minval == 10.0
    assert data.mask.loc["g1"].to_dict() == {
        "zero": True,
        "minimum": False,
        "na": True,
        "higher": False,
    }
    assert data.zeros.loc["g1"].to_dict() == {
        "zero": True,
        "minimum": False,
        "na": False,
        "higher": False,
    }
    assert np.isclose(data.areas_log.loc["g1", "minimum"], expected_minimum)
    assert np.isclose(data.areas_log.loc["g1", "higher"], expected_higher)
    assert data.areas_log.loc["g1", "minimum"] > 0

    missing_values = data.areas_log.loc["g1", ["zero", "na"]]
    if fill_na_zero:
        assert np.allclose(missing_values, 0.0)
    else:
        assert missing_values.isna().all()
