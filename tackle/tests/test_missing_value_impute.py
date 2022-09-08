import numpy as np
import pytest
from pathlib import Path

import pandas as pd
from tackle.utils import impute_missing

@pytest.fixture
def rectangular_mat():
    nrows = 100

    nrows = 8
    ncols = 10
    mat = np.random.normal(size=(nrows, ncols))
    bool_mat = np.random.randint(0, 2, size=mat.shape)
    bool_mat = bool_mat.astype(bool)

    mat[ bool_mat ] = np.nan

    df = pd.DataFrame(mat)
    return df

def test_old_impute_missing(rectangular_mat):


    res = impute_missing(rectangular_mat, make_plot=False)
    assert res.isna().stack().all() == False
    return