import numpy as np
import pytest
from pathlib import Path

import pandas as pd
import tackle.utils as utils
from tackle.utils import impute_missing


# how to make a teardown
# @pytest.fixture(scope='session')
# def clear_files_teardown():
#     yield None
#     os.system("rm -rf logs json")

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

def test_impute_missing_orig(rectangular_mat):

    res = impute_missing(rectangular_mat, make_plot=False)
    assert res.isna().stack().all() == False
    return


def test_impute_missing_gaussian_dispatches_legacy(monkeypatch, rectangular_mat):
    calls = {}

    def fake_old(frame, **kwargs):
        calls["legacy"] = kwargs
        return frame.fillna(0)

    def fake_mqish(frame, **kwargs):
        calls["mqish"] = kwargs
        return frame.fillna(1)

    monkeypatch.setattr(utils, "impute_missing_old", fake_old)
    monkeypatch.setattr(utils, "impute_missing_mqish", fake_mqish)

    res = utils.impute_missing_gaussian(
        rectangular_mat,
        method="legacy",
        make_plot=False,
    )

    assert "legacy" in calls
    assert "mqish" not in calls
    assert res.isna().sum().sum() == 0


def test_impute_missing_gaussian_dispatches_mqish(monkeypatch, rectangular_mat):
    calls = {}

    def fake_old(frame, **kwargs):
        calls["legacy"] = kwargs
        return frame.fillna(0)

    def fake_mqish(frame, **kwargs):
        calls["mqish"] = kwargs
        return frame.fillna(1)

    monkeypatch.setattr(utils, "impute_missing_old", fake_old)
    monkeypatch.setattr(utils, "impute_missing_mqish", fake_mqish)

    res = utils.impute_missing_gaussian(
        rectangular_mat,
        method="mqish",
        scale=0.1,
        make_plot=False,
    )

    assert "legacy" not in calls
    assert calls["mqish"]["effective_width"] == 0.1
    assert res.isna().sum().sum() == 0
