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


def test_legacy_imputation_reports_final_sd_after_averaging_draws():
    rng = np.random.default_rng(7)
    frame = pd.DataFrame(
        {
            "observed": rng.normal(size=4000),
            "missing": np.nan,
        }
    )
    missing = frame.isna()
    observed = frame.stack().dropna()

    imputed = utils.impute_missing_old(
        frame,
        downshift=1.8,
        scale=0.8,
        n_draws=8,
        make_plot=False,
    )
    summary = utils.summarize_imputation_distribution(imputed, observed, missing)

    expected_relative_sd = 0.8 / np.sqrt(8)
    assert summary["imputed_count"] == 4000
    assert summary["final_imputed_sd"] == pytest.approx(
        imputed.loc[missing["missing"], "missing"].std()
    )
    assert summary["final_imputed_sd_fraction_of_observed"] == pytest.approx(
        expected_relative_sd, abs=0.02
    )


def test_imputation_plot_reports_realized_sd_instead_of_draw_scale():
    frame = pd.DataFrame(
        {"S1": [1.0, 3.0], "S2": [10.0, 14.0]},
        index=["g1", "g2"],
    )
    missing = pd.DataFrame(
        {"S1": [False, False], "S2": [True, True]},
        index=frame.index,
    )

    summary = utils.plot_imputed(
        frame,
        frame["S1"],
        missing,
        downshift=1.8,
        scale=0.8,
    )
    title = utils.plt.gca().get_title()
    utils.plt.close(utils.plt.gcf())

    assert summary["final_imputed_sd"] == pytest.approx(np.sqrt(8))
    assert "final imputed SD: 2.83" in title
    assert "scale" not in title
    assert utils.imputation_distribution_plot_stem(1.8, np.sqrt(8)) == (
        "distribution_ds_1.8_imputed_sd_2.83"
    )
