import numpy as np
import pytest

import pandas as pd

from tackle.utils import impute_missing_lupine


@pytest.fixture
def small_rectangular_mat():
    # Small matrix with some NaNs to keep Lupine runtime reasonable
    rng = np.random.default_rng(42)
    mat = rng.normal(size=(6, 4))
    mask = rng.random(size=mat.shape) < 0.3
    mat[mask] = np.nan
    return pd.DataFrame(mat)


@pytest.mark.skipif(
    pytest.importorskip("lupine", reason="lupine package is required for this test")
    is None,
    reason="lupine package not available",
)
def test_impute_missing_lupine_basic(small_rectangular_mat):
    df = small_rectangular_mat
    # Use a single model on CPU to keep tests fast and deterministic-ish
    res = impute_missing_lupine(df, n_models=1, device="cpu", biased=False)

    # Shape is preserved
    assert res.shape == df.shape
    # No NaNs remain after imputation
    assert not res.isna().any().any()


@pytest.mark.skipif(
    pytest.importorskip("lupine", reason="lupine package is required for this test")
    is None,
    reason="lupine package not available",
)
def test_impute_missing_lupine_joint_mode_falls_back(small_rectangular_mat):
    df = small_rectangular_mat
    res = impute_missing_lupine(df, n_models=1, device="cpu", biased=False, mode="joint")
    assert res.shape == df.shape
    assert not res.isna().any().any()
