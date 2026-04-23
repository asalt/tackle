from types import SimpleNamespace

import pandas as pd

import tackle.containers as containers
from tackle.containers import Data
from tackle.impute_cache import (
    compute_imputation_cache_key,
    load_imputation_cache,
    save_imputation_cache,
)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {"S1": [1.0, None, 3.0], "S2": [4.0, 5.0, None]},
        index=["101", "202", "303"],
    )


def test_imputation_cache_key_distinguishes_methods():
    frame = _frame()
    legacy = compute_imputation_cache_key(
        frame,
        backend="gaussian",
        gaussian_method="legacy",
        params={"downshift": 1.8, "scale": 0.8},
    )
    mqish = compute_imputation_cache_key(
        frame,
        backend="gaussian",
        gaussian_method="mqish",
        params={"downshift": 1.8, "effective_width": 0.3},
    )
    lupine = compute_imputation_cache_key(
        frame,
        backend="lupine",
        lupine_mode="local",
        params={"mode": "local"},
    )

    assert legacy != mqish
    assert legacy != lupine
    assert mqish != lupine


def test_imputation_cache_roundtrip(tmp_path):
    frame = _frame()
    key = "abc123"

    save_imputation_cache(
        tmp_path,
        key,
        frame,
        metadata={"backend": "gaussian", "gaussian_method": "legacy"},
        overwrite=True,
    )

    loaded = load_imputation_cache(tmp_path, key)
    assert loaded is not None
    assert loaded.equals(frame)


def test_data_impute_missing_reuses_cache(monkeypatch, tmp_path):
    calls = {"gaussian": 0}

    def fake_gaussian(frame, **kwargs):
        calls["gaussian"] += 1
        return frame.fillna(float(calls["gaussian"]))

    monkeypatch.setattr(containers, "impute_missing_gaussian", fake_gaussian)

    dummy = SimpleNamespace(
        imputation_backend="gaussian",
        gaussian_method="legacy",
        lupine_mode="local",
        cache_impute=True,
        cache_overwrite=False,
        impute_missing_values=True,
        outpath=str(tmp_path),
    )

    frame = _frame()
    first = Data.impute_missing(dummy, frame, downshift=1.8, scale=0.8, make_plot=False)
    second = Data.impute_missing(dummy, frame, downshift=1.8, scale=0.8, make_plot=False)

    assert calls["gaussian"] == 1
    assert first.equals(second)


def test_data_impute_missing_overwrites_cache(monkeypatch, tmp_path):
    calls = {"gaussian": 0}

    def fake_gaussian(frame, **kwargs):
        calls["gaussian"] += 1
        return frame.fillna(float(calls["gaussian"]))

    monkeypatch.setattr(containers, "impute_missing_gaussian", fake_gaussian)

    dummy = SimpleNamespace(
        imputation_backend="gaussian",
        gaussian_method="legacy",
        lupine_mode="local",
        cache_impute=True,
        cache_overwrite=True,
        impute_missing_values=True,
        outpath=str(tmp_path),
    )

    frame = _frame()
    first = Data.impute_missing(dummy, frame, downshift=1.8, scale=0.8, make_plot=False)
    second = Data.impute_missing(dummy, frame, downshift=1.8, scale=0.8, make_plot=False)

    assert calls["gaussian"] == 2
    assert not first.equals(second)
