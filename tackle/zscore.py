from __future__ import annotations

import numpy as np
import pandas as pd


def my_zscore(
    values: pd.Series,
    minval: float | None = None,
    remask: bool = True,
    fillna: bool = True,
) -> pd.Series:
    """Apply tackle's detection-aware feature z-score.

    Nondetections temporarily receive ``observed minimum - observed SD`` so
    sparse detections contribute relative to absence during scaling. A unit
    offset handles constant or singleton detections; its magnitude cancels
    during standardization. The original nonfinite mask is restored by default.
    """

    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    missing = pd.Series(
        ~np.isfinite(numeric.to_numpy()),
        index=numeric.index,
    )
    numeric = numeric.mask(missing)
    observed = numeric.loc[~missing]
    if observed.empty:
        return pd.Series(np.nan, index=numeric.index, name=numeric.name, dtype=float)

    if minval is None or pd.isna(minval):
        observed_sd = float(observed.std(ddof=1))
        offset = (
            observed_sd if np.isfinite(observed_sd) and observed_sd > 0 else 1.0
        )
        minval = float(observed.min()) - offset
    elif not np.isfinite(float(minval)):
        minval = 0.0

    scaled_input = numeric.fillna(float(minval)) if fillna else numeric
    scaled_sd = float(scaled_input.std(ddof=1))
    if not np.isfinite(scaled_sd) or scaled_sd <= 0:
        result = pd.Series(0.0, index=numeric.index, name=numeric.name)
        result.loc[scaled_input.isna()] = np.nan
    else:
        result = (scaled_input - float(scaled_input.mean())) / scaled_sd
    if remask:
        result.loc[missing] = np.nan
    return result
