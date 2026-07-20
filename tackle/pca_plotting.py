from __future__ import annotations

import math
from collections.abc import Sequence

import pandas as pd


def _round_up_quarter(value: float) -> float:
    return math.ceil(float(value) * 4.0) / 4.0


def resolve_pca2_figsize(
    requested: Sequence[float] | None,
    *,
    sample_metadata: pd.DataFrame,
    color: str | None,
    marker: str | None,
) -> tuple[float, float]:
    """Resolve a compact PCA size, expanding only for a demanding legend.

    Six by seven inches leaves a useful portrait-oriented plotting panel with a
    normal right-side legend.  Extra legend entries primarily require height;
    unusually long labels and a second legend require some additional width.
    An explicit ``--figsize`` always wins unchanged.
    """

    if requested is not None:
        width, height = requested
        return float(width), float(height)

    legend_fields = []
    for field in (color, marker):
        if field and field not in legend_fields and field in sample_metadata.columns:
            legend_fields.append(field)

    level_labels: list[str] = []
    for field in legend_fields:
        values = sample_metadata[field].dropna().astype(str)
        level_labels.extend(dict.fromkeys(values).keys())

    total_levels = len(level_labels)
    longest_label = max((len(value) for value in level_labels), default=0)

    width = 6.0
    height = 7.0
    if len(legend_fields) > 1:
        width += 0.5
    if longest_label > 18:
        width += min(3.0, 0.08 * (longest_label - 18))
    if total_levels > 12:
        height += min(7.0, 0.2 * (total_levels - 12))

    return _round_up_quarter(width), _round_up_quarter(height)
