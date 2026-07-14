from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from scipy.stats import f as f_distribution


@dataclass(frozen=True)
class PcaTestScope:
    """A named set of PCA score columns tested as one multivariate response."""

    name: str
    pcs: tuple[str, ...]
    plot_key: str | None = None
    selection: str = "explicit"


@dataclass(frozen=True)
class WelchJamesResult:
    statistic: float
    numerator_df: float
    denominator_df: float
    p_value: float
    status: str = "ok"
    message: str = ""


def _pc_number(column: object) -> int | None:
    match = re.fullmatch(r"PC(\d+)", str(column), flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def _scope_from_pcs(
    pcs: Sequence[str],
    *,
    name: str | None = None,
    selection: str = "explicit",
) -> PcaTestScope:
    normalized = tuple(f"PC{int(str(pc).removeprefix('PC').removeprefix('pc'))}" for pc in pcs)
    if len(normalized) < 2:
        raise ValueError("A PCA separation scope needs at least two components.")
    if len(set(normalized)) != len(normalized):
        raise ValueError("A PCA separation scope cannot repeat a component.")
    plot_key = None
    if len(normalized) == 2:
        plot_key = f"{normalized[0].lower()}_vs_{normalized[1].lower()}"
    return PcaTestScope(
        name=name or "_".join(normalized),
        pcs=normalized,
        plot_key=plot_key,
        selection=selection,
    )


def resolve_pca_test_scopes(
    requested: Iterable[str] | None,
    *,
    scores: pd.DataFrame,
    max_pc: int,
) -> list[PcaTestScope]:
    """Resolve displayed planes or explicit PC sets against available score columns.

    ``displayed`` expands to every two-dimensional plane produced by ``pca2`` up
    to ``max_pc``. ``leading`` (and its backwards-compatible alias ``all``)
    supplies all available components as candidates; analysis then selects the
    largest leading PC1..PCk block for which the omnibus and all required
    pairwise Welch-James tests are estimable.
    """

    available = {
        number: str(column)
        for column in scores.columns
        if (number := _pc_number(column)) is not None
    }
    if len(available) < 2:
        raise ValueError("PCA separation testing requires at least two score columns.")

    tokens = [str(value).strip() for value in (requested or ()) if str(value).strip()]
    if not tokens:
        tokens = ["displayed", "leading"]

    scopes: list[PcaTestScope] = []
    for token in tokens:
        lowered = token.lower()
        if lowered == "displayed":
            displayed = [number for number in sorted(available) if number <= int(max_pc)]
            if len(displayed) < 2:
                raise ValueError("--max-pc must be at least 2 for displayed-plane tests.")
            for left, right in combinations(displayed, 2):
                scopes.append(
                    _scope_from_pcs(
                        (available[left], available[right]),
                        selection="displayed",
                    )
                )
            continue

        if lowered in {"all", "leading"}:
            all_components = [available[number] for number in sorted(available)]
            scopes.append(
                _scope_from_pcs(
                    all_components,
                    name="leading_estimable_pcs",
                    selection="leading_estimable",
                )
            )
            continue

        parts = [part for part in re.split(r"[,+:]", token) if part.strip()]
        numbers: list[int] = []
        for part in parts:
            match = re.fullmatch(r"(?:PC)?(\d+)", part.strip(), flags=re.IGNORECASE)
            if not match:
                raise ValueError(
                    f"Invalid PCA test scope {token!r}; use displayed, leading, or e.g. 1,2."
                )
            numbers.append(int(match.group(1)))
        missing = [number for number in numbers if number not in available]
        if missing:
            raise ValueError(
                "Requested PCA component(s) are unavailable: "
                + ", ".join(f"PC{number}" for number in missing)
            )
        scopes.append(_scope_from_pcs([available[number] for number in numbers]))

    unique: list[PcaTestScope] = []
    seen: set[tuple[str, tuple[str, ...]]] = set()
    for scope in scopes:
        scope_kind = "leading" if scope.selection == "leading_estimable" else "fixed"
        key = (scope_kind, scope.pcs)
        if key in seen:
            continue
        seen.add(key)
        unique.append(scope)
    return unique


def euclidean_r_squared(values: np.ndarray, groups: Sequence[object]) -> float:
    """Return centroid R2 = 1 - SSwithin / SStotal for Euclidean coordinates."""

    matrix = np.asarray(values, dtype=float)
    labels = np.asarray(groups, dtype=object)
    if matrix.ndim != 2 or matrix.shape[0] != labels.shape[0]:
        raise ValueError("PCA values and grouping labels have incompatible dimensions.")
    if matrix.shape[0] == 0 or not np.isfinite(matrix).all():
        return np.nan

    grand = matrix.mean(axis=0)
    ss_total = float(np.square(matrix - grand).sum())
    if not np.isfinite(ss_total) or ss_total <= 0:
        return np.nan

    ss_within = 0.0
    for level in pd.unique(labels):
        subset = matrix[labels == level]
        centroid = subset.mean(axis=0)
        ss_within += float(np.square(subset - centroid).sum())
    return float(np.clip(1.0 - ss_within / ss_total, 0.0, 1.0))


def pairwise_centroid_geometry(
    values: np.ndarray,
    groups: Sequence[object],
) -> dict[str, object]:
    """Summarize two-group separation and within-group score-space radii."""

    matrix = np.asarray(values, dtype=float)
    labels = np.asarray(groups, dtype=object)
    if matrix.ndim != 2 or matrix.shape[0] != labels.shape[0]:
        raise ValueError("PCA values and grouping labels have incompatible dimensions.")
    levels = list(pd.unique(labels))
    empty = {
        "geometry_group_a": None,
        "geometry_group_b": None,
        "centroid_distance": np.nan,
        "rms_radius_a": np.nan,
        "rms_radius_b": np.nan,
        "pooled_rms_radius": np.nan,
        "standardized_separation": np.nan,
    }
    if len(levels) != 2 or matrix.shape[0] == 0 or not np.isfinite(matrix).all():
        return empty

    left, right = levels
    left_values = matrix[labels == left]
    right_values = matrix[labels == right]
    if left_values.size == 0 or right_values.size == 0:
        return empty
    left_centroid = left_values.mean(axis=0)
    right_centroid = right_values.mean(axis=0)
    centroid_distance = float(np.linalg.norm(left_centroid - right_centroid))
    radius_a = float(
        np.sqrt(np.square(left_values - left_centroid).sum(axis=1).mean())
    )
    radius_b = float(
        np.sqrt(np.square(right_values - right_centroid).sum(axis=1).mean())
    )
    pooled_radius = float(np.sqrt((radius_a**2 + radius_b**2) / 2.0))
    standardized = (
        centroid_distance / pooled_radius
        if np.isfinite(pooled_radius) and pooled_radius > 0
        else np.nan
    )
    return {
        "geometry_group_a": str(left),
        "geometry_group_b": str(right),
        "centroid_distance": centroid_distance,
        "rms_radius_a": radius_a,
        "rms_radius_b": radius_b,
        "pooled_rms_radius": pooled_radius,
        "standardized_separation": float(standardized),
    }


def _failed_welch_james(status: str, message: str) -> WelchJamesResult:
    return WelchJamesResult(
        statistic=np.nan,
        numerator_df=np.nan,
        denominator_df=np.nan,
        p_value=np.nan,
        status=status,
        message=message,
    )


def welch_james_test(values: np.ndarray, groups: Sequence[object]) -> WelchJamesResult:
    """Johansen's heteroscedastic Welch-James ADF test for group mean vectors.

    This implements the one-way, between-sample specialization of the general
    matrix statistic. It deliberately requires an estimable covariance contrast
    and never substitutes a pseudoinverse or regularized covariance silently.
    """

    matrix = np.asarray(values, dtype=float)
    labels = np.asarray(groups, dtype=object)
    if matrix.ndim != 2 or matrix.shape[0] != labels.shape[0]:
        raise ValueError("PCA values and grouping labels have incompatible dimensions.")
    if matrix.shape[1] < 1 or not np.isfinite(matrix).all():
        return _failed_welch_james(
            "invalid_values", "The selected PCA score matrix contains non-finite values."
        )

    levels = list(pd.unique(labels))
    if len(levels) < 2:
        return _failed_welch_james(
            "insufficient_groups", "At least two nonempty groups are required."
        )

    group_values = [matrix[labels == level] for level in levels]
    group_sizes = [int(group.shape[0]) for group in group_values]
    if any(size < 2 for size in group_sizes):
        return _failed_welch_james(
            "insufficient_group_size",
            "Every group needs at least two samples for a covariance estimate.",
        )

    dimensions = matrix.shape[1]
    means = np.concatenate([group.mean(axis=0) for group in group_values])
    covariance_blocks = []
    for group, size in zip(group_values, group_sizes):
        covariance = np.atleast_2d(np.cov(group, rowvar=False, ddof=1)).astype(float)
        if covariance.shape != (dimensions, dimensions):
            covariance = covariance.reshape(dimensions, dimensions)
        covariance_blocks.append(covariance / size)
    sigma = block_diag(*covariance_blocks)

    contrasts = np.zeros((len(levels) - 1, len(levels)), dtype=float)
    contrasts[:, 0] = 1.0
    for row in range(len(levels) - 1):
        contrasts[row, row + 1] = -1.0
    contrast = np.kron(contrasts, np.eye(dimensions))
    numerator_df = int(np.linalg.matrix_rank(contrast))
    contrast_covariance = contrast @ sigma @ contrast.T
    contrast_covariance = (contrast_covariance + contrast_covariance.T) / 2.0

    covariance_rank = int(np.linalg.matrix_rank(contrast_covariance))
    if covariance_rank < numerator_df:
        return _failed_welch_james(
            "singular_covariance",
            f"The covariance contrast has rank {covariance_rank}, below the required {numerator_df}.",
        )
    condition = float(np.linalg.cond(contrast_covariance))
    if not np.isfinite(condition) or condition > 1e12:
        return _failed_welch_james(
            "ill_conditioned_covariance",
            f"The covariance contrast condition number is {condition:.3g}.",
        )

    inverse_covariance = np.linalg.solve(
        contrast_covariance, np.eye(numerator_df, dtype=float)
    )
    hypothesis = contrast @ means
    twj = float(hypothesis.T @ inverse_covariance @ hypothesis)

    projection = sigma @ contrast.T @ inverse_covariance @ contrast
    correction_sum = 0.0
    block_size = dimensions
    for group_index, size in enumerate(group_sizes):
        selector = np.zeros_like(sigma)
        start = group_index * block_size
        selector[start : start + block_size, start : start + block_size] = np.eye(
            block_size
        )
        component = projection @ selector
        trace = float(np.trace(component))
        correction_sum += (trace**2 + float(np.trace(component @ component))) / (
            size - 1
        )
    correction = 0.5 * correction_sum
    if not np.isfinite(correction) or correction <= 0:
        return _failed_welch_james(
            "invalid_df_correction",
            "The Welch-James approximate-df correction is not positive and finite.",
        )

    denominator_df = numerator_df * (numerator_df + 2.0) / (3.0 * correction)
    divisor = numerator_df + 2.0 * correction - (
        6.0 * correction / (numerator_df + 2.0)
    )
    statistic = twj / divisor
    p_value = float(
        f_distribution.sf(statistic, numerator_df, denominator_df)
    )
    if not all(np.isfinite(value) for value in (statistic, denominator_df, p_value)):
        return _failed_welch_james(
            "nonfinite_result", "The Welch-James approximation produced a non-finite result."
        )
    return WelchJamesResult(
        statistic=float(statistic),
        numerator_df=float(numerator_df),
        denominator_df=float(denominator_df),
        p_value=p_value,
    )


def adjust_pvalues(values: Sequence[float], method: str = "holm") -> np.ndarray:
    """Adjust finite p-values without requiring statsmodels."""

    pvalues = np.asarray(values, dtype=float)
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    finite_positions = np.flatnonzero(np.isfinite(pvalues))
    if finite_positions.size == 0:
        return adjusted
    finite = pvalues[finite_positions]
    method_key = str(method).strip().lower().replace("-", "_")
    if method_key in {"none", "raw"}:
        result = finite
    elif method_key == "bonferroni":
        result = np.minimum(1.0, finite * len(finite))
    else:
        order = np.argsort(finite, kind="stable")
        ordered = finite[order]
        count = len(ordered)
        if method_key == "holm":
            result_ordered = np.maximum.accumulate(
                (count - np.arange(count)) * ordered
            )
        elif method_key == "hochberg":
            result_ordered = np.minimum.accumulate(
                ((count - np.arange(count)) * ordered)[::-1]
            )[::-1]
        elif method_key in {"bh", "fdr_bh"}:
            result_ordered = np.minimum.accumulate(
                (ordered * count / np.arange(1, count + 1))[::-1]
            )[::-1]
        else:
            raise ValueError(
                f"Unknown p-value adjustment {method!r}; use holm, hochberg, "
                "bonferroni, BH, or none."
            )
        result_ordered = np.clip(result_ordered, 0.0, 1.0)
        result = np.empty_like(result_ordered)
        result[order] = result_ordered
    adjusted[finite_positions] = result
    return adjusted


def _clean_group_series(series: pd.Series) -> pd.Series:
    result = series.copy().astype(object)
    missing = result.isna()
    as_text = result.astype(str).str.strip()
    missing |= as_text.str.lower().isin({"", "na", "nan", "<na>", "null"})
    result = as_text.mask(missing)
    return result


def _group_sizes_text(groups: pd.Series) -> str:
    counts = groups.value_counts(sort=False)
    return ";".join(f"{level}={int(count)}" for level, count in counts.items())


def _empty_pairwise_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "group_field",
            "scope",
            "pcs",
            "n_pcs",
            "explained_variance_pct",
            "group_a",
            "group_b",
            "n_a",
            "n_b",
            "centroid_distance",
            "rms_radius_a",
            "rms_radius_b",
            "pooled_rms_radius",
            "standardized_separation",
            "r2",
            "welch_james_f",
            "numerator_df",
            "denominator_df",
            "p_value",
            "p_adjust_method",
            "p_adj",
            "p_adj_all_scopes",
            "status",
            "message",
            "method",
        ]
    )


def _select_scope_data(
    score_frame: pd.DataFrame,
    group_series: pd.Series,
    pcs: Sequence[str],
) -> tuple[pd.DataFrame, pd.Series, list[object]]:
    selected = score_frame.loc[:, list(pcs)].apply(pd.to_numeric, errors="coerce")
    keep = group_series.notna() & np.isfinite(selected).all(axis=1)
    selected = selected.loc[keep]
    groups = group_series.loc[keep]
    return selected, groups, list(pd.unique(groups))


def _largest_common_estimable_leading_scope(
    scope: PcaTestScope,
    *,
    score_frame: pd.DataFrame,
    group_series: pd.Series,
) -> PcaTestScope | None:
    """Choose the largest PC1..PCk block valid globally and for every pair."""

    for component_count in range(len(scope.pcs), 0, -1):
        pcs = scope.pcs[:component_count]
        selected, groups, levels = _select_scope_data(
            score_frame, group_series, pcs
        )
        values = selected.to_numpy(dtype=float)
        if welch_james_test(values, groups.to_numpy(dtype=object)).status != "ok":
            continue
        pairwise_ok = True
        if len(levels) > 2:
            for left, right in combinations(levels, 2):
                pair_keep = groups.isin([left, right])
                pair_result = welch_james_test(
                    selected.loc[pair_keep].to_numpy(dtype=float),
                    groups.loc[pair_keep].to_numpy(dtype=object),
                )
                if pair_result.status != "ok":
                    pairwise_ok = False
                    break
        if pairwise_ok:
            return PcaTestScope(
                name="leading_estimable_pcs",
                pcs=tuple(pcs),
                plot_key=None,
                selection="leading_estimable",
            )
    return None


def analyze_pca_separation(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    group_fields: Iterable[str],
    scopes: Sequence[PcaTestScope],
    p_adjust_method: str = "holm",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run plane/scope omnibus tests and post-hoc pairwise group tests."""

    score_frame = scores.copy()
    score_frame.index = score_frame.index.astype(str)
    metadata_frame = metadata.copy()
    metadata_frame.index = metadata_frame.index.astype(str)
    metadata_frame = metadata_frame.reindex(score_frame.index)
    component_variances = score_frame.apply(pd.to_numeric, errors="coerce").var(
        axis=0, ddof=1
    )
    total_variance = float(component_variances.sum())

    def explained_variance_pct(pcs: Sequence[str]) -> float:
        if not np.isfinite(total_variance) or total_variance <= 0:
            return np.nan
        return float(component_variances.loc[list(pcs)].sum() / total_variance * 100.0)

    omnibus_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []
    for field in group_fields:
        if field not in metadata_frame.columns:
            raise ValueError(f"PCA test metadata field {field!r} was not found.")
        group_series = _clean_group_series(metadata_frame[field])

        field_omnibus_start = len(omnibus_rows)
        field_pairwise_start = len(pairwise_rows)
        ordered_scopes = [
            scope for scope in scopes if scope.selection != "leading_estimable"
        ] + [scope for scope in scopes if scope.selection == "leading_estimable"]
        for requested_scope in ordered_scopes:
            scope = requested_scope
            if scope.selection == "leading_estimable":
                resolved = _largest_common_estimable_leading_scope(
                    scope,
                    score_frame=score_frame,
                    group_series=group_series,
                )
                if resolved is None:
                    continue
                scope = resolved
                pcs_text = ",".join(scope.pcs)
                duplicate_rows = [
                    row
                    for row in omnibus_rows[field_omnibus_start:]
                    if row["pcs"] == pcs_text
                ]
                if duplicate_rows:
                    continue

            selected, groups, levels = _select_scope_data(
                score_frame, group_series, scope.pcs
            )
            values = selected.to_numpy(dtype=float)
            result = welch_james_test(values, groups.to_numpy(dtype=object))
            row = {
                "group_field": str(field),
                "scope": scope.name,
                "pcs": ",".join(scope.pcs),
                "n_pcs": len(scope.pcs),
                "explained_variance_pct": explained_variance_pct(scope.pcs),
                "n_samples": int(len(selected)),
                "n_groups": int(len(levels)),
                "group_sizes": _group_sizes_text(groups),
                "r2": euclidean_r_squared(values, groups.to_numpy(dtype=object)),
                "welch_james_f": result.statistic,
                "numerator_df": result.numerator_df,
                "denominator_df": result.denominator_df,
                "p_value": result.p_value,
                "p_adjust_method": p_adjust_method,
                "p_adj": np.nan,
                "status": result.status,
                "message": result.message,
                "method": "Johansen Welch-James ADF",
            }
            omnibus_rows.append(row)

            if len(levels) < 2:
                continue
            scope_pairwise_start = len(pairwise_rows)
            for left, right in combinations(levels, 2):
                pair_keep = groups.isin([left, right])
                pair_values = selected.loc[pair_keep].to_numpy(dtype=float)
                pair_groups = groups.loc[pair_keep]
                pair_result = welch_james_test(
                    pair_values, pair_groups.to_numpy(dtype=object)
                )
                pair_geometry = pairwise_centroid_geometry(
                    pair_values, pair_groups.to_numpy(dtype=object)
                )
                counts = pair_groups.value_counts(sort=False)
                pairwise_rows.append(
                    {
                        "group_field": str(field),
                        "scope": scope.name,
                        "pcs": ",".join(scope.pcs),
                        "n_pcs": len(scope.pcs),
                        "explained_variance_pct": explained_variance_pct(scope.pcs),
                        "group_a": str(left),
                        "group_b": str(right),
                        "n_a": int(counts.get(left, 0)),
                        "n_b": int(counts.get(right, 0)),
                        "centroid_distance": pair_geometry["centroid_distance"],
                        "rms_radius_a": pair_geometry["rms_radius_a"],
                        "rms_radius_b": pair_geometry["rms_radius_b"],
                        "pooled_rms_radius": pair_geometry["pooled_rms_radius"],
                        "standardized_separation": pair_geometry[
                            "standardized_separation"
                        ],
                        "r2": euclidean_r_squared(
                            pair_values, pair_groups.to_numpy(dtype=object)
                        ),
                        "welch_james_f": pair_result.statistic,
                        "numerator_df": pair_result.numerator_df,
                        "denominator_df": pair_result.denominator_df,
                        "p_value": pair_result.p_value,
                        "p_adjust_method": p_adjust_method,
                        "p_adj": np.nan,
                        "p_adj_all_scopes": np.nan,
                        "status": pair_result.status,
                        "message": pair_result.message,
                        "method": "Johansen Welch-James ADF",
                    }
                )
            scope_slice = slice(scope_pairwise_start, len(pairwise_rows))
            adjusted = adjust_pvalues(
                [row["p_value"] for row in pairwise_rows[scope_slice]],
                p_adjust_method,
            )
            for pair_row, p_adj in zip(pairwise_rows[scope_slice], adjusted):
                pair_row["p_adj"] = p_adj

        omnibus_slice = slice(field_omnibus_start, len(omnibus_rows))
        omnibus_adjusted = adjust_pvalues(
            [row["p_value"] for row in omnibus_rows[omnibus_slice]],
            p_adjust_method,
        )
        for omnibus_row, p_adj in zip(omnibus_rows[omnibus_slice], omnibus_adjusted):
            omnibus_row["p_adj"] = p_adj

        pairwise_slice = slice(field_pairwise_start, len(pairwise_rows))
        all_pairwise_adjusted = adjust_pvalues(
            [row["p_value"] for row in pairwise_rows[pairwise_slice]],
            p_adjust_method,
        )
        for pair_row, p_adj in zip(
            pairwise_rows[pairwise_slice], all_pairwise_adjusted
        ):
            pair_row["p_adj_all_scopes"] = p_adj

    pairwise_frame = (
        pd.DataFrame(pairwise_rows) if pairwise_rows else _empty_pairwise_frame()
    )
    return pd.DataFrame(omnibus_rows), pairwise_frame


def _p_adjustment_caption_label(method: object) -> str:
    key = str(method).strip().lower().replace("-", "_")
    labels = {
        "none": "p",
        "raw": "p",
        "holm": "Holm-adjusted p",
        "hochberg": "Hochberg-adjusted p",
        "bonferroni": "Bonferroni-adjusted p",
        "bh": "BH-adjusted p",
        "fdr_bh": "BH-adjusted p",
    }
    return labels.get(key, f"{str(method).strip()}-adjusted p")


def format_pca_test_caption(
    rows: pd.DataFrame,
    pairwise_rows: pd.DataFrame | None = None,
) -> str:
    """Format compact plot-caption lines for one displayed PCA plane."""

    lines: list[str] = []
    for row in rows.to_dict(orient="records"):
        r2 = row.get("r2")
        prefix = f"{row['group_field']}: R²={r2:.3f}" if np.isfinite(r2) else f"{row['group_field']}: R²=NA"
        if row.get("status") == "ok":
            p_label = _p_adjustment_caption_label(row.get("p_adjust_method", "none"))
            line = (
                prefix
                + "; WJ "
                + f"F*({row['numerator_df']:.0f}, {row['denominator_df']:.2f})="
                + f"{row['welch_james_f']:.2f}; {p_label}={row['p_adj']:.3g}"
            )
        else:
            line = prefix + f"; WJ not estimable ({row.get('status')})"
        geometry_row: dict[str, object] | None = None
        if pairwise_rows is not None and not pairwise_rows.empty:
            candidates = pairwise_rows
            for key in ("group_field", "scope"):
                if key in candidates.columns and key in row:
                    candidates = candidates.loc[candidates[key] == row[key]]
            if len(candidates) == 1:
                geometry_row = candidates.iloc[0].to_dict()
        if geometry_row is not None and np.isfinite(
            geometry_row.get("standardized_separation", np.nan)
        ):
            line += (
                f"\nCentroid distance = {geometry_row['centroid_distance']:.2f}"
                f"\nRMS radii ({geometry_row['group_a']}, {geometry_row['group_b']}) = "
                f"{geometry_row['rms_radius_a']:.2f}, {geometry_row['rms_radius_b']:.2f}"
                "\nStandardized separation = "
                f"{geometry_row['standardized_separation']:.2f}"
            )
        lines.append(line)
    return "\n".join(lines)
