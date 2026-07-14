from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import utils
from .gct_io import write_gctx
from .zscore import my_zscore


def correlation_contract() -> dict[str, Any]:
    """Return the metric and clustering choices shared by the CLI and plotter."""

    return {
        "metrics": ("l2", "l1", "pearson", "spearman"),
        "metric_aliases": {"euclidean": "l2", "manhattan": "l1"},
        "linkages": (
            "auto",
            "single",
            "complete",
            "average",
            "weighted",
            "centroid",
            "median",
            "ward.D2",
        ),
        "default_linkage": "auto",
        "automatic_linkages": {
            "l2": "ward.D2",
            "l1": "average",
            "pearson": "ward.D2",
            "spearman": "ward.D2",
        },
        "linkage_aliases": {"weighted": "mcquitty"},
    }


def normalize_metric(metric: str) -> str:
    value = str(metric).strip().lower()
    contract = correlation_contract()
    value = contract["metric_aliases"].get(value, value)
    metrics = contract["metrics"]
    if value not in metrics:
        allowed = ", ".join(metrics)
        raise ValueError(
            f"Unsupported correlation metric {metric!r}; choose {allowed}."
        )
    return value


def normalize_linkage(linkage: str) -> str:
    value = str(linkage).strip()
    linkages = correlation_contract()["linkages"]
    if value not in linkages:
        allowed = ", ".join(linkages)
        raise ValueError(
            f"Unsupported correlation linkage {linkage!r}; choose {allowed}."
        )
    return value


def resolve_linkage(metric: str, linkage: str) -> str:
    """Resolve the metric-aware automatic linkage to an explicit R method."""

    metric = normalize_metric(metric)
    linkage = normalize_linkage(linkage)
    if linkage == "auto":
        return str(correlation_contract()["automatic_linkages"][metric])
    return linkage


@dataclass(frozen=True)
class SampleMetricResult:
    """Displayed metric, clustering dissimilarity, and source overlap counts."""

    metric: str
    metric_matrix: pd.DataFrame
    clustering_distance_matrix: pd.DataFrame
    overlap_counts: pd.DataFrame


def compute_pairwise_counts(matrix: pd.DataFrame) -> pd.DataFrame:
    """Count jointly observed features for every pair of samples."""

    finite = np.isfinite(matrix.to_numpy(dtype=float))
    values = finite.astype(np.int64).T @ finite.astype(np.int64)
    counts = pd.DataFrame(
        values,
        index=matrix.columns.copy(),
        columns=matrix.columns.copy(),
    )
    counts.index.name = "sample"
    counts.columns.name = "sample"
    return counts


def exclude_correlation_samples(
    matrix: pd.DataFrame,
    sample_exclude: Iterable[str] | None,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Drop named samples in input order and return a compact selection record."""

    requested = list(dict.fromkeys(map(str, sample_exclude or ())))
    available = set(matrix.columns.astype(str))
    applied_columns = [
        sample for sample in matrix.columns if str(sample) in requested
    ]
    applied = [str(sample) for sample in applied_columns]
    unknown = [sample for sample in requested if sample not in available]
    selected = matrix.drop(columns=applied_columns) if applied_columns else matrix.copy()
    if selected.shape[1] < 2:
        raise ValueError(
            "Correlation requires at least two samples after --sample-exclude."
        )
    return selected, {
        "requested": requested,
        "applied": applied,
        "unknown": unknown,
    }


def flatten_metadata_fields(values: Iterable[str] | None) -> list[str]:
    """Normalize repeatable and colon-delimited metadata options."""

    result: list[str] = []
    for value in values or ():
        for field in str(value).split(":"):
            field = field.strip()
            if field and field not in result:
                result.append(field)
    return result


def prepare_logged_matrix(data_obj) -> pd.DataFrame:
    """Return the masked, log-shifted feature-by-sample matrix."""

    matrix = data_obj.areas_log_shifted.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    mask = getattr(data_obj, "mask", None)
    if mask is not None:
        mask = mask.copy()
        mask.index = mask.index.astype(str)
        mask.columns = mask.columns.astype(str)
        mask = mask.reindex(
            index=matrix.index, columns=matrix.columns, fill_value=False
        )
        matrix = matrix.mask(mask.astype(bool))

    if matrix.shape[1] < 2:
        raise ValueError("Correlation requires at least two samples.")
    if not np.isfinite(matrix.to_numpy(dtype=float)).any():
        raise ValueError("The masked log-expression matrix has no finite values.")

    matrix.index.name = matrix.index.name or "GeneID"
    return matrix


def _zscore_block(block: pd.DataFrame) -> pd.DataFrame:
    """Apply cluster2-style detection-aware z-scoring to feature rows."""

    result = block.apply(my_zscore, axis=1, remask=True, fillna=True)
    result.index = block.index.copy()
    result.columns = block.columns.copy()
    return result


def zscore_feature_rows(
    matrix: pd.DataFrame,
    *,
    metadata: pd.DataFrame | None = None,
    by: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Z-score each feature globally or within metadata-defined sample groups."""

    by_fields = flatten_metadata_fields(by)
    if not by_fields:
        return _zscore_block(matrix)
    if metadata is None:
        raise ValueError("Metadata is required when --z-score-by is used.")

    metadata = metadata.copy()
    metadata.index = metadata.index.astype(str)
    missing_fields = [field for field in by_fields if field not in metadata.columns]
    if missing_fields:
        raise ValueError(
            "Z-score metadata field(s) not found: " + ", ".join(missing_fields)
        )
    missing_samples = [
        sample for sample in matrix.columns if sample not in metadata.index
    ]
    if missing_samples:
        raise ValueError(
            "Metadata is missing sample(s): " + ", ".join(missing_samples[:10])
        )

    group_data = metadata.loc[matrix.columns, by_fields].copy()
    for field in by_fields:
        group_data[field] = utils.normalize_metadata_str_values(group_data[field])

    result = pd.DataFrame(np.nan, index=matrix.index, columns=matrix.columns)
    grouped = group_data.groupby(by_fields, sort=False, dropna=False, observed=True)
    for _key, rows in grouped:
        samples = list(rows.index)
        result.loc[:, samples] = _zscore_block(matrix.loc[:, samples])
    result.index.name = matrix.index.name
    return result


def compute_sample_metric(
    matrix: pd.DataFrame,
    metric: str,
    *,
    overlap_counts: pd.DataFrame | None = None,
) -> SampleMetricResult:
    """Compute the displayed metric and its clustering dissimilarity.

    L2 values are pairwise-complete RMS distances and L1 values are
    pairwise-complete mean absolute distances, so sample pairs with different
    overlap counts remain comparable.
    Pearson and Spearman values are displayed as correlation coefficients and
    are converted to chord distances ``sqrt(2 * (1 - r))`` for clustering.
    """

    metric = normalize_metric(metric)
    values = matrix.to_numpy(dtype=float)
    if overlap_counts is None:
        counts = compute_pairwise_counts(matrix)
    else:
        counts = overlap_counts.copy()
        counts.index = counts.index.astype(str)
        counts.columns = counts.columns.astype(str)
        sample_ids = list(matrix.columns.astype(str))
        if list(counts.index) != sample_ids or list(counts.columns) != sample_ids:
            raise ValueError(
                "Pairwise counts must have the same sample identifiers and ordering "
                "as the metric matrix."
            )

    sample_ids = matrix.columns.copy()
    if metric in {"pearson", "spearman"}:
        metric_matrix = matrix.corr(method=metric, min_periods=2)
        np.fill_diagonal(metric_matrix.values, 1.0)
        clipped = metric_matrix.clip(lower=-1.0, upper=1.0)
        clustering_distance_matrix = np.sqrt(2.0 * (1.0 - clipped))
        np.fill_diagonal(clustering_distance_matrix.values, 0.0)
    else:
        finite = np.isfinite(values)
        if not finite.any():
            raise ValueError(
                f"{metric.upper()} distance cannot be computed without finite values."
            )
        n_samples = values.shape[1]
        rms = np.full((n_samples, n_samples), np.nan, dtype=float)
        for left in range(n_samples):
            for right in range(left, n_samples):
                jointly_observed = finite[:, left] & finite[:, right]
                n_joint = int(jointly_observed.sum())
                if n_joint == 0:
                    continue
                delta = (
                    values[jointly_observed, left]
                    - values[jointly_observed, right]
                )
                value = (
                    float(np.sqrt(np.mean(np.square(delta))))
                    if metric == "l2"
                    else float(np.mean(np.abs(delta)))
                )
                rms[left, right] = value
                rms[right, left] = value
        metric_matrix = pd.DataFrame(rms, index=sample_ids, columns=sample_ids)
        clustering_distance_matrix = metric_matrix.copy()

    for frame in (metric_matrix, clustering_distance_matrix):
        frame.index.name = "sample"
        frame.columns.name = "sample"
    counts.index.name = "sample"
    counts.columns.name = "sample"
    return SampleMetricResult(
        metric=metric,
        metric_matrix=metric_matrix,
        clustering_distance_matrix=clustering_distance_matrix,
        overlap_counts=counts,
    )


def _string_keyed_mapping(values: Mapping[Any, Any] | None) -> dict[str, Any]:
    return {str(key): value for key, value in (values or {}).items()}


def _mapper_metadata_frame(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    if "GeneID" not in frame.columns:
        frame.insert(0, "GeneID", frame.index)
    frame["GeneID"] = frame["GeneID"].astype(str)
    return frame.drop_duplicates("GeneID", keep="first").set_index("GeneID")


def _missing_metadata(values: pd.Series) -> pd.Series:
    text = values.astype("string")
    return values.isna() | text.fillna("").str.strip().eq("")


def build_correlation_rdesc(data_obj, gene_ids: Sequence[Any]) -> pd.DataFrame:
    """Build feature metadata using the same mappers as tackle's GCT export."""

    from .containers import get_annotation_mapper, get_gene_mapper

    ids = pd.Index([str(gene_id) for gene_id in gene_ids], name="GeneID")
    rdesc = pd.DataFrame(index=ids)
    rdesc["GeneID"] = ids

    # The gene mapper is authoritative for gene identity. The broader annotation
    # mapper then contributes the optional localization and pathway descriptors.
    mapper_frames = (
        _mapper_metadata_frame(get_gene_mapper().df),
        _mapper_metadata_frame(get_annotation_mapper().df),
    )
    for source in mapper_frames:
        source = source.reindex(ids)
        for field in source.columns:
            if field == "GeneID":
                continue
            if field not in rdesc.columns:
                rdesc[field] = source[field]
            else:
                missing = _missing_metadata(rdesc[field])
                rdesc.loc[missing, field] = source.loc[missing, field]

    # Prefer the symbols and functional categories actually carried by this
    # Data object over a generic mapper value.
    preferred = {
        "GeneSymbol": _string_keyed_mapping(getattr(data_obj, "gid_symbol", None)),
        "FunCats": _string_keyed_mapping(getattr(data_obj, "gid_funcat_mapping", None)),
    }
    for field, mapping in preferred.items():
        if field not in rdesc.columns:
            rdesc[field] = ""
        carried = pd.Series(ids.map(mapping), index=ids)
        present = ~_missing_metadata(carried)
        rdesc.loc[present, field] = carried.loc[present]

    front = [
        field
        for field in ("GeneID", "GeneSymbol", "GeneDescription", "FunCats")
        if field in rdesc.columns
    ]
    remainder = [field for field in rdesc.columns if field not in front]
    rdesc = rdesc.loc[:, front + remainder]
    # GCT row descriptors are textual fields. Normalizing mixed annotation
    # columns here avoids lossy/fallback type guessing during pandas-to-R
    # conversion while retaining the exact displayed metadata values.
    for field in rdesc.columns:
        rdesc[field] = rdesc[field].where(rdesc[field].notna(), "").astype(str)
    rdesc["id"] = ids
    return rdesc


def build_correlation_cdesc(data_obj, samples: Sequence[Any]) -> pd.DataFrame:
    """Return every carried sample metadata field, aligned to matrix columns."""

    sample_ids = [str(sample) for sample in samples]
    cdesc = data_obj.col_metadata.copy()
    cdesc.index = cdesc.index.astype(str)
    missing = [sample for sample in sample_ids if sample not in cdesc.index]
    if missing:
        raise ValueError("Sample metadata missing for: " + ", ".join(missing[:10]))
    cdesc = cdesc.loc[sample_ids].copy()
    cdesc["id"] = sample_ids
    return cdesc


def write_correlation_input_gctx(
    data_obj,
    matrix: pd.DataFrame,
    out_path: str | Path,
) -> Path:
    """Write the exact transformed feature-by-sample input as float64 GCTX."""

    rdesc = build_correlation_rdesc(data_obj, matrix.index)
    cdesc = build_correlation_cdesc(data_obj, matrix.columns)
    return write_gctx(
        matrix,
        Path(out_path),
        row_metadata=rdesc,
        col_metadata=cdesc,
        matrix_dtype="float64",
        content_addressed=True,
    )


def write_correlation_sample_gctx(
    data_obj,
    matrix: pd.DataFrame,
    out_path: str | Path,
) -> Path:
    """Write a sample-by-sample matrix with original sample metadata on both axes."""

    row_ids = list(matrix.index.astype(str))
    col_ids = list(matrix.columns.astype(str))
    if row_ids != col_ids:
        raise ValueError(
            "Correlation summary GCTX matrices must be square with identical "
            "sample identifiers and ordering on both axes."
        )
    sample_metadata = build_correlation_cdesc(data_obj, col_ids)
    return write_gctx(
        matrix,
        Path(out_path),
        row_metadata=sample_metadata,
        col_metadata=sample_metadata,
        matrix_dtype="float64",
        content_addressed=True,
    )


# Compatibility name for callers from the first correlation implementation.
write_correlation_input_gct = write_correlation_input_gctx


def select_plot_metadata(
    metadata: pd.DataFrame,
    samples: Iterable[str],
    *,
    include: Iterable[str] | None = None,
    exclude: Iterable[str] | None = None,
    cut_by: Iterable[str] | None = None,
    max_categories: int = 20,
) -> tuple[pd.DataFrame, list[str]]:
    """Select and normalize metadata annotations while retaining all cut fields."""

    samples = [str(sample) for sample in samples]
    metadata = metadata.copy()
    metadata.index = metadata.index.astype(str)
    missing_samples = [sample for sample in samples if sample not in metadata.index]
    if missing_samples:
        raise ValueError(
            "Metadata is missing sample(s): " + ", ".join(missing_samples[:10])
        )
    metadata = metadata.loc[samples]

    include_fields = flatten_metadata_fields(include)
    exclude_fields = set(flatten_metadata_fields(exclude))
    cut_fields = flatten_metadata_fields(cut_by)
    requested = include_fields + [x for x in cut_fields if x not in include_fields]
    missing_fields = [field for field in requested if field not in metadata.columns]
    if missing_fields:
        raise ValueError("Metadata field(s) not found: " + ", ".join(missing_fields))

    if include_fields:
        fields = list(include_fields)
    else:
        technical = {"recno", "runno", "searchno", "label"}
        fields = []
        for field in metadata.columns:
            if field in technical or field in exclude_fields:
                continue
            series = metadata[field]
            nunique = int(series.nunique(dropna=True))
            if nunique <= 1:
                continue
            if pd.api.types.is_numeric_dtype(series):
                fields.append(field)
            elif nunique < len(samples) and nunique <= int(max_categories):
                fields.append(field)

    fields = [field for field in fields if field not in exclude_fields]
    for field in cut_fields:
        if field not in fields:
            fields.append(field)

    selected = metadata.loc[:, fields].copy() if fields else pd.DataFrame(index=samples)
    for field in selected.columns:
        if not pd.api.types.is_numeric_dtype(selected[field]):
            selected[field] = utils.normalize_metadata_str_values(selected[field])
    selected.index.name = "sample"
    return selected, cut_fields


def metadata_color_payload(
    metadata: pd.DataFrame,
    configured_colors: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, dict[str, str]]:
    payload: dict[str, dict[str, str]] = {}
    configured_colors = configured_colors or {}
    for field in metadata.columns:
        series = metadata[field]
        if pd.api.types.is_numeric_dtype(series):
            continue
        overrides = configured_colors.get(field)
        mapping = utils.get_color_mapping_with_defaults(
            series,
            overrides if isinstance(overrides, dict) else None,
        )
        if mapping:
            payload[str(field)] = {
                str(key): str(value) for key, value in mapping.items()
            }
    return payload


def plot_correlation_heatmap(
    metric_matrix: pd.DataFrame,
    clustering_distance_matrix: pd.DataFrame,
    overlap_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    outname: str,
    metric: str,
    linkage: str,
    z_score: bool,
    z_score_by: Iterable[str] | None,
    cut_by: Iterable[str] | None,
    cluster: bool,
    annotate: bool,
    file_fmts: Iterable[str],
    metadata_colors: Mapping[str, Mapping[str, str]] | None = None,
    fig_width: float | None = None,
    fig_height: float | None = None,
    png_res: int = 300,
) -> list[str]:
    """Render the precomputed sample metric with ComplexHeatmap."""

    import rpy2.robjects as robjects
    from rpy2.robjects import conversion, pandas2ri
    from rpy2.robjects.conversion import localconverter

    metric = normalize_metric(metric)
    linkage = resolve_linkage(metric, linkage)
    cut_fields = flatten_metadata_fields(cut_by)
    z_by_fields = flatten_metadata_fields(z_score_by)
    metadata_for_r = metadata.copy()
    metadata_for_r.insert(0, "sample", metadata_for_r.index.astype(str))
    colors_json = json.dumps(
        metadata_color_payload(metadata, metadata_colors), sort_keys=True
    )

    with localconverter(robjects.default_converter + pandas2ri.converter):
        metric_r = conversion.py2rpy(metric_matrix)
        distance_r = conversion.py2rpy(clustering_distance_matrix)
        overlap_r = conversion.py2rpy(overlap_counts)
        metadata_r = conversion.py2rpy(metadata_for_r)

    r_file = Path(__file__).resolve().parent / "R" / "correlationplot.R"
    robjects.r["source"](str(r_file))
    r_plot = robjects.r["correlation_heatmap"]

    z_label = "feature z-scored log expression" if z_score else "log expression"
    if z_by_fields:
        z_label += " within " + " + ".join(z_by_fields)
    if metric == "l2":
        measure_label = "RMS distance"
    elif metric == "l1":
        measure_label = "mean absolute distance"
    else:
        measure_label = "correlation"
    metric_label = metric.upper() if metric in {"l2", "l1"} else metric.title()
    title = f"{metric_label} sample {measure_label} ({z_label})"

    result = r_plot(
        metric_r,
        distance_matrix=distance_r,
        overlap_counts=overlap_r,
        metadata=metadata_r,
        metric=metric,
        linkage=linkage,
        cut_by=(
            robjects.vectors.StrVector(cut_fields) if cut_fields else robjects.NULL
        ),
        cluster=bool(cluster),
        annotate=bool(annotate),
        outname=os.path.abspath(outname),
        outfiletypes=robjects.vectors.StrVector(list(file_fmts)),
        metadata_colors_json=colors_json,
        title=title,
        fig_width=fig_width if fig_width is not None else robjects.NULL,
        fig_height=fig_height if fig_height is not None else robjects.NULL,
        png_res=int(png_res),
    )
    return [str(path) for path in result]
