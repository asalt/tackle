"""Local, QUAL-only PSM quality-control summaries for tackle metrics."""

from __future__ import annotations

import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping
from warnings import warn

import pandas as pd


class PsmQcDependencyError(RuntimeError):
    """Raised when the optional PSM QC dependency is unavailable."""


@dataclass(frozen=True)
class PsmQcArtifacts:
    summary_tsv: Path
    ids_over_time_tsv: Path
    charge_distribution_tsv: Path
    missed_cleavage_distribution_tsv: Path
    plot_paths: tuple[Path, ...]


def _require_polars():
    try:
        import polars as pl
    except ImportError as exc:  # pragma: no cover - depends on test environment
        raise PsmQcDependencyError(
            "PSM QC requires Polars; install it with `pip install 'tackle[psm-qc]'`."
        ) from exc
    return pl


def require_psm_qc_dependencies() -> None:
    """Fail early with the installation hint when PSM QC is requested."""
    _require_polars()


def _collect_streaming(lazy_frame):
    try:
        return lazy_frame.collect(engine="streaming")
    except TypeError as exc:
        if "engine" not in str(exc):
            raise
        return lazy_frame.collect(streaming=True)


def _selected_columns() -> tuple[str, ...]:
    return (
        "SpectrumFile",
        "FirstScan",
        "Precursor.Id",
        "RTmin",
        "Sequence",
        "SequenceModi",
        "Charge",
        "MissedCleavages",
        "q_value",
        "PEP",
        "DeltaMassPPM",
        "PrecursorArea",
        "PSM_UseFLAG",
        "oriFLAG",
    )


def _summary_columns() -> tuple[str, ...]:
    return (
        "recno",
        "runno",
        "searchno",
        "run_id",
        "config_sections",
        "source_file",
        "status",
        "message",
        "psm_use_filter",
        "ori_filter",
        "identity_basis",
        "rows_raw",
        "rows_psm_use",
        "rows_canonical",
        "unique_sequences",
        "unique_modified_sequences",
        "rt_min",
        "rt_max",
        "q_value_median",
        "q_value_p95",
        "pep_median",
        "pep_p95",
        "delta_mass_ppm_median",
        "delta_mass_ppm_abs_p95",
        "precursor_area_median",
        "precursor_area_p95",
    )


def _long_columns(value_name: str) -> tuple[str, ...]:
    return (
        "recno",
        "runno",
        "searchno",
        "run_id",
        "config_sections",
        value_name,
        "count",
        "fraction",
    )


def _canonical_number(value, field: str) -> str:
    if value is None or str(value).strip() == "":
        raise ValueError(f"missing {field}")
    number = float(str(value).strip())
    if not number.is_integer() or number <= 0:
        raise ValueError(f"invalid {field}: {value!r}")
    return str(int(number))


def _config_triplets(config: Mapping) -> tuple[OrderedDict, list[dict]]:
    grouped = OrderedDict()
    errors = []
    for section, record in config.items():
        if str(section).startswith("__"):
            continue
        try:
            triplet = tuple(
                _canonical_number(
                    record.get(field, 1 if field != "recno" else None), field
                )
                for field in ("recno", "runno", "searchno")
            )
        except (TypeError, ValueError) as exc:
            errors.append(
                {
                    "recno": record.get("recno"),
                    "runno": record.get("runno", 1),
                    "searchno": record.get("searchno", 1),
                    "run_id": "",
                    "config_sections": str(section),
                    "source_file": "",
                    "status": "invalid_config",
                    "message": str(exc),
                }
            )
            continue
        grouped.setdefault(triplet, []).append(str(section))
    return grouped, errors


def _qual_file_index(data_dir: Path) -> list[Path]:
    return sorted(path.resolve() for path in data_dir.rglob("*_psms_QUAL.tsv"))


def _resolve_qual_file(index: Iterable[Path], triplet: tuple[str, str, str]):
    prefix = "_".join(triplet) + "_"
    matches = [
        path
        for path in index
        if path.name.startswith(prefix) and path.name.endswith("psms_QUAL.tsv")
    ]
    if len(matches) == 1:
        return matches[0], "ok", ""
    if not matches:
        return None, "missing", f"No local QUAL file matched {prefix}*psms_QUAL.tsv"
    joined = "; ".join(str(path) for path in matches)
    return None, "ambiguous", f"Multiple local QUAL files matched: {joined}"


def _numeric_series(frame, column: str):
    pl = _require_polars()
    if column not in frame.columns:
        return None
    return frame.get_column(column).cast(pl.Float64, strict=False).drop_nulls()


def _stat(frame, column: str, operation: str):
    series = _numeric_series(frame, column)
    if series is None or len(series) == 0:
        return None
    if operation == "median":
        return float(series.median())
    if operation == "p95":
        return float(series.quantile(0.95, interpolation="linear"))
    if operation == "abs_p95":
        return float(series.abs().quantile(0.95, interpolation="linear"))
    if operation == "min":
        return float(series.min())
    if operation == "max":
        return float(series.max())
    raise ValueError(f"Unknown summary operation: {operation}")


def _n_unique(frame, column: str):
    if column not in frame.columns:
        return None
    return int(frame.get_column(column).drop_nulls().n_unique())


def _filter_flag(frame, column: str):
    pl = _require_polars()
    return frame.filter(pl.col(column).cast(pl.Float64, strict=False) == 1)


def _canonicalize_psms(frame):
    if "oriFLAG" in frame.columns:
        return _filter_flag(frame, "oriFLAG"), "applied", "oriFLAG"

    candidates = (
        ("SpectrumFile", "FirstScan"),
        ("SpectrumFile", "Precursor.Id"),
        ("SpectrumFile", "SequenceModi", "Charge", "RTmin"),
    )
    for columns in candidates:
        if all(column in frame.columns for column in columns):
            return (
                frame.unique(subset=list(columns), maintain_order=True),
                "not_available",
                "+".join(columns),
            )
    return frame, "not_available", "row"


def _run_fields(triplet, sections) -> dict:
    recno, runno, searchno = triplet
    return {
        "recno": recno,
        "runno": runno,
        "searchno": searchno,
        "run_id": "_".join(triplet),
        "config_sections": ";".join(sections),
    }


def _ids_over_time(frame, run_fields: dict):
    pl = _require_polars()
    columns = (*run_fields.keys(), "rt_minute", "ids")
    if "RTmin" not in frame.columns:
        return pd.DataFrame(columns=columns)
    result = (
        frame.with_columns(
            pl.col("RTmin")
            .cast(pl.Float64, strict=False)
            .floor()
            .cast(pl.Int64)
            .alias("rt_minute")
        )
        .drop_nulls("rt_minute")
        .group_by("rt_minute")
        .agg(pl.len().alias("ids"))
        .sort("rt_minute")
    )
    result = pd.DataFrame(result.to_dicts())
    for column, value in reversed(tuple(run_fields.items())):
        result.insert(0, column, value)
    return result.loc[:, columns]


def _distribution(frame, column: str, value_name: str, run_fields: dict):
    pl = _require_polars()
    columns = _long_columns(value_name)
    if column not in frame.columns:
        return pd.DataFrame(columns=columns)
    result = (
        frame.with_columns(
            pl.col(column)
            .cast(pl.Float64, strict=False)
            .cast(pl.Int64, strict=False)
            .alias(value_name)
        )
        .drop_nulls(value_name)
        .group_by(value_name)
        .agg(pl.len().alias("count"))
        .sort(value_name)
        .with_columns((pl.col("count") / pl.col("count").sum()).alias("fraction"))
    )
    result = pd.DataFrame(result.to_dicts())
    for column_name, value in reversed(tuple(run_fields.items())):
        result.insert(0, column_name, value)
    return result.loc[:, columns]


def _summarize_qual(path: Path, triplet, sections):
    pl = _require_polars()
    run_fields = _run_fields(triplet, sections)
    lazy = pl.scan_csv(
        path,
        separator="\t",
        infer_schema=False,
        null_values=["", "NA", "NaN"],
    )
    available = set(lazy.collect_schema().names())
    selected = [column for column in _selected_columns() if column in available]
    if not selected:
        raise ValueError("QUAL file contains none of the supported PSM QC columns")
    frame = _collect_streaming(lazy.select(selected))

    rows_raw = frame.height
    if "PSM_UseFLAG" in frame.columns:
        frame = _filter_flag(frame, "PSM_UseFLAG")
        psm_use_filter = "applied"
    else:
        psm_use_filter = "not_available"
    rows_psm_use = frame.height

    frame, ori_filter, identity_basis = _canonicalize_psms(frame)
    summary = {
        **run_fields,
        "source_file": str(path),
        "status": "ok",
        "message": "",
        "psm_use_filter": psm_use_filter,
        "ori_filter": ori_filter,
        "identity_basis": identity_basis,
        "rows_raw": rows_raw,
        "rows_psm_use": rows_psm_use,
        "rows_canonical": frame.height,
        "unique_sequences": _n_unique(frame, "Sequence"),
        "unique_modified_sequences": _n_unique(frame, "SequenceModi"),
        "rt_min": _stat(frame, "RTmin", "min"),
        "rt_max": _stat(frame, "RTmin", "max"),
        "q_value_median": _stat(frame, "q_value", "median"),
        "q_value_p95": _stat(frame, "q_value", "p95"),
        "pep_median": _stat(frame, "PEP", "median"),
        "pep_p95": _stat(frame, "PEP", "p95"),
        "delta_mass_ppm_median": _stat(frame, "DeltaMassPPM", "median"),
        "delta_mass_ppm_abs_p95": _stat(frame, "DeltaMassPPM", "abs_p95"),
        "precursor_area_median": _stat(frame, "PrecursorArea", "median"),
        "precursor_area_p95": _stat(frame, "PrecursorArea", "p95"),
    }
    return (
        summary,
        _ids_over_time(frame, run_fields),
        _distribution(frame, "Charge", "charge", run_fields),
        _distribution(frame, "MissedCleavages", "missed_cleavages", run_fields),
    )


def _with_summary_columns(frame: pd.DataFrame) -> pd.DataFrame:
    for column in _summary_columns():
        if column not in frame.columns:
            frame[column] = None
    return frame.loc[:, _summary_columns()]


def _concat_frames(frames: list[pd.DataFrame], columns) -> pd.DataFrame:
    nonempty = [frame for frame in frames if not frame.empty]
    if not nonempty:
        return pd.DataFrame(columns=columns)
    return pd.concat(nonempty, ignore_index=True).loc[:, columns]


def _normalized_suffix(file_fmt: str) -> str:
    return file_fmt if str(file_fmt).startswith(".") else f".{file_fmt}"


def _save_figure(fig, output_root: Path, file_fmts, png_res: int) -> list[Path]:
    import matplotlib.pyplot as plt

    paths = []
    for file_fmt in file_fmts:
        suffix = _normalized_suffix(file_fmt)
        path = Path(f"{output_root}{suffix}")
        save_kws = {"bbox_inches": "tight"}
        if suffix.lower() == ".png":
            save_kws["dpi"] = png_res
        fig.savefig(path, **save_kws)
        paths.append(path)
    plt.close(fig)
    return paths


def _plot_ids_over_time(frame: pd.DataFrame, output_root, file_fmts, png_res):
    if frame.empty:
        return []
    import matplotlib.pyplot as plt

    run_ids = frame["run_id"].drop_duplicates().tolist()
    ncols = min(3, len(run_ids))
    nrows = math.ceil(len(run_ids) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.0 * ncols, 2.8 * nrows),
        squeeze=False,
        sharey=False,
    )
    for ax, run_id in zip(axes.flat, run_ids):
        subset = frame.loc[frame["run_id"] == run_id]
        ax.plot(subset["rt_minute"], subset["ids"], color="#2066b2", linewidth=1.2)
        ax.fill_between(subset["rt_minute"], subset["ids"], color="#76a9e0", alpha=0.35)
        ax.set_title(run_id)
        ax.set_xlabel("Acquisition time (min)")
        ax.set_ylabel("Canonical PSMs / min")
        ax.grid(axis="y", color="#d5d5d5", linewidth=0.6)
    for ax in axes.flat[len(run_ids) :]:
        ax.set_visible(False)
    fig.suptitle("PSM identifications over acquisition time", fontweight="bold")
    fig.tight_layout()
    return _save_figure(fig, output_root, file_fmts, png_res)


def _plot_composition(
    charge: pd.DataFrame,
    missed: pd.DataFrame,
    output_root,
    file_fmts,
    png_res,
):
    if charge.empty and missed.empty:
        return []
    import matplotlib.pyplot as plt
    import numpy as np

    frames = (
        (charge, "charge", "Charge state"),
        (missed, "missed_cleavages", "Missed cleavages"),
    )
    run_ids = list(
        dict.fromkeys(
            [
                *charge.get("run_id", pd.Series(dtype=str)).tolist(),
                *missed.get("run_id", pd.Series(dtype=str)).tolist(),
            ]
        )
    )
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.0, max(3.5, 1.8 + (0.45 * len(run_ids)))),
        squeeze=False,
    )
    image = None
    for ax, (frame, value_name, title) in zip(axes.flat, frames):
        if frame.empty:
            ax.set_visible(False)
            continue
        pivot = frame.pivot(
            index="run_id", columns=value_name, values="fraction"
        ).reindex(run_ids)
        pivot = pivot.reindex(sorted(pivot.columns), axis=1).fillna(0)
        image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=1)
        ax.set_title(title)
        ax.set_xticks(
            range(len(pivot.columns)), [str(value) for value in pivot.columns]
        )
        ax.set_yticks(range(len(pivot.index)), pivot.index)
        ax.set_xlabel(value_name.replace("_", " ").title())
        ax.set_ylabel("Run")
        for row_index, column_index in np.ndindex(pivot.shape):
            value = pivot.iat[row_index, column_index]
            if value >= 0.02:
                ax.text(
                    column_index,
                    row_index,
                    f"{value:.0%}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if value >= 0.55 else "#222222",
                )
    if image is not None:
        colorbar_ax = fig.add_axes((0.93, 0.22, 0.014, 0.52))
        fig.colorbar(
            image,
            cax=colorbar_ax,
            label="Fraction of canonical PSMs",
        )
    fig.suptitle("PSM composition", fontweight="bold")
    fig.subplots_adjust(left=0.16, right=0.91, bottom=0.16, top=0.82, wspace=0.45)
    return _save_figure(fig, output_root, file_fmts, png_res)


def run_psm_qc(
    *,
    config: Mapping,
    data_dir,
    output_stem,
    file_fmts=(".png",),
    png_res: int = 300,
) -> PsmQcArtifacts:
    """Summarize each unique configured search from one local QUAL scan."""
    _require_polars()
    data_dir = Path(data_dir).expanduser().resolve()
    if not data_dir.is_dir():
        raise ValueError(f"PSM data directory does not exist: {data_dir}")
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    grouped, config_errors = _config_triplets(config)
    qual_index = _qual_file_index(data_dir)
    summaries = list(config_errors)
    ids_frames = []
    charge_frames = []
    missed_frames = []

    for triplet, sections in grouped.items():
        run_fields = _run_fields(triplet, sections)
        path, status, message = _resolve_qual_file(qual_index, triplet)
        if path is None:
            summaries.append(
                {
                    **run_fields,
                    "source_file": "",
                    "status": status,
                    "message": message,
                }
            )
            warn(f"PSM QC {run_fields['run_id']}: {message}")
            continue
        try:
            summary, ids, charge, missed = _summarize_qual(path, triplet, sections)
        except Exception as exc:
            summaries.append(
                {
                    **run_fields,
                    "source_file": str(path),
                    "status": "error",
                    "message": str(exc),
                }
            )
            warn(f"PSM QC {run_fields['run_id']} failed: {exc}")
            continue
        summaries.append(summary)
        ids_frames.append(ids)
        charge_frames.append(charge)
        missed_frames.append(missed)

    summary_frame = _with_summary_columns(pd.DataFrame(summaries))
    ids_frame = _concat_frames(
        ids_frames,
        (*_run_fields(("", "", ""), []).keys(), "rt_minute", "ids"),
    )
    charge_frame = _concat_frames(charge_frames, _long_columns("charge"))
    missed_frame = _concat_frames(missed_frames, _long_columns("missed_cleavages"))

    summary_tsv = Path(f"{output_stem}_summary.tsv")
    ids_tsv = Path(f"{output_stem}_ids_over_time.tsv")
    charge_tsv = Path(f"{output_stem}_charge_distribution.tsv")
    missed_tsv = Path(f"{output_stem}_missed_cleavage_distribution.tsv")
    summary_frame.to_csv(summary_tsv, sep="\t", index=False)
    ids_frame.to_csv(ids_tsv, sep="\t", index=False)
    charge_frame.to_csv(charge_tsv, sep="\t", index=False)
    missed_frame.to_csv(missed_tsv, sep="\t", index=False)

    plot_paths = [
        *_plot_ids_over_time(
            ids_frame, Path(f"{output_stem}_ids_over_time"), file_fmts, png_res
        ),
        *_plot_composition(
            charge_frame,
            missed_frame,
            Path(f"{output_stem}_composition"),
            file_fmts,
            png_res,
        ),
    ]
    return PsmQcArtifacts(
        summary_tsv=summary_tsv,
        ids_over_time_tsv=ids_tsv,
        charge_distribution_tsv=charge_tsv,
        missed_cleavage_distribution_tsv=missed_tsv,
        plot_paths=tuple(plot_paths),
    )
