from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd


def _safe_float(value: Any, default: float) -> float:
    try:
        return float(value)
    except Exception:
        return float(default)


@dataclass(frozen=True)
class Cluster2FigsizeEnv:
    width_scale: float = 1.0
    min_figwidth: float = 5.4
    max_figwidth: float = 48.0
    width_margin_overhead: float = 1.2
    legend_col_width: float = 0.1
    legend_row_height: float = 0.35

    ENV_WIDTH_SCALE: str = "TACKLE_WIDTH_SCALE"
    ENV_MIN_FIGWIDTH: str = "TACKLE_MIN_FIGWIDTH"
    ENV_MAX_FIGWIDTH: str = "TACKLE_MAX_FIGWIDTH"
    ENV_WIDTH_MARGIN_OVERHEAD: str = "TACKLE_WIDTH_MARGIN_OVERHEAD"
    ENV_LEGEND_COL_WIDTH: str = "TACKLE_LEGEND_COL_WIDTH"
    ENV_LEGEND_ROW_HEIGHT: str = "TACKLE_LEGEND_ROW_HEIGHT"

    @classmethod
    def from_env(cls) -> "Cluster2FigsizeEnv":
        return cls(
            width_scale=_safe_float(os.getenv(cls.ENV_WIDTH_SCALE), cls.width_scale),
            min_figwidth=_safe_float(os.getenv(cls.ENV_MIN_FIGWIDTH), cls.min_figwidth),
            max_figwidth=_safe_float(os.getenv(cls.ENV_MAX_FIGWIDTH), cls.max_figwidth),
            width_margin_overhead=_safe_float(
                os.getenv(cls.ENV_WIDTH_MARGIN_OVERHEAD), cls.width_margin_overhead
            ),
            legend_col_width=_safe_float(
                os.getenv(cls.ENV_LEGEND_COL_WIDTH), cls.legend_col_width
            ),
            legend_row_height=_safe_float(
                os.getenv(cls.ENV_LEGEND_ROW_HEIGHT), cls.legend_row_height
            ),
        )


@dataclass(frozen=True)
class Cluster2FigsizeResult:
    figwidth: float
    figheight: float
    debug: dict[str, float]


def _base_dims(n_rows: int, n_cols: int, *, add_title: bool) -> tuple[float, float]:
    row_slope = 0.24 if n_rows <= 250 else 0.18
    col_slope = 0.24 if n_cols <= 35 else 0.18
    base_h = 4.2 + (n_rows * row_slope)
    if add_title:
        base_h += 0.36
    base_w = 7.6 + (n_cols * col_slope)
    return float(base_h), float(base_w)


def _ensure_min_w(w: float, min_w: float) -> float:
    return max(float(w), float(min_w))


def _clamp_w(w: float, min_w: float, max_w: float) -> float:
    return max(min(float(w), float(max_w)), float(min_w))


def _apply_extras(
    h: float,
    w: float,
    *,
    col_cluster: bool,
    row_annot_cols: Optional[int],
    add_description: bool,
    margins_overhead: float,
    row_annot_side: str,
    row_names_side: str,
    show_gene_symbols: bool,
) -> tuple[float, float, dict[str, float]]:
    info: dict[str, float] = {
        "col_cluster": 0.0,
        "row_annot_w": 0.0,
        "row_annot_h": 0.0,
        "desc": 0.0,
        "margins": 0.0,
        "left_row_annot": 0.0,
        "left_row_names": 0.0,
    }
    if col_cluster:
        h += 3.6
        info["col_cluster"] = 3.6
    if row_annot_cols is not None and row_annot_cols > 0:
        _w = 0.4 + (0.26 * row_annot_cols)
        _h = 0.2 + (0.40 * row_annot_cols)
        w += _w
        h += _h
        info["row_annot_w"] = float(_w)
        info["row_annot_h"] = float(_h)
    if add_description:
        w += 1.8
        info["desc"] = 1.8

    w += margins_overhead
    info["margins"] = float(margins_overhead)

    if row_annot_side == "left" and row_annot_cols is not None and row_annot_cols > 0:
        w += 1.0
        info["left_row_annot"] = 1.0
    if row_names_side == "left" and show_gene_symbols:
        w += 0.8
        info["left_row_names"] = 0.8

    return float(h), float(w), info


def _legend_layout(
    *,
    col_meta: Optional[pd.DataFrame],
    row_annot_df: Optional[pd.DataFrame],
    legend_col_width: float,
    legend_row_height: float,
) -> tuple[float, float, int, int]:
    legend_groups = 1
    dense_units = 0

    if col_meta is not None:
        meta_cols = [c for c in col_meta.columns if c != "name"]
        legend_groups += len(meta_cols)
        for c in meta_cols:
            try:
                dense_units += max(int(col_meta[c].nunique()) - 6, 0)
            except Exception:
                pass

    if row_annot_df is not None:
        legend_groups += int(len(row_annot_df.columns))
        for c in list(row_annot_df.columns):
            try:
                dense_units += max(int(pd.Series(row_annot_df[c]).nunique()) - 6, 0)
            except Exception:
                pass

    legend_width_extra = (legend_groups * legend_col_width) + (0.08 * dense_units)
    legend_height_extra = (legend_groups * legend_row_height) + (0.04 * dense_units)
    return (
        float(legend_width_extra),
        float(legend_height_extra),
        int(legend_groups),
        int(dense_units),
    )


def compute_cluster2_figsize(
    *,
    n_rows: int,
    n_cols: int,
    figsize: Optional[tuple[Any, Any]] = None,
    optimal_figsize: bool,
    has_title: bool,
    col_cluster: bool,
    row_annot_df: Optional[pd.DataFrame],
    col_meta: Optional[pd.DataFrame],
    add_description: bool,
    row_annot_side: str,
    row_names_side: str,
    show_gene_symbols: bool,
    gene_symbol_fontsize: int,
    has_cut_by: bool,
    env: Optional[Cluster2FigsizeEnv] = None,
) -> Cluster2FigsizeResult:
    if env is None:
        env = Cluster2FigsizeEnv.from_env()

    if figsize is not None and len(figsize) != 2:
        raise ValueError(f"figsize must be a 2-tuple, got {figsize!r}")

    annot_cols: Optional[int] = None
    if row_annot_df is not None:
        try:
            annot_cols = int(len(row_annot_df.columns))
        except Exception:
            annot_cols = None

    legend_width_extra, legend_height_extra, legend_groups, dense_units = _legend_layout(
        col_meta=col_meta,
        row_annot_df=row_annot_df,
        legend_col_width=env.legend_col_width,
        legend_row_height=env.legend_row_height,
    )

    cut_by_extra = 0.4 if has_cut_by else 0.0

    base_h = float("nan")
    base_w = float("nan")
    extra_info: dict[str, float] = {
        "margins": 0.0,
        "row_annot_w": 0.0,
        "row_annot_h": 0.0,
        "desc": 0.0,
        "left_row_annot": 0.0,
        "left_row_names": 0.0,
    }
    w_pre = float("nan")

    # Case A: width missing
    if figsize is not None and (figsize[0] is None) and (figsize[1] is not None):
        if optimal_figsize:
            base_h, base_w = _base_dims(n_rows, n_cols, add_title=has_title)
        else:
            base_h = 10.56
            base_w = _ensure_min_w(min(n_cols / 2.0, 16.0), env.min_figwidth)
        h, w, extra_info = _apply_extras(
            base_h,
            base_w,
            col_cluster=col_cluster,
            row_annot_cols=annot_cols,
            add_description=add_description,
            margins_overhead=env.width_margin_overhead,
            row_annot_side=row_annot_side,
            row_names_side=row_names_side,
            show_gene_symbols=show_gene_symbols,
        )
        w_pre = (w + legend_width_extra + cut_by_extra) * env.width_scale
        figwidth = _clamp_w(w_pre, env.min_figwidth, env.max_figwidth)
        figheight = float(figsize[1])

    # Case A2: height missing
    elif figsize is not None and (figsize[1] is None) and (figsize[0] is not None):
        if optimal_figsize:
            base_h, base_w = _base_dims(n_rows, n_cols, add_title=has_title)
        else:
            base_h = 10.56
            base_w = _ensure_min_w(min(n_cols / 2.0, 16.0), env.min_figwidth)
        h, w, extra_info = _apply_extras(
            base_h,
            base_w,
            col_cluster=col_cluster,
            row_annot_cols=annot_cols,
            add_description=add_description,
            margins_overhead=env.width_margin_overhead,
            row_annot_side=row_annot_side,
            row_names_side=row_names_side,
            show_gene_symbols=show_gene_symbols,
        )
        h = h + legend_height_extra
        w_pre = (w + legend_width_extra + cut_by_extra) * env.width_scale
        figwidth = _clamp_w(_safe_float(figsize[0], env.min_figwidth), env.min_figwidth, env.max_figwidth)
        figheight = float(h)

    # Case B: neither provided (or both None)
    elif (figsize is None) or (figsize[0] is None and figsize[1] is None):
        if optimal_figsize:
            base_h, base_w = _base_dims(n_rows, n_cols, add_title=has_title)
        else:
            base_h = 12.0
            base_w = _ensure_min_w(min(n_cols / 2.0, 16.0), env.min_figwidth)
        h, w, extra_info = _apply_extras(
            base_h,
            base_w,
            col_cluster=col_cluster,
            row_annot_cols=annot_cols,
            add_description=add_description,
            margins_overhead=env.width_margin_overhead,
            row_annot_side=row_annot_side,
            row_names_side=row_names_side,
            show_gene_symbols=show_gene_symbols,
        )
        h = h + legend_height_extra
        w_pre = (w + legend_width_extra + cut_by_extra) * env.width_scale
        figwidth = _clamp_w(w_pre, env.min_figwidth, env.max_figwidth)
        figheight = float(h)

    # Case C: both provided
    else:
        figwidth = float(figsize[0])
        figheight = float(figsize[1])

    if show_gene_symbols and not optimal_figsize:
        figheight = max(((gene_symbol_fontsize + 2) / 72) * n_rows, 12)
        if figheight > 218:
            figheight = 218

    debug = {
        "base_w": float(base_w),
        "base_h": float(base_h),
        "legend_width_extra": float(legend_width_extra),
        "legend_height_extra": float(legend_height_extra),
        "legend_groups": float(legend_groups),
        "dense_units": float(dense_units),
        "cut_by_extra": float(cut_by_extra),
        "width_scale": float(env.width_scale),
        "pre_clamp_width": float(w_pre),
        "final_w": float(figwidth),
        "final_h": float(figheight),
        "margins_overhead": float(env.width_margin_overhead),
        "row_annot_w": float(extra_info.get("row_annot_w", 0.0)),
        "row_annot_h": float(extra_info.get("row_annot_h", 0.0)),
        "desc": float(extra_info.get("desc", 0.0)),
        "margins": float(extra_info.get("margins", 0.0)),
        "left_row_annot": float(extra_info.get("left_row_annot", 0.0)),
        "left_row_names": float(extra_info.get("left_row_names", 0.0)),
    }

    return Cluster2FigsizeResult(
        figwidth=float(figwidth),
        figheight=float(figheight),
        debug=debug,
    )
