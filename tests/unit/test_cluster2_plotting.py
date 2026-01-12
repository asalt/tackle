import pandas as pd
import pytest

from tackle.cluster2.plotting import Cluster2FigsizeEnv, compute_cluster2_figsize


def _env() -> Cluster2FigsizeEnv:
    return Cluster2FigsizeEnv(
        width_scale=1.0,
        min_figwidth=5.4,
        max_figwidth=48.0,
        width_margin_overhead=1.2,
        legend_col_width=0.1,
        legend_row_height=0.35,
    )


def test_compute_cluster2_figsize_respects_figsize_when_no_gene_symbols():
    res = compute_cluster2_figsize(
        n_rows=25,
        n_cols=10,
        figsize=(11, 9),
        optimal_figsize=True,
        has_title=False,
        col_cluster=False,
        row_annot_df=None,
        col_meta=None,
        add_description=False,
        row_annot_side="right",
        row_names_side="right",
        show_gene_symbols=False,
        gene_symbol_fontsize=12,
        has_cut_by=False,
        env=_env(),
    )
    assert res.figwidth == 11.0
    assert res.figheight == 9.0


def test_compute_cluster2_figsize_overrides_height_for_gene_symbols_when_not_optimal():
    n_rows = 100
    font = 12
    res = compute_cluster2_figsize(
        n_rows=n_rows,
        n_cols=10,
        figsize=(11, 5),
        optimal_figsize=False,
        has_title=False,
        col_cluster=False,
        row_annot_df=None,
        col_meta=None,
        add_description=False,
        row_annot_side="right",
        row_names_side="right",
        show_gene_symbols=True,
        gene_symbol_fontsize=font,
        has_cut_by=False,
        env=_env(),
    )
    expected = max(((font + 2) / 72) * n_rows, 12)
    assert res.figheight == pytest.approx(expected)


def test_compute_cluster2_figsize_fills_missing_width_and_clamps():
    env = _env()
    res = compute_cluster2_figsize(
        n_rows=75,
        n_cols=20,
        figsize=(None, 10),
        optimal_figsize=False,
        has_title=True,
        col_cluster=True,
        row_annot_df=None,
        col_meta=None,
        add_description=True,
        row_annot_side="right",
        row_names_side="right",
        show_gene_symbols=False,
        gene_symbol_fontsize=12,
        has_cut_by=False,
        env=env,
    )
    assert res.figheight == 10.0
    assert env.min_figwidth <= res.figwidth <= env.max_figwidth


def test_compute_cluster2_figsize_counts_legend_groups_and_cut_by():
    env = _env()
    col_meta = pd.DataFrame(
        {
            "name": ["s1", "s2"],
            "treatment": ["A", "B"],
            "batch": ["1", "1"],
        }
    )
    row_annot_df = pd.DataFrame({"annot": ["x", "y", "z"]})
    res = compute_cluster2_figsize(
        n_rows=50,
        n_cols=10,
        figsize=None,
        optimal_figsize=True,
        has_title=False,
        col_cluster=False,
        row_annot_df=row_annot_df,
        col_meta=col_meta,
        add_description=False,
        row_annot_side="right",
        row_names_side="right",
        show_gene_symbols=False,
        gene_symbol_fontsize=12,
        has_cut_by=True,
        env=env,
    )
    assert int(res.debug["legend_groups"]) == 4
    assert res.debug["cut_by_extra"] == 0.4

