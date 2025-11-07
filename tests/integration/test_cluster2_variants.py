import os
from pathlib import Path

import pytest

pytest.importorskip("rpy2")
from rpy2.rinterface_lib.embedded import RRuntimeError

import tackle.clusterplot_dispatcher as cdisp


def test_cluster2_with_scaling_and_zscore(ctx):
    try:
        cdisp.run(
            ctx=ctx,
            add_description=False,
            annotate=None,
            annotate_genes=False,
            cmap=None,
            cut_by=None,
            color_low="blue",
            color_mid="white",
            color_high="red",
            col_cluster=True,
            row_cluster=True,
            cluster_row_slices=True,
            cluster_col_slices=True,
            figwidth=6.0,
            figheight=6.0,
            figsize=None,
            force_plot_genes=False,
            genefile=None,
            genefile_sheet=0,
            gene_symbols=True,
            genesymbols=True,
            gene_symbol_fontsize=8.0,
            gene_annot=None,
            gsea_input=None,
            highlight_geneids=(),
            highlight_geneids_table=None,
            linear=False,
            legend_include=(),
            legend_exclude=(),
            optimal_figsize=False,
            sample_reference=None,
            sample_include=None,
            sample_exclude=None,
            linkage="ward",
            max_autoclusters=10,
            nclusters=3,
            cluster_func="Kmeans",
            main_title=None,
            order_by_abundance=False,
            volcano_file=None,
            volcano_filter_params=(2.0, 0.05, "pAdj"),
            volcano_direction="both",
            volcano_sortby="pValue",
            cluster_file=(None, None),
            row_annot_side="left",
            row_dend_side="left",
            row_names_side="right",
            seed=42,
            show_metadata=True,
            standard_scale="1",  # row scale
            show_missing_values=False,
            cluster_fillna="min",
            z_score="0",  # column z-score
            z_score_by=None,
            z_score_fillna=False,
            add_human_ratios=False,
            volcano_topn=50,
        )
    except RRuntimeError as err:
        pytest.skip(f"R runtime prerequisites missing: {err}")

    cluster_dir = Path(ctx.obj["data_obj"].outpath) / "cluster2"
    images = list(cluster_dir.rglob("*.png"))
    assert images, "Expected a clustermap image with scaling/z-score options"


def test_cluster2_sample_include_reduces_columns(ctx):
    try:
        cdisp.run(
            ctx=ctx,
            add_description=False,
            annotate=None,
            annotate_genes=False,
            cmap=None,
            cut_by=None,
            color_low="blue",
            color_mid="white",
            color_high="red",
            col_cluster=True,
            row_cluster=True,
            cluster_row_slices=True,
            cluster_col_slices=True,
            figwidth=6.0,
            figheight=6.0,
            figsize=None,
            force_plot_genes=False,
            genefile=None,
            genefile_sheet=0,
            gene_symbols=True,
            genesymbols=True,
            gene_symbol_fontsize=8.0,
            gene_annot=None,
            gsea_input=None,
            highlight_geneids=(),
            highlight_geneids_table=None,
            linear=False,
            legend_include=(),
            legend_exclude=(),
            optimal_figsize=False,
            sample_reference="condition",
            sample_include=["A"],
            sample_exclude=None,
            linkage="ward",
            max_autoclusters=10,
            nclusters=3,
            cluster_func="Kmeans",
            main_title=None,
            order_by_abundance=False,
            volcano_file=None,
            volcano_filter_params=(2.0, 0.05, "pAdj"),
            volcano_direction="both",
            volcano_sortby="pValue",
            cluster_file=(None, None),
            row_annot_side="left",
            row_dend_side="left",
            row_names_side="right",
            seed=7,
            show_metadata=True,
            standard_scale="None",
            show_missing_values=False,
            cluster_fillna="min",
            z_score="None",
            z_score_by=None,
            z_score_fillna=False,
            add_human_ratios=False,
            volcano_topn=50,
        )
    except RRuntimeError as err:
        pytest.skip(f"R runtime prerequisites missing: {err}")

    cluster_dir = Path(ctx.obj["data_obj"].outpath) / "cluster2"
    imgs = list(cluster_dir.rglob("*.png"))
    assert imgs, "Expected at least one clustermap image"
    # Assert filename hints reduced columns (look for 'x4' ending since 8 samples -> 4 after include)
    assert any(img.name.endswith("x4.png") for img in imgs), "Expected column count suffix x4 in output name"

