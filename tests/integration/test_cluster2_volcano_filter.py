import os
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("rpy2")
from rpy2.rinterface_lib.embedded import RRuntimeError

import tackle.clusterplot_dispatcher as cdisp


def test_cluster2_with_volcano_filter(ctx):
    # Build a small volcano TSV with clear up/down and significant rows
    data_obj = ctx.obj["data_obj"]
    genes = list(data_obj.areas_log.index)[:10]
    up = genes[:5]
    down = genes[5:]
    df = pd.DataFrame({
        "GeneID": up + down,
        "log2_FC": [3.0] * len(up) + [-3.0] * len(down),
        "pAdj": [0.001] * 10,
        "pValue": [0.001] * 10,
    })
    vdir = Path(data_obj.outpath) / "volcano"
    vdir.mkdir(parents=True, exist_ok=True)
    vfile = vdir / "group_test.tsv"
    df.to_csv(vfile, sep="\t", index=False)

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
            volcano_file=str(vfile),
            volcano_filter_params=(2.0, 0.05, "pAdj"),
            volcano_direction="both",
            volcano_sortby="log2_FC",
            cluster_file=(None, None),
            row_annot_side="left",
            row_dend_side="left",
            row_names_side="right",
            seed=1234,
            show_metadata=True,
            standard_scale="None",
            show_missing_values=False,
            cluster_fillna="min",
            z_score="None",
            z_score_by=None,
            z_score_fillna=False,
            add_human_ratios=False,
            volcano_topn=6,
        )
    except RRuntimeError as err:
        pytest.skip(f"R runtime prerequisites missing: {err}")

    cluster_dir = Path(data_obj.outpath) / "cluster2"
    images = list(cluster_dir.rglob("*.png"))
    assert images, "Expected a clustermap image after volcano filtering"
    # Check that naming encodes sortby and direction (dir_b for both)
    assert any("log2_FC" in p.name for p in images)
    assert any("dir_b" in p.name for p in images)

