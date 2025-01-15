# cluster2.py

# from . import rutils
import os
import re
from copy import deepcopy
from functools import partial
import numpy as np
import pandas as pd

import rpy2
from rpy2.robjects import r
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

grdevices = importr("grDevices")

from . import utils
from . import grdevice_helper
from . import containers
from .utils import (
    fix_name,
    _get_logger,
    validate_in_config,
    parse_gid_file,
    get_file_name,
    get_default_color_mapping,
    get_outname,
)

from .containers import (
    get_annotation_mapper,
    get_gene_mapper,
    get_hgene_mapper,
)

logger = _get_logger(__name__)


def run(
    ctx,
    add_description,
    annotate,
    annotate_genes,
    cmap,
    cut_by,
    color_low,
    color_mid,
    color_high,
    col_cluster,
    row_cluster,
    cluster_row_slices,
    cluster_col_slices,
    figsize,
    force_plot_genes,
    genefile,
    genefile_sheet,
    gene_symbols,
    genesymbols,
    gene_symbol_fontsize,
    gene_annot,
    gsea_input,
    highlight_geneids,
    highlight_geneids_table,
    linear,
    legend_include,
    legend_exclude,
    optimal_figsize,
    sample_reference,
    sample_include,
    linkage,
    max_autoclusters,
    nclusters,
    cluster_func,
    main_title,
    order_by_abundance,
    volcano_file,
    volcano_filter_params,
    volcano_topn,
    volcano_direction,
    volcano_sortby,
    cluster_file,
    row_annot_side,
    row_dend_side,
    row_names_side,
    seed,
    show_metadata,
    standard_scale,
    show_missing_values,
    z_score,
    z_score_by,
    add_human_ratios,
):
    outname_kws = dict()

    outname_kws["rds"] = "l" if row_dend_side == "left" else "r"
    outname_kws["cc"] = "T" if col_cluster else "F"
    outname_kws["rc"] = "T" if row_cluster else "F"

    if genesymbols is True and gene_symbols is False:
        gene_symbols = True
    if z_score == "None":
        z_score = None
    if cluster_func == "none":
        cluster_func = None

    # =================================================================
    if gsea_input is not None:
        raise NotImplementedError()

    pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(
        os.path.split(os.path.abspath(__file__))[0], "R", "clusterplot.R"
    )

    r_source(r_file)
    cluster2 = robjects.r["cluster2"]

    # =================================================================

    data_obj = ctx.obj["data_obj"]
    col_meta = data_obj.col_metadata.copy().astype(str).fillna("")
    col_meta = data_obj.col_metadata

    if order_by_abundance and row_cluster:
        raise NotImplementedError("Not Implemented !")

    if not figsize:  # returns empty tuple if not specified
        figsize = None

    missing_values = "masked" if show_missing_values else "unmasked"

    X = data_obj.areas_log
    X.index = X.index.astype(str)
    if show_missing_values:
        _mask = data_obj.mask
        X[_mask] = np.nan
    # filter if sample include is mentioned
    if sample_reference is not None:
        # we take a subset of the data
        col_meta = col_meta.loc[col_meta[sample_reference].isin(sample_include)]
        X = X[col_meta.index]

    genes = None
    column_title = None
    if genefile and not isinstance(genefile, dict):
        genes = parse_gid_file(genefile, sheet=genefile_sheet)  # default 0
        outname_kws["genefile"] = fix_name(os.path.splitext(genefile)[0])
        column_title = outname_kws["genefile"]
    elif genefile and isinstance(genefile, dict):  # this is if already processed
        _key = list(genefile.keys())[0]
        genes = genefile[_key]
        outname_kws["genefile"] = _key
    if genefile:
        _tokeep = [x for x in genes if x in X.index]  # preserves order
        # X = X.loc[set(X.index) & set(genes)]
        X = X.loc[_tokeep]
    # print(len(X))

    ## filter by volcano output
    # ================================ volcano file parsing =================================
    _tmp = list()
    _fc, _pval, _ptype = volcano_filter_params
    # for f in volcano_file:
    if volcano_file is not None:
        logger.info(f"Loading volcano file: {volcano_file}")
        volcanofile_basename = os.path.split(volcano_file)[-1]
        if "Batch" in volcanofile_basename:
            volcanofile_basename = volcanofile_basename[
                volcanofile_basename.find("Batch") + 12 :
            ]
            # name_group = volcanofile_basename[ volcanofile_basename.find("_group"): ]
        name_group = re.search(r"(?<=group)[_]?(.*)(?=\.tsv)", volcanofile_basename)
        if name_group is None:
            name_group = volcanofile_basename
        else:
            name_group = name_group.group(1)
        logger.info(f"volcanofile name group is {name_group}")
        outname_kws["volcano_file"] = name_group
        outname_kws["direction"] = volcano_direction
        column_title = name_group
        _df = pd.read_table(volcano_file)
        if "pValue" not in _df and "p-value" in _df:
            _df = _df.rename(columns={"p-value": "pValue"})
        if "log2_FC" not in _df and "Value" in _df:  #
            _df = _df.rename(columns={"Value": "log2_FC"})

        if not np.isfinite(volcano_topn):
            _df = _df[(abs(_df["log2_FC"]) > np.log2(_fc)) & (_df[_ptype] < _pval)]
        else:  # np.isfinite(volcano_topn):
            #
            _n = (
                int(volcano_topn / 2)
                if volcano_direction == "both"
                else int(volcano_topn)
            )
            if volcano_sortby == "log2_FC":
                _up = _df.sort_values(by="log2_FC", ascending=False).head(_n)
                _dn = _df.sort_values(by="log2_FC", ascending=False).tail(_n)
            elif volcano_sortby == "pValue":
                _up = (
                    _df.query("log2_FC>0")
                    .sort_values(by=["pValue", "log2_FC"], ascending=[True, False])
                    .head(_n)
                )
                _dn = (
                    _df.query("log2_FC<0")
                    .sort_values(by=["pValue", "log2_FC"], ascending=[True, True])
                    .head(_n)
                )
                _dn = _dn[::-1]

            if volcano_direction == "both":
                _df = pd.concat([_up, _dn])
            elif volcano_direction == "up":
                _df = _up
            elif volcano_direction == "down":
                _df = _dn
        _tmp.append(_df)
    if _tmp:
        _dfs = pd.concat(_tmp)
        _genes = _dfs.GeneID.astype(str).unique()
        _tokeep = [x for x in _genes if x in X.index]  # preserve order
        # X = X.loc[set(X.index) & set(genes)]
        X = X.loc[_tokeep]

    # ================================ end of volcano file parsing =================================
    if cluster_file[0] is not None:
        if cluster_file[0].endswith("xlsx"):
            _df = pd.read_excel(cluster_file[0])
        else:
            _df = pd.read_table(cluster_file[0])

        _df = _df.rename(
            columns={
                "Cluster": "cluster",
                "Kmeans_Cluster": "cluster",
                "kmeans_cluster": "cluster",
                "kmeans_Cluster": "cluster",
                "Kmeans_cluster": "cluster",
                "PAM_cluster": "cluster",
            }
        )

        _thecluster = cluster_file[1]
        if "cluster" not in _df:
            raise ValueError("improper cluster file, does not have column `cluster`")
        if _thecluster not in _df["cluster"]:
            raise ValueError(f"cluster {_thecluster} not in {cluster_file[0]}")
        _genes = _df[_df["cluster"] == _thecluster]["GeneID"].astype(str)
        _tokeep = [x for x in _genes if x in X.index]
        X = X.loc[_tokeep]
        outname_kws["subcluster"] = _thecluster
        main_title += f"\nCluster {_thecluster}"
    # ========================= end of special ifle parsing ==========================

    gids_to_annotate = None
    if gene_annot:
        gids_to_annotate = parse_gid_file(gene_annot)

    if linear:
        X = 10**X

    if standard_scale is not None and standard_scale != "None":
        if standard_scale == 1 or standard_scale == "1":
            X = ((X.T / X.max(1)).T).fillna(0)
        elif standard_scale == 0 or standard_scale == "0":
            X = (X / X.max(0)).fillna(0)

    symbols = [data_obj.gid_symbol.get(x, "?") for x in X.index]

    genemapper = get_gene_mapper()
    symbols = [
        data_obj.gid_symbol.get(x, genemapper.symbol.get(str(x), str(x)))
        for x in X.index
    ]
    X = X.reset_index().assign(GeneSymbol=symbols)

    # X = data_obj.areas_log_shifted.copy()
    # X[data_obj.areas == 0] = 0 # fill the zeros back
    # X[data_obj.mask] = np.NaN

    # data_obj.set_highlight_gids(highlight_geneids)

    row_annot_track = list()
    # TODO parse and look for geneids
    for file_ in highlight_geneids:

        if file_.endswith("xlsx") or file_.endswith("xls"):
            df_ = pd.read_excel(file_, header=None, dtype=str)
        else:
            df_ = pd.read_csv(file_, sep="\t", header=None, dtype=str)
        title_ = get_file_name(file_)
        if len(df_.columns) != 2:
            raise ValueError(
                """File {} should have two columns, first being GeneID,
            second being corresponding value. No headers""".format(
                    file_
                )
            )
        df_.columns = ["GeneID", title_]
        df_ = df_.set_index("GeneID")
        row_annot_track.append(df_)

    if highlight_geneids_table:
        df_ = parse_gid_file(highlight_geneids_table)
        # import ipdb; ipdb.set_trace()
        # row_annot_track.append(df_)
    # #xx = containers.add_annotations(X)
    # this needs rewriting
    if annotate_genes:
        aa = containers.get_annotation_mapper()
        row_annot_track.append(
            aa.df[
                [x for x in aa.df if x != "GeneSymbol" and x != "MitoCarta_Pathways"]
            ].set_index("GeneID")
        )
        outname_kws["geneanno"] = "T"

    row_annot_df = None
    if row_annot_track:
        row_annot_df = pd.concat(row_annot_track, axis=1)
        # make a dataframe that spans all genes about to be plotted
        ixs_ = X.GeneID.astype(str)
        missing_ = list(set(ixs_) - set(row_annot_df.index))
        intersect_ixs_ = list(set(ixs_) & set(row_annot_df.index))
        missing_df_ = pd.DataFrame(index=missing_, columns=row_annot_df.columns)
        row_annot_df = pd.concat(
            [row_annot_df.loc[intersect_ixs_], missing_df_], axis=0
        ).fillna("")

    # data_obj.standard_scale    = data_obj.clean_input(standard_scale)
    # data_obj.z_score           = data_obj.clean_input(z_score)

    # col_meta = data_obj.col_metadata.copy().pipe(clean_categorical)

    if add_human_ratios:
        col_meta["HS_ratio"] = data_obj.taxon_ratios["9606"]
    # _expids = ('recno', 'runno', 'searchno', 'label')
    _expids = [
        "recno",
        "runno",
        "searchno",
    ]
    if ("label" in col_meta.columns) and (col_meta["label"].nunique() == 1):
        _expids.append("label")

    # valid_entries = validate_in_config(*legend_include, *legend_exclude, valid_entries=col_meta.index)
    valid_entries = validate_in_config(
        *legend_include, *legend_exclude, valid_entries=col_meta.columns
    )
    legend_include = [
        x for x in legend_include if x in (set(legend_include) & set(valid_entries))
    ]
    legend_exclude = set(legend_exclude) & set(valid_entries)
    to_exclude = set(_expids) | set(legend_exclude)
    if legend_include:
        to_include = legend_include
    else:
        to_include = set(col_meta.columns) - to_exclude
    # col_meta = col_meta.loc[[x for x in col_meta.index if x not in _expids]]
    col_meta = col_meta[[x for x in col_meta.columns if x in to_include]]
    if col_meta.empty:
        col_meta = None
    else:
        col_meta.index.name = "name"
        col_meta = col_meta.reset_index()
        col_meta = utils.set_pandas_datatypes(col_meta)

    ## ============================================================
    # experimental, best to be put in a separate function
    if cut_by is not None:
        if ":" in cut_by:
            cut_by = cut_by.split(":")
        else:
            cut_by = [cut_by]
        cut_by = np.array(cut_by)
        for c in cut_by:
            if col_meta is not None and c not in col_meta.columns:
                print("{} not in column metadata".format(cut_by))
                cut_by = None
                break
    if cut_by is not None:
        if hasattr(cut_by, "__iter__"):
            outname_kws["cut_by"] = str.join("_", cut_by)
        else:
            outname_kws["cut_by"] = cut_by
    elif cut_by is None:
        cut_by = robjects.NULL
    ## ============================================================
    annot_mat = None
    if annotate:
        # annot_mat = (data_obj.data.loc[ idx[X.index.tolist(), annotate], : ]
        _cols = [x for x in data_obj.data.columns if x not in ("Metric")]
        annot_mat = (
            data_obj.data.loc[
                (data_obj.data.GeneID.isin(X.GeneID.tolist()))
                & (data_obj.data.Metric == annotate)
            ][_cols]
            # .reset_index(level=1, drop=True)
            .set_index("GeneID")
            .fillna(0)
            .astype(int)
        )
        annot_mat.columns.name = annotate
        outname_kws["annot"] = annotate

        annot_mat["GeneID"] = annot_mat.index
        annot_mat = annot_mat[["GeneID"] + [x for x in annot_mat if x != "GeneID"]]

    # ============================================================
    if z_score_by is not None:
        outname_kws["zby"] = z_score_by

    if linear:
        outname_kws["linear"] = "linear"
    if standard_scale is not None:
        if standard_scale == 1 or standard_scale == "1":
            outname_kws["scale"] = "row"
        if standard_scale == 0 or standard_scale == "0":
            outname_kws["scale"] = "col"

    #
    # data_obj.do_cluster()
    outname_func = partial(
        get_outname,
        # name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        batch=data_obj.batch_applied,
        batch_method="param" if not data_obj.batch_nonparametric else "nonparam",
        link=linkage,
        missing_values=missing_values,
        normtype=data_obj.normtype,
        outpath=data_obj.outpath,  # , "cluster2")
    )
    # =================================================================
    if gene_symbols and add_description:
        # this could definitely be improved
        description_frame = data_obj.data.query("Metric == 'GeneDescription'")
        _cols = [x for x in description_frame if x not in ("GeneID", "Metric")]
        _descriptions = description_frame[_cols].stack().unique()
        if all(x == "" for x in _descriptions):
            description_frame = data_obj.data.query(
                "Metric == 'Description'"
            )  # try another
            _descriptions = description_frame[_cols].stack().unique()
        if len(description_frame) == 0 or all(x == "" for x in _descriptions):
            pass  # fail
        else:
            description_frame = description_frame.bfill(axis=1).ffill(axis=1)
            cols_for_merge = ["GeneID", description_frame.columns[-1]]
            to_merge = description_frame[cols_for_merge].rename(
                columns={description_frame.columns[-1]: "Description"}
            )
            X = pd.merge(X, to_merge, how="left")
            # try to remove excess info from description that are not useful for display
            # this is uniprot based removal of identifiers
            X["Description"] = X.Description.str.replace(r"OS=.*", "", regex=True)
            X["GeneSymbol"] = X.GeneSymbol.astype(str) + " " + X.Description.astype(str)
            X = X.drop("Description", axis=1)
        outname_kws["descr"] = "T"

    metadata_colorsR = None
    metadata_color_lists = None

    # columns to assign colors:
    # TODO prevent crash if col_meta is empty
    # metacats = col_meta.dtypes[col_meta.dtypes.isin(["string", "object"])].index
    # metacats = col_meta.dtypes[
    #     (col_meta.dtypes == "string") | (col_meta.dtypes == "object")
    # ].index
    # gene_metacats = None
    metacats = col_meta.columns
    if row_annot_df is not None:
        metacats = set(metacats) | set(row_annot_df.columns)

    metadata_color_list = list()
    for metacat in metacats:
        if data_obj.metadata_colors is not None and metacat in data_obj.metadata_colors:
            # metadata_colorsR = robjects.ListVector([])
            # TODO handle when there is not a color for a category
            mapping = data_obj.metadata_colors[metacat]
            entry = robjects.vectors.ListVector(
                {metacat: robjects.vectors.ListVector(mapping)}
            )
            metadata_color_list.append(entry)
            # for value in :  # value is a key-value mapping
            # returns a numpy array with loss of names..
            # entry = robjects.vectors.ListVector({key: robjects.vectors.ListVector(x)})
            # if key != 'response':
            #     continue
            # it works!
            # actually needs to be an "atomic" vector and not a "list vector"
            # this is taken care of in R..
            # Unclear if this can be converted here. Calling R's `unlist` function on the resulting ListVector
        else:
            # if metacat not in col_meta:
            #     continue
            if metacat in col_meta:
                themapping = get_default_color_mapping(
                    col_meta[metacat]
                )  # set everything but float
            if row_annot_df is not None and metacat in row_annot_df:
                themapping = get_default_color_mapping(
                    row_annot_df[metacat]
                )  # set everything but float
            # import ipdb; ipdb.set_trace()
            if themapping is not None:
                entry = robjects.vectors.ListVector(
                    {metacat: robjects.vectors.ListVector(themapping)}
                )
                metadata_color_list.append(entry)
            else:
                pass  # for numeric

                # metadata_color_list.append(
                #     robjects.vectors.ListVector({metacat: robjects.NULL})
                # )

    if metadata_color_list:
        metadata_colorsR = metadata_color_list

    # print(metadata_color_list)
    # ========================================================

    if gids_to_annotate:
        # need to convert to a flat R vector
        gids_to_annotate = robjects.vectors.StrVector(
            [str(x) for x in gids_to_annotate]
        )

    # can we make sure this is always the case everywhere, it makes things much simpler and should
    # avoid bugs
    X["GeneID"] = X.GeneID.astype(str)

    if X.empty:
        print("No data!")
        return

    col_data = robjects.NULL
    if show_metadata and not col_meta is None:
        col_data = col_meta
    #     # cannot convert Categorical column of Integers to Category in py2r
    #     col_data = col_meta.pipe(clean_categorical)  # does this fix the problem?
    # col_data = utils.set_pandas_datatypes(col_data)

    # rpy2 does not map None to robjects.NULL
    if row_annot_df is None:
        row_annot_df = robjects.NULL

    pandas2ri.activate()

    def plot_and_save(
        X,
        out,
        # gr_kws=None,
        grdevice=None,  # depreciated
        annot_mat=None,
        main_title=None,
        column_title=column_title,  # defined in outer scope
        file_fmt=".pdf",
        figsize=None,
        **kws,
    ):
        """
        # gr_kws is redefined here (update fighiehgt /width), no need to take it as input
        """
        print("Saving", out, "...", end="", flush=True)
        # have to put this here to preserve the layout set by ComplexHeatmap::draw
        if annot_mat is None:
            annot_mat = robjects.NULL
        if main_title is None:
            main_title = robjects.NULL

        min_figwidth = 5
        if figsize is None:  # either None or length 2 tuple
            if optimal_figsize:
                figheight = 4 + (X.shape[0] * 0.22)
                figwidth = 8 + (X.shape[1] * 0.26)
                if col_cluster:
                    figheight += 3.2
            else:
                figheight = 12
                figwidth = max(min(len(X.columns) / 2, 16), min_figwidth)
            if row_annot_df is not None and row_annot_df is not robjects.NULL:
                figwidth += 0.15 * len(row_annot_df.columns)
                figheight += 0.4 * len(row_annot_df.columns)
        else:
            figwidth, figheight = figsize
            if (
                gene_symbols and not optimal_figsize
            ):  # make sure there is enough room for the symbols
                figheight = max(((gene_symbol_fontsize + 2) / 72) * len(X), 12)
                if figheight > 218:  # maximum figheight in inches
                    FONTSIZE = max(218 / figheight, 6)
                    figheight = 218
        print(figwidth, figheight)

        call_kws = dict(
            data=X,
            color_low=color_low,
            color_mid=color_mid,
            color_high=color_high,
            cut_by=cut_by,
            annot_mat=annot_mat,
            the_annotation=annotate or robjects.NULL,
            z_score=(
                z_score if z_score != "None" and z_score is not None else robjects.NULL
            ),
            z_score_by=z_score_by or robjects.NULL,
            # data_obj.normtype, # can add normtype info here to label on cbar, perhaps
            row_annot_df=row_annot_df,
            col_data=col_data,
            # cmap_name=cmap or np.nan,
            gids_to_annotate=gids_to_annotate or robjects.NULL,
            force_plot_genes=force_plot_genes,
            # genes=genes or robjects.NULL, # this has been filtered above
            show_gene_symbols=gene_symbols,
            standard_scale=standard_scale or robjects.NULL,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            # metadata=data_obj.config if show_metadata else None,
            nclusters=nclusters or robjects.NULL,
            cluster_func=cluster_func or robjects.NULL,
            max_autoclusters=max_autoclusters,
            show_missing_values=show_missing_values,
            main_title=main_title or column_title or robjects.NULL,
            # mask=data_obj.mask,
            # figsize=figsize,
            linear=linear,  # linear transformation not done in R yet
            normed=data_obj.normed,
            linkage=linkage,
            gene_symbol_fontsize=gene_symbol_fontsize,
            # legend_include=legend_include,
            # legend_exclude=legend_exclude,
            row_dend_side=row_dend_side,
            row_names_side=row_names_side,
            cluster_row_slices=cluster_row_slices,
            cluster_col_slices=cluster_col_slices,
            order_by_abundance=order_by_abundance,
            seed=seed or robjects.NULL,
            metadata_colors=metadata_colorsR or robjects.NULL,
            savedir=os.path.abspath(
                os.path.join(data_obj.outpath, "cluster2")
            ),  # this is for saving things within the r function
            fixed_size=True if optimal_figsize else False,
        )
        call_kws.update(kws)
        # logger.info(f"call_kws: {call_kws}")

        if len(out) > 299:  # quick fix
            out_path, out_name = os.path.split(out)
            out_name = (
                out_name.replace("treatment", "treat")
                .replace("clustermap", "cmap")
                .replace("normtype", "norm")
                .replace("volcano_file", "vfile")
                .replace("direction_", "")
                .replace("genotype", "geno")
            )
            out = os.path.join(out_path, out_name)

        logger.info(f"figheight: {figheight}, figwidth: {figwidth}")

        print(figwidth, figheight)
        grdevice = grdevice_helper.get_device(
            filetype=file_fmt, width=figwidth, height=figheight
        )
        grdevice(file=out)  # open file for saving
        ret = cluster2(**call_kws)  # draw
        # print(ret[0])
        grdevices.dev_off()  # close file
        print("done.", flush=True)

        sil_df = None
        try:
            sil_df = rpy2.robjects.pandas2ri.rpy2py_dataframe(ret[1])
        except Exception as e:
            pass

        if sil_df is not None:
            row_orders = ret[2]
            the_orders = [
                row_orders.rx2(str(n)) - 1 for n in row_orders.names
            ]  # subtract 1  for zero indexing
            the_orders = [x for y in the_orders for x in y]

            cluster_metrics = sil_df.iloc[the_orders]

            out = outname_func("clustermap_clusters", **outname_kws)
            print("saving", out)
            cluster_metrics.to_csv(out, index=False, sep="\t")

    # ==============================================================================================

    if cluster_func is not None:
        outname_kws[cluster_func] = nclusters
    # outname = outname_func("clustermap", **outname_kws)

    for file_fmt in ctx.obj["file_fmts"]:
        # grdevice = gr_devices[file_fmt] # this is now done within plot_and_save
        # gr_kw = gr_kws[file_fmt]
        annot_mat_to_pass = annot_mat
        # if len(X) > 300 and annotate:
        #     annot_mat_to_pass = None
        #     logger.info(f"number of genes is {len(X)} >> 300, skipping annotation")
        #     if "annotate" in outname_kws:
        #         outname_kws.pop("annotate")

        #################################################################
        ##                          make plot                          ##
        #################################################################

        _ncol = len([x for x in X.columns if x not in ("GeneID", "GeneSymbol")])
        _nrow = X.shape[0]
        nrow_ncol = "_{0}x{1}".format(_nrow, _ncol)
        # import ipdb; ipdb.set_trace()
        outname_base_name = "clustermap"
        extra_outname_kws = dict()
        if "volcano_file" in outname_kws:  # 'hack' to nest 1 folder deeper
            volcano_file = outname_kws.pop("volcano_file")
            outname_base_name = os.path.join(outname_base_name, volcano_file)
            extra_outname_kws["volcano_file"] = volcano_file

        out = outname_func(outname_base_name, **outname_kws) + nrow_ncol + file_fmt

        outname_kws.update(extra_outname_kws)

        plot_and_save(
            X,
            out,
            # grdevice,
            annot_mat=annot_mat_to_pass,
            main_title=main_title,
            file_fmt=file_fmt,
            figsize=figsize,
        )

        ##################################################################
        ##                     plot any annotations                     ##
        ##################################################################

        if data_obj.annotations is None:
            continue

        for annotation in data_obj.annotations:
            annotator = get_annotation_mapper()
            # annot_df = annotator.get_annot(annotation)
            annot_df = annotator.df

            subX = X[X.GeneID.isin(annot_df.GeneID)]
            if subX.empty:  # try mapping to homologene
                logger.info(f"Trying to map genes to hs through homologene")
                hgene_mapper = get_hgene_mapper()
                hg_gene_dict = hgene_mapper.map_to_human(X.GeneID)
                hs_genes = set(hg_gene_dict.values())
                _annot_genes = set(annot_df.GeneID) & set(hg_gene_dict.values())
                _tokeep = [k for k, v in hg_gene_dict.items() if v in _annot_genes]
                _tokeep = set(_tokeep)
                logger.info(
                    f"{annotation}: Successfully remapped {len(_tokeep)} genes to hs genes"
                )
                subX = X[X.GeneID.isin(_tokeep)]
            if subX.empty:  # if still empty, continue
                continue

            sub_annot_mat = None
            if annotate and annot_mat is not None:
                sub_annot_mat = annot_mat[annot_mat.GeneID.isin(subX.GeneID)]
            _show_gene_symbols = False
            if len(subX) < 141 and bool(gene_symbols) == False:
                _show_gene_symbols = True
                logger.info(
                    f"number of genes is {len(subX)} << 101, adding symbols. Specify --gene-symbols to override>."
                )
            if gene_symbols is True:  # force override
                _show_gene_symbols = True

            nrow_col = "_{0}x{1}".format(*subX.shape)
            out = (
                outname_func("clustermap", geneset=fix_name(annotation), **outname_kws)
                + nrow_ncol
                + file_fmt
            )
            out = fix_name(out)  # make it shorter
            plot_and_save(
                subX,
                out,
                grdevice,
                main_title=annotation,
                annot_mat=sub_annot_mat,
                show_gene_symbols=_show_gene_symbols,
                file_fmt=file_fmt,
            )

        # heatmap = ret.rx('heatmap')
        # print(heatmap)

    # ==============================================================================================

    # ==============================================================================================
    # this doesn't work yet, fix it later
    # discrete_clusters = ret.rx2('discrete_clusters')

    # if discrete_clusters is not robjects.NULL:
    #     the_clusters = discrete_clusters[0]
    #     X['Cluster'] = the_clusters

    # name = outname_func('{}_{}'.format(cluster_func, nclusters)) + '.tsv'
    # X[['GeneID', 'GeneSymbol', 'Cluster']].sort_values(by='Cluster').to_csv(name, sep='\t', index=False)
