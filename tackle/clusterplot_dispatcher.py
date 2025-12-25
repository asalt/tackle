# cluster2.py

# from . import rutils
import os
import re
from copy import deepcopy
from functools import partial
import numpy as np
import pandas as pd


from . import utils
from . import containers
from . import statfile_sorter
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


def _series_matches_includes(series: pd.Series, include_values) -> pd.Series:
    """Return boolean mask matching include_values against a metadata series.

    Supports numeric metadata columns even when include_values are provided as strings
    (e.g. "1" will match 1.0).
    """
    if include_values is None:
        raise ValueError("include_values cannot be None")
    include_list = list(include_values)
    if not include_list:
        return pd.Series(False, index=series.index)

    as_str = [str(v) for v in include_list]
    numeric_includes = []
    string_includes = []
    for v in as_str:
        try:
            numeric_includes.append(float(v))
        except ValueError:
            string_includes.append(v)

    mask = pd.Series(False, index=series.index)
    if numeric_includes:
        numeric_series = pd.to_numeric(series, errors="coerce")
        mask |= numeric_series.isin(numeric_includes)
    if string_includes:
        mask |= series.astype(str).isin(string_includes)
    return mask


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
    figwidth,
    figheight,
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
    sample_exclude,  # list
    linkage,
    max_autoclusters,
    nclusters,
    cluster_func,
    main_title,
    order_by_abundance,
    volcano_file,
    volcano_filter_params,
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
    cluster_fillna,
    z_score,
    z_score_by,
    z_score_fillna,
    add_human_ratios,
    volcano_topn=50,
):

    from . import grdevice_helper
    import rpy2
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter

    grdevices = importr("grDevices")

    outname_kws = dict()
    print(volcano_direction)

    # hacky sloppy name assignment
    outname_kws["rds" + "l" if row_dend_side == "left" else "rds" + "r"] = ""
    outname_kws["rc" + "T" if row_cluster else "rc" + "F"] = ""
    outname_kws["cc" + "T" if col_cluster else "cc" + "F"] = ""

    if genesymbols is True and gene_symbols is False:
        gene_symbols = True
    if z_score == "None":
        z_score = None
    if cluster_func == "none":
        cluster_func = None

    # =================================================================
    if gsea_input is not None:
        raise NotImplementedError()

    # pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(
        os.path.split(os.path.abspath(__file__))[0], "R", "clusterplot.R"
    )

    r_source(r_file)
    cluster2 = robjects.r["cluster2"]

    # =================================================================

    data_obj = ctx.obj["data_obj"]
    # col_meta = data_obj.col_metadata.copy().astype(str).fillna("")
    col_meta = data_obj.col_metadata.fillna("") # compromise

    if order_by_abundance and row_cluster:
        raise NotImplementedError("Not Implemented !")

    if not figsize:  # returns empty tuple if not specified
        figsize = (figwidth, figheight)
    missing_values = "masked" if show_missing_values else "unmasked"
    if cluster_fillna == "avg":
        missing_values += "_avg"

    X = data_obj.areas_log
    X.index = X.index.astype(str)
    if show_missing_values:
        _mask = data_obj.mask
        X[_mask] = np.nan

    # filter if sample include is mentioned
    if sample_reference is not None:
        if not sample_include:
            raise ValueError(
                "`--sample-reference` requires at least one `--sample-include` value"
            )
        if sample_reference not in col_meta.columns:
            raise ValueError(
                f"Unknown `--sample-reference` {sample_reference!r}; available: {list(col_meta.columns)}"
            )
        # we take a subset of the data
        mask = _series_matches_includes(col_meta[sample_reference], sample_include)
        col_meta = col_meta.loc[mask]
        X = X[col_meta.index]

    if sample_exclude is not None:
        to_exclude = set(sample_exclude) & set(col_meta.index)
        to_keep = [
            x for x in col_meta.index if x not in to_exclude
        ]  # use a list instead of set to preserve order
        col_meta = col_meta.loc[to_keep]
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
        X = X[~X.index.duplicated(keep="first")]  # just in case there are duplicate ids
        if force_plot_genes:  # add these in and fill as NA
            # Account for codepaths where GeneID column is not created yet
            existing_gene_ids = X.index
            if "GeneID" in X.columns:
                existing_gene_ids = X["GeneID"].astype(str)
            missing = set(genes) - set(existing_gene_ids)
            if missing:
                Xmissing = pd.DataFrame(index=list(missing), columns=X.columns)
                Xmissing.index.name = X.index.name
                X = pd.concat([X, Xmissing])

    # print(len(X))

    ## filter by volcano output
    # ================================ volcano file parsing =================================

    _tmp = list()
    _fc, _pval, _ptype = volcano_filter_params
    # for f in volcano_file:
    if volcano_file is not None:
        X = statfile_sorter.sort_files(
            [volcano_file],
            X,
            sort_by=volcano_sortby,
            direction=volcano_direction,
            topn=volcano_topn,
            fc=_fc,
            pval_cutoff=_pval,
            pval_type=_ptype
        )
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
        outname_kws["vfile"] = (
            name_group.replace(":", "_")
            .replace("+", "_")
            .replace("?", "qmk")
            .replace("|", "or")
            .replace(r"/", "_")
            .replace(r"\\", "_")
        )
        outname_kws[volcano_sortby] = ""
        outname_kws[f"dir_{volcano_direction[0]}"] = ""
        column_title = name_group  # _df = pd.read_table(volcano_file) # if "pValue" not in _df and "p-value" in _df:
        #     _df = _df.rename(columns={"p-value": "pValue"})
        # if "log2_FC" not in _df and "Value" in _df:  #
        #     _df = _df.rename(columns={"Value": "log2_FC"})

        # if not np.isfinite(volcano_topn):
        #     _df = _df[(abs(_df["log2_FC"]) > np.log2(_fc)) & (_df[_ptype] < _pval)]
        # else:  # np.isfinite(volcano_topn):
        #     #
        #     _n = (
        #         int(volcano_topn / 2)
        #         if volcano_direction == "both"
        #         else int(volcano_topn)
        #     )
        #     if volcano_sortby == "log2_FC":
        #         _up = _df.sort_values(by="log2_FC", ascending=False).head(_n)
        #         _dn = _df.sort_values(by="log2_FC", ascending=False).tail(_n)
        #     elif volcano_sortby == "pValue":
        #         _up = (
        #             _df.query("log2_FC>0")
        #             .sort_values(by=["pValue", "log2_FC"], ascending=[True, False])
        #             .head(_n)
        #         )
        #         _dn = (
        #             _df.query("log2_FC<0")
        #             .sort_values(by=["pValue", "log2_FC"], ascending=[True, True])
        #             .head(_n)
        #         )
        #         _dn = _dn[::-1]

        #     if volcano_direction == "both":
        #         _df = pd.concat([_up, _dn])
        #     elif volcano_direction == "up":
        #         _df = _up
        #     elif volcano_direction == "down":
        #         _df = _dn
        # _tmp.append(_df)
    # if _tmp:
    #     _dfs = pd.concat(_tmp)
    #     _genes = _dfs.GeneID.astype(str).unique()
    #     _tokeep = [x for x in _genes if x in X.index]  # preserve order
    #     # X = X.loc[set(X.index) & set(genes)]
    #     X = X.loc[_tokeep]

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
    # ========================= end of special file parsing ==========================

    gids_to_annotate = None
    if gene_annot:
        gids_to_annotate = parse_gid_file(gene_annot)  # info

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
        df_[title_] = df_[title_].fillna("")
        df_ = df_.set_index("GeneID")
        # checking for uniqueness of the index would help avoid crashes
        row_annot_track.append(df_)

    if highlight_geneids_table:
        df_ = parse_gid_file(highlight_geneids_table)
        # row_annot_track.append(df_)
    # #xx = containers.add_annotations(X)
    # this needs rewriting
    if annotate_genes:
        aa = containers.get_annotation_mapper()
        annot_mapping = aa.map_gene_ids(X.GeneID.tolist(), taxon="10090")
        row_annot_track.append(
            annot_mapping[
                [x for x in aa.df if x != "GeneSymbol" and x != "MitoCarta_Pathways"]
            ].set_index("GeneID")
        )
        outname_kws["ganno"] = ""

    row_annot_df = None
    if row_annot_track:
        if len(row_annot_track) > 1:  # more than 1 track to add
            row_annot_df = pd.concat(row_annot_track, axis=1)
        else:  # if length one it duplicate the whole thing?
            row_annot_df = row_annot_track[0]
        # make a dataframe that spans all genes about to be plotted
        ixs_ = X.GeneID.astype(str)
        missing_ = list(set(ixs_) - set(row_annot_df.index))
        intersect_ixs_ = list(set(ixs_) & set(row_annot_df.index))
        missing_df_ = pd.DataFrame(index=missing_, columns=row_annot_df.columns)
        row_annot_df = row_annot_df.loc[intersect_ixs_]
        row_annot_df = pd.concat([row_annot_df, missing_df_], axis=0).fillna("")

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
        # Normalize to a list of strings
        if isinstance(cut_by, str):
            cut_by = cut_by.split(":") if ":" in cut_by else [cut_by]
        else:
            cut_by = list(cut_by)
        # Validate columns exist in col_meta

        if col_meta is not None:
            for c in cut_by:
                if c not in col_meta.columns:
                    print(f"{c} not in column metadata")
                    cut_by = None
                    break
            if hasattr(cut_by, "__iter__") and cut_by is not None:
                outname_kws["cut"] = str.join("_", cut_by)
            elif cut_by is not None:
                outname_kws["cut"] = cut_by
            if cut_by is not None:
                cut_by = robjects.vectors.StrVector(
                    cut_by
                )  # actually make it the correct r object type
    elif cut_by is None:
        cut_by = None  # robjects.NULL # this is fine we make it robjects.NULL later
    if not cut_by:
        cut_by = None
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

        missing = set(X.GeneID) - set(annot_mat.index)
        if missing:
            _missing = pd.DataFrame(index=list(missing), columns=X.columns).fillna(0)
            _missing.index.name = X.index.name
            annot_mat = pd.concat([annot_mat, _missing])

        annot_mat["GeneID"] = annot_mat.index.astype(str)
        annot_mat = annot_mat[["GeneID"] + [x for x in annot_mat if x != "GeneID"]]
        # fill

    # ============================================================

    if z_score != "None" and z_score is not None:
        outname_kws["z"] = z_score
    else:
        outname_kws["z"] = "none"
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
    _d = {
        # "fna"+"T" if z_score_fillna else "F": "",
        # "mval"+"T" if missing_values else "F": "",
        "l"
        + linkage: ""
    }
    outname_func = partial(
        get_outname,
        # name=data_obj.outpath_name,
        tx=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        batch=data_obj.batch_applied,
        batch_method="param" if not data_obj.batch_nonparametric else "nonparam",
        # fna="T" if z_score_fillna else "F",
        # l=linkage,
        # mval="T" if missing_values else "F",
        norm=data_obj.normtype,
        outpath=data_obj.outpath,  # , "cluster2")
        **_d,
    )
    # =================================================================
    if gene_symbols and add_description:
        # this could definitely be improved
        description_frame = data_obj.data.query("Metric == 'GeneDescription'")
        _cols = [x for x in description_frame if x not in ("GeneID", "Metric")]
        # _descriptions = description_frame[_cols].stack().unique()
        _descriptions = description_frame[_cols].fillna("").stack().unique()
        if all(x == "" for x in _descriptions):
            description_frame = data_obj.data.query(
                "Metric == 'Description'"
            )  # try another
            _descriptions = description_frame[_cols].fillna("").stack().unique()
        if len(description_frame) == 0 or all(x == "" for x in _descriptions):
            pass  # fail
        else:
            description_frame = description_frame.bfill(axis=1).ffill(axis=1).fillna("")
            cols_for_merge = ["GeneID", description_frame.columns[-1]]
            to_merge = description_frame[cols_for_merge].rename(
                columns={description_frame.columns[-1]: "Description"}
            )
            X = pd.merge(X, to_merge, how="left")
            # try to remove excess info from description that are not useful for display
            # this is uniprot based removal of identifiers
            X["Description"] = X.Description.fillna("").str.replace(
                r"OS=.*", "", regex=True
            )
            X["GeneSymbol"] = (
                X.GeneSymbol.astype(str) + " " + X.Description.fillna("").astype(str)
            )
            X = X.drop("Description", axis=1)
        outname_kws["desc"] = "T"

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
        col_data = col_meta# .fillna("NA")
    #     # cannot convert Categorical column of Integers to Category in py2r
    #     col_data = col_meta.pipe(clean_categorical)  # does this fix the problem?
    # col_data = utils.set_pandas_datatypes(col_data)

    # rpy2 does not map None to robjects.NULL
    if row_annot_df is None:
        row_annot_df = robjects.NULL

    # pandas2ri.activate()

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
        nonlocal row_cluster, col_cluster
        """
        # gr_kws is redefined here (update fighiehgt /width), no need to take it as input
        """
        print("Saving", out, "...", end="", flush=True)
        # have to put this here to preserve the layout set by ComplexHeatmap::draw
        if annot_mat is None:
            annot_mat = robjects.NULL
        if main_title is None:
            main_title = robjects.NULL

        # Centralized sizing logic for clarity and easier troubleshooting
        def _base_dims(n_rows, n_cols, add_title):
            """Return baseline height/width in inches for the heatmap area."""
            # More generous scaling for large N
            row_slope = 0.24 if n_rows <= 250 else 0.18
            # Dial back width slope and baseline a bit further
            col_slope = 0.24 if n_cols <= 35 else 0.18
            base_h = 4.2 + (n_rows * row_slope)
            if add_title:
                base_h += 0.36
            base_w = 7.6 + (n_cols * col_slope)
            return base_h, base_w

        def _apply_extras(h, w, *, col_cluster_flag, row_annot_cols, add_desc,
                          margins_overhead, row_annot_side, row_names_side,
                          gene_symbols):
            info = {"col_cluster": 0.0, "row_annot_w": 0.0, "row_annot_h": 0.0,
                    "desc": 0.0, "margins": 0.0, "left_row_annot": 0.0,
                    "left_row_names": 0.0}
            if col_cluster_flag:
                h += 3.6  # dendrogram/labels overhead when clustering columns
                info["col_cluster"] = 3.6
            if row_annot_cols is not None and row_annot_cols > 0:
                # Allowance for legends/labels from row annotations
                _w = 0.4 + (0.26 * row_annot_cols)
                _h = 0.2 + (0.40 * row_annot_cols)
                w += _w
                h += _h
                info["row_annot_w"] = _w
                info["row_annot_h"] = _h
            if add_desc:
                w += 1.8
                info["desc"] = 1.8
            # Fixed overhead for left/right margins, dendrograms, etc.
            w += margins_overhead
            info["margins"] = margins_overhead
            # If row annotations sit on the left, reserve extra space
            if row_annot_side == 'left' and row_annot_cols and row_annot_cols > 0:
                w += 1.0
                info["left_row_annot"] = 1.0
            # If row names are on the left and gene symbols are shown, add a bit more
            if (row_names_side == 'left') and gene_symbols:
                w += 0.8
                info["left_row_names"] = 0.8
            return h, w, info

        def _ensure_min_w(w, min_w=5.4):
            return max(w, min_w)

        def _clamp_w(w, min_w=5.4, max_w=48.0):
            return max(min(w, max_w), min_w)

        # Determine figure size based on provided figsize tuple
        n_rows = int(X.shape[0])
        n_cols = int(len([c for c in X.columns if c not in ("GeneID", "GeneSymbol")]))
        has_title = not (main_title == robjects.NULL)
        annot_cols = None
        if row_annot_df is not None and row_annot_df is not robjects.NULL:
            try:
                annot_cols = int(len(row_annot_df.columns))
            except Exception:
                annot_cols = None

        # Configurable knobs via environment (for quick tuning without code changes)
        try:
            width_scale = float(os.getenv("TACKLE_WIDTH_SCALE", "1.0"))
        except Exception:
            width_scale = 1.0
        try:
            min_w_env = float(os.getenv("TACKLE_MIN_FIGWIDTH", "5.4"))
        except Exception:
            min_w_env = 5.4
        try:
            max_w_env = float(os.getenv("TACKLE_MAX_FIGWIDTH", "48.0"))
        except Exception:
            max_w_env = 48.0
        try:
            margins_overhead = float(os.getenv("TACKLE_WIDTH_MARGIN_OVERHEAD", "1.2"))
        except Exception:
            margins_overhead = 1.2

        # Legend-aware width bump (bottom legends are horizontal)
        legend_groups = 1  # include the heatmap legend itself
        dense_units = 0
        if 'col_meta' in locals() and (col_meta is not None):
            meta_cols = [c for c in col_meta.columns if c != 'name']
            legend_groups += len(meta_cols)
            # penalize heavy cardinality
            for c in meta_cols:
                try:
                    dense_units += max(int(col_meta[c].nunique()) - 6, 0)
                except Exception:
                    pass
        if row_annot_df is not None and row_annot_df is not robjects.NULL:
            try:
                legend_groups += int(len(row_annot_df.columns))
                for c in list(row_annot_df.columns):
                    try:
                        dense_units += max(int(pd.Series(row_annot_df[c]).nunique()) - 6, 0)
                    except Exception:
                        pass
            except Exception:
                pass

        # Reserve width per legend column to avoid clipping (dialed back by default)
        try:
            LEGEND_COL_WIDTH_IN = float(os.getenv("TACKLE_LEGEND_COL_WIDTH", "0.1"))
        except Exception:
            LEGEND_COL_WIDTH_IN = 0.1
        legend_width_extra = (legend_groups * LEGEND_COL_WIDTH_IN) + (0.08 * dense_units)

        # Reserve some height for bottom legends (they sit below the plot)
        try:
            LEGEND_ROW_HEIGHT_IN = float(os.getenv("TACKLE_LEGEND_ROW_HEIGHT", "0.35"))
        except Exception:
            LEGEND_ROW_HEIGHT_IN = 0.35
        legend_height_extra = (legend_groups * LEGEND_ROW_HEIGHT_IN) + (0.04 * dense_units)

        # cut_by splits create extra headers/labels in R, give a small bump
        has_cut_by = 'cut_by' in locals() and (cut_by is not None and cut_by is not robjects.NULL)
        cut_by_extra = 0.4 if has_cut_by else 0.0

        # Case A: figsize provided with one dimension missing
        if figsize is not None and (figsize[0] is None) and (figsize[1] is not None):
            # width missing, compute width only
            if optimal_figsize:
                base_h, base_w = _base_dims(n_rows, n_cols, has_title)
            else:
                # Non-optimal: conservative defaults
                base_h, base_w = 10.56, _ensure_min_w(min(n_cols / 2.0, 16.0), min_w_env)
            h, w, extra_info = _apply_extras(base_h, base_w,
                                             col_cluster_flag=col_cluster,
                                             row_annot_cols=annot_cols,
                                             add_desc=add_description,
                                             margins_overhead=margins_overhead,
                                             row_annot_side=row_annot_side,
                                             row_names_side=row_names_side,
                                             gene_symbols=gene_symbols)
            w_pre = (w + legend_width_extra + cut_by_extra) * width_scale
            figwidth = _clamp_w(w_pre, min_w_env, max_w_env)
            figheight = float(figsize[1])
            figsize = (figwidth, figheight)

        elif figsize is not None and (figsize[1] is None) and (figsize[0] is not None):
            # height missing, compute height only
            if optimal_figsize:
                base_h, base_w = _base_dims(n_rows, n_cols, has_title)
            else:
                base_h, base_w = 10.56, _ensure_min_w(min(n_cols / 2.0, 16.0), min_w_env)
            h, w, extra_info = _apply_extras(base_h, base_w,
                                             col_cluster_flag=col_cluster,
                                             row_annot_cols=annot_cols,
                                             add_desc=add_description,
                                             margins_overhead=margins_overhead,
                                             row_annot_side=row_annot_side,
                                             row_names_side=row_names_side,
                                             gene_symbols=gene_symbols)
            # add legend height here since we are computing height
            h = h + legend_height_extra
            w_pre = (w + legend_width_extra + cut_by_extra) * width_scale
            figwidth = _clamp_w(float(figsize[0]), min_w_env, max_w_env)
            figheight = h
            figsize = (figwidth, figheight)

        # Case B: figsize not provided or both set to None
        elif (figsize is None) or (figsize[0] is None and figsize[1] is None):
            if optimal_figsize:
                base_h, base_w = _base_dims(n_rows, n_cols, has_title)
            else:
                # Historical defaults that tended to look balanced
                base_h, base_w = 12.0, _ensure_min_w(min(n_cols / 2.0, 16.0), min_w_env)
            h, w, extra_info = _apply_extras(base_h, base_w,
                                             col_cluster_flag=col_cluster,
                                             row_annot_cols=annot_cols,
                                             add_desc=add_description,
                                             margins_overhead=margins_overhead,
                                             row_annot_side=row_annot_side,
                                             row_names_side=row_names_side,
                                             gene_symbols=gene_symbols)
            # add legend height here since we are computing height
            h = h + legend_height_extra
            w_pre = (w + legend_width_extra + cut_by_extra) * width_scale
            figwidth = _clamp_w(w_pre, min_w_env, max_w_env)
            figheight = h
            figsize = (figwidth, figheight)

        # Case C: both width and height provided, respect them as-is
        figwidth, figheight = figsize

        # Ensure room for gene symbols only when not using optimal sizing
        if gene_symbols and not optimal_figsize:
            figheight = max(((gene_symbol_fontsize + 2) / 72) * n_rows, 12)
            if figheight > 218:  # cap figheight and scale font in R
                FONTSIZE = max(218 / figheight, 6)
                figheight = 218

        if X.shape[0] < 2:
            row_cluster, col_cluster = False, False
        call_kws = dict(
            data=X,
            color_low=color_low,
            color_mid=color_mid,
            color_high=color_high,
            annot_mat=annot_mat,
            the_annotation=annotate or robjects.NULL,
            z_score=(
                z_score if z_score != "None" and z_score is not None else robjects.NULL
            ),
            z_score_by=z_score_by or robjects.NULL,
            z_score_fillna=z_score_fillna,
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
            cluster_fillna=cluster_fillna,
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
            figwidth=figwidth or robjects.NULL,
        )
        if cut_by is not None and cut_by is not robjects.NULL:
            kws["cut_by"] = cut_by
        call_kws.update(kws)
        # logger.info(f"call_kws: {call_kws}")

        if len(out) > 299:  # quick fix
            out_path, out_name = os.path.split(out)
            out_name = (
                out_name.replace("treatment", "treat")
                .replace("clustermap", "cmap")
                .replace("normtype", "norm")
                .replace("genotype", "geno")
            )
            out = os.path.join(out_path, out_name)

        # Detailed breakdown to simplify troubleshooting
        try:
            logger.info(
                "figsize: base_w=%.2f, legends_w=%.2f, margins=%.2f, row_annot_w=%.2f, desc=%.2f, cut_by=%.2f, left_row_annot=%.2f, left_row_names=%.2f, width_scale=%.2f, pre_clamp=%.2f, final_w=%.2f, base_h=%.2f, legends_h=%.2f, final_h=%.2f"
                % (
                    base_w if 'base_w' in locals() else -1,
                    legend_width_extra,
                    extra_info.get("margins", 0.0) if 'extra_info' in locals() else 0.0,
                    extra_info.get("row_annot_w", 0.0) if 'extra_info' in locals() else 0.0,
                    extra_info.get("desc", 0.0) if 'extra_info' in locals() else 0.0,
                    cut_by_extra,
                    extra_info.get("left_row_annot", 0.0) if 'extra_info' in locals() else 0.0,
                    extra_info.get("left_row_names", 0.0) if 'extra_info' in locals() else 0.0,
                    width_scale,
                    w_pre if 'w_pre' in locals() else -1,
                    figwidth,
                    base_h if 'base_h' in locals() else -1,
                    legend_height_extra,
                    figheight,
                )
            )
        except Exception:
            pass

        with localconverter(
            robjects.default_converter + pandas2ri.converter
        ):  # no tuples
            grdevice = grdevice_helper.get_device(
                filetype=file_fmt, width=figwidth, height=figheight
            )
            grdevice(file=out)  # open file for saving
            ret = cluster2(**call_kws)  # draw
            # print(ret[0])
            grdevices.dev_off()  # close file
            print(".done", flush=True)

        sil_df = None
        try:
            sil_df = rpy2.robjects.pandas2ri.rpy2py_dataframe(ret[1])
        except Exception as e:
            pass

        if sil_df is not None:
            row_orders = ret[2]
            the_orders = [
                [x - 1 for x in row_orders.rx2(int(str(n)))] for n in row_orders.names
            ]
            the_orders = [x for y in the_orders for x in y]

            # [0]

            # the_orders = [
            #     row_orders.rx2(int(str(n))) - 1 for n in row_orders.names
            # ]  # subtract 1  for zero indexing

            cluster_metrics = sil_df.iloc[the_orders]

            out = outname_func("clustermap", **outname_kws) + ".tsv"
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
        outname_base_name = "clustermap"
        this_outname_kws = outname_kws.copy()
        if "vfile" in outname_kws:  # 'hack' to nest 1 folder deeper
            volcano_file = this_outname_kws.pop("vfile")
            outname_base_name = os.path.join(outname_base_name, volcano_file)
            # extra_outname_kws["vf"] = volcano_file
        if "genefile" in outname_kws:
            _genefile = this_outname_kws.pop("genefile")
            outname_base_name = os.path.join(outname_base_name, _genefile)

        # outname_kws.update(extra_outname_kws)
        # this_outname_kws = outname_kws.copy()  # <- defensive copy
        # this_outname_kws.update(extra_outname_kws)
        out = outname_func(outname_base_name, **this_outname_kws) + nrow_ncol + file_fmt

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
                # grdevice,
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
