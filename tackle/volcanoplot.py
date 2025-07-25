import os
import itertools

import numpy as np
import pandas as pd

try:
    pd.NA
except AttributeError:
    pd.NA = "NA"  # not sure if using this, could be in main.py

from .utils import get_outname, parse_gid_file, fix_name

# from .containers import GeneMapper


# try:
#     from rpy2.robjects import r
#     import rpy2.robjects as robjects
#     from rpy2.robjects import pandas2ri
#     from rpy2.robjects.packages import importr
#     _viaR = True
# except ModuleNotFoundError:
#     _viaR = False
#     print("Must install rpy2")


def volcanoplot(
    ctx,
    foldchange,
    expression_data,
    number,
    only_sig=False,
    sig=0.05,
    genes=None,
    direction="both",
    sig_metric="pAdj",
    number_by="log2_FC",
    yaxis="pAdj",
    label_scale=1.4,
    marker_scale=1.2,
    highlight_geneids=None,
    force_highlight_geneids=False,
    formula=None,
    contrasts=None,
    impute_missing_values=False,
    width=5,
    height=5,
    annot_scale=1.0,
    bg_marker_color="#22222288",
    pch=16,
    alpha=1.0,
    fill_na_zero=False,
    extra_outname_info=None,
    color_down="blue",
    color_up="red",
    global_xmax=None,
    global_ymax=None,
):
    data_obj = ctx.obj["data_obj"]


    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter


    # gm = GeneMapper()
    gm = get_gene_mapper()

    if yaxis not in ("pValue", "pAdj"):
        raise ValueError("Must choose between `pValue` and `pAdj`")

    group = data_obj.group  #
    if group is None and formula is None:
        print("Must supply a group value.")
        return

    # if group and data_obj.col_metadata[group].nunique() < 2:
    #     print("Error in volcanoplot, number of groups must be at least 2.")
    #     return
    # if data_obj.col_metadata.loc[group].nunique() != 2:
    #     print('Error in volcanoplot, number of groups must be exactly 2.')
    #     return
    if group and data_obj.col_metadata[group].value_counts().min() < 2:
        print("Each group must have at least 2 replicates.")
        return

    if group:
        groups = dict()
        for grp in data_obj.col_metadata[group].unique():
            samples = (
                (data_obj.col_metadata[group] == grp).where(lambda x: x).dropna().index
            )
            groups[grp] = samples

    # this is where the missing value imputation results will come from
    results = data_obj.stat_model(
        formula=formula,
        contrasts_str=contrasts,
        impute_missing_values=impute_missing_values,
        fill_na_zero=fill_na_zero,
    )

    meta = data_obj.col_metadata

    def fix_group_name(group, entries):
        # for entry in meta.index:
        # print(group, entries)
        group = group.split("_")
        entries = sorted(entries, key=lambda x: len(x))
        for ix, entry in enumerate(entries):
            groupres = list()
            for g in group:
                if (
                    g.startswith(entry)
                    and not any(g.startswith(e) for e in entries[ix + 1 :])
                    and False
                ):  # ignore this
                    i = g.find(entry)
                    if i != 0:
                        continue
                    res = g[len(entry) :]
                    # res = g.lstrip(entry)
                else:
                    res = g
                groupres.append(res)
            group = groupres
            # group = [x.lstrip(entry) if x.startswith(entry) else x for x in group]
            # print(entry, group)
        return ":".join(group)

        # start = group.find(entry)
        # end = len(entry) + start
        # if start == -1:
        #     continue
        # group_lst = list(group)
        # group = ''.join(group_lst[:start] + group_lst[end:])
        # return group

    max_fc = float(max(x.log2_FC.abs().max() for x in results.values()))

    def _clean(x):
        return x.strip("\( ").strip(" \)").strip()

    # prepare table for volcanoplot and export
    if genes is not None:  # only plot these select genes
        _genes = set(genes) & set(df.index)
        df = df.loc[_genes]

    if global_xmax is None or global_xmax is True:
        global_xmax = float(abs(max(x["log2_FC"].abs().max() for x in results.values())))
    elif global_xmax is False:
        global_xmax = robjects.NA_Real
    if global_ymax is None or global_ymax is True:
        global_ymax = -np.log10(
            min(x["pValue"][x["pValue"] > 0].min() for x in results.values())
        )
        global_ymax = float(global_ymax)
        # Use min to find the smallest p-value, then take the -log10. Exclude non-positive p-values.
    elif global_ymax is False:
        global_ymax = robjects.NA_Real

    # TODO add check to ensure xmax and ymax  not smaller than all actual x and y values
    # now this
    for comparison, df in results.items():  # results contains dataf
        # df already contains the expression data by default now, as returned by the data_obj stat running method

        _comparison_name = None
        _comparison = comparison
        if "=" in _comparison:
            _comparison_name, _comparison = _comparison.split("=", maxsplit=1)

        # too nested
        # group0, group1 establish
        groups = _comparison.split(" - ")
        if len(groups) == 2:
            # group0, group1 = [x.strip() for x in comparison.split('-')]
            group1, group0 = [_clean(x) for x in _comparison.split(" - ")]
            # print(group0, group1)
            group0_fix, group1_fix = (
                fix_group_name(group0, meta.columns),
                fix_group_name(group1, meta.columns),
            )
            # print(group0, group1)
        else:
            # more complex formula
            # TODO handle this, maybe with user-config?
            # pass
            # group1, group0 = [x.strip() for x in comparison.split("-")]
            group1, group0 = [x.strip() for x in _comparison.split(" - ")]
            group0_fix, group1_fix = (
                fix_group_name(group0, meta.columns),
                fix_group_name(group1, meta.columns),
            )
            group0_fix, group1_fix = "Down", "Up"

        # if we do this here, the columns get added before the expression values
        df["GeneSymbol"] = df.index.map(
            lambda x: data_obj.gid_symbol.get(x, gm.symbol.get(str(x), x))
        )
        df["FunCats"] = df.index.map(lambda x: data_obj.gid_funcat_mapping.get(x, ""))
        df["GeneDescription"] = df.index.map(lambda x: gm.description.get(str(x), ""))
        df.index.name = "GeneID"
        df["highlight"] = False
        if highlight_geneids is not None:
            df.loc[df.index.intersection(highlight_geneids), "highlight"] = True
        #     df.
        # if highlight_geneids

        df["signedlogP"] = df.apply(
            lambda x: -np.log10(x["pValue"]) * (1 if x["log2_FC"] > 0 else -1), axis=1
        ).fillna(0)

        # # df['GeneSymbol'] = df.index.map(lambda x: data_obj.gid_symbol.get(x, '?'))
        # df["GeneSymbol"] = df.index.map(
        #     lambda x: data_obj.gid_symbol.get(x, gm.symbol.get(str(x), x))
        # )
        # df["FunCats"] = df.index.map(lambda x: data_obj.gid_funcat_mapping.get(x, ""))
        # df["GeneDescription"] = df.index.map(lambda x: gm.description.get(str(x), ""))
        # df.index.name = "GeneID"
        # df["highlight"] = False
        # df["signedlogP"] = df.apply(
        #     lambda x: -np.log10(x["pValue"]) * (1 if x["log2_FC"] > 0 else -1), axis=1
        # )
        # if highlight_geneids is not None:
        #     df.loc[set(highlight_geneids) & set(df.index), "highlight"] = True
        _xtra = {}
        direction_codes = {"up": "U", "down": "D", "both": "B"}
        # if direction != "both":
        _xtra["dir"] = direction_codes.get(direction, "?")

        _b = None
        if not data_obj.batch_nonparametric and data_obj.batch_applied == True:
            _b = "parametric"
        elif data_obj.batch_applied == True:
            _b = ("nonparam",)

        # impute_missing_values = (impute_missing_values,)
        # fill_na_zero = (fill_na_zero,)
        if _comparison_name is None:
            _groupname = "{}_vs_{}".format(group0_fix, group1_fix)
        else:
            _groupname = _comparison_name

        outname = get_outname(
            "volcano",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            # batch_me=_b,
            #sort=number_by,
            colors_only=data_obj.colors_only,
            batch=data_obj.batch_applied,
            normtype=data_obj.normtype,
            imv="T" if impute_missing_values else "F",
            fna="T" if fill_na_zero else "F",
            group=_groupname.replace(":", "_").replace(r"(", "").replace(")", "").replace("+", "_").replace(r"/", "dv"),
            outpath=os.path.join(data_obj.outpath, "volcano"),
            **_xtra,
        )

        slicepoint = 170
        space = min(10, 255 - slicepoint)  # Adjust space to fit within max length
        if len(outname) > 255:
            outname = outname.replace("timepoint", "T")
            outname = outname.replace("time", "T")
            outname = outname.replace("genotype", "geno")

        while len(outname) > 255:
            # Ensure you're slicing within the bounds of the string
            if slicepoint + space >= len(outname):
                space = len(outname) - slicepoint - 1  # Avoid out-of-bounds slicing

            outname = outname[:slicepoint] + ".." + outname[slicepoint + space :]
            slicepoint += 20
            space = min(10, 255 - slicepoint)  # Adjust space to fit within max length

        out = outname + ".tsv"
        print("Saving", out, "...", end="", flush=True)
        export_data = df  # this is a "results" dataframe from the stat running method
        if only_sig:
            _log2_cutoff = np.sqrt(foldchange)
            export_data = df.query("pAdj < @sig & abs(log2_FC) > @_log2_cutoff")
        if expression_data:
            try:
                group_entries = [group0.split(group)[-1], group1.split(group)[-1]]
                exps = data_obj.col_metadata[
                    data_obj.col_metadata[group].isin(group_entries)
                ].index
                if exps.empty:
                    raise ValueError()
                export_data = export_data.join(data_obj.areas_log_shifted[exps])
            except Exception as e:
                print(
                    "Error trying to subselect data for export. Try specifying group if you have not."
                )
                export_data = export_data.join(data_obj.areas_log_shifted / np.log10(2))

        export_cols = [x for x in export_data.columns if x not in ("highlight",)]
        export_data[export_cols].to_csv(out, sep="\t")
        print("done", flush=True)

        # pandas2ri.activate()
        r_source = robjects.r["source"]
        r_file = os.path.join(
            os.path.split(os.path.abspath(__file__))[0], "R", "volcanoplot.R"
        )
        r_source(r_file)
        Rvolcanoplot = robjects.r["volcanoplot"]

        file_fmts = ctx.obj["file_fmts"]
        grdevices = importr("grDevices")
        gr_devices = {
            ".png": grdevices.png,
            ".pdf": grdevices.pdf,
            ".svg": grdevices.svg,
        }
        gr_kws = {
            ".png": dict(width=width, height=height, units="in", res=300),
            ".pdf": dict(
                width=width,
                height=height,
            ),
            ".svg": dict(
                width=width,
                height=height,
            ),
        }

        df["FunCats"] = df.FunCats.fillna("").astype(str)
        df["GeneSymbol"] = df["GeneSymbol"].fillna("").astype(str)
        df.index = df.index.astype(
            "str"
        )  # make uniform dtype so rpy2 does not crash in conversion
        df = df.reset_index()
        df = df[~df["t"].isna()]

        # make new outname

        if extra_outname_info is not None:
            _xtra["n"] = extra_outname_info

        outname = get_outname(
            "volcano",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            # batch_me=_b,
            sort=number_by,
            colors_only=data_obj.colors_only,
            batch=data_obj.batch_applied,
            normtype=data_obj.normtype,
            imv="T" if impute_missing_values else "F",
            fna="T" if fill_na_zero else "F",
            group=_groupname.replace(":", "_").replace(r"(", "").replace(")", "").replace("+", "_").replace(r"/", "dv"),
            outpath=os.path.join(data_obj.outpath, "volcano"),
            **_xtra,
        )
        for file_fmt in file_fmts:
            with localconverter(robjects.default_converter + pandas2ri.converter): # no tuples
                grdevice = gr_devices[file_fmt]
                gr_kw = gr_kws[file_fmt]
                out = outname + file_fmt
                print("Saving", out, "...", end="", flush=True)

                grdevice(file=out, **gr_kw)

                # Rvolcanoplot(pandas2ri.py2ri(df.reset_index()), max_labels=number, fc_cutoff=foldchange,
                # _data = df.reset_index()
                # _data['FunCats'] = _data.FunCats.fillna('')
                Rvolcanoplot(
                    df,
                    max_labels=number,
                    fc_cutoff=foldchange,
                    number_by=number_by,
                    direction=direction,
                    force_highlight_geneids=force_highlight_geneids,
                    bg_marker_color=bg_marker_color,
                    color_down=color_down,
                    color_up=color_up,
                    sig=sig,
                    sig_metric=sig_metric,
                    yax=yaxis,
                    label_cex=label_scale,
                    annot_cex=annot_scale,
                    marker_cex=marker_scale,
                    max_fc=max_fc,
                    point_size=1.4,
                    group0=group0,
                    group1=group1,
                    alpha=alpha,
                    pch=pch,
                    global_xmax=global_xmax,
                    global_ymax=global_ymax,
                    # **kws,
                )
                grdevices.dev_off()
                print("done.", flush=True)
        ## end for comparison, df in results.items():

    return


from .containers import get_gene_mapper

# =============================================================================================
# end
# =============================================================================================
