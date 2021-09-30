import os
import itertools


import numpy as np
import pandas as pd

try:
    pd.NA
except AttributeError:
    pd.NA = "NA"

from .utils import get_outname, parse_gid_file, fix_name

from .containers import GeneMapper


def volcanoplot(
    ctx,
    foldchange,
    expression_data,
    number,
    only_sig=False,
    sig=0.05,
    genes=None,
    sig_metric="pAdj",
    number_by="log2_FC",
    yaxis="pAdj",
    scale=1.2,
    highlight_geneids=None,
    formula=None,
    contrasts=None,
    impute_missing_values=False,
    width=5,
    height=5,
):
    print(sig_metric)

    data_obj = ctx.obj["data_obj"]
    gm = GeneMapper()

    if yaxis not in ("pValue", "pAdj"):
        raise ValueError("Must choose between `pValue` and `pAdj`")

    group = data_obj.group  #
    if group is None and formula is None:
        print("Must supply a group value.")
        return

    if group and data_obj.col_metadata[group].nunique() < 2:
        print("Error in volcanoplot, number of groups must be at least 2.")
        return
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

    results = data_obj.stat_model(
        formula=formula,
        contrasts_str=contrasts,
        impute_missing_values=impute_missing_values,
    )

    meta = data_obj.col_metadata
    # import ipdb; ipdb.set_trace()

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

    max_fc = max(x.log2_FC.abs().max() for x in results.values())

    for comparison, df in results.items():
        # group0, group1 establish
        groups = comparison.split(" - ")
        if len(groups) == 2:
            # group0, group1 = [x.strip() for x in comparison.split('-')]
            group1, group0 = [x.strip() for x in comparison.split(" - ")]
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
            # import ipdb; ipdb.set_trace()
            group1, group0 = [x.strip() for x in comparison.split(" - ")]
            group0_fix, group1_fix = (
                fix_group_name(group0, meta.columns),
                fix_group_name(group1, meta.columns),
            )
            group0_fix, group1_fix = "Down", "Up"

        # df['GeneSymbol'] = df.index.map(lambda x: data_obj.gid_symbol.get(x, '?'))
        df["GeneSymbol"] = df.index.map(
            lambda x: data_obj.gid_symbol.get(x, gm.symbol.get(str(x), x))
        )
        df["FunCats"] = df.index.map(lambda x: data_obj.gid_funcat_mapping.get(x, ""))
        df["GeneDescription"] = df.index.map(lambda x: gm.description.get(str(x), ""))
        df.index.name = "GeneID"
        df["highlight"] = False
        df["signedlogP"] = df.apply(
            lambda x: -np.log10(x["pValue"]) * (1 if x["log2_FC"] > 0 else -1), axis=1
        )
        if highlight_geneids is not None:
            df.loc[set(highlight_geneids) & set(df.index), "highlight"] = True

        if genes is not None:  # only plot these select genes
            _genes = set(genes) & set(df.index)
            df = df.loc[_genes]

        outname = get_outname(
            "volcanoplot",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            colors_only=data_obj.colors_only,
            batch=data_obj.batch_applied,
            batch_method="parametric"
            if not data_obj.batch_nonparametric
            else "nonparametric",
            outpath=data_obj.outpath,
            group="{}_vs_{}".format(fix_name(group0_fix), fix_name(group1_fix)),
        )

        if len(outname) > 255:
            outname = outname[:255]
        out = outname + ".tsv"

        print("Saving", out, "...", end="", flush=True)
        export_data = df
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
                export_data = export_data.join(data_obj.areas_log_shifted)
        export_cols = [x for x in export_data.columns if x not in ("highlight",)]
        export_data[export_cols].to_csv(out, sep="\t")
        print("done", flush=True)

        try:
            from rpy2.robjects import r
            import rpy2.robjects as robjects
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr

            _viaR = True
        except ModuleNotFoundError:
            _viaR = False

        if not _viaR:
            print("Must install rpy2")
            return

        pandas2ri.activate()
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
        df["GeneSymbol"] = df.GeneSymbol.fillna("").astype(str)
        df.index = df.index.astype(
            "str"
        )  # make uniform dtype so rpy2 does not crash in conversion
        df = df.reset_index()

        for file_fmt in file_fmts:

            grdevice = gr_devices[file_fmt]
            gr_kw = gr_kws[file_fmt]
            out = outname + file_fmt
            print("Saving", out, "...", end="", flush=True)

            grdevice(file=out, **gr_kw)

            # Rvolcanoplot(pandas2ri.py2ri(df.reset_index()), max_labels=number, fc_cutoff=foldchange,
            Rvolcanoplot(
                df.reset_index(),
                max_labels=number,
                fc_cutoff=foldchange,
                number_by=number_by,
                sig=sig,
                sig_metric=sig_metric,
                yaxis=yaxis,
                label_cex=scale,
                max_fc=max_fc,
                group0=group0,
                group1=group1,
            )

            grdevices.dev_off()
            print("done.", flush=True)

    return

    # =============================================================================================
    # end
    # =============================================================================================
