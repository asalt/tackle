import os
from collections import OrderedDict

import matplotlib

from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib import markers
from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import seaborn as sb


from .utils import get_outname, save_multiple, genefilter, filter_sra, filter_taxon

from .containers import TAXON_MAPPER

idx = pd.IndexSlice

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"


def make_metrics(
    data_obj, file_fmts, before_filter=False, before_norm=False, full=False
):

    rc = {
        "font.family": "sans-serif",
        "font.sans-serif": [
            "DejaVu Sans",
            "Arial",
            "Liberation Sans",
            "Bitstream Vera Sans",
            "sans-serif",
        ],
        "legend.frameon": True,
    }

    sb.set_context("talk")
    sb.set_palette("muted")
    sb.set_color_codes()
    sb.set_style("white", rc)

    data = data_obj.metric_values

    if before_filter:
        kws = dict(before="filter")
        taxon = "all"
    else:
        kws = dict(after="filter")
        taxon = data_obj.taxon

    outname = get_outname(
        "metrics",
        name=data_obj.outpath_name,
        taxon=taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        outpath=data_obj.outpath,
        **kws
    )

    sra = pd.DataFrame(data=[data[n]["SRA"] for n in data.keys()], index=data.keys())
    gpg = pd.DataFrame(
        data=[data[n]["GPGroups"] for n in data.keys()], index=data.keys()
    )
    psms = pd.DataFrame(data=[data[n]["PSMs"] for n in data.keys()], index=data.keys())
    peptides = pd.DataFrame(
        data=[data[n]["Peptides"] for n in data.keys()], index=data.keys()
    )
    area = OrderedDict((n, data[n]["Area"]) for n in data.keys())

    # ========================================================================

    trypsin = OrderedDict((n, data[n]["Trypsin"]) for n in data.keys())
    trypsin = pd.DataFrame(trypsin)
    trypsin_df = trypsin.apply(lambda x: x / sum(x))

    trypsinP = OrderedDict((n, data[n]["Trypsin/P"]) for n in data.keys())
    trypsinP = pd.DataFrame(trypsinP).fillna(0)
    trypsinP_df = trypsinP.apply(lambda x: x / sum(x))
    # trypsinP_df.index.name = 'miscuts'
    # trypsinP_df = trypsin_df.reset_index()

    trypsin_dfr = trypsin_df.melt(
        var_name="name", value_name="Trypsin", ignore_index=False
    )
    trypsin_dfr.index.name = "miscuts"
    trypsin_dfr = trypsin_dfr.reset_index()
    trypsinP_dfr = trypsinP_df.melt(
        var_name="name", value_name="Trypsin/P", ignore_index=False
    )
    trypsinP_dfr.index.name = "miscuts"
    trypsinP_dfr = trypsinP_dfr.reset_index()

    miscut_frame = pd.DataFrame.merge(
        trypsin_dfr, trypsinP_dfr, on=["miscuts", "name"], how="outer"
    ).sort_values(by=["miscuts", "name"])

    export_name = get_outname(
        "miscut_ratio",
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        outpath=data_obj.outpath,
    )
    miscut_frame.to_csv(export_name + ".tsv", sep="\t", index=False)

    # ========================================================================

    frames = list()
    area_name = "AreaSum_dstrAdj" if before_norm else "AreaSum_dstrAdj_normed"
    for name, value in area.items():
        frame = pd.Series(value).to_frame(area_name).multiply(1e9).apply(np.log10)
        frame["Name"] = name
        frames.append(frame)
    area_df = pd.concat(frames)

    to_export = (
        pd.concat(
            [
                sra[["S", "R", "A"]],
                gpg.rename(columns={0: "GPGroups"}),
                psms.rename(columns={x: x + "_psms" for x in psms.columns}).fillna(0),
                peptides.rename(columns={x: x + "_peptides" for x in psms.columns}),
            ],
            axis=1,
        )
        .fillna(0)
        .astype(int)
    )

    export_name = get_outname(
        "metrics",
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        outpath=data_obj.outpath,
    )
    to_export.to_csv(export_name + ".tab", sep="\t", index=True)

    # ==================================================================
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], "R", "metrics.R")
    r_source(r_file)
    Rmetrics = robjects.r["metrics"]
    Rmetrics(to_export, savename=export_name, exts=[x.lstrip(".") for x in file_fmts])
    # ==================================================================

    ggridges = importr("ggridges")
    rboxplot = r["boxplot"]
    rboxplot

    import rpy2.robjects.lib.ggplot2 as gg

    area_df['Name'] = pd.Categorical(area_df['Name'], ordered=True)

    plot = (
        gg.ggplot(area_df)
        + gg.aes_string(y="Name", x=area_name)
        + ggridges.stat_density_ridges(quantile_lines=True, alpha=0.8)
        + gg.theme_classic(base_size=12)
    )
    plot.plot()

    outname = get_outname(
        "metrics_dist",
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        outpath=data_obj.outpath,
        after="filter"
        # **kws
    )
    for ffmt in file_fmts:
        plot.save(outname + ffmt)

    # area = pd.DataFrame(data=[data[n]['area'] for n in data.keys()], index=data.keys())

    green = "darkgreen"
    yellow = "gold"
    red = "firebrick"

    # fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,8), sharex=True, sharey=False)
    # fig = plt.figure(figsize=(18,10))
    # gs = gridspec.GridSpec(2,3)

    # ax_sra  = fig.add_subplot(gs[0, 0])
    # ax_gpg  = fig.add_subplot(gs[0, 1])
    # ax_psms = fig.add_subplot(gs[1, 0])
    # ax_pept = fig.add_subplot(gs[1, 1])
    # ax_area = fig.add_subplot(gs[0:, 2])

    # sra[['S', 'R', 'A']].plot.bar(stacked=True, ax=ax_sra, color=[green, yellow, red], title='SRA')
    # gpg.plot.bar(ax=ax_gpg, legend=False, title='Gene Product Groups')
    # psms[['Total', 'u2g']].plot.bar(stacked=False, ax=ax_psms, title='PSMs', width=.75)
    # peptides[['Total', 'u2g', 'Strict', 'Strict_u2g']].plot.bar(stacked=False, ax=ax_pept, title='Peptides',
    #                                                             width=.8)
    # plt.setp(ax_sra.get_xticklabels(), visible=False)
    # plt.setp(ax_gpg.get_xticklabels(), visible=False)

    # # sb.violinplot(y='Name', x=area_name, data=area_df, ax=ax_area)
    # ### ????? This is not working for some reason...... n
    # ### *** TypeError: ufunc 'isfinite' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''

    # area_df['Name'] = area_df['Name'].astype(str)
    # area_df[area_name] = area_df[area_name].astype(float)
    # # sb.boxenplot(y='Name', x=area_name, data=area_df, ax=ax_area)
    # sb.boxenplot(y=area_df['Name'].values, x=area_df[area_name].values, ax=ax_area)
    # ax_area.yaxis.tick_right()
    # ax_area.set_ylabel('')
    # if before_norm:
    #     ax_area.set_xlabel('log$_{10}$ iBAQ dstrAdj')
    # else:
    #     ax_area.set_xlabel('log$_{10}$ iBAQ dstrAdj normed')
    # # plt.setp( ax_area.xaxis.get_majorticklabels(), rotation=90 )
    # plt.setp( ax_area.xaxis.get_majorticklabels(), rotation=0 )

    # FONTSIZE = 10
    # ncols = {0: 1, 2:2, 3:1, 4:0}
    # for ix, ax in enumerate((ax_sra, ax_gpg, ax_psms, ax_pept, ax_area)):
    #     ax.yaxis.grid(True, lw=.25, color='grey', ls=':')
    #     if ix == 4:
    #         ticklabels = ax.yaxis.get_ticklabels()
    #     else:
    #         ticklabels = ax.xaxis.get_ticklabels()
    #     for tick in ticklabels:
    #         txt = tick.get_text()
    #         newsize = FONTSIZE + 3 if ix == 4 else FONTSIZE
    #         textlen = len(txt)

    #         # resizing so labels fit, not best approach
    #         if textlen > 7 and textlen <= 9:
    #             newsize -= 1
    #         if textlen >= 10 and textlen < 13:
    #             newsize -= 1
    #         elif textlen >= 13:
    #             newsize -= 1
    #         if len(data.keys()) >= 17 and ix != 4:
    #             newsize -= 2
    #         tick.set_size(newsize)
    #         # if ix == 4:
    #         #     print(txt, textlen, newsize)
    #     sb.despine(ax=ax)
    #     if ix == 1:
    #         continue  #  no legend for gpgroups
    #     ax.legend(loc='lower left', ncol=ncols[ix], fontsize=8)

    # sb.despine(ax=ax_area, right=False, top=True, left=True, bottom=False)
    # ax_area.xaxis.grid(True, lw=.25, color='grey', ls=':')

    # fig.subplots_adjust(hspace=.15, top=.95, left=.1, right=.9)
    # save_multiple(fig, outname, *file_fmts)

    if full:  # also plot info from PSMs
        raise NotImplementedError("Not yet")
        return  # not done
        config = data_obj.config
        data_dir = data_obj.data_dir

        psms = OrderedDict()
        for name, record in config.items():
            if name.startswith("__"):
                continue
            recno = record.get("recno")
            runno = record.get("runno")
            searchno = record.get("searcno")
            df = ispec.PSMs(recno, runno, searchno, data_dir=data_dir).df.query(
                "oriFLAG==1"
            )

            df.index = pd.to_timedelta(df.RTmin, unit="m")

            psms[name] = df

        nrows = len(psms)

        m = max([pd.Timedelta(df.RTmin.max(), unit="m") for df in psms.values()])

        fig = plt.figure(figsize=(12, 8))
        gs = gridspec.GridSpec(
            nrows, 3, width_ratios=(3, 1, 3), left=0.1, right=0.9, bottom=0.1, top=0.9
        )
        # plot IDs/min over time
        prev_ax = None
        for ix, (name, df) in enumerate(psms.items()):
            if prev_ax:
                # ax = fig.add_subplot(gs[ix, 0], sharex=prev_ax, sharey=prev_ax)
                ax = fig.add_subplot(gs[ix, 0], sharey=prev_ax)
            else:
                ax = fig.add_subplot(gs[ix, 0])
            # df.index = pd.to_timedelta( df.RTmin, unit='m' )
            df.loc[m] = np.NaN
            g = df.groupby(
                [
                    pd.Grouper(freq="Min"),
                ]
            )
            g.size().plot(ax=ax)
            label_ax = fig.add_subplot(gs[ix, 1])
            label_ax.annotate(
                name,
                xy=(0.5, 0.5),
                xycoords="axes fraction",
                va="center",
                ha="center",
                size=12,
            )
            sb.despine(ax=label_ax, left=True, bottom=True)
            plt.setp(label_ax.get_xticklabels(), visible=False)
            plt.setp(label_ax.get_yticklabels(), visible=False)
            ax.xaxis.grid(True, lw=0.25, color="grey", ls=":")
            if ix < nrows - 1:
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.set_xlabel("")
            sb.despine(ax=ax)
            prev_ax = ax

    # ==========================  gene count overlap ============================================
    bins = (
        0.001,
        0.10,
        0.15,
        0.20,
        0.25,
        0.3,
        0.40,
        0.50,
        0.60,
        0.75,
        0.80,
        0.90,
        1.00,
    )
    ncols = len(data_obj.data.columns)
    fracs = [x / ncols for x in range(1, ncols + 1)]
    bin_indices = sorted(set(np.searchsorted(bins, fracs)))
    bins_kept = [0, *[bins[x] for x in bin_indices], 1.1]

    df = data_obj.data.pipe(filter_sra, SRA="S")
    ## Below doesn't work, always returns after filter since before_filter data not saved
    ## for this purpose.
    ## will require explicit saving like the other data.
    # if before_filter:
    #     df = data_obj.data.pipe(filter_sra, SRA='S')
    # elif not before_filter:  #after funcat filter
    #     dfs = dict()
    #     for c in data_obj.data.columns:
    #         frame = (data_obj.data[[c]].reset_index().pivot('level_0', 'level_1', c)
    #                  .pipe(genefilter,
    #                        funcats=data_obj.funcats,
    #                        funcats_inverse=data_obj.funcats_inverse,
    #                        geneid_subset=data_obj.geneid_subset,
    #                        ignored_geneid_subset=data_obj.ignore_geneid_subset
    #                  )
    #         )
    #         dfs[c] = frame
    #     stacked_data = [ df.stack() for df in dfs.values() ]
    #     df = (pd.concat( stacked_data, axis=1, keys=dfs.keys() ).pipe(filter_sra, SRA='S')
    #                  .pipe(filter_taxon, taxon=TAXON_MAPPER.get(data_obj.taxon))
    #     )

    # counts = (df.loc[ idx[:, 'SRA'], :]
    counts = df.loc[df.Metric == "SRA"].where(lambda x: x == "S").count(1)

    count_frac = counts / ncols
    bin_labels = ["≥ {:.0%}".format(x) for x in bins_kept[:-1]]
    bin_labels[0] = "> 0%"
    bin_labels[-1] = "100%"
    # bin_labels.append('100%')
    count_frac_binned = pd.cut(
        count_frac, bins_kept, include_lowest=True, right=False, labels=bin_labels
    )
    summary = (
        count_frac_binned.value_counts()
        .sort_index(ascending=False)
        .where(lambda x: x > 0)
        .dropna()
        .cumsum()
    )

    new_rc = {"xtick.labelsize": 14, "ytick.labelsize": 14, "font.size": 14}
    matplotlib.rcParams.update(new_rc)
    fig, ax = plt.subplots()
    summary.plot.bar(ax=ax, color="b")
    ax.set_ylabel("Gene Counts (Strict)")
    ax.grid(True, axis="y", linewidth=0.5, color="#444444", ls="--")
    # sb.set()
    # matplotlib.rc('xtick', labelsize=14)
    # matplotlib.rc('xtick', labelsize=14)
    # matplotlib.rc('ytick', labelsize=14)
    # matplotlib.rc('font', size=14)
    fig.autofmt_xdate(ha="center")
    fig.tight_layout()

    # Note change this to taxon=taxon when updated
    # also right now always returns after filter
    outname = get_outname(
        "metrics_genecounts",
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        outpath=data_obj.outpath,
        after="filter"
        # **kws
    )
    save_multiple(fig, outname, *file_fmts)

    # # ==========================  gene count overlap
    # # this is just not practical when the number of experiments grows.
    # # Need smarter filtering of inters_size_bounds and inters_degree_bounds
    # import pyupset
    # from pyupset.visualisation import plot as pu_plot
    # from pyupset.visualisation import UpSetPlot
    # data_dict = OrderedDict()
    # for col in data_obj.data:
    #     print(col)
    #     boolean = data_obj.data.loc[ idx[:, 'SRA'], col ] == 'S'
    #     gids = boolean.where(lambda x: x).dropna().index.get_level_values(0)
    #     data_dict[col] = pd.DataFrame({'GeneID': gids})

    # print('making fig')
    # fig_dict = pu_plot(data_dict, unique_keys=['GeneID'], inters_size_bounds=(1, np.inf),
    # )
    # fig = fig_dict['figure']

    # outname = get_outname('metrics_upset', name=data_obj.outpath_name, taxon=data_obj.taxon,
    #                       non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
    #                       batch=data_obj.batch_applied,
    #                       batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
    #                       outpath=data_obj.outpath,
    #                       **kws
    # )
    # save_multiple(fig, outname, *file_fmts)
