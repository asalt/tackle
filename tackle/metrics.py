import os
from collections import OrderedDict
from pathlib import Path
from warnings import warn

from functools import partial
import matplotlib

from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib import markers
from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import seaborn as sb


from .utils import (
    get_outname,
    save_multiple,
    genefilter,
    filter_sra,
    filter_taxon,
    iter_named_items,
)
from . import grdevice_helper


# from .containers import TAXON_MAPPER

idx = pd.IndexSlice

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"


def _build_metrics_from_filtered(df_filtered: pd.DataFrame, *, before_norm: bool):
    if df_filtered is None or df_filtered.empty:
        return OrderedDict()

    df_long = df_filtered.reset_index()
    sample_cols = [c for c in df_long.columns if c not in ("GeneID", "Metric")]
    if not sample_cols:
        return OrderedDict()

    available = set(df_long["Metric"].unique())

    def _metric_frame(metric, numeric=True):
        if metric not in available:
            return None
        frame = df_long.loc[df_long["Metric"] == metric, sample_cols]
        if numeric:
            frame = frame.apply(pd.to_numeric, errors="coerce")
        return frame

    def _sum_series(series):
        return pd.to_numeric(series, errors="coerce").fillna(0).sum()

    sra_frame = _metric_frame("SRA", numeric=False)
    gpg_frame = _metric_frame("GPGroup", numeric=False)
    psms_frame = _metric_frame("PSMs")
    psms_u2g_frame = _metric_frame("PSMs_u2g")
    pept_frame = _metric_frame("PeptideCount")
    pept_u2g_frame = _metric_frame("PeptideCount_u2g")
    pept_s_frame = _metric_frame("PeptideCount_S")
    pept_s_u2g_frame = _metric_frame("PeptideCount_S_u2g")
    area_metric = "iBAQ_dstrAdj" if before_norm else "area"
    area_frame = _metric_frame(area_metric)

    data = OrderedDict()
    for sample in sample_cols:
        entry = dict()
        entry["SRA"] = (
            sra_frame[sample].value_counts(dropna=True).to_dict()
            if sra_frame is not None
            else {}
        )
        entry["GPGroups"] = (
            gpg_frame[sample].dropna().nunique() if gpg_frame is not None else 0
        )
        entry["PSMs"] = {
            "Total": _sum_series(psms_frame[sample]) if psms_frame is not None else 0,
            "u2g": _sum_series(psms_u2g_frame[sample])
            if psms_u2g_frame is not None
            else 0,
        }
        entry["Peptides"] = {
            "Total": _sum_series(pept_frame[sample]) if pept_frame is not None else 0,
            "u2g": _sum_series(pept_u2g_frame[sample])
            if pept_u2g_frame is not None
            else 0,
            "Strict": _sum_series(pept_s_frame[sample])
            if pept_s_frame is not None
            else 0,
            "Strict_u2g": _sum_series(pept_s_u2g_frame[sample])
            if pept_s_u2g_frame is not None
            else 0,
        }
        entry["Area"] = (
            area_frame[sample]
            .where(lambda x: x > 0)
            .dropna()
            .values
            if area_frame is not None
            else np.array([])
        )
        data[sample] = entry

    return data


def _normalize_ratio_frame(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.fillna(0)
    return frame.apply(lambda col: col / col.sum() if col.sum() else col * 0)


def _miscut_long_frame(frame: pd.DataFrame, *, value_name: str) -> pd.DataFrame:
    normalized = _normalize_ratio_frame(frame)
    reset = normalized.reset_index()
    first_col = reset.columns[0]
    if first_col != "miscuts":
        reset = reset.rename(columns={first_col: "miscuts"})
    return reset.melt(id_vars="miscuts", var_name="name", value_name=value_name)


def _append_plot_png(plot_png_paths, plot_root: str) -> None:
    plot_path = Path(f"{plot_root}.png")
    if plot_path.exists():
        plot_png_paths.append(plot_path)


def make_metrics(
    data_obj, file_fmts, png_res=300, before_filter=False, before_norm=False, full=False
):
    from .metrics_html import build_metrics_html_report

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

    if before_filter:
        data = data_obj.metric_values
    else:
        data = _build_metrics_from_filtered(data_obj.df_filtered, before_norm=before_norm)
        if not data:
            warn("Filtered metrics unavailable; falling back to pre-filter metrics.")
            data = data_obj.metric_values

    kws = dict()
    if before_filter:
        kws["filter"] = "before"
        taxon = "all"
    else:
        kws["filter"] = "after"
        taxon = data_obj.taxon
    #
    if before_norm:
        kws["norm"] = "before"
    else:
        kws["norm"] = "after"

    namegen = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method=(
            "parametric" if not data_obj.batch_nonparametric else "nonparametric"
        ),
        normtype=data_obj.normtype,
        outpath=os.path.join(data_obj.outpath, "metrics"),
        # **kws,
    )
    # outname = namegen("metrics")

    sra = pd.DataFrame(data=[data[n]["SRA"] for n in data.keys()], index=data.keys())
    # check if all SRA in sra
    for _k in "S", "R", "A":
        if _k not in sra:
            sra[_k] = 0

    gpg = pd.DataFrame(
        data=[data[n]["GPGroups"] for n in data.keys()], index=data.keys()
    )
    psms = pd.DataFrame(data=[data[n]["PSMs"] for n in data.keys()], index=data.keys())
    peptides = pd.DataFrame(
        data=[data[n]["Peptides"] for n in data.keys()], index=data.keys()
    )
    area = OrderedDict((n, data[n]["Area"]) for n in data.keys())

    # ========================================================================

    miscut_frame = None
    trypsin_ready = bool(data) and all("Trypsin" in data[n] for n in data.keys())
    trypsinp_ready = bool(data) and all("Trypsin/P" in data[n] for n in data.keys())
    if trypsin_ready and trypsinp_ready:
        trypsin = OrderedDict((n, data[n]["Trypsin"]) for n in data.keys())
        trypsin_df = _normalize_ratio_frame(pd.DataFrame(trypsin))

        trypsinP = OrderedDict((n, data[n]["Trypsin/P"]) for n in data.keys())
        trypsinP_df = _normalize_ratio_frame(pd.DataFrame(trypsinP))

        trypsin_dfr = _miscut_long_frame(trypsin_df, value_name="Trypsin")
        trypsinP_dfr = _miscut_long_frame(trypsinP_df, value_name="Trypsin/P")

        miscut_frame = pd.DataFrame.merge(
            trypsin_dfr, trypsinP_dfr, on=["miscuts", "name"], how="outer"
        ).sort_values(by=["miscuts", "name"])

        export_name = namegen("miscut_ratio")
        miscut_frame.to_csv(export_name + ".tsv", sep="\t", index=False)
    else:
        warn("Skipping miscut_ratio: Trypsin metrics unavailable.")

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

    export_name = namegen("metrics")
    to_export.to_csv(export_name + ".tsv", sep="\t", index=True)

    # ==================================================================
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, conversion
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter

    # with localconverter(robjects.default_converter + pandas2ri.converter): # no tuples

    # pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], "R", "metrics.R")
    r_source(r_file)
    Rmetrics = robjects.r["metrics"]

    with localconverter(robjects.default_converter + pandas2ri.converter):  # no tuples
        to_export_r = conversion.py2rpy(to_export)

    metrics_plots = Rmetrics(to_export_r, return_plots=True)

    r_print = robjects.r["print"]
    grdevices = importr("grDevices")
    num_samples = to_export.shape[0]
    plot_width = min(24, max(9, num_samples // 2))
    plot_height = 9
    plot_png_paths = []

    for plot_name, plot_obj in iter_named_items(metrics_plots):
        plot_root = f"{export_name}_{plot_name}"
        for file_fmt in file_fmts:
            out = f"{plot_root}{file_fmt}"
            grdevice = grdevice_helper.get_device(
                filetype=file_fmt,
                width=plot_width,
                height=plot_height,
                res=png_res,
            )
            grdevice(file=out)
            try:
                r_print(plot_obj)
            finally:
                grdevices.dev_off()
        _append_plot_png(plot_png_paths, plot_root)
    # ==================================================================

    # ggridges = importr("ggridges")
    # rboxplot = r["boxplot"]
    # rboxplot

    area_df["Name"] = pd.Categorical(area_df["Name"], ordered=True)

    # with localconverter(robjects.default_converter + pandas2ri.converter): # no tuples
    # this is not working with new rpy2
    # plot = (
    #     gg.ggplot(area_df)
    #     + gg.aes_string(y="Name", x=area_name)
    #     + ggridges.stat_density_ridges(quantile_lines=True, alpha=0.8)
    #     + gg.theme_classic(base_size=12)
    # )
    # plot.plot()

    outname = namegen("metrics_dist", **kws)

    grdevices = importr("grDevices")



    num_samples = len(area_df)
    num_rows = (num_samples // 6) + 1
    num_rows = 1
    # ncol = num_facets # ? needs to be number of samples ber column
    # panel_width = 2.2
    # panel_height = 3
    # # figwidth = max(7, ncol * panel_width)
    # # figheight = panel_height
    figheight = 4 * num_rows
    figwidth = 2 + (0.75 * area_df.Name.nunique()//num_rows)



    print(area_df.head())
    _code = f"""
      library(rlang)
      ycol="{area_name}"
      p <- ggplot2::ggplot(df, aes(x=Name, y=!!sym(ycol))) +
        geom_violin() +
        geom_boxplot(width=.12) +
        labs(y = ycol) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust=0.5,  size = 10)) 

        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        # facet_wrap(~Name, nrow={num_rows}) +
      print(p)
    
    """
    with localconverter(robjects.default_converter + pandas2ri.converter):  # no tuples

        robjects.r.assign("df", area_df)
        robjects.r(_code)

        for file_fmt in file_fmts:
            out = outname + file_fmt
            grdevice = grdevice_helper.get_device(
                filetype=file_fmt,
                width=figwidth,
                height=figheight,
                res=png_res,
            )  # returns a partial grdevice func
            grdevice(file=out)  # open file for saving
            robjects.r(_code)
            # print(ret[0])
            grdevices.dev_off()  # close file
            print(".done", flush=True)
        _append_plot_png(plot_png_paths, outname)

    # area = pd.DataFrame(data=[data[n]['area'] for n in data.keys()], index=data.keys())

    green = "darkgreen"
    yellow = "gold"
    red = "firebrick"

    # outname = get_outname(
    #    "metrics_dist",
    #    name=data_obj.outpath_name,
    #    taxon=data_obj.taxon,
    #    non_zeros=data_obj.non_zeros,
    #    colors_only=data_obj.colors_only,
    #    batch=data_obj.batch_applied,
    #    batch_method="parametric"
    #    if not data_obj.batch_nonparametric
    #    else "nonparametric",
    #    outpath=os.path.join(data_obj.outpath, "metrics"),
    #    **kws,
    # )

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
    # with localconverter(robjects.default_converter + pandas2ri.converter): # no tuples

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
    # outname = get_outname(
    #     "metrics_genecounts",
    #     name=data_obj.outpath_name,
    #     taxon=data_obj.taxon,
    #     non_zeros=data_obj.non_zeros,
    #     colors_only=data_obj.colors_only,
    #     batch=data_obj.batch_applied,
    #     batch_method="parametric"
    #     if not data_obj.batch_nonparametric
    #     else "nonparametric",
    #     outpath=data_obj.outpath,
    #     after="filter",
    #     **kws,
    # )

    outname = namegen("metrics_genecounts")
    save_multiple(fig, outname, *file_fmts)
    _append_plot_png(plot_png_paths, outname)

    metrics_html_path = None
    try:
        metrics_html_path = build_metrics_html_report(
            out_html=export_name + ".html",
            metrics_df=to_export,
            metadata_df=data_obj.col_metadata,
            miscut_df=miscut_frame,
            plot_paths=plot_png_paths,
            export_stem=export_name,
            title=f"Tackle Metrics Report: {data_obj.outpath_name}",
            analysis_label=data_obj.outpath_name,
            before_filter=before_filter,
            before_norm=before_norm,
        )
    except Exception as exc:
        warn(f"Skipping metrics HTML report: {exc}")

    return {
        "metrics_tsv": Path(export_name + ".tsv"),
        "metrics_html": metrics_html_path,
        "miscut_tsv": None
        if miscut_frame is None
        else Path(namegen("miscut_ratio") + ".tsv"),
        "plot_pngs": plot_png_paths,
    }

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
