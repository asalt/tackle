import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import seaborn as sb

from .utils import *

rc = {"font.family": "serif", "font.serif": ["Times", "Palatino", "serif"]}
sb.set_context("paper", font_scale=0.8)
sb.set_style("white", rc)


def scatterplot(
    ibaqs_log_shifted,
    mask=None,
    stat="pearson",
    colors_only=False,
    shade_correlation=True,
    outname="name",
    file_fmts=None,
):
    try:
        from rpy2.robjects import r
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr

        _viaR = True
    except ModuleNotFoundError:
        _viaR = False

    cor_mat = ibaqs_log_shifted.T.corr(method=stat)
    f = outname + "_corr.tsv"
    cor_mat.to_csv(f, sep="\t")

    if _viaR:
        pandas2ri.activate()
        r_source = robjects.r["source"]
        r_file = os.path.join(
            os.path.split(os.path.abspath(__file__))[0], "R", "scatter.R"
        )
        r_source(r_file)
        grdevices = importr("grDevices")
        Rscatterplot = robjects.r["scattermat"]

        plt_size = ibaqs_log_shifted.shape[1] * 0.75
        if file_fmts is None:
            file_fmts = (".png",)
        # file_fmts = (".png",)
        gr_devices = {
            ".png": grdevices.png,
            ".pdf": grdevices.pdf,
            ".svg": grdevices.svg,
        }
        gr_kws = {
            ".png": dict(width=plt_size, height=plt_size, units="in", res=300),
            ".pdf": dict(
                width=plt_size,
                height=plt_size,
            ),
            ".svg": dict(
                width=plt_size,
                height=plt_size,
            ),
        }
        for file_fmt in file_fmts:
            out = outname + file_fmt
            grdevice = gr_devices[file_fmt]
            gr_kw = gr_kws[file_fmt]
            print("Saving", out, "...", end="", flush=True)
            grdevice(file=out, **gr_kw)
            Rscatterplot(
                ibaqs_log_shifted,
                method=stat,
                interactive=False,
                colors_only=colors_only,
            )
            grdevices.dev_off()
            print("done.", flush=True)

        return

    # this is all old
    # xymin = 0
    xymin = np.floor(ibaqs_log_shifted[ibaqs_log_shifted > 0].min().min())
    xymax = np.ceil(ibaqs_log_shifted.max().max())

    g = sb.PairGrid(ibaqs_log_shifted)
    g.map_upper(
        plot_delegator,
        stat=stat,
        filter_zeros=True,
        upper_or_lower="upper",
        colors_only=colors_only,
        shade_correlation=shade_correlation,
    )
    g.map_lower(
        plot_delegator,
        stat=stat,
        filter_zeros=True,
        upper_or_lower="lower",
        xymin=xymin,
        xymax=xymax,
        colors_only=colors_only,
        shade_correlation=shade_correlation,
    )
    g.map_diag(hist, xmin=xymin, xmax=xymax, colors_only=colors_only)
    if shade_correlation:
        color_diag(g)
    sb.despine(fig=g.fig, left=True, bottom=True)
    remove_ticklabels(fig=g.fig)
    #  adjust the spacing between subplots
    hspace = g.fig.subplotpars.hspace
    wspace = g.fig.subplotpars.wspace
    g.fig.subplots_adjust(
        hspace=hspace * 0.1,
        wspace=wspace * 0.1,
        right=0.9,
        bottom=0.1,
        left=0.1,
        top=0.9,
    )
    # g.fig.suptitle(outpath_name.replace('_', ' '))
    # cbar_ax = g.fig.add_axes([.85, .15, .05, .7])
    # cbar_ax = g.fig.add_axes([.95, .15, .05, .85])
    # plot_cbar(cbar_ax)
    range_ax = g.fig.add_axes([0.25, 0.05, 0.65, 0.05])
    # range_ax = g.fig.add_axes([.10, .06, .50, .06])
    range_ax.set_xlim((xymin, xymax))
    props = dict(color="black", linewidth=2, markeredgewidth=2)
    with sb.plotting_context(context="notebook", font_scale=1.4):
        font_size = mpl.rcParams["font.size"] * 0.9
        make_xaxis(range_ax, fmt_str="%1.0f", **props)
        range_ax.set_xlabel(
            "log10 iBAQ dstrAdj",
            labelpad=16,
            fontdict={"size": font_size, "weight": "normal"},
        )
        sb.despine(ax=range_ax, left=True, bottom=True)
        remove_ticklabels(ax=range_ax)

    return g
