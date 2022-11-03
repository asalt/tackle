# plot_gsea_results.py
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

import toolz
import pandas as pd


from itertools import tee
from .utils import save_multiple


def pairwise(iterable):
    """
    itertools.pairwise(iterable)¶
    Return successive overlapping pairs taken from the input iterable.

    The number of 2-tuples in the output iterator will be one fewer than the number of inputs. It will be empty if the input iterable has fewer than two values.

    Roughly equivalent to:
    """
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


import click

STAT_METRICS = click.Choice(["FWER p-val", "FDR q-val"])


def filter_results(
    down: pd.DataFrame, up: pd.DataFrame, cutoff=0.25, stat_metric="FWER p-val"
) -> pd.DataFrame:

    # if stat_metric is None or stat_metric is STAT_METRICS:
    # if stat_metric is STAT_METRICS:
    #     stat_metric = "FWER p-val"
    if stat_metric not in STAT_METRICS.choices:
        raise ValueError(f"Must be one of {STAT_METRICS.choices}")

    # this is not great code
    gsea_sig = pd.DataFrame()
    while len(gsea_sig) < 10 or (
        "NES" in gsea_sig
        and len(gsea_sig) < 5
        and (
            gsea_sig.query("NES>0").pipe(len) > 1
            and gsea_sig.query("NES<0").pipe(len) > 1
        )
    ):
        if cutoff > 1:
            break
        gsea_sig = pd.concat(
            [
                down[down[stat_metric] < cutoff],
                up[up[stat_metric] < cutoff],
            ]
        )
        cutoff += 0.1

    if gsea_sig.empty:
        print("No gene sets to plot!")
        return
    number = 999
    cutoff_val = min(abs(gsea_sig.head(number)["NES"].dropna()))
    tokeep1 = (
        gsea_sig.query("NES>0 & abs(NES)>=@cutoff_val")
        .NES.sort_values(ascending=False)
        .head(number)
        .index
    )
    tokeep2 = (
        gsea_sig.query("NES<0 & abs(NES)>=@cutoff_val")
        .NES.sort_values(ascending=False)
        .head(number)
        .index
    )
    idx = [x for y in [tokeep1, tokeep2] for x in y]
    gsea_sig = gsea_sig.loc[idx]
    # gesa_sig = gsea_sig.sort_values(by='NES', ascending=False)

    # gsea_sig["color"] = gsea_sig["FWER p-val"].apply(

    return gsea_sig

    # m = map(lambda x: set(x.parent), pair  )
    # assert len()
    # results = [*pair]
    # assert

    # parse result
    # GSEA outputs the summary files of the form:
    # gsea_report_for_[groupname]_[digit_timestamp].xls
    # group0 = glob.glob(
    #     os.path.join(new_folder, "gsea_report_for_{}_[0-9]*.xls".format(groups[0]))
    # )
    # group1 = glob.glob(
    #     os.path.join(new_folder, "gsea_report_for_{}_[0-9]*.xls".format(groups[1]))
    # )
    # assert len(group0) == len(group1) == 1
    # group0_df = pd.read_table(group0[0], index_col="NAME")
    # group1_df = pd.read_table(group1[0], index_col="NAME")


def plot(
    gsea_sig: pd.DataFrame,
    stat_metric="FWER p-val",
    group0="group0",
    group1="group1",
    outname="gseaout",
    **kws,
):
    if stat_metric not in STAT_METRICS.choices:
        raise ValueError(f"Must be one of {STAT_METRICS}")

    stat_metric_cutoff = 0.25 if stat_metric == "FWER p-val" else 0.05

    cmap = mpl.cm.Reds_r
    bounds = np.linspace(0, 1, 21)
    cnorm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    powernorm = mpl.colors.PowerNorm(0.5, vmin=0, vmax=1)

    cmap = mpl.cm.Reds_r
    bounds = np.linspace(0, 1, 21)
    cnorm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    powernorm = mpl.colors.PowerNorm(0.5, vmin=0, vmax=1)

    gsea_sig["color"] = gsea_sig[stat_metric].apply(
        lambda x: mpl.colors.to_hex(cmap(cnorm(powernorm(x))))
    )
    import textwrap

    gsea_sig.index = [x.replace("HALLMARK_", "") for x in gsea_sig.index]
    gsea_sig.index = [x.replace("REACTOME_", "") for x in gsea_sig.index]
    gsea_sig.index = gsea_sig.index + [
        "*" if x < stat_metric_cutoff else "" for x in gsea_sig[stat_metric]
    ]
    gsea_sig.index = gsea_sig.index.map(
        lambda x: textwrap.fill(x.replace("_", " "), 34, break_long_words=False)
    )

    plt.rc(
        "font",
        **{
            "family": "sans-serif",
            "sans-serif": [
                "DejaVu Sans",
                "Arial",
                "Liberation Sans",
                "Bitstream Vera Sans",
                "sans-serif",
            ],
        },
    )

    nes_range = gsea_sig.NES.abs().max() * 2
    figwidth = 8
    # figwidth = np.round(nes_range * 2.5, decimals=1)
    figheight = max(5, min(gsea_sig.pipe(len) // 1.25, 24))
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = mpl.gridspec.GridSpec(
        2,
        1,
        width_ratios=[
            1,
        ],
        height_ratios=[19, 1],
        hspace=0.4,
    )
    # ax0 = plt.subplot(gs[0])
    ax0 = fig.add_subplot(gs[0])
    # ax1 = plt.subplot(gs[1])
    # cax = plt.subplot(gs[1:])
    cax = fig.add_subplot(gs[1:])
    # should always work now
    try:
        gsea_sig["NES"].fillna(0).plot.barh(
            ax=ax0, color=gsea_sig.color, edgecolor="#222222", linewidth=2
        )
    except Exception as e:
        raise Exception(e)  ## ?? unforseen error
    ax0.axvline(color="#222222")
    # gsea_sig['NES_group1'].fillna(0).plot.barh(colormap=cmap, ax=ax1)
    # ax1.set_yticklabels(())

    gradient = np.apply_along_axis(
        lambda x: cnorm(powernorm(x)), 0, np.linspace(0, 1, 256)
    )
    gradient = np.vstack((gradient, gradient))
    cax.imshow(gradient, aspect="auto", cmap=cmap)
    cax.set_yticklabels(())

    start, end = cax.get_xlim()
    cax.xaxis.set_ticks(np.linspace(start, end, 11))
    cax.set_xticklabels(
        ["{:.2f}".format(x) for x in np.linspace(0, 1, 11)],
        fontdict=dict(size=10),
    )
    # cax.set_xlabel("FWER p-val")
    cax.set_xlabel(stat_metric)
    ax0.grid(axis="x")

    ax0.set_ylabel("")
    ax0.set_xlabel("NES")
    max_nes = gsea_sig.NES.abs().max() + 0.25
    ax0.set_xlim(-max_nes, max_nes)

    for tick in ax0.yaxis.get_ticklabels():
        txt = tick.get_text()
        if len(txt) > 20:
            size = 8
        elif len(txt) > 30:
            size = 6
        else:
            size = 9
        tick.set_size(size)

    # plt.tight_layout()
    # fig = plt.gcf()
    ax0.text(0, 1.04, group0, transform=ax0.transAxes)
    ax0.text(1, 1.04, group1, transform=ax0.transAxes, ha="right")
    ax0.text(
        # -0.04,
        # -0.12,
        -2.2,
        -2.5,
        f"* {stat_metric} < {stat_metric_cutoff}",
        fontsize=12,
        # transform=ax0.transAxes,
    )
    gs.tight_layout(fig, rect=(0, 0, 1, 0.96))
    # fig.subplots_adjust(left=.4)
    # fig.tight_layout()
    outname = str(outname)
    save_multiple(fig, outname, *(".png", ".pdf"))


def main(p: Path):
    print(f"looking at directory {p}")
    g = p.glob("**/*gsea_report*xls")
    from collections import defaultdict


    d = defaultdict(list)
    for entry in g:
        d[entry.parent].append(entry)

    print(d)
    # for pair in pairwise(g):
    for parent, pair in d.items():
        group1 = re.search("(?<=cls_)(.*)(?=_versus)", parent.name)
        group0 = re.search("(?<=versus_).*(?=_pathway)", parent.name)
        group1 = group1.group()
        group0 = group0.group()

        m = map(
            lambda path: re.search("(?<=gsea_report_for_)(.*)(?=_\d+)", path.name), pair
        )
        f = filter(None, m)
        # group0, group1 = [x.group() for x in f]
        assert len(set(x.parent for x in pair)) == 1
        down = [x for x in pair if f"_{group0}" in x.name][0]
        up = [x for x in pair if f"_{group1}" in x.name][0]
        assert down != up
        results = [pd.read_table(x, index_col="NAME") for x in (down, up)]

        # results = [pd.read_table(x, index_col="NAME") for x in pair]
        filtered_results = filter_results(*results)
        parent = pair[0].parent

        outname = parent.name
        # if any(x=='0um72h' for x in (group0, group1)):
        #     if any(x=='250um72h' for x in (group0, group1)):
        #         import ipdb; ipdb.set_trace()
        #    print(parent, '\n', pair, '\n', group0, group1)
        plot(filtered_results, outname=outname, group0=group0, group1=group1)
        # xx = [x.name for x in pair]

        # import ipdb

        # ipdb.set_trace()

        # outname = parent.joinpath(parent.name)
        # outname = Path(".")
