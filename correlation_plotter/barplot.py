from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib import gridspec
import seaborn as sb

from .utils import save_multiple

def plotter(data, linear=False, colors=None, legend_colors=None, legend_name=None, title=None):

    w = len(data.columns) * .75
    w = max(w, 4)
    ylabel = 'log$_{10}$ iBAQ' if not linear else 'iBAQ'

    kwargs = dict()
    if colors:
        barcolors = [colors[x] for x in data.columns]
        kwargs['color'] = barcolors

    # fig, ax = plt.subplots(figsize=(w,5))
    fig = plt.figure(figsize=(w,5))
    gs = gridspec.GridSpec(1,10)
    ax = fig.add_subplot(gs[0, 0:9])
    ax_leg = fig.add_subplot(gs[0, 9])

    data.loc['mean'].plot.bar(yerr=data.loc['std'], rot=0, ax=ax, **kwargs)
    fig.autofmt_xdate()
    ax.set_ylabel(ylabel)
    sb.despine(ax=ax_leg, bottom=True, left=True)
    ax_leg.set_xticks([])
    ax_leg.set_yticks([])

    if legend_colors is not None:
        handles, labels = list(), list()
        for key, value in legend_colors.items():
            handles.append(mpl.patches.Patch(color=value))
            labels.append(key)
        if len(labels) <= 20:
            ncols = max(len(labels) // 4, 1)
        else:
            ncols = max(len(labels) // 8, 1)

        leg = ax_leg.legend( handles, labels,
                             loc='center', ncol=ncols,
                             title=legend_name
        )
        ax_leg.add_artist(leg)

    if title:
        ax.set_title(title)

    fig.tight_layout(w_pad=.1)
    return fig, ax

def barplot(X, genes, metadata, average=None, color=None, cmap=None, linear=False,
            base_outfunc=None, file_fmts=('.png',), gid_symbol=None):
    """
    Plot barplots for each gene in a list of genes
    :file_fmts: iterable of valid file formats for matplotlib to save figures
    """
    if cmap is None:
        cmap = 'tab10'
    if gid_symbol is None:
        gid_symbol = dict()

    if color:
        color_groups = metadata.T.groupby(color).groups
        cpalette = iter(sb.color_palette(cmap, n_colors=len(color_groups)))
        colors = dict()
        legend_colors = OrderedDict()
        for g, cols in color_groups.items():
            thecolor = to_hex(next(cpalette))
            legend_colors[g] = thecolor
            for c in cols:
                colors[c] = thecolor
    else:
        colors = None
        legend_colors = None

        # colors is a dict with sample_name : hex_color



    if linear:
        X = X.apply(lambda x: np.power(10, x))

    if average is not None:
        groupcols = metadata.T.groupby(average).groups # dict with groupname : [cols]
    else:
        groupcols = {x: x for x in X.columns}
    for gene in genes:


        edata = X.loc[gene]
        data = OrderedDict()
        for k, cols in groupcols.items():
            mean_ = edata.loc[cols].mean()
            std_ = edata.loc[cols].std()
            data[k] = {'mean': mean_, 'std': std_}

        df = pd.DataFrame(data)

        symbol = gid_symbol.get(gene, '')
        title = "{} {}".format(gene, symbol)
        fig, ax = plotter(df, colors=colors, linear=linear, legend_name=color, legend_colors=legend_colors,
                          title=title)

        if average:
            name = '{}_{}'.format(str(gene), average)
        elif not average:
            name = str(gene)
        if linear:
            name += '_linear'
        outname = base_outfunc(name)
        save_multiple(fig, outname, *file_fmts, dpi=300)
