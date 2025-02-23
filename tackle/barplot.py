from collections import OrderedDict

from  warnings import warn
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib import gridspec
import seaborn as sb

from .utils import save_multiple

def plotter(data, linear=False, z_score=False, colors=None, legend_colors=None, legend_name=None, title=None):

    sb.set_context('notebook')

    w = len(data.columns) * .75
    w = max(w, 4)
    ylabel = 'log$_{10}$ iBAQ' if not linear else 'iBAQ'
    if z_score:
        ylabel += ' z scored'

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
    # fig.autofmt_xdate()
    # plt.setp( ax[1].xaxis.get_majorticklabels(), rotation=90 )

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

def barplot(X, genes, metadata, average=None, color=None, color_order=None, cmap=None, linear=False,
            figsize=None, xtickrotation=None, xticksize=None,
            z_score=False, base_outfunc=None, file_fmts=('.png',), gid_symbol=None, metadata_colors=None,
            retain_order=False
):
    """
    Plot barplots for each gene in a list of genes
    :file_fmts: iterable of valid file formats for matplotlib to save figures
    """

    orig_rc = mpl.rcParams
    # sb.set_context('notebook', font_scale=.8)


    rc = {'font.size': 12.0, 'axes.labelsize': 12.0, 'axes.titlesize': 12.0, 'xtick.labelsize': 11.0,
          'ytick.labelsize': 11.0, 'legend.fontsize': 11.0, 'axes.linewidth': 1.25, 'grid.linewidth': 1.0,
          'lines.linewidth': 1.5, 'lines.markersize': 6.0, 'patch.linewidth': 1.0, 'xtick.major.width': 1.25,
          'ytick.major.width': 1.25, 'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0, 'xtick.major.size':
          6.0, 'ytick.major.size': 6.0, 'xtick.minor.size': 4.0, 'ytick.minor.size': 4.0}
    sb.set_context('notebook', rc=rc)
    mpl.rcParams.update(rc)
    # why does none of this work??

    if xtickrotation is None:
        xtickrotation=30


    if cmap is None:
        cmap = 'tab10'
    if gid_symbol is None:
        gid_symbol = dict()

    if color:
        color_groups = metadata.groupby(color).groups
        cpalette = iter(sb.color_palette(cmap, n_colors=len(color_groups)))
        colors = dict()
        legend_colors = OrderedDict()
        if average:
            color_entries = {x:k for k,v in metadata.groupby([average]).groups.items() for x in v}
        else:
            color_entries = None
        for g, cols in color_groups.items():
            if metadata_colors and g in metadata_colors:
               thecolor = metadata_colors[g]
            elif g == 'NA':
                thecolor = 'grey'
            else:
                thecolor = to_hex(next(cpalette))
            legend_colors[g] = thecolor
            for c in cols:
                if color_entries:
                    key = color_entries.get(c)
                    if key not in colors:
                        colors[key] = thecolor
                else:
                    colors[c] = thecolor
        # if 'NA' in colors:
        #     color['NA'] = 'grey'
        # if 'NA' in legend_colors:
        #     legend_colors['NA'] = 'grey'

        # the_order = [x for y in color_groups.values() for x in y]
        # the_colors = [colors[x] for x in the_order]
        # the_order = [x for y in legend_colors.values() for x in y]
        # the_colors = [colors[x] for x in the_order]

        the_order = list()

        if color_order is not None:
            _color_order = color_order.split('|')
            _missing = set(legend_colors.keys()) - set(_color_order)
            for m in _missing:
                _color_order.append(m)
            _iter = ((x, legend_colors.get(x)) for x in _color_order)

        else:
            _iter = legend_colors.items()
            the_order = list()

        if not retain_order:
            for entry, _color in _iter:
                samples = [x for x, y in colors.items() if y == _color]
                for s in samples:
                    the_order.append(s)


    else:
        colors = None
        legend_colors = None
        the_order, the_colors = None, None

        # colors is a dict with sample_name : hex_color

    if linear:
        X = X.apply(lambda x: np.power(10, x))
    elif z_score:
        X = sb.matrix.ClusterGrid.z_score(X, 0)

    if average is not None:
        groupcols = metadata.groupby(average).groups # dict with groupname : [cols]
    else:
        groupcols = OrderedDict([(x, x) for x in X.columns])

    for gene in genes:

        if gene not in X.index:
            warn('GeneID {} not in dataset, skipping..'.format(gene))
            continue
        edata = X.loc[gene].to_frame('Expression').join(metadata).reset_index()
        if the_order:
            if not average:
                new_ix_order = list()
                for x in the_order:
                    query = edata[edata['index'] == x]
                    if not query.empty:
                        new_ix_order.append(query.index[0])
                # new_ix_order = [edata[edata['index'] == x].index[0] for x in the_order]
            elif average: # then need to extract each individual entry
                new_ix_order = list()
                for x in the_order:
                    for y in groupcols[x]:
                        new_ix_order.append(
                            edata[edata['index'] == y].index[0]
                        )

            edata = edata.loc[new_ix_order]
        # edata = X.loc[gene]

        # data = OrderedDict()
        # for k, cols in groupcols.items():
        #     mean_ = edata.loc[cols].mean()
        #     std_ = edata.loc[cols].std()
        #     data[k] = {'mean': mean_, 'std': std_}

        # df = pd.DataFrame(data)

        symbol = gid_symbol.get(gene, '')
        title = "{} {}".format(gene, symbol)


        if not average:
            data_size = len(edata)
        else:
            data_size = metadata[average].nunique()

        MAX_ROW = 18

        nrow = int( np.ceil(data_size / MAX_ROW) )
        # w = data_size * .75
        w = min(data_size, MAX_ROW) * .75
        w = max(w, 4)

        # fig, ax = plt.subplots(figsize=(w,5))
        if figsize:
            w, h = figsize
        else:
            h = 5 * (nrow*.75)
            h = min(h, 15)
        print(w,h)

        # plt.rcParams['figure.constrained_layout.use'] = True
        fig = plt.figure(figsize=(w,h))
        # gs = gridspec.GridSpec(1,10)
        gs = gridspec.GridSpec(nrow, 10, wspace=.4, figure=fig)

        # ax = fig.add_subplot(gs[0, 0:9])


        x = 'index' if not average else average
        chunks = np.array_split(np.arange(data_size), nrow)

        axs = list()
        for chunk, row in zip(chunks, range(nrow)):
            if axs:
                ax = fig.add_subplot(gs[row, 0:9], sharey=axs[row-1])
            else:
                ax = fig.add_subplot(gs[row, 0:9])

            if not average:
                sb.barplot(x=x, y='Expression', data=edata.iloc[chunk], ax=ax, ci='sd', capsize=.2,
                        # palette=colors, order=the_order,
                        palette=colors, order=None,
                )

            if average:

                # if the_order is not None:
                if the_order:
                    group_chunk = the_order[ chunk[0]:chunk[-1]+1 ]
                else:
                    group_chunk = edata[average].unique()[chunk]

                if metadata_colors and average in metadata_colors:
                    color_mapping = metadata_colors[average]
                    # group_chunk is the collection of items to plot, a subset of the
                    # "average" indicator category, in the correct order
                    colors = [color_mapping.get(x) for x in group_chunk]

                sb.barplot(x=x, y='Expression', data=edata[edata[average].isin(group_chunk)],
                           ax=ax, ci='sd', capsize=.2,
                        # palette=colors, order=the_order,
                        palette=colors, order=None,
                )

                sb.swarmplot(x=x, y='Expression', hue=None, data=edata[edata[average].isin(group_chunk)],
                             order=None,
                             # order=the_order,
                             hue_order=None, dodge=False, orient=None, color=None, palette=colors,
                             size=5, edgecolor='gray', linewidth=1., ax=ax, )


            yticks = ax.get_yticks()
            print(yticks)
            if len(yticks) < 3:
                ax.yaxis.set_major_locator(mpl.ticker.LinearLocator(3))
                ax.yaxis.set_minor_locator(MultipleLocator(1))


            ax.grid(axis='y', which='both')
            axs.append(ax)


            ylabel = 'log$_{10}$ iBAQ' if not linear else 'iBAQ'
            if z_score:
                ylabel += ' z scored'
            ax.set_ylabel(ylabel)
            ax.set_xlabel('')

            ax.legend([], frameon=False)
            # if xtickrotation is None:
            #     xtickrotation=45
            #     # fig.autofmt_xdate()
            # else:
            for tick in ax.get_xticklabels():
                tick.set_rotation(xtickrotation)
                if xticksize is not None:
                    tick.set_fontsize(xticksize)

        axs[0].set_title(title)
        ax_leg = fig.add_subplot(gs[:, 9])


        if legend_colors is not None:
            handles, labels = list(), list()
            if color_order:
                _iter = (( x, legend_colors.get(x) ) for x in _color_order if x in legend_colors )
            else:
                _iter = legend_colors.items()
            for key, value in _iter:
                handles.append(mpl.patches.Patch(color=value))
                labels.append(key)
            if len(labels) <= 20:
                ncols = max(len(labels) // 4, 1)
            else:
                ncols = max(len(labels) // 8, 1)

            leg = ax_leg.legend( handles, labels,
                                 loc='center', ncol=ncols,
                                 title=color, frameon=False,
                                 handleheight=.8, handlelength=.8,
            )
            ax_leg.add_artist(leg)

        # ax.set_ylabel(ylabel)
        sb.despine(ax=ax_leg, bottom=True, left=True)
        ax_leg.set_xticks([])
        ax_leg.set_yticks([])

        # if len(chunks) == 1:
        #     fig.subplots_adjust(left=.10, bottom=.20, top=.90, hspace=.8)
        # else:
        #     fig.subplots_adjust(left=.10, bottom=.08, top=.96, hspace=.8)

        # do this at the end, after all plots are drawn
        for ax in axs:
            yticks = ax.get_yticks()
            print(yticks)
            if len(yticks) < 4:

                ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(
                    np.ceil(ax.get_ylim()[1] / 3)
                )
                )

                # ax.yaxis.set_major_locator(mpl.ticker.LinearLocator(3))
                # ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))

        # ax.set_xlabel('')

        # ax_leg.legend

        # fig, ax = plotter(df, colors=colors, linear=linear, z_score=z_score, legend_name=color, legend_colors=legend_colors,
        #                   title=title)

        if average:
            name = '{}_{}'.format(str(gene), average)
        elif not average:
            name = str(gene)
        if linear:
            name += '_linear'
        if z_score:
            name += '_zscore'
        outname = base_outfunc(name)
        save_multiple(fig, outname, *file_fmts, dpi=300)
        plt.close(fig)

    mpl.rcParams.update(orig_rc)
