"""

"""
import sys
import os
import re
import operator as op
from collections import OrderedDict
from functools import partial

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sb
from seaborn.distributions import _freedman_diaconis_bins
import click


from bcmproteomics_ext import ispec
sb.set_context('notebook', font_scale=1.4)

from utils import *

sys.setrecursionlimit(10000)

TAXON_MAPPER = {'human': 9606,
                'mouse': 10090}

def run(records, ZEROS=0, stat='spearman', taxon='all', data_dir=None, OUTPATH='../figures',
        funcat=None, geneid_subset=None, highlight_gids=None, highlight_gid_names=None):

    if stat not in ('pearson', 'spearman'):
        raise ValueError('Must select from `pearson` or `spearman`')

    exps = OrderedDict()
    for name, record in records.items():
        exp = ispec.E2G(data_dir=data_dir, **record)
        if len(exp) == 0:
            print('No data in {!r}'.format(exp))
            continue
        if funcat:
            df = exp.df[ exp.df['FunCats'].fillna('').str.contains(funcat)]
        else:
            df = exp.df
        if geneid_subset:
            df = df.loc[geneid_subset]
        exps[name] = df

    if ZEROS == 'max':
        ZEROS = len(exps)

    panel = pd.Panel(exps)

    dummy_filter = lambda x, *args, **kwargs: x
    taxon_filter = TAXON_MAPPER.get(taxon)
    if taxon_filter is None:
        filter_func = dummy_filter
    else:
        filter_func = partial(filter_taxon, taxon=taxon_filter)

    panel_filtered = (panel.pipe(filter_observations, 'iBAQ_dstrAdj', ZEROS)
                      .pipe(filter_func)
    )

    ibaqs = panel_filtered.minor_xs('iBAQ_dstrAdj').astype(float)

    ibaqs_log = np.log10(ibaqs.fillna(0)+1e-10)

    minval = ibaqs_log.min().min()
    shift_val = np.ceil(np.abs(minval))

    ibaqs_log_shifted = ibaqs_log + shift_val

    # xymin = 0
    xymin = np.floor(ibaqs_log_shifted[ibaqs_log_shifted > 0].min().min())
    xymax = np.ceil(ibaqs_log_shifted.max().max())

    g = sb.PairGrid(ibaqs_log_shifted)
    g.map_upper(plot_delegator, stat=stat, filter_zeros=True, upper_or_lower='upper')
    g.map_lower(plot_delegator, stat=stat, filter_zeros=True, upper_or_lower='lower',
                xymin=xymin, xymax=xymax)
    g.map_diag(hist, xmin=xymin, xmax=xymax)

    sb.despine(fig=g.fig, left=True, bottom=True)
    remove_ticklabels(fig=g.fig)

    #  adjust the spacing between subplots
    hspace = g.fig.subplotpars.hspace
    wspace = g.fig.subplotpars.wspace
    g.fig.subplots_adjust(hspace=hspace*.1, wspace=wspace*.1, right=.8, bottom=.2,
                          left=.2, top=.8)

    outpath_name = os.path.split(OUTPATH)[-1]
    # g.fig.suptitle(outpath_name.replace('_', ' '))
    # cbar_ax = g.fig.add_axes([.85, .15, .05, .7])
    cbar_ax = g.fig.add_axes([.85, .15, .05, .65])
    plot_cbar(cbar_ax)

    # range_ax = g.fig.add_axes([.25, .05, .65, .05])
    range_ax = g.fig.add_axes([.20, .06, .50, .06])
    range_ax.set_xlim((xymin, xymax))
    props = dict(color='black', linewidth=2, markeredgewidth=2)
    with sb.plotting_context(context='notebook', font_scale=1.4):
        font_size = mpl.rcParams['font.size'] * .9
        make_xaxis(range_ax, fmt_str='%1.0f', **props)
        range_ax.set_xlabel('log10 iBAQ dstrAdj', labelpad=16,
                            fontdict={'size': font_size, 'weight': 'normal'},)
        sb.despine(ax=range_ax, left=True, bottom=True)
        remove_ticklabels(ax=range_ax)


    outname = os.path.join(OUTPATH,
                           '{}_scatter_{}_{}less_zeros'.format(outpath_name, taxon, ZEROS))
    save_multiple(g, outname, '.png', dpi=96)

    ibaqs_zscore = (ibaqs_log_shifted - ibaqs_log_shifted.mean()) / ibaqs_log_shifted.std()


    row_colors = None
    if highlight_gids:

        cmap = sb.color_palette(n_colors=max(8, len(highlight_gids)))
        colors_dfs = list()
        for ix, hgid in enumerate(highlight_gids):
            color = cmap[ix]
            highlights = {gid: color for gid in hgid}
            colors = ibaqs_zscore.index.map( lambda x: highlights.get(x, 'white') )
            colors_df = pd.DataFrame(colors, index=ibaqs_zscore.index)
            colors_dfs.append(colors_df)
        row_colors = pd.concat(colors_dfs,axis=1)
        row_colors.columns = highlight_gid_names

    g = sb.clustermap(ibaqs_zscore, row_colors=row_colors, yticklabels=False)
    # g.ax_heatmap.set_yticklabels([], rotation=0)
    outname = os.path.join(OUTPATH,
                           '{}_clustermap_{}_{}less_zeros'.format(outpath_name, taxon, ZEROS))
    save_multiple(g, outname, '.png',)

@click.command()
@click.option('--data-dir', type=click.Path(exists=True, file_okay=False),
              default='../data/raw',
              help='optional location to store and read e2g files')
@click.option('--geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True,
              help="""Optional list of geneids to subset by.
              Should have 1 geneid per line. """)
@click.option('--highlight-geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=True,
              help="""Optional list of geneids to highlight by.
              Should have 1 geneid per line. """)
@click.option('--funcat', type=str, default=None, show_default=True,
              help="""Optional gene subset based on funcat or funcats,
              regular expression allowed. """)
@click.option('--stat', type=click.Choice(['pearson', 'spearman']),
              default='spearman', show_default=True)
@click.option('--taxon', type=click.Choice(['human', 'mouse', 'all']),
              default='all', show_default=True)
@click.option('-z', '--zeros', default=0, show_default=True,
              help='Number of zeros tolerated across all samples.')
@click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
def main(data_dir, geneids, highlight_geneids, funcat, stat, taxon, zeros, experiment_file):
    fname, ext = os.path.splitext(experiment_file)

    geneid_subset=None
    if geneids:
        geneid_subset = parse_gid_file(geneids)
        if len(geneid_subset) == 0:
            print('Non geneids found in file {}'.format(geneids))

    highlight_gids = None
    highlight_gid_names = None
    if highlight_geneids:
        highlight_gids = list()
        highlight_gid_names = list()

        for ix, h_gid in enumerate(highlight_geneids):
            highlight_gid = parse_gid_file(h_gid)
            highlight_gids.append(highlight_gid)
            h_gid_name = get_file_name(h_gid)
            if h_gid_name:
                highlight_gid_names.append(h_gid_name)
            else:
                highlight_gid_names.append(ix)

        if len(highlight_gids) == 0:
            print('Non geneids found in file {}'.format(highlight_geneids))

    try:
        analysis_name = re.findall('\w+', fname)[-1]
    except IndexError:
        print('Error parsing configfile name.')
        analysis_name = 'Unnamed'

    OUTPATH = os.path.join('../figures/', analysis_name)
    if not os.path.exists(OUTPATH):
        os.mkdir(OUTPATH)

    data = read_config(experiment_file)
    if len(data) == 0:
        raise ValueError('No items in configfile.')


    # experiment_file
    run(data, ZEROS=zeros, stat=stat, taxon=taxon, data_dir=data_dir, OUTPATH=OUTPATH,
        funcat=funcat, geneid_subset=geneid_subset, highlight_gids=highlight_gids,
        highlight_gid_names=highlight_gid_names
    )


if __name__ == '__main__':

    main()


    # EXPS = [
    #     dict(
    #         recno=32918,
    #         runno=1,
    #         searchno=1),
    #     dict(
    #         recno=32919,
    #         runno=1,
    #         searchno=1
    #     ),

    # ]
    # ZEROS = 0 # a number or 'max'
    # stat='spearman'

    # main(EXPS, ZEROS, stat)

    # experiment_file = './configtest.ini'
    # fname, ext = os.path.splitext(experiment_file)
    # geneids = './geneid_subset.txt'
    # highlight_geneids = './highlight_subset.txt'
