"""

"""
import sys
import os
import re
import operator as op
from collections import OrderedDict

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

def run(records, ZEROS=0, stat='spearman', data_dir=None):

    if stat not in ('pearson', 'spearman'):
        raise ValueError('Must select from `pearson` or `spearman`')



    exps = OrderedDict()
    for name, record in records.items():
        exp = ispec.E2G(**record, data_dir=data_dir)
        if len(exp) == 0:
            print('No data in {!r}'.format(exp))
            continue
        exps[name] = exp.df

    if ZEROS == 'max':
        ZEROS = len(exps)

    panel = pd.Panel(exps)
    panel_filtered = (panel.pipe(filter_observations, 'iBAQ_dstrAdj', ZEROS)
                      .pipe(filter_taxon)
    )

    ibaqs = panel_filtered.minor_xs('iBAQ_dstrAdj').astype(float)

    ibaqs_log = ibaqs.apply(np.log10)

    minval = ibaqs_log.min().min()
    shift_val = np.ceil(np.abs(minval))

    ibaqs_log_shifted = ibaqs_log + shift_val

    xymin = np.floor(ibaqs_log_shifted.min().min())
    xymax = np.ceil(ibaqs_log_shifted.max().max())

    # min_nonzero = ibaqs.where(lambda x: x > 0).min().min()
    # ibaqs = ((ibaqs.fillna(0) + min_nonzero*.1)
    #         .apply(np.log10)
    # )


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
    g.fig.subplots_adjust(hspace=hspace*.25, wspace=wspace*.25, right=.8, bottom=.2)

    cbar_ax = g.fig.add_axes([.85, .15, .05, .7])
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


    # save_multiple(g, '../figures/correlationplot2/scatter_human_3less_zeros', '.png', '.pdf',
    #               dpi=96)
    save_multiple(g, '../figures/scatter_human_{}less_zeros'.format(ZEROS), '.png',
                  dpi=96)


    min_nonzero = ibaqs.where(lambda x: x > 0).min().min()
    ibaqs_nonzero = (ibaqs.fillna(0) + min_nonzero*.1).apply(np.log10)
    ibaqs_zscore = (ibaqs_nonzero - ibaqs_nonzero.mean()) / ibaqs_nonzero.std()
    g = sb.clustermap(ibaqs_zscore)
    g.ax_heatmap.set_yticklabels([])
    save_multiple(g, '../figures/clustermap_human_{}less_zeros'.format(ZEROS), '.png',)

@click.command()
@click.option('--data-dir', type=click.Path(exists=True, file_okay=False),
              default='../data/raw',
              help='optional location to store and read e2g files')
@click.option('--stat', type=click.Choice(['pearson', 'spearman']),
              default='spearman', show_default=True)
@click.option('-z', '--zeros', default=0, show_default=True,
              help='Number of zeros tolerated across all samples.')
@click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
def main(data_dir, stat, zeros, experiment_file):
    data = read_config(experiment_file)
    if len(data) == 0:
        raise ValueError('No items in configfile.')


    # experiment_file
    run(data, ZEROS=zeros, stat=stat, data_dir=data_dir)


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
