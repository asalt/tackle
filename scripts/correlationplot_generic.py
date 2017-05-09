"""

"""
import sys
import os
import re
import operator as op

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
sb.set_context('notebook', font_scale=1.2)

from utils import *

sys.setrecursionlimit(10000)

def run(records, ZEROS=0, stat='spearman', data_dir=None):

    if stat not in ('pearson', 'spearman'):
        raise ValueError('Must select from `pearson` or `spearman`')



    exps = dict()
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

    ibaqs = panel_filtered.minor_xs('iBAQ_dstrAdj')
    min_nonzero = ibaqs.where(lambda x: x > 0).min().min()
    ibaqs = ((ibaqs.fillna(0) + min_nonzero*.1)
            .apply(np.log10)
    )




    g = sb.PairGrid(ibaqs)
    g.map_offdiag(scatter, stat='spearman')
    g.map_diag(plt.hist)
    # save_multiple(g, '../figures/correlationplot2/scatter_human_3less_zeros', '.png', '.pdf',
    #               dpi=96)
    save_multiple(g, '../figures/scatter_human_{}less_zeros'.format(ZEROS), '.png',
                dpi=96)


    ibaqs_zscore = (ibaqs - ibaqs.mean()) / ibaqs.std()
    g = sb.clustermap(ibaqs_zscore)
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
