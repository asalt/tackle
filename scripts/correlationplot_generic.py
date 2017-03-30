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


from bcmproteomics_ext import ispec
sb.set_context('notebook', font_scale=1.2)

from utils import *

sys.setrecursionlimit(10000)

def main(records, ZEROS=0, stat='spearman'):

    if stat not in ('pearson', 'spearman'):
        raise ValueError('Must select from `pearson` or `spearman`')



    exps = dict()
    for record in records:
        exp = ispec.E2G(data_dir='../data/raw/', **record)
        if len(exp) == 0:
            print('No data in {!r}'.format(exp))
            continue
        exps[repr(exp)] = exp.df

    if ZEROS == 'max':
        ZEROS = len(exps)

    panel = pd.Panel(exps)
    panel_filtered = (panel.pipe(filter_observations, 'iBAQ_dstrAdj', ZEROS)
                      .pipe(filter_taxon)
    )

    ibaqs = panel_filtered.minor_xs('iBAQ_dstrAdj').astype(float)

    ibaqs_log = ibaqs.apply(np.log10)
    # min_nonzero = ibaqs.where(lambda x: x > 0).min().min()
    # ibaqs = ((ibaqs.fillna(0) + min_nonzero*.1)
    #         .apply(np.log10)
    # )


    g = sb.PairGrid(ibaqs_log)
    g.map_offdiag(scatter, stat='spearman', filter_zeros=True)
    g.map_diag(hist)
    # save_multiple(g, '../figures/correlationplot2/scatter_human_3less_zeros', '.png', '.pdf',
    #               dpi=96)
    save_multiple(g, '../figures/scatter_human_{}less_zeros'.format(ZEROS), '.png',
                  dpi=96)


    min_nonzero = ibaqs.where(lambda x: x > 0).min().min()
    ibaqs_nonzero = (ibaqs.fillna(0) + min_nonzero*.1).apply(np.log10)
    ibaqs_zscore = (ibaqs_nonzero - ibaqs_nonzero.mean()) / ibaqs_nonzero.std()
    g = sb.clustermap(ibaqs_zscore)
    g.ax_heatmap.setyticklabels([])
    save_multiple(g, '../figures/clustermap_human_{}less_zeros'.format(ZEROS), '.png',)

if __name__ == '__main__':

    EXPS = [
        dict(
            recno=32918,
            runno=1,
            searchno=1),
        dict(
            recno=32919,
            runno=1,
            searchno=1
        ),

    ]
    ZEROS = 'max'  # a number or 'max', for the clustermap
    stat='spearman'

    main(EXPS, ZEROS, stat)
