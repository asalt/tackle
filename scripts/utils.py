"""

"""
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


def filter_observations(panel, column, threshold):
    """
    Filter by less than or equal to threshold of 0 observations
    """

    indices = (panel.minor_xs(column)
               .fillna(0)
               .where(lambda x: x == 0)
               .count(1)
               .where(lambda x: x <= threshold)
               .dropna()
               .index
    )
    return panel.loc[:, indices, :]


def filter_taxon(panel, taxon=9606):
    indices = ((panel.minor_xs('TaxonID') == taxon)
               .any(1)
               .where(lambda x : x == True)
               .dropna()
               .index
    )
    return panel.loc[:, indices, :]


def pearson_r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

def spearman_r2(x, y):
    return stats.spearmanr(x, y)[0] ** 2


rsquare_colors = sb.color_palette("coolwarm", n_colors=100)
def scatter(x, y, stat='pearson', **kwargs):


    kwargs['alpha'] = .2
    ax = plt.gca()
    ax.scatter(x, y, **kwargs)
    if stat == 'pearson':
        rsquared = pearson_r2(x,y)
        text = 'Pearson'
    elif stat == 'spearman':
        rsquared = spearman_r2(x,y)
        text = 'Spearman'
    text += ' r$^2$ {:.2f}'.format(rsquared)
    ax.annotate(text, xy=(0.05, 0.95), xycoords='axes fraction')


    ax_bg_ix = int(round(rsquared, 2) * 100 )
    ax_bg = rsquare_colors[ax_bg_ix]
    ax.patch.set_facecolor(ax_bg)
    # anchored_text = AnchoredText(text, loc=2)
    # ax.add_artist(anchored_text)


def save_multiple(fig, filename, *exts, verbose=True, dpi=300, **save_kwargs):
    """save a figure to a specific file multiple
    times with different extensions"""
    for ext in exts:
        out = filename+ext
        if verbose:
            print("Saving", out, '...', end='', flush=True)
        fig.savefig(out, dpi=dpi, **save_kwargs)
        if verbose:
            print('done.', flush=True)
