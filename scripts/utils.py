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

def hist(x, **kwargs):
    X = x[ (~x.isnull()) & ~(x.abs() == np.inf) ]
    plt.hist(X.values, **kwargs)

rsquare_colors = sb.color_palette("coolwarm", n_colors=100)

def plot_delegator(x, y, stat='pearson', filter_zeros=True,
                    upper_or_lower='upper', **kwargs):
    if upper_or_lower == 'upper':
        func = annotate_stat
    elif upper_or_lower == 'lower':
        func = scatter

    x_nonzero = x[ (~x.isnull()) & ~(x.abs() == np.inf) ].index
    y_nonzero = y[ (~y.isnull()) & ~(y.abs() == np.inf) ].index

    nonzeros = list(set(x_nonzero) & set(y_nonzero))

    X, Y = x, y

    if filter_zeros:
        X = x.loc[nonzeros]
        Y = y.loc[nonzeros]

    kwargs['alpha'] = .2
    ax = plt.gca()

    if stat == 'pearson':
        rsquared = pearson_r2(X,Y)
        text = 'Pearson'
    elif stat == 'spearman':
        rsquared = spearman_r2(X,Y)
        text = 'Spearman'
    text = 'n = {:,}\nr$^2$ = {:.2f}'.format(len(X), rsquared)


    ax_bg_ix = int(round(rsquared, 2) * 100 )
    ax_bg = rsquare_colors[ax_bg_ix]
    ax.patch.set_facecolor(ax_bg)
    # kwargs['text'] = text

    func(X, Y, ax, text=text, **kwargs)


def annotate_stat(x, y, ax, text, **kwargs):

    ax.annotate(text, xy=(0.5, 0.5), xycoords='axes fraction',
                va='center', ha='center'
    )

def scatter(x, y, ax, **kwargs):
    try:
        kwargs.pop('text')
    except KeyError:
        pass

    ax.scatter(x, y, c='k', **kwargs)

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
