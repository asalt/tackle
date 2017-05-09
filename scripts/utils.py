"""

"""
import os
import re
import configparser
import operator as op
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sb
from seaborn.distributions import _freedman_diaconis_bins
import click


from bcmproteomics_ext import ispec
sb.set_context('notebook', font_scale=1.8)


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


def pearson_r(x, y):
    return stats.pearsonr(x, y)[0]

def spearman_r(x, y):
    return stats.spearmanr(x, y)[0]

def hist(x, xmin=None, xmax=None, **kwargs):
    if 'color' in kwargs:
        color = kwargs.pop('color')
    X = x[ (~x.isnull()) & ~(x.abs() == np.inf) ]
    plt.hist(X.values, color='grey', **kwargs)

    # sb.despine(ax=ax, left=True, bottom=True)
    ax = plt.gca()
    if xmin and xmax:
        ax.set_xlim((xmin, xmax))

N_COLORS = 100
r_colors = sb.color_palette("coolwarm", n_colors=N_COLORS+1)
STEP = .2

def remove_ticklabels(fig=None, ax=None):
    # idea via seaborn/utils.py :: despine

    # Get references to the axes we want
    if fig is None and ax is None:
        axes = plt.gcf().axes
    elif fig is not None:
        axes = fig.axes
    elif ax is not None:
        axes = [ax]

    for ax_i in axes:
        ax_i.set_xticklabels([])
        ax_i.set_yticklabels([])

def plot_cbar(ax):

    labels = list(reversed(['{:.1f}'.format(x) for x in np.arange(1, -1.1, -STEP)]))
    cmap = mpl.colors.ListedColormap(r_colors)
    cbar = mpl.colorbar.ColorbarBase( ax, cmap=cmap)
    cbar.set_ticks(np.arange(0, 2.1, STEP/2))
    cbar.set_ticklabels(labels)
    ax.set_ylabel('r value')

def make_xaxis(ax, yloc=0, offset=0.05, fmt_str='%1.1f', **props):
    xmin, xmax = ax.get_xlim()
    locs = [loc for loc in ax.xaxis.get_majorticklocs()
            if loc >= xmin and loc <= xmax]
    tickline, = ax.plot(locs, [yloc]*len(locs), linestyle='',
                        marker=mpl.lines.TICKDOWN, **props)
    axline, = ax.plot([xmin, xmax], [yloc, yloc], **props)
    tickline.set_clip_on(False)
    axline.set_clip_on(False)
    for loc in locs:
        ax.text(loc, yloc - offset, fmt_str % loc,
                horizontalalignment='center',
                verticalalignment='top')



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
        r = pearson_r(X,Y)
        text = 'Pearson'
    elif stat == 'spearman':
        r = spearman_r(X,Y)
        text = 'Spearman'
    text = 'n = {:,}\nr = {:.2f}'.format(len(X), r)


    ax_bg_ix = int(round(r+1, 2) * N_COLORS/2 )  # add 1 to shift from -1 - 1 to 0 - 2 for indexing
    ax_bg = r_colors[ax_bg_ix]
    ax.patch.set_facecolor(ax_bg)
    # kwargs['text'] = text

    func(X, Y, ax, text=text, **kwargs)


def annotate_stat(x, y, ax, text, **kwargs):

    ax.annotate(text, xy=(0.5, 0.5), xycoords='axes fraction',
                va='center', ha='center', size=30
    )
    # sb.despine(ax=ax, left=True, bottom=True)

def scatter(x, y, ax, xymin, xymax, **kwargs):
    try:
        kwargs.pop('text')
    except KeyError:
        pass

    ax.scatter(x, y, c='k', **kwargs)
    # sb.despine(ax=ax, left=True, bottom=True)

    if xymin and xymax:
        ax.set_xlim((xymin, xymax))
        ax.set_ylim((xymin, xymax))

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

def make_config(path='.'):
    config = configparser.ConfigParser()
    config['Name'] = OrderedDict((('recno', 12345),
                                  ('runno', 1),
                                  ('searchno', 1)))
    file_ = os.path.join(path, 'generic_config.ini')
    with open(file_, 'w') as cf:
        print('Creating ', file_)
        config.write(cf)
        cf.write('#runno and searchno are optional, default to 1\n')

def read_config(configfile):
    """reads config file and returns the data"""
    config = configparser.ConfigParser()
    with open(configfile, 'r') as f:
        config.read_file(f)

    sections = config.sections()  # retains order
    FIELDS = ('recno', 'runno', 'searchno',)

    data = defaultdict(lambda : dict(runno=1, searchno=1))  # does not retain order (no guarantee)
    for section_key in sections:
        section = config[section_key]
        for field in FIELDS:
            value = section.get(field)
            if value is None:
                continue
            data[section_key][field] = value
        if 'recno' not in data[section_key]:  # record number is required
            print(section_key, 'does not have recno defined, skipping')
            data.pop(section_key)

    ordered_data = OrderedDict()
    for key in sections:
        ordered_data[key] = data[key]



    return ordered_data
