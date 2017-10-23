import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import seaborn as sb

from utils import *

rc = {'font.family': 'serif',
      'font.serif': ['Times', 'Palatino', 'serif']}
sb.set_context('paper')
sb.set_style('white', rc)

def scatterplot(ibaqs_log_shifted, stat='pearson', colors_only=False, shade_correlation=True):

    # xymin = 0
    xymin = np.floor(ibaqs_log_shifted[ibaqs_log_shifted > 0].min().min())
    xymax = np.ceil(ibaqs_log_shifted.max().max())

    g = sb.PairGrid(ibaqs_log_shifted)
    g.map_upper(plot_delegator, stat=stat, filter_zeros=True, upper_or_lower='upper', colors_only=colors_only, shade_correlation=shade_correlation)
    g.map_lower(plot_delegator, stat=stat, filter_zeros=True, upper_or_lower='lower',
                xymin=xymin, xymax=xymax, colors_only=colors_only, shade_correlation=shade_correlation)
    g.map_diag(hist, xmin=xymin, xmax=xymax, colors_only=colors_only)
    if shade_correlation:
        color_diag(g)
    sb.despine(fig=g.fig, left=True, bottom=True)
    remove_ticklabels(fig=g.fig)
    #  adjust the spacing between subplots
    hspace = g.fig.subplotpars.hspace
    wspace = g.fig.subplotpars.wspace
    g.fig.subplots_adjust(hspace=hspace*.1, wspace=wspace*.1, right=.8, bottom=.2,
                          left=.2, top=.8)
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

    return g
