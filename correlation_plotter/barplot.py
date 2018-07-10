from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import seaborn as sb

from .utils import save_multiple

def plotter(data, color=None, cmap='tab10'):

    w = len(data.columns) * .75
    w = max(w, 4)


    fig, ax = plt.subplots(figsize=(w,5))
    data.loc['mean'].plot.bar(yerr=data.loc['std'], rot=0, ax=ax)
    fig.autofmt_xdate()
    return fig, ax

def barplot(X, genes, metadata, average=None, color=None, cmap=None, linear=False, base_outfunc=None, file_fmts=('.png',)):
    """
    Plot barplots for each gene in a list of genes
    :file_fmts: iterable of valid file formats for matplotlib to save figures
    """
    if cmap is None:
        cmap = 'tab10'

    if color:
        color_groups = metadata.T.groupby(color).groups
        colors = dict()
        sb.color_palette(color, n_colors=len(color_groups))




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
        fig, ax = plotter(df, color=None, cmap=cmap)
        if average:
            name = '{}_{}'.format(str(gene), average)
        else:
            name = str(gene)
        outname = base_outfunc(name)
        save_multiple(fig, outname, *file_fmts, dpi=300)
