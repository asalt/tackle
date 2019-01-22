import warnings
import itertools

import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import markers
from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import seaborn as sb

from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.pipeline import Pipeline, make_pipeline

from adjustText import adjust_text

from .utils import *


class PCAplot:

    markers = ('o', 'v', 's', 'd', '*', 'X', 'P', 'h', '<', 'H', 'D', '>', 'p', '^', )

    def __init__(self, X, color_label=None, marker_label=None, col_data=None, metadata=None, annotate=False):
        """
        :X: DataFrame with columns as sample names and rows as GeneIDs. Index consists of GeneIDs
        :metadata: Dictionary of dictionaries of this format:
                   {'__PCA__': {'color': color_field, 'marker': marker_field}
        :col_data: DataFrame with columns as sample names and rows as metadata. Index consists of
                   metadata names:
                                A     B      C      D
                    treatment  ctrl  ctrl  treat  treat
                    batch         0     1      0      1
        """
        self.X = X
        self.metadata = metadata
        self.annotate = annotate
        if metadata.get('__PCA__') is not None:
            warnings.warn("""
            Specifying __PCA__ in config file is depreciated. Specify color, marker, and annot at the command line
            """
            )

        if col_data is not None:
            col_data = col_data.T
            for col in col_data:
                col_data[col] = pd.Categorical(col_data[col])
        to_drop = [x for x in col_data.columns if x.startswith('_')]
        self.col_data = col_data.drop(to_drop, axis=1)
        params_annotate = None

        if annotate == False:
            pca_params = metadata.get('__PCA__')
            if pca_params is not None:
                warnings.warn("""Specifying annotation in config file is depreciated,
                it can now be specified as a command line option `--annotate`""", DeprecationWarning)
                params_annotate = pca_params.get('annot')
            if params_annotate is not None:
                annotate = True if 'true' in params_annotate.lower() else False
            self.annotate = annotate

        X_centered = X.sub(X.mean(1), axis='index')


        U, s, V = np.linalg.svd(X_centered)
        self.s = s
        self.vars = s**2 / (s**2).sum()

        eigen = s**2
        sumvariance = np.cumsum(eigen)
        sumvariance /= sumvariance[-1]

        # no current support for multiple multiplexed samples
        components = pd.DataFrame(data=V, columns=X.columns)  # should be same as col_data.index except removal of (any) experiments  with no data

        try:
            # components = pd.DataFrame(data=V, columns=col_data.index)
            df = col_data.join(components.T)
        except ValueError:
            components = pd.DataFrame(data=V, columns=col_data.columns)
            mapping = col_data.to_dict(orient='records')[0]
            df = components.rename(columns=mapping).T
            df['color'] = pd.Categorical(df.index)
            # df = col_data.T.join(components.T)
            pca_params = dict(color = 'color')

        # var1, var2, *rest = [x/s.sum() for x in s]
        # var1, var2, *rest = s**2 / (s**2).sum()

        if self.metadata:
            pca_params = metadata.get('__PCA__')
            if pca_params is not None and 'color' in pca_params:
                color_label = pca_params.get('color')
                try:
                    n_colors = df[color_label].nunique()
                except KeyError:
                    warnings.warn('The label {} is not in the metadata'.format(color_label))
                    n_colors = 1
            if pca_params is not None and 'marker' in pca_params:
                marker_label = pca_params.get('marker')

        self.color_label = color_label
        if color_label:
            n_colors = df[color_label].nunique()
        else:
            n_colors = 1

        if n_colors <= 10:
            colors = sb.color_palette('tab10', n_colors=10)
        else:
            colors = sb.color_palette('cubehelix', n_colors=n_colors)

        if color_label:
            color_mapper = [ colors[ix] for ix in df[color_label].cat.codes ]
        else:
            color_mapper = [ colors[0] for _ in df.index ]
        df['_color'] = color_mapper

        self.marker_label = marker_label

        if self.marker_label:
            try:
                marker_mapper = [ self.markers[ix]
                                  for ix in df[marker_label].cat.codes ]
            except IndexError:
                raise TooManyCategories('Not enough marker shapes')
        else:
            marker_mapper = [ 'o' for _ in df.index ]
        df['_marker'] = marker_mapper

        self.df = df

    def plot_pc(self, x=1, y=2):

        df           = self.df
        color_label  = self.color_label
        marker_label = self.marker_label
        annotate     = self.annotate

        color_handles, color_labels   = list(), list()
        marker_handles, marker_labels = list(), list()
        texts = list()
        maxcols = 1
        for label in (color_label, marker_label):
            if not label:
                continue
            try:
                ncols = np.ceil( df[label].nunique()/5 )
                maxcols = max([ maxcols, ncols ])
            except ValueError:
                pass

        figsize = [6.4, 4.8]
        # figsize[1] += 1.4*(maxcols-1)


        fig, ax = plt.subplots(figsize=figsize)
        for name, row in df.iterrows():  # plot first and second components
            color  = row['_color']
            marker = row['_marker']
            ax.scatter( row[x-1], row[y-1], color=color, marker=marker )

            if annotate:
                texts.append( ax.text(row[x-1], row[y-1], row.name, size=6 ) )

            if color_label:
                name = row[color_label]
                if name not in color_labels:
                    color_labels.append(name)
                    # color_handle = matplotlib.patches.Patch(color=color)
                    color_handle = matplotlib.patches.Rectangle((0., 0.),   # (x,y)
                                                                0.5,          # width
                                                                0.5, color=color
                    )
                    color_handles.append(color_handle)
            if marker_label:
                name = row[marker_label]
                if name not in marker_labels:
                    marker_labels.append(name)
                    marker_handle = plt.Line2D((0, 1), (0, 0), color='gray', marker=marker,
                                            linestyle='')
                    marker_handles.append(marker_handle)

        if annotate:
            res = adjust_text(texts, arrowprops=dict(arrowstyle="-", color='grey', lw=0.5),
                            force_points=0.1, expand_text=(1.2, 1.2), expand_points=(1.4, 1.4)
            )

        legends = list()
        maxcol = 1
        if color_label:
            ncol = 1
            # if len(color_labels) > 5:
            #     ncol += 1
            fontsize = 12
            longest_label = max(len(x) for x in color_labels)
            if longest_label > 12:
                fontsize=8
            leg1 = ax.legend(color_handles, color_labels, title=color_label, bbox_to_anchor=(1.04, 1),
                             loc='upper left', ncol=ncol, labelspacing=.2, fontsize=fontsize,
                            borderaxespad=0)
            legends.append(leg1)
            maxcol = max([maxcol, ncol])
        if marker_label:
            ncol = 1
            # if len(marker_labels) > 5:
            #     ncol += 1
            leg2 = ax.legend(marker_handles, marker_labels, loc='lower left', title=marker_label,
                             bbox_to_anchor=(1.04, 0), borderaxespad=0, ncol=ncol,
                             labelspacing=.2
            )
            legends.append(leg2)
            maxcol = max([maxcol, ncol])

        var1, var2 = self.vars[x-1], self.vars[y-1]
        ax.set_xlabel('PC{} ({:.2%})'.format(x, var1))
        ax.set_ylabel('PC{} ({:.2%})'.format(y, var2))

        for leg in legends:
            ax.add_artist(leg)

        ax.axhline(color='grey')
        ax.axvline(color='grey')
        right = .75
        if maxcol < 1:
            right -= .25
        fig.tight_layout(rect=[0,0,.75,1])


        return fig, ax

    def plot_vars(self):

        sb.set_style("ticks")

        fig, ax = plt.subplots()

        VARMAX = 12  # don't plot more than 12 here, not necessary
        varmax = min([VARMAX, len(self.vars)])

        # ax.scatter( range( 1, len(self.vars)+1 ), self.vars, c='#333333' )
        ax.scatter( range( 1, varmax+1 ), self.vars[:varmax], c='#333333' )
        ax.grid(axis='y')
        ax.set_xlabel('Principal Component')
        for f in (ax.set_xticks, ax.set_xticklabels):
            # f( range(1, len(self.vars)+1) )
            f( range(1, varmax+1) )
        ax.set_ylabel('Variance')

        sb.despine(ax=ax, trim=True)
        fig.tight_layout()
        return fig, ax

def pcaplot(X, metadata=None, col_data=None, annotate=False, max_pc=2, color_label=None, marker_label=None, genes=None):
    rc = {'font.family': 'sans-serif',
          'text.usetex': False,
        # "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans",
        #                     "Bitstream Vera Sans", "sans-serif"],
        'legend.frameon': False,
    }
    sb.set_context('talk')
    sb.set_style('white', rc)

    if genes is not None:  # only plot these select genes
        _genes = set(genes) & set(X.index)
        X = X.loc[_genes]

    pca = PCAplot(X, color_label=color_label, marker_label=marker_label, metadata=metadata,
                  col_data=col_data, annotate=annotate)

    figs = dict()

    for a,b in itertools.combinations(range(1, max_pc+1), 2):

        try:
            fig, ax = pca.plot_pc(a, b)
            figs['pcaplot_{}_{}'.format(a, b)] = fig
        except IndexError:
            print('PC {} doesn\'t exist!')


    fig, ax = pca.plot_vars()
    figs['pcaplot_variance'] = fig

    return figs

pcaplot.__doc__ = PCAplot.__init__.__doc__ + """
        :max_pc: maximum number of principal components to plot

        returns
        :figs: dict with keys indicating plotted principal components and values each a matplotlib.figure.Figure
            which can be saved
"""
