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




# def plot_pc(df, annot_text=False):

#     texts = list()

#     fig, ax = plt.subplots()
#     for name, row in df.iterrows():  # plot first and second components
#         color  = row['_color']
#         marker = row['_marker']
#         ax.scatter( row[0], row[1], color=color, marker=marker )

#         if annot_text:
#             texts.append( ax.text(row[0], row[1], row.name, size=6 ) )

#         if color_label:
#             name = row[color_label]
#             if name not in color_labels:
#                 color_labels.append(name)
#                 color_handle = matplotlib.patches.Patch(color=color)
#                 color_handles.append(color_handle)
#         if marker_label:
#             name = row[marker_label]
#             if name not in marker_labels:
#                 marker_labels.append(name)
#                 marker_handle = plt.Line2D((0, 1), (0, 0), color='gray', marker=marker,
#                                            linestyle='')
#                 marker_handles.append(marker_handle)

#     if annot_text:
#         res = adjust_text(texts, arrowprops=dict(arrowstyle="-", color='grey', lw=0.5),
#                           force_points=0.1, expand_text=(1.2, 1.2), expand_points=(1.4, 1.4)
#         )

#     legends = list()
#     if color_label:
#         leg1 = ax.legend(color_handles, color_labels, title=color_label, bbox_to_anchor=(1.04, 1),
#                          loc='upper left',
#                          borderaxespad=0)
#         legends.append(leg1)
#     if marker_label:
#         leg2 = ax.legend(marker_handles, marker_labels, loc='lower left', title=marker_label,
#                          bbox_to_anchor=(1.04, 0), borderaxespad=0)
#         legends.append(leg2)

#     ax.set_xlabel('PC1 ({:.2%})'.format(var1))
#     ax.set_ylabel('PC2 ({:.2%})'.format(var2))

#     for leg in legends:
#         ax.add_artist(leg)

#     ax.axhline(color='grey')
#     ax.axvline(color='grey')
#     fig.tight_layout(rect=[0,0,0.75,1])


#     return fig, ax


# def pcaplot(X, metadata=None, col_data=None):
#     """
#     looks for __PCA__ section in input config.ini file
#     """

#     rc = {'font.family': 'sans-serif',
#           "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans",
#                               "Bitstream Vera Sans", "sans-serif"],
#     }

#     sb.set_context('notebook')
#     sb.set_palette('muted')
#     sb.set_color_codes()
#     sb.set_style('white', rc)

#     if col_data is not None:
#         col_data = col_data.T
#         for col in col_data:
#             col_data[col] = pd.Categorical(col_data[col])
#     to_drop = [x for x in col_data.columns if x.startswith('_')]
#     col_data = col_data.drop(to_drop, axis=1)

#     pca_params = metadata.get('__PCA__')
#     annot_text = pca_params.get('annot')
#     if annot_text is not None:
#         annot_text = True if 'true' in annot_text.lower() else False

#     X_centered = X.sub(X.mean(1), axis='index')

#     U, s, V = np.linalg.svd(X_centered)

#     eigen = s**2
#     sumvariance = np.cumsum(eigen)
#     sumvariance /= sumvariance[-1]

#     # no current support for multiple multiplexed samples
#     components = pd.DataFrame(data=V, columns=X.columns)  # should be same as col_data.index except removal of (any) experiments  with no data
#     try:
#         # components = pd.DataFrame(data=V, columns=col_data.index)
#         df = col_data.join(components.T)
#     except ValueError:
#         components = pd.DataFrame(data=V, columns=col_data.columns)
#         mapping = col_data.to_dict(orient='records')[0]
#         df = components.rename(columns=mapping).T
#         df['color'] = pd.Categorical(df.index)
#         # df = col_data.T.join(components.T)
#         pca_params = dict(color = 'color')

#     # var1, var2, *rest = [x/s.sum() for x in s]
#     var1, var2, *rest = s**2 / (s**2).sum()


#     if pca_params is not None and 'color' in pca_params:
#         color_label = pca_params.get('color')
#         try:
#             n_colors = df[color_label].nunique()
#         except KeyError:
#             warnings.warn('The label {} is not in the metadata'.format(color_label))
#             n_colors = 1
#     else:
#         n_colors = 1
#         color_label = None

#     if n_colors <= 10:
#         colors = sb.color_palette('tab10', n_colors=10)
#     else:
#         colors = sb.color_palette('cubehelix', n_colors=n_colors)

#     if color_label is not None:
#         color_mapper = [ colors[ix] for ix in df[color_label].cat.codes ]
#     else:
#         color_mapper = [ colors[0] for _ in df.index ]
#     df['_color'] = color_mapper

#     if pca_params is not None:
#         marker_label = pca_params.get('marker')
#     else:
#         marker_label = None

#     my_markers = ('o', 'v', 's', 'd', '*', 'X', 'P', 'h', '<', 'H', 'D', '>', 'p', '^', )

#     if marker_label:
#         try:
#             marker_mapper = [ my_markers[ix]
#                               for ix in df[marker_label].cat.codes ]
#         except IndexError:
#             raise TooManyCategories('Not enough marker shapes')
#     else:
#         marker_mapper = [ 'o' for _ in df.index ]
#     df['_marker'] = marker_mapper

#     # shade_label  = pca_params.get('shade')

#     color_handles, color_labels   = list(), list()
#     marker_handles, marker_labels = list(), list()

#     import ipdb; ipdb.set_trace()
#     figs = dict()
#     fig, ax = plot_pc(df, annot_text=annot_text)
#     figs['pcaplot_1_2'] = fig

#     return figs

class PCAplot:

    markers = ('o', 'v', 's', 'd', '*', 'X', 'P', 'h', '<', 'H', 'D', '>', 'p', '^', )

    def __init__(self, X, metadata, col_data, annotate=False):
        self.X = X
        self.metadata = metadata
        self.annotate = annotate

        if col_data is not None:
            col_data = col_data.T
            for col in col_data:
                col_data[col] = pd.Categorical(col_data[col])
        to_drop = [x for x in col_data.columns if x.startswith('_')]
        self.col_data = col_data.drop(to_drop, axis=1)

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

        pca_params = metadata.get('__PCA__')
        if pca_params is not None and 'color' in pca_params:
            color_label = pca_params.get('color')
            try:
                n_colors = df[color_label].nunique()
            except KeyError:
                warnings.warn('The label {} is not in the metadata'.format(color_label))
                n_colors = 1
        else:
            n_colors = 1
            color_label = None
        self.color_label = color_label

        if n_colors <= 10:
            colors = sb.color_palette('tab10', n_colors=10)
        else:
            colors = sb.color_palette('cubehelix', n_colors=n_colors)

        if color_label is not None:
            color_mapper = [ colors[ix] for ix in df[color_label].cat.codes ]
        else:
            color_mapper = [ colors[0] for _ in df.index ]
        df['_color'] = color_mapper

        if pca_params is not None:
            marker_label = pca_params.get('marker')
        else:
            marker_label = None
        self.marker_label = marker_label

        if marker_label:
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

        fig, ax = plt.subplots()
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
                    color_handle = matplotlib.patches.Patch(color=color)
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
        if color_label:
            leg1 = ax.legend(color_handles, color_labels, title=color_label, bbox_to_anchor=(1.04, 1),
                            loc='upper left',
                            borderaxespad=0)
            legends.append(leg1)
        if marker_label:
            leg2 = ax.legend(marker_handles, marker_labels, loc='lower left', title=marker_label,
                            bbox_to_anchor=(1.04, 0), borderaxespad=0)
            legends.append(leg2)

        var1, var2 = self.vars[x-1], self.vars[y-1]
        ax.set_xlabel('PC{} ({:.2%})'.format(x, var1))
        ax.set_ylabel('PC{} ({:.2%})'.format(y, var2))

        for leg in legends:
            ax.add_artist(leg)

        ax.axhline(color='grey')
        ax.axvline(color='grey')
        fig.tight_layout(rect=[0,0,0.75,1])


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

def pcaplot(X, metadata=None, col_data=None, annotate=False, max_pc=2):

    rc = {'font.family': 'sans-serif',
        "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans",
                            "Bitstream Vera Sans", "sans-serif"],
        'legend.frameon': False,
    }
    sb.set_context('talk')
    sb.set_style('white', rc)

    pca = PCAplot(X, metadata, col_data, annotate=annotate)

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
