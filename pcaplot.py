import warnings

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import markers
from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import seaborn as sb

from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.pipeline import Pipeline, make_pipeline


from utils import *

rc = {'font.family': 'serif',
      'font.serif': ['Times', 'Palatino', 'serif']}
sb.set_context('paper')
sb.set_style('white', rc)

def pcaplot(X, metadata=None, col_data=None):
    """
    looks for __PCA__ section in input config.ini file
    """
    # print(len(X))

    # if metadata is not None:
    if col_data is not None:
        # col_data = parse_metadata(metadata).T
        col_data = col_data.T
        for col in col_data:
            col_data[col] = pd.Categorical(col_data[col])
    to_drop = [x for x in col_data.columns if x.startswith('_')]
    col_data = col_data.drop(to_drop, axis=1)

    pca_params = metadata.get('__PCA__')

    # l2 normalization down genes
    # X_normed = pd.DataFrame(data=normalize(X, axis=0), columns=X.columns, index=X.index)


    X_centered = X.sub(X.mean(1), axis='index')
    U, s, V = np.linalg.svd(X_centered)

    eigen = s**2
    sumvariance = np.cumsum(eigen)
    sumvariance /= sumvariance[-1]

    # X_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(X.T).T

    # X_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(X_normed.T).T

    # pca = PCA(n_components=min(100, len(X))).fit(X.T)
    # pca = PCA(n_components=min(100, len(X))).fit(X_scaled.T)

    # pca = TruncatedSVD(min(100, len(X))).fit(X.T)
    # pca = TruncatedSVD(min(100, len(X))).fit(X_scaled.T)

    # pca = PCA().fit(X)

    # components_ is array of [n_components, n_features]
    # components = pd.DataFrame(data=pca.components_[:, 0:8], columns=col_data.index)
    # var1, var2, *rest = pca.explained_variance_ratio_

    # no current support for multiple multiplexed samples
    try:
        components = pd.DataFrame(data=V, columns=col_data.index)
        df = col_data.join(components.T)
    except ValueError:
        components = pd.DataFrame(data=V, columns=col_data.columns)
        mapping = col_data.to_dict(orient='records')[0]
        df = components.rename(columns=mapping).T
        df['color'] = pd.Categorical(df.index)
        # df = col_data.T.join(components.T)
        pca_params = dict(color = 'color')

    var1, var2, *rest = [x/s.sum() for x in s]


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

    colors = sb.color_palette('tab10', n_colors=max(10, n_colors))
    if color_label is not None:
        color_mapper = [ colors[ix] for ix in df[color_label].cat.codes ]
    else:
        color_mapper = [ colors[0] for _ in df.index ]
    df['_color'] = color_mapper

    if pca_params is not None:
        marker_label = pca_params.get('marker')
    else:
        marker_label = None

    my_markers = ('o', 'v', 's', 'd', '*', 'X', 'P', 'h', '<', 'H', 'D', '>', 'p', '^', )

    if marker_label:
        try:
            marker_mapper = [ my_markers[ix]
                              for ix in df[marker_label].cat.codes ]
        except IndexError:
            raise TooManyCategories('Not enough marker shapes')
    else:
        marker_mapper = [ 'o' for _ in df.index ]
    df['_marker'] = marker_mapper

    # shade_label  = pca_params.get('shade')

    color_handles, color_labels   = list(), list()
    marker_handles, marker_labels = list(), list()

    fig, ax = plt.subplots()
    for name, row in df.iterrows():  # plot first and second components
        color  = row['_color']
        marker = row['_marker']
        ax.scatter( row[0], row[1], color=color, marker=marker )

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

    ax.set_xlabel('PC1 ({:.2%})'.format(var1))
    ax.set_ylabel('PC2 ({:.2%})'.format(var2))

    for leg in legends:
        ax.add_artist(leg)

    ax.axhline(color='grey')
    ax.axvline(color='grey')
    fig.tight_layout(rect=[0,0,0.75,1])


    return fig, ax


    # comp1, comp2, *rest = pca.components_
