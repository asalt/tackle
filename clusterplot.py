import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import rgb2hex

import numpy as np
import pandas as pd
import seaborn as sb
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

from utils import *

rc = {'font.family': 'serif',
      'font.serif': ['Times', 'Palatino', 'serif']}
sb.set_context('paper')
sb.set_style('white', rc)
sb.set_palette('muted')
sb.set_color_codes()

def _calculate_box_sizes(size_vector):
    """
    returns optimal box widths for
    different number of labels
    input is size vector such as
    [1, 5, 1]
    Return optimized sizes adding up to 1:
    [0.14, 0.71, 0.14]
    """
    size_vector = np.array(size_vector)
    sizes = size_vector / size_vector.sum()
    cumsum = np.cumsum(sizes)
    start = [0, *cumsum][:-1]
    return start

def plot_silhouette_scores(scores, start, end):

    fig, ax = plt.subplots()
    ax.plot( range(start, end+1), scores, marker='o' )

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel('Silhouette Score')
    ax.set_xlabel('Number of Clusters')

    return fig, ax


def calc_optimal_clusters(data, start=2, end=20):
    best_score = -np.inf
    best_cluster = None

    scores = list()

    for i in range(start, end+1):
        kmeans = KMeans(n_clusters=i).fit(data)

        score = silhouette_score(data, kmeans.labels_)

        scores.append(score)

        if score > best_score:
            best_score = score
            best_cluster = i

    fig, ax = plot_silhouette_scores(scores, start, end)

    return best_cluster, fig, ax

def silhouette_plot(data, labels):

    fig, ax = plt.subplots()
    n_clusters = len(set(labels))

    # ax.set_xlim([-0.1, 1])
    ax.set_ylim([0, len(data) + (n_clusters + 1) * 10])

    silhouette_avg = silhouette_score(data, labels)
    sample_silhouette_values = silhouette_samples(data, labels)

    y_lower = 10

    cmap = iter(sb.color_palette('hls', n_colors=max(6, n_clusters)))
    cmap_mapping = {val : rgb2hex(next(cmap)) for val in range(n_clusters)}


    for i in range(n_clusters):
        ith_cluster_silhouette_values = sample_silhouette_values[labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cmap_mapping[i]
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title('Silhouette Plot for the Various Clusters.')
    ax.set_xlabel("Silhouette Coefficient Values")
    ax.set_ylabel("Cluster Label")
    ax.set_yticks([])  # Clear the yaxis labels / ticks

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    return fig, ax


def calc_kmeans(data, nclusters, seed=None):

    autofig, autoax = None, None
    if nclusters == 'auto':
        nclusters, autofig, autoax = calc_optimal_clusters(data)

    kmeans = KMeans(n_clusters=nclusters, random_state=seed).fit(data)

    fig, ax = silhouette_plot(data, kmeans.labels_)


    ret = {'nclusters': nclusters, 'auto': {'fig': autofig, 'ax': autoax},
           'silhouette': {'fig': fig, 'ax': ax},
           'kmeans': kmeans,
           'nclusters': nclusters
    }

    return ret



def clusterplot(data, highlight_gids=None, highlight_gid_names=None,
                gid_symbol=None, nclusters=None, gene_symbols=None, z_score=None, standard_scale=None, row_cluster=True, seed=None,
                col_cluster=True, metadata=None, col_data=None):
    """
    :nclusters: None, 'auto', or positive integer

    """
    retval = dict()

    row_colors = None
    extra_artists = None
    if highlight_gids:

        cmap = sb.color_palette('hls', n_colors=max(6, len(highlight_gids)))
        colors_dfs = list()
        for ix, hgid in enumerate(highlight_gids):
            color = cmap[ix]
            highlights = {gid: color for gid in hgid}
            colors = [highlights.get(x, (1., 1., 1.)) for x in data.index]
            # colors = data.index.map( lambda x: highlights.get(x, 'white') )
            colors_df = pd.Series(colors, index=data.index).to_frame(highlight_gid_names[ix])
            colors_dfs.append(colors_df)
        row_colors = pd.concat(colors_dfs,axis=1)
        # row_colors.columns = highlight_gid_names

    col_colors = None
    # if metadata is not None:
    if col_data is not None:
        # col_data = parse_metadata(metadata)
        col_colors = pd.DataFrame(columns=col_data.columns,
                                  index=col_data.index)

        for info in col_data.index:
            cmap = iter(sb.color_palette('hls', n_colors=max(6, col_data.loc[info].nunique())))
            mapping = {val : next(cmap) for val in col_data.loc[info].unique()}
            colors = col_data.loc[info].map(mapping)
            col_colors.loc[info] = colors
        col_colors = col_colors.T


    if gene_symbols:
        clustermap_symbols = [gid_symbol.get(x, '?') for x in data.index]
        data.index = clustermap_symbols
        if row_colors is not None:
            row_colors.index = clustermap_symbols

    # figheight = min(len(data) / 6, 100)
    figheight = 12
    figwidth  = len(data.columns) / 2
    # if col_colors is not None:
    #     figwidth -= (max(len(x) for x in col_colors.columns) * .16667)

    if nclusters is not None:

        kmeans_result = calc_kmeans(data, nclusters, seed)
        kmeans = kmeans_result['kmeans']
        clusters = pd.Series(data=kmeans.labels_, index=data.index)
        cluster_order = clusters.sort_values().index

        plot_data = data.loc[cluster_order]

        cmap = iter(sb.color_palette('hls', n_colors=max(6, kmeans.n_clusters)))
        cmap_mapping = {val : rgb2hex(next(cmap)) for val in range(kmeans.n_clusters)}
        cluster_colors = clusters.map(cmap_mapping).to_frame('Cluster')

        cluster_data = data.assign(Cluster=clusters)
        cluster_data['silhouette_score'] = silhouette_samples(data, kmeans.labels_)


        if row_colors is None:
            row_colors = cluster_colors
        else:
            row_colors = pd.concat([row_colors, cluster_colors])

        row_cluster = False
        kmeans_result['data'] = cluster_data
        retval['kmeans'] = kmeans_result

    else:
        plot_data = data



    g = sb.clustermap(plot_data,
                      row_colors=row_colors if row_colors is not None and not row_colors.empty else None,
                      col_colors=col_colors if col_colors is not None and not col_colors.empty else None,
                      yticklabels=False if not gene_symbols else True,
                      z_score=z_score, standard_scale=standard_scale,
                      figsize=(figwidth, figheight),
                      row_cluster=row_cluster, col_cluster=col_cluster,
                      cmap = 'YlOrRd' if z_score is None else 'RdBu_r',
                      center = 0 if z_score is not None else None
    )
    if gene_symbols:
        for tick in g.ax_heatmap.yaxis.get_ticklabels():
            tick.set_rotation(0)
            tick.set_size(tick.get_size()*.4)
    for tick in g.ax_heatmap.xaxis.get_ticklabels():
        tick.set_rotation(90)
    if g.ax_row_colors:
        ticks = g.ax_row_colors.xaxis.get_ticklabels()
        if len(ticks) > 1:
            scale = .9/ (np.log1p(len(ticks)) *.9)
            for tick in ticks:
                tick.set_size(tick.get_size()*scale)


    if col_colors is not None:
        widths = _calculate_box_sizes( col_colors.nunique() )
        col_colors_t = col_colors.T
        bboxes = [(x, 1.02, 1, .2) for x in widths]
        # bboxes = [(x, 1.02, 1, .2) for x in np.arange(0, 1, 1/len(col_colors_t.index))]
        legends = list()
        for bbox, ix in zip(bboxes, col_colors_t.index):
            col_name        = ix
            col_labels      = col_data.loc[ix].sort_values().drop_duplicates()
            col_names       = col_labels.values
            label_colors    = col_colors_t.loc[ix, col_labels.index].values
            handles, labels = list(), list()
            for n, c in zip(col_names, label_colors):
                handle = matplotlib.patches.Patch(color=c,)
                handles.append(handle)
                labels.append(n)
            leg = g.ax_col_dendrogram.legend( handles, labels, bbox_to_anchor=bbox,
                                              loc='lower left', ncol=len(col_names) // 2,
                                              title=col_name
            )
            legends.append(leg)
            g.ax_col_dendrogram.add_artist(leg)
            extra_artists = legends

    retval['clustermap'] = dict(clustergrid=g, extra_artists=extra_artists)

    return retval
    # return g, extra_artists
