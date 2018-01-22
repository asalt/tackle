import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import rgb2hex

import numpy as np
import pandas as pd
import seaborn as sb
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_samples, silhouette_score

from .utils import *
from .containers import MyClusterGrid

rc = {'font.family': 'sans-serif',}
      # 'font.serif': ['Times', 'Palatino', 'serif']}
sb.set_context('notebook')
sb.set_palette('muted')
sb.set_color_codes()
sb.set_style('white', rc)
# mpl.rcParams.update(rc)

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


    for i in reversed(range(n_clusters)):
        ith_cluster_silhouette_values = sample_silhouette_values[labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cmap_mapping[i]
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i+1))

        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title('Silhouette Plot for the Various Clusters')
    ax.set_xlabel("Silhouette Coefficient Values")
    ax.set_ylabel("Cluster Label")
    ax.set_yticks([])  # Clear the yaxis labels / ticks

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")
    sb.despine(fig=fig, left=True, bottom=True)

    return fig, ax


def calc_kmeans(data, nclusters, seed=None, max_autoclusters=30):

    autofig, autoax = None, None
    if nclusters == 'auto':
        nclusters, autofig, autoax = calc_optimal_clusters(data, end=max_autoclusters)

    kmeans = KMeans(n_clusters=nclusters, random_state=seed).fit(data)

    fig, ax = silhouette_plot(data, kmeans.labels_)


    ret = {'nclusters': nclusters, 'auto': {'fig': autofig, 'ax': autoax},
           'silhouette': {'fig': fig, 'ax': ax},
           'kmeans': kmeans,
           'nclusters': nclusters
    }

    return ret



def clusterplot(data, dbscan=False, highlight_gids=None, highlight_gid_names=None, gid_symbol=None,
                nclusters=None, gene_symbols=None, z_score=None, standard_scale=None, mask=None,
                show_missing_values=True, max_autoclusters=30, row_cluster=True,
                seed=None, col_cluster=True, metadata=None, col_data=None, figsize=None,
                normed=False, linkage='average'):
    """
    :nclusters: None, 'auto', or positive integer

    """

    rc = {'font.family': 'sans-serif',}
    sb.set_context('notebook')
    sb.set_palette('muted')
    sb.set_color_codes()
    sb.set_style('white', rc)

    retval = dict()


    if dbscan or nclusters:  # do not perform hierarchical clustering and KMeans (or DBSCAN)
        row_cluster = False

    row_colors = None
    extra_artists = None
    if highlight_gids:

        cmap = [rgb2hex(x) for x in
                sb.color_palette('hls', n_colors=max(6, len(highlight_gids)))
        ]
        colors_dfs = list()
        for ix, hgid in enumerate(highlight_gids):
            color = cmap[ix]
            highlights = {gid: color for gid in hgid}
            colors = [highlights.get(x, '#ffffff') for x in data.index]
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
            cmap = iter(rgb2hex(x) for x in
                        sb.color_palette('hls', n_colors=max(6, col_data.loc[info].nunique())))
            mapping = {val : next(cmap) for val in col_data.loc[info].unique()}
            colors = col_data.loc[info].map(mapping)
            col_colors.loc[info] = colors
        col_colors = col_colors.T

    if gene_symbols:  # change index to symbols
        assert all(data.index == mask.index)
        clustermap_symbols = [gid_symbol.get(x, '?') for x in data.index]
        data.index = clustermap_symbols
        mask.index = clustermap_symbols
        if row_colors is not None:
            row_colors.index = clustermap_symbols


    # if nclusters is not None or dbscan:
    if z_score is not None:
        data_t = sb.matrix.ClusterGrid.z_score(data, z_score)
    else:
        data_t = data


    if nclusters is not None:

        kmeans_result = calc_kmeans(data_t, nclusters, seed, max_autoclusters)
        kmeans = kmeans_result['kmeans']
        clusters = pd.Series(data=kmeans.labels_, index=data.index)
        cluster_order = clusters.sort_values().index

        # plot_data = data.loc[cluster_order]
        plot_data = data_t.loc[cluster_order]

        cmap = iter(rgb2hex(x) for x in sb.color_palette('hls', n_colors=max(6, kmeans.n_clusters)))
        cmap_mapping = {val : next(cmap) for val in range(kmeans.n_clusters)}
        cluster_colors = clusters.map(cmap_mapping).to_frame('Cluster')

        cluster_data = data_t.copy()
        if show_missing_values:
            # cluster_data[mask] = np.NAN
            cluster_data[mask] = 0
        cluster_data = cluster_data.assign(Cluster=clusters+1)
        cluster_data['silhouette_score'] = silhouette_samples(data_t, kmeans.labels_)

        if row_colors is None:
            row_colors = cluster_colors
        else:
            row_colors = pd.concat([cluster_colors, row_colors], axis=1)

        kmeans_result['data'] = cluster_data
        retval['kmeans'] = kmeans_result


    elif dbscan:
        db = DBSCAN().fit(data_t)
        clusters = pd.Series(data=db.labels_, index=data.index)
        n_clusters = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0)
        if n_clusters == 0:
            raise ValueError('No clusters found!')
        cluster_order = clusters.sort_values().index

        plot_data = data_t.loc[cluster_order]

        cmap = iter(sb.color_palette('hls', n_colors=max(6, n_clusters)))
        cmap_mapping = {val : rgb2hex(next(cmap)) for val in range(n_clusters)}
        cmap_mapping[-1] = 'k'
        cluster_colors = clusters.map(cmap_mapping).to_frame('Cluster')

        cluster_data = data_t.assign(Cluster=clusters+1)
        cluster_data['silhouette_score'] = silhouette_samples(data_t, db.labels_)
        cluster_data.loc[cluster_data.Cluster == -1, 'silhouette_score'] = np.nan

        if row_colors is None:
            row_colors = cluster_colors
        else:
            row_colors = pd.concat([row_colors, cluster_colors])

        valid_gids = clusters.where(lambda x: x != -1).dropna().index
        fig, ax = silhouette_plot(data.loc[valid_gids],
                                  clusters.loc[valid_gids].values)

        retval['dbscan'] = {'nclusters': n_clusters,
                            'silhouette': {'fig': fig, 'ax': ax},
                            'data' : cluster_data
        }

    else:
        plot_data = data_t

    cmap_name = 'YlOrRd' if (z_score is None and normed == False) else 'RdBu_r'
    cmap = cmap_name
    # cmap = mpl.cm.get_cmap(cmap_name)
    # robust = True if normed == True else False
    robust = False
    if (z_score is not None or normed==True):  # adapted from seaborn.matrix.ClusterMap
        center = 0
        # vmin = np.percentile(plot_data, 2) if robust else plot_data.min().min()
        # vmax = np.percentile(plot_data, 98) if robust else plot_data.max().max()
    #     low, high = 2, 98
        if normed:
            robust = True
    #         low, high = 30, 70
    #     vmin = np.percentile(plot_data, low) if robust else plot_data.min().min()
    #     vmax = np.percentile(plot_data, high) if robust else plot_data.max().max()
    #     # maxval = max(abs(vmin), vmax)
    #     # vmax = maxval
    #     # vmin = -maxval

    #     vrange = max(vmax - center, center - vmin)

    #     normalize = mpl.colors.Normalize(center - vrange, center + vrange)
    #     # cmin, cmax = normlize([vmin, vmax])
    #     cmin, cmax = normalize([-vrange, vrange])
    #     cc = np.linspace(cmin, cmax, 256)
    #     cmap = mpl.colors.ListedColormap(cmap(cc))

    # cmap.set_bad(color='gray')

    # figheight = min(len(data) / 6, 100)
    heatmap_height_ratio = .8  # this is the default (seaborn). Needs to be increased when figs are very long

    FONTSIZE = 8

    if figsize is None:

        figheight = 12
        heatmap_height_ratio = .8  # this is the default from Seaborn

        if gene_symbols:  # make sure there is enough room for the symbols
            figheight = max(((FONTSIZE+2)/72) * len(plot_data), 12)
            if figheight > 218: # maximum figheight in inches
                FONTSIZE = max(218 / figheight, 6)
                figheight = 218


        min_figwidth = 4
        if col_colors is not None and not col_colors.empty:
            for _ in range(1, len(col_colors.columns)):
                min_figwidth += 1  # make room for more labels in legend after 2 labels
        if gene_symbols:
            min_figwidth += 1 # for labels?

        figwidth  = max( min( len(data.columns) / 2, 8), min_figwidth )
        # if col_colors is not None:
        #     figwidth -= (max(len(x) for x in col_colors.columns) * .16667)
        figsize = (figwidth, figheight)

    dendrogram_width_ratio = None
    heatmap_width_ratio = .8  # default
    if gene_symbols:  # make sure there is enough room for the symbols
        heatmap_height_ratio = max(.8, .8 * (.2*(figheight-12)) )
        # heatmap_width_ratio = 1.7
        # dendrogram_width_ratio = .06
        # print(figheight, heatmap_height_ratio)
        pass

    if figheight > 12:
        dendrogram_width_ratio  = .16 + max(0,
                                            (.1*np.log10(figwidth))*(figwidth-12)
        )

    # a minimal subclass of seaborn ClusterGrid for scaling.
    plotter = MyClusterGrid(plot_data,
                            figsize=figsize,
                            row_colors=row_colors if row_colors is not None and not row_colors.empty else None,
                            col_colors=col_colors if col_colors is not None and not col_colors.empty else None,
                            mask=mask.loc[plot_data.index] if show_missing_values else None,
                            heatmap_height_ratio=heatmap_height_ratio,
                            dendrogram_width_ratio=dendrogram_width_ratio,
                            heatmap_width_ratio=heatmap_width_ratio
    )

    g = plotter.plot(method=linkage, metric='euclidean',
                     row_cluster=row_cluster, col_cluster=col_cluster,
                     row_linkage=None, col_linkage=None,
                     colorbar_kws=None,
                     cmap=cmap,
                     # cmap='RdBu_r',
                     robust=True,
                     yticklabels=False if not gene_symbols else True,
                     center=0 if cmap != 'YlOrRd' else None,
                     vmax=plot_data.max().max() if cmap == 'YlOrRd' else None,
                     rasterized=True,
    )
    if figheight <= 12:
        hspace =.01
        wspace = .01
    else:
        hspace = .01 / (22*figheight)
        wspace = .01

    g.fig.subplots_adjust(hspace=hspace, wspace=wspace,
                          left=.5/figwidth, right=1-1./figwidth,
                          bottom=1/figheight, top=1-(1.4/figheight)
    )


    # g = sb.clustermap(plot_data,
    #                   row_colors=row_colors if row_colors is not None and not row_colors.empty else None,
    #                   col_colors=col_colors if col_colors is not None and not col_colors.empty else None,
    #                   yticklabels=False if not gene_symbols else True,
    #                   # z_score=z_score, standard_scale=standard_scale,
    #                   figsize=figsize,
    #                   row_cluster=row_cluster, col_cluster=col_cluster,
    #                   cmap=cmap,
    #                   mask=mask.loc[plot_data.index] if show_missing_values else None,
    #                   # center = 0 if z_score is not None else None
    # )

    if gene_symbols:
        for tick in g.ax_heatmap.yaxis.get_ticklabels():
            tick.set_rotation(0)
            # tick.set_size(tick.get_size()*.4)
            txt = tick.get_text()
            if len(txt) > 7 and len(txt) < 9:
                tick.set_size(FONTSIZE-1)
            elif len(txt) >= 9:
                tick.set_size(FONTSIZE-3)
            else:
                tick.set_size(FONTSIZE)
        g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
    for tick in g.ax_heatmap.xaxis.get_ticklabels():
        tick.set_rotation(90)
    if g.ax_row_colors:
        ticks = g.ax_row_colors.xaxis.get_ticklabels()
        if len(ticks) > 1:
            scale = .9/ (np.log1p(len(ticks)) *.9)
            for tick in ticks:
                tick.set_size(tick.get_size()*scale)

    if row_colors is not None and 'Cluster' in row_colors.columns:
        # annotate cluster numbers
        _positions = dict()
        _lastpos = 0
        for cluster, color in sorted(cmap_mapping.items(), reverse=False):
            length = row_colors[ row_colors['Cluster'] == color ].pipe(len)
            # _positions[cluster] = ( length+_lastpos ) // 2;
            _positions[  _lastpos + ( length // 2 ) ] = cluster
            _lastpos += length

        g.ax_row_colors.set_yticks(np.arange(0, row_colors.pipe(len), 1))
        _yticklabels = ['' if ix not in _positions
                        else _positions[ix] + 1
                        for ix, _ in enumerate(g.ax_row_colors.get_yticklabels())]
        g.ax_row_colors.set_yticklabels(_yticklabels, fontsize=12)
        g.ax_row_colors.tick_params(axis='y', pad=8)  # so not too close to spine
        g.ax_row_colors.spines["left"].set_position(("axes", 0.0)) # green one

    if col_colors is not None:
        col_label_lengths = col_data.applymap(len).max(1) + col_colors.nunique()
        # widths = _calculate_box_sizes( col_colors.nunique() )
        widths = _calculate_box_sizes( col_label_lengths )
        col_colors_t = col_colors.T
        bbox_y0 = 1.44 if col_cluster else .8
        bboxes = [(x, bbox_y0, 1, .2) for x in widths]  # (x0, y0, width, height)
        # bboxes = [(x, 1.02, 1, .2) for x in np.arange(0, 1, 1/len(col_colors_t.index))]
        legends = list()
        for bbox, ix in zip(bboxes, col_colors_t.index):
            col_name        = ix
            col_labels      = col_data.loc[ix].drop_duplicates()
            col_names       = col_labels.values
            label_colors    = col_colors_t.loc[ix, col_labels.index].values
            handles, labels = list(), list()
            for n, c in zip(col_names, label_colors):
                handle = mpl.patches.Patch(color=c,)
                handles.append(handle)
                labels.append(n)
            leg = g.ax_col_dendrogram.legend( handles, labels, bbox_to_anchor=bbox,
                                              loc='upper left', ncol=max(len(col_names) // 3, 1),
                                              title=col_name
            )
            legends.append(leg)
            g.ax_col_dendrogram.add_artist(leg)
            extra_artists = legends


    # make sure there is enough room on the right side for labels
    if col_colors is not None and not col_colors.empty:

        width, height = g.fig.get_size_inches()

        longest_label = max(len(x) for x in col_colors.columns) + 6  # add a little padding
        char_width = (430/1000) # approx via from https://www.math.utah.edu/~beebe/fonts/afm-widths.html
        longest_length = longest_label * char_width
        inch_shift = longest_length * 12/72  # 72 pts in an inch

        shift = 1 - (inch_shift / width)

        g.gs.update(right=shift)  # add some room on the right so everything fits

    g.ax_heatmap.yaxis.set_label('')


    retval['clustermap'] = dict(clustergrid=g, extra_artists=extra_artists)

    return retval
    # return g, extra_artists
