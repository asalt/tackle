import textwrap
from collections import OrderedDict
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import rgb2hex
import matplotlib.ticker as ticker


import numpy as np
import pandas as pd
import seaborn as sb
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_samples, silhouette_score
from scipy.cluster import hierarchy

from .utils import *
from .containers import MyClusterGrid, GeneMapper, MyDendrogramPlotter
sb.matrix._DendrogramPlotter = MyDendrogramPlotter


rc = {'font.family': 'sans-serif',
      'text.usetex': False,
      'svg.fonttype': 'none'
}
      # 'font.serif': ['Times', 'Palatino', 'serif']}
sb.set_context('notebook')
sb.set_palette('muted')
sb.set_color_codes()
sb.set_style('white', rc)
# mpl.rcParams.update(rc)

def _calculate_box_sizes(size_vector, start_pos=0, end_pos=1):
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
    cumsum = np.cumsum(sizes) * end_pos
    start = [start_pos, *cumsum][:-1]
    return start

def plot_silhouette_scores(scores, start, end):

    fig, ax = plt.subplots()
    ax.plot( range(start, end+1), scores, marker='o' )

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel('Silhouette Score')
    ax.set_xlabel('Number of Clusters')

    return fig, ax


def calc_optimal_clusters(data, start=4, end=20, random_state=None):
    best_score = -np.inf
    best_cluster = None

    scores = list()

    for i in range(start, end+1):
        kmeans = KMeans(n_clusters=i, random_state=random_state).fit(data)

        score = silhouette_score(data, kmeans.labels_, random_state=random_state)

        scores.append(score)

        if score > best_score:
            best_score = score
            best_cluster = i

    fig, ax = plot_silhouette_scores(scores, start, end)

    return best_cluster, fig, ax

def silhouette_plot(data, labels, random_state=None):

    fig, ax = plt.subplots()
    n_clusters = len(set(labels))

    # ax.set_xlim([-0.1, 1])
    ax.set_ylim([0, len(data) + (n_clusters + 1) * 10])

    silhouette_avg = silhouette_score(data, labels, random_state=random_state)
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
        nclusters, autofig, autoax = calc_optimal_clusters(data, end=max_autoclusters,
                                                           random_state=seed)

    kmeans = KMeans(n_clusters=nclusters, random_state=seed).fit(data)

    clusters = pd.Series(data=kmeans.labels_, index=data.index)
    cluster_order = clusters.sort_values().index

    # order the kmeans clusters via hierachical clustering for visual display
    _d = dict()
    for c in clusters.unique():
        sel = data.loc[ clusters[clusters==c].index ]
        _d[c] = sel.mean() # or median, or mean/std, or some variation
    _df = pd.DataFrame(_d).T
    _linkage = hierarchy.linkage(_df, method='ward', optimal_ordering=True)
    _dend = hierarchy.dendrogram(_linkage, no_plot=True)
    _order = _dend['leaves']
    cluster_order = _df.iloc[_order].index
    # cluster_order = _df.index
    clusters_categorical = pd.Series(pd.Categorical(clusters, categories=cluster_order, ordered=True),
                                        index=clusters.index).sort_values()
    cluster_order_remapping = {c: i for i, c in enumerate(cluster_order)}
    clusters_categorical = clusters_categorical.map(cluster_order_remapping)
    clusters = clusters.map(cluster_order_remapping)

    cluster_centers = kmeans.cluster_centers_[ cluster_order ]
    cluster_order_remapping

    suborder = list()
    for c in clusters_categorical.cat.categories:
        ids = clusters_categorical.where(lambda x: x == c).dropna().index

        _linkage = hierarchy.linkage(data.loc[ids], method='ward', optimal_ordering=True)
        _dend = hierarchy.dendrogram(_linkage, no_plot=True)
        _order = _dend['leaves']
        suborder.append(data.loc[ids].iloc[_order].index)

    suborder = np.concatenate(suborder)

    clusters_categorical = clusters_categorical.loc[suborder]


    silhouette_scores = silhouette_samples(data.loc[clusters_categorical.index],
                                           clusters_categorical.values)

    # fig, ax = silhouette_plot(data, kmeans.labels_)
    fig, ax = silhouette_plot(data, clusters.values, random_state=seed)

    # ==============================================================
    # cluster centers plot
    # cluster_centers = kmeans_result['cluster_centers']

    # cmap = iter(rgb2hex(x) for x in sb.color_palette('hls', n_colors=nclusters))
    # cmap_mapping = OrderedDict(((val,  next(cmap)) for val in range(nclusters)))
    cmap_mapping = {val: rgb2hex(c) for val, c in enumerate(sb.color_palette('hls', n_colors=nclusters))}
    print(list(rgb2hex(x) for x in sb.color_palette('hls', nclusters)))
    print('tackle cmap', cmap_mapping)
    cluster_colors = clusters.loc[clusters_categorical.index].map(cmap_mapping).to_frame('Cluster')

    # print(cluster_centers.shape)
    figsize = [max(cluster_centers.shape[1]*.75, max(cluster_centers.shape[0]*.5, 16)), cluster_centers.shape[0]]
    # print(figsize)
    ccenter_fig, ccenter_axs = plt.subplots(nrows=cluster_centers.shape[0], sharex=True,
                                            sharey=True, figsize=figsize)


    ymin, ymax = cluster_centers.min()*1.1, cluster_centers.max()*1.1
    for ax, category in zip(ccenter_axs.flat, clusters_categorical.cat.categories):
        color=cmap_mapping[category]
        ax.plot(cluster_centers[category], color=color, marker='o')
        ax.axhline(y=0, ls='--', color='grey', alpha=.6)
        bbox = dict(boxstyle="round", ec=color, fc=color+'aa', alpha=0.5)
        # ax.set_ylabel(category+1, rotation=0, backgroundcolor=cmap_mapping[category]+'aa', boxstyle='round')
        ax.set_ylabel(category+1, rotation=0, bbox=bbox)
        ax.set_yticklabels([])
        # ax.set_ylim(-ymax, ymax)
        # ax.yaxis.set_major_locator(ticker.LinearLocator(5))
        ax.grid(axis='x', zorder=9)
        ax.patch.set_alpha(0)

    ax.set_xticks(range(len(data.columns)))
    ticks = list()
    for tick in data.columns:
        txt = textwrap.fill(tick.replace('_', ' ', 1), 6, break_long_words=False)
        ticks.append(txt)
    ax.set_xticklabels(ticks, rotation=90)
    ccenter_axs[0].set_title('Cluster Centers')
    ymax = max(np.abs(ax.get_ylim()))*1.14
    ax.set_ylim(-ymax, ymax)
    # ccenter_fig.subplots_adjust(hspace=-.2)
    ccenter_fig.tight_layout()
    sb.despine(fig=ccenter_fig, left=True, bottom=True)


    # ==============================================================

    ret = {'nclusters': nclusters, 'auto': {'fig': autofig, 'ax': autoax},
           'silhouette': {'fig': fig, 'ax': ax},
           'clusters': clusters,
           # 'kmeans': kmeans,
           'cluster_centers': cluster_centers,
           'nclusters': clusters.nunique(),
           'clusters_categorical': clusters_categorical,
           'silhouette_scores': silhouette_scores,
           'cluster_center_plot': dict(fig=ccenter_fig, axs=ccenter_axs)
    }

    return ret



def clusterplot(data, annot_mat=None,
                cmap_name=None, dbscan=False, genes=None, highlight_gids=None, highlight_gid_names=None, gid_symbol=None,
                nclusters=None, gene_symbols=None, z_score=None, z_score_by=None, standard_scale=None, mask=None,
                show_missing_values=True, max_autoclusters=30, row_cluster=True,
                seed=None, col_cluster=True, metadata=None, col_data=None, figsize=None,
                normed=False, linkage='average',
                gene_symbol_fontsize=8, legend_include=None, legend_exclude=None,
                metadata_colors=None, circle_col_markers=False, circle_col_marker_size=12,
                force_optimal_ordering=False,
                square=False, force_plot_genes=False,
                main_title=None, order_by_abundance=False
):
    """
    :nclusters: None, 'auto', or positive integer
    :force_plot_genes: bool, default False. When True, will force display of all genes in `genes` as missing values

    """

    rc = {'font.family': 'sans-serif',
          "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans",
                              "Bitstream Vera Sans", "sans-serif"],
    }

    if order_by_abundance:
        gene_order = data.mean(1).sort_values(ascending=False).index
        data = data.loc[gene_order]
        mask = mask.loc[gene_order]
        annot_mat = annot_mat.loc[gene_order]


    sb.set_context('notebook')
    sb.set_palette('muted')
    sb.set_color_codes()
    sb.set_style('white', rc)

    retval = dict()
    data = data.copy()
    if mask is not None:
        mask = mask.copy()
    if annot_mat is not None:
        annot_mat = annot_mat.copy()

    if genes is not None:  # only plot these select genes
        if not force_plot_genes:
            _genes = [x for x in genes if x in (set(genes) & set(data.index))]
        else:
            _genes = genes
        # data = data.loc[_genes].fillna(0)
        # mask = mask.loc[_genes].fillna(True)
        data = data.reindex(_genes).fillna(0)
        mask = mask.reindex(_genes).fillna(True)
        if annot_mat is not None:
            annot_mat = annot_mat.reindex(_genes).fillna(0).astype(int)

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
    if not legend_exclude:
        legend_exclude = tuple()
    if not legend_include and col_data is not None:
        # legend_include = col_data.index.tolist()
        legend_include = col_data.columns.tolist()


    if col_data is not None:
        # col_data = parse_metadata(metadata)
        # _legend_data = set(legend_include) - set(legend_exclude) if legend_include else col_data.index
        _legend_data = set(legend_include) - set(legend_exclude) if legend_include else col_data
        # _legend_data = sorted(_legend_data)
        _legend_data = [x for x in legend_include if x in _legend_data] # preserve order

        # col_data = col_data.loc[_legend_data].copy()
        col_data = col_data[_legend_data].copy()
        # if 'label' in col_data.index:
        #     if col_data.loc['label'].nunique() == 1:
        #         col_data = col_data.loc[col_data.index.drop('label')]
        if 'label' in col_data:
            if col_data['label'].nunique() == 1:
                # col_data = col_data[col_data.index.drop('label')]
                col_data = col_data.drop('label', axis=1)

        col_colors = pd.DataFrame(columns=col_data.columns,
                                  index=col_data.index)

        # for info in col_data.index:
        for info in col_data.columns:
            # n_colors=max(10, col_data.loc[info].dropna().nunique())
            n_colors=max(10, col_data[info].dropna().nunique())
            if n_colors <= 10:
                colors = sb.color_palette('tab10', n_colors=n_colors)
            else:
                # colors = sb.color_palette('cubehelix', n_colors=n_colors)
                colors = sb.cubehelix_palette(n_colors, start=.1, rot=.75)
            cmap = iter(rgb2hex(x) for x in
                        colors)

            # mapping = {val : next(cmap) if not pd.isna(val) else 'grey'
            #            for val in col_data.loc[info].unique()}
            # colors = col_data.loc[info].map(mapping)
            # col_colors.loc[info] = colors

            mapping = dict()
            for val in sorted(col_data[info].unique()):
                if pd.isna(val):
                    continue
                if val == 'NA':
                    mapping['NA'] = 'grey'
                    continue
                if metadata_colors and info in metadata_colors:
                    c = metadata_colors[info].get(val, next(cmap))
                else:
                    c = next(cmap)
                mapping[val] = c
#             mapping = {val : next(cmap) if not pd.isna(val) else 'grey'
#                        for val in col_data[info].unique()}

            colors = col_data[info].map(mapping) # cannot use np.nan as a dictionary key!
            if colors.isna().any():
                colors.loc[colors.isna()] = 'grey'  # cannot use np.nan as a dictionary key!
            colors.loc['NA'] = 'grey'  # cannot use np.nan as a dictionary key!
            col_colors[info] = colors

        # force NANs grey
        # col_colors.loc[ col_data.loc[col_colors.index].isna() ] = '#444444'
        # col_colors = col_colors.T

    _genemapper = GeneMapper()

    _geneids = data.index.copy()
    if gene_symbols:  # change index to symbols
        assert all(data.index == mask.index)
        if annot_mat is not None:
            assert all(data.index == annot_mat.index)
        clustermap_symbols = [gid_symbol.get(str(x), _genemapper.symbol.get(str(x), str(x)))
                              for x in data.index]
        data.index = clustermap_symbols
        mask.index = clustermap_symbols
        if annot_mat is not None:
            annot_mat.index = clustermap_symbols
        if row_colors is not None:
            row_colors.index = clustermap_symbols

    # if nclusters is not None or dbscan:
    if standard_scale is not None:
        if standard_scale == 1:
            data_t = (data.T / data.max(1)).T
        elif standard_scale == 0:
            data_t = (data / data.max(1))
    elif z_score is not None:
        if z_score_by:
            data_t = pd.concat([
                sb.matrix.ClusterGrid.z_score(data[col_data[col_data[z_score_by]==x].index],
                                                   z_score)
                                              for x in col_data[z_score_by].unique()
            ], sort=False, axis=1)
        else:
            data_t = sb.matrix.ClusterGrid.z_score(data, z_score)
        _minval = data_t.min().min()
        data_t = data_t.fillna(_minval)
        # data[ data < 0] = 0
        # data_t = sb.matrix.ClusterGrid.z_score(data.replace(0, np.NAN), z_score)
        # data_t = data_t.fillna(data_t.min().min())
    else:
        data_t = data


    if nclusters is not None:

        kmeans_result = calc_kmeans(data_t, nclusters, seed, max_autoclusters)
        # kmeans = kmeans_result['kmeans']
        clusters = kmeans_result['clusters']
        nclusters = kmeans_result['nclusters']
        clusters_categorical = kmeans_result['clusters_categorical'] #optimally ordered for visualization

        plot_data = data_t.loc[clusters_categorical.index]

        cmap = iter(rgb2hex(x) for x in sb.color_palette('hls', n_colors=max(6, nclusters)))
        cmap_mapping = {val : next(cmap) for val in range(nclusters)}
        cluster_colors = clusters.loc[clusters_categorical.index].map(cmap_mapping).to_frame('Cluster')

        # this is all done in the calc_kmeans function now
        # # ==============================================================
        # # cluster centers plot
        # cluster_centers = kmeans_result['cluster_centers']
        # print(cluster_centers.shape)
        # figsize = [max(cluster_centers.shape[1]*.75, 16), cluster_centers.shape[0]]
        # print(figsize)
        # ccenter_fig, axs = plt.subplots(nrows=cluster_centers.shape[0], sharex=True, sharey=True, figsize=figsize)
        # for ax, category in zip(axs.flat, clusters_categorical.cat.categories):
        #     color=cmap_mapping[category]
        #     ax.plot(cluster_centers[category], color=color, marker='o')
        #     ax.axhline(y=0, ls='--', color='grey', alpha=.6)
        #     bbox = dict(boxstyle="round", ec=color, fc=color+'aa', alpha=0.5)
        #     # ax.set_ylabel(category+1, rotation=0, backgroundcolor=cmap_mapping[category]+'aa', boxstyle='round')
        #     ax.set_ylabel(category+1, rotation=0, bbox=bbox)
        #     ax.set_yticklabels([])
        # ax.set_xticks(range(len(plot_data.columns)))
        # ax.set_xticklabels(plot_data.columns)
        # axs[0].set_title('Cluster Centers')
        # ymax = max(np.abs(ax.get_ylim()))*1.14
        # ax.set_ylim(-ymax, ymax)
        # ccenter_fig.tight_layout()
        # sb.despine(fig=ccenter_fig, left=True, bottom=True)

        # kmeans_result['cluster_center_plot'] = dict(fig=ccenter_fig)

        # # ==============================================================



        cluster_data = data_t.copy()
        # ?? this is not right
        # if show_missing_values:
            # cluster_data[mask] = np.NAN
            # cluster_data[mask] = 0
        cluster_data = cluster_data.assign(Cluster=clusters+1)
        cluster_data['silhouette_score'] = kmeans_result['silhouette_scores']

        if not gene_symbols:
            cluster_data['GeneSymbol'] = cluster_data.index.map(lambda x: gid_symbol.get(x,
                                                                        _genemapper.symbol.get(x, x)
            ))
        elif gene_symbols:
            cluster_data['GeneID'] = _geneids

        if row_colors is None:
            row_colors = cluster_colors
        else:
            row_colors = pd.concat([cluster_colors, row_colors], axis=1)

        kmeans_result['data'] = cluster_data.loc[clusters_categorical.index]
        if gene_symbols:
            kmeans_result['data'].index.name = 'GeneSymbol'
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

    if cmap_name is None:
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

    FONTSIZE = gene_symbol_fontsize # default 8

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

        # figwidth  = max( min( len(data.columns) / 2, 8), min_figwidth )
        figwidth  = max( min( len(data.columns) / 2, 16), min_figwidth )
        # if col_colors is not None:
        #     figwidth -= (max(len(x) for x in col_colors.columns) * .16667)
        figsize = (figwidth, figheight)
    else:
        figwidth, figheight = figsize


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

    # # TODO this should be done to make unique
    # index_counts = plot_data.index.value_counts().where(lambda x : x > 1).dropna()
    # for ix, entry in index_counts.items():
    #     break

    # simple workaround, fix this:
    if len(set(plot_data.index)) < len(plot_data.index):
        plot_data = plot_data.groupby(plot_data.index).mean()
        mask = mask.groupby(mask.index).mean()
        if order_by_abundance:
            plot_data = plot_data.loc[gene_order]


    if nclusters is not None: # put back in order again
        plot_data = plot_data.loc[clusters_categorical.index]
        mask = mask.loc[clusters_categorical.index]

    # a minimal subclass of seaborn ClusterGrid for scaling.
    plotter = MyClusterGrid(plot_data,
                            figsize=figsize,
                            row_colors=row_colors if row_colors is not None and not row_colors.empty else None,
                            col_colors=col_colors if col_colors is not None and not col_colors.empty else None,
                            mask=mask.loc[plot_data.index] if show_missing_values else None,
                            circle_col_markers=circle_col_markers,
                            circle_col_marker_size=circle_col_marker_size,
                            force_optimal_ordering=force_optimal_ordering,
                            # heatmap_height_ratio=heatmap_height_ratio,
                            # dendrogram_width_ratio=dendrogram_width_ratio,
                            # heatmap_width_ratio=heatmap_width_ratio,
    )

    col_color_kws=dict(fontsize=FONTSIZE)
    if circle_col_markers:
        col_color_kws['linewidths'] = .5
        col_color_kws['linecolor'] = '#111111'


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
                     xticklabels=plot_data.columns,
                     # col_color_kws=dict(fontsize=12),
                     col_color_kws=col_color_kws,
                     annot=annot_mat.loc[plot_data.index] if annot_mat is not None else None,
                     annot_kws=dict(size=8),
                     square=square,


    )

    if annot_mat is not None:
        _txt = 'Annotation : {}'.format(annot_mat.columns.name)
        g.ax_row_dendrogram.text(0, -.14, _txt, transform=g.ax_row_dendrogram.transAxes, fontsize=8)


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
        if len(tick.get_text()) > 10:
            tick.set_size(tick.get_size()*.75)
    if g.ax_row_colors:
        ticks = g.ax_row_colors.xaxis.get_ticklabels()
        if len(ticks) > 1:
            scale = .9/ (np.log1p(len(ticks)) *.9)
            for tick in ticks:
                tick.set_size(tick.get_size()*scale)

    g.ax_heatmap.set_ylabel('') # don't need the GeneID label
    # plt.rcParams['figure.constrained_layout.use'] = True

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
        # col_label_lengths = col_data.fillna('').applymap(len).max(1) + col_colors.nunique()
        # col_label_lengths = col_data.astype(str).applymap(len).max(0) + col_colors.nunique()
        col_label_lengths = col_data.applymap(len).max(0) + col_colors.nunique()
        # widths = _calculate_box_sizes( col_colors.nunique() )
        # widths = _calculate_box_sizes( col_label_lengths, start_pos=-.2, end_pos=1.2 )
        widths = _calculate_box_sizes( col_label_lengths, start_pos=0.0, end_pos=1.1 )
        col_colors_t = col_colors.T

        # bbox_y0 = 1.44 if col_cluster else .8
        # bbox_y0 = 2.16 if col_cluster else 1.2
        bbox_y0 = 2.46 if col_cluster else 1.8

        bboxes = [(x, bbox_y0, 1, .2) for x in widths]  # (x0, y0, width, height)
        # bboxes = [(x, 1.02, 1, .2) for x in np.arange(0, 1, 1/len(col_colors_t.index))]
        legends = list()
        for bbox, ix in zip(bboxes, col_colors_t.index):
            col_name        = ix
            # col_labels      = col_data.loc[ix].drop_duplicates()
            col_labels      = col_data[ix].drop_duplicates()
            col_names       = col_labels.values
            label_colors    = col_colors_t.loc[ix, col_labels.index].values
            handles, labels = list(), list()
            for n, c in zip(col_names, label_colors):
                if circle_col_markers:
                    handle = mpl.lines.Line2D(range(1), range(1), color="none", marker='o', markerfacecolor=c,
                                              markeredgewidth=.5, markeredgecolor='#111111',
                                              markersize=12)
                else:
                    handle = mpl.patches.Patch(color=c,)
                handles.append(handle)
                labels.append(n)
            ncols = max(len(col_names) // 3, 1)
            # if len(col_names) <= 10:
            #     ncols = max(len(col_names) // 3, 1)
            # else:
            #     ncols = max(len(col_names) // 5, 1)
            # ncols = min(ncols, 20)
            leg = g.ax_col_dendrogram.legend( handles, labels, bbox_to_anchor=bbox,
                                              loc='upper left', ncol=ncols,
                                              frameon=False,
                                              title=textwrap.fill(col_name.replace('_',' '),
                                                                  8,
                                                                  break_long_words=False)
            )
            legends.append(leg)
            g.ax_col_dendrogram.add_artist(leg)
            extra_artists = legends


    if figheight <= 12:
        hspace =.01
        wspace = .01
    else:
        hspace = .01 / (22*figheight)
        wspace = .01


    if main_title is not None:
        g.fig.suptitle(main_title)

    g.fig.subplots_adjust(hspace=hspace, wspace=wspace,
                          # left=.5/figwidth, right=1-1./figwidth,
                          left=.04, right=.88,
                          bottom=1/figheight, top=1-(1.4/figheight)
    )


    # make sure there is enough room on the right side for labels
    if col_colors is not None and not col_colors.empty:

        width, height = g.fig.get_size_inches()

        longest_label = max(len(x) for x in col_colors.columns) + 6  # add a little padding
        char_width = (430/1000) # approx via from https://www.math.utah.edu/~beebe/fonts/afm-widths.html
        longest_length = longest_label * char_width
        # inch_shift = longest_length * 12/72  # 72 pts in an inch
        # inch_shift = longest_length * FONTSIZE/72  # 72 pts in an inch

        inch_shift = longest_length * (0.0393701*FONTSIZE/12) # 1 mm, 0.039 inches


        # shift = 1 - (inch_shift / width)
        shift = 1 - (inch_shift)

        g.gs.update(right=shift)  # add some room on the right so everything fits

    g.ax_heatmap.yaxis.set_label('')


    retval['clustermap'] = dict(clustergrid=g, extra_artists=extra_artists)

    return retval
    # return g, extra_artists
