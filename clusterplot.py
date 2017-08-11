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


def clusterplot(ibaqs_log_shifted, highlight_gids=None, highlight_gid_names=None,
                gid_symbol=None, gene_symbols=None, z_score=None, standard_scale=None, row_cluster=True,
                col_cluster=True, metadata=None, col_data=None):
    row_colors = None
    extra_artists = None
    if highlight_gids:

        cmap = sb.color_palette('hls', n_colors=max(6, len(highlight_gids)))
        colors_dfs = list()
        for ix, hgid in enumerate(highlight_gids):
            color = cmap[ix]
            highlights = {gid: color for gid in hgid}
            colors = [highlights.get(x, (1., 1., 1.)) for x in ibaqs_log_shifted.index]
            # colors = ibaqs_log_shifted.index.map( lambda x: highlights.get(x, 'white') )
            colors_df = pd.Series(colors, index=ibaqs_log_shifted.index).to_frame(highlight_gid_names[ix])
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
        clustermap_symbols = [gid_symbol.get(x, '?') for x in ibaqs_log_shifted.index]
        ibaqs_log_shifted.index = clustermap_symbols
        if row_colors is not None:
            row_colors.index = clustermap_symbols

    # figheight = min(len(ibaqs_log_shifted) / 6, 100)
    figheight = 12
    figwidth  = len(ibaqs_log_shifted.columns) / 2
    # if col_colors is not None:
    #     figwidth -= (max(len(x) for x in col_colors.columns) * .16667)
    g = sb.clustermap(ibaqs_log_shifted,
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
    return g, extra_artists
