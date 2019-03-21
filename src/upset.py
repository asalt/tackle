"""
"Fork" from pyupset UpSetPlot to produce better figure
"""
import os
from glob import glob
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pyupset
from pyupset.visualisation import plot as pu_plot
from pyupset.visualisation import UpSetPlot
import seaborn as sb

sb.set_context('talk')
sb.set_style('ticks')



from matplotlib.patches import Rectangle, Circle
class MyUpSetPlot(UpSetPlot):

    def __init__(self, *args, bound_min=1, h_ratio=4, v_ratio=1, dot_size=None, figsize=None, **kwargs):
        self.bound_min = bound_min
        self.h_ratio = h_ratio
        self.v_ratio = v_ratio
        self.dot_size = dot_size
        self.figsize = figsize
        super().__init__(*args, **kwargs)

    def _base_sets_plot(self, sorted_sets, sorted_set_names):
        """
        Plots horizontal bar plot for base set sizes.

        :param sorted_sets: list of data frames, sorted according to user's directives.
        :param sorted_set_names: list of names for the data frames.
        :return: Axes.
        """
        ax = self.ax_setsize
        ax.invert_xaxis()
        height = .6
        # bar_bottoms = self.y_values - height / 2
        bar_bottoms = self.y_values

        ax.barh(bar_bottoms, [len(x) for x in sorted_sets], height=height, color=self.greys[1])

        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 4))

        self._strip_axes(ax, keep_spines=['bottom'], keep_ticklabels=['bottom'])

        ax.set_ylim((height / 2, ax.get_ylim()[1] + height / 2))
        xlim = ax.get_xlim()
        gap = max(xlim) / 500.0 * 20
        ax.set_xlim(xlim[0] + gap, xlim[1] - gap)
        xlim = ax.get_xlim()
        ax.spines['bottom'].set_bounds(xlim[0], xlim[1])

        # ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(3))
        # ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(3))
        # ax.xaxis.set_major_locator(mpl.ticker.AutoLocator(integer=True, min_n_ticks=3, prune='both'))
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, nbins=4))

        # bracket_height = ax.transData.inverted().transform([(0, 0), (0, ax.get_ylim()[1])])
        # bracket_height = np.abs(bracket_height[1, 1] - bracket_height[0, 1])
        # for i, (x, y) in enumerate(zip([len(x) for x in sorted_sets], self.y_values)):
        #     ax.annotate(sorted_set_names[i], rotation=90, ha='right', va='bottom', fontsize=15,
        #                     xy=(x, y), xycoords='data',
        #                     xytext=(-30, 0), textcoords='offset points',
        #                     arrowprops=dict(arrowstyle="-[, widthB=%s"%(bracket_height,),
        #                                     shrinkA=1,
        #                                     shrinkB=3,
        #                                     connectionstyle='arc,angleA=-180, angleB=180, armB=30',
        #                                     ),
        #                     )

        ax.set_xlabel("Set size", fontweight='bold', fontsize=13)

        return ax.get_ylim()


    def _inters_sizes_plot(self, ordered_in_sets, inters_sizes):
        """
        Plots bar plot for intersection sizes.
        :param ordered_in_sets: array of tuples. Each tuple represents an intersection. The array is sorted according
        to the user's directives

        :param inters_sizes: array of ints. Sorted, likewise.

        :return: Axes
        """
        ax = self.ax_intbars
        width = .5
        self._strip_axes(ax, keep_spines=['left'], keep_ticklabels=['left'])

        # bar_bottom_left = self.x_values - width / 2
        bar_bottom_left = self.x_values

        bar_colors = [self._color_for_query(frozenset(inter)) for inter in ordered_in_sets]

        ax.bar(bar_bottom_left, inters_sizes, width=width, color=bar_colors, linewidth=0)
        ax.set_ylim(bottom=0)
        ylim = ax.get_ylim()
        label_vertical_gap = (ylim[1] - ylim[0]) / 60

        for x, y in zip(self.x_values, inters_sizes):
            ax.text(x, y + label_vertical_gap, "%.4g" % y,
                    # rotation=90, ha='center', va='bottom')
                    rotation=0, ha='center', va='bottom')

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 4))

        gap = max(ylim) / 500.0 * 20
        ax.set_ylim(ylim[0] - gap, ylim[1] + gap)
        ylim = ax.get_ylim()
        ax.spines['left'].set_bounds(ylim[0], ylim[1])

        ax.yaxis.grid(True, lw=.25, color='grey', ls=':')
        ax.set_axisbelow(True)
        ax.set_ylabel("Intersection size", labelpad=6, fontweight='bold', fontsize=13)

        return ax.get_xlim()


    def _inters_matrix(self, ordered_in_sets, ordered_out_sets, xlims, ylims, set_row_map):
        """
        Plots intersection matrix.

        :param ordered_in_sets: Array of tuples representing sets included in an intersection. Sorted according to
        the user's directives.

        :param ordered_out_sets: Array of tuples representing sets excluded from an intersection. Sorted likewise.

        :param xlims: tuple. x limits for the intersection matrix plot.

        :param ylims: tuple. y limits for the intersection matrix plot.

        :param set_row_map: dict. Maps data frames (base sets) names to a row of the intersection matrix

        :return: Axes
        """
        ax = self.ax_intmatrix
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        if len(self.x_values) > 1:
            row_width = self.x_values[1] - self.x_values[0]
        else:
            row_width = self.x_values[0]

        self._strip_axes(ax)

        background = plt.cm.Greys([.09])[0]
        if not self.dot_size:
            dot_size = min(300, max(300 - (50*(len(self.x_values)-4)), 50))
        else:
            dot_size = self.dot_size

        for r, y in enumerate(self.y_values):
            if r % 2 == 0:
                ax.add_patch(Rectangle((xlims[0], y - row_width / 2), height=row_width,
                                       width=xlims[1],
                                       color=background, zorder=0))

        for col_num, (in_sets, out_sets) in enumerate(zip(ordered_in_sets, ordered_out_sets)):
            in_y = [set_row_map[s] for s in in_sets]
            out_y = [set_row_map[s] for s in out_sets]
            # in_circles = [Circle((self.x_values[col_num], y), radius=dot_size, color=self.greys[1]) for y in in_y]
            # out_circles = [Circle((self.x_values[col_num], y), radius=dot_size, color=self.greys[0]) for y in out_y]
            # for c in chain.from_iterable([in_circles, out_circles]):
            # ax.add_patch(c)
            ax.scatter(np.repeat(self.x_values[col_num], len(in_y)), in_y,
                       color=np.tile(self._color_for_query(frozenset(in_sets)), (len(in_y), 1)),s=dot_size)
            ax.scatter(np.repeat(self.x_values[col_num], len(out_y)), out_y, color=self.greys[0], s=dot_size)
            ax.vlines(self.x_values[col_num], min(in_y), max(in_y), lw=3.5, color=self._color_for_query(frozenset(in_sets)))


    def _prepare_figure(self, additional_plots):
        """
        Prepares the figure, axes (and their grid) taking into account the additional plots.

        :param additional_plots: list of dictionaries as specified in plot()
        :return: references to the newly created figure and axes
        """
        # fig = plt.figure(figsize=(17, 10.5))
        if self.figsize is not None:
            h, w = self.figsize
        else:
            h, w = 12, 10.5
        fig = plt.figure(figsize=(h, w))
        if additional_plots:
            main_gs = gridspec.GridSpec(3, 1, hspace=.4)
            # main_gs = gridspec.GridSpec(3, 1, hspace=.4, height_ratios=(1, 3, 1))
            topgs = main_gs[:2, 0]
            botgs = main_gs[2, 0]
        else:
            topgs = gridspec.GridSpec(1, 1)[0, 0]
        fig_cols = self.cols + 5
        # fig_rows = self.rows + self.rows * 4
        fig_rows = self.rows + self.rows * self.h_ratio

        # gs_top = gridspec.GridSpecFromSubplotSpec(fig_rows, fig_cols, subplot_spec=topgs, wspace=.1, hspace=.2)

        # setsize_w, setsize_h = 3, self.rows
        setsize_w, setsize_h = int(self.v_ratio*3), self.rows
        tablesize_w, tablesize_h = setsize_w + 2, self.rows
        # tablesize_w, tablesize_h = setsize_w + 1, self.rows
        intmatrix_w, intmatrix_h = tablesize_w + self.cols, self.rows
        # intbars_w, intbars_h = tablesize_w + self.cols, self.rows * 4
        intbars_w, intbars_h = tablesize_w + self.cols, self.rows * self.h_ratio
        # print(tablesize_w, intmatrix_w, intbars_w)
        # print(tablesize_h, intmatrix_h, intbars_h)
        gs_top = gridspec.GridSpecFromSubplotSpec(fig_rows, fig_cols, subplot_spec=topgs, wspace=.1, hspace=.2,
                                                  # height_ratios = [*[1]*32, *[5]*4, *[1]*4]
                                                  # height_ratios = [*[1]*intbars_h, *[4]*(intmatrix_h//2), *[1]*(tablesize_h//2)]
                                                  )
        ax_setsize = plt.subplot(gs_top[-1:-setsize_h, 0:setsize_w])
        ax_tablenames = plt.subplot(gs_top[-1:-tablesize_h, setsize_w:tablesize_w])
        ax_intmatrix = plt.subplot(gs_top[-1:-intmatrix_h, tablesize_w:intmatrix_w])
        # ax_intbars = plt.subplot(gs_top[:self.rows * 4 - 1, tablesize_w:intbars_w])
        ax_intbars = plt.subplot(gs_top[:self.rows * self.h_ratio - 1, tablesize_w:intbars_w])

        add_ax = []
        if additional_plots:
            num_plots = len(additional_plots)
            num_bot_rows, num_bot_cols = int(np.ceil(num_plots / 2)), 2
            # num_bot_rows, num_bot_cols = int(np.ceil(num_plots / 2)), int(np.floor(num_plots / 2 )) + 1
            gs_bottom = gridspec.GridSpecFromSubplotSpec(num_bot_rows, num_bot_cols,
                                                         subplot_spec=botgs, wspace=.15, hspace=.2)
            from itertools import product

            for r, c in product(range(num_bot_rows), range(num_bot_cols)):
               if r+c+1>num_plots: break
               new_plot = plt.subplot(gs_bottom[r, c])
               add_ax.append(new_plotL)

        return fig, ax_intbars, ax_intmatrix, ax_setsize, ax_tablenames, tuple(add_ax)

# from pyupset.visualization
def strip_axes(ax, keep_spines=None, keep_ticklabels=None):
    """
    Removes spines and tick labels from ax, except those specified by the user.

    :param ax: Axes on which to operate.
    :param keep_spines: Names of spines to keep.
    :param keep_ticklabels: Names of tick labels to keep.

    Possible names are 'left'|'right'|'top'|'bottom'.
    """
    tick_params_dict = {'which': 'both',
                        'bottom': 'off',
                        'top': 'off',
                        'left': 'off',
                        'right': 'off',
                        'labelbottom': 'off',
                        'labeltop': 'off',
                        'labelleft': 'off',
                        'labelright': 'off'}
    if keep_ticklabels is None:
        keep_ticklabels = []
    if keep_spines is None:
        keep_spines = []
    lab_keys = [(k, "".join(["label", k])) for k in keep_ticklabels]
    for k in lab_keys:
        tick_params_dict[k[0]] = 'on'
        tick_params_dict[k[1]] = 'on'
    ax.tick_params(**tick_params_dict)
    for sname, spine in ax.spines.items():
        if sname not in keep_spines:
            spine.set_visible(False)

# pyupset.visualisation.UpSetPlot = MyUpSetPlot
greys = plt.cm.Greys([.22, .8])
from pyupset.visualisation import DataExtractor
def plot(data_dict, *, unique_keys=None, sort_by='size', inters_size_bounds=(0, np.inf),
         inters_degree_bounds=(1, np.inf), additional_plots=None, query=None, h_ratio=4, v_ratio=1,
         figsize=(12, 10.5), dot_size=None):
    """
    Plots a main set of graph showing intersection size, intersection matrix and the size of base sets. If given,
    additional plots are placed below the main graph.

    :param data_dict: dictionary like {data_frame_name: data_frame}

    :param unique_keys: list. Specifies the names of the columns that, together, can uniquely identify a row. If left
    empty, pyUpSet will try to use all common columns in the data frames and may possibly raise an exception (no
    common columns) or produce unexpected results (columns in different data frames with same name but different
    meanings/data).

    :param sort_by: 'size' or 'degree'. The order in which to sort the intersection bar chart and matrix in the main
    graph

    :param inters_size_bounds: tuple. Specifies the size limits of the intersections that will be displayed.
    Intersections (and relative data) whose size is outside the interval will not be plotted. Defaults to (0, np.inf).

    :param inters_degree_bounds: tuple. Specified the degree limits of the intersections that will be displayed.
    Intersections (and relative data) whose degree is outside the interval will not be plotted. Defaults to (0, np.inf).

    :param additional_plots: list of dictionaries. See below for details.

    :param query: list of tuples. See below for details.

    :return: dictionary of matplotlib objects, namely the figure and the axes.

    :raise ValueError: if no unique_keys are specified and the data frames have no common column names.

    The syntax to specify additional plots follows the signature of the corresponding matplotlib method in an Axes
    class. For each additional plot one specifies a dictionary with the kind of plot, the columns name to retrieve
    relevant data and the kwargs to pass to the plot function, as in `{'kind':'scatter', 'data':{'x':'col_1',
    'y':'col_2'}, 'kwargs':{'s':50}}`.

    Currently supported additional plots: scatter.

    It is also possible to highlight intersections. This is done through the `query` argument, where the
    intersections to highligh must be specified with the names used as keys in the data_dict.

    """
    query = [] if query is None else query
    ap = [] if additional_plots is None else additional_plots
    all_columns = unique_keys if unique_keys is not None else __get_all_common_columns(data_dict)
    all_columns = list(all_columns)

    plot_data = DataExtractor(data_dict, all_columns)
    ordered_inters_sizes, ordered_in_sets, ordered_out_sets = \
        plot_data.get_filtered_intersections(sort_by,inters_size_bounds,inters_degree_bounds)
    ordered_dfs, ordered_df_names = plot_data.ordered_dfs, plot_data.ordered_df_names

    upset = MyUpSetPlot(len(ordered_dfs), len(ordered_in_sets), additional_plots, query, h_ratio=h_ratio, v_ratio=v_ratio, dot_size=dot_size, figsize=figsize)
    fig_dict = upset.main_plot(ordered_dfs, ordered_df_names, ordered_in_sets, ordered_out_sets,
                               ordered_inters_sizes)
    fig_dict['additional'] = []


    # ap = [{kind:'', data:{x:'', y:''}, s:'', ..., kwargs:''}]
    for i, graph_settings in enumerate(ap):
        plot_kind = graph_settings.pop('kind')
        data_vars = graph_settings.pop('data_quantities')
        graph_properties = graph_settings.get('graph_properties', {})
        data_values = plot_data.extract_data_for(data_vars, query)
        ax = upset.additional_plot(i, plot_kind, data_values, graph_properties, labels=data_vars)
        fig_dict['additional'].append(ax)

    return fig_dict


def make_plot(data, unique_keys, suptitle='', bound_min=1, h_ratio=4, v_ratio=1, query=None, dot_size=None, figsize=(12, 10.5)):

    # res = pyupset.plot({k: data[k]['sig'].to_frame() for k in data.keys()},
    res = plot(data,
               inters_size_bounds=(bound_min, np.inf),
               unique_keys=unique_keys,
               h_ratio=h_ratio,
               v_ratio=v_ratio,
               # additional_plots=[{'kind':'hist', 'data_quantities': {'x': 'GeneID'}}]
               # additional_plots=[{'kind':'scatter', 'data_quantities': {'x': 'GeneID', 'y':'GeneID'}}]
               query=query,
               dot_size=dot_size, figsize=figsize
    )
    fig = res['figure']
    bot_ax = res.get('additional')
    if bot_ax:
        bot_ax = bot_ax[0]
        # bot_ax = res['additional'][0]
        bot_ax.cla()

    fig.subplots_adjust(top=.92)

    # total_data = np.array(sorted( ((k, data[k]['total']) for k in data.keys() )) )

    if bot_ax:
        x, y = total_data[:, 0], list(map( int, total_data[:, 1] ))
        bot_ax.bar(x, y, color=greys[1], width=.5)

        ylim, xlim = (0, max(y)*1.1), (-.5, len(x))
        bot_ax.set_ylim(ylim[0], ylim[1])
        bot_ax.set_xlim(xlim[0], xlim[1])

        bot_ax.yaxis.set_major_locator(mpl.ticker.LinearLocator(5))
        bot_ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))


        bot_ax.spines['left'].set_bounds(ylim[0], ylim[1])
        bot_ax.spines['bottom'].set_bounds(xlim[0], xlim[1])

        bot_ax.yaxis.grid(True, lw=.25, color='grey', ls=':')
        bot_ax.set_xlabel('Number of nonzeros')
        bot_ax.set_ylabel('Total GeneIDs')

    if suptitle:
        fig.suptitle(suptitle)

    return res
