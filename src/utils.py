"""

"""
import os
import re
import configparser
import operator as op
from collections import OrderedDict, defaultdict
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sb
from seaborn.distributions import _freedman_diaconis_bins as seaborn_bin_calc
import click


from bcmproteomics_ext import ispec
sb.set_context('notebook', font_scale=1.8)

idx = pd.IndexSlice

N_COLORS = 100
# r_colors = sb.color_palette("coolwarm", n_colors=N_COLORS+1)
r_colors = sb.color_palette("RdBu_r", n_colors=N_COLORS+1)
STEP = .2


# def filter_observations(panel, column, threshold):
#     """
#     Filter by less than or equal to threshold of 0 observations
#     """
#     indices = (panel.minor_xs(column)
#                .fillna(0)
#                .where(lambda x: x != 0)
#                .count(1)
#                .where(lambda x: x >= threshold)
#                .dropna()
#                .index
#     )
#     return panel.loc[:, indices, :]

def filter_observations(df, column, nonzero_value, subgroup=None, metadata=None):
    """
    format is:
                         Sample1  Sample2  Sample3 ..
    GeneID column_name
    a      PSMs          1        2         3
    a      FunCats       DBTF     DBTF      DBTF
    a      iBAQ_dstrAdj  1        1         NA
    """
    if isinstance(nonzero_value, int):
        threshold = nonzero_value

    if subgroup is not None and metadata is None:
        raise ValueError('Must provide metadata if specifying subgroup')

    if subgroup is None:
        if isinstance(nonzero_value, float): #then ratio of total
            threshold = len(df.columns) * nonzero_value

        mask = (df.loc[ idx[:, column], :].fillna(0)
                .where(lambda x : x != 0)
                .count(1)
                .where(lambda x: x >= threshold)
                .dropna())

        gids = mask.index.get_level_values(0)

        return df.loc[ gids.values ]

    else:

        all_gids = set()

        for sample, grp in metadata.T.groupby(subgroup):
            columns = grp.index

            if isinstance(nonzero_value, float): #then ratio of total
                threshold = len(columns) * nonzero_value

            mask = (df.loc[ idx[:, column], columns].fillna(0)
                    .where(lambda x : x != 0)
                    .count(1)
                    .where(lambda x: x >= threshold)
                    .dropna())

            gids = mask.index.get_level_values(0)
            all_gids |= set(gids)

        return df.loc[ idx[tuple(all_gids), :], : ]



def filter_sra(df, SRA='S'):

    mask = ((df.loc[ idx[:, 'SRA'], :] == 'S')
            .any(1)
            .where(lambda x: x)
            .dropna()
    )
    gids = mask.index.get_level_values(0)

    return df.loc[ gids.values ]

def filter_funcats(df_long, funcats):

    mask = (df_long.loc[ idx[:, 'FunCats'], df_long.columns[0] ].str.contains(funcats)
            .where(lambda x: x)
            .dropna()
    )
    gids = mask.index.get_level_values(0)
    return df_long.loc[gids.values]



# def filter_taxon(panel, taxon=9606):
#     indices = ((panel.minor_xs('TaxonID') == taxon)
#                .any(1)
#                .where(lambda x : x == True)
#                .dropna()
#                .index
#     )
#     return panel.loc[:, indices, :]

def filter_taxon(df, taxon=9606):
    mask = df.loc[ idx[:, ('TaxonID')], : ] == 9606
    gids = mask[mask].dropna().index.get_level_values(0)
    return df.loc[ gids.values ]

def pearson_r(x, y):
    return stats.pearsonr(x, y)[0]

def spearman_r(x, y):
    return stats.spearmanr(x, y)[0]

def color_diag(g):
    for ax in np.diag(g.axes):
        ax.set_facecolor(r_colors[-1])

def hist(x, xmin=None, xmax=None, colors_only=False, **kwargs):
    if colors_only:
        return

    ax = plt.gca()
    if 'color' in kwargs:
        color = kwargs.pop('color')
    if 'bins' in kwargs:
        kwargs.pop('bins')
    if 'edgecolor' in kwargs:
        kwargs.pop('edgecolor')
    X = x[ x>0 ]
    try:
        nbins = seaborn_bin_calc(X)
    except ZeroDivisionError:
        nbins = 10
    # print(nbins)
    ax.hist(X.values, color='k', bins=nbins, edgecolor='none', **kwargs)

    # sb.despine(ax=ax, left=True, bottom=True)
    if xmin and xmax:
        ax.set_xlim((xmin, xmax))


def remove_ticklabels(fig=None, ax=None):
    # idea via seaborn/utils.py :: despine

    # Get references to the axes we want
    if fig is None and ax is None:
        axes = plt.gcf().axes
    elif fig is not None:
        axes = fig.axes
    elif ax is not None:
        axes = [ax]

    for ax_i in axes:
        ax_i.set_xticklabels([])
        ax_i.set_yticklabels([])

def plot_cbar(ax):

    labels = list(reversed(['{:.1f}'.format(x) for x in np.arange(1, -1.1, -STEP)]))
    cmap = mpl.colors.ListedColormap(r_colors)
    cbar = mpl.colorbar.ColorbarBase( ax, cmap=cmap)
    cbar.set_ticks(np.arange(0, 2.1, STEP/2))
    cbar.set_ticklabels(labels)
    ax.set_ylabel('r value')

def make_xaxis(ax, yloc=0, offset=0.05, fmt_str='%1.1f', **props):

    xmin, xmax = ax.get_xlim()
    locs = [loc for loc in ax.xaxis.get_majorticklocs()
            if loc >= xmin and loc <= xmax]
    tickline, = ax.plot(locs, [yloc]*len(locs), linestyle='',
                        marker=mpl.lines.TICKDOWN, **props)
    axline, = ax.plot([xmin, xmax], [yloc, yloc], **props)
    tickline.set_clip_on(False)
    axline.set_clip_on(False)
    for loc in locs:
        ax.text(loc, yloc - offset, fmt_str % loc,
                horizontalalignment='center',
                verticalalignment='top')



def plot_delegator(x, y, stat='pearson', filter_zeros=True,
                   upper_or_lower='upper', colors_only=False,
                   shade_correlation=True, **kwargs):
    if upper_or_lower == 'upper':
        func = annotate_stat
    elif upper_or_lower == 'lower':
        func = scatter

    # x_nonzero = x[ (~x.isnull()) & ~(x.abs() == np.inf) ].index
    # y_nonzero = y[ (~y.isnull()) & ~(y.abs() == np.inf) ].index
    x_nonzero = x[ x > 0 ].index
    y_nonzero = y[ y > 0 ].index

    nonzeros = list(set(x_nonzero) & set(y_nonzero))
    # nonzeros = list(set(x_nonzero) | set(y_nonzero))

    X, Y = x, y

    if filter_zeros:
        X = x.loc[nonzeros]
        Y = y.loc[nonzeros]

    kwargs['alpha'] = .4
    ax = plt.gca()

    if stat == 'pearson':
        r = pearson_r(X,Y)
        text = 'Pearson'
    elif stat == 'spearman':
        r = spearman_r(X,Y)
        text = 'Spearman'
    text = 'n = {:,}\nr = {:.2f}'.format(len(X), r)

    if not shade_correlation:
        pass
    elif np.isnan(r):
        print("Could not calculate r")
    else:
        ax_bg_ix = int(round(r+1, 2) * N_COLORS/2 )  # add 1 to shift from -1 - 1 to 0 - 2 for indexing
        ax_bg = r_colors[ax_bg_ix]
        ax.patch.set_facecolor(ax_bg)
        ax.patch.set_alpha(.5)
        # kwargs['text'] = text
    if colors_only:
        return

    func(X, Y, ax, text=text, **kwargs)


def annotate_stat(x, y, ax, text, **kwargs):

    # textsplit = text.split('\n')
    # maxlen = max(len(x) for x in textsplit)
    # size = 30
    # if maxlen >= 9:
    #     size -= 6
    # if maxlen >= 15:
    #     size -= 6

    size = mpl.rcParams['font.size']
    ax.annotate(text, xy=(0.5, 0.5), xycoords='axes fraction',
                va='center', ha='center', size=size
    )
    # sb.despine(ax=ax, left=True, bottom=True)

def scatter(x, y, ax, xymin, xymax, **kwargs):

    s = 4
    marker = '.'
    if 'text' in kwargs:
        kwargs.pop('text')
    if 'color' in kwargs:
        kwargs.pop('color')
    if 'alpha' in kwargs:
        alpha = kwargs.pop('alpha')
    if 's' in kwargs:
        s = kwargs.pop('s')
    if 'marker' in kwargs:
        marker = kwargs.pop('marker')


    ax.scatter(x, y, color='#222222', alpha=alpha, s=s, marker=marker, **kwargs)
    # sb.despine(ax=ax, left=True, bottom=True)

    if xymin and xymax:
        ax.set_xlim((xymin, xymax))
        ax.set_ylim((xymin, xymax))

    # anchored_text = AnchoredText(text, loc=2)
    # ax.add_artist(anchored_text)


def save_multiple(fig, filename, *exts, verbose=True, dpi=300, **save_kwargs):
    """save a figure to a specific file multiple
    times with different extensions"""
    # rasterized = True
    rasterized = False
    # if 'rasterized' in save_kwargs:
    #     rasterized = save_kwargs.pop('rasterized')

    for ext in exts:
        out = os.path.abspath(filename+ext)

        if ext == '.pdf':
            rasterized=True
        else:
            rasterized=False

        if verbose:
            print("Saving", out, '...', end='', flush=True)
        fig.savefig(out, dpi=dpi, rasterized=rasterized, **save_kwargs)
        if verbose:
            print('done.', flush=True)

def make_config(path='.'):
    config = configparser.ConfigParser()
    config.optionxform = str
    config['Name'] = OrderedDict((('recno', 12345),
                                  ('runno', 1),
                                  ('searchno', 1)))
    file_ = os.path.join(path, 'generic_config.ini')
    with open(file_, 'w') as cf:
        print('Creating ', file_)
        config.write(cf)
        cf.write('#runno and searchno are optional, default to 1\n')

@lru_cache()
def read_config(configfile, enforce=True):
    """reads config file and returns the data.
    Needs to have recno at a minimum.
    recno, runno, searchno are used for getting data.
    Other fields are used for PCA plots and (optionally) clustermaps.
    """
    config = configparser.ConfigParser()
    config.optionxform = str
    with open(configfile, 'r') as f:
        config.read_file(f)

    sections = config.sections()  # retains order
    FIELDS = ('recno', 'runno', 'searchno', 'label')


    data = defaultdict(lambda : dict(runno=1, searchno=1, label='none'))  # does not retain order (no guarantee)
    for section_key in sections:
        section = config[section_key]
        other_fields = set(section.keys()) - set(FIELDS)
        for field in FIELDS:
            value = section.get(field)
            if value is None:
                continue
            data[section_key][field] = value
        for field in other_fields:
            value = section.get(field)
            data[section_key][field] = value
        if section_key.startswith('__'):
            pass
        elif 'recno' not in data[section_key] and enforce:  # record number is required
            print(section_key, 'does not have recno defined, skipping')
            data.pop(section_key)

    ordered_data = OrderedDict()
    for key in sections:
        if key not in data.keys():
            continue
        ordered_data[key] = data[key]

    return ordered_data

def parse_gid_file(gids):

    gid_out = list()
    with open(gids, 'r') as f:
        for line in f:
            try:
                gid = float(line.strip())
                gid_out.append(gid)
            except ValueError:
                pass
    return set(gid_out)

def get_file_name(full_file):
    fname, ext = os.path.splitext(full_file)
    grp = re.search('\w+', fname)
    if grp:
        return grp.group()
    else:
        return None


# def fillna_meta(panel, col):
#     df = panel.loc[:, :, col]
#     panel.loc[:, :, col] = (df.fillna(method='ffill', axis=1)
#                             .fillna(method='bfill', axis=1)
#     )

def fillna_meta(df, index_col):
    """
    Fill NANs across rows
    """
    if index_col not in df.index.get_level_values(1).unique():
        return
    selection = df.loc[idx[:, index_col], :]
    df.loc[idx[:, index_col], :] = (selection.fillna(method='ffill', axis=1)
                                    .fillna(method='bfill', axis=1)
    )



def parse_metadata(metadata):
    # expids = ('recno', 'runno', 'searchno')
    expids = tuple()
    metadata_filtered = OrderedDict([(k,v) for k, v in metadata.items()
                                     if not k.startswith('__')
    ])
    # col_data = pd.DataFrame.from_dict(metadata, orient='columns').filter(regex='^(?!__)')
    # col_data = pd.DataFrame.from_dict(metadata_filtered, orient='columns')
    col_data = pd.DataFrame(metadata_filtered, columns=metadata_filtered.keys())
    col_data = col_data.loc[[x for x in col_data.index if x not in expids]]
    return col_data


class iFOT:
    def __init__(self):
        self.file = os.path.join( os.path.split(os.path.abspath(__file__))[0],
                                  'data', 'geneignore.txt'
        )
        self._to_ignore = None
    @property
    def to_ignore(self):
        if self._to_ignore is None:
            self._to_ignore = pd.read_table(self.file, comment='#', header=None, names=['GeneID'])['GeneID']
            # a pandas.Series
        return self._to_ignore
    def filter(self, genes):
        return [g for g in genes if (self.bool_real(g) and g not in self.to_ignore)]
    @staticmethod
    def bool_real(x):
        flag = True
        try:
            flag = np.isfinite(x)
        except TypeError:
            pass
        return flag and bool(x)
ifot_normalizer = iFOT()

UNANNOTATED_TIDS = (6239,)

def normalize(df, name='name', ifot=False, ifot_ki=False, ifot_tf=False, median=False,
              outcol=None, taxon=None):

    if ifot: # normalize by ifot but without keratins
        norm_ = df.loc[ifot_normalizer.filter(df.index), 'iBAQ_dstrAdj'].sum()
    elif median:
        nonzero_ix = df.query('iBAQ_dstrAdj > 0').index
        norm_ = df.loc[ifot_normalizer.filter(nonzero_ix), 'iBAQ_dstrAdj'].median()
    elif ifot_ki:
        if taxon and taxon in UNANNOTATED_TIDS:
            norm_ = 1
        else:
            norm_ = df.loc[df['FunCats'].fillna('').str.contains('KI'), 'iBAQ_dstrAdj'].sum()
    elif ifot_tf:
        if taxon and taxon in UNANNOTATED_TIDS:
            norm_ = 1
        else:
            norm_ = df.loc[df['FunCats'].fillna('').str.contains('TF'), 'iBAQ_dstrAdj'].sum()
    else:
        norm_ = 1
    if norm_ == 0:
        # error = '{} has a sum of 0 when trying to normalize, aborting'.format(name)
        error = '{} has a sum of 0 when trying to normalize, skipping normalization'.format(name)
        print(error)
        print(kwargs)
        # raise click.Abort()
        sum_ = 1
    if outcol is None:
        outcol = 'area'
    return df['iBAQ_dstrAdj'] / norm_
    # df[outcol] = df['iBAQ_dstrAdj'] / norm_  # use generic 'area' name for all normalization procedures
    # return df

def genefilter(df, funcats=None, funcats_inverse=None, geneid_subset=None, ignored_geneid_subset=None):

    if funcats:  # do this after possible normalization
        df = df[df['FunCats'].fillna('').str.contains(funcats, case=False)]
    if funcats_inverse:  # do this after possible normalization
        df = df[~df['FunCats'].fillna('').str.contains(funcats_inverse, case=False)]
    if geneid_subset:  # do this at the end
        df = df.loc[geneid_subset]
    if ignored_geneid_subset:
        tokeep = set(df.index) - set(ignored_geneid_subset)
        df = df.loc[tokeep]
    valid_ixs = (x for x in df.index if not np.isnan(x))
    return df.loc[valid_ixs]

def filter_and_assign(df, name, funcats=None, funcats_inverse=None, geneid_subset=None,
                      ignored_geneid_subset=None, ifot=False, ifot_ki=False, ifot_tf=False, median=False):
    """Filter by funcats and geneid_subset (if given)
       remove NAN GeneIDs"""

    if ifot: # normalize by ifot but without keratins
        norm_ = df.loc[ifot_normalizer.filter(df.index), 'iBAQ_dstrAdj'].sum()
    elif median:
        norm_ = df.loc[ifot_normalizer.filter(df.index), 'iBAQ_dstrAdj'].median()
    elif ifot_ki:
        norm_ = df.loc[df['FunCats'].fillna('').str.contains('KI'), 'iBAQ_dstrAdj'].sum()
    elif ifot_tf:
        norm_ = df.loc[df['FunCats'].fillna('').str.contains('TF'), 'iBAQ_dstrAdj'].sum()
    else:
        norm_ = 1
    if norm_ == 0:
        error = '{} has a sum of 0 when trying to normalize, aborting'.format(name)
        print(error)
        raise click.Abort()
        # sum_ = 1
    df['area'] = df['iBAQ_dstrAdj'] / norm_  # use generic 'area' name for all normalization procedures

    if funcats:  # do this after possible normalization
        df = df[df['FunCats'].fillna('').str.contains(funcats, case=False)]
    if funcats_inverse:  # do this after possible normalization
        df = df[~df['FunCats'].fillna('').str.contains(funcats_inverse, case=False)]
    if geneid_subset:  # do this at the end
        df = df.loc[geneid_subset]
    if ignored_geneid_subset:
        tokeep = set(df.index) - set(ignored_geneid_subset)
        df = df.loc[tokeep]
    valid_ixs = (x for x in df.index if not np.isnan(x))
    return df.loc[valid_ixs]

def assign_cols(df, name):
    """Filter by funcats and geneid_subset (if given)
       remove NAN GeneIDs"""
    # if funcats:
    #     df = df[df['FunCats'].fillna('').str.contains(funcats, case=False)]
    # if geneid_subset:
    #     df = df.loc[geneid_subset]
    valid_ixs = (x for x in df.index if not np.isnan(x))
    return df.loc[valid_ixs]

def get_outname(plottype: str, name, taxon, non_zeros, colors_only, batch=None,
                batch_method='parametric',
                outpath='.', **kwargs):
    colors = 'colors_only' if colors_only else 'annotated'

    kwarg_values = list()
    for key, value in kwargs.items():
        s = '{}_{}'.format(key, value)
        kwarg_values.append(s)
    kwarg_string = '_'.join(kwarg_values) if kwarg_values else ''

    batch_str = '{}Batch_{}_'.format(batch_method, batch) if batch else ''

    '{}'.format(kwargs)
    fname = '{}_{}_{}_{}_{}more_nonzero_{}{}'.format(name, plottype, taxon, colors, non_zeros,
                                                     batch_str,
                                                     kwarg_string).strip('_')
    return os.path.join(outpath, fname)

class TooManyCategories(Exception):
    pass

from collections import OrderedDict, Callable

class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))
