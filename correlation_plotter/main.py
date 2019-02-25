"""

"""
import sys
import os
import re
import json
import glob
from datetime import datetime
import operator as op
from collections import OrderedDict
from functools import partial
import copy
from tempfile import NamedTemporaryFile
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sb
from seaborn.distributions import _freedman_diaconis_bins
import click
# try:
#     from quick import gui_option
# except ModuleNotFoundError:
#     pass

# rc = {'font.family': 'serif',
#       'font.serif': ['Times', 'Palatino', 'serif']}
# sb.set_context('paper')
# sb.set_style('white', rc)

rc = {'font.family': 'sans-serif',}
sb.set_context('notebook')
sb.set_style('white', rc)
sb.set_palette('muted')
sb.set_color_codes()

__version__ = '0.39'

from bcmproteomics_ext import ispec
sb.set_context('notebook', font_scale=1.4)

from .scatterplot import scatterplot
from .clusterplot import clusterplot
from .metrics import make_metrics
from .pcaplot import pcaplot
from .utils import *
from .containers import Data
from .barplot import barplot
# from cluster_to_plotly import cluster_to_plotly

sys.setrecursionlimit(10000)

TAXON_MAPPER = {'human': 9606,
                'mouse': 10090}

from contextlib import contextmanager
# workaround for windows
@contextmanager
def named_temp(*args, **kwargs):
    f = NamedTemporaryFile(*args, delete=False, **kwargs)
    try:
        yield f
    finally:
        try:
            os.unlink(f.name)
        except OSError:
            pass



def run(data_obj):


    # if 'scatter' in plots or 'all' in data_obj.plots:
    if data_obj.make_plot('scatter'):
        g = scatterplot(data_obj.areas_log_shifted, stat=data_obj.stat,
                        colors_only=data_obj.colors_only,
                        shade_correlation=data_obj.shade_correlation)
        outname = get_outname('scatter', name=data_obj.outpath_name, taxon=data_obj.taxon, non_zeros=data_obj.non_zeros,
                              colors_only=data_obj.colors_only, outpath=data_obj.outpath)
        save_multiple(g, outname, '.png', dpi=96)

    # ibaqs_zscore = (ibaqs_log_shifted - ibaqs_log_shifted.mean(axis=1)) / ibaqs_log_shifted.std(axis=1)
    # ibaqs_log_shifted
    # if 'cluster' in plots or 'all' in plots:
    if data_obj.make_plot('cluster'):
        g, extra_artists = clusterplot(data_obj.areas_log_shifted, highlight_gids=data_obj.highlight_gids,
                                       highlight_gid_names=data_obj.highlight_gid_names,
                                       gid_symbol=data_obj.gid_symbol,
                                       gene_symbols=data_obj.gene_symbols, z_score=data_obj.z_score,
                                       standard_scale=data_obj.standard_scale,
                                       row_cluster=data_obj.row_cluster, col_cluster=data_obj.col_cluster,
                                       metadata=data_obj.config if data_obj.show_metadata else None,
                                       col_data = data_obj.col_metadata
        )
        outname = get_outname('clustermap', name=data_obj.outpath_name, taxon=data_obj.taxon,
                              non_zeros=data_obj.non_zeros,
                              colors_only=data_obj.colors_only, outpath=data_obj.outpath)
        # print(outname)

        bbox_inches='tight'
        if extra_artists is not None:
            bbox_inches=None
        # extra_artists=None
        save_multiple(g, outname, '.png', bbox_extra_artists=extra_artists, bbox_inches=bbox_inches)

    # if 'pca' in plots or 'all' in plots:
    if data_obj.make_plot('pca'):
        fig, ax = pcaplot(data_obj.areas_log_shifted, data_obj.config, col_data = data_obj.col_metadata)
        outname = get_outname('pcaplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
                              non_zeros=data_obj.non_zeros,
                              colors_only=data_obj.colors_only, outpath=data_obj.outpath,
        )
        save_multiple(fig, outname, '.png')

# def file_or_subcommand(ctx, param, value):
#     pass
class Path_or_Subcommand(click.Path):

    def __init__(self, *args, **kwargs ):
        super(Path_or_Subcommand, self).__init__(*args, **kwargs)

    def convert(self, value, param, ctx):

        commands = ctx.command.commands.keys()

        if value in commands:
            help_txt = globals()[value].get_help(ctx)
            click.echo(help_txt)
            sys.exit(0)

        return click.Path.convert(self, value, param, ctx)

class int_or_ratio(click.ParamType):
    name = 'integer'

    def convert(self, value, param, ctx):
        try:
            float(value)
        except ValueError:
            self.fail('%s is not a valid integer or float' % value, param, ctx)
        try:
            return int(value)
        except ValueError: # then float
            val = float(value)
            if val > 1:
                self.fail('%s cannot be greater than 1 if specified as a float' % value, param, ctx)
            return val


def validate_cluster_number(ctx, param, value):

    if value == 'auto':
        return 'auto'

    elif value == 'None' or value is None:
        return None

    elif value.isdigit():
        retval = int(value)
        if retval < 1:
            raise click.BadParameter('must be set to at least 1')
        return retval

    else:
        raise click.BadParameter('must be one of `None`, `auto` or an integer')

def validate_in_config(*entries, valid_entries=None):
    import difflib
    correct = list()
    for entry in entries:
        if entry not in  valid_entries:
            closest_match = difflib.get_close_matches(entry, valid_entries, n=1, cutoff=.6)
            if closest_match:
                response2 = 'Did you  mean {}'.format(closest_match[0])
            else:
                response2 = ""
            click.echo("""{} is not in the config file. {}
            """.format(entry, response2))
        else:
            correct.append(entry)
    return correct




def validate_seed(ctx, param, value):

    if value == 'None' or value is None:
        return None

    elif value.isdigit():
        return int(value)

    else:
        raise click.BadParameter('Must be an integer or `None`')

# def validate_configfile(experiment_file, nonzero_subgroup=None, batch=None, group=None, covariate=None, pairs=None):
def validate_configfile(experiment_file, **kwargs):
    """
    check to ensure all metadata is present in all samples
    :experiment_file: config file to validate metadata entries
    :**kwargs: metadata entries to test
    """
    config = read_config(experiment_file)
    config = copy.deepcopy(config)  # because we're about to change it

    config_len = len([x for x in config.keys()
                      if not x.startswith('__')])

    dunder_fields = dict(__PCA__   = ('color', 'marker'),
                         __norm__  = ('label', 'group'),
                         __batch__ = ('batch', ),
    )

    for entry in config.keys():
        if entry.startswith('__'):
            dunder = config.get(entry)
            fields = dunder_fields.get(entry)
            if fields is None:
                continue
            # validate each field
            for field in fields:
                if field not in config[entry]: continue #
                count = 0
                value = config[entry][field]  # value to tally

                for k, v in config.items():
                    if k in dunder_fields:
                        continue
                    ret = v.get(value)
                    if ret is not None:
                        count += 1

                if count == 0:
                    raise click.BadParameter("""Invalid entry for {}:
                    {} is not annotated in the config file""".format(entry, value))
                if count < config_len:
                    raise click.BadParameter("""Invalid entry for {}:
                    {} is not annotated for all experiments in the config file
                    ({} / {} )""".format(entry, value, count, config_len
                    ))

    if '__norm__' in config.keys():
        # need to do additional validation
        norm_info = config.get('__norm__')
        control = norm_info['control']
        group   = norm_info['group']
        label   = norm_info['label']
        groups = set()
        counter = 0
        for k, v in config.items():
            if k in dunder_fields:
                continue
            subgroup = v.get(group)
            groups.add(subgroup)
            ret = v.get(label)
            if ret == control:
                counter += 1
        if counter != len(groups):
            raise click.BadParameter("""
            Invalid specification for `control` in the __norm__ specification.
            Expected {} to be specified by {} a total of {} times
            but is observed {} time(s).
            """.format(control, label, len(groups), counter))

    def check_group(name, name_str):
        count = 0
        for k, v in config.items():
            if k in dunder_fields:
                continue
            ret = v.get(name)
            if ret is not None:
                count += 1

        if count == 0:
            raise click.BadParameter("""{} is specified for {} but
            is not annotated in the config file""".format(name, name_str))
        if count < config_len:
            raise click.BadParameter("""{} is specified for {} but
            is not annotated for all experiments in the config file ({} / {})""".format(name, name_str,
                                                                                        count, config_len
            ))

    # if nonzero_subgroup is not None:
    #     check_group(nonzero_subgroup, 'nonzero_subgroup')
    # if batch is not None:
    #     check_group(batch, 'batch')
    # if group is not None:
    #     check_group(group, 'group')
    # if covariate is not None:
    #     check_group(covariate, 'covariate')

    # check to ensure all metadata is present in all samples
    for kw_name, value in kwargs.items():
        if value is not None:
            check_group(value, kw_name)

    return  # all passed

# @gui_option
@click.group(chain=True)
@click.option('--additional-info', type=click.Path(exists=True, dir_okay=False), default=None,
              help='.ini file with metadata for isobaric data used for scatter and PCA plots')
@click.option('--batch', type=str, default=None, help='Metadata entry to group experiments for batch correction via ComBat (Requires rpy2, R, and sva installations)')
@click.option('--batch-nonparametric', is_flag=True, default=False, help='Use nonparametric method for batch correction with ComBat (only used if --batch is also specified)')
@click.option('--batch-noimputation', is_flag=True, default=False, help='Leave original missing values after batch correction')
@click.option('--covariate', type=str, default=None, help='Metadata entry to use as covariate for batch correction via ComBat')
@click.option('--data-dir', type=click.Path(exists=False, file_okay=False),
              default='./data/', show_default=True,
              help='location to store and read e2g files')
@click.option('--file-format', type=click.Choice(('.png', '.pdf', '.svg')), default=('.png',),
              show_default=True, multiple=True, help='File format for any plots')
@click.option('--funcats', type=str, default=None, show_default=True,
              help="""Optional gene subset based on funcat or funcats,
              regular expression allowed. """)
@click.option('--funcats-inverse', type=str, default=None, show_default=True,
              help="""Optional gene filtering based on funcat or funcats, regular expression allowed. """)
@click.option('--geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True,
              help="""Optional list of geneids to subset by.
              Should have 1 geneid per line. """)
@click.option('--group', type=str, default=None, help='Metadata entry to calculate p-values for differential across (Requires rpy2, R, and sva installations)')
@click.option('--limma', default=False, is_flag=True, help='Use limma moderated t test.')
@click.option('--block', default=None, is_flag=False, help='Bio or Tech rep for block design for limma')
@click.option('--pairs', type=str, default=None, help='Metadata entry that indicates sample pairs for running pairwise statistical tests.')
@click.option('--ignore-geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True,
              help="""Optional list of geneids to ignore. Should have 1 geneid per line. """)
@click.option('--iFOT', default=False, show_default=True, is_flag=True,
              help="""Calculate iFOT (divide by total input per experiment)""")
@click.option('--iFOT-KI', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT based on kinases')
@click.option('--iFOT-TF', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT based on transcription factors')
@click.option('--median', default=False, show_default=True, is_flag=True,
              help='Normalize based on median sample expression value')
@click.option('-n', '--name', type=str, default='',
              help='An optional name for the analysis that will place all results in a subfolder.')
@click.option('--result-dir', type=click.Path(exists=False, file_okay=False),
              default='./results', show_default=True,
              help='Base directory to store results. Will be created if does not exist.')
@click.option('--taxon', type=click.Choice(['human', 'mouse', 'all']),
              default='all', show_default=True)
@click.option('--non-zeros', default=1, show_default=True, type=int_or_ratio(),
              help="""Minimum number of non zeros OR fraction of nonzeros allowed for each gene product
              across samples. If a decimal is specified (e.g. 1.0), this indicates 100% of values are nonzero.
              If an integer is specified (1), this indicates that 1 value is nonzero.
              """)
@click.option('--nonzero-subgroup', type=str, default=None, help='')
# @click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
@click.argument('experiment_file', type=Path_or_Subcommand(exists=True, dir_okay=False))
@click.pass_context
def main(ctx, additional_info, batch, batch_nonparametric, batch_noimputation, covariate, data_dir,
         file_format, funcats, funcats_inverse, geneids, group, limma, block, pairs, ignore_geneids, ifot,
         ifot_ki, ifot_tf, median, name, result_dir, taxon, non_zeros, nonzero_subgroup,
         experiment_file):
    """
    """
         # name, taxon, non_zeros, experiment_file):

    if ifot + ifot_ki + ifot_tf + median> 1:
        raise click.BadParameter('Cannot specify a combination of `iFOT`, `iFOT-KI`, `iFOT-TF`, `median`')

    # validate_subgroup(nonzero_subgroup, experiment_file)
    validate_configfile(experiment_file, nonzero_subgroup=nonzero_subgroup, batch=batch, group=group, covariate=covariate, pairs=pairs)

    metrics, metrics_after_filter, metrics_unnormed_area = False, False, True
    if 'metrics' in sys.argv:  #know this ahead of time, accumulate metrics during data load
        metrics = True
        if '--after-filter' in sys.argv:
            metrics_after_filter = True #
        if '--after-norm' in sys.argv:
            metrics_unnormed_area = False

    export_all = False
    if all(x in sys.argv for x in ('export', '--level')) and any(x in sys.argv for x in ('all', 'align')):
        #know this ahead of time, calculate more things during data load
        export_all = True

    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    if ctx.obj is None:
        ctx.obj = dict()

    ctx.obj['file_fmts'] = file_format


    analysis_name = get_file_name(experiment_file)
    if analysis_name is None:
        print('Error parsing configfile name.')
        analysis_name = 'Unnamed'

    now = datetime.now()
    context = click.get_current_context()
    params = context.params

    data_obj = Data(additional_info=additional_info, batch=batch,
                    batch_nonparametric=batch_nonparametric, batch_noimputation=batch_noimputation,
                    covariate=covariate, data_dir=data_dir, base_dir=result_dir, funcats=funcats,
                    funcats_inverse=funcats_inverse, geneids=geneids, group=group, pairs=pairs,
                    ifot=ifot, ifot_ki=ifot_ki, ifot_tf=ifot_tf, median=median, name=name,
                    non_zeros=non_zeros, nonzero_subgroup=nonzero_subgroup, taxon=taxon,
                    experiment_file=experiment_file, metrics=metrics, limma=limma,
                    block=block,
                    metrics_after_filter=metrics_after_filter,
                    metrics_unnormed_area=metrics_unnormed_area, ignore_geneids=ignore_geneids,
                    export_all=export_all
    )

    outname = get_outname('metadata', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                          batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath)+'.tab'
    data_obj.col_metadata.to_csv(outname, sep='\t')

    # cf = 'correlatioplot_args_{}.json'.format(now.strftime('%Y_%m_%d_%H_%M_%S'))
    # with open(os.path.join(data_obj.outpath, cf), 'w') as f:
    #     json.dump(params, f)
    outname = get_outname('context', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                          batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath)+'.tab'
    params = dict(ctx.params)
    params['file_format'] = ' | '.join(params['file_format'])
    param_df = pd.DataFrame(params, index=[analysis_name]).T
    param_df.to_csv(outname, sep='\t')

    ctx.obj['data_obj'] = data_obj

@main.command('scatter')
@click.option('--colors-only', default=False, is_flag=True, show_default=True,
              help="Only plot colors on correlationplot, no datapoints")
@click.option('--shade-correlation/--no-shade-correlation', default=True, is_flag=True,
              show_default=True, help="")
@click.option('--stat', type=click.Choice(['pearson', 'spearman']),
              default='pearson', show_default=True)
@click.pass_context
def scatter(ctx, colors_only, shade_correlation, stat):

    data_obj = ctx.obj['data_obj']
    file_fmts = ctx.obj['file_fmts']

    _ = data_obj.areas_log_shifted
    outname = get_outname('scatter', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=colors_only,
                          batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath)

    X = data_obj.areas_log_shifted.copy()
    # to_mask = (data_obj.mask | (data_obj.areas_log==-10))
    to_mask = (data_obj.mask | (data_obj.areas_log==data_obj.minval_log)).dropna(1)
    X[to_mask] = np.NaN
    # g = scatterplot(X.replace(0, np.NaN), stat=stat,
    g = scatterplot(X, stat=stat,
                    colors_only=colors_only, shade_correlation=shade_correlation,
                    outname=outname,
                    file_fmts=file_fmts,
                    mask=data_obj.mask
    )
    if g is None: # plotted and saved via R
        return

    save_multiple(g, outname, *file_fmts, dpi=96)

@main.command('export')
@click.option('--level', type=click.Choice(['all', 'align', 'area']), default='area',
              help="""Export data table of the filtered list of gene products used for plotting
              `all` returns all the data in long format
              `align` returns all data formatted for import into align!
              """)
@click.option('--genesymbols', default=False, is_flag=True, show_default=True,
              help='Include GeneSymbols in data export when `level` is set to `area`')
@click.pass_context
def export(ctx, level, genesymbols):

    data_obj = ctx.obj['data_obj']
    data_obj.perform_data_export(level, genesymbols=genesymbols)

@main.command('cluster')
@click.option('--col-cluster/--no-col-cluster', default=True, is_flag=True, show_default=True,
              help="""Cluster columns via hierarchical clustering.
              Note this is overridden by specifying `nclusters`""")
@click.option('--cmap', default=None, show_default=True)
@click.option('--figsize', nargs=2, type=float, default=None, show_default=True,
              help='''Optionally specify the figuresize (width, height) in inches
              If not specified, tries to use a reasonable default depending on the number of
              samples.
              ''')
@click.option('--genefile', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=False,
              help="""File of geneids to plot.
              Should have 1 geneid per line. """)
@click.option('--gene-symbols', default=False, is_flag=True, show_default=True,
              help="Show Gene Symbols on clustermap")
@click.option('--gene-symbol-fontsize', default=8, show_default=True,
              help="Gene Symbol font size")
@click.option('--highlight-geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=True,
              help="""Optional list of geneids to highlight by.
              Should have 1 geneid per line. """)
@click.option('--legend-include', type=str, multiple=True,
              help="""Specific entries in the config file to include in the legend.
              (Default all are included)""")
@click.option('--legend-exclude', type=str, multiple=True,
              help="""Specific entries in the config file to ignore for the legend.
              (Default all are included)""")
@click.option('--linkage', type=click.Choice(['single', 'complete', 'average', 'weighted', 'centroid',
                                              'median', 'ward']),
              default='ward', show_default=True,
              help='linkage method for hierarchical clustering'
)
@click.option('--max-autoclusters', default=30, show_default=True, help="""Max number of clusters to try
when `auto` is set for `--nclusters`""")
@click.option('--nclusters', default=None, callback=validate_cluster_number, show_default=True,
              help="""If specified by an integer, use that number of clusters via k-means clustering. If specified as `auto`, will try to find the optimal number of clusters""")
@click.option('--dbscan', default=False, is_flag=True, show_default=True,
              help="""Use DBSCAN algorithm to find and cluster data. Cannot be used with nclusters specification""")
@click.option('--row-cluster/--no-row-cluster', default=True, is_flag=True, show_default=True,
              help="Cluster rows via hierarchical clustering")
@click.option('--seed', default=None, help='seed for kmeans clustering', callback=validate_seed,
              show_default=True)
@click.option('--show-metadata/--hide-metadata', default=True, show_default=True,
              is_flag=True,
              help="""Show metadata on clustermap if present""")
@click.option('--standard-scale', type=click.Choice(['None', '0', '1']),
              default='None', show_default=True)
@click.option('--show-missing-values/--hide-missing-values', default=True, is_flag=True, show_default=True,
              help="""Whether or not to show missing values on the cluster plot and missing values""")
@click.option('--z-score', type=click.Choice(['None', '0', '1']),
              default='0', show_default=True)
@click.pass_context
def cluster(ctx, cmap, col_cluster, dbscan, figsize, genefile, gene_symbols, gene_symbol_fontsize,
            highlight_geneids, legend_include, legend_exclude, linkage, max_autoclusters, nclusters,
            row_cluster, seed, show_metadata, standard_scale, show_missing_values, z_score):

    if not figsize:  # returns empty tuple if not specified
        figsize = None

    if nclusters is not None and dbscan:
        raise click.BadOptionUsage('Cannot specify `nclusters` and use DBSCAN')

    genes = None
    if genefile:
        genes = parse_gid_file(genefile)

    data_obj = ctx.obj['data_obj']
    data_obj.set_highlight_gids(highlight_geneids)
    data_obj.standard_scale    = data_obj.clean_input(standard_scale)
    data_obj.z_score           = data_obj.clean_input(z_score)

    col_meta = data_obj.col_metadata
    _expids = ('recno', 'runno', 'searchno')

    valid_entries = validate_in_config(*legend_include, *legend_exclude, valid_entries=col_meta.index)
    legend_include = set(legend_include) & set(valid_entries)
    legend_exclude = set(legend_exclude) & set(valid_entries)
    col_meta = col_meta.loc[[x for x in col_meta.index if x not in _expids]]
    result = clusterplot(data_obj.areas_log_shifted,
                         cmap_name=cmap,
                         highlight_gids=data_obj.highlight_gids,
                         highlight_gid_names=data_obj.highlight_gid_names,
                         genes=genes,
                         gid_symbol=data_obj.gid_symbol,
                         gene_symbols=gene_symbols, z_score=data_obj.z_score,
                         standard_scale=data_obj.standard_scale,
                         row_cluster=row_cluster, col_cluster=col_cluster,
                         # metadata=data_obj.config if show_metadata else None,
                         col_data = col_meta if show_metadata else None,
                         nclusters=nclusters,
                         dbscan=dbscan,
                         max_autoclusters=max_autoclusters,
                         show_missing_values=show_missing_values,
                         mask=data_obj.mask,
                         figsize=figsize,
                         normed=data_obj.normed,
                         linkage=linkage,
                         gene_symbol_fontsize=gene_symbol_fontsize,
                         legend_include=legend_include,
                         legend_exclude=legend_exclude
    )

    g = result['clustermap']['clustergrid']
    extra_artists = result['clustermap']['extra_artists']
    missing_values = 'masked' if show_missing_values else 'unmasked'
    outname_func = partial(get_outname, name=data_obj.outpath_name,
                           taxon=data_obj.taxon, non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                           batch=data_obj.batch_applied,
                           batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                           outpath=data_obj.outpath, missing_values=missing_values)

    kmeans_res = result.get('kmeans')
    dbscan_res = result.get('dbscan')
    if kmeans_res is not None:
        kmeans_clusters = kmeans_res['nclusters']
        outname = outname_func('clustermap_{}clusters'.format(kmeans_clusters))
    elif dbscan_res is not None:
        dbscan_clusters = dbscan_res['nclusters']
        outname = outname_func('clustermap_{}clusters'.format(dbscan_clusters))
    else:
        outname = outname_func('clustermap')

    bbox_inches='tight'
    if extra_artists is not None:
        bbox_inches=None
        # extra_artists=None
    # save_multiple(g, outname, '.png', bbox_extra_artists=extra_artists, bbox_inches=bbox_inches)
    file_fmts = ctx.obj['file_fmts']
    save_multiple(g, outname, *file_fmts, bbox_extra_artists=extra_artists, bbox_inches=bbox_inches)


    if kmeans_res is not None:
        kmeans_data = kmeans_res['data']
        outname = os.path.abspath(outname_func('{}clusters_labels'.format(kmeans_clusters))+'.tab')
        kmeans_data.to_csv(outname, index=True, sep='\t', na_rep='NaN')
        print('Saved:', outname)

        fig = kmeans_res['auto'].get('fig')
        if fig is not None:
            outname = outname_func('cluster_optimized_results')
            save_multiple(fig, outname, *file_fmts)

        fig = kmeans_res['silhouette'].get('fig')
        if fig is not None:
            outname = outname_func('{}clusters_silhouette'.format(kmeans_clusters))
            save_multiple(fig, outname, *file_fmts)

    if dbscan_res is not None:
        dbscan_data = dbscan_res['data']

        outname = os.path.abspath(outname_func('{}clusters_labels'.format(dbscan_clusters))+'.tab')
        dbscan_data.to_csv(outname, index=True, sep='\t')
        print('Saved:', outname)

        fig = dbscan_res['silhouette'].get('fig')
        if fig is not None:
            outname = outname_func('{}clusters_silhouette'.format(dbscan_clusters))
            save_multiple(fig, outname, *file_fmts)



@main.command('pca')
@click.option('--annotate', default=False, show_default=True, is_flag=True,
              help="Annotate points on PC plot")
@click.option('--max-pc', default=2, show_default=True,
              help='Maximum PC to plot. Plots all combinations up to this number.')
@click.option('--color', default='', show_default=True, is_flag=False,
              help="What meta entry to color PCA")
@click.option('--marker', default='', show_default=True, is_flag=False,
              help="What meta entry to mark PCA")
@click.option('--genefile', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=False,
              help="""File of geneids to plot.
              Should have 1 geneid per line. """)
@click.pass_context
def pca(ctx, annotate, max_pc, color, marker, genefile):

    data_obj = ctx.obj['data_obj']

    genes = None
    if genefile:
        genes = parse_gid_file(genefile)

    # # fig, ax = pcaplot(data_obj.areas_log_shifted, data_obj.config, col_data = data_obj.col_metadata)

    figs = pcaplot(data_obj.areas_log_shifted, metadata=data_obj.config, col_data=data_obj.col_metadata,
                   annotate=annotate, max_pc=max_pc, color_label=color, marker_label=marker, genes=genes)

    # outname_func = get_outname('pcaplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
    #                            non_zeros=data_obj.non_zeros, batch=data_obj.batch_applied,
    #                            batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
    #                            colors_only=data_obj.colors_only, outpath=data_obj.outpath,
    # )

    outname_func = partial(get_outname, name=data_obj.outpath_name, taxon=data_obj.taxon,
                           non_zeros=data_obj.non_zeros, batch=data_obj.batch_applied,
                           batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                           colors_only=data_obj.colors_only, outpath=data_obj.outpath,
    )

    file_fmts = ctx.obj['file_fmts']
    for name, fig in figs.items():
        outname = outname_func(name)
        save_multiple(fig, outname, *file_fmts)

@main.command('metrics')
@click.option('--full', is_flag=True, default=False,
              help='Calculate more metrics, requires PSMs data'
)
@click.option('--before-filter/--after-filter', default=True, is_flag=True, show_default=True,
              help="Whether or not to show metrics before or after filtering")
@click.option('--before-norm/--after-norm', default=True, is_flag=True, show_default=True,
              help="Whether or not to show area before or after normalization")
@click.pass_context
def metrics(ctx, full, before_filter, before_norm):

    data_obj = ctx.obj['data_obj']
    file_fmts = ctx.obj['file_fmts']

    make_metrics(data_obj, file_fmts, before_filter=before_filter, before_norm=before_norm, full=full)

from .overlap import make_overlap
@main.command('overlap')
@click.option('--group', type=str, default=None, help='Metadata entry to group samples for assessing overlap')
@click.option('--non-zeros', default=1., show_default=True, type=int_or_ratio(),
              help="""Minimum number of non zeros OR fraction of nonzeros allowed for each sample
              (or sample group. If a decimal is specified (e.g. 1.0), this indicates 100% of values are nonzero.
              If an integer is specified (1), this indicates that 1 value is nonzero.
              """)
@click.pass_context
def overlap(ctx, group, non_zeros):
    """
    Plot gene product overlap across experiments
    """
    data_obj = ctx.obj['data_obj']
    file_fmts = ctx.obj['file_fmts']

    if group:
        validate_configfile(data_obj.experiment_file, group=group)

    make_overlap(data_obj, group=group, non_zeros=non_zeros)


@main.command('volcano')
@click.option('-f', '--foldchange', type=int, default=4, show_default=True,
              help='fold change cutoff'
)
@click.option('-e', '--expression-data', default=False, is_flag=True, show_default=True,
              help='Include expression data for each sample in tabular output.'
)
@click.option('-n', '--number', type=int, default=35, show_default=True,
              help='Maximum number of significant genes to highlight (annotate) in plot'
)
@click.option('-o', '--only-sig', default=False, is_flag=True, show_default=True,
              help='Only export genes that are significantly different (based on set cutoff)'
)
@click.option('-s', '--scale', type=float, default=1.2, show_default=True,
              help='To what extent to scale the labels'
)
@click.option('--sig', type=float, default=.05, show_default=True,
              help='Significance cutoff for (B.H. adjusted) pvalue'
)
@click.option('--p-adj/--p-value', default=True, is_flag=True, show_default=True,
              help="Whether to plot padj or pvalue on volcano plot (does not change underlying data)")
@click.option('--highlight-geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=False,
              help="""Optional list of geneids to also highlight. Should have 1 geneid per line. """)
@click.pass_context
def volcano(ctx, foldchange, expression_data, number, only_sig, sig, scale, p_adj, highlight_geneids):
    """
    Draw volcanoplot and highlight significant (FDR corrected pvalue < .05 and > 2 fold change)
    """
    from .volcanoplot import volcanoplot
    yaxis = 'pAdj' if p_adj else 'pValue'
    volcanoplot(ctx, foldchange, expression_data, number=number, only_sig=only_sig, sig=sig,
                yaxis=yaxis, scale=scale, highlight_geneids=highlight_geneids)



@main.command('gsea')
@click.option('--show-result/--no-show-result', default=True, show_default=True)
@click.option('--collapse', type=bool, default=False, show_default=True)
@click.option('--geneset', type=click.Choice(('hallmark', 'go_biological', 'curated.CP.all',
                                              'curated.CP.KEGG', 'oncogenic', 'curated.CP.BioCarta',
                                              'curated.CP.Reactome', 'curated.CGP',
                                              'go.All', 'go.Bio', 'go.Cell', 'go.Molecular'
)),
              default=('hallmark',), show_default=True, multiple=True
)
@click.option('--metric', type=click.Choice(('Signal2Noise', 'tTest', 'Cosine', 'Euclidean', 'Manhatten',
                                             'Pearson', 'Ratio_of_Classes', 'Diff_of_Classes')),
              default='Signal2Noise', show_default=True)
@click.option('--mode', type=click.Choice(('Max_probe', 'Median_of_probes')),
              default='Max_probe', show_default=True)
@click.option('--norm', type=click.Choice(('meandiv', 'None')), default='meandiv', show_default=True)
@click.option('--number-of-permutations', default=1000, show_default=True)
@click.option('--permute', type=click.Choice(('phenotype', 'gene_set')), default='phenotype',
              show_default=True,
)
@click.option('--plot-top-x', default=20, show_default=True)
@click.option('--rnd-type', type=click.Choice(('no_balance', 'equalize_and_balance')),
              default='no_balance', show_default=True)
@click.option('--scoring-scheme', type=click.Choice(('weighted', 'classic', 'weighted_p2', 'weighted_p1.5')),
              default='weighted', show_default=True)
@click.option('--set-max', default=500, show_default=True)
@click.option('--set-min', default=15, show_default=True)
@click.option('--sort', type=click.Choice(('real', 'abs')), default='real', show_default=True)
@click.option('-n', '--number', default=9999, help='Number of pathways to plot in output', show_default=True)
@click.pass_context
def gsea(ctx, show_result, collapse, geneset, metric, mode, number_of_permutations, norm, permute,
         plot_top_x, rnd_type, scoring_scheme, sort, set_max, set_min, number):
    """
    Run GSEA on specified groups
    """
    data_obj = ctx.obj['data_obj']
    file_fmts = ctx.obj['file_fmts']

    plt.rc('font',**{'family':'sans-serif','sans-serif':["DejaVu Sans", "Arial", "Liberation Sans",
                            "Bitstream Vera Sans", "sans-serif"]})
    group = data_obj.group #
    if group is None:
        raise ValueError('Must specify Group')


    # expression = data_obj.areas_log_shifted.copy().fillna(0)
    expression = data_obj.areas_log_shifted.copy().fillna('na')
    expression.index.name = 'NAME'

    pheno = data_obj.col_metadata


    nsamples = len(pheno.columns)
    ngroups  = pheno.loc[group].nunique()
    groups   = pheno.loc[group].unique()
    # pheno_indicator = dict()
    # for ix, grp in enumerate(groups):
    #     pheno_indicator[grp] = ix
    # classes  = list(map(str, [pheno_indicator[grp] for grp in pheno.loc[group]]))
    classes = [grp.replace(' ', '_') for grp in pheno.loc[group]]

    cls_comparison = ''
    if ngroups == 2:  # reverse it
        cls_comparison = '#{1}_versus_{0}'.format(*groups)
    elif ngroups != 2:
        raise ValueError('Must have 2 groups')


    namegen = partial(get_outname, name=data_obj.outpath_name, taxon=data_obj.taxon,
                      non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                      batch=data_obj.batch_applied,
                      batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                      outpath=data_obj.outpath)


    # param_file = os.path.abspath(namegen('gsea_params') + '.txt')
    param_file = namegen('gsea_params') + '.txt'

    gsea_jar = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                            'GSEA', 'gsea-3.0.jar')

    geneset_mapping = {'hallmark': 'h.all.v6.1.entrez.gmt',
                       'go_biological': 'c5.bp.v6.1.entrez.gmt',
                       'curated.CGP': 'c2.all.v6.1.entrez.gmt',
                       'curated.CP.all': 'c2.cp.v6.1.entrez.gmt',
                       'curated.CP.BioCarta': 'c2.cp.biocarta.v6.1.entrez.gmt',
                       'curated.CP.KEGG': 'c2.cp.kegg.v6.1.entrez.gmt',
                       'curated.CP.Reactome': 'c2.cp.reactome.v6.1.entrez.gmt',
                       'oncogenic': 'c6.all.v6.1.entrez.gmt',
                       'go.All': 'c5.all.v6.1.entrez.gmt',
                       'go.Bio': 'c5.bp.v6.1.entrez.gmt',
                       'go.Cell': 'c5.cc.v6.1.entrez.gmt',
                       'go.Molecular': 'c5.mf.v6.1.entrez.gmt'
                       # : 'msigdb.v6.1.entrez.gmt',
    }
    # get most recent, sort by name and take last
    homologene_f = sorted(glob.glob(os.path.join(os.path.split(os.path.abspath(__file__))[0],
                                                 'data', 'homologene*data')),
                          reverse=True)[0]


    homologene = (pd.read_table(homologene_f, header=None,
                                names=('Homologene', 'TaxonID', 'GeneID',
                                       'Symbol', 'ProteinGI', 'ProteinAccession'))
    )
    # check if we have non-human GeneIDs
    hgene_query = homologene[ homologene.GeneID.isin(expression.index) ]
    if hgene_query.TaxonID.nunique() > 1:
        raise ValueError('No support for multi-species GSEA')
    if hgene_query.TaxonID.nunique() == 1 and hgene_query.TaxonID.unique()[0] != 9606:
        # remap
        print('Remapping {} GeneIDs to human'.format(hgene_query.TaxonID.unique()[0]))
        gid_hgene = hgene_query[['GeneID', 'Homologene']].set_index('GeneID')['Homologene'].to_dict()
        hgene_hugid = (homologene.query('TaxonID==9606') [['GeneID', 'Homologene']]
                       .set_index('Homologene')['GeneID'].to_dict()
        )
        expression.index = expression.index.map( lambda x: hgene_hugid.get( gid_hgene.get(x) ))

    expression = expression.loc[ expression.index.dropna() ]
    expression.index = expression.index.astype(int)
    if expression.index.nunique() < len(expression.index):
        expression = expression.groupby(expression.index).mean()


    collapse = 'true' if collapse else 'false'

    for gs in geneset:

        outname = namegen('gsea', pathway=gs)
        outdir = os.path.abspath(os.path.join(data_obj.outpath, 'gsea'))
        report_name = os.path.split(outname)[-1]

        f = geneset_mapping.get(gs)

        geneset_file = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                                    'GSEA', 'genesets', f)

        # with NamedTemporaryFile(suffix='.txt') as f, NamedTemporaryFile(mode='w', suffix='.cls') as g:
        with named_temp(suffix='.txt') as f, named_temp(mode='w', suffix='.cls') as g:
            expression.to_csv(f.name, sep='\t')
            f.close() # windows compat
        # with open('./results/gsea/kip_dda_kinases_pheno.cls', 'w') as g:
            g.write('{} {} 1\n'.format(nsamples, ngroups))
            g.write('# {}\n'.format(' '.join(groups)))
            # g.write('{}\n'.format(' '.join(classes)))
            g.write('{}\n'.format(' '.join(pheno.loc[group])))
            g.file.flush()
            g.close() # windows compat

            # rpt_name\t{rpt_name}
            params = """
            collapse\t{collapse}
            metric\t{metric}
            mode\t{mode}
            norm\t{norm}
            order\tdescending
            include_only_symbols\tfalse
            permute\t{permute}
            plot_top_x\t{plot_top_x}
            rnd_type\t{rnd_type}
            set_max\t{set_max}
            set_min\t{set_min}
            nperm\t{nperm}
            res\t{res}
            cls\t{cls}{cls_comparison}
            gmx\t{gmx}
            out\t{outdir}
            rpt_label\t{rpt_label}
            gui\tfalse
            """.format(collapse=collapse, metric=metric, mode=mode, norm=norm, permute=permute,
                    rnd_type=rnd_type, scoring_scheme=scoring_scheme, sort=sort, set_max=set_max,
                    set_min=set_min, nperm=number_of_permutations, plot_top_x=plot_top_x,
                    cls_comparison=cls_comparison,
                    # rpt_name=report_name,
                    res=os.path.abspath(f.name), cls=os.path.abspath(g.name),
                    gmx=os.path.abspath(geneset_file), outdir=outdir, rpt_label=report_name )
            with open(param_file, 'w') as f:
                f.write(params)


            res = subprocess.run(['java', '-Xmx8192m', '-cp', gsea_jar, 'xtools.gsea.Gsea',
                                '-param_file', param_file,
            ],
                                # stdout=subprocess.PIPE
            )

            res.check_returncode()

        folders = [os.path.abspath(os.path.join(outdir, f)) for f in os.listdir(outdir) if ('Gsea' in f)]
        folders.sort(key=os.path.getmtime)

        new_folder = folders[-1]
        index = os.path.join(new_folder, 'index.html')
        if sys.platform == "darwin":  # check if on OSX
            index = 'file://' + index

        if show_result:
            import webbrowser
            webbrowser.open(index)

        # parse result
        # GSEA outputs the summary files of the form:
        # gsea_report_for_[groupname]_[digit_timestamp].xls
        group0 = glob.glob( os.path.join(new_folder, 'gsea_report_for_{}_[0-9]*.xls'.format(groups[0])) )
        group1 = glob.glob( os.path.join(new_folder, 'gsea_report_for_{}_[0-9]*.xls'.format(groups[1])) )
        assert len(group0) == len(group1) == 1
        group0_df = pd.read_table(group0[0], index_col='NAME')
        group1_df = pd.read_table(group1[0], index_col='NAME')

        # gsea_sig = group0_df[ group0_df['FWER p-val' ] < .55 ].join(
        #     group1_df[ group1_df['FWER p-val'] < .25 ],
        #     lsuffix='_group0', rsuffix='_group1',
        #     how='outer'
        # )

        cmap = mpl.cm.Reds_r
        bounds = np.linspace(0, 1, 21)
        cnorm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        powernorm = mpl.colors.PowerNorm(.5, vmin=0, vmax=1)

        gsea_sig = pd.DataFrame()
        cutoff = .25
        while len(gsea_sig) < 5:
            if cutoff > 1:
                break
            gsea_sig = pd.concat([group0_df[ group0_df['FWER p-val'] < cutoff ],
                                  group1_df[ group1_df['FWER p-val'] < cutoff ],
            ])
            cutoff += .1

        if gsea_sig.empty:
            print('No gene sets to plot!')
            # return
        tokeep = gsea_sig.NES.sort_values(ascending=False).head(number).index
        idx = [x for x in gsea_sig.index if x in tokeep]
        gsea_sig = gsea_sig.loc[idx]

        gsea_sig['color'] = gsea_sig['FWER p-val'].apply(lambda x:
                                mpl.colors.to_hex( cmap(cnorm(powernorm(x))) )
        )
        import textwrap
        gsea_sig.index = gsea_sig.index +  ['*' if x<.25 else '' for x in gsea_sig['FWER p-val'] ]
        gsea_sig.index = gsea_sig.index.map(lambda x: textwrap.fill(x.replace('_', ' '),
                                                                    24,
                                                                    break_long_words=False) )
        # print(gsea_sig[['FWER p-val', 'NES', 'color']])

        # mpl.colors.Normalize(vmin=1.,vmax=1.)
        if gsea_sig.NES.max() * gsea_sig.NES.min() < 0:
            nes_range = gsea_sig.NES.max() + abs(gsea_sig.NES.min())
        else:
            nes_range = max(abs(gsea_sig.NES.max()), 2.5)

        figwidth = np.round(nes_range*2.5, decimals=1)


        figheight = max(4, min(gsea_sig.pipe(len) // 1.5, 14))
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = mpl.gridspec.GridSpec(2, 1,
                                   width_ratios=[1,],
                                   height_ratios=[19, 1],
                                   hspace=.4,
        )
        # ax0 = plt.subplot(gs[0])
        ax0 = fig.add_subplot(gs[0])
        # ax1 = plt.subplot(gs[1])
        # cax = plt.subplot(gs[1:])
        cax = fig.add_subplot(gs[1:])
        gsea_sig['NES'].fillna(0).plot.barh(ax=ax0, color=gsea_sig.color, edgecolor='#222222',
                                            linewidth=2)
        ax0.axvline(color='#222222')
        # gsea_sig['NES_group1'].fillna(0).plot.barh(colormap=cmap, ax=ax1)
        # ax1.set_yticklabels(())

        gradient = np.apply_along_axis(lambda x: cnorm(powernorm(x)), 0, np.linspace(0, 1, 256))
        gradient = np.vstack((gradient, gradient))
        cax.imshow(gradient, aspect='auto', cmap=cmap)
        cax.set_yticklabels(())

        start, end = cax.get_xlim()
        cax.xaxis.set_ticks(np.linspace(start, end, 11))
        cax.set_xticklabels( ['{:.2f}'.format(x) for x in np.linspace(0, 1, 11)] )
        cax.set_xlabel('FWER p-val')
        ax0.grid(axis='x')

        ax0.set_ylabel('')
        ax0.set_xlabel('NES')

        for tick in ax0.yaxis.get_ticklabels():
            txt = tick.get_text()
            if len(txt) > 20:
                size = 8
            elif len(txt) > 30:
                size = 6
            else:
                size = 9
            tick.set_size(size)


        # plt.tight_layout()
        # fig = plt.gcf()
        groups[0], groups[1]
        ax0.text(0, 1.04, groups[0], transform=ax0.transAxes)
        ax0.text(1, 1.04, groups[1], transform=ax0.transAxes, ha='right')
        gs.tight_layout(fig, rect=(0, 0, 1, .96))
        # fig.subplots_adjust(left=.4)
        # fig.tight_layout()
        save_multiple(fig, outname, *file_fmts)

@main.command('bar')
@click.option('--gene', type=int,
              default=None, show_default=True, multiple=True,
              help="Gene to plot. Multiple allowed")
@click.option('--genefile', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=False,
              help="""File of geneids to plot.
              Should have 1 geneid per line. """)
@click.option('--average', type=str, default=None, help='Metadata entry to group data and plot average. Cannot be used with color')
@click.option('--color', type=str, default=None, help='Metadata entry to group color bars. Cannot be used with average')
@click.option('--cmap', default=None, show_default=True, help="""
Any valid, qualitative, matplotlib colormap. See https://matplotlib.org/examples/color/colormaps_reference.html.
""")
@click.option('--linear', default=False, is_flag=True, help='Plot linear values (default log10)')
@click.pass_context
def bar(ctx, average, color, cmap, gene, genefile, linear):

    data_obj = ctx.obj['data_obj']
    col_meta = data_obj.col_metadata
    if genefile:
        gene = gene + tuple(parse_gid_file(genefile))
    if len(gene) == 0:
        raise ValueError("Must supply at least 1 gene")

    outpath = os.path.join(data_obj.outpath, 'bar')
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    outfunc = partial(get_outname, name=data_obj.outpath_name, taxon=data_obj.taxon,
                      non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                      batch=data_obj.batch_applied,
                      batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                      outpath=outpath,
    )
    if color is not None and average is not None:
        raise ValueError('Cannot specify color and average at the same time.')

    barplot(data_obj.areas_log_shifted, genes=gene, color=color, cmap=cmap, metadata=col_meta, average=average,
            linear=linear, base_outfunc=outfunc, file_fmts=ctx.obj['file_fmts'], gid_symbol=data_obj.gid_symbol)

if __name__ == '__main__':

    main(obj={})
