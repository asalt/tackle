"""

"""
import sys
import os
import re
import json
from datetime import datetime
import operator as op
from collections import OrderedDict
from functools import partial
import copy

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
from .pcaplot import pcaplot
from .utils import *
from .containers import Data
# from cluster_to_plotly import cluster_to_plotly

sys.setrecursionlimit(10000)

TAXON_MAPPER = {'human': 9606,
                'mouse': 10090}


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
        print(outname)

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

def validate_seed(ctx, param, value):

    if value == 'None' or value is None:
        return None

    elif value.isdigit():
        return int(value)

    else:
        raise click.BadParameter('Must be an integer or `None`')

def validate_configfile(experiment_file, nonzero_subgroup=None, batch=None, group=None):
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

    if nonzero_subgroup is not None:
        check_group(nonzero_subgroup, 'nonzero_subgroup')
    if batch is not None:
        check_group(batch, 'batch')
    if group is not None:
        check_group(group, 'group')

    return


@click.group(chain=True)
@click.option('--additional-info', type=click.Path(exists=True, dir_okay=False), default=None,
              help='.ini file with metadata for isobaric data used for scatter and PCA plots')
@click.option('--batch', type=str, default=None, help='Metadata entry to group experiments for batch correction via ComBat (Requires rpy2, R, and sva installations)')
@click.option('--batch-nonparametric', is_flag=True, default=False, help='Use nonparametric method for batch correction with ComBat (only used if --batch is also specified)')
@click.option('--data-dir', type=click.Path(exists=False, file_okay=False),
              default='./data/', show_default=True,
              help='location to store and read e2g files')
@click.option('--file-format', type=click.Choice(('.png', '.pdf', '.svg')), default=('.png',),
              show_default=True, multiple=True, help='File format for any plots')
@click.option('--funcats', type=str, default=None, show_default=True,
              help="""Optional gene subset based on funcat or funcats,
              regular expression allowed. """)
@click.option('--geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True,
              help="""Optional list of geneids to subset by.
              Should have 1 geneid per line. """)
@click.option('--group', type=str, default=None, help='Metadata entry to calculate p-values for differential across (Requires rpy2, R, and sva installations)')
@click.option('--iFOT', default=False, show_default=True, is_flag=True,
              help="""Calculate iFOT (divide by total input per experiment)""")
@click.option('--iFOT-KI', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT based on kinases')
@click.option('--iFOT-TF', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT based on transcription factors')
@click.option('-n', '--name', type=str, default='',
              help='An optional name for the analysis that will place all results in a subfolder.')
@click.option('--result-dir', type=click.Path(exists=False, file_okay=False),
              default='./results', show_default=True,
              help='Base directory to store results. Will be created if does not exist.')
@click.option('--taxon', type=click.Choice(['human', 'mouse', 'all']),
              default='all', show_default=True)
@click.option('--non-zeros', default=0, show_default=True,
              help='Minimum number of non zeros allowed for each gene product across samples.')
@click.option('--nonzero-subgroup', type=str, default=None, help='')
# @click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
@click.argument('experiment_file', type=Path_or_Subcommand(exists=True, dir_okay=False))
@click.pass_context
def main(ctx, additional_info, batch, batch_nonparametric, data_dir, file_format, funcats, geneids,
         group, ifot, ifot_ki, ifot_tf, name, result_dir, taxon, non_zeros, nonzero_subgroup, experiment_file):
    """
    """
         # name, taxon, non_zeros, experiment_file):

    if ifot + ifot_ki + ifot_tf > 1:
        raise click.BadParameter('Cannot specify a combination of `iFOT`, `iFOT-KI` and `iFOT-TF`')

    # validate_subgroup(nonzero_subgroup, experiment_file)
    validate_configfile(experiment_file, nonzero_subgroup=nonzero_subgroup, batch=batch, group=group)

    metrics, metrics_after_filter = False, False
    if 'metrics' in sys.argv:  #know this ahead of time, accumulate metrics during data load
        metrics = True
        if '--after-filter' in sys.argv:
            metrics_after_filter = True #

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

    data_obj = Data(additional_info=additional_info, batch=batch, batch_nonparametric=batch_nonparametric,
                    data_dir=data_dir, base_dir=result_dir,
                    funcats=funcats, geneids=geneids, group=group, ifot=ifot, ifot_ki=ifot_ki,
                    ifot_tf=ifot_tf, name=name, non_zeros=non_zeros,
                    nonzero_subgroup=nonzero_subgroup, taxon=taxon, experiment_file=experiment_file,
                    metrics=metrics, metrics_after_filter=metrics_after_filter
    )

    # cf = 'correlatioplot_args_{}.json'.format(now.strftime('%Y_%m_%d_%H_%M_%S'))
    # with open(os.path.join(data_obj.outpath, cf), 'w') as f:
    #     json.dump(params, f)

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
    to_mask = (data_obj.mask | (data_obj.areas_log==-10))
    X[to_mask] = np.NaN
    g = scatterplot(X.replace(0, np.NaN), stat=stat,
                    colors_only=colors_only, shade_correlation=shade_correlation,
                    outname=outname,
                    file_fmts=file_fmts,
                    mask=data_obj.mask
    )
    if g is None: # plotted and saved via R
        return

    save_multiple(g, outname, *file_fmts, dpi=96)

@main.command('export')
@click.option('--level', type=click.Choice(['all', 'area']), default='area',
              help="""Export data table of the filtered list of gene products used for plotting""")
@click.pass_context
def export(ctx, level):

    data_obj = ctx.obj['data_obj']
    data_obj.perform_data_export(level)

@main.command('cluster')
@click.option('--col-cluster/--no-col-cluster', default=True, is_flag=True, show_default=True,
              help="""Cluster columns via hierarchical clustering.
              Note this is overridden by specifying `nclusters`""")
@click.option('--figsize', nargs=2, type=float, default=None, show_default=True,
              help='''Optionally specify the figuresize (width, height) in inches
              If not specified, tries to use a reasonable default depending on the number of
              samples.
              ''')
@click.option('--gene-symbols', default=False, is_flag=True, show_default=True,
              help="Show Gene Symbols on clustermap")
@click.option('--highlight-geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=True,
              help="""Optional list of geneids to highlight by.
              Should have 1 geneid per line. """)
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
def cluster(ctx, col_cluster, dbscan, figsize, gene_symbols, highlight_geneids, linkage, max_autoclusters,
            nclusters, row_cluster, seed, show_metadata, standard_scale, show_missing_values,
            z_score):

    if not figsize:  # returns empty tuple if not specified
        figsize = None

    if nclusters is not None and dbscan:
        raise click.BadOptionUsage('Cannot specify `nclusters` and use DBSCAN')

    data_obj = ctx.obj['data_obj']
    data_obj.set_highlight_gids(highlight_geneids)
    data_obj.standard_scale    = data_obj.clean_input(standard_scale)
    data_obj.z_score           = data_obj.clean_input(z_score)
    result = clusterplot(data_obj.areas_log_shifted, highlight_gids=data_obj.highlight_gids,
                         highlight_gid_names=data_obj.highlight_gid_names,
                         gid_symbol=data_obj.gid_symbol,
                         gene_symbols=gene_symbols, z_score=data_obj.z_score,
                         standard_scale=data_obj.standard_scale,
                         row_cluster=row_cluster, col_cluster=col_cluster,
                         # metadata=data_obj.config if show_metadata else None,
                         col_data = data_obj.col_metadata if show_metadata else None,
                         nclusters=nclusters,
                         dbscan=dbscan,
                         max_autoclusters=max_autoclusters,
                         show_missing_values=show_missing_values,
                         mask=data_obj.mask,
                         figsize=figsize,
                         normed=data_obj.normed,
                         linkage=linkage
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
@click.pass_context
def pca(ctx):

    data_obj = ctx.obj['data_obj']

    fig, ax = pcaplot(data_obj.areas_log_shifted, data_obj.config, col_data = data_obj.col_metadata)
    outname = get_outname('pcaplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          colors_only=data_obj.colors_only, outpath=data_obj.outpath,
    )
    file_fmts = ctx.obj['file_fmts']
    save_multiple(fig, outname, *file_fmts)

@main.command('metrics')
@click.option('--full', is_flag=True, default=False,
              help='Calculate more metrics, requires PSMs data'
)
@click.option('--before-filter/--after-filter', default=True, is_flag=True, show_default=True,
              help="Whether or not to show metrics before or after filtering")
@click.pass_context
def metrics(ctx, full, before_filter):

    rc = {'font.family': 'sans-serif',
          "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans",
                              "Bitstream Vera Sans", "sans-serif"],
          'legend.frameon': True,
    }

    sb.set_context('notebook')
    sb.set_palette('muted')
    sb.set_color_codes()
    sb.set_style('white', rc)

    data_obj = ctx.obj['data_obj']
    file_fmts = ctx.obj['file_fmts']

    data = data_obj.metric_values

    if before_filter:
        kws = dict(before='filter')
    else:
        kws = dict(after='filter')

    outname = get_outname('metrics', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                          batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath,
                          **kws
    )

    sra = pd.DataFrame(data=[data[n]['SRA'] for n in data.keys()], index=data.keys())
    gpg = pd.DataFrame(data=[data[n]['GPGroups'] for n in data.keys()], index=data.keys())
    psms = pd.DataFrame(data=[data[n]['PSMs'] for n in data.keys()], index=data.keys())
    peptides = pd.DataFrame(data=[data[n]['Peptides'] for n in data.keys()], index=data.keys())
    area = OrderedDict((n, data[n]['Area']) for n in data.keys())

    frames = list()
    for name, value in area.items():
        frame = pd.Series(value).to_frame('AreaSum_dstrAdj').multiply(1e9).apply(np.log10)
        frame['Name'] = name
        frames.append(frame)
    area_df = pd.concat(frames)

    # area = pd.DataFrame(data=[data[n]['area'] for n in data.keys()], index=data.keys())

    green = 'darkgreen'; yellow = 'gold'; red ='firebrick'

    # fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,8), sharex=True, sharey=False)
    fig = plt.figure(figsize=(18,10))
    gs = gridspec.GridSpec(2,3)

    ax_sra  = fig.add_subplot(gs[0, 0])
    ax_gpg  = fig.add_subplot(gs[0, 1])
    ax_psms = fig.add_subplot(gs[1, 0])
    ax_pept = fig.add_subplot(gs[1, 1])
    ax_area = fig.add_subplot(gs[0:, 2])

    sra[['S', 'R', 'A']].plot.bar(stacked=True, ax=ax_sra, color=[green, yellow, red], title='SRA')
    gpg.plot.bar(ax=ax_gpg, legend=False, title='Gene Product Groups')
    psms[['Total', 'u2g']].plot.bar(stacked=False, ax=ax_psms, title='PSMs', width=.75)
    peptides[['Total', 'u2g', 'Strict', 'Strict_u2g']].plot.bar(stacked=False, ax=ax_pept, title='Peptides',
                                                                width=.8)
    plt.setp(ax_sra.get_xticklabels(), visible=False)
    plt.setp(ax_gpg.get_xticklabels(), visible=False)

    sb.violinplot(y='Name', x='AreaSum_dstrAdj', data=area_df, ax=ax_area)
    ax_area.yaxis.tick_right()
    ax_area.set_ylabel('')
    ax_area.set_xlabel('log$_{10}$ AreaSum dstrAdj')
    # plt.setp( ax_area.xaxis.get_majorticklabels(), rotation=90 )
    plt.setp( ax_area.xaxis.get_majorticklabels(), rotation=0 )

    FONTSIZE = 10
    ncols = {0: 1, 2:2, 3:1}
    for ix, ax in enumerate((ax_sra, ax_gpg, ax_psms, ax_pept,)):
        ax.yaxis.grid(True, lw=.25, color='grey', ls=':')
        for tick in ax.xaxis.get_ticklabels():
            txt = tick.get_text()
            newsize = FONTSIZE

            if len(txt) > 7 and len(txt) <= 9:
                newsize -= 1
            if len(txt) >= 10 and len(txt) < 13:
                newsize -= 1
            elif len(txt) >= 13:
                newsize -= 1
            tick.set_size(newsize)
        sb.despine(ax=ax)
        if ix == 1:
            continue  #  no legend for gpgroups
        ax.legend(loc='best', ncol=ncols[ix])

    sb.despine(ax=ax_area, right=False, top=True, left=True, bottom=False)
    ax_area.xaxis.grid(True, lw=.25, color='grey', ls=':')

    fig.subplots_adjust(hspace=.15, top=.95, left=.1, right=.9)
    save_multiple(fig, outname, *file_fmts)

    if full:  # also plot info from PSMs
        return # not done
        config = data_obj.config
        data_dir = data_obj.data_dir

        psms = OrderedDict()
        for name, record in config.items():
            if name.startswith('__'):
                continue
            recno = record.get('recno')
            runno = record.get('runno')
            searchno = record.get('searcno')
            df = ispec.PSMs(recno, runno, searchno, data_dir=data_dir).df.query('oriFLAG==1')

            df.index = pd.to_timedelta( df.RTmin, unit='m' )

            psms[name] = df



        nrows = len(psms)

        m = max([ pd.Timedelta(df.RTmin.max(), unit='m') for df in psms.values() ])

        fig = plt.figure(figsize=(12, 8))
        gs = gridspec.GridSpec(nrows, 3, width_ratios=(3,1,3), left=.1, right=.9, bottom=.1, top=.9)
        # plot IDs/min over time
        prev_ax = None
        for ix, (name, df) in enumerate(psms.items()):
            if prev_ax:
                # ax = fig.add_subplot(gs[ix, 0], sharex=prev_ax, sharey=prev_ax)
                ax = fig.add_subplot(gs[ix, 0],  sharey=prev_ax)
            else:
                ax = fig.add_subplot(gs[ix, 0])
            # df.index = pd.to_timedelta( df.RTmin, unit='m' )
            df.loc[m] = np.NaN
            g = df.groupby([ pd.Grouper(freq='Min'), ])
            g.size().plot(ax=ax)
            label_ax = fig.add_subplot(gs[ix, 1])
            label_ax.annotate(name, xy=(0.5, 0.5), xycoords='axes fraction',
                              va='center', ha='center', size=12
            )
            sb.despine(ax=label_ax, left=True, bottom=True)
            plt.setp(label_ax.get_xticklabels(), visible=False)
            plt.setp(label_ax.get_yticklabels(), visible=False)
            ax.xaxis.grid(True, lw=.25, color='grey', ls=':')
            if ix < nrows-1:
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.set_xlabel('')
            sb.despine(ax=ax)
            prev_ax = ax




@main.command('volcano')
@click.option('-n', '--number', type=int, default=35, show_default=True,
              help='Maximum number of significant genes to highlight (annotate) in plot'
)
@click.pass_context
def volcano(ctx, number):
    """
    Draw volcanoplot and highlight significant (FDR corrected pvalue < .05 and > 2 fold change)
    """

    data_obj = ctx.obj['data_obj']

    group = data_obj.group #
    if group is None:
        print('Must supply a group value.')
        return

    if data_obj.col_metadata.loc[group].nunique() != 2:
        print('Error in volcanoplot, number of groups must be exactly 2.')
        return
    if data_obj.col_metadata.loc[group].value_counts().min() < 3:
        print('Each group must have at least 3 replicates.')
        return

    groups = dict()
    for grp in data_obj.col_metadata.loc[group].unique():
        samples = ((data_obj.col_metadata.loc[group] == grp)
                   .where(lambda x: x)
                   .dropna()
                   .index
        )
        groups[grp] = samples

    group0, group1 = data_obj.col_metadata.loc[group].values[[0, -1]]
    samples0, samples1 = groups[group0], groups[group1]


    values = data_obj.areas_log_shifted

    qvals   = data_obj.qvalues

    log2_fc = ((values[samples1].mean(1) - values[samples0].mean(1))
               .apply(lambda x: np.power(10, x))
               # .pipe(np.power, 10)
               .pipe(np.log2)
    )
    log2_fc.name = 'log2_Fold_Change'

    df = qvals.join(log2_fc.to_frame())
    df['GeneSymbol'] = df.index.map(lambda x: data_obj.gid_symbol.get(x, '?'))
    df['FunCats']    = df.index.map(lambda x: data_obj.gid_funcat_mapping.get(x, ''))
    df.index.name = 'GeneID'

    outname = get_outname('volcanoplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                          batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath)

    out = outname + '.tab'
    print("Saving", out, '...', end='', flush=True)
    df.to_csv(out, sep='\t')
    print('done', flush=True)

    try:
        from rpy2.robjects import r
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        _viaR = True
    except ModuleNotFoundError:
        _viaR = False

    if _viaR:

        pandas2ri.activate()
        r_source = robjects.r['source']
        r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                              'R', 'volcanoplot.R')
        r_source(r_file)
        Rvolcanoplot = robjects.r['volcanoplot']

        file_fmts = ctx.obj['file_fmts']
        grdevices = importr('grDevices')
        gr_devices = {'.png': grdevices.png,
                      '.pdf': grdevices.pdf,
                      '.svg': grdevices.svg}
        gr_kws = {'.png': dict(width=5, height=5, units='in', res=300),
                  '.pdf': dict(width=5, height=5,),
                  '.svg': dict(width=5, height=5,)
        }
        for file_fmt in file_fmts:

            grdevice = gr_devices[file_fmt]
            gr_kw = gr_kws[file_fmt]
            out = outname + file_fmt
            print("Saving", out, '...', end='', flush=True)

            grdevice(file=out, **gr_kw)

            Rvolcanoplot(pandas2ri.py2ri(df.reset_index()), max_labels=number)

            grdevices.dev_off()
            print('done.', flush=True)

    else:
        print('Must install rpy2')


if __name__ == '__main__':

    main(obj={})
