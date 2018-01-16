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

__version__ = '0.38'

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

def validate_configfile(experiment_file, nonzero_subgroup=None, batch=None):
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

    def check_group(name):
        count = 0
        for k, v in config.items():
            if k in dunder_fields:
                continue
            ret = v.get(nonzero_subgroup)
            if ret is not None:
                count += 1

        if count == 0:
            raise click.BadParameter("""{} is specified for nonzero subgroups but
            is not annotated in the config file""".format(nonzero_subgroup))
        if count < config_len:
            raise click.BadParameter("""{} is specified for nonzero subgroups but
            is not annotated for all experiments in the config file ({} / {})""".format(nonzero_subgroup,
                                                                                        count, config_len
            ))
    if nonzero_subgroup is not None:
        check_group(nonzero_subgroup)
    if batch is not None:
        check_group(batch)

    return


@click.group(chain=True)
@click.option('--additional-info', type=click.Path(exists=True, dir_okay=False), default=None,
              help='.ini file with metadata for isobaric data used for scatter and PCA plots')
@click.option('--batch', type=str, default=None, help='Metadata entry to group experiments for batch correction via ComBat (Requires rpy2, R, and sva installations)')
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
@click.option('--iFOT', default=False, show_default=True, is_flag=True,
              help="""Calculate iFOT (divide by total input per experiment)""")
@click.option('--iFOT-KI', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT based on kinases')
@click.option('--iFOT-TF', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT based on transcription factors')
@click.option('-n', '--name', type=str, default='',
              help='An optional name for the analysis that will place all results in a subfolder.')
@click.option('--taxon', type=click.Choice(['human', 'mouse', 'all']),
              default='all', show_default=True)
@click.option('--non-zeros', default=0, show_default=True,
              help='Minimum number of non zeros allowed for each gene product across samples.')
@click.option('--nonzero-subgroup', type=str, default=None, help='')
# @click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
@click.argument('experiment_file', type=Path_or_Subcommand(exists=True, dir_okay=False))
@click.pass_context
def main(ctx, additional_info, batch, data_dir, file_format, funcats, geneids, ifot, ifot_ki, ifot_tf,
         name, taxon, non_zeros, nonzero_subgroup, experiment_file):
    """
    """
         # name, taxon, non_zeros, experiment_file):

    if ifot + ifot_ki + ifot_tf > 1:
        raise click.BadParameter('Cannot specify a combination of `iFOT`, `iFOT-KI` and `iFOT-TF`')

    # validate_subgroup(nonzero_subgroup, experiment_file)
    validate_configfile(experiment_file, nonzero_subgroup=nonzero_subgroup, batch=batch)


    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

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

    data_obj = Data(additional_info=additional_info, batch=batch, data_dir=data_dir, funcats=funcats,
                    geneids=geneids, ifot=ifot, ifot_ki=ifot_ki, ifot_tf=ifot_tf, name=name,
                    non_zeros=non_zeros, nonzero_subgroup=nonzero_subgroup, taxon=taxon,
                    experiment_file=experiment_file)

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
                          outpath=data_obj.outpath)

    X = data_obj.areas_log_shifted.copy()
    X[data_obj.mask] = np.NaN
    g = scatterplot(X, stat=stat,
                    colors_only=colors_only, shade_correlation=shade_correlation,
                    outname=outname+file_fmts[0])
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
def cluster(ctx, col_cluster, dbscan, figsize, gene_symbols, highlight_geneids, max_autoclusters,
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
    )

    g = result['clustermap']['clustergrid']
    extra_artists = result['clustermap']['extra_artists']
    missing_values = 'masked' if show_missing_values else 'unmasked'
    outname_func = partial(get_outname, name=data_obj.outpath_name,
                           taxon=data_obj.taxon, non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                           outpath=data_obj.outpath, missing_values=missing_values, batch=data_obj.batch_applied)

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
                            colors_only=data_obj.colors_only, outpath=data_obj.outpath,
    )
    file_fmts = ctx.obj['file_fmts']
    save_multiple(fig, outname, *file_fmts)



# @click.command()
# @click.option('--additional-info', type=click.Path(exists=True, dir_okay=False), default=None,
#               help='.ini file with metadata for isobaric data used for scatter and PCA plots')
# @click.option('--data-dir', type=click.Path(exists=True, file_okay=False),
#               default='./data/', show_default=True,
#               help='optional location to store and read e2g files')
# @click.option('--col-cluster/--no-col-cluster', default=True, is_flag=True, show_default=True,
#               help="Cluster columns")
# @click.option('--export-data', type=click.Choice(['None', 'all', 'area']), default=None,
#               help="""Export data table of the filtered list of gene products used for plotting""")
# @click.option('--shade-correlation/--no-shade-correlation', default=True, is_flag=True,
#               show_default=True, help="")
# @click.option('--colors-only', default=False, is_flag=True, show_default=True,
#               help="Only plot colors on correlationplot, no datapoints")
# @click.option('--gene-symbols', default=False, is_flag=True, show_default=True,
#               help="Show Gene Symbols on clustermap")
# @click.option('--geneids', type=click.Path(exists=True, dir_okay=False),
#               default=None, show_default=True,
#               help="""Optional list of geneids to subset by.
#               Should have 1 geneid per line. """)
# @click.option('--highlight-geneids', type=click.Path(exists=True, dir_okay=False),
#               default=None, show_default=True, multiple=True,
#               help="""Optional list of geneids to highlight by.
#               Should have 1 geneid per line. """)
# @click.option('--funcats', type=str, default=None, show_default=True,
#               help="""Optional gene subset based on funcat or funcats,
#               regular expression allowed. """)
# @click.option('--iFOT', default=False, show_default=True, is_flag=True,
#               help='Calculate iFOT (divide by total input per experiment)')
# @click.option('--plots', type=click.Choice(['none', 'scatter', 'cluster', 'pca', 'all']), multiple=True,
#               default=('all',), show_default=True)
# @click.option('-n', '--name', type=str, default='',
#               help='An optional name for the analysis that will place all results in a subfolder.')
# @click.option('--non-non_zeros', default=0, show_default=True,
#               help='Minimum number of non-non_zeros allowed across samples.')
# @click.option('--row-cluster/--no-row-cluster', default=True, is_flag=True, show_default=True,
#               help="Cluster rows")
# @click.option('--stat', type=click.Choice(['pearson', 'spearman']),
#               default='spearman', show_default=True)
# @click.option('--standard-scale', type=click.Choice(['None', '0', '1']),
#               default='None', show_default=True)
# @click.option('--show-metadata', default=False, show_default=True,
#               is_flag=True,
#               help="""Show metadata on clustermap if present""")
# @click.option('--taxon', type=click.Choice(['human', 'mouse', 'all']),
#               default='all', show_default=True)
# @click.option('--z-score', type=click.Choice(['None', '0', '1']),
#               default='0', show_default=True)
# @click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
# def main(additional_info, data_dir, colors_only, export_data, gene_symbols, geneids,
#          highlight_geneids, funcats, ifot, stat, taxon, plots, name, non_zeros, experiment_file,
#          col_cluster, row_cluster, standard_scale, show_metadata, shade_correlation, z_score):

#     analysis_name = get_file_name(experiment_file)
#     if analysis_name is None:
#         print('Error parsing configfile name.')
#         analysis_name = 'Unnamed'

#     # OUTPATH = os.path.join('../figures/', analysis_name)
#     # if not os.path.exists(OUTPATH):
#     #     os.mkdir(OUTPATH)
#     # if name:
#     #     OUTPATH = os.path.join(OUTPATH, name)
#     #     if not os.path.exists(OUTPATH):
#     #         os.mkdir(OUTPATH)

#     now = datetime.now()
#     context = click.get_current_context()
#     params = context.params

#     data_obj = Data(**params)

#     cf = 'correlatioplot_args_{}.json'.format(now.strftime('%Y_%m_%d_%H_%M_%S'))
#     with open(os.path.join(data_obj.outpath, cf), 'w') as f:
#         json.dump(params, f)


#     fname, ext = os.path.splitext(experiment_file)

#     geneid_subset=None
#     if geneids:
#         geneid_subset = parse_gid_file(geneids)
#         if len(geneid_subset) == 0:
#             print('Non geneids found in file {}'.format(geneids))

#     highlight_gids = None
#     highlight_gid_names = None
#     if highlight_geneids:
#         highlight_gids = list()
#         highlight_gid_names = list()

#         for ix, h_gid in enumerate(highlight_geneids):
#             highlight_gid = parse_gid_file(h_gid)
#             highlight_gids.append(highlight_gid)
#             h_gid_name = get_file_name(h_gid)
#             if h_gid_name:
#                 highlight_gid_names.append(h_gid_name)
#             else:
#                 highlight_gid_names.append(ix)

#         if len(highlight_gids) == 0:
#             print('Non geneids found in file {}'.format(highlight_geneids))



#     data = read_config(experiment_file)
#     if additional_info:
#         labeled_meta = read_config(additional_info, enforce=False)
#     else:
#         labeled_meta = None

#     if len(data) == 0:
#         raise ValueError('No items in configfile.')

#     if z_score == 'None':
#         z_score = None
#     else:
#         z_score = int(z_score)

#     # experiment_file
#     run(data_obj)

#     # run(data, NON_ZEROS=non_zeros, stat=stat, taxon=taxon, data_dir=data_dir, OUTPATH=OUTPATH,
#     #     funcats=funcats, geneid_subset=geneid_subset, highlight_gids=highlight_gids, plots=plots,
#     #     highlight_gid_names=highlight_gid_names, colors_only=colors_only, gene_symbols=gene_symbols,
#     #     col_cluster=col_cluster, row_cluster=row_cluster, show_metadata=show_metadata,
#     #     shade_correlation=shade_correlation, z_score=z_score, labeled_meta=labeled_meta
#     # )


if __name__ == '__main__':

    main(obj={})
