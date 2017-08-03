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

rc = {'font.family': 'serif',
      'font.serif': ['Times', 'Palatino', 'serif']}
sb.set_context('paper')
sb.set_style('white', rc)

__version__ = '0.2'

from bcmproteomics_ext import ispec
sb.set_context('notebook', font_scale=1.4)

from scatterplot import scatterplot
from clusterplot import clusterplot
from pcaplot import pcaplot
from utils import *
from containers import Data
# from cluster_to_plotly import cluster_to_plotly

sys.setrecursionlimit(10000)

TAXON_MAPPER = {'human': 9606,
                'mouse': 10090}


# def run(records, NON_ZEROS=0, stat='spearman', taxon='all', data_dir=None, OUTPATH='../figures',
#         funcats=None, geneid_subset=None, highlight_gids=None, highlight_gid_names=None,
#         colors_only=False, gene_symbols=False, col_cluster=True, row_cluster=True,
#         show_metadata=False, shade_correlation=True, z_score=0, plots=('all',),
#         labeled_meta=None):
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
        save_multiple(g, outname, '.png', bbox_extra_artists=extra_artists, bbox_inches=bbox_inches)

    # if 'pca' in plots or 'all' in plots:
    if data_obj.make_plot('pca'):
        fig, ax = pcaplot(data_obj.areas_log_shifted, data_obj.config, col_data = data_obj.col_metadata)
        outname = get_outname('pcaplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
                              non_zeros=data_obj.non_zeros,
                              colors_only=data_obj.colors_only, outpath=data_obj.outpath,
        )
        save_multiple(fig, outname, '.png')

@click.command()
@click.option('--additional-info', type=click.Path(exists=True, dir_okay=False), default=None,
              help='.ini file with metadata for isobaric data used for scatter and PCA plots')
@click.option('--data-dir', type=click.Path(exists=True, file_okay=False),
              default='../data/raw', show_default=True,
              help='optional location to store and read e2g files')
@click.option('--col-cluster/--no-col-cluster', default=True, is_flag=True, show_default=True,
              help="Cluster columns")
@click.option('--export-data', type=click.Choice(['None', 'all', 'area']), default=None,
              help="""Export data table of the filtered list of gene products used for plotting""")
@click.option('--shade-correlation/--no-shade-correlation', default=True, is_flag=True,
              show_default=True, help="")
@click.option('--colors-only', default=False, is_flag=True, show_default=True,
              help="Only plot colors on correlationplot, no datapoints")
@click.option('--gene-symbols', default=False, is_flag=True, show_default=True,
              help="Show Gene Symbols on clustermap")
@click.option('--geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True,
              help="""Optional list of geneids to subset by.
              Should have 1 geneid per line. """)
@click.option('--highlight-geneids', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True, multiple=True,
              help="""Optional list of geneids to highlight by.
              Should have 1 geneid per line. """)
@click.option('--funcats', type=str, default=None, show_default=True,
              help="""Optional gene subset based on funcat or funcats,
              regular expression allowed. """)
@click.option('--iFOT', default=False, show_default=True, is_flag=True,
              help='Calculate iFOT (divide by total input per experiment)')
@click.option('--plots', type=click.Choice(['none', 'scatter', 'cluster', 'pca', 'all']), multiple=True,
              default=('all',), show_default=True)
@click.option('-n', '--name', type=str, default='',
              help='An optional name for the analysis that will place all results in a subfolder.')
@click.option('--non-zeros', default=0, show_default=True,
              help='Minimum number of non-non_zeros allowed across samples.')
@click.option('--row-cluster/--no-row-cluster', default=True, is_flag=True, show_default=True,
              help="Cluster rows")
@click.option('--stat', type=click.Choice(['pearson', 'spearman']),
              default='spearman', show_default=True)
@click.option('--standard-scale', type=click.Choice(['None', '0', '1']),
              default='None', show_default=True)
@click.option('--show-metadata', default=False, show_default=True,
              is_flag=True,
              help="""Show metadata on clustermap if present""")
@click.option('--taxon', type=click.Choice(['human', 'mouse', 'all']),
              default='all', show_default=True)
@click.option('--z-score', type=click.Choice(['None', '0', '1']),
              default='0', show_default=True)
@click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
def main(additional_info, data_dir, colors_only, export_data,
         gene_symbols, geneids, highlight_geneids, funcats, ifot,
         stat, taxon, plots, name, non_zeros, experiment_file, col_cluster, row_cluster,
         standard_scale,
         show_metadata, shade_correlation, z_score):

    analysis_name = get_file_name(experiment_file)
    if analysis_name is None:
        print('Error parsing configfile name.')
        analysis_name = 'Unnamed'

    OUTPATH = os.path.join('../figures/', analysis_name)
    if not os.path.exists(OUTPATH):
        os.mkdir(OUTPATH)
    if name:
        OUTPATH = os.path.join(OUTPATH, name)
        if not os.path.exists(OUTPATH):
            os.mkdir(OUTPATH)

    now = datetime.now()
    context = click.get_current_context()
    params = context.params

    data_obj = Data(**params)

    cf = 'correlatioplot_args_{}.json'.format(now.strftime('%Y_%m_%d_%H_%M_%S'))
    with open(os.path.join(data_obj.outpath, cf), 'w') as f:
        json.dump(params, f)


    fname, ext = os.path.splitext(experiment_file)

    geneid_subset=None
    if geneids:
        geneid_subset = parse_gid_file(geneids)
        if len(geneid_subset) == 0:
            print('Non geneids found in file {}'.format(geneids))

    highlight_gids = None
    highlight_gid_names = None
    if highlight_geneids:
        highlight_gids = list()
        highlight_gid_names = list()

        for ix, h_gid in enumerate(highlight_geneids):
            highlight_gid = parse_gid_file(h_gid)
            highlight_gids.append(highlight_gid)
            h_gid_name = get_file_name(h_gid)
            if h_gid_name:
                highlight_gid_names.append(h_gid_name)
            else:
                highlight_gid_names.append(ix)

        if len(highlight_gids) == 0:
            print('Non geneids found in file {}'.format(highlight_geneids))



    data = read_config(experiment_file)
    if additional_info:
        labeled_meta = read_config(additional_info, enforce=False)
    else:
        labeled_meta = None

    if len(data) == 0:
        raise ValueError('No items in configfile.')

    if z_score == 'None':
        z_score = None
    else:
        z_score = int(z_score)

    # experiment_file
    run(data_obj)

    # run(data, NON_ZEROS=non_zeros, stat=stat, taxon=taxon, data_dir=data_dir, OUTPATH=OUTPATH,
    #     funcats=funcats, geneid_subset=geneid_subset, highlight_gids=highlight_gids, plots=plots,
    #     highlight_gid_names=highlight_gid_names, colors_only=colors_only, gene_symbols=gene_symbols,
    #     col_cluster=col_cluster, row_cluster=row_cluster, show_metadata=show_metadata,
    #     shade_correlation=shade_correlation, z_score=z_score, labeled_meta=labeled_meta
    # )


if __name__ == '__main__':

    main()


    # EXPS = [
    #     dict(
    #         recno=32918,
    #         runno=1,
    #         searchno=1),
    #     dict(
    #         recno=32919,
    #         runno=1,
    #         searchno=1
    #     ),

    # ]
    # NON_ZEROS = 0 # a number or 'max'
    # stat='spearman'

    # main(EXPS, NON_ZEROS, stat)

    # experiment_file = './configtest.ini'
    # fname, ext = os.path.splitext(experiment_file)
    # geneids = './geneid_subset.txt'
    # highlight_geneids = './highlight_subset.txt'
