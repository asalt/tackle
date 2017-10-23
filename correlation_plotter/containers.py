import sys
import os
import re
import json
from datetime import datetime
import operator as op
from collections import OrderedDict
from functools import partial
from warnings import warn

import numpy as np
import pandas as pd
from scipy import stats
from seaborn.matrix import ClusterGrid

from utils import *


idx = pd.IndexSlice

TAXON_MAPPER = {'human': 9606,
                'mouse': 10090}

def join_and_create_path(*strings, verbose=True):
    """Joins strings and returns resulting path.
    Creates path if does not exist
    """
    path = os.path.join(*strings)
    if not os.path.exists(path):
        os.mkdir(path)
        if verbose:
            print('Created new directory: {}'.format(os.path.abspath(path)))
    return path


class Data:

    def __init__(self, additional_info=None, col_cluster=True,
                 colors_only=False, data_dir='./data',
                 base_dir='./results',
                 experiment_file=None, funcats=None,
                 gene_symbols=False, geneids=None,
                 highlight_geneids=None,
                 name=None, non_zeros=0,
                 nonzero_subgroup=None,
                 plots=('all', ),
                 row_cluster=True, shade_correlation=True,
                 show_metadata=False,
                 standard_scale='None',
                 stat='pearson',
                 taxon='all',
                 z_score='0',
                 export_data=None,
                 ifot=False, ifot_ki=False, ifot_tf=False):
        "docstring"

        if experiment_file is None:
            raise ValueError('Must specify valid experiment_file')

        self.additional_info   = additional_info
        self.col_cluster       = col_cluster
        self.colors_only       = colors_only
        self.data_dir          = data_dir
        self.experiment_file   = experiment_file
        self.funcats           = funcats
        self.gene_symbols      = gene_symbols
        self.geneids           = geneids
        self.highlight_geneids = highlight_geneids
        self.non_zeros         = non_zeros
        self.nonzero_subgroup  = nonzero_subgroup
        self.plots             = plots
        self.row_cluster       = row_cluster
        self.shade_correlation = shade_correlation
        self.show_metadata     = show_metadata
        self.stat              = stat
        self.taxon             = taxon
        self.standard_scale    = self.clean_input(standard_scale)
        self.z_score           = self.clean_input(z_score)
        self.export_data       = None if export_data == 'None' else export_data
        self.ifot              = ifot
        self.ifot_ki           = ifot_ki
        self.ifot_tf           = ifot_tf
        self.base_dir          = base_dir

        self.outpath           = None
        self.analysis_name     = None


        self.set_analysis_name(experiment_file)
        self.set_outpath(base_dir, self.analysis_name, name)
        self.outpath_name = os.path.split(self.outpath)[-1]

        self.geneid_subset = None
        self.set_geneid_subset(geneids)

        self.highlight_gids, self.highlight_gid_names = None, None
        self.set_highlight_gids(highlight_geneids)

        self.config = read_config(experiment_file)
        if len(self.config) == 0:
            raise ValueError('No items in configfile.')
        self.labeled_meta = None
        if additional_info:
            self.labeled_meta = read_config(additional_info, enforce=False)

        self.panel = None
        self.load_data()

        self._gid_symbol = None

        self.panel_filtered = None
        self.filter_data()

        # self.ibaqs, self.ibaqs_log, self.ibaqs_log_shifted = (None, ) * 3
        self._areas, self._areas_log, self._areas_log_shifted = (None, ) * 3

        # self.perform_data_export()

    @property
    def areas(self):
        if self._areas is None:
            self.set_area_dfs()
        return self._areas

    @property
    def areas_log(self):
        if self._areas_log is None:
            self.set_area_dfs()
        return self._areas_log

    @property
    def areas_log_shifted(self):
        if self._areas_log_shifted is None:
            self.set_area_dfs()
        return self._areas_log_shifted

    @property
    def mask(self):
        if self._mask is None:
            self.set_area_dfs()
        return self._mask


    @staticmethod
    def clean_input(raw):
        if raw is None:
            return raw
        if raw.isdigit():
            return int(raw)
        elif raw == 'None':
            return None
        elif raw is None:
            return None
        else:  # invalid input
            warn('''Invalid input for z_score: `{}`.
            Setting to `None`.
            Choose from `None`, `0` or `1`.
            '''.format(raw))
            return None


    def set_outpath(self, path, *args):
        outpath = join_and_create_path(path)
        for arg in args:
            if arg:
                outpath = join_and_create_path(outpath, arg)
        self.outpath = outpath

    def set_analysis_name(self, experiment_file):
        ini_file = os.path.basename(experiment_file)
        analysis_name = os.path.splitext(ini_file)[0]
        self.analysis_name = analysis_name


    def set_geneid_subset(self, geneids):
        if geneids is None:
            self.geneid_subset = None
            return
        self.geneid_subset = parse_gid_file(geneids)
        if len(self.geneid_subset) == 0:
            warn('No geneids found in file {}'.format(geneids))

    def set_highlight_gids(self, highlight_geneids):
        if highlight_geneids is None:
            self.highlight_gids = None
            return
        highlight_gids = list()
        highlight_gid_names = list()
        for ix, h_gid in enumerate(highlight_geneids):
            highlight_gid = parse_gid_file(h_gid)
            if len(highlight_gid) == 0:
                warn('Non geneids found in file {}'.format(highlight_geneids))

            highlight_gids.append(highlight_gid)
            h_gid_name = get_file_name(h_gid)
            if h_gid_name:
                highlight_gid_names.append(h_gid_name)
            else:
                highlight_gid_names.append(ix)

        # self.highlight_gids = highlight_gid_names
        self.highlight_gids = highlight_gids
        self.highlight_gid_names = highlight_gid_names


    def _assign_labeled(self, record, exp, exps, name,
                        funcats=None, geneid_subset=None):
        labels = dict()
        for key, value in self.config[name].items():
            if key.isdigit():
                df = (exp.df[ exp.df.EXPLabelFLAG == int(key) ]
                      .set_index(['GeneID'])
                      # .pipe(filter_and_assign, value, funcats, geneid_subset)
                      .pipe(assign_cols, value)
                )
                labels[value] = df

        if self.labeled_meta and not all(v in self.labeled_meta for v in labels.keys()):
            warn('Mismatch between `labeled_meta` and input labels')

        else:
            pca = self.config.get('__PCA__')
            for key in self.labeled_meta.keys():
                if key not in labels:
                    continue
                newkey = '{}|{}'.format(name, key)
                exps[newkey] = labels[key]

            # exps = OrderedDict((key, exps[key])
            #                    for key in self.labeled_meta.keys())

            # re-assign config to labeled_meta
            # since this is a split of the each experiment ID
            # within each isobaric experiment
            # self.config = self.labeled_meta
            # self.config['__PCA__'] = pca

        return exps

    def load_data(self):

        exps = OrderedDict()
        config = self.config
        col_metadata = None

        for name, record in config.items():
            # print('Loading', record)
            if name.startswith('__'):
                continue
            labeltype = config[name].get('__LABELTYPE__', 'LF')
            recno = record.get('recno')
            runno = record.get('runno')
            searchno = record.get('searcno')
            exp = ispec.E2G(recno, runno, searchno, data_dir=self.data_dir)
            if len(exp) == 0:
                print('No data in {!r}'.format(exp))
                continue
            if labeltype == 'TMT' or labeltype == 'iTRAQ':
                exps = self._assign_labeled(record, exp, exps, name, self.funcats, self.geneid_subset)
            else:
                df = filter_and_assign(exp.df, name, self.funcats, self.geneid_subset,
                                       self.ifot, self.ifot_ki, self.ifot_tf
                )
                # df = assign_cols(exp.df, name)
                exps[name] = df.set_index(df.index.astype(int))

        # self.multi = pd.concat(exps.values(), keys=exps.keys())
        self.exps = exps
        stacked_data = [ df.stack() for df in exps.values() ]
        self.data = pd.concat( stacked_data, axis=1, keys=exps.keys() )
        # self.panel = pd.Panel(exps)
        for ax in ('GeneCapacity', 'GeneSymbol', 'GeneDescription',
                   'FunCats', 'TaxonID'):
            fillna_meta(self.data, ax)

        if self.additional_info:
            labeled_meta = read_config(self.additional_info, enforce=False)
            additional_metadata = parse_metadata(labeled_meta)
            metadata_dict = dict()
            for col in self.data.columns:
                exp, label = col.split('|')
                try:
                    metadata = additional_metadata[label]
                except KeyError:
                    warn('Mismatch between `labeled_meta` and input labels')
                    continue
                metadata_dict[col] = metadata
            col_metadata = pd.DataFrame.from_dict(metadata_dict)

        else:
            col_metadata = parse_metadata(self.config)

        self.col_metadata = col_metadata

    @property
    def gid_symbol(self):
        if self._gid_symbol is None:
            sel = self.data.loc[ idx[:, 'GeneSymbol'], self.data.columns[0] ]
            sel.index = sel.index.droplevel(1)
            self._gid_symbol = sel.to_dict()
            # self._gid_symbol = self.data.loc[ idx[:, 'GeneSymbol'], self.data.columns[0] ]
            # self._gid_symbol = self.panel.iloc[0]['GeneSymbol'].to_dict()
        return self._gid_symbol

    def filter_data(self):
        dummy_filter = lambda x, *args, **kwargs: x
        taxon_filter = TAXON_MAPPER.get(self.taxon)
        if taxon_filter is None:
            filter_func = dummy_filter
        else:
            filter_func = partial(filter_taxon, taxon=taxon_filter)

        df_filtered = (self.data.pipe(filter_observations, 'iBAQ_dstrAdj',
                                      self.non_zeros, self.nonzero_subgroup, self.col_metadata)
                       .pipe(filter_sra, SRA='S')
                       .pipe(filter_func)
        )
        df_filtered.index.names = ['GeneID', 'Metric']
        self.df_filtered = df_filtered


    def set_area_dfs(self):
        # self.areas = self.panel_filtered.minor_xs('iBAQ_dstrAdj').astype(float)
        self._areas = self.df_filtered.loc[ idx[:, 'iBAQ_dstrAdj'], : ]
        self._areas.index = self._areas.index.droplevel(1)  # don't need the second index
        self._mask = self._areas.applymap(np.isnan)
        if len(self.areas) == 0:
            raise ValueError('No data')

        # if self.ifot is True:
        #     sum_ = self._areas.sum(0)
        #     self._areas = self._areas / sum_
        #     # self.areas /= self.areas.sum(0)

        # elif self.ifot_ki is True:
        #     gids = self.df_filtered.pipe(filter_funcats, 'KI').index.get_level_values(0)
        #     sum_ = self._areas.loc[gids].sum(0)
        #     self._areas = self._areas / sum_

        # elif self.ifot_tf is True:
        #     gids = self.df_filtered.pipe(filter_funcats, 'TF')
        #     sum_ = self._areas.loc[gids].sum(0)
        #     self._areas = self._areas / sum_

        gids = set(self._areas.index)
        if self.funcats:
            gids &= set(self.df_filtered.pipe(filter_funcats, self.funcats).index.get_level_values(0))
        if self.geneid_subset:
            gids &= set(self.geneid_subset)

        gids = tuple(gids)

        self._areas_log = np.log10(self._areas.fillna(0)+1e-10)
        # fillna with the mean value. This prevents skewing of normalization such as
        # z score. The NAN values are held in the self.mask dataframe
        # self._areas_log = np.log10(self._areas.T.fillna(self._areas.mean(axis=1)).T + 1e-8)
        minval = self._areas_log.min().min()
        shift_val = np.ceil(np.abs(minval))

        self._areas_log_shifted = self._areas_log + shift_val

    def make_plot(self, pltname):
        if 'all' in self.plots:
            return True
        if pltname in self.plots:
            return True
        return False

    def perform_data_export(self, level='all'):
        # if self.export_data is None:
        #     return

        fname = '{}_data_{}_{}_more_zeros.tab'.format(level,
                                                      self.outpath_name,
                                                      self.non_zeros)
        outname = os.path.abspath(os.path.join(self.outpath, fname))
        # if self.export_data == 'all':
        if level == 'all':
            self.df_filtered.to_csv(outname, sep='\t')
        elif level == 'area':
            self.areas_log_shifted.to_csv(outname, sep='\t')
        print('Exported', outname)


class MyClusterGrid(ClusterGrid):

    def __init__(self, *args, heatmap_height_ratio=.8,
                 dendrogram_width_ratio=.16, **kwargs):
        self.heatmap_height_ratio = heatmap_height_ratio
        self.dendrogram_width_ratio = dendrogram_width_ratio
        super().__init__(*args, **kwargs)

    def dim_ratios(self, side_colors, axis, figsize, side_colors_ratio=0.05):
        # need to adjust the heatmap height ratio for long figures
        # such that it fills up more of the given room
        ratios = ClusterGrid.dim_ratios(self, side_colors, axis, figsize, side_colors_ratio=0.05)

        if axis == 0:  # calculating height ratios
            ratios[-1] = self.heatmap_height_ratio

            if self.dendrogram_width_ratio: #
                ratios[0] = self.dendrogram_width_ratio

        elif axis == 1 and self.dendrogram_width_ratio:  # calculating width ratios
            ratios[1] = self.dendrogram_width_ratio
            ratios[0] = self.dendrogram_width_ratio * .4

        # print(axis, ':', ratios)
        return ratios
