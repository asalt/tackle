import sys
import os
import re
import json
from datetime import datetime
import operator as op
from collections import OrderedDict
from functools import partial
from warnings import warn
from matplotlib import cm, gridspec

import numpy as np
import pandas as pd
from scipy import stats
from seaborn.matrix import ClusterGrid
from seaborn.matrix import _HeatMapper as HeatMapper
from seaborn.matrix import _matrix_mask

from .utils import *


idx = pd.IndexSlice

TAXON_MAPPER = {'human': 9606,
                'mouse': 10090}


LABEL_MAPPER = {'none': 0,  # hard coded number IDs for labels
                '126': 1260,
                '127C': 1270,
                '127N': 1271,
                '128C': 1280,
                '128N': 1281,
                '129C': 1290,
                '129N': 1291,
                '130C': 1300,
                '130N': 1301,
                '131': 1310,
                '113': 113,
                '114': 114,
                '115': 115,
                '116': 116,
                '117': 117,
                '118': 118,
                '119': 119,
                '121': 121,
}

def maybe_int(x):
    try:
        return int(x)
    except ValueError:
        warn('Value {} cannot be converted to int'.format(x))
        return x

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

class GeneMapper:
    def __init__(self):
        self.file = os.path.join( os.path.split(os.path.abspath(__file__))[0],
                                  'data', 'genetable_hsmmgg.tab'
        )
        self._df = None
        self._symbol = None
        self._funcat = None
        self._description = None

    @property
    def df(self):
        if self._df is None:
            self._df = pd.read_table(self.file, index_col='GeneID')
        return self._df

    @property
    def symbol(self):
        if self._symbol is None:
            self._symbol = self.df['GeneSymbol'].to_dict()
        return self._symbol

    @property
    def funcat(self):
        if self._funcat is None:
            self._funcat = self.df['FunCats'].to_dict()
        return self._funcat

    @property
    def description(self):
        if self._description is None:
            self._description = self.df['GeneDescription'].to_dict()
        return self._description

_genemapper = GeneMapper()

class Data:

    def __init__(self, additional_info=None, batch=None,
                 batch_nonparametric=False,
                 batch_noimputation=False,
                 covariate=None,
                 col_cluster=True,
                 colors_only=False, data_dir='./data',
                 base_dir='./results',
                 experiment_file=None, funcats=None,
                 funcats_inverse=None,
                 gene_symbols=False, geneids=None,
                 group=None, pairs=None,
                 highlight_geneids=None,
                 ignore_geneids=None,
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
                 ifot=False, ifot_ki=False, ifot_tf=False, median=False,
                 set_outpath=True,
                 metrics=False, metrics_after_filter=True,

    ):
        "docstring"

        if experiment_file is None:
            raise ValueError('Must specify valid experiment_file')

        self.additional_info      = additional_info
        self.batch                = batch
        self.batch_nonparametric  = batch_nonparametric
        self.batch_noimputation   = batch_noimputation
        self.covariate            = covariate
        self.col_cluster          = col_cluster
        self.colors_only          = colors_only
        self.data_dir             = data_dir
        self.experiment_file      = experiment_file
        self.funcats              = funcats
        self.funcats_inverse      = funcats_inverse
        self.gene_symbols         = gene_symbols
        self.geneids              = geneids
        self.ignore_geneids       = ignore_geneids
        self.group                = group
        self.pairs                = pairs
        self.highlight_geneids    = highlight_geneids
        self.non_zeros            = non_zeros
        self.nonzero_subgroup     = nonzero_subgroup
        self.plots                = plots
        self.row_cluster          = row_cluster
        self.shade_correlation    = shade_correlation
        self.show_metadata        = show_metadata
        self.stat                 = stat
        self.taxon                = taxon
        self.standard_scale       = self.clean_input(standard_scale)
        self.z_score              = self.clean_input(z_score)
        self.export_data          = None if export_data == 'None' else export_data
        self.ifot                 = ifot
        self.ifot_ki              = ifot_ki
        self.ifot_tf              = ifot_tf
        self.median               = median
        self.base_dir             = base_dir
        self.metrics              = metrics
        self.metrics_after_filter = metrics_after_filter
        self._metric_values       = None

        self.outpath              = None
        self.analysis_name        = None

        self.normed               = False
        self.batch_applied        = None  # set to the batch (only) upon successful batch correction
        self.gid_funcat_mapping   = None  # loaded with data

        self.set_analysis_name(experiment_file)
        if set_outpath:
            self.set_outpath(base_dir, self.analysis_name, name)
            self.outpath_name = os.path.split(self.outpath)[-1]
        else:
            self.outpath = None
            self.outpath_name = None


        # self.geneid_subset = None
        self.geneid_subset = self.set_geneid_subset(geneids)
        self.ignore_geneid_subset = self.set_geneid_subset(ignore_geneids)

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

        self._qvalues = None
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

    @property
    def zeros(self):
        if self._zeros is None:
            self.set_area_dfs()
        return self._zeros

    @property
    def qvalues(self):
        if self._qvalues is None:
            self.calc_qvals()
        return self._qvalues

    @property
    def metric_values(self):
        return self._metric_values


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
            return None
            # self.geneid_subset = None
            # return
        # self.geneid_subset = parse_gid_file(geneids)
        geneid_subset = parse_gid_file(geneids)
        if len(geneid_subset) == 0:
            warn('No geneids found in file {}'.format(geneids))
        return geneid_subset

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
                if self.metrics and not self.metrics_after_filter:
                    self._update_metrics(df, name)

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

    def _update_metrics(self, df, name):
        if self.metric_values is None:
            self._metric_values = DefaultOrderedDict(lambda : defaultdict(None) )

        sra = df.SRA.value_counts().to_dict()
        gpg = df.GPGroup.nunique()
        psms = {'Total': df.PSMs.sum(),
                'u2g':   df.PSMs_u2g.sum(),
        }
        peptides = {'Total':     df.PeptideCount.sum(),
                   'u2g':        df.PeptideCount_u2g.sum(),
                   'Strict':     df.PeptideCount_S.sum(),
                   'Strict_u2g': df.PeptideCount_S_u2g.sum(),
        }

        self._metric_values[name]['SRA']      = sra
        self._metric_values[name]['GPGroups'] = gpg
        self._metric_values[name]['PSMs']     = psms
        self._metric_values[name]['Peptides'] = peptides
        self._metric_values[name]['Area']     = df.AreaSum_dstrAdj.where(lambda x: x > 0 ).dropna().values
        # self._metric_values[name]['GeneIDs']  = exp.df.GeneID.unique()



    def load_data(self):

        exps = OrderedDict()
        gid_funcat_mapping = dict()
        config = self.config
        col_metadata = None

        for name, record in config.items():
            # print('Loading', record)
            if name.startswith('__'):
                continue
            labeltype = config[name].get('__LABELTYPE__', 'LF') # depreciated
            recno = record.get('recno')
            runno = record.get('runno')
            searchno = record.get('searcno')
            label = record.get('label')
            labelquery = LABEL_MAPPER.get(label, 0)
            exp = ispec.E2G(recno, runno, searchno, data_dir=self.data_dir)
            df = exp.df.query('EXPLabelFLAG==@labelquery').copy()

            if not df.index.name == 'GeneID':
                df.index = df.GeneID

            if df.empty:
                warn('No data for {!r}, skipping'.format(exp))
                continue
                # raise ValueError('No data for {!r}'.format(exp))


            if self.metrics and not self.metrics_after_filter:
                self._update_metrics(df, name)

            # exp.df['GeneID'] = exp.df['GeneID'].astype(int)
            df['GeneID'] = df['GeneID'].apply(maybe_int)
            funcats_dict = df.drop_duplicates('GeneID').set_index('GeneID')['FunCats'].to_dict()
            gid_funcat_mapping.update(funcats_dict)

            if labeltype == 'TMT' or labeltype == 'iTRAQ': # depreciated
                exps = self._assign_labeled(record, exp, exps, name, self.funcats, self.geneid_subset)
            else:
                df = filter_and_assign(df, name, self.funcats, self.funcats_inverse,
                                       self.geneid_subset, self.ignore_geneid_subset, self.ifot,
                                       self.ifot_ki, self.ifot_tf, self.median)
                # df = assign_cols(exp.df, name)
                if self.metrics and self.metrics_after_filter:
                    self._update_metrics(df, name)
                # exps[name] = df.set_index(df.index.astype(int))
                df.index = [maybe_int(x) for x in df.index]
                exps[name] = df

        self.gid_funcat_mapping = gid_funcat_mapping

        # self.multi = pd.concat(exps.values(), keys=exps.keys())
        self.exps = exps
        # _cols = ['TaxonID', 'IDSet', 'GeneSymbol', 'iBAQ_dstrAdj', 'FunCats', 'SRA']
        # stacked_data = [ df[_cols].stack() for df in exps.values() ]
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
            try:
                sel = self.data.loc[ idx[:, 'GeneSymbol'], self.data.columns[0] ]
                sel.index = sel.index.droplevel(1)
                self._gid_symbol = sel.to_dict()
            except KeyError:
                return dict()
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


        batch_info = self.config.get('__batch__')
        if batch_info:
            batch = batch_info.get('batch')
            metadata = self.col_metadata.T  # rows are experiments, cols are metadata
            areas = self._areas.copy()
            for ix, g in metadata.groupby(batch):
                meanvals = areas[g.index].mean(1)
                areas.loc[:, g.index] = (self._areas[g.index]
                                       .div(meanvals.fillna(0)+1e-20,
                                            axis='index')
                )

            self._areas = areas


        self._mask = self._areas.applymap(np.isnan)
        self._zeros = self._areas == 0
        if len(self.areas) == 0:
            raise ValueError('No data')

        gids = set(self._areas.index)
        if self.funcats:
            gids &= set(self.df_filtered.pipe(filter_funcats, self.funcats).index.get_level_values(0))
        if self.geneid_subset:
            gids &= set(self.geneid_subset)

        gids = tuple(gids)
        self.minval = self._areas.replace(0, np.NAN).stack().dropna().min()

        # self._areas_log = np.log10(self._areas.fillna(0)+1e-10)
        self._areas_log = np.log10(self._areas.replace(0, np.NAN).fillna(self.minval/2))
        self._areas_log.index.name = 'GeneID'
        # fillna with the mean value. This prevents skewing of normalization such as
        # z score. The NAN values are held in the self.mask dataframe
        # self._areas_log = np.log10(self._areas.T.fillna(self._areas.mean(axis=1)).T + 1e-8)

        # if norm_info is not None:  # do not shift the values
        #     self._areas_log_shifted = self._areas_log
        #     return

        minval = self._areas_log.min().min()
        shift_val = np.ceil(np.abs(minval))
        self.minval_log = minval

        self._areas_log_shifted = self._areas_log + shift_val
        self._areas_log_shifted.index.name = 'GeneID'


        # if specified, normalize by a specified control group
        norm_info = self.config.get('__norm__')
        if norm_info is not None:
            control = norm_info['control']
            group   = norm_info['group']
            label   = norm_info['label']
            metadata = self.col_metadata.T  # rows are experiments, cols are metadata
            areas = self._areas_log_shifted.copy()
            ctrl_exps = list()
            for ix, g in metadata.groupby(group):
                ctrl_exp = g[ g[label] == control ].index[0] # should only be one
                ctrl_exps.append(ctrl_exp)
                # to_normalize = list(set(g.index) - set([ctrl_exp]))
                to_normalize = g.index
                areas.loc[:, to_normalize] = (self._areas_log_shifted[to_normalize]
                                              .sub(self._areas_log_shifted[ctrl_exp]
                                                   .fillna(0) + 1e-20,
                                                   axis='index')
                )
            self._areas = (areas.drop(ctrl_exps, axis=1)
                           .where(lambda x: x!= 0)
                           .dropna(how='all')
            )
            self.minval_log = self._areas.replace(0, np.NAN).stack().dropna().min()
            finite = self._areas.pipe(np.isfinite)
            maxval_log = self._areas[ finite ].stack().dropna().max()
            self._areas_log_shifted = (self._areas.fillna(self.minval_log/2)
                                       .replace(np.inf, maxval_log*1.5)
            )
            sample_cols = [x for x in self.col_metadata.columns if x not in ctrl_exps]
            sample_ixs  = [x for x in self.col_metadata.index if x != label]
            self.col_metadata = self.col_metadata.loc[sample_ixs, sample_cols]
            self._mask = self.mask[sample_cols]
            self.normed = True

        if self.batch is not None:
            # try batch normalization via ComBat
            self.batch_normalize()

        # if self.group is not None:
        #     self.calc_qvals()


    def batch_normalize(self):

        try:
            from rpy2.rinterface import RRuntimeError
        except Exception as e:
            RRuntimeError = type('RRuntimeError', (Exception,), {})
        try:
            from rpy2.robjects import r
            from rpy2.robjects.packages import importr
            sva = importr('sva')
        except ModuleNotFoundError as e:
            print('Failure in rpy2 import and load', e, sep='\n')
            return
        except RRuntimeError as e:
            print('sva is not installed', e, sep='\n')
            return

        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        grdevices = importr('grDevices')

        pheno = self.col_metadata.T
        r.assign('pheno', pheno)

        if self.covariate is not None:
            r('mod  <- model.matrix(~as.factor({}), pheno)'.format(self.covariate))
            mod = r['mod']
            self.batch_applied = self.batch + '_Cov_{}'.format(self.group)
        else:
            mod = r('model.matrix(~1, pheno)')
            self.batch_applied = self.batch + '_noCov'
        # r.assign('batch', 'pheno${}'.format(self.batch))
        batch = pheno[self.batch]
        # res = sva.ComBat(dat=self._areas_log_shifted.fillna(0), batch=batch,
        #                  mod=mod, par_prior=True, mean_only=False)
        if not self.batch_nonparametric:
            plot_prior = True
            outname = get_outname('Combat_prior_plots', name=self.outpath_name, taxon=self.taxon,
                                  non_zeros=self.non_zeros, colors_only=self.colors_only,
                                  batch=self.batch_applied,
                                  batch_method = 'parametric' if not self.batch_nonparametric else 'nonparametric',
            outpath=self.outpath)
            grdevices.png(file=outname+'.png', width=5, height=5, units='in', res=300)
        else:
            plot_prior = False
        res = sva.ComBat(dat=self._areas_log_shifted.fillna(0).values, batch=batch,
                         mod=mod, par_prior=not self.batch_nonparametric, mean_only=False, prior_plots=plot_prior)

        if plot_prior:
            grdevices.dev_off()

        df = pd.DataFrame(index=self.areas_log_shifted.index,
                          columns=self.areas_log_shifted.columns,
                          data=pandas2ri.ri2py(res)
        )
        # df = pandas2ri.ri2py(res)
        nas = sum(df.isnull().any(1))
        if nas > 0:
            print('{} Gene Product(s) became NAN after batch normalization, dropping'.format(nas))

        # df.index = df.index.astype(int)
        df.index = self._areas_log_shifted.index
        # df.index = [maybe_int(x) for x in df.index]
        df.columns = self._areas_log_shifted.columns
        # reassign mask - ComBat can impute some NA values
        # TODO: resolve this for normed data
        if not self.normed:
            if not self.batch_noimputation:  # else leave old mask
                thresh = self.areas_log_shifted[ (~self.mask) & (self._areas_log_shifted > 0)].min().min() * .9
                new_mask = (df[ self.mask ] < thresh)
                new_mask.columns = self._areas_log_shifted.columns
                self._mask = new_mask


        self._areas_log_shifted = df.dropna(how='any')
        # self._areas_log_shifted.columns = self._areas_log.columns  # r changes certain characters in column names
        self._areas_log_shifted.columns = self._areas.columns  # r changes certain characters in column names
        # self._areas_log_shifted.index = self._areas_log_shifted.index.astype(int)  #r converts rownames to str

        self._areas_log_shifted.index = [maybe_int(x) for x in self._areas_log_shifted.index]

        self._areas_log_shifted.index.name = 'GeneID'

    def calc_qvals(self):
        try: ModuleNotFoundError
        except: ModuleNotFoundError = type('ModuleNotFoundError', (Exception,), {})
        try:
            from rpy2.rinterface import RRuntimeError
        except Exception as e:
            RRuntimeError = type('RRuntimeError', (Exception,), {})
        try:
            from rpy2.robjects import r
            from rpy2.robjects.packages import importr
            sva = importr('sva')
        except ModuleNotFoundError as e:
            print('Failure in rpy2 import and load', e, sep='\n')
            return
        except RRuntimeError as e:
            print('sva is not installed', e, sep='\n')
            return
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        r_source = r['source']
        r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                              'R', 'pvalue_cov.R')
        r_source(r_file)

        pheno = self.col_metadata.T
        r.assign('pheno', pheno)
        r('mod0 <- model.matrix(~1, pheno)')
        r('mod  <- model.matrix(~as.factor({}), pheno)'.format(self.group))
        r.assign('edata', self.areas_log_shifted.fillna(0))

        if self.covariate is not None:
            ncov = pheno[self.covariate].nunique()
            r.assign('ncov', pheno[self.covariate].nunique())
        else:
            r.assign('ncov', 0)
            ncov = 0

        if not self.pairs:
            pvalues = r('pvalue.batch(as.matrix(edata), mod, mod0, ncov)')
        else: # ttest rel
            from .ttest_rel_cov import ttest_rel_cov
            groups = self.col_metadata.loc[self.group].unique()
            # only have t-test implemented here
            if len(groups) > 2:
                raise NotImplementedError('Only have support for 2 groups')
            group0, group1 = groups
            meta_ordered = self.col_metadata.sort_index(by=[self.group, self.pairs], axis=1)
            # better way to do this?
            cols0 = (meta_ordered.T[self.group] == group0).apply(lambda x: x if x else np.nan).dropna().index
            cols1 = (meta_ordered.T[self.group] == group1).apply(lambda x: x if x else np.nan).dropna().index
            # self.
            # t test on every gene set
            def ttest_func(row):
                t, p = ttest_rel_cov(row[cols0], row[cols1], ncov=ncov)
                return p
            pvalues = self.areas_log_shifted.apply(ttest_func, axis=1)



        p_adjust = r['p.adjust']

        # pvalues = f_pvalue(self.areas_log_shifted.fillna(0), mod, mod0)
        qvalues = p_adjust(pvalues, method='BH')
        qvalues = pd.DataFrame(index=self.areas_log_shifted.index,
                               data=np.array([pvalues, pandas2ri.ri2py(qvalues)]).T,
                               columns=['pValue', 'qValue']
        ).sort_values(by='pValue')
        # qvalues.name = 'q_value'

        self._qvalues = qvalues

    def make_plot(self, pltname):
        if 'all' in self.plots:
            return True
        if pltname in self.plots:
            return True
        return False

    def perform_data_export(self, level='all', genesymbols=False):
        # if self.export_data is None:
        #     return

        # fname = '{}_data_{}_{}_more_zeros.tab'.format(level,
        #                                               self.outpath_name,
        #                                               self.non_zeros)

        outname = get_outname('data_{}'.format(level), name=self.outpath_name, taxon=self.taxon,
                              non_zeros=self.non_zeros, colors_only=self.colors_only,
                              batch=self.batch_applied,
                              batch_method = 'parametric' if not self.batch_nonparametric else 'nonparametric',
                              outpath=self.outpath) + '.tab'

        # outname = os.path.abspath(os.path.join(self.outpath, fname))
        # if self.export_data == 'all':
        if level == 'all':
            self.df_filtered.to_csv(outname, sep='\t')
        elif level == 'area':
            export = self.areas_log_shifted.copy()
            export[self.mask] = np.NaN
            order = export.columns
            if genesymbols:
                # export['GeneSymbol'] = export.index.map(lambda x: self.gid_symbol.get(x, '?'))
                export['GeneSymbol'] = export.index.map(lambda x: self.gid_symbol.get(x,
                                                                      _genemapper.symbol.get(x, '?')
                ))
                order = ['GeneSymbol']  # index column is GeneID, add GeneSymbol
                order += [x for x in export.columns if x not in order]
            export[order].to_csv(outname, sep='\t')
        print('Exported', outname)


from six import string_types
class MyHeatMapper(HeatMapper):

    def _determine_cmap_params(self, plot_data, vmin, vmax,
                               cmap, center, robust):
        """Use some heuristics to set good defaults for colorbar and range."""
        calc_data = plot_data.data[~np.isnan(plot_data.data)]
        if vmin is None:
            vmin = np.percentile(calc_data, 2) if robust else calc_data.min()
            # vmin = np.percentile(calc_data, 20) if robust else calc_data.min()
        if vmax is None:
            vmax = np.percentile(calc_data, 98) if robust else calc_data.max()
            # vmax = np.percentile(calc_data, 75) if robust else calc_data.max()
        self.vmin, self.vmax = vmin, vmax

        # Choose default colormaps if not provided
        if cmap is None:
            if center is None:
                self.cmap = cm.rocket
            else:
                self.cmap = cm.icefire
        elif isinstance(cmap, string_types):
            self.cmap = mpl.cm.get_cmap(cmap)
        elif isinstance(cmap, list):
            self.cmap = mpl.colors.ListedColormap(cmap)
        else:
            self.cmap = cmap

        # Recenter a divergent colormap
        if center is not None:
            vrange = max(vmax - center, center - vmin)
            normlize = mpl.colors.Normalize(center - vrange, center + vrange)
            cmin, cmax = normlize([vmin, vmax])
            cc = np.linspace(cmin, cmax, 256)
            self.cmap = mpl.colors.ListedColormap(self.cmap(cc))
        self.cmap.set_bad(color='gray')

sb.matrix._HeatMapper = MyHeatMapper

class MyClusterGrid(ClusterGrid):

    # def __init__(self, *args, heatmap_height_ratio=.8,
    #              dendrogram_width_ratio=.16, heatmap_width_ratio=.8,
    #              expected_size_dendrogram=1.0, expected_size_colors=0.25:
    #              **kwargs):
    def __init__(self, data, pivot_kws=None, z_score=None, standard_scale=None,
                 figsize=None, row_colors=None, col_colors=None, mask=None,
                 expected_size_dendrogram=1.0,
                 expected_size_colors=0.25):
        """Grid object for organizing clustered heatmap input on to axes"""

        if isinstance(data, pd.DataFrame):
            self.data = data
        else:
            self.data = pd.DataFrame(data)

        self.data2d = self.format_data(self.data, pivot_kws, z_score,
                                       standard_scale)

        self.mask = _matrix_mask(self.data2d, mask)

        self.expected_size_dendrogram = expected_size_dendrogram
        self.expected_size_side_colors = expected_size_colors

        if figsize is None:
            width, height = 10, 10
            figsize = (width, height)
        self.fig = plt.figure(figsize=figsize)

        self.row_colors, self.row_color_labels = \
            self._preprocess_colors(data, row_colors, axis=0)
        self.col_colors, self.col_color_labels = \
            self._preprocess_colors(data, col_colors, axis=1)

        width_ratios = self.dim_ratios(self.row_colors,
                                       figsize=figsize,
                                       axis=0)
        height_ratios = self.dim_ratios(self.col_colors,
                                        figsize=figsize,
                                        axis=1)
        nrows = 3 if self.col_colors is None else 4
        ncols = 3 if self.row_colors is None else 4
        self.gs = gridspec.GridSpec(nrows, ncols, wspace=0.01, hspace=0.01,
                                    width_ratios=width_ratios,
                                    height_ratios=height_ratios)

        self.ax_row_dendrogram = self.fig.add_subplot(self.gs[nrows - 1, 0:2])
        self.ax_col_dendrogram = self.fig.add_subplot(self.gs[0:2, ncols - 1])
        self.ax_row_dendrogram.set_axis_off()
        self.ax_col_dendrogram.set_axis_off()

        self.ax_row_colors = None
        self.ax_col_colors = None

        if self.row_colors is not None:
            self.ax_row_colors = self.fig.add_subplot(
                self.gs[nrows - 1, ncols - 2])
        if self.col_colors is not None:
            self.ax_col_colors = self.fig.add_subplot(
                self.gs[nrows - 2, ncols - 1])

        self.ax_heatmap = self.fig.add_subplot(self.gs[nrows - 1, ncols - 1])

        # colorbar for scale to left corner
        if self.col_colors is not None:
            cbar_max = 3
        else:
            cbar_max = 2

        self.cax = self.fig.add_subplot(self.gs[0:cbar_max, 0])

        self.dendrogram_row = None
        self.dendrogram_col = None
        # self.heatmap_height_ratio = heatmap_height_ratio
        # self.dendrogram_width_ratio = dendrogram_width_ratio
        # self.heatmap_width_ratio=heatmap_width_ratio

        # super().__init__(*args, **kwargs)
    def dim_ratios(self, side_colors, axis, figsize):
        """Get the proportions of the figure taken up by each axes
        """
        figdim = figsize[axis]

        expected_size_for_dendrogram = self.expected_size_dendrogram  # Inches
        expected_size_for_side_colors = self.expected_size_side_colors  # Inches

        # Get resizing proportion of this figure for the dendrogram and
        # colorbar, so only the heatmap gets bigger but the dendrogram stays
        # the same size.
        dendrogram = expected_size_for_dendrogram / figdim

        # add the colorbar
        colorbar_width = .8 * dendrogram
        colorbar_height = .2 * dendrogram
        if axis == 1:
            ratios = [colorbar_width, colorbar_height]
        else:
            ratios = [colorbar_height, colorbar_width]

        if side_colors is not None:
            colors_shape = np.asarray(side_colors).shape
            # This happens when a series or a list is passed
            if len(colors_shape) <= 2:
                n_colors = 1
            # And this happens when a dataframe is passed, the first dimension is number of colors
            else:
                n_colors = colors_shape[0]

            # Multiply side colors size by the number of colors
            expected_size_for_side_colors = n_colors * expected_size_for_side_colors

            side_colors_ratio = expected_size_for_side_colors / figdim

            # Add room for the colors
            ratios += [side_colors_ratio]

        # Add the ratio for the heatmap itself
        ratios.append(1 - sum(ratios))

        return ratios

    # def dim_ratios(self, side_colors, axis, figsize, side_colors_ratio=0.05):
    #     """need to adjust the heatmap height ratio for long figures
    #     such that it fills up more of the given room heatmap.
    #     heatmap_width_ratio is set at .8, default. Tweak to add room for labels
    #     """
    #     ratios = ClusterGrid.dim_ratios(self, side_colors, axis, figsize, side_colors_ratio=side_colors_ratio)

    #     if axis == 0:  # calculating height ratios
    #         ratios[-1] = self.heatmap_height_ratio

    #         if self.dendrogram_width_ratio: #
    #             ratios[0] = self.dendrogram_width_ratio

    #     elif axis == 1 and self.dendrogram_width_ratio:  # calculating width ratios
    #         ratios[0] = self.dendrogram_width_ratio * .4
    #         ratios[1] = self.dendrogram_width_ratio
    #     elif axis ==  1:
    #         ratios[-1] = self.heatmap_width_ratio

    #     # print(axis, ':', ratios)
    #     return ratios
