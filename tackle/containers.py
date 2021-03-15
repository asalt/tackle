import sys
import os
import re
import json
from functools import lru_cache
from datetime import datetime
import operator as op
from collections import OrderedDict, Counter
import itertools
from functools import partial
from warnings import warn
from matplotlib import cm, gridspec

import numpy as np
import pandas as pd
from scipy import stats
from seaborn.matrix import ClusterGrid
from seaborn.matrix import _HeatMapper as HeatMapper
from seaborn.matrix import _matrix_mask

try:
    pd.NA
except AttributeError:
    pd.NA = "NA"

z_score = ClusterGrid.z_score
# python implementation of my R implementation
def my_zscore(x, minval=None, remask=True):
    mask = x.isna()

    if minval is None:
        minval = x.dropna().min()
    x = x.fillna(minval)
    ret = z_score(x)
    if remask:
        ret.loc[mask] = np.nan
    return ret

import logging

def _get_logger():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler("tackle.log")
    # fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

logger = _get_logger()


idx = pd.IndexSlice

TAXON_MAPPER = {
    "human": 9606,
    "mouse": 10090,
    "celegans": 6239,
}


LABEL_MAPPER = {
    "none": 0,  # hard coded number IDs for labels
    "126": 1260,
    "127C": 1270,
    "127N": 1271,
    "128C": 1280,
    "128N": 1281,
    "129C": 1290,
    "129N": 1291,
    "130C": 1300,
    "130N": 1301,
    "131N": 1310,
    "131C": 1311,
    "132N": 1321,
    "132C": 1320,
    "133N": 1331,
    "133C": 1330,
    "134": 1340,
    "1260": 1260,
    "1270": 1270,
    "1271": 1271,
    "1280": 1280,
    "1281": 1281,
    "1290": 1290,
    "1291": 1291,
    "1300": 1300,
    "1301": 1301,
    "1310": 1310,
    "1311": 1311,
    "1320": 1320,
    "1321": 1321,
    "1330": 1330,
    "1331": 1331,
    "131" : 1310,
    1260: "TMT_126",
    1270: "TMT_127_C",
    1271: "TMT_127_N",
    1280: "TMT_128_C",
    1281: "TMT_128_N",
    1290: "TMT_129_C",
    1291: "TMT_129_N",
    1300: "TMT_130_C",
    1301: "TMT_130_N",
    1310: "TMT_131_C",
    1311: "TMT_131_N",
    1321: "TMT_132_C",
    1320: "TMT_132_N",
    1331: "TMT_133_C",
    1330: "TMT_133_N",
    1340: "TMT_134",
    "113": 113,
    "114": 114,
    "115": 115,
    "116": 116,
    "117": 117,
    "118": 118,
    "119": 119,
    "121": 121,
}


def join_and_create_path(*strings, verbose=True):
    """Joins strings and returns resulting path.
    Creates path if does not exist
    """
    path = os.path.join(*strings)
    if not os.path.exists(path):
        os.mkdir(path)
        if verbose:
            print("Created new directory: {}".format(os.path.abspath(path)))
    return path

PWD = os.path.split(os.path.abspath(__file__))[0]

class GeneMapper:
    def __init__(self):
        self.file = os.path.join(
            PWD,
            # 'data', 'genetable_hsmmgg.tab'
            # 'data', 'genetable20200501.tsv'
            "data",
            "genetable20201208.tsv",
        )
        self._df = None
        self._symbol = None
        self._funcat = None
        self._description = None
        self._taxon = None

    @property
    def df(self):
        if self._df is None:
            self._df = pd.read_table(self.file, index_col="GeneID", dtype=str)
            self._df.index = self._df.index.astype(str)
            self._df["FunCats"] = self._df["FunCats"].fillna("")
        return self._df

    @property
    def symbol(self):
        if self._symbol is None:
            self._symbol = self.df["GeneSymbol"].to_dict()
        return self._symbol

    @property
    def funcat(self):
        if self._funcat is None:
            self._funcat = self.df["FunCats"].to_dict()
        return self._funcat

    @property
    def description(self):
        if self._description is None:
            self._description = self.df["GeneDescription"].to_dict()
        return self._description

    @property
    def taxon(self):
        if self._taxon is None:
            self._taxon = self.df["TaxonID"].to_dict()
        return self._taxon

class Annotations:
    def __init__(self):
        self.file = os.path.join(
            PWD, "data", "combined_annotations.tsv"
        )
        if not os.path.exists(self.file):
            logger.error(f"Could not find {self.file}")
            return
        logger.info(f"Loading annotations file {self.file}")
        self.df = pd.read_table(self.file, dtype=str)

        self._categories = [x for x in self.df if x not in ("GeneID", "GeneSymbol")]

    @property
    def categories(self):
        return self._categories

    def get_annot(self, cat):

        df = self.df
        if cat not in df:
            raise ValueError("{cat} must be one of {self.categories}")
        return df[ ~df[cat].isna() ]


_genemapper = None
def get_gene_mapper(cls=GeneMapper):
    global _genemapper
    if _genemapper is None:
        _genemapper = GeneMapper()
    return _genemapper

_annotmapper = None
def get_annotation_mapper(cls=Annotations):
    global _annotmapper
    if _annotmapper is None:
        _annotmapper = Annotations()
    return _annotmapper



def assign_sra(df):

    # df['SRA'] = 'A'
    # cat_type = CategoricalDtype(categories=['S', 'R', 'A'],
    #                             ordered=True)
    # df['SRA'] = df['SRA'].astype(cat_type)
    df["SRA"] = pd.Categorical(
        ["A"] * len(df), categories=["S", "R", "A"], ordered=True
    )
    # df['SRA'] = df['SRA'].astype('category', categories=['S', 'R', 'A'],
    #                              ordered=True)
    df.loc[(df["IDSet"] == 1) & (df["IDGroup_u2g"] <= 3), "SRA"] = "S"

    df.loc[(df["IDSet"] == 2) & (df["IDGroup"] <= 3), "SRA"] = "S"

    df.loc[
        (df["IDSet"] == 1) & (df["SRA"] != "S") & (df["IDGroup_u2g"] <= 5), "SRA"
    ] = "R"

    df.loc[(df["IDSet"] == 2) & (df["SRA"] != "S") & (df["IDGroup"] <= 5), "SRA"] = "R"

    return df


class Data:
    def __init__(
        self,
            annotations=None,
        additional_info=None,
        batch=None,
        batch_nonparametric=False,
        batch_noimputation=False,
        covariate=None,
        cmap_file=None,
        col_cluster=True,
        colors_only=False,
        data_dir="./data",
        base_dir="./results",
        experiment_file=None,
        funcats=None,
        funcats_inverse=None,
        gene_symbols=False,
        geneids=None,
        group=None,
        pairs=None,
        limma=False,
        block=None,
        highlight_geneids=None,
        ignore_geneids=None,
        name=None,
        non_zeros=0,
        nonzero_subgroup=None,
        unique_pepts=0,
        plots=("all",),
        row_cluster=True,
        shade_correlation=True,
        show_metadata=False,
        standard_scale="None",
        stat="pearson",
        taxon="all",
        z_score="0",
        export_all=False,
        SRA="S",
        number_sra=1,
        ifot=False,
        ifot_ki=False,
        ifot_tf=False,
        median=False,
        genefile_norm=None,
        normalize_across_species=False,
        set_outpath=True,
        outpath=None,
        outpath_name=None,
        metrics=False,
        metrics_after_filter=True,
        metrics_unnormed_area=True,
        cluster_annotate_cols=None,
        impute_missing_values=False,
        only_local=False,
    ):
        "docstring"

        if experiment_file is None:
            raise ValueError("Must specify valid experiment_file")

        self.annotations = annotations
        self.additional_info = additional_info
        self.batch = batch
        self.batch_nonparametric = batch_nonparametric
        self.batch_noimputation = batch_noimputation
        self.covariate = covariate
        self.col_cluster = col_cluster
        self.colors_only = colors_only
        self.data_dir = data_dir
        self.experiment_file = experiment_file
        self.funcats = funcats
        self.funcats_inverse = funcats_inverse
        self.gene_symbols = gene_symbols
        self.geneids = geneids
        self.ignore_geneids = ignore_geneids
        self.group = group
        self.pairs = pairs
        self.limma = limma
        self.block = block
        self.highlight_geneids = highlight_geneids
        self.non_zeros = non_zeros
        self.nonzero_subgroup = nonzero_subgroup
        self.unique_pepts = unique_pepts
        self.plots = plots
        self.row_cluster = row_cluster
        self.shade_correlation = shade_correlation
        self.show_metadata = show_metadata
        self.stat = stat
        self.taxon = taxon
        self.standard_scale = self.clean_input(standard_scale)
        self.z_score = self.clean_input(z_score)
        self.export_all = export_all
        self.ifot = ifot
        self.ifot_ki = ifot_ki
        self.ifot_tf = ifot_tf
        self.median = median
        self.genefile_norm = genefile_norm
        self.normalize_across_species = normalize_across_species
        self.base_dir = base_dir
        self.metrics = metrics
        self.metrics_after_filter = metrics_after_filter
        self.metrics_unnormed_area = metrics_unnormed_area
        self._metric_values = None
        self.SRA = SRA
        self.number_sra = number_sra
        self.cluster_annotate_cols = cluster_annotate_cols

        self.outpath = None
        self.analysis_name = None

        self.normed = False
        self.batch_applied = (
            None  # set to the batch (only) upon successful batch correction
        )
        self.gid_funcat_mapping = None  # loaded with data

        self.impute_missing_values = impute_missing_values  # use gaussian distribution 2 sd down from mean of data

        self.set_analysis_name(experiment_file)
        if set_outpath:
            self.set_outpath(base_dir, self.analysis_name, name)
            self.outpath_name = os.path.split(self.outpath)[-1]
        else:
            self.outpath = outpath
            self.outpath_name = outpath_name

        # self.geneid_subset = None
        self.geneid_subset = self.set_geneid_subset(geneids)
        self.ignore_geneid_subset = self.set_geneid_subset(ignore_geneids)

        self.highlight_gids, self.highlight_gid_names = None, None
        self.set_highlight_gids(highlight_geneids)

        self.config = read_config(experiment_file)
        if len(self.config) == 0:
            raise ValueError("No items in configfile.")
        self.labeled_meta = None
        if additional_info:
            self.labeled_meta = read_config(additional_info, enforce=False)

        self.panel = None
        self.load_data(only_local=only_local)

        self._gid_symbol = None

        self.panel_filtered = None
        self.filter_data()

        # self.ibaqs, self.ibaqs_log, self.ibaqs_log_shifted = (None, ) * 3
        self._areas, self._areas_log, self._areas_log_shifted = (None,) * 3

        self._padj = None
        # self.perform_data_export()

        if cmap_file:
            with open(cmap_file) as f:
                jdata = json.load(f)

            self.metadata_colors = jdata
        else:
            self.metadata_colors = None

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
    def padj(self):
        if self._padj is None:
            self.calc_padj()
        return self._padj

    @property
    def metric_values(self):
        return self._metric_values

    @staticmethod
    def clean_input(raw):
        if raw is None:
            return raw
        if raw.isdigit():
            return int(raw)
        elif raw == "None":
            return None
        elif raw is None:
            return None
        else:  # invalid input
            warn(
                """Invalid input for z_score: `{}`.
            Setting to `None`.
            Choose from `None`, `0` or `1`.
            """.format(
                    raw
                )
            )
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
            warn("No geneids found in file {}".format(geneids))
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
                warn("Non geneids found in file {}".format(highlight_geneids))

            highlight_gids.append(highlight_gid)
            h_gid_name = get_file_name(h_gid)
            if h_gid_name:
                highlight_gid_names.append(h_gid_name)
            else:
                highlight_gid_names.append(ix)

        # self.highlight_gids = highlight_gid_names
        self.highlight_gids = highlight_gids
        self.highlight_gid_names = highlight_gid_names

    def _assign_labeled(
        self, record, exp, exps, name, funcats=None, geneid_subset=None
    ):
        labels = dict()
        for key, value in self.config[name].items():
            if key.isdigit():
                df = (
                    exp.df[exp.df.EXPLabelFLAG == int(key)].set_index(["GeneID"])
                    # .pipe(filter_and_assign, value, funcats, geneid_subset)
                    .pipe(assign_cols, value)
                )
                labels[value] = df
                if self.metrics and not self.metrics_after_filter:
                    self._update_metrics(df, name)

        if self.labeled_meta and not all(v in self.labeled_meta for v in labels.keys()):
            warn("Mismatch between `labeled_meta` and input labels")

        else:
            pca = self.config.get("__PCA__")
            for key in self.labeled_meta.keys():
                if key not in labels:
                    continue
                newkey = "{}|{}".format(name, key)
                exps[newkey] = labels[key]

            # exps = OrderedDict((key, exps[key])
            #                    for key in self.labeled_meta.keys())

            # re-assign config to labeled_meta
            # since this is a split of the each experiment ID
            # within each isobaric experiment
            # self.config = self.labeled_meta
            # self.config['__PCA__'] = pca

        return exps

    def _update_metrics(self, df, name, area_column="iBAQ_dstrAdj"):
        if self.metric_values is None:
            self._metric_values = DefaultOrderedDict(lambda: defaultdict(None))

        sra = df.SRA.value_counts().to_dict()
        gpg = df.GPGroup.nunique()
        psms = {
            "Total": df.PSMs.sum(),
            "u2g": df.PSMs_u2g.sum(),
        }
        peptides = {
            "Total": df.PeptideCount.sum(),
            "u2g": df.PeptideCount_u2g.sum(),
            "Strict": df.PeptideCount_S.sum(),
            "Strict_u2g": df.PeptideCount_S_u2g.sum(),
        }
        self._metric_values[name]["SRA"] = sra
        self._metric_values[name]["GPGroups"] = gpg
        self._metric_values[name]["PSMs"] = psms
        self._metric_values[name]["Peptides"] = peptides
        self._metric_values[name]["Area"] = (
            df[area_column].where(lambda x: x > 0).dropna().values
        )

        # digestion efficiency
        allpepts = [
            y for x in df.PeptidePrint.apply(lambda x: x.split("_")).values for y in x
        ]
        allpepts = set(allpepts)
        ## trypsin/P
        miscut_regex = re.compile("[kr].*[kr]", re.I)

        # trypsin
        cutsite_regex = re.compile("([KR](?=[^P]))", re.I)
        miscuts = [len(cutsite_regex.findall(x)) for x in allpepts]
        counter = Counter(miscuts)
        self._metric_values[name]["Trypsin"] = counter

        # trypsin/P
        cutsite_regex = re.compile("[kr](?=.)", re.I)
        miscuts = [len(cutsite_regex.findall(x)) for x in allpepts]
        counter = Counter(miscuts)
        self._metric_values[name]["Trypsin/P"] = counter

        # if self.metrics_unnormed_area:
        #     self._metric_values[name]['Area']     = df.AreaSum_dstrAdj.where(lambda x: x > 0 ).dropna().values
        # elif not self.metrics_unnormed_area:
        #     self._metric_values[name]['Area']     = df.iBAQ_dstrAdj.where(lambda x: x > 0 ).dropna().values
        # self._metric_values[name]['GeneIDs']  = exp.df.GeneID.unique()

    @staticmethod
    @lru_cache()
    def get_e2g(recno, runno, searchno, data_dir, only_local=False):
        e2g = ispec.E2G(
            recno, runno, searchno, data_dir=data_dir, only_local=only_local
        )
        e2g.df["GeneID"] = e2g.df["GeneID"].apply(maybe_int)
        e2g.df.index = e2g.df.GeneID
        standardize_meta(e2g.df)
        return e2g

    def load_data(self, only_local=False):

        exps = OrderedDict()
        gid_funcat_mapping = dict()
        config = self.config

        col_metadata = parse_metadata(self.config)
        taxon_ratios = col_metadata.copy()

        for name, record in config.items():
            # print('Loading', name)
            if name.startswith("__"):
                continue
            labeltype = config[name].get("__LABELTYPE__", "LF")  # depreciated
            recno = record.get("recno")
            runno = record.get("runno")
            searchno = record.get("searchno")
            label = record.get("label")
            #if label and label not in LABEL_MAPPER:

            labelquery = LABEL_MAPPER.get(label, label)

            print("Getting", recno, runno, searchno, label, "to/from", self.data_dir)
            exp = self.get_e2g(
                recno, runno, searchno, data_dir=self.data_dir, only_local=only_local
            )

            if "EXPLabelFLAG" not in exp.df and "LabelFLAG" in exp.df:
                exp.df.rename(columns={"LabelFLAG": "EXPLabelFLAG"}, inplace=True)
            df = exp.df.query("EXPLabelFLAG==@labelquery").copy()
            if df.empty:
                warn(
                    "\n"
                    + "#" * 80
                    + "\n"
                    + "Could not find label {}, in mapping, check if correct label specified. Skipping to avoid error".format(
                        label
                    )
                    + "\n"
                    + "#" * 80
                )
                continue

            # TODO fix this in bcmproteomics
            df.index = df.index.astype("int").astype("str")

            if (
                "SRA" not in df or df.SRA.isna().any()
            ):  # some old grouper experiments don't have it
                df = assign_sra(df)

            if df.empty:
                warn("No data for {!r}, label {}, skipping".format(exp, label))
                continue
            if "9606" not in exp.taxon_ratios:
                # should always be here, or None?
                pass
            # taxon_ratios.loc[name, '9606'] = exp.taxon_ratios['9606']
            # taxon_ratios.loc[name, '10090'] = exp.taxon_ratios['10090']
            # taxon_ratios.loc[name, '9031'] = exp.taxon_ratios['9031']
            taxon_ratios.loc[name, "9606"] = exp.taxon_ratios.get("9606")
            taxon_ratios.loc[name, "10090"] = exp.taxon_ratios.get("10090")
            taxon_ratios.loc[name, "9031"] = exp.taxon_ratios.get("9031")

            # df = exp.df.query('LabelFLAG==@labelquery').copy()

            if not df.index.name == "GeneID":
                df.index = df.GeneID

                # raise ValueError('No data for {!r}'.format(exp))

            if (
                self.metrics
                and not self.metrics_after_filter
                and self.metrics_unnormed_area
            ):
                self._update_metrics(df, name, area_column="iBAQ_dstrAdj")

            # exp.df['GeneID'] = exp.df['GeneID'].astype(int)
            df["GeneID"] = df["GeneID"].apply(maybe_int)

            funcats_dict = (
                df.drop_duplicates("GeneID").set_index("GeneID")["FunCats"].to_dict()
            )
            gid_funcat_mapping.update(funcats_dict)

            if (df.FunCats == "").all():
                df["FunCats"] = df.GeneID.map(_genemapper.funcat).fillna("")
            if "TaxonID" not in df or df.TaxonID.isna().any():
                if "TaxonID" in df:
                    loc = df[df.TaxonID.isna()].index
                    # TODO read taxon in as a "string" to avoid upcast to float with missing values...
                    df.loc[loc, "TaxonID"] = [
                        _genemapper.taxon.get(str(int(x))) for x in loc
                    ]
                else:
                    df.loc[:, "TaxonID"] = [
                        _genemapper.taxon.get(str(int(x))) for x in df.index
                    ]

            if labeltype == "TMT" or labeltype == "iTRAQ":  # depreciated
                exps = self._assign_labeled(
                    record, exp, exps, name, self.funcats, self.geneid_subset
                )
            if self.normalize_across_species:

                df["area"] = normalize(
                    df,
                    name,
                    ifot=self.ifot,
                    ifot_ki=self.ifot_ki,
                    ifot_tf=self.ifot_tf,
                    median=self.median,
                    genefile_norm=self.genefile_norm,
                )

                if self.export_all:
                    df.loc[:, "iBAQ_dstrAdj_FOT"] = normalize(df, ifot=True)
                    df.loc[:, "iBAQ_dstrAdj_MED"] = normalize(df, median=True)
            else:
                # df = filter_and_assign(df, name, self.funcats, self.funcats_inverse,
                #                        self.geneid_subset, self.ignore_geneid_subset, self.ifot,
                #                        self.ifot_ki, self.ifot_tf, self.median)
                for taxonid in df.TaxonID.unique():
                    if taxonid == 0.0:
                        continue  # invalid
                    df.loc[df.TaxonID == taxonid, "area"] = normalize(
                        df.loc[df.TaxonID == taxonid],
                        name,
                        ifot=self.ifot,
                        ifot_ki=self.ifot_ki,
                        ifot_tf=self.ifot_tf,
                        median=self.median,
                        genefile_norm=self.genefile_norm,
                        taxon=taxonid,
                    )
                # df = normalize(df, name, ifot=self.ifot, ifot_ki=self.ifot_ki, ifot_tf=self.ifot_tf,
                #                median=self.median)
                if self.export_all:  # have to calculate more columns

                    for taxonid in df.TaxonID.unique():

                        df.loc[df.TaxonID == taxonid, "iBAQ_dstrAdj_FOT"] = normalize(
                            df.loc[df.TaxonID == taxonid], ifot=True
                        )

                        df.loc[df.TaxonID == taxonid, "iBAQ_dstrAdj_MED"] = normalize(
                            df.loc[df.TaxonID == taxonid], median=True
                        )

                    # df = (df.pipe(normalize, ifot=True, outcol='iBAQ_dstrAdj_FOT')
                    #       .pipe(normalize, median=True, outcol='iBAQ_dstrAdj_MED')
                    # )

            dummy_filter = lambda x, *args, **kwargs: x
            taxon_filter = TAXON_MAPPER.get(self.taxon)
            if taxon_filter is None:
                filter_func = dummy_filter
            else:
                filter_func = lambda x: x[x["TaxonID"] == taxon_filter]
            if (
                self.metrics
                and not self.metrics_after_filter
                and not self.metrics_unnormed_area
            ):
                self._update_metrics(df, name, area_column="area")

            df = genefilter(
                df,
                funcats=self.funcats,
                funcats_inverse=self.funcats_inverse,
                geneid_subset=self.geneid_subset,
                ignored_geneid_subset=self.ignore_geneid_subset,
            ).pipe(filter_func)

            # unique peptide filter
            df = df[df.PeptideCount_u2g >= self.unique_pepts]

            # df = assign_cols(exp.df, name)
            if (
                self.metrics
                and self.metrics_after_filter
                and not self.metrics_unnormed_area
            ):
                self._update_metrics(df, name, area_column="area")
            elif (
                self.metrics
                and self.metrics_after_filter
                and self.metrics_unnormed_area
            ):
                raise NotImplementedError(
                    "Not available. Use --after-norm in conjunction with --after-filter"
                )
            # exps[name] = df.set_index(df.index.astype(int))
            df.index = [maybe_int(x) for x in df.index]
            exps[name] = df

        self.taxon_ratios = taxon_ratios

        self.gid_funcat_mapping = gid_funcat_mapping

        # filter down by exps that actually have data and were loaded
        self.col_metadata = col_metadata.loc[exps.keys()]

        # self.multi = pd.concat(exps.values(), keys=exps.keys())
        self.exps = exps
        ## TODO can check to ensure not exporting all data and stack this smaller amount of data

        if not self.export_all:
            _cols = [
                "TaxonID",
                "IDSet",
                "GeneSymbol",
                "iBAQ_dstrAdj",
                "FunCats",
                "SRA",
                "area",
                "PeptideCount_u2g",
            ]
            if self.cluster_annotate_cols:
                for x in self.cluster_annotate_cols:
                    if x not in _cols:
                        _cols.append(x)
            stacked_data = [df[_cols].stack() for df in exps.values()]
        else:
            stacked_data = [df.stack() for df in exps.values()]
        print("stacking...", flush=True, end="")
        self.data = pd.concat(stacked_data, axis=1, keys=exps.keys())
        self.data.index.names = ["GeneID", "Metric"]
        self.data = self.data.reset_index()
        print("done", flush=True)
        # self.panel = pd.Panel(exps)
        for ax in (
            "GeneCapacity",
            "GeneSymbol",
            "GeneDescription",
            "FunCats",
            "TaxonID",
        ):
            fillna_meta(self.data, ax)

        # if self.additional_info: # depreciated, remove
        #     labeled_meta = read_config(self.additional_info, enforce=False)
        #     additional_metadata = parse_metadata(labeled_meta)
        #     metadata_dict = dict()
        #     for col in self.data.columns:
        #         exp, label = col.split('|')
        #         try:
        #             metadata = additional_metadata[label]
        #         except KeyError:
        #             warn('Mismatch between `labeled_meta` and input labels')
        #             continue
        #         metadata_dict[col] = metadata
        #     col_metadata = pd.DataFrame.from_dict(metadata_dict)

        # else:
        #     col_metadata = parse_metadata(self.config)

        # self.col_metadata = col_metadata

    @property
    def gid_symbol(self):
        if self._gid_symbol is None:
            try:
                # sel = self.data.loc[ idx[:, 'GeneSymbol'], self.data.columns[0] ]
                sel = self.data.loc[
                    self.data.Metric == "GeneSymbol", self.data.columns[0:3]
                ]
                # sel.index = sel.index.droplevel(1)
                self._gid_symbol = sel.set_index("GeneID")[sel.columns[-1]].to_dict()
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
        df_filtered = (
            self.data.pipe(
                filter_observations,
                "area",
                self.non_zeros,
                self.nonzero_subgroup,
                self.col_metadata,
            ).pipe(filter_sra, SRA=self.SRA, number_sra=self.number_sra)
            # .pipe(filter_upept, number=1)
            .pipe(filter_func)
        )
        self.df_filtered = df_filtered.set_index(["GeneID", "Metric"])

    # def impute_missing(self, frame):
    #     _norm_notna = frame.replace(0, np.NAN).stack().apply(np.log10)
    #     _norm_notna += np.abs(_norm_notna.min())
    #     _mean = _norm_notna.mean()
    #     _sd = _norm_notna.std()
    #     _norm = stats.norm(loc=_mean-(_sd*2), scale=_sd)
    #     _number_na = self._areas.replace(0, np.NAN).isna().sum().sum()
    #     # print(frame.replace(0, np.NAN).isna().sum())
    #     random_values = _norm.rvs(size=_number_na, random_state=1234)

    #     _areas_log = np.log10(frame.replace(0, np.NAN))
    #     _areas_log += np.abs(_areas_log.min().min())

    #     start_ix = 0
    #     for col in _areas_log:
    #         last_ix = _areas_log[col].isna().sum()
    #         # print(_areas_log[col].isna().sum())
    #         _areas_log.loc[_areas_log[col].isna(), col] = random_values[start_ix: start_ix+last_ix]
    #         start_ix += last_ix
    #     return _areas_log

    def set_area_dfs(self):
        # self.areas = self.panel_filtered.minor_xs('iBAQ_dstrAdj').astype(float)
        self._areas = self.df_filtered.loc[idx[:, "area"], :]
        self._areas.index = self._areas.index.droplevel(
            1
        )  # don't need the second index

        batch_info = self.config.get("__batch__")  # depreciated
        if batch_info:
            batch = batch_info.get("batch")
            # metadata = self.col_metadata.T  # rows are experiments, cols are metadata
            metadata = self.col_metadata  # rows are experiments, cols are metadata
            areas = self._areas.copy()
            for ix, g in metadata.groupby(batch):
                meanvals = areas[g.index].mean(1)
                areas.loc[:, g.index] = self._areas[g.index].div(
                    meanvals.fillna(0) + 1e-20, axis="index"
                )

            self._areas = areas

        self._mask = self._areas.applymap(np.isnan)
        self._zeros = self._areas == 0
        if len(self.areas) == 0:
            raise ValueError("No data")

        gids = set(self._areas.index)
        if self.funcats:
            gids &= set(
                self.df_filtered.pipe(
                    filter_funcats, self.funcats
                ).index.get_level_values(0)
            )
        if self.geneid_subset:
            gids &= set(self.geneid_subset)

        gids = tuple(gids)
        self.minval = self._areas.replace(0, np.NAN).stack().dropna().min()

        # if self.impute_missing_values or 1:
        if self.impute_missing_values :
            import ipdb; ipdb.set_trace()
            # downshift=2.
            # scale = .5
            downshift, scale = 1.8, 0.8

            inpute_plotname = get_outname(
                "distribution_ds_{:.2g}_scale_{:.2g}".format(downshift, scale),
                name=self.outpath_name,
                taxon=self.taxon,
                non_zeros=self.non_zeros,
                colors_only=self.colors_only,
                batch=self.batch_applied,
                batch_method="parametric"
                if not self.batch_nonparametric
                else "nonparametric",
                outpath=self.outpath,
            )

            to_impute = (
                self._areas.replace(0, np.NAN)
                .divide(self.minval)
                .add(1)
                .applymap(np.log10)
            )
            imputed = utils.impute_missing(
                to_impute, downshift=downshift, scale=scale, make_plot=True
            )
            plt.savefig(inpute_plotname + ".png", dpi=90)
            plt.close(plt.gcf())
            self._areas_log = imputed

        # self._areas_log = np.log10(self._areas.fillna(0)+1e-10)
        else:
            # self._areas_log = np.log10(
            #     self._areas.replace(0, np.NAN)
            #     .fillna(self.minval / 2)
            #     .divide(self.minval / 2)
            # )
            self._areas_log = np.log10(
                self._areas.replace(0, np.NAN)
                .fillna(self.minval )
                .divide(self.minval )
            )

        self._areas_log.index.name = "GeneID"
        # fillna with the mean value. This prevents skewing of normalization such as
        # z score. The NAN values are held in the self.mask dataframe
        # self._areas_log = np.log10(self._areas.T.fillna(self._areas.mean(axis=1)).T + 1e-8)

        # if norm_info is not None:  # do not shift the values
        #     self._areas_log_shifted = self._areas_log
        #     return

        # don't need this section anymore since now minval is 0
        minval = self._areas_log.min().min()
        shift_val = np.ceil(np.abs(minval))
        self.minval_log = minval

        self._areas_log_shifted = self._areas_log + shift_val
        self._areas_log_shifted.index.name = "GeneID"
        # if self.geneid_subset: # don't do this here, already done
        #     self._areas_log_shifted = self._areas_log_shifted.loc[self.geneid_subset]

        # if specified, normalize by a specified control group
        norm_info = self.config.get("__norm__")
        if norm_info is not None:
            control = norm_info["control"]
            group = norm_info["group"]
            label = norm_info["label"]
            # metadata = self.col_metadata.T  # rows are experiments, cols are metadat
            metadata = self.col_metadata  # rows are experiments, cols are metadata
            areas = self._areas_log_shifted.copy()
            ctrl_exps = list()

            for ix, g in metadata.groupby(group):
                ctrl_exp = g[g[label] == control].index[0]  # should only be one
                ctrl_exps.append(ctrl_exp)
                # to_normalize = list(set(g.index) - set([ctrl_exp]))
                to_normalize = g.index
                areas.loc[:, to_normalize] = self._areas_log_shifted[to_normalize].sub(
                    self._areas_log_shifted[ctrl_exp].fillna(0) + 1e-20, axis="index"
                )
            self._areas = (
                areas.drop(ctrl_exps, axis=1).where(lambda x: x != 0).dropna(how="all")
            )
            self.minval_log = self._areas.replace(0, np.NAN).stack().dropna().min()
            finite = self._areas.pipe(np.isfinite)

            # have to take care of number/0 case
            maxval_log = self._areas[finite].stack().dropna().max()
            self._areas_log_shifted = self._areas.fillna(self.minval_log / 2).replace(
                np.inf, maxval_log * 1.5
            )
            # not sure if this is right, can check:
            # sample_cols = [x for x in self.col_metadata.columns if x not in ctrl_exps]
            sample_cols = self.col_metadata.columns
            sample_ixs = [x for x in self.col_metadata.index if x not in ctrl_exps]

            self.col_metadata = self.col_metadata.loc[sample_ixs, sample_cols]
            self._mask = self.mask[sample_ixs]
            self.normed = True

        if self.export_all and not self.normed:
            new_cols = list()
            for col in "iBAQ_dstrAdj", "iBAQ_dstrAdj_FOT", "iBAQ_dstrAdj_MED":
                frame = self.df_filtered.loc[idx[:, col], :]
                if self.impute_missing_values:
                    frame_log = self.impute_missing(frame).reset_index()
                else:
                    minval = frame.replace(0, np.NAN).stack().dropna().min()
                    frame_log = np.log10(frame.replace(0, np.NAN).fillna(minval / 2))
                    minval_log = frame_log.replace(0, np.NAN).stack().dropna().min()
                    shift_val = np.ceil(np.abs(minval_log))
                    frame_log = (
                        frame_log.fillna(minval_log / 2) + shift_val
                    ).reset_index()
                frame_log["Metric"] = col + "_log10"
                new_cols.append(frame_log.set_index(["GeneID", "Metric"]))

                # self.df_filtered.loc[ idx[:, col+'_log'], :] = np.NAN
                # self.data.loc[ idx[:, col+'_log'], :] = frame_log.values
                # self.data.loc[ idx[frame_log.index.levels[0].values, col+'_log'], :] = frame_log.values
            self.df_filtered = pd.concat([self.df_filtered, *new_cols])

        if self.batch is not None:
            self._areas_log_shifted = self.batch_normalize(self.areas_log_shifted)
            # try batch normalization via ComBat
            if (
                self.export_all and not self.normed
            ):  # batch normalize the other requested area columns
                for col in "iBAQ_dstrAdj", "iBAQ_dstrAdj_FOT", "iBAQ_dstrAdj_MED":
                    frame = (
                        self.df_filtered.loc[idx[:, col + "_log10"], :]
                        .reset_index(level=1, drop=True)
                        .astype(float)
                    )
                    # get rid of "Metric" index level, batch_normalize not expecting it
                    res = self.batch_normalize(frame, prior_plot=False)
                    self.df_filtered.loc[idx[:, col + "_log10"], :] = res.values
                    self.df_filtered.loc[idx[:, col], :] = np.power(res, 10).values
            elif self.export_all and self.normed:
                warn(
                    """No support for full export with norm channel. Please ignore
                        iBAQ_dstrAdj, iBAQ_dstrAdj_FOT, and iBAQ_dstrAdj_MED
                        as they have not been properly batch corrected and
                        normalized
                """
                )

        if (
            self.export_all and not self.normed
        ):  # now calculate z scores for all these extra columns
            new_cols = list()
            for col in (
                "iBAQ_dstrAdj_log10",
                "iBAQ_dstrAdj_FOT_log10",
                "iBAQ_dstrAdj_MED_log10",
            ):
                frame = (
                    self.df_filtered.loc[idx[:, col], :]
                    .apply(z_score, axis=1)
                    .reset_index()
                )
                frame["Metric"] = col + "_zscore"
                new_cols.append(frame.set_index(["GeneID", "Metric"]))
            self.df_filtered = pd.concat([self.df_filtered, *new_cols])

            # new_cols.append(frame_log.set_index(['GeneID', 'Metric']))
            # self.data

        # if self.group is not None:
        #     self.calc_qvals()

    def batch_normalize(self, data, prior_plot=True):
        """"""

        try:
            from rpy2.rinterface import RRuntimeError
        except Exception as e:
            RRuntimeError = type("RRuntimeError", (Exception,), {})
        try:
            from rpy2.robjects import r
            from rpy2.robjects.packages import importr

            sva = importr("sva")
        except ModuleNotFoundError as e:
            print("Failure in rpy2 import and load", e, sep="\n")
            return
        except RRuntimeError as e:
            print("sva is not installed", e, sep="\n")
            return

        from rpy2.robjects import pandas2ri

        pandas2ri.activate()
        grdevices = importr("grDevices")

        # pheno = self.col_metadata.T
        pheno = self.col_metadata.copy()
        for col in ("recno", "runno", "searchno"):
            pheno[col] = pheno[col].astype(str)
        try:
            pd.NA
        except AttributeError:
            pd.NA = "NA"  # bugfix for pandas < 1.0
        r.assign("pheno", pheno)

        if self.covariate is not None:
            r("mod  <- model.matrix(~as.factor({}), pheno)".format(self.covariate))
            mod = r["mod"]
            self.batch_applied = self.batch + "_Cov_{}".format(self.group)
        else:
            mod = r("model.matrix(~1, pheno)")
            self.batch_applied = self.batch + "_noCov"
        # r.assign('batch', 'pheno${}'.format(self.batch))
        batch = pheno[self.batch]
        # res = sva.ComBat(dat=self._areas_log_shifted.fillna(0), batch=batch,
        #                  mod=mod, par_prior=True, mean_only=False)
        if not self.batch_nonparametric and prior_plot:
            plot_prior = True
            outname = get_outname(
                "Combat_prior_plots",
                name=self.outpath_name,
                taxon=self.taxon,
                non_zeros=self.non_zeros,
                colors_only=self.colors_only,
                batch=self.batch_applied,
                batch_method="parametric"
                if not self.batch_nonparametric
                else "nonparametric",
                outpath=self.outpath,
            )
            grdevices.png(file=outname + ".png", width=5, height=5, units="in", res=300)
        else:
            plot_prior = False
        res = sva.ComBat(
            dat=data.values,
            batch=batch,
            mod=mod,
            par_prior=not self.batch_nonparametric,
            mean_only=False,
            prior_plots=plot_prior,
        )

        if plot_prior:
            grdevices.dev_off()

        try:
            df = pd.DataFrame(
                index=data.index,
                columns=data.columns,
                # data=pandas2ri.ri2py(res)
                data=res,
            )
        except ValueError:  # rpy2.9.5 support
            df = pd.DataFrame(
                index=data.index, columns=data.columns, data=pandas2ri.ri2py(res)
            )

        # df = pandas2ri.ri2py(res)
        nas = sum(df.isnull().any(1))
        if nas > 0:
            print(
                "{} Gene Product(s) became NAN after batch normalization, dropping".format(
                    nas
                )
            )

        # df.index = df.index.astype(int)
        df.index = data.index
        # df.index = [maybe_int(x) for x in df.index]
        df.columns = data.columns

        # ===============================================================================
        # TODO:why do this?
        # reassign mask - ComBat can impute some NA values
        # : resolve this for normed data
        # if not self.normed: # ??

        # if not self.batch_noimputation:  # else leave old mask
        #     thresh = (
        #         self.areas_log_shifted[(~self.mask) & (self._areas_log_shifted > 0)]
        #         .min()
        #         .min()
        #     )
        #     new_mask = df[self.mask] <= thresh
        #     new_mask.columns = self._areas_log_shifted.columns
        #     self._mask = new_mask

        # fill in combat funny values back to zero
        # df[self.mask] = 0

        # self._areas_log_shifted = df.dropna(how='any')

        df = df.dropna(how="any")
        # self._areas_log_shifted.columns = self._areas_log.columns  # r changes certain characters in column names
        # self._areas_log_shifted.columns = self._areas.columns  # r changes certain characters in column names
        df.columns = data.columns  # r changes certain characters in column names
        # self._areas_log_shifted.index = self._areas_log_shifted.index.astype(int)  #r converts rownames to str

        df.index = [maybe_int(x) for x in df.index]

        df.index.name = "GeneID"
        return df

    # def calc_padj(self):
    def stat_model(self, formula=None, contrasts_str=None, impute_missing_values=False):
        """
        Still work in progress.

        Should return a dictionary of form:
        {comparison: DataFrame of results}
        """
        try:
            ModuleNotFoundError
        except:
            ModuleNotFoundError = type("ModuleNotFoundError", (Exception,), {})
        try:
            from rpy2.rinterface import RRuntimeError
        except Exception as e:
            RRuntimeError = type("RRuntimeError", (Exception,), {})
        try:
            from rpy2 import robjects
            from rpy2.robjects import r
            from rpy2.robjects.packages import importr

            sva = importr("sva")
        except ModuleNotFoundError as e:
            print("Failure in rpy2 import and load", e, sep="\n")
            return
        except RRuntimeError as e:
            print("sva is not installed", e, sep="\n")
            return
        from rpy2.robjects import pandas2ri

        pandas2ri.activate()
        r_source = r["source"]
        r_file = os.path.join(
            os.path.split(os.path.abspath(__file__))[0], "R", "pvalue_cov.R"
        )
        r_source(r_file)

        mat = self.areas_log_shifted
        if impute_missing_values:

            downshift, scale = 1.8, 0.8

            inpute_plotname = get_outname(
                "distribution_ds_{:.2g}_scale_{:.2g}".format(downshift, scale),
                name=self.outpath_name,
                taxon=self.taxon,
                non_zeros=self.non_zeros,
                colors_only=self.colors_only,
                batch=self.batch_applied,
                batch_method="parametric"
                if not self.batch_nonparametric
                else "nonparametric",
                outpath=self.outpath,
            )

            mask_values = mat[self.mask].copy()

            mat[self.mask] = np.nan
            # mat[mat == 0] = np.nan
            mat = utils.impute_missing(
                mat, downshift=downshift, scale=scale, make_plot=True
            )
            # now add back the mask values?
            #- this might be important when batch correction is performed, otherwise should have no effect?
            # mat = mat + mask_values.fillna(0)
            plt.savefig(inpute_plotname + ".png", dpi=90)
            plt.close(plt.gcf())

        r.assign("edata", mat)

        # pheno = self.col_metadata.T
        pheno = self.col_metadata.copy()
        for col in ("recno", "runno", "searchno"):
            pheno[col] = pheno[col].astype(str)
        robjects.r.assign("pheno", pheno)
        r("mod0 <- model.matrix(~1, pheno)")

        if self.group and not formula:
            mod = r("mod  <- model.matrix(~0+{}, pheno)".format(self.group))
        elif formula:
            mod = r("mod <- model.matrix({}, pheno)".format(formula))
        else:
            raise ValueError("Must specify 1 of `group` or `formula`")

        if self.covariate is not None:
            ncov = pheno[self.covariate].nunique()
            r.assign("ncov", pheno[self.covariate].nunique())
        else:
            r.assign("ncov", 0)
            ncov = 0

        if not self.pairs and not self.limma:  # standard t test
            pvalues = r("pvalue.batch(as.matrix(edata), mod, mod0, ncov)")
        elif not self.pairs and self.limma:
            importr("limma")
            # r('suppressMessages(library(dplyr))')
            importr("dplyr", on_conflict="warn")
            r("block <- NULL")
            r("cor <- NULL")
            if self.block:

                r('block <- as.factor(pheno[["{}"]])'.format(self.block))
                r("corfit <- duplicateCorrelation(edata, design = mod,  block = block)")
                r("cor <- corfit$consensus")

            fit = r("""fit <- lmFit(as.matrix(edata), mod, block = block, cor = cor)""")
            # need to make valid R colnames
            variables = robjects.r("colnames(mod)")
            fixed_vars = [
                x.replace(
                    ":",
                    "_",
                )
                .replace(" ", "_")
                .replace("-", "_")
                .replace("+", "_")
                for x in variables
            ]
            robjects.r.assign("fixed_vars", fixed_vars)
            robjects.r("colnames(mod) <- fixed_vars")

            if contrasts_str is None:
                contrasts_array = list()
                for group in itertools.combinations(
                    (x for x in fixed_vars if "Intercept" not in x), 2
                ):
                    contrast = "{} - {}".format(*group)
                    contrasts_array.append(contrast)

            elif contrasts_str:
                contrasts_array = [
                    x.strip() for x in contrasts_str.split(",") if x.strip()
                ]

            robjects.r("print(mod)")
            # robjects.r('print(fit)')

            robjects.r.assign("contrasts_array", contrasts_array)

            robjects.r(
                """contrasts_matrix <- makeContrasts(contrasts = contrasts_array,
                                                            levels = mod
            )"""
            )

            contrast_fit = robjects.r(
                """
            fit2 <- contrasts.fit(fit, contrasts_matrix) %>% eBayes(robust=TRUE, trend=TRUE)
            """
            )

            results = dict()
            for coef, name in enumerate(contrasts_array, 1):

                result = r(
                    """ topTable(fit2, n=Inf, sort.by='none', coef={}, confint=TRUE) """.format(
                        coef
                    )
                ).rename(
                    columns={
                        "logFC": "log2_FC",
                        "adj.P.Val": "pAdj",
                        "P.Value": "pValue",
                    }
                )
                result["log2_FC"] = result["log2_FC"].apply(lambda x: x / np.log10(2))
                result["CI.L"] = result["CI.L"].apply(lambda x: x / np.log10(2))
                result["CI.R"] = result["CI.R"].apply(lambda x: x / np.log10(2))

                # DON'T NEED TO DO THIS ANYMORE
                # we ensure the order of result is equal to order of areas_log_shifted
                # to preserve GeneID order
                # result.index = self.areas_log_shifted.index

                results[name] = result

            return results

            results = r(
                """lmFit(as.matrix(edata), mod, block = block, cor = cor) %>%
                       eBayes(robust=TRUE, trend=TRUE) %>%
                       topTable(n=Inf, sort.by='none')
            """.format(
                    self.group
                )
            )
            pvalues = results["P.Value"]
            padj = results["adj.P.Val"]
        elif self.pairs and self.limma:
            importr("limma")
            r("library(dplyr)")
            r("mod  <- model.matrix(~{}+{}, pheno)".format(self.group, self.pairs))
            results = r(
                """lmFit(as.matrix(edata), mod) %>%
                       eBayes(robust=TRUE, trend=TRUE) %>%
                       topTable(n=Inf, sort.by='none', coef="{}")
            """.format(
                    self.group + pheno[self.group].iloc[-1]
                )
            )
            pvalues = results["P.Value"]
            padj = results["adj.P.Val"]

        else:  # ttest rel
            from .ttest_rel_cov import ttest_rel_cov

            # groups = self.col_metadata.loc[self.group].unique()
            groups = self.col_metadata[self.group].unique()
            # only have t-test implemented here
            if len(groups) > 2:
                raise NotImplementedError("Only have support for 2 groups")
            group0, group1 = groups
            # meta_ordered = self.col_metadata.sort_index(by=[self.group, self.pairs], axis=1)
            meta_ordered = self.col_metadata.sort_values(
                by=[self.group, self.pairs], axis=1
            )
            # better way to do this?
            # cols0 = (meta_ordered.T[self.group] == group0).apply(lambda x: x if x else np.nan).dropna().index
            # cols1 = (meta_ordered.T[self.group] == group1).apply(lambda x: x if x else np.nan).dropna().index
            cols0 = (
                (meta_ordered[self.group] == group0)
                .apply(lambda x: x if x else np.nan)
                .dropna()
                .index
            )
            cols1 = (
                (meta_ordered[self.group] == group1)
                .apply(lambda x: x if x else np.nan)
                .dropna()
                .index
            )
            # self.
            # t test on every gene set
            def ttest_func(row):
                t, p = ttest_rel_cov(row[cols0], row[cols1], ncov=ncov)
                return p

            pvalues = self.areas_log_shifted.apply(ttest_func, axis=1)

        p_adjust = r["p.adjust"]

        # pvalues = f_pvalue(self.areas_log_shifted.fillna(0), mod, mod0)
        try:
            padj  # limma calculates this
        except NameError:
            padj = pandas2ri.ri2py(p_adjust(pvalues, method="BH"))

        padj = pd.DataFrame(
            index=self.areas_log_shifted.index,
            data=np.array([pvalues, padj]).T,
            columns=["pValue", "pAdj"],
        ).sort_values(by="pValue")
        # qvalues.name = 'q_value'

        self._padj = padj

    def make_plot(self, pltname):
        if "all" in self.plots:
            return True
        if pltname in self.plots:
            return True
        return False

    def perform_data_export(self, level="all", genesymbols=False, linear=False):

        # fname = '{}_data_{}_{}_more_zeros.tab'.format(level,
        #                                               self.outpath_name,
        #                                               self.non_zeros)

        # outname = os.path.abspath(os.path.join(self.outpath, fname))
        # if self.export_data == 'all':
        self.areas_log_shifted  # make sure it's created

        level_formatter = level
        if level == "area" and linear:
            level_formatter = level + "_linear"

        outname = (
            get_outname(
                "data_{}".format(level_formatter),
                name=self.outpath_name,
                taxon=self.taxon,
                non_zeros=self.non_zeros,
                colors_only=self.colors_only,
                batch=self.batch_applied,
                batch_method="parametric"
                if not self.batch_nonparametric
                else "nonparametric",
                outpath=self.outpath,
            )
            + ".tsv"
        )

        if level == "all":
            self.df_filtered.sort_index(level=[0, 1]).to_csv(outname, sep="\t")

        elif level == "MSPC":
            # TODO: export ALL data, not just filtered
            export = self.df_filtered.sort_index(level=[0, 1])
            data = list()
            # cols = export.index.get_level_values(1).unique()
            gene_metadata_cols = [
                "GeneID",
                "TaxonID",
                "GeneSymbol",
                "GeneDescription",
                "FunCats",
                "GeneCapacity",
            ]
            # TODO fix this
            for c in gene_metadata_cols:
                try:
                    export.loc[idx[:, c], :] = (
                        export.loc[idx[:, c], :]
                        .fillna(
                            method="ffill",
                            axis=1,
                        )
                        .fillna(
                            method="bfill",
                            axis=1,
                        )
                    )
                except KeyError:  # BAD
                    print("KeyError!")
                    # pass

            gene_metadata_cols = [
                "GeneID",
                "TaxonID",
                "GeneSymbol",
                "GeneDescription",
                "FunCats",
                "GeneCapacity",
            ]

            cols = [
                "SRA",
                "PeptideCount",
                "PeptideCount_u2g",
                "PSMs",
                "AreaSum_dstrAdj",
                "iBAQ_dstrAdj",
                "iBAQ_dstrAdj_log10",
                "iBAQ_dstrAdj_MED",
                "iBAQ_dstrAdj_MED_log10",
                "iBAQ_dstrAdj_MED_log10_zscore",
            ]

            for col in export.columns:
                renamer = {x: "{}_{}".format(x, col) for x in cols}

                subdf = (
                    export.loc[idx[:, :], col]
                    .reset_index()
                    .pivot(index="GeneID", columns="Metric")
                )
                if "GeneID" in subdf:
                    subdf = subdf.drop("GeneID", 1)
                subdf.columns = subdf.columns.droplevel(0)
                subdf.index = subdf.index.map(maybe_int).astype(str)
                subdf["GeneID"] = subdf["GeneID"].apply(maybe_int).astype(str)
                # this makes sure we don't crash if any columns are missing
                _cols = [x for x in cols if x in subdf]
                subdf = subdf.set_index([x for x in gene_metadata_cols if x in subdf])[
                    _cols
                ].rename(columns=renamer)
                data.append(subdf)

            for_export = pd.concat(data, axis=1).reset_index()
            for_export["GeneID"] = for_export["GeneID"].apply(maybe_int)
            for_export["TaxonID"] = for_export["TaxonID"].apply(maybe_int)
            print("Writing {}".format(outname))
            for_export.to_csv(outname, sep="\t", index=False)

        elif level == "align":

            export = self.df_filtered.sort_index(level=[0, 1])
            column_number_mapping = dict()
            data = list()
            # cols = export.index.get_level_values(1).unique()
            gene_metadata_cols = [
                "GeneID",
                "TaxonID",
                "GeneSymbol",
                "GeneDescription",
                "FunCats",
                "GeneCapacity",
            ]
            for c in gene_metadata_cols:
                try:
                    export.loc[idx[:, c], :] = (
                        export.loc[idx[:, c], :]
                        .fillna(
                            method="ffill",
                            axis=1,
                        )
                        .fillna(
                            method="bfill",
                            axis=1,
                        )
                    )
                except KeyError:
                    pass

            cols = [
                "SRA",
                "IDSet",
                "IDGroup",
                "IDGroup_u2g",
                "GPGroup",
                "GPGroups_All",
                "ProteinGI_GIDGroups",
                "ProteinGI_GIDGroupCount",
                "PeptidePrint",
                "Coverage",
                "Coverage_u2g",
                "PeptideCount",
                "PeptideCount_u2g",
                "PeptideCount_S",
                "PeptideCount_S_u2g",
                "PSMs",
                "PSMs_u2g",
                "AreaSum_u2g_0",
                "AreaSum_u2g_all",
                "AreaSum_max",
                "AreaSum_dstrAdj",
                "iBAQ_dstrAdj",
                "iBAQ_dstrAdj_log10",
                "iBAQ_dstrAdj_log10_zscore",
                "iBAQ_dstrAdj_FOT",
                "iBAQ_dstrAdj_FOT_log10",
                "iBAQ_dstrAdj_FOT_log10_zscore",
                "iBAQ_dstrAdj_MED",
                "iBAQ_dstrAdj_MED_log10",
                "iBAQ_dstrAdj_MED_log10_zscore",
            ]

            if "GeneCapacity" in export.index.get_level_values(1):
                if not all(
                    export.loc[idx[:, "GeneCapacity"], :].T.nunique() == 1
                ):  # different genecapacities...
                    gene_metadata_cols = [
                        x for x in gene_metadata_cols if x != "GeneCapacity"
                    ]
                    cols.append("GeneCapacity")

            gene_metadata = dict()

            for ix, col in enumerate(export.columns, 1):

                identifier = self.col_metadata.loc[col][
                    ["recno", "runno", "searchno", "label"]
                ].to_dict()

                renamer = {
                    x: "{}_{}_{}_{recno}_{runno}_{searchno}_{label}".format(
                        x, ix, col, **identifier
                    )
                    for x in cols
                }

                subdf = (
                    export.loc[idx[:, :], col]
                    .reset_index()
                    .pivot(index="GeneID", columns="Metric")
                )
                if "GeneID" in subdf:
                    subdf = subdf.drop("GeneID", 1)
                subdf.columns = subdf.columns.droplevel(0)
                subdf.index = subdf.index.map(maybe_int).astype(str)
                subdf["GeneID"] = subdf["GeneID"].apply(maybe_int).astype(str)
                # subdf['GeneID'] = subdf.index #hack
                # subdf['Description'] =
                # this makes sure we don't crash if any columns are missing
                _cols = [x for x in cols if x in subdf]

                subdf = subdf.set_index([x for x in gene_metadata_cols if x in subdf])[
                    _cols
                ].rename(columns=renamer)

                metadata = dict(self.config[col])
                metadata["name"] = col
                to_pop = [
                    x
                    for x in metadata
                    if x not in ("recno", "runno", "searchno", "label", "name")
                ]
                for p in to_pop:
                    metadata.pop(p)
                column_number_mapping[ix] = metadata
                data.append(subdf)
            for_export = pd.concat(data, axis=1).reset_index()
            for_export["GeneID"] = for_export["GeneID"].apply(maybe_int)
            for_export["TaxonID"] = for_export["TaxonID"].apply(maybe_int)
            for_export.to_csv(outname, sep="\t", index=False)
            meta_df = pd.DataFrame(column_number_mapping).T
            outname = (
                get_outname(
                    "metadata_{}".format(level),
                    name=self.outpath_name,
                    taxon=self.taxon,
                    non_zeros=self.non_zeros,
                    colors_only=self.colors_only,
                    batch=self.batch_applied,
                    batch_method="parametric"
                    if not self.batch_nonparametric
                    else "nonparametric",
                    outpath=self.outpath,
                )
                + ".tsv"
            )
            meta_df.to_csv(outname, sep="\t")

        elif level == "area":
            export = self.areas_log_shifted.copy()
            if linear:
                export = export.apply(lambda x: 10 ** x)
            if not self.impute_missing_values:
                export[self.areas == 0] = 0  # fill the zeros back
                export[self.mask] = np.NaN
            order = export.columns
            if genesymbols:
                # export['GeneSymbol'] = export.index.map(lambda x: self.gid_symbol.get(x, '?'))
                export["GeneSymbol"] = export.index.map(
                    lambda x: self.gid_symbol.get(
                        x,
                        # _genemapper.symbol.get(x, '?')
                        _genemapper.symbol.get(str(int(x)), x),
                    )
                )
                order = ["GeneSymbol"]  # index column is GeneID, add GeneSymbol
                order += [x for x in export.columns if x not in order]
            export[order].to_csv(outname, sep="\t")

        # elif level == 'SRA':
        #     export = self.data.loc[ self.data.Metric=='SRA' ]
        elif level == 'zscore':  # export zscore of the data
            # do sturf
            export = (
                self.df_filtered.loc[idx[:, 'area'], :]
                .apply(my_zscore, axis=1)
                .reset_index()
                )

            #     .apply(z_score, axis=1)
            #     .reset_index()
            # )


        else:
            # export some other column of data
            export = self.data.loc[self.data.Metric == level]

            export["GeneSymbol"] = export.GeneID.map(
                lambda x: self.gid_symbol.get(
                    x,
                    # _genemapper.symbol.get(x, '?')
                    _genemapper.symbol.get(str(int(x)), x),
                )
            )
            order = ["GeneID", "GeneSymbol"]  # add GeneSymbol
            order += [x for x in export.columns if x not in order and x != "Metric"]
            export[order].to_csv(outname, sep="\t", index=False)

            print("Exported", outname)


# ========================================================================================================= #

from six import string_types


class MyHeatMapper(HeatMapper):
    def _draw_data(self, ax, **kws):
        return ax.pcolormesh(self.plot_data, **kws)

    def _determine_cmap_params(self, plot_data, vmin, vmax, cmap, center, robust):
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
        self.cmap.set_bad(color="gray")
        # self.cmap.set_bad(color='white')

    def plot(self, ax, cax, kws):
        """Draw the heatmap on the provided Axes."""
        # Remove all the Axes spines
        despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        mesh = self._draw_data(
            ax, vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, **kws
        )

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Possibly add a colorbar
        if self.cbar:
            fontsize = None
            # doesn't work...
            if "fontsize" in self.cbar_kws:
                fontsize = self.cbar_kws.pop("fontsize")
            cb = ax.figure.colorbar(mesh, cax, ax, **self.cbar_kws)
            cb.outline.set_linewidth(0)
            if fontsize:
                for tick in ax.get_yticklabels():
                    tick.set_size(fontsize)
                    print(tick.get_text(), tick.get_size())
                # ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize)
            # If rasterized is passed to pcolormesh, also rasterize the
            # colorbar to avoid white lines on the PDF rendering
            if kws.get("rasterized", False):
                cb.solids.set_rasterized(True)

        # Add row and column labels
        if isinstance(self.xticks, string_types) and self.xticks == "auto":
            xticks, xticklabels = self._auto_ticks(ax, self.xticklabels, 0)
        else:
            xticks, xticklabels = self.xticks, self.xticklabels
        if isinstance(self.yticks, string_types) and self.yticks == "auto":
            yticks, yticklabels = self._auto_ticks(ax, self.yticklabels, 1)
        else:
            yticks, yticklabels = self.yticks, self.yticklabels

        ax.set(xticks=xticks, yticks=yticks)
        xtl = ax.set_xticklabels(xticklabels)
        ytl = ax.set_yticklabels(yticklabels, rotation="vertical")

        # Possibly rotate them if they overlap
        ax.figure.draw(ax.figure.canvas.get_renderer())
        if axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, mesh)
        # Invert the y axis to show the plot in matrix form
        ax.invert_yaxis()


import seaborn as sb
sb.matrix._HeatMapper = MyHeatMapper
from seaborn import heatmap
from seaborn import despine
from seaborn.matrix import _matrix_mask, axis_ticklabels_overlap


class _ScatterMapper(MyHeatMapper):
    """
    Draw a scattermap plot, similar to heatmap plot, but use scatter dots instead of heatmap
    """

    def __init__(
        self,
        data,
        marker,
        marker_size,
        vmin,
        vmax,
        cmap,
        center,
        robust,
        cbar,
        cbar_kws,
        xticklabels=True,
        yticklabels=True,
        mask=None,
    ):

        super(_ScatterMapper, self).__init__(
            data,
            vmin,
            vmax,
            cmap,
            center,
            robust,
            cbar=cbar,
            cbar_kws=cbar_kws,
            xticklabels=xticklabels,
            yticklabels=yticklabels,
            mask=mask,
            # Don't support annotation
            annot=False,
            fmt=None,
            annot_kws=None,
        )

        self.marker = marker

        if isinstance(marker_size, float) or isinstance(marker_size, int):
            self.marker_size = marker_size
        elif isinstance(marker_size, pd.DataFrame):
            self.marker_size = marker_size.loc[
                self.data.index, self.data.columns
            ].values
        else:
            self.marker_size = marker_size

    def _draw_data(self, ax, **kws):

        data = self.plot_data

        range_y = np.arange(data.shape[0], dtype=int) + 0.5
        range_x = np.arange(data.shape[1], dtype=int) + 0.45
        x, y = np.meshgrid(range_x, range_y)
        kws["rasterized"] = False
        return ax.scatter(x, y, c=data, marker=self.marker, s=self.marker_size, **kws)


def scattermap(
    data,
    marker="o",
    marker_size=100,
    vmin=None,
    vmax=None,
    cmap=None,
    center=None,
    robust=False,
    linewidths=0,
    linecolor="white",
    cbar=True,
    cbar_kws=None,
    cbar_ax=None,
    square=False,
    xticklabels="auto",
    yticklabels="auto",
    mask=None,
    ax=None,
    **kwargs
):

    plotter = _ScatterMapper(
        data,
        marker,
        marker_size,
        vmin,
        vmax,
        cmap,
        center,
        robust,
        cbar,
        cbar_kws,
        xticklabels,
        yticklabels,
        mask,
    )

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax


class MyClusterGrid(ClusterGrid):

    # def __init__(self, *args, heatmap_height_ratio=.8,
    #              dendrogram_width_ratio=.16, heatmap_width_ratio=.8,
    #              expected_size_dendrogram=1.0, expected_size_colors=0.25:
    #              **kwargs):
    def __init__(
        self,
        data,
        pivot_kws=None,
        z_score=None,
        standard_scale=None,
        figsize=None,
        row_colors=None,
        col_colors=None,
        mask=None,
        expected_size_dendrogram=1.0,
        circle_col_markers=False,
        circle_col_marker_size=60,
        force_optimal_ordering=False,
        expected_size_colors=0.25,
    ):
        """Grid object for organizing clustered heatmap input on to axes"""

        if isinstance(data, pd.DataFrame):
            self.data = data
        else:
            self.data = pd.DataFrame(data)

        self.circle_col_markers = circle_col_markers
        self.circle_col_marker_size = circle_col_marker_size

        self.data2d = self.format_data(self.data, pivot_kws, z_score, standard_scale)

        self.mask = _matrix_mask(self.data2d, mask)

        self.expected_size_dendrogram = expected_size_dendrogram
        self.expected_size_side_colors = expected_size_colors

        self.force_optimal_ordering = force_optimal_ordering

        if figsize is None:
            width, height = 10, 10
            figsize = (width, height)
        self.fig = plt.figure(figsize=figsize)

        self.row_colors, self.row_color_labels = self._preprocess_colors(
            data, row_colors, axis=0
        )
        self.col_colors, self.col_color_labels = self._preprocess_colors(
            data, col_colors, axis=1
        )

        width_ratios = self.dim_ratios(self.row_colors, figsize=figsize, axis=0)
        height_ratios = self.dim_ratios(self.col_colors, figsize=figsize, axis=1)
        nrows = 3 if self.col_colors is None else 4
        ncols = 3 if self.row_colors is None else 4
        self.gs = gridspec.GridSpec(
            nrows,
            ncols,
            wspace=0.01,
            hspace=0.01,
            width_ratios=width_ratios,
            height_ratios=height_ratios,
        )

        self.ax_row_dendrogram = self.fig.add_subplot(self.gs[nrows - 1, 0:2])
        self.ax_col_dendrogram = self.fig.add_subplot(self.gs[0:2, ncols - 1])
        self.ax_row_dendrogram.set_axis_off()
        self.ax_col_dendrogram.set_axis_off()

        self.ax_row_colors = None
        self.ax_col_colors = None

        if self.row_colors is not None:
            self.ax_row_colors = self.fig.add_subplot(self.gs[nrows - 1, ncols - 2])
        if self.col_colors is not None:
            self.ax_col_colors = self.fig.add_subplot(self.gs[nrows - 2, ncols - 1])

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
        """Get the proportions of the figure taken up by each axes"""
        figdim = figsize[axis]

        expected_size_for_dendrogram = self.expected_size_dendrogram  # Inches
        expected_size_for_side_colors = self.expected_size_side_colors  # Inches

        # Get resizing proportion of this figure for the dendrogram and
        # colorbar, so only the heatmap gets bigger but the dendrogram stays
        # the same size.
        dendrogram = expected_size_for_dendrogram / figdim

        # add the colorbar
        colorbar_width = 0.8 * dendrogram
        colorbar_height = 0.2 * dendrogram
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

    def plot(
        self,
        metric,
        method,
        colorbar_kws,
        row_cluster,
        col_cluster,
        row_linkage,
        col_linkage,
        row_color_kws=None,
        col_color_kws=None,
        annot_kws=None,
        tree_kws=None,
        **kws
    ):
        colorbar_kws = {} if colorbar_kws is None else colorbar_kws
        row_color_kws = {} if row_color_kws is None else row_color_kws
        col_color_kws = {} if col_color_kws is None else col_color_kws
        annot_kws = {} if annot_kws is None else annot_kws
        tree_kws = {} if tree_kws is None else tree_kws
        self.plot_dendrograms(
            row_cluster,
            col_cluster,
            metric,
            method,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            force_optimal_ordering=self.force_optimal_ordering,
            tree_kws=tree_kws,
        )
        try:
            xind = self.dendrogram_col.reordered_ind
        except AttributeError:
            xind = np.arange(self.data2d.shape[1])
        try:
            yind = self.dendrogram_row.reordered_ind
        except AttributeError:
            yind = np.arange(self.data2d.shape[0])

        annot = None
        if "annot" in kws:
            annot = kws.pop("annot")

        self.plot_colors(xind, yind, row_color_kws, col_color_kws, **kws)
        self.plot_matrix(
            colorbar_kws, xind, yind, annot=annot, annot_kws=annot_kws, **kws
        )
        return self

    def plot_matrix(self, colorbar_kws, xind, yind, **kws):
        self.data2d = self.data2d.iloc[yind, xind]
        self.mask = self.mask.iloc[yind, xind]

        # Try to reorganize specified tick labels, if provided
        xtl = kws.pop("xticklabels", "auto")
        try:
            xtl = np.asarray(xtl)[xind]
        except (TypeError, IndexError):
            pass
        ytl = kws.pop("yticklabels", "auto")
        try:
            ytl = np.asarray(ytl)[yind]
        except (TypeError, IndexError):
            pass

        annot = None
        if "annot" in kws and isinstance(kws["annot"], pd.DataFrame):
            ## not working for some reason:
            annot = kws.pop("annot").iloc[yind, xind]
            # annot = kws.pop('annot').loc[self.data2d.index, self.data2d.columns]
        elif "annot" in kws:
            annot = kws.pop("annot")

        heatmap(
            self.data2d,
            ax=self.ax_heatmap,
            cbar_ax=self.cax,
            cbar_kws=colorbar_kws,
            mask=self.mask,
            annot=annot,
            fmt="",
            xticklabels=xtl,
            yticklabels=ytl,
            **kws
        )

        # need to remove..
        if kws.get("cbar") == False:
            self.cax.set_xticks([])
            self.cax.set_yticks([])

        ytl = self.ax_heatmap.get_yticklabels()
        ytl_rot = None if not ytl else ytl[0].get_rotation()
        self.ax_heatmap.yaxis.set_ticks_position("right")
        self.ax_heatmap.yaxis.set_label_position("right")
        if ytl_rot is not None:
            ytl = self.ax_heatmap.get_yticklabels()
            plt.setp(ytl, rotation=ytl_rot)

    def plot_dendrograms(
        self,
        row_cluster,
        col_cluster,
        metric,
        method,
        row_linkage,
        col_linkage,
        force_optimal_ordering,
        tree_kws=None,
    ):
        # Plot the row dendrogram
        if row_cluster:
            self.dendrogram_row = dendrogram(
                self.data2d,
                metric=metric,
                method=method,
                label=False,
                axis=0,
                ax=self.ax_row_dendrogram,
                rotate=True,
                linkage=row_linkage,
                force_optimal_ordering=force_optimal_ordering,
                tree_kws=tree_kws,
            )
        else:
            self.ax_row_dendrogram.set_xticks([])
            self.ax_row_dendrogram.set_yticks([])
        # PLot the column dendrogram
        if col_cluster:
            self.dendrogram_col = dendrogram(
                self.data2d,
                metric=metric,
                method=method,
                label=False,
                axis=1,
                ax=self.ax_col_dendrogram,
                linkage=col_linkage,
                force_optimal_ordering=force_optimal_ordering,
                tree_kws=tree_kws,
            )
        else:
            self.ax_col_dendrogram.set_xticks([])
            self.ax_col_dendrogram.set_yticks([])
        despine(ax=self.ax_row_dendrogram, bottom=True, left=True)
        despine(ax=self.ax_col_dendrogram, bottom=True, left=True)

    def plot_colors(self, xind, yind, row_color_kws=None, col_color_kws=None, **kws):
        """Plots color labels between the dendrogram and the heatmap

        Parameters
        ----------
        heatmap_kws : dict
            Keyword arguments heatmap
        """
        # Remove any custom colormap and centering
        kws = kws.copy()
        kws.pop("cmap", None)
        kws.pop("center", None)
        kws.pop("vmin", None)
        kws.pop("vmax", None)
        kws.pop("robust", None)
        kws.pop("xticklabels", None)
        kws.pop("yticklabels", None)

        # Plot the row colors
        if self.row_colors is not None:
            matrix, cmap = self.color_list_to_matrix_and_cmap(
                self.row_colors, yind, axis=0
            )

            # Get row_color labels
            if self.row_color_labels is not None:
                row_color_labels = self.row_color_labels
            else:
                row_color_labels = False
            row_color_kws = row_color_kws.copy()
            for x in (
                "cmap",
                "center",
                "vmin",
                "vmax",
                "robust",
                "xticklabels",
                "yticklabels",
            ):
                row_color_kws.pop(x, None)
            full_kws = kws.copy()
            full_kws.update(row_color_kws)
            # TODO: ???
            if "cbar" in full_kws:
                full_kws.pop("cbar")

            heatmap(
                matrix,
                cmap=cmap,
                cbar=False,
                ax=self.ax_row_colors,
                xticklabels=row_color_labels,
                yticklabels=False,
                **full_kws
            )

            # Adjust rotation of labels
            if row_color_labels is not False:
                plt.setp(self.ax_row_colors.get_xticklabels(), rotation=90)
        else:
            despine(self.ax_row_colors, left=True, bottom=True)

        # Plot the column colors
        if self.col_colors is not None:
            matrix, cmap = self.color_list_to_matrix_and_cmap(
                self.col_colors, xind, axis=1
            )

            # Get col_color labels
            if self.col_color_labels is not None:
                col_color_labels = self.col_color_labels
            else:
                col_color_labels = False

            col_color_kws = col_color_kws.copy()
            for x in (
                "cmap",
                "center",
                "vmin",
                "vmax",
                "robust",
                "xticklabels",
                "yticklabels",
            ):
                col_color_kws.pop(x, None)
            full_kws = kws.copy()
            full_kws.update(col_color_kws)

            fontsize = 12
            if "fontsize" in col_color_kws:
                fontsize = col_color_kws.pop("fontsize")

            # TODO: ???
            if "cbar" in kws:
                kws.pop("cbar")
            if self.circle_col_markers:
                scattermap(
                    matrix,
                    cmap=cmap,
                    cbar=False,
                    ax=self.ax_col_colors,
                    marker_size=self.circle_col_marker_size,
                    xticklabels=False,
                    yticklabels=col_color_labels,
                    **col_color_kws,
                    **kws
                )
            else:
                heatmap(
                    matrix,
                    cmap=cmap,
                    cbar=False,
                    ax=self.ax_col_colors,
                    xticklabels=False,
                    yticklabels=col_color_labels,
                    **col_color_kws,
                    **kws
                )

            # scattermap(matrix, cmap=cmap, cbar=False, ax=self.ax_col_colors, marker_size=100,
            #            xticklabels=False, yticklabels=col_color_labels, **kws)

            # Adjust rotation of labels, place on right side
            if col_color_labels is not False:
                self.ax_col_colors.yaxis.tick_right()
                plt.setp(
                    self.ax_col_colors.get_yticklabels(), rotation=0, fontsize=fontsize
                )
        else:
            despine(self.ax_col_colors, left=True, bottom=True)


from scipy.cluster import hierarchy


class MyDendrogramPlotter(sb.matrix._DendrogramPlotter):
    def __init__(
        self,
        data,
        linkage,
        metric,
        method,
        axis,
        label,
        rotate,
        force_optimal_ordering=False,
    ):

        self.force_optimal_ordering = force_optimal_ordering
        super(sb.matrix._DendrogramPlotter, self).__init__(
            data, linkage, metric, method, axis, label, rotate
        )

    def _calculate_linkage_scipy(self):
        if np.product(self.shape) >= 10000:
            UserWarning(
                "This will be slow... (gentle suggestion: " '"pip install fastcluster")'
            )
        optimal_ordering = False
        if self.array.shape[0] < 100:
            optimal_ordering = True

        if self.force_optimal_ordering:
            optimal_ordering = True

        if optimal_ordering:
            print("optimal ordering true")

        linkage = hierarchy.linkage(
            self.array,
            method=self.method,
            metric=self.metric,
            optimal_ordering=optimal_ordering,
        )
        return linkage


sb.matrix._DendrogramPlotter = MyDendrogramPlotter


def dendrogram(
    data,
    linkage=None,
    axis=1,
    label=True,
    metric="euclidean",
    method="average",
    rotate=False,
    ax=None,
    force_optimal_ordering=False,
    tree_kws=None,
):
    """Draw a tree diagram of relationships within a matrix

    Parameters
    ----------
    data : pandas.DataFrame
        Rectangular data
    linkage : numpy.array, optional
        Linkage matrix
    axis : int, optional
        Which axis to use to calculate linkage. 0 is rows, 1 is columns.
    label : bool, optional
        If True, label the dendrogram at leaves with column or row names
    metric : str, optional
        Distance metric. Anything valid for scipy.spatial.distance.pdist
    method : str, optional
        Linkage method to use. Anything valid for
        scipy.cluster.hierarchy.linkage
    rotate : bool, optional
        When plotting the matrix, whether to rotate it 90 degrees
        counter-clockwise, so the leaves face right
    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis

    Returns
    -------
    dendrogramplotter : _DendrogramPlotter
        A Dendrogram plotter object.

    Notes
    -----
    Access the reordered dendrogram indices with
    dendrogramplotter.reordered_ind

    """
    plotter = MyDendrogramPlotter(
        data,
        linkage=linkage,
        axis=axis,
        metric=metric,
        method=method,
        label=label,
        rotate=rotate,
        force_optimal_ordering=force_optimal_ordering,
    )

    if ax is None:
        ax = plt.gca()
    return plotter.plot(ax=ax, tree_kws=tree_kws)

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


from . import utils
from .utils import *
