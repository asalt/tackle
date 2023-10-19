"""

"""
from dis import show_code
import sys
import os
from pathlib import Path
import itertools
import gseapy as gp

try:
    import rpy2
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    robjects.pandas2ri.activate()
except ModuleNotFoundError:
    print("rpy2 needs to be installed")

#  import re
#  import json
import glob
from datetime import datetime

#  import operator as op
#  from collections import OrderedDict
from functools import partial
import copy
import subprocess
from warnings import warn

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
#  # from matplotlib import gridspec
#  import matplotlib.pyplot as plt
#  from matplotlib.offsetbox import AnchoredText
import seaborn as sb

#  from seaborn.distributions import _freedman_diaconis_bins
#  import click

# try:
#     from quick import gui_option
# except ModuleNotFoundError:
#     pass

rc = {"font.family": "serif", "font.serif": ["Times", "Palatino", "serif"]}
sb.set_context("paper")
sb.set_style("white", rc)

rc = {
    "font.family": "sans-serif",
}
sb.set_context("notebook")
sb.set_style("white", rc)
sb.set_palette("muted")
sb.set_color_codes()

__version__ = "0.5.01"


GENESET_MAPPING = {
    "hallmark": "h.all.v7.0.entrez.gmt",
    "go_biological": "c5.bp.v7.0.entrez.gmt",
    "curated.all": "c2.all.v7.0.entrez.gmt",
    "curated.CGP": "c2.cgp.v7.0.entrez.gmt",
    "curated.CP.all": "c2.cp.v7.0.entrez.gmt",
    "curated.CP.BioCarta": "c2.cp.biocarta.v7.0.entrez.gmt",
    "curated.CP.KEGG": "c2.cp.kegg.v7.0.entrez.gmt",
    "curated.CP.Reactome": "c2.cp.reactome.v7.0.entrez.gmt",
    "curated.CP.PID": "c2.cp.pid.v7.0.entrez.gmt",
    "oncogenic.C6": "c6.all.v7.0.entrez.gmt",
    "go.All": "c5.all.v7.0.entrez.gmt",
    "go.Bio": "c5.bp.v7.0.entrez.gmt",
    "go.Cell": "c5.cc.v7.0.entrez.gmt",
    "go.Molecular": "c5.mf.v7.0.entrez.gmt",
    "motif.gene.sets": "c3.all.v7.0.entrez.gmt",
    "hallmark.Mm": "mh.all.v2022.1.Mm.entrez.gmt",
    "go.Bio.Mm": "m5.go.bp.v2022.1.Mm.entrez.gmt",
    "go.Cell.Mm": "m5.go.cc.v2022.1.Mm.entrez.gmt",
    "go.Molecular.Mm": "m5.go.mf.v2022.1.Mm.entrez.gmt",
}

from bcmproteomics_ext import ispec

sb.set_context("notebook", font_scale=1.4)

from .scatterplot import scatterplot
from .clusterplot import clusterplot
from .metrics import make_metrics
from .pcaplot import pcaplot
from .utils import fix_name, _get_logger
from . import utils
from .utils import *
from .containers import (
    Data,
    GeneMapper,
    get_annotation_mapper,
    get_gene_mapper,
    get_hgene_mapper,
)
from .barplot import barplot

# from cluster_to_plotly import cluster_to_plotly

sys.setrecursionlimit(10000)

TAXON_MAPPER = {"human": 9606, "mouse": 10090, "celegans": 6239}


## what is the smarter way to do this?
import logging


# def _get_logger():
#     logger = logging.getLogger(__name__)
#     logger.setLevel(logging.DEBUG)
#
#     ch = logging.StreamHandler()
#     ch.setLevel(logging.INFO)
#     # create formatter and add it to the handlers
#     formatter = logging.Formatter(
#         "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
#     )
#     ch.setFormatter(formatter)
#     # add the handlers to the logger
#     logger.addHandler(ch)
#     try:
#         fh = logging.FileHandler("tackle.log")
#         fh.setFormatter(formatter)
#         logger.addHandler(fh)
#     except PermissionError:
#         pass
#
#     # fh.setLevel(logging.DEBUG)
#     # create console handler with a higher log level
#     return logger


logger = _get_logger(__name__)


def run(data_obj):
    # if 'scatter' in plots or 'all' in data_obj.plots:
    if data_obj.make_plot("scatter"):
        g = scatterplot(
            data_obj.areas_log_shifted,
            stat=data_obj.stat,
            colors_only=data_obj.colors_only,
            shade_correlation=data_obj.shade_correlation,
        )
        outname = get_outname(
            "scatter",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            colors_only=data_obj.colors_only,
            outpath=data_obj.outpath,
        )
        save_multiple(g, outname, ".png", dpi=96)

    # ibaqs_zscore = (ibaqs_log_shifted - ibaqs_log_shifted.mean(axis=1)) / ibaqs_log_shifted.std(axis=1)
    # ibaqs_log_shifted
    # if 'cluster' in plots or 'all' in plots:
    if data_obj.make_plot("cluster"):
        g, extra_artists = clusterplot(
            data_obj.areas_log_shifted,
            highlight_gids=data_obj.highlight_gids,
            highlight_gid_names=data_obj.highlight_gid_names,
            gid_symbol=data_obj.gid_symbol,
            gene_symbols=data_obj.gene_symbols,
            z_score=data_obj.z_score,
            standard_scale=data_obj.standard_scale,
            row_cluster=data_obj.row_cluster,
            col_cluster=data_obj.col_cluster,
            metadata=data_obj.config if data_obj.show_metadata else None,
            col_data=data_obj.col_metadata,
        )
        outname = get_outname(
            "clustermap",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            colors_only=data_obj.colors_only,
            outpath=data_obj.outpath,
        )
        # print(outname)

        bbox_inches = "tight"
        if extra_artists is not None:
            bbox_inches = None
        # extra_artists=None
        save_multiple(
            g,
            outname,
            ".png",
            bbox_extra_artists=extra_artists,
            bbox_inches=bbox_inches,
        )

    # if 'pca' in plots or 'all' in plots:
    if data_obj.make_plot("pca"):
        fig, ax = pcaplot(
            data_obj.areas_log_shifted, data_obj.config, col_data=data_obj.col_metadata
        )
        outname = get_outname(
            "pcaplot",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            colors_only=data_obj.colors_only,
            outpath=data_obj.outpath,
        )
        save_multiple(fig, outname, ".png")


# def file_or_subcommand(ctx, param, value):
#     pass
class Path_or_Geneset(click.Path):
    "a real file path or the name of a valid gene set"

    def __init__(self, *args, **kwargs):
        super(Path_or_Geneset, self).__init__(*args, **kwargs)

    @classmethod
    def look_for_genes(self, value: str, geneset_str: str):
        f = GENESET_MAPPING.get(geneset_str)
        geneset_file = os.path.join(
            os.path.split(os.path.abspath(__file__))[0], "GSEA", "genesets", f
        )
        gs = gp.parser.read_gmt(geneset_file)
        matches = [x for x in gs.keys() if value in x or value in x.lower()]
        if len(matches) == 0:
            print("No matches found")
            return
        if len(matches) > 1:
            print("more than 1 match, just returning one at a time")
            return
            # todo more than one at the same time
        # return {m, gs['m'] for m in matches}
        ##
        return {matches[0]: gs[matches[0]]}

    def convert(self, value, param, ctx):
        try:
            res = click.Path.convert(self, value, param, ctx)
            return res
        except:
            # f = GENESET_MAPPING.get("hallmark")
            choices = (
                "hallmark",
                "go.Bio",
                "go.Cell",
                "go.Molecular",
                "go.All",
                "go.Cell.Mm",
                "go.Bio.Mm",
                "go.Molecular.Mm",
            )

            for choice in choices:
                out = self.look_for_genes(value, choice)
                if out is not None:
                    break
            return out
            # try hallmark first


class Path_or_Subcommand(click.Path):
    EXCEPTIONS = ("make_config", "replot_gsea")  # this one we just run

    def __init__(self, *args, **kwargs):
        super(Path_or_Subcommand, self).__init__(*args, **kwargs)

    def convert(self, value, param, ctx):
        commands = ctx.command.commands.keys()

        if value in self.EXCEPTIONS:
            return value

        if value in commands:
            help_txt = globals()[value].get_help(ctx)
            click.echo(help_txt)
            sys.exit(0)

        return click.Path.convert(self, value, param, ctx)


class int_or_ratio(click.ParamType):
    name = "integer"

    def convert(self, value, param, ctx):
        try:
            float(value)
        except ValueError:
            self.fail("%s is not a valid integer or float" % value, param, ctx)
        try:
            return int(value)
        except ValueError:  # then float
            val = float(value)
            if val > 1:
                self.fail(
                    "%s cannot be greater than 1 if specified as a float" % value,
                    param,
                    ctx,
                )
            return val


class Path_or_Glob(click.Path):
    def __init__(self, *args, **kwargs):
        super(Path_or_Glob, self).__init__(*args, **kwargs)

    def convert(self, value, param, ctx):
        try:
            return click.Path.convert(self, value, param, ctx)
        except click.BadParameter as e:
            globres = glob.glob(value)
            if not globres:
                e.show()
                sys.exit(1)


def validate_cluster_number(ctx, param, value):
    if value == "auto":
        return "auto"

    elif value == "None" or value is None:
        return None

    elif value.isdigit():
        retval = int(value)
        if retval < 1:
            raise click.BadParameter("must be set to at least 1")
        return retval

    else:
        raise click.BadParameter("must be one of `None`, `auto` or an integer")


def validate_in_config(*entries, valid_entries=None):
    import difflib

    correct = list()
    for entry in entries:
        if entry not in valid_entries:
            closest_match = difflib.get_close_matches(
                entry, valid_entries, n=1, cutoff=0.6
            )
            if closest_match:
                response2 = "Did you  mean {}".format(closest_match[0])
            else:
                response2 = ""
            click.echo(
                """{} is not in the config file. {}
            """.format(
                    entry, response2
                )
            )
        else:
            correct.append(entry)
    return correct


def validate_seed(ctx, param, value):
    if value == "None" or value is None:
        return None

    if isinstance(value, (int, float)):
        return int(value)

    elif value.isdigit():
        return int(value)

    else:
        raise click.BadParameter("Must be an integer or `None`")


# def validate_configfile(experiment_file, nonzero_subgroup=None, batch=None, group=None, covariate=None, pairs=None):
def validate_configfile(experiment_file, **kwargs):
    """
    check to ensure all metadata is present in all samples
    :experiment_file: config file to validate metadata entries
    :**kwargs: metadata entries to test
    """
    config = read_config(experiment_file)
    config = copy.deepcopy(config)  # because we're about to change it

    config_len = len([x for x in config.keys() if not x.startswith("__")])

    dunder_fields = dict(
        __PCA__=("color", "marker"),
        __norm__=("label", "group"),
        __batch__=("batch",),
    )

    for entry in config.keys():
        if entry.startswith("__"):
            dunder = config.get(entry)
            fields = dunder_fields.get(entry)
            if fields is None:
                continue
            # validate each field
            for field in fields:
                if field not in config[entry]:
                    continue  #
                count = 0
                value = config[entry][field]  # value to tally

                for k, v in config.items():
                    if k in dunder_fields:
                        continue
                    ret = v.get(value)
                    if ret is not None:
                        count += 1

                if count == 0:
                    raise click.BadParameter(
                        """Invalid entry for {}:
                    {} is not annotated in the config file""".format(
                            entry, value
                        )
                    )
                if count < config_len:
                    raise click.BadParameter(
                        """Invalid entry for {}:
                    {} is not annotated for all experiments in the config file
                    ({} / {} )""".format(
                            entry, value, count, config_len
                        )
                    )

    if "__norm__" in config.keys():
        # need to do additional validation
        norm_info = config.get("__norm__")
        control = norm_info["control"]
        group = norm_info["group"]
        label = norm_info["label"]
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
            raise click.BadParameter(
                """
            Invalid specification for `control` in the __norm__ specification.
            Expected {} to be specified by {} a total of {} times
            but is observed {} time(s).
            """.format(
                    control, label, len(groups), counter
                )
            )

    def check_group(name, name_str):
        count = 0
        missing = list()
        for k, v in config.items():
            if k in dunder_fields:
                continue
            ret = v.get(name)
            if ret is not None:
                count += 1
            else:
                missing.append(k)

        if count == 0:
            raise click.BadParameter(
                """{} is specified for {} but
            is not annotated in the config file""".format(
                    name, name_str
                )
            )
        if count < config_len:
            raise click.BadParameter(
                """{} is specified for {} but
            is not annotated for all experiments in the config file ({} / {}).
            Check entries {}
            """.format(
                    name, name_str, count, config_len, ", ".join(missing)
                )
            )

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


# ANNOTATION_CHOICES = (
#     "IDG",
#     "IO",
#     "CYTO_NUC",
#     "ER_GOLGI",
#     "SECRETED",
#     "DBTF",
#     "NUCLEUS",
#     "RTK",
#     "MATRISOME",
#     "SurfaceLabel",
#     "CellMembrane",
#     "Secreted",
#     "glycomineN",
#     "glycomineO",
#     "glycomineN/O",
#     "Membrane_Secreted",
#     "_all",
# )

ANNOTATION_CHOICES = get_annotation_mapper().categories or tuple()
ANNOTATION_CHOICES = [*ANNOTATION_CHOICES, "_all"]
# ANNOTATION_CHOICES = ["_all"]


# @gui_option
@click.group(chain=True)
@click.option(
    "--additional-info",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help=".ini file with metadata for isobaric data used for scatter and PCA plots",
)
@click.option(
    "-a",
    "--annotations",
    type=click.Choice(ANNOTATION_CHOICES),
    multiple=True,
    default=None,
    help="analyses to be performed on subsets of genes",
)
@click.option(
    "--batch",
    type=str,
    default=None,
    help="Metadata entry to group experiments for batch correction via ComBat (Requires rpy2, R, and sva installations)",
)
@click.option(
    "--batch-nonparametric",
    is_flag=True,
    default=False,
    help="Use nonparametric method for batch correction with ComBat (only used if --batch is also specified)",
)
@click.option(
    "--batch-noimputation",
    is_flag=True,
    default=False,
    help="Leave original missing values after batch correction",
)
@click.option(
    "--covariate",
    type=str,
    default=None,
    help="Metadata entry to use as covariate for batch correction via ComBat",
)
@click.option(
    "--cmap-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    help="JSON file associating specific metadata entries to pre-defined, valid colors",
)
@click.option(
    "--data-dir",
    type=click.Path(exists=False, file_okay=False),
    default="./data/",
    show_default=True,
    help="location to store and read e2g files",
)
@click.option(
    "--fill-na-zero / --no-fill-na-zero",
    show_default=True,
    default=True,
    is_flag=True,
    help="""
""",
)
@click.option(
    "--only-load-local",
    default=False,
    is_flag=True,
    show_default=True,
    help="Only try to load data locally, skip calls to iSPEC",
)
@click.option(
    "--file-format",
    type=click.Choice((".png", ".pdf", ".svg")),
    default=(".png",),
    show_default=True,
    multiple=True,
    help="File format for any plots",
)
@click.option(
    "--funcats",
    type=str,
    default=None,
    show_default=True,
    help="""Optional gene subset based on funcat or funcats,
              regular expression allowed. """,
)
@click.option(
    "--funcats-inverse",
    type=str,
    default=None,
    show_default=True,
    help="""Optional gene filtering based on funcat or funcats, regular expression allowed. """,
)
@click.option(
    "--geneids",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    help="""Optional list of geneids to subset by.
              Should have 1 geneid per line. """,
)
@click.option(
    "--group",
    type=str,
    default=None,
    help="Metadata entry to calculate p-values for differential across (Requires rpy2, R, and sva installations)",
)
@click.option(
    "--limma/ --no-limma",
    default=True,
    is_flag=True,
    help="Use limma moderated t test.",
)
@click.option(
    "--block",
    default=None,
    is_flag=False,
    help="Bio or Tech rep for block design for limma",
)
@click.option(
    "--pairs",
    type=str,
    default=None,
    help="Metadata entry that indicates sample pairs for running pairwise statistical tests.",
)
@click.option(
    "--ignore-geneids",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    help="""Optional list of geneids to ignore. Should have 1 geneid per line. """,
)
@click.option(
    "--iFOT",
    default=False,
    show_default=True,
    is_flag=True,
    help="""Calculate iFOT (divide by total input per experiment)""",
)
@click.option(
    "--iFOT-KI",
    default=False,
    show_default=True,
    is_flag=True,
    help="Calculate iFOT based on kinases",
)
@click.option(
    "--iFOT-TF",
    default=False,
    show_default=True,
    is_flag=True,
    help="Calculate iFOT based on transcription factors",
)
@click.option(
    "--genefile-norm",
    type=click.Path(exists=True, dir_okay=False),
    show_default=True,
    is_flag=False,
    help="Normalize by a list of genes",
)
@click.option(
    "--median",
    default=False,
    show_default=True,
    is_flag=True,
    help="Normalize based on median sample expression value",
)
@click.option(
    "--quantile75",
    default=False,
    show_default=True,
    is_flag=True,
    help="Normalize based on 75 pct quantile sample expression value",
)
@click.option(
    "--quantile90",
    default=False,
    show_default=True,
    is_flag=True,
    help="Normalize based on 90 pct quantile sample expression value",
)
@click.option(
    "--normalize-across-species", is_flag=True, default=False, show_default=True
)
@click.option(
    "--impute-missing-values / --no-impute-missing-values",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "-n",
    "--name",
    type=str,
    default="",
    help="An optional name for the analysis that will place all results in a subfolder.",
)
@click.option(
    "--result-dir",
    type=click.Path(exists=False, file_okay=False),
    default="./results",
    show_default=True,
    help="Base directory to store results. Will be created if does not exist.",
)
@click.option("--SRA", type=click.Choice(["S", "R", "A"]))
@click.option(
    "--number-sra",
    default=1,
    help="number of `SRA` level gene products required",
    show_default=True,
)
@click.option(
    "--taxon",
    type=click.Choice(["human", "mouse", "celegans", "all"]),
    default="all",
    show_default=True,
)
@click.option(
    "--non-zeros",
    default=1,
    show_default=True,
    type=int_or_ratio(),
    help="""Minimum number of non zeros OR fraction of nonzeros allowed for each gene product
              across samples. If a decimal is specified (e.g. 1.0), this indicates 100% of values are nonzero.
              If an integer is specified (1), this indicates that 1 value is nonzero.
              """,
)
@click.option(
    "--unique-pepts",
    default=0,
    show_default=True,
    type=int,
    help="number of unique peptides required per gene. Also must satisfy the nonzero argument.",
)
@click.option("--ref-group-name", default=None)
@click.option("--ref-control-channel", default=None)
@click.option("--nonzero-subgroup", type=str, default=None, help="")
# @click.argument('experiment_file', type=click.Path(exists=True, dir_okay=False))
@click.argument("experiment_file", type=Path_or_Subcommand(exists=True, dir_okay=False))
@click.pass_context
def main(
    ctx,
    additional_info,
    annotations,
    batch,
    batch_nonparametric,
    batch_noimputation,
    covariate,
    cmap_file,
    data_dir,
    only_load_local,
    file_format,
    fill_na_zero,
    funcats,
    funcats_inverse,
    geneids,
    group,
    impute_missing_values,
    limma,
    block,
    pairs,
    ignore_geneids,
    ifot,
    ifot_ki,
    ifot_tf,
    genefile_norm,
    median,
    quantile75,
    quantile90,
    name,
    normalize_across_species,
    result_dir,
    taxon,
    non_zeros,
    nonzero_subgroup,
    unique_pepts,
    sra,
    number_sra,
    experiment_file,
    ref_group_name,
    ref_control_channel,
):
    """"""
    # name, taxon, non_zeros, experiment_file):
    norm_info = None
    if ref_control_channel and ref_group_name:
        norm_info = dict(
            control=ref_control_channel, group=ref_group_name, label="label"
        )

    if experiment_file in Path_or_Subcommand.EXCEPTIONS:
        # then is it actually a subcommand (only make_config right now)
        # we wish to run on its own, without loading data

        # hacky
        if experiment_file == "make_config":
            parser = make_config.make_parser(ctx)
            _cmd = make_config
        elif experiment_file == "replot_gsea":
            parser = replot_gsea.make_parser(ctx)
            _cmd = replot_gsea
        the_kws = parser.parse_args(sys.argv[2:])[0]
        ctx.invoke(_cmd, **the_kws)
        ctx.exit(0)

    if not limma:
        raise click.BadOptionUsage(
            "limma", "At the moment, only use of `limma` is supported"
        )

    sumed_norm_flags = ifot + ifot_ki + ifot_tf + median + quantile75 + quantile90
    if (sumed_norm_flags > 1) or (sumed_norm_flags > 0 and genefile_norm):
        raise click.BadParameter(
            "Cannot specify a combination of `iFOT`, `iFOT-KI`, `iFOT-TF`, `median`, `genefile_norm`"
        )

    # validate_subgroup(nonzero_subgroup, experiment_file)
    validate_configfile(
        experiment_file,
        nonzero_subgroup=nonzero_subgroup,
        batch=batch,
        group=group,
        covariate=covariate,
        pairs=pairs,
    )

    metrics, metrics_after_filter, metrics_unnormed_area = False, False, True
    if (
        "metrics" in sys.argv
    ):  # know this ahead of time, accumulate metrics during data load
        metrics = True
        if "--after-filter" in sys.argv:
            metrics_after_filter = True  #
        if "--after-norm" in sys.argv:
            metrics_unnormed_area = False

    # for keeping track of what to stack from data
    cluster_annotate_cols = None

    # extract cluster annotation column for cluster if present
    # and not if --annotate is specified for pca
    # this is needed now so we keep the proper column when loading data
    # the logic should work!

    if "--annotate" in sys.argv and any(
        x in sys.argv for x in ("cluster", "cluster2", "gsea")
    ):
        cluster_annotate_cols = list()
        # _annot_args = [i for i, x in enumerate(sys.argv) if x.strip() == "--annotate"]
        _annot_args = [i for i, x in enumerate(sys.argv) if x.strip() == "--annotate"]
        for i in _annot_args:
            try:
                _q = sys.argv[i + 1]
            except IndexError:
                continue
            if _q in [
                "PSMs",
                "PSMs_u2g",
                "PeptideCount",
                "PeptideCount_S",
                "PeptideCount_S_u2g",
                "PeptideCount_u2g",
                "SRA",
            ]:
                cluster_annotate_cols.append(_q)
    if "--level" in sys.argv:  # for data export
        if cluster_annotate_cols is None:
            cluster_annotate_cols = list()
        _level_args = [i for i, x in enumerate(sys.argv) if x.strip() == "--level"]
        for i in _level_args:
            try:
                _q = sys.argv[i + 1]
            except IndexError:
                continue
            if _q in [
                "PSMs",
                "PSMs_u2g",
                "PeptideCount",
                "PeptideCount_S",
                "PeptideCount_S_u2g",
                "PeptideCount_u2g",
                "SRA",
            ]:
                cluster_annotate_cols.append(_q)
    if cluster_annotate_cols is not None:
        cluster_annotate_cols = list(set(cluster_annotate_cols))

    if annotations is not None and "_all" in annotations:
        annotations = [x for x in ANNOTATION_CHOICES if x not in ("_all")]

    export_all = False
    if all(x in sys.argv for x in ("export", "--level")) and any(
        x in sys.argv
        for x in ("all", "align", "MSPC")
        # x in sys.argv for x in ("all", "align")
    ):
        # know this ahead of time, calculate more things during data load
        export_all = True

    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    if ctx.obj is None:
        ctx.obj = dict()

    ctx.obj["file_fmts"] = file_format

    analysis_name = get_file_name(experiment_file)
    if analysis_name is None:
        print("Error parsing configfile name.")
        analysis_name = "Unnamed"

    now = datetime.now()
    # context = click.get_current_context() # same thing as ctx
    params = ctx.params

    data_obj = Data(
        additional_info=additional_info,
        annotations=annotations,
        batch=batch,
        batch_nonparametric=batch_nonparametric,
        batch_noimputation=batch_noimputation,
        covariate=covariate,
        cmap_file=cmap_file,
        data_dir=data_dir,
        base_dir=result_dir,
        fill_na_zero=fill_na_zero,
        funcats=funcats,
        funcats_inverse=funcats_inverse,
        geneids=geneids,
        group=group,
        pairs=pairs,
        ifot=ifot,
        ifot_ki=ifot_ki,
        ifot_tf=ifot_tf,
        median=median,
        quantile75=quantile75,
        quantile90=quantile90,
        genefile_norm=genefile_norm,
        name=name,
        non_zeros=non_zeros,
        nonzero_subgroup=nonzero_subgroup,
        unique_pepts=unique_pepts,
        taxon=taxon,
        impute_missing_values=impute_missing_values,
        normalize_across_species=normalize_across_species,
        experiment_file=experiment_file,
        metrics=metrics,
        limma=limma,
        block=block,
        SRA=sra,
        number_sra=number_sra,
        metrics_after_filter=metrics_after_filter,
        metrics_unnormed_area=metrics_unnormed_area,
        ignore_geneids=ignore_geneids,
        export_all=export_all,
        cluster_annotate_cols=cluster_annotate_cols,
        only_local=only_load_local,
        norm_info=norm_info,
    )

    outname = (
        get_outname(
            "metadata",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            colors_only=data_obj.colors_only,
            batch=data_obj.batch_applied,
            batch_method="parametric"
            if not data_obj.batch_nonparametric
            else "nonparametric",
            outpath=data_obj.outpath,
        )
        + ".tab"
    )
    data_obj.col_metadata.to_csv(outname, sep="\t")

    taxon_ratios = data_obj.taxon_ratios
    if (
        not taxon_ratios["9606"].nunique()
        == taxon_ratios["10090"].nunique()
        == taxon_ratios["9031"].nunique()
    ):
        outname = (
            get_outname(
                "taxon_ratios",
                name=data_obj.outpath_name,
                taxon=data_obj.taxon,
                non_zeros=data_obj.non_zeros,
                colors_only=data_obj.colors_only,
                batch=data_obj.batch_applied,
                batch_method="parametric"
                if not data_obj.batch_nonparametric
                else "nonparametric",
                outpath=data_obj.outpath,
            )
            + ".tab"
        )
        taxon_ratios.to_csv(outname, sep="\t")

    # cf = 'correlatioplot_args_{}.json'.format(now.strftime('%Y_%m_%d_%H_%M_%S'))
    # with open(os.path.join(data_obj.outpath, cf), 'w') as f:
    #     json.dump(params, f)
    outname = (
        get_outname(
            "context",
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            colors_only=data_obj.colors_only,
            batch=data_obj.batch_applied,
            batch_method="parametric"
            if not data_obj.batch_nonparametric
            else "nonparametric",
            outpath=data_obj.outpath,
        )
        + ".tab"
    )
    params = dict(ctx.params)
    params["file_format"] = " | ".join(params["file_format"])
    params["annotations"] = " | ".join(params["annotations"])
    param_df = pd.DataFrame(params, index=[analysis_name]).T
    param_df.to_csv(outname, sep="\t")

    ctx.obj["data_obj"] = data_obj


# @main.command("make-jk")
# @click.pass_context
# def convert(context,):


@main.command("make_config")
@click.option(
    "--delimiter",
    type=str,
    default=None,
    show_default=True,
    help="Delimiter used in the input file.",
)
@click.option(
    "--excel",
    default=False,
    is_flag=True,
    show_default=True,
    help="""Input file is an excel sheet.""",
)
@click.option(
    "--excel-sheetnumber",
    default=0,
    type=int,
    show_default=True,
    help="""Sheet number to use if input is excel file. For use with `--excel`.""",
)
@click.option(
    "--infer-inputfile/--no--infer-inputfile",
    default=True,
    is_flag=True,
    show_default=True,
    help="""Automatically try to infer input file and delimiter based on file extension.
              Note explicit use of `--delimiter` and `--excel` will overwrite this behavior.
              """,
)
@click.option(
    "--output",
    default=None,
    show_default=True,
    help="Name and location to save output. Default is same location as input with same basename.",
)
@click.argument("inputfile", type=click.Path(exists=True, dir_okay=False))
@click.pass_context
def make_config(
    ctx,
    delimiter,
    excel,
    excel_sheetnumber,
    infer_inputfile,
    output,
    inputfile,
    help=False,
):
    """
    Convert a delimited (tab, csv, excel, etc) file into a config file.
    """
    if help is True or (help is False and inputfile is None):
        help_txt = globals()["make_config"].get_help(ctx)
        print(help_txt)
        sys.exit(0)

    from .Ini2Csv.main import Csv2Conf

    SEP_DICT = {
        "tab": "\t",
        "tsv": "\t",
        ".ssv": ";",
        "csv": ",",
        "xls": "excel",
        "xlsx": "excel",
    }

    filefront, ext = os.path.splitext(inputfile)
    ext = ext[1:]

    if infer_inputfile and delimiter is None:
        sep = SEP_DICT.get(ext)
        if sep is None:
            warn(
                "Separator not known for extension {}, falling back to default tab".format(
                    ext
                )
            )
            sep = "\t"

    elif delimiter and not excel:
        sep = delimiter

    elif excel:
        sep = "excel"

    if output is not None:
        out = output
    else:
        out = filefront + ".conf"

    c = Csv2Conf(inputfile, sep=sep, excel_sheetnumber=excel_sheetnumber)
    print("Writing {}".format(out), end=" ...", flush=True)
    c.export(out)
    print("done")


@main.command("scatter")
@click.option(
    "--colors-only",
    default=False,
    is_flag=True,
    show_default=True,
    help="Only plot colors on correlationplot, no datapoints",
)
@click.option(
    "--shade-correlation/--no-shade-correlation",
    default=True,
    is_flag=True,
    show_default=True,
    help="",
)
@click.option(
    "--stat",
    type=click.Choice(["pearson", "spearman"]),
    default="pearson",
    show_default=True,
)
@click.pass_context
def scatter(ctx, colors_only, shade_correlation, stat):
    data_obj = ctx.obj["data_obj"]
    file_fmts = ctx.obj["file_fmts"]

    _ = data_obj.areas_log_shifted
    outname = get_outname(
        "scatter",
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=colors_only,
        batch=data_obj.batch_applied,
        stat=stat,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        outpath=data_obj.outpath,
    )

    X = data_obj.areas_log_shifted.copy()
    # to_mask = (data_obj.mask | (data_obj.areas_log==-10))
    # ??
    # to_mask = (data_obj.mask | (data_obj.areas_log == data_obj.minval_log)).dropna(1)
    to_mask = (data_obj.mask | (data_obj.areas_log == data_obj.areas_log.min())).dropna(
        1
    )
    X[to_mask] = np.NaN
    # g = scatterplot(X.replace(0, np.NaN), stat=stat,
    g = scatterplot(
        X,
        stat=stat,
        colors_only=colors_only,
        shade_correlation=shade_correlation,
        outname=outname,
        file_fmts=file_fmts,
        mask=data_obj.mask,
    )
    if g is None:  # plotted and saved via R
        return

    save_multiple(g, outname, *file_fmts, dpi=96)


@main.command("export")
@click.option(
    "--level",
    type=click.Choice(
        [
            "all",
            "align",
            "area",
            "MSPC",
            "SRA",
            "PeptideCount",
            "PeptideCount_u2g",
            "AreaSum_dstrAdj",
            "AreaSum_u2g_0",
            "AreaSum_u2g_max",
            "zscore",
            "gct",
        ]
    ),
    default=("area",),
    multiple=True,
    help="""Export data table of the filtered list of gene products used for plotting
              `all` returns all the data in long format
              `align` returns all data formatted for import into align!
              `SRA` returns data matrix with SRA values per gene product
              """,
)
@click.option(
    "--linear",
    default=False,
    show_default=True,
    is_flag=True,
    help="Export linear (not logged) values when exporting as area",
)
@click.option(
    "--genesymbols",
    default=False,
    is_flag=True,
    show_default=True,
    help="alias for --gene-symbols",
)
@click.option("--gene-symbols", default=False, is_flag=True, show_default=True)
@click.pass_context
def export(ctx, level, genesymbols, gene_symbols, linear):
    # =====

    # =====

    data_obj = ctx.obj["data_obj"]
    for l in level:
        data_obj.perform_data_export(
            l, genesymbols=genesymbols or gene_symbols, linear=linear
        )


@main.command("cluster")
@click.option(
    "--annotate",
    type=click.Choice(
        [
            "PSMs",
            "PSMs_u2g",
            "PeptideCount",
            "PeptideCount_S",
            "PeptideCount_S_u2g",
            "PeptideCount_u2g",
            "SRA",
        ]
    ),
    default=None,
    show_default=True,
)
@click.option(
    "--col-cluster/--no-col-cluster",
    default=True,
    is_flag=True,
    show_default=True,
    help="""Cluster columns via hierarchical clustering.
              Note this is overridden by specifying `nclusters`""",
)
@click.option("--cmap", default=None, show_default=True)
@click.option("--circle-col-markers", is_flag=True, default=False, show_default=True)
@click.option("--circle-col-marker-size", default=60, show_default=True)
@click.option(
    "--figsize",
    nargs=2,
    type=float,
    default=None,
    show_default=True,
    help="""Optionally specify the figuresize (width, height) in inches
              If not specified, tries to use a reasonable default depending on the number of
              samples.
              """,
)
@click.option(
    "--force-optimal-ordering",
    is_flag=True,
    default=False,
)
@click.option(
    "--genefile",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.option(
    "--force-plot-genes",
    is_flag=True,
    default=False,
    help="""
              Will force display of all genes in `genefile` as missing values
              """,
)
@click.option(
    "--gene-symbols",
    default=False,
    is_flag=True,
    show_default=True,
    help="Show Gene Symbols on clustermap",
)
@click.option(
    "--gene-symbol-fontsize", default=8, show_default=True, help="Gene Symbol font size"
)
@click.option(
    "--highlight-geneids",
    # type=click.Path(exists=True, dir_okay=False),
    type=Path_or_Geneset(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=True,
    help="""Optional list of geneids to highlight by.
              Should have 1 geneid per line. """,
)
@click.option(
    "--linear", default=False, is_flag=True, help="Plot linear values (default log10)"
)
@click.option(
    "--legend-include",
    type=str,
    multiple=True,
    help="""Specific entries in the config file to include in the legend.
              (Default all are included)""",
)
@click.option(
    "--legend-exclude",
    type=str,
    multiple=True,
    help="""Specific entries in the config file to ignore for the legend.
              (Default all are included)""",
)
@click.option(
    "--linkage",
    type=click.Choice(
        ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
    ),
    default="ward",
    show_default=True,
    help="linkage method for hierarchical clustering",
)
@click.option(
    "--max-autoclusters",
    default=30,
    show_default=True,
    help="""Max number of clusters to try
when `auto` is set for `--nclusters`""",
)
@click.option(
    "--nclusters",
    default=None,
    callback=validate_cluster_number,
    show_default=True,
    help="""If specified by an integer, use that number of clusters via k-means clustering. If specified as `auto`, will try to find the optimal number of clusters""",
)
@click.option(
    "--dbscan",
    default=False,
    is_flag=True,
    show_default=True,
    help="""Use DBSCAN algorithm to find and cluster data. Cannot be used with nclusters specification""",
)
@click.option(
    "--row-cluster/--no-row-cluster",
    default=True,
    is_flag=True,
    show_default=True,
    help="Cluster rows via hierarchical clustering",
)
@click.option("--order-by-abundance", default=False, is_flag=True, show_default=True)
@click.option(
    "--seed",
    default=1234,
    help="seed for kmeans clustering",
    callback=validate_seed,
    show_default=True,
)
@click.option(
    "--show-metadata/--hide-metadata",
    default=True,
    show_default=True,
    is_flag=True,
    help="""Show metadata on clustermap if present""",
)
@click.option(
    "--standard-scale",
    type=click.Choice(["None", "0", "1"]),
    default="None",
    show_default=True,
)
@click.option(
    "--show-missing-values/--hide-missing-values",
    default=True,
    is_flag=True,
    show_default=True,
    help="""Whether or not to show missing values on the cluster plot and missing values""",
)
@click.option("--square", default=False, is_flag=True, help="Force square heatmap")
@click.option(
    "--z-score", type=click.Choice(["None", "0", "1"]), default="0", show_default=True
)
@click.option("--z-score-by", type=str, default=None, show_default=True)
@click.pass_context
def cluster(
    ctx,
    annotate,
    cmap,
    circle_col_markers,
    circle_col_marker_size,
    col_cluster,
    dbscan,
    figsize,
    force_optimal_ordering,
    force_plot_genes,
    genefile,
    gene_symbols,
    gene_symbol_fontsize,
    highlight_geneids,
    linear,
    legend_include,
    legend_exclude,
    linkage,
    max_autoclusters,
    nclusters,
    order_by_abundance,
    row_cluster,
    seed,
    show_metadata,
    standard_scale,
    show_missing_values,
    square,
    z_score,
    z_score_by,
):
    if order_by_abundance and row_cluster:
        print("Not ")

    if not figsize:  # returns empty tuple if not specified
        figsize = None

    if nclusters is not None and dbscan:
        raise click.BadOptionUsage("Cannot specify `nclusters` and use DBSCAN")

    genes = None
    if genefile:
        genes = parse_gid_file(genefile)

    data_obj = ctx.obj["data_obj"]

    X = data_obj.areas_log_shifted
    if linear:
        X = 10**X

    # X = data_obj.areas_log_shifted.copy()
    # X[data_obj.areas == 0] = 0 # fill the zeros back
    # X[data_obj.mask] = np.NaN

    data_obj.set_highlight_gids(highlight_geneids)
    data_obj.standard_scale = data_obj.clean_input(standard_scale)
    data_obj.z_score = data_obj.clean_input(z_score)

    col_meta = data_obj.col_metadata
    _expids = ("recno", "runno", "searchno")

    # valid_entries = validate_in_config(*legend_include, *legend_exclude, valid_entries=col_meta.index)
    valid_entries = validate_in_config(
        *legend_include, *legend_exclude, valid_entries=col_meta.columns
    )
    legend_include = [
        x for x in legend_include if x in (set(legend_include) & set(valid_entries))
    ]
    legend_exclude = set(legend_exclude) & set(valid_entries)
    # col_meta = col_meta.loc[[x for x in col_meta.index if x not in _expids]]
    col_meta = col_meta[[x for x in col_meta.columns if x not in _expids]]

    annot_mat = None
    if annotate:
        # annot_mat = (data_obj.data.loc[ idx[X.index.tolist(), annotate], : ]
        _cols = [x for x in data_obj.data.columns if x not in ("Metric")]
        annot_mat = (
            data_obj.data.loc[
                (data_obj.data.GeneID.isin(X.index.tolist()))
                & (data_obj.data.Metric == annotate)
            ][_cols]
            # .reset_index(level=1, drop=True)
            .set_index("GeneID")
            .fillna(0)
            .astype(int)
        )
        annot_mat.columns.name = annotate

    result = clusterplot(
        X,
        annot_mat=annot_mat or robjects.NULL,
        cmap_name=cmap,
        highlight_gids=data_obj.highlight_gids,
        highlight_gid_names=data_obj.highlight_gid_names,
        force_plot_genes=force_plot_genes,
        genes=genes,
        gid_symbol=data_obj.gid_symbol,
        gene_symbols=gene_symbols,
        z_score=data_obj.z_score,
        standard_scale=data_obj.standard_scale,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        # metadata=data_obj.config if show_metadata else None,
        col_data=col_meta if show_metadata else None,
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
        legend_exclude=legend_exclude,
        order_by_abundance=order_by_abundance,
        seed=seed,
        metadata_colors=data_obj.metadata_colors,
        circle_col_markers=circle_col_markers,
        circle_col_marker_size=circle_col_marker_size,
        square=square,
        z_score_by=z_score_by,
    )

    g = result["clustermap"]["clustergrid"]
    extra_artists = result["clustermap"]["extra_artists"]
    missing_values = "masked" if show_missing_values else "unmasked"
    outname_func = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        outpath=data_obj.outpath,
        missing_values=missing_values,
    )

    kmeans_res = result.get("kmeans")
    dbscan_res = result.get("dbscan")
    if kmeans_res is not None:
        kmeans_clusters = kmeans_res["nclusters"]
        outname = outname_func("clustermap_{}clusters".format(kmeans_clusters))
    elif dbscan_res is not None:
        dbscan_clusters = dbscan_res["nclusters"]
        outname = outname_func("clustermap_{}clusters".format(dbscan_clusters))
    else:
        outname = outname_func("clustermap")

    bbox_inches = "tight"
    if extra_artists is not None:
        bbox_inches = None
        # extra_artists=None
    # save_multiple(g, outname, '.png', bbox_extra_artists=extra_artists, bbox_inches=bbox_inches)
    file_fmts = ctx.obj["file_fmts"]
    save_multiple(
        g,
        outname,
        *file_fmts,
        bbox_extra_artists=extra_artists,
        bbox_inches=bbox_inches,
    )

    if kmeans_res is not None:
        fig = kmeans_res["cluster_center_plot"]["fig"]
        outname = outname_func("{}clusters_centers".format(kmeans_clusters))
        save_multiple(fig, outname, *file_fmts)
        plt.close(fig)

        kmeans_data = kmeans_res["data"]
        outname = os.path.abspath(
            outname_func("{}clusters_labels".format(kmeans_clusters)) + ".tsv"
        )
        kmeans_data.to_csv(outname, index=True, sep="\t", na_rep="NaN")
        print("Saved:", outname)

        fig = kmeans_res["auto"].get("fig")
        if fig is not None:
            outname = outname_func("cluster_optimized_results")
            save_multiple(fig, outname, *file_fmts)
            plt.close(fig)

        fig = kmeans_res["silhouette"].get("fig")
        if fig is not None:
            outname = outname_func("{}clusters_silhouette".format(kmeans_clusters))
            save_multiple(fig, outname, *file_fmts)
            plt.close(fig)

    if dbscan_res is not None:
        dbscan_data = dbscan_res["data"]

        outname = os.path.abspath(
            outname_func("{}clusters_labels".format(dbscan_clusters)) + ".tab"
        )
        dbscan_data.to_csv(outname, index=True, sep="\t")
        print("Saved:", outname)

        fig = dbscan_res["silhouette"].get("fig")
        if fig is not None:
            outname = outname_func("{}clusters_silhouette".format(dbscan_clusters))
            save_multiple(fig, outname, *file_fmts)
            plt.close(fig)


@main.command("pca")
@click.option(
    "--annotate",
    default=False,
    show_default=True,
    is_flag=True,
    help="Annotate points on PC plot",
)
@click.option(
    "--max-pc",
    default=2,
    show_default=True,
    help="Maximum PC to plot. Plots all combinations up to this number.",
)
@click.option(
    "--color",
    default="",
    show_default=True,
    is_flag=False,
    help="What meta entry to color PCA",
)
@click.option(
    "--marker",
    default="",
    show_default=True,
    is_flag=False,
    help="What meta entry to mark PCA",
)
@click.option(
    "--genefile",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.pass_context
def pca(ctx, annotate, max_pc, color, marker, genefile):
    data_obj = ctx.obj["data_obj"]

    genes = None
    if genefile:
        genes = parse_gid_file(genefile)

    # # fig, ax = pcaplot(data_obj.areas_log_shifted, data_obj.config, col_data = data_obj.col_metadata)

    figs = pcaplot(
        data_obj.areas_log_shifted,
        metadata=data_obj.config,
        col_data=data_obj.col_metadata,
        annotate=annotate,
        max_pc=max_pc,
        color_label=color,
        marker_label=marker,
        genes=genes,
        metadata_colors=data_obj.metadata_colors,
    )

    # outname_func = get_outname('pcaplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
    #                            non_zeros=data_obj.non_zeros, batch=data_obj.batch_applied,
    #                            batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
    #                            colors_only=data_obj.colors_only, outpath=data_obj.outpath,
    # )

    outname_func = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        colors_only=data_obj.colors_only,
        outpath=data_obj.outpath,
    )

    file_fmts = ctx.obj["file_fmts"]
    for name, fig in figs.items():
        outname = outname_func(name)
        save_multiple(fig, outname, *file_fmts)


@main.command("pca2")
@click.option(
    "--annotate",
    default=False,
    show_default=True,
    is_flag=True,
    help="Annotate points on PC plot",
)
@click.option("--frame", is_flag=True, default=False, show_default=True)
@click.option(
    "--normalize-by",
    default=None,
    show_default=True,
    help="Metadata group to normalize by",
)
@click.option("--center / --no-center", is_flag=True, default=True, show_default=True)
@click.option("--scale / --no-scale", is_flag=True, default=False, show_default=True)
# @click.option('--annotate', type=click.Choice(['t', 'norm', 'k'])
@click.option(
    "--max-pc",
    default=2,
    show_default=True,
    help="Maximum PC to plot. Plots all combinations up to this number.",
)
@click.option(
    "--color",
    default=None,
    show_default=True,
    is_flag=False,
    help="What meta entry to color PCA",
)
@click.option(
    "--marker",
    default=None,
    show_default=True,
    is_flag=False,
    help="What meta entry to mark PCA",
)
@click.option(
    "--genefile",
    type=Path_or_Geneset(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.pass_context
def pca2(
    ctx, annotate, frame, normalize_by, center, scale, max_pc, color, marker, genefile
):
    outname_kws = dict()

    # try:
    #     import rpy2.robjects as robjects
    #     from rpy2.robjects import pandas2ri
    #     from rpy2.robjects.packages import importr
    # except ModuleNotFoundError:
    #     print("rpy2 needs to be installed")
    #     return

    data_obj = ctx.obj["data_obj"]

    X = data_obj.areas_log_shifted
    X.index = X.index.astype(str)
    col_meta = data_obj.col_metadata.copy()

    genes = None
    if genefile:
        genes = parse_gid_file(genefile)
        _tokeep = [x for x in genes if x in X.index]  # preserve order
        X = X.loc[_tokeep]
        outname_kws["genefile"] = fix_name(os.path.splitext(genefile)[0])

    # ======================================
    outname_func = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        center=center,
        scale=scale,
        norm_by=normalize_by,
        # outpath=data_obj.outpath,
        outpath=os.path.join(data_obj.outpath, "pca"),
        **outname_kws,
    )
    # ======================================

    # now convert to correct orientation
    # variable <color> <shape> gene1 gene2 gene3 ...

    # dfm = df.melt(id_vars=['GeneID', 'GeneSymbol'])
    dfm = (
        X.fillna(0)
        .reset_index()
        .melt(id_vars=["GeneID"])
        .merge(col_meta, left_on="variable", right_index=True)
    )
    # dfm.to_csv(outname_func('pca_input')+'.tsv', sep='\t', index=False)
    if color in dfm:
        dfm[color] = dfm[color].astype(str)
    if marker in dfm:
        dfm[marker] = dfm[marker].astype(str)

    robjects.pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], "R", "pcaplot.R")
    r_source(r_file)
    pca2 = robjects.r["pca2"]

    # # columns to assign colors:
    # metadata_color_list = list()
    # for metacat in (color, marker):
    #     if not metacat:
    #         continue
    #     if data_obj.metadata_colors is not None and metacat in data_obj.metadata_colors:
    #         mapping = data_obj.metadata_colors[metacat]
    #         # entry = robjects.vectors.ListVector({metacat: robjects.vectors.ListVector(mapping)})
    #         metadata_color_list.append(mapping)
    #     else:
    #         ncolors = col_meta[metacat].nunique()
    #         cmap = sb.color_palette(n_colors=ncolors)
    #         color_iter = map(mpl.colors.rgb2hex, cmap)
    #         themapping = {x: c for x,c in zip(col_meta[metacat].unique(), color_iter)}
    #         # entry = robjects.vectors.ListVector({metacat: robjects.vectors.ListVector(themapping)})
    #         metadata_color_list.append(themapping)
    # metadata_colorframe = pd.DataFrame(metadata_color_list)

    # columns to assign colors:
    metadata_color_list = list()
    # for metacat in (color, marker):
    if not color:
        logger.warning(f"Color not specified. Using recno as color variable.")
        # raise ValueError("must specify color")
        # print(data_obj.col_metadata.columns)
        # color = data_obj.col_metadata.columns[0]
        color = "recno"

    if color:
        if data_obj.metadata_colors is not None and color in data_obj.metadata_colors:
            mapping = data_obj.metadata_colors[color]
            # entry = robjects.vectors.ListVector({metacat: robjects.vectors.ListVector(mapping)})
            keys = col_meta[color].unique()
            mapping = {k: mapping[k] for k in keys}
            color_list = robjects.vectors.ListVector(mapping)
        else:
            ncolors = col_meta[color].nunique()
            cmap = sb.color_palette(palette="bright", n_colors=ncolors)
            color_iter = map(mpl.colors.rgb2hex, cmap)
            themapping = {x: c for x, c in zip(col_meta[color].unique(), color_iter)}
            # entry = robjects.vectors.ListVector({metacat: robjects.vectors.ListVector(themapping)})
            color_list = robjects.vectors.ListVector(themapping)

    # if marker:
    #     if data_obj.metadata_colors is not None and marker in data_obj.metadata_colors:
    #         mapping = data_obj.metadata_colors[marker]
    #         # entry = robjects.vectors.ListVector({metacat: robjects.vectors.ListVector(mapping)})
    #         marker_list = robjects.vectors.ListVector(mapping)
    #     else:
    #         ncolors = col_meta[marker].nunique()
    #         cmap = sb.color_palette(n_colors=ncolors)
    #         color_iter = map(mpl.colors.rgb2hex, cmap)
    #         themapping = {x: c for x,c in zip(col_meta[marker].unique(), color_iter)}
    #         # entry = robjects.vectors.ListVector({metacat: robjects.vectors.ListVector(themapping)})
    #         marker_list = robjects.vectors.ListVector(mapping)

    # metadata_colorframe = pd.DataFrame(metadata_color_list)

    file_fmts = ctx.obj["file_fmts"]
    pca2(
        dfm,
        color=color or robjects.NULL,
        shape=marker or robjects.NULL,
        outfiletypes=np.array(file_fmts),
        outname=outname_func("pca2"),
        center=center,
        scale=scale,
        normalize_by=normalize_by or robjects.NULL,
        label=annotate,
        showframe=frame,
        max_pc=max_pc,
        color_list=color_list,
        marker_list=robjects.NULL,
        title=outname_kws.get("genefile") or robjects.NULL,
        annot_str=outname_kws.get("genefile") or robjects.NULL,
    )


@main.command("cluster2")
@click.option(
    "--annotate",
    type=click.Choice(
        [
            "PSMs",
            "PSMs_u2g",
            "PeptideCount",
            "PeptideCount_S",
            "PeptideCount_S_u2g",
            "PeptideCount_u2g",
            "SRA",
        ]
    ),
    default=None,
    show_default=True,
)
@click.option(
    "--cluster-col-slices/--no-cluster-col-slices",
    default=True,
    is_flag=True,
    show_default=True,
    help="whether or not to cluster col slices",
)
@click.option(
    "--cluster-row-slices/--no-cluster-row-slices",
    default=True,
    is_flag=True,
    show_default=True,
    help="whether or not to cluster row slices",
)
@click.option(
    "--col-cluster/--no-col-cluster",
    default=True,
    is_flag=True,
    show_default=True,
    help="""Cluster columns via hierarchical clustering.
              Note this is overridden by specifying `nclusters`""",
)
@click.option("--cmap", default=None, show_default=True)
@click.option(
    "--cut-by",
    default=None,
    show_default=True,
    help="""metadata variable to cut the heatmap
              into segments""",
)
# @click.option('--circle-col-markers', is_flag=True, default=False, show_default=True)
# @click.option('--circle-col-marker-size', default=60, show_default=True)
@click.option("--color-low", default="blue", show_default=True)
@click.option("--color-mid", default="white", show_default=True)
@click.option("--color-high", default="red", show_default=True)
@click.option(
    "--figsize",
    nargs=2,
    type=float,
    default=None,
    show_default=True,
    help="""Optionally specify the figuresize (width, height) in inches
              If not specified, tries to use a reasonable default depending on the number of
              samples.
              """,
)
@click.option("--add-human-ratios", default=False, is_flag=True, show_default=True)
@click.option(
    "--genefile",
    # type=click.Path(exists=True, dir_okay=False),
    type=Path_or_Geneset(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.option(
    "--genefile-sheet",
    type=int,
    default=0,
    show_default=True,
    help="Use specified sheet number for genefile with multiple excel sheets",
)
@click.option(
    "--force-plot-genes",
    is_flag=True,
    default=False,
    help="""
              Will force display of all genes in `genefile` as missing values
              """,
)
@click.option(
    "--gene-symbols",
    default=False,
    is_flag=True,
    show_default=True,
    help="Show Gene Symbols on clustermap",
)
@click.option(
    "--genesymbols", default=False, is_flag=True, help="Alias for --gene-symbols"
)
@click.option(
    "--gene-symbol-fontsize",
    default=8,
    show_default=True,
    type=float,
    help="Gene Symbol font size",
)
@click.option(
    "--gene-annot", type=click.Path(exists=True, dir_okay=False), help="annotate"
)
# @click.option("--colorbar-direction", click.Choice(["horizontal", "vertical"]))
@click.option(
    "--highlight-geneids",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=True,
    help="""Optional list of geneids to make into a track
              on the left side of the heatmap
              Should have 1 geneid per line.
              """,
)
@click.option(
    "--volcano-file",
    type=Path_or_Glob(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""tackle volcanoplot output tsv file (...)
              to subselect genes meeting certain thresholds
              """,
)
@click.option(
    "--volcano-filter-params",
    type=(float, float, click.Choice(["pValue", "pAdj"])),
    default=(0, 1, "pAdj"),
    show_default=True,
)
@click.option(
    "--volcano-topn", default=100, show_default=True, help="Top n to plot total"
)
@click.option(
    "--volcano-sortby",
    default="log2_FC",
    type=click.Choice(["log2_FC", "pValue", "pAdj"]),
    show_default=True,
)
@click.option(
    "--volcano-direction",
    type=click.Choice(["up", "down", "both"]),
    default="both",
    show_default=True,
)
@click.option(
    "--cluster-file",
    type=(click.Path(exists=True, dir_okay=False), int),
    default=(None, 0),
    show_default=True,
    multiple=False,
    help="""make clusterplot for a specific subcluster""",
)
@click.option(
    "--linear", default=False, is_flag=True, help="Plot linear values (default log10)"
)
@click.option(
    "--legend-include",
    type=str,
    multiple=True,
    help="""Specific entries in the config file to include in the legend.
              (Default all are included)""",
)
@click.option(
    "--legend-exclude",
    type=str,
    multiple=True,
    help="""Specific entries in the config file to ignore for the legend.
              (Default all are included)""",
)
@click.option(
    "--sample-reference",
    type=str,
    help="metadata key to choose what samples to keep as specified by `--sample-include`",
)
@click.option(
    "--sample-include",
    type=str,
    multiple=True,
    help="samples to include b ased on metadata entry `--sample-reference`",
)
@click.option(
    "--linkage",
    type=click.Choice(
        ["single", "complete", "average", "weighted", "centroid", "median", "ward.D2"]
    ),
    default="ward.D2",
    show_default=True,
    help="linkage method for hierarchical clustering",
)
@click.option("--main-title", default="", show_default=True)
@click.option(
    "--max-autoclusters",
    default=30,
    show_default=True,
    help="""Max number of clusters to try
when `auto` is set for `--nclusters`""",
)
@click.option(
    "--nclusters",
    default=None,
    callback=validate_cluster_number,
    show_default=True,
    help="""If specified by an integer, use that number of clusters via k-means clustering.
             NOT_IMPLEMENTED: If specified as `auto`, will try to find the optimal number of clusters""",
)
@click.option(
    "--cluster-func",
    default="none",
    show_default=True,
    type=click.Choice(["none", "Kmeans", "kmeans", "PAM"]),
    help="""Function used to break data into k distinct clusters,
              k is specified by `n-clusters`
              """,
)
@click.option(
    "--row-annot-side",
    default="left",
    type=click.Choice(["left", "right"]),
    show_default=True,
)
@click.option(
    "--row-cluster/--no-row-cluster",
    default=True,
    is_flag=True,
    show_default=True,
    help="Cluster rows via hierarchical clustering",
)
@click.option(
    "--row-dend-side",
    default="left",
    type=click.Choice(["left", "right"]),
    show_default=True,
)
@click.option(
    "--row-names-side",
    default="right",
    type=click.Choice(["left", "right"]),
    show_default=True,
)
@click.option("--order-by-abundance", default=False, is_flag=True, show_default=True)
@click.option(
    "--seed",
    default="1234",
    help="seed for kmeans clustering",
    callback=validate_seed,
    show_default=True,
)
@click.option(
    "--show-metadata/--hide-metadata",
    default=True,
    show_default=True,
    is_flag=True,
    help="""Show metadata on clustermap if present""",
)
@click.option(
    "--standard-scale",
    type=click.Choice(["None", "0", "1"]),
    default="None",
    show_default=True,
)
@click.option(
    "--show-missing-values/--hide-missing-values",
    default=True,
    is_flag=True,
    show_default=True,
    help="""Whether or not to show missing values on the cluster plot and missing values""",
)
@click.option(
    "--z-score", type=click.Choice(["None", "0", "1"]), default="0", show_default=True
)
@click.option("--z-score-by", type=str, default=None, show_default=True)
@click.pass_context
def cluster2(
    ctx,
    annotate,
    cmap,
    cut_by,
    color_low,
    color_mid,
    color_high,
    col_cluster,
    row_cluster,
    cluster_row_slices,
    cluster_col_slices,
    figsize,
    force_plot_genes,
    genefile,
    genefile_sheet,
    gene_symbols,
    genesymbols,
    gene_symbol_fontsize,
    gene_annot,
    highlight_geneids,
    linear,
    legend_include,
    legend_exclude,
    sample_reference,
    sample_include,
    linkage,
    max_autoclusters,
    nclusters,
    cluster_func,
    main_title,
    order_by_abundance,
    volcano_file,
    volcano_filter_params,
    volcano_topn,
    volcano_direction,
    volcano_sortby,
    cluster_file,
    row_annot_side,
    row_dend_side,
    row_names_side,
    seed,
    show_metadata,
    standard_scale,
    show_missing_values,
    z_score,
    z_score_by,
    add_human_ratios,
):
    outname_kws = dict()

    outname_kws["rds"] = "l" if row_dend_side == "left" else "r"
    outname_kws["cc"] = "T" if col_cluster else "F"
    outname_kws["rc"] = "T" if row_cluster else "F"

    if genesymbols is True and gene_symbols is False:
        gene_symbols = True
    if z_score == "None":
        z_score = None
    if cluster_func == "none":
        cluster_func = None

    # =================================================================

    pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(
        os.path.split(os.path.abspath(__file__))[0], "R", "clusterplot.R"
    )

    r_source(r_file)
    grdevices = importr("grDevices")
    cluster2 = robjects.r["cluster2"]

    # =================================================================

    data_obj = ctx.obj["data_obj"]
    col_meta = data_obj.col_metadata.copy().astype(str).fillna("")
    col_meta = data_obj.col_metadata

    if order_by_abundance and row_cluster:
        raise NotImplementedError("Not Implemented !")

    if not figsize:  # returns empty tuple if not specified
        figsize = None

    missing_values = "masked" if show_missing_values else "unmasked"

    X = data_obj.areas_log
    X.index = X.index.astype(str)
    if show_missing_values:
        _mask = data_obj.mask
        X[_mask] = np.nan
    # filter if sample include is mentioned
    if sample_reference is not None:
        # we take a subset of the data
        col_meta = col_meta.loc[col_meta[sample_reference].isin(sample_include)]
        X = X[col_meta.index]

    genes = None
    if genefile and not isinstance(genefile, dict):
        genes = parse_gid_file(genefile, sheet=genefile_sheet)  # default 0
        outname_kws["genefile"] = fix_name(os.path.splitext(genefile)[0])
    elif genefile and isinstance(genefile, dict):
        _key = list(genefile.keys())[0]
        genes = genefile[_key]
        outname_kws["genefile"] = _key
    if genefile:
        _tokeep = [x for x in genes if x in X.index]  # preserves order
        # X = X.loc[set(X.index) & set(genes)]
        X = X.loc[_tokeep]
    # print(len(X))

    ## filter by volcano output
    _tmp = list()
    _fc, _pval, _ptype = volcano_filter_params
    # for f in volcano_file:
    if volcano_file is not None:
        logger.info(f"Loading volcano file: {volcano_file}")
        volcanofile_basename = os.path.split(volcano_file)[-1]
        name_group = re.search("(?<=group)[_]?(.*)(?=\.tsv)", volcanofile_basename)
        if name_group is None:
            name_group = volcanofile_basename
        else:
            name_group = name_group.group(1)
        logger.info(f"volcanofile name group is {name_group}")
        outname_kws["volcano_file"] = name_group
        outname_kws["direction"] = volcano_direction
        _df = pd.read_table(volcano_file)
        if not np.isfinite(volcano_topn):
            _df = _df[(abs(_df["log2_FC"]) > np.log2(_fc)) & (_df[_ptype] < _pval)]
        else:  # np.isfinite(volcano_topn):
            #
            _n = (
                int(volcano_topn / 2)
                if volcano_direction == "both"
                else int(volcano_topn)
            )
            if volcano_sortby == "log2_FC":
                _up = _df.sort_values(by="log2_FC", ascending=False).head(_n)
                _dn = _df.sort_values(by="log2_FC", ascending=False).tail(_n)
            elif volcano_sortby == "pValue":
                _up = (
                    _df.query("log2_FC>0")
                    .sort_values(by="pValue", ascending=True)
                    .head(_n)
                )
                _dn = (
                    _df.query("log2_FC<0")
                    .sort_values(by="pValue", ascending=True)
                    .head(_n)
                )

            if volcano_direction == "both":
                _df = pd.concat([_up, _dn])
            elif volcano_direction == "up":
                _df = _up
            elif volcano_direction == "down":
                _df = _dn
        _tmp.append(_df)
    if _tmp:
        _dfs = pd.concat(_tmp)
        _genes = _dfs.GeneID.astype(str).unique()
        _tokeep = [x for x in _genes if x in X.index]  # preserve order
        # X = X.loc[set(X.index) & set(genes)]
        X = X.loc[_tokeep]

    if cluster_file[0] is not None:
        if cluster_file[0].endswith("xlsx"):
            _df = pd.read_excel(cluster_file[0])
        else:
            _df = pd.read_table(cluster_file[0])

        _df = _df.rename(
            columns={
                "Cluster": "cluster",
                "Kmeans_Cluster": "cluster",
                "kmeans_cluster": "cluster",
                "kmeans_Cluster": "cluster",
                "Kmeans_cluster": "cluster",
                "PAM_cluster": "cluster",
            }
        )

        _thecluster = cluster_file[1]
        if "cluster" not in _df:
            raise ValueError("improper cluster file, does not have column `cluster`")
        if _thecluster not in _df["cluster"]:
            raise ValueError(f"cluster {_thecluster} not in {cluster_file[0]}")
        _genes = _df[_df["cluster"] == _thecluster]["GeneID"].astype(str)
        _tokeep = [x for x in _genes if x in X.index]
        X = X.loc[_tokeep]
        outname_kws["subcluster"] = _thecluster
        main_title += f"\nCluster {_thecluster}"

    gids_to_annotate = None
    if gene_annot:
        gids_to_annotate = parse_gid_file(gene_annot)

    if linear:
        X = 10**X

    if standard_scale is not None:
        if standard_scale == 1 or standard_scale == "1":
            X = ((X.T / X.max(1)).T).fillna(0)
        elif standard_scale == 0 or standard_scale == "0":
            X = (X / X.max(1)).fillna(0)

    symbols = [data_obj.gid_symbol.get(x, "?") for x in X.index]

    genemapper = GeneMapper()
    symbols = [
        data_obj.gid_symbol.get(x, genemapper.symbol.get(str(x), str(x)))
        for x in X.index
    ]
    X = X.reset_index().assign(GeneSymbol=symbols)

    # X = data_obj.areas_log_shifted.copy()
    # X[data_obj.areas == 0] = 0 # fill the zeros back
    # X[data_obj.mask] = np.NaN

    # data_obj.set_highlight_gids(highlight_geneids)

    row_annot_track = list()
    # TODO parse and look for geneids
    for file_ in highlight_geneids:
        if file_.endswith("xlsx") or file_.endswith("xls"):
            df_ = pd.read_excel(file_, header=None, dtype=str)
        else:
            df_ = pd.read_csv(file_, sep="\t", header=None, dtype=str)
        title_ = get_file_name(file_)
        if len(df_.columns) != 2:
            raise ValueError(
                """File {} should have two columns, first being GeneID,
            second being corresponding value. No headers""".format(
                    file_
                )
            )
        df_.columns = ["GeneID", title_]
        df_ = df_.set_index("GeneID")
        row_annot_track.append(df_)

    row_annot_df = None
    if row_annot_track:
        row_annot_df = pd.concat(row_annot_track, axis=1)
        # make a dataframe that spans all genes about to be plotted
        ixs_ = X.GeneID.astype(str)
        missing_ = set(ixs_) - set(row_annot_df.index)
        intersect_ixs_ = set(ixs_) & set(row_annot_df.index)
        missing_df_ = pd.DataFrame(index=missing_, columns=row_annot_df.columns)
        row_annot_df = pd.concat(
            [row_annot_df.loc[intersect_ixs_], missing_df_], axis=0
        ).fillna("")

    # data_obj.standard_scale    = data_obj.clean_input(standard_scale)
    # data_obj.z_score           = data_obj.clean_input(z_score)

    # col_meta = data_obj.col_metadata.copy().pipe(clean_categorical)
    if add_human_ratios:
        col_meta["HS_ratio"] = data_obj.taxon_ratios["9606"]
    # _expids = ('recno', 'runno', 'searchno', 'label')
    _expids = [
        "recno",
        "runno",
        "searchno",
    ]
    if ("label" in col_meta.columns) and (col_meta["label"].nunique() == 1):
        _expids.append("label")

    # valid_entries = validate_in_config(*legend_include, *legend_exclude, valid_entries=col_meta.index)
    valid_entries = validate_in_config(
        *legend_include, *legend_exclude, valid_entries=col_meta.columns
    )
    legend_include = [
        x for x in legend_include if x in (set(legend_include) & set(valid_entries))
    ]
    legend_exclude = set(legend_exclude) & set(valid_entries)
    to_exclude = set(_expids) | set(legend_exclude)
    if legend_include:
        to_include = legend_include
    else:
        to_include = set(col_meta.columns) - to_exclude
    # col_meta = col_meta.loc[[x for x in col_meta.index if x not in _expids]]
    col_meta = col_meta[[x for x in col_meta.columns if x in to_include]]
    if col_meta.empty:
        col_meta = None
    else:
        col_meta.index.name = "name"
        col_meta = col_meta.reset_index()
        col_meta = utils.set_pandas_datatypes(col_meta)

    ## ============================================================
    # experimental, best to be put in a separate function
    if cut_by is not None:
        _to_none = True
        if ":" in cut_by:
            cut_by = cut_by.split(":")
        else:
            cut_by = [cut_by]
        cut_by = np.array(cut_by)
        for c in cut_by:
            if col_meta is not None and c in col_meta.columns:
                _to_none = False
            if _to_none:
                print("{} not in column metadata".format(cut_by))
                cut_by = robjects.NULL
    if cut_by is None:
        cut_by = robjects.NULL
    elif cut_by is not None:
        outname_kws["cut_by"] = str.join("_", cut_by)

    ## ============================================================
    annot_mat = None
    if annotate:
        # annot_mat = (data_obj.data.loc[ idx[X.index.tolist(), annotate], : ]
        _cols = [x for x in data_obj.data.columns if x not in ("Metric")]
        annot_mat = (
            data_obj.data.loc[
                (data_obj.data.GeneID.isin(X.GeneID.tolist()))
                & (data_obj.data.Metric == annotate)
            ][_cols]
            # .reset_index(level=1, drop=True)
            .set_index("GeneID")
            .fillna(0)
            .astype(int)
        )
        annot_mat.columns.name = annotate
        outname_kws["annotate"] = annotate

        annot_mat["GeneID"] = annot_mat.index
        annot_mat = annot_mat[["GeneID"] + [x for x in annot_mat if x != "GeneID"]]

    # ============================================================
    if z_score_by is not None:
        outname_kws["zscore_by"] = z_score_by

    if linear:
        outname_kws["linear"] = "linear"
    if standard_scale is not None:
        if standard_scale == 1 or standard_scale == "1":
            outname_kws["scale"] = "row"
        if standard_scale == 0 or standard_scale == "0":
            outname_kws["scale"] = "col"

    #
    # data_obj.do_cluster()
    outname_func = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        batch=data_obj.batch_applied,
        batch_method="param" if not data_obj.batch_nonparametric else "nonparam",
        link=linkage,
        missing_values=missing_values,
        normtype=data_obj.normtype,
        outpath=os.path.join(data_obj.outpath, "cluster2"),
    )
    # =================================================================

    min_figwidth = 5
    if figsize is None:  # either None or length 2 tuple
        figheight = 12
        figwidth = max(min(len(X.columns) / 2, 16), min_figwidth)
    else:
        figwidth, figheight = figsize
    if gene_symbols:  # make sure there is enough room for the symbols
        figheight = max(((gene_symbol_fontsize + 2) / 72) * len(X), 12)
        if figheight > 218:  # maximum figheight in inches
            FONTSIZE = max(218 / figheight, 6)
            figheight = 218
    logger.info(f"figsize: {figwidth} x {figheight}")

    # gr_devices = {".png": grdevices.png, ".pdf": grdevices.pdf, ".svg": grdevices.svg}
    gr_devices = {
        ".png": grdevices.png,
        ".pdf": grdevices.cairo_pdf,
        ".svg": grdevices.svg,
    }
    gr_kws = {
        ".png": dict(width=figwidth, height=figheight, units="in", res=300),
        ".pdf": dict(
            width=figwidth,
            height=figheight,
        ),
        ".svg": dict(
            width=figwidth,
            height=figheight,
        ),
    }

    # DONE: IDEA - for all entries in the column and row data that do not have a predefined colormap,
    # assign defaults here (otherwise ComplexHeatmap will assign random colors, which are not always good choices)

    metadata_colorsR = None
    metadata_color_lists = None

    # columns to assign colors:
    # TODO prevent crash if col_meta is empty
    # metacats = col_meta.dtypes[col_meta.dtypes.isin(["string", "object"])].index
    # metacats = col_meta.dtypes[
    #     (col_meta.dtypes == "string") | (col_meta.dtypes == "object")
    # ].index
    # gene_metacats = None
    metacats = col_meta.columns
    if row_annot_df is not None:
        metacats = set(metacats) | set(row_annot_df.columns)

    metadata_color_list = list()
    for metacat in metacats:
        if data_obj.metadata_colors is not None and metacat in data_obj.metadata_colors:
            # metadata_colorsR = robjects.ListVector([])
            # TODO handle when there is not a color for a category
            mapping = data_obj.metadata_colors[metacat]
            entry = robjects.vectors.ListVector(
                {metacat: robjects.vectors.ListVector(mapping)}
            )
            metadata_color_list.append(entry)
            # for value in :  # value is a key-value mapping
            # returns a numpy array with loss of names..
            # entry = robjects.vectors.ListVector({key: robjects.vectors.ListVector(x)})
            # if key != 'response':
            #     continue
            # it works!
            # actually needs to be an "atomic" vector and not a "list vector"
            # this is taken care of in R..
            # Unclear if this can be converted here. Calling R's `unlist` function on the resulting ListVector
        else:
            if metacat not in col_meta:
                continue
            themapping = get_default_color_mapping(
                col_meta[metacat]
            )  # set everything but float
            if themapping is not None:
                entry = robjects.vectors.ListVector(
                    {metacat: robjects.vectors.ListVector(themapping)}
                )
                metadata_color_list.append(entry)

    if metadata_color_list:
        metadata_colorsR = metadata_color_list

    # print(metadata_color_list)
    # ========================================================

    if gids_to_annotate:
        # need to convert to a flat R vector
        gids_to_annotate = robjects.vectors.StrVector(
            [str(x) for x in gids_to_annotate]
        )

    # can we make sure this is always the case everywhere, it makes things much simpler and should
    # avoid bugs
    X["GeneID"] = X.GeneID.astype(str)

    if X.empty:
        print("No data!")
        return

    col_data = robjects.NULL
    if show_metadata and not col_meta is None:
        col_data = col_meta
    #     # cannot convert Categorical column of Integers to Category in py2r
    #     col_data = col_meta.pipe(clean_categorical)  # does this fix the problem?
    # col_data = utils.set_pandas_datatypes(col_data)

    # rpy2 does not map None to robjects.NULL
    if row_annot_df is None:
        row_annot_df = robjects.NULL

    pandas2ri.activate()

    def plot_and_save(
        X,
        out,
        grdevice,
        gr_kws=None,
        annot_mat=None,
        main_title=None,
        file_fmt=".pdf",
        **kws,
    ):
        """
        # gr_kws is redefined here (update fighiehgt /width), no need to take it as input
        """
        if gr_kws is not None:
            logger.info(f"gr_kws: {gr_kws}")
        if gr_kws is None:
            gr_kws = dict()
        # grdevice(file=out, **gr_kw)  # grDevices::png or pdf, etc
        print("Saving", out, "...", end="", flush=True)
        # have to put this here to preserve the layout set by ComplexHeatmap::draw
        if annot_mat is None:
            annot_mat = robjects.NULL
        if main_title is None:
            main_title = robjects.NULL

        figheight = 12
        if gene_symbols:
            figheight = max(((gene_symbol_fontsize + 2) / 72) * len(X), 12)
            # figheight = max(
            #     ((gene_symbol_fontsize + 2) / 72) * (0.7 * (len(X) + 4)), 12
            # )
            if figheight > 218:  # maximum figheight in inches
                figheight = 218
        min_figwidth = 3
        if figsize is None:  # either None or length 2 tuple
            figheight = 12
            figwidth = max(min(len(X.columns) / 2, 16), min_figwidth)
        else:
            figwidth, figheight = figsize

        logger.info(f"figheight: {figheight}, figwidth: {figwidth}")
        gr_kws = {
            ".png": dict(width=figwidth, height=figheight, units="in", res=300),
            ".pdf": dict(
                width=figwidth,
                height=figheight,
            ),
            ".svg": dict(
                width=figwidth,
                height=figheight,
            ),
        }
        gr_kw = gr_kws[file_fmt]

        call_kws = dict(
            data=X,
            color_low=color_low,
            color_mid=color_mid,
            color_high=color_high,
            cut_by=cut_by,
            annot_mat=annot_mat,
            the_annotation=annotate or robjects.NULL,
            z_score=z_score
            if z_score != "None" and z_score is not None
            else robjects.NULL,
            z_score_by=z_score_by or robjects.NULL,
            row_annot_df=row_annot_df,
            col_data=col_data,
            # cmap_name=cmap or np.nan,
            gids_to_annotate=gids_to_annotate or robjects.NULL,
            force_plot_genes=force_plot_genes,
            # genes=genes or robjects.NULL, # this has been filtered above
            show_gene_symbols=gene_symbols,
            standard_scale=True if standard_scale == True else robjects.NULL,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            # metadata=data_obj.config if show_metadata else None,
            nclusters=nclusters or robjects.NULL,
            cluster_func=cluster_func or robjects.NULL,
            max_autoclusters=max_autoclusters,
            show_missing_values=show_missing_values,
            main_title=main_title,
            # mask=data_obj.mask,
            # figsize=figsize,
            linear=linear,  # linear transformation not done in R yet
            normed=data_obj.normed,
            linkage=linkage,
            gene_symbol_fontsize=gene_symbol_fontsize,
            # legend_include=legend_include,
            # legend_exclude=legend_exclude,
            row_dend_side=row_dend_side,
            row_names_side=row_names_side,
            cluster_row_slices=cluster_row_slices,
            cluster_col_slices=cluster_col_slices,
            order_by_abundance=order_by_abundance,
            seed=seed or robjects.NULL,
            metadata_colors=metadata_colorsR or robjects.NULL,
            savedir=os.path.abspath(
                os.path.join(data_obj.outpath, "cluster2")
            ),  # this is for saving things within the r function
        )
        call_kws.update(kws)
        logger.info(f"call_kws: {call_kws}")

        grdevice(file=out, **gr_kw)  # grDevices::png or pdf, etc
        ret = cluster2(**call_kws)
        # print(ret[0])
        grdevices.dev_off()
        print("done.", flush=True)

        sil_df = None
        try:
            sil_df = rpy2.robjects.pandas2ri.rpy2py_dataframe(ret[1])
        except Exception as e:
            pass

        if sil_df is not None:
            row_orders = ret[2]
            the_orders = [
                row_orders.rx2(str(n)) - 1 for n in row_orders.names
            ]  # subtract 1  for zero indexing
            the_orders = [x for y in the_orders for x in y]

            cluster_metrics = sil_df.iloc[the_orders]

            out = outname_func("clustermap_clusters", **outname_kws)
            print("saving", out)
            cluster_metrics.to_csv(out, index=False, sep="\t")

    # ==============================================================================================

    if cluster_func is not None:
        outname_kws[cluster_func] = nclusters
    outname = outname_func("clustermap", **outname_kws)

    for file_fmt in ctx.obj["file_fmts"]:
        grdevice = gr_devices[file_fmt]
        gr_kw = gr_kws[file_fmt]
        out = outname + file_fmt
        annot_mat_to_pass = annot_mat
        if len(X) > 300 and annotate:
            annot_mat_to_pass = None
            logger.info(f"number of genes is {len(X)} >> 300, skipping annotation")
            if "annotate" in outname_kws:
                outname_kws.pop("annotate")

        #################################################################
        ##                          make plot                          ##
        #################################################################

        plot_and_save(
            X,
            out,
            grdevice,
            annot_mat=annot_mat_to_pass,
            main_title=main_title,
            file_fmt=file_fmt,
        )

        ##################################################################
        ##                     plot any annotations                     ##
        ##################################################################

        if data_obj.annotations is None:
            continue

        for annotation in data_obj.annotations:
            annotator = get_annotation_mapper()
            annot_df = annotator.get_annot(annotation)

            subX = X[X.GeneID.isin(annot_df.GeneID)]
            if subX.empty:  # try mapping to homologene
                logger.info(f"Trying to map genes to hs through homologene")
                hgene_mapper = get_hgene_mapper()
                hg_gene_dict = hgene_mapper.map_to_human(X.GeneID)
                hs_genes = set(hg_gene_dict.values())
                _annot_genes = set(annot_df.GeneID) & set(hg_gene_dict.values())
                _tokeep = [k for k, v in hg_gene_dict.items() if v in _annot_genes]
                _tokeep = set(_tokeep)
                logger.info(
                    f"{annotation}: Successfully remapped {len(_tokeep)} genes to hs genes"
                )
                subX = X[X.GeneID.isin(_tokeep)]
            if subX.empty:  # if still empty, continue
                continue

            sub_annot_mat = None
            if annotate and annot_mat is not None:
                sub_annot_mat = annot_mat[annot_mat.GeneID.isin(subX.GeneID)]
            _show_gene_symbols = False
            if len(subX) < 141 and bool(gene_symbols) == False:
                _show_gene_symbols = True
                logger.info(
                    f"number of genes is {len(subX)} << 101, adding symbols. Specify --gene-symbols to override>."
                )
            if gene_symbols is True:  # force override
                _show_gene_symbols = True

            out = (
                outname_func("clustermap", geneset=fix_name(annotation), **outname_kws)
                + file_fmt
            )
            plot_and_save(
                subX,
                out,
                grdevice,
                main_title=annotation,
                annot_mat=sub_annot_mat,
                show_gene_symbols=_show_gene_symbols,
                file_fmt=file_fmt,
            )

        # heatmap = ret.rx('heatmap')
        # print(heatmap)

    # ==============================================================================================

    # ==============================================================================================
    # this doesn't work yet, fix it later
    # discrete_clusters = ret.rx2('discrete_clusters')

    # if discrete_clusters is not robjects.NULL:
    #     the_clusters = discrete_clusters[0]
    #     X['Cluster'] = the_clusters

    # name = outname_func('{}_{}'.format(cluster_func, nclusters)) + '.tsv'
    # X[['GeneID', 'GeneSymbol', 'Cluster']].sort_values(by='Cluster').to_csv(name, sep='\t', index=False)


@main.command("metrics")
@click.option(
    "--full",
    is_flag=True,
    default=False,
    help="Calculate more metrics, requires PSMs data",
)
@click.option(
    "--before-filter/--after-filter",
    default=True,
    is_flag=True,
    show_default=True,
    help="Whether or not to show metrics before or after filtering",
)
@click.option(
    "--before-norm/--after-norm",
    default=True,
    is_flag=True,
    show_default=True,
    help="Whether or not to show area before or after normalization",
)
@click.pass_context
def metrics(ctx, full, before_filter, before_norm):
    data_obj = ctx.obj["data_obj"]
    file_fmts = ctx.obj["file_fmts"]

    make_metrics(
        data_obj,
        file_fmts,
        before_filter=before_filter,
        before_norm=before_norm,
        full=full,
    )


from .overlap import make_overlap


@main.command("overlap")
@click.option(
    "--figsize",
    nargs=2,
    type=float,
    default=(12, 10.5),
    show_default=True,
    help="""Optionally specify the figuresize (width, height) in inches
              """,
)
@click.option(
    "--group",
    type=str,
    default=None,
    help="Metadata entry to group samples for assessing overlap",
)
@click.option(
    "--maxsize", default=15, help="Max number of overlaps to plot", show_default=True
)
@click.option(
    "--non-zeros",
    default=1.0,
    show_default=True,
    type=int_or_ratio(),
    help="""Minimum number of non zeros OR fraction of nonzeros allowed for each sample
              (or sample group. If a decimal is specified (e.g. 1.0), this indicates 100% of values are nonzero.
              If an integer is specified (1), this indicates that 1 value is nonzero.
              """,
)
@click.pass_context
def overlap(ctx, figsize, group, maxsize, non_zeros):
    """
    Plot gene product overlap across experiments
    """
    data_obj = ctx.obj["data_obj"]
    file_fmts = ctx.obj["file_fmts"]

    if group:
        validate_configfile(data_obj.experiment_file, group=group)

    make_overlap(
        data_obj,
        figsize=figsize,
        group=group,
        maxsize=maxsize,
        non_zeros=non_zeros,
        file_fmts=file_fmts,
    )


@main.command("volcano")
@click.option(
    "--bg-marker-color",
    default="#888888",
    help='"rgb(a)" for background marker color',
)
@click.option(
    "--color-down",
    default="blue",
    help='"rgb(a)" for downregulated',
)
@click.option(
    "--color-up",
    default="red",
    help='"rgb(a)" for upregulated',
)
@click.option(
    "-f",
    "--foldchange",
    type=float,
    default=4,
    show_default=True,
    help="fold change cutoff",
)
@click.option(
    "-d",
    "--direction",
    type=click.Choice(("up", "down", "both")),
    default="both",
    show_default=True,
    help="direction of fold change",
    # no suggestion appeared
)
@click.option(
    "-e",
    "--expression-data",
    default=False,
    is_flag=True,
    show_default=True,
    help="Include expression data for each sample in tabular output.",
)
@click.option("--annot-scale", default=1.0, show_default=True)
@click.option(
    "-n",
    "--number",
    type=int,
    default=35,
    show_default=True,
    help="Maximum number of significant genes to highlight (annotate) in plot",
)
@click.option(
    "--number-by",
    type=click.Choice(("abs_log2_FC", "log2_FC", "pValue")),
    default="log2_FC",
    show_default=True,
    help="""How to determine the top n genes to label on plot.
              `log2_Fold_Change` takes top n/2 genes that are up and down""",
)
@click.option(
    "-o",
    "--only-sig",
    default=False,
    is_flag=True,
    show_default=True,
    help="Only export genes that are significantly different (based on set cutoff)",
)
@click.option(
    "-s",
    "--label-scale",
    type=float,
    default=1.5,
    show_default=True,
    help="To what extent to scale the labels",
)
@click.option(
    "--scale",
    type=float,
    help="alias for --label-scale",
    default=-1,
)
@click.option(
    "--marker-scale",
    default=1.4,
    show_default=True,
)
@click.option(
    "--sig",
    type=float,
    default=0.05,
    show_default=True,
    help="Significance cutoff for (B.H. adjusted) pvalue",
)
@click.option(
    "--sig-metric",
    type=click.Choice(("pValue", "pAdj")),
    default="pAdj",
    show_default=True,
    help="Whether to use pValue or B.H. pAdj value for gene highlighting cutoff",
)
@click.option(
    "--p-value/--p-adj",
    default=True,
    is_flag=True,
    show_default=True,
    help="Whether to plot padj or pvalue on volcano plot (does not change underlying data)",
)
@click.option(
    "--highlight-geneids",
    # type=click.Path(exists=True, dir_okay=False),
    type=Path_or_Geneset(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""Optional list of geneids to also highlight. Should have 1 geneid per line. """,
)
@click.option(
    "--force-highlight-geneids",
    is_flag=True,
    default=False,
    show_default=True,
    help="plot all genes specified in `--highlight-geneids` regardless of significance value",
)
@click.option(
    "--formula",
    default=None,
    show_default=True,
    help="""more complex linear regression formula for use with limma.
              Supersedes `group` option""",
)
@click.option(
    "--impute-missing-values / --no-impute-missing-values",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--contrasts",
    default=None,
    show_default=True,
)
@click.option(
    "--genefile",
    # type=click.Path(exists=True, dir_okay=False),
    type=Path_or_Geneset(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.option(
    "--figsize",
    nargs=2,
    type=float,
    default=(5, 5),
    show_default=True,
    help="""Optionally specify the figuresize (width, height) in inches
              """,
)
@click.option(
    "--pch",
    default=16,
    show_default=True,
    help="<experimental> use 21 for circles with borders",
)
@click.option(
    "--alpha",
    default=1.0,
    show_default=True,
    help="<experimental> alpha for top genes that are not highlight_selected_genes",
)
@click.option(
    "--fill-na-zero / --no-fill-na-zero",
    show_default=True,
    default=True,
    is_flag=True,
    help="""
""",
)
@click.pass_context
def volcano(
    ctx,
    annot_scale,
    bg_marker_color,
    color_down,
    color_up,
    # point_size,
    direction,
    foldchange,
    expression_data,
    number,
    number_by,
    only_sig,
    sig,
    sig_metric,
    label_scale,
    marker_scale,
    scale,
    impute_missing_values,
    p_value,
    highlight_geneids,
    force_highlight_geneids,
    formula,
    contrasts,
    genefile,
    figsize,
    pch,
    alpha,
    fill_na_zero,
):
    """
    Draw volcanoplot and highlight significant (FDR corrected pvalue < .05 and > 2 fold change)
    """

    if scale != -1:  # check if scale was changed, then use depreciated option
        warn("--scale is depreciated, please use --marker-scale in the future")
        marker_scale = scale
    from .volcanoplot import volcanoplot

    # yaxis = 'pAdj' if p_adj else 'pValue'
    yaxis = "pValue" if p_value else "pAdj"

    genes = None
    #
    if genefile and isinstance(genefile, Path):
        genes = parse_gid_file(genefile)

    gids_to_highlight = None
    name_for_gids_to_highlight = None
    if highlight_geneids is not None and not isinstance(highlight_geneids, dict):
        gids_to_highlight = parse_gid_file(highlight_geneids)
        name_for_gids_to_highlight = Path(highlight_geneids).name
    elif isinstance(highlight_geneids, dict):
        name_for_gids_to_highlight = list(highlight_geneids.keys())[0]
        gids_to_highlight = highlight_geneids[name_for_gids_to_highlight]

    width, height = figsize

    volcanoplot(
        ctx,
        foldchange,
        expression_data,
        number=number,
        number_by=number_by,
        direction=direction,
        only_sig=only_sig,
        sig=sig,
        sig_metric=sig_metric,
        yaxis=yaxis,
        label_scale=label_scale,
        marker_scale=marker_scale,
        formula=formula,
        contrasts=contrasts,
        genes=genes,
        width=width,
        height=height,
        highlight_geneids=gids_to_highlight,
        impute_missing_values=impute_missing_values,
        bg_marker_color=bg_marker_color,
        annot_scale=annot_scale,
        pch=pch,
        alpha=alpha,
        force_highlight_geneids=force_highlight_geneids,
        fill_na_zero=fill_na_zero,
        extra_outname_info=name_for_gids_to_highlight,
        color_down=color_down,
        color_up=color_up,
    )


@main.command("gsea")
@click.option("--show-result/--no-show-result", default=True, show_default=True)
@click.option("--collapse", type=bool, default=False, show_default=True)
@click.option(
    "--gmt", type=click.Path(exists=True, dir_okay=False), help="Custom GMT geneset"
)
@click.option(
    "--only-human", is_flag=True, help="Only use human genes in a mixed species dataset"
)
@click.option(
    "--use-cluster2",
    is_flag=True,
    help="Use the new cluster2 heatmap routine to plot interesting pathways",
)
@click.option(
    "--geneset",
    type=click.Choice(
        GENESET_MAPPING.keys()
        # (
        #    "hallmark",
        #    "go_biological",
        #    "curated.CP.all",
        #    "curated.CP.KEGG",
        #    "curated.CGP",
        #    "oncogenic.C6",
        #    "curated.CP.BioCarta",
        #    "curated.CP.Reactome",
        #    "curated.CGP",
        #    "curated.CP.PID",
        #    "go.All",
        #    "go.Bio",
        #    "go.Cell",
        #    "go.Molecular",
        #    "motif.gene.sets",
        # )
    ),
    default=("hallmark",),
    show_default=True,
    multiple=True,
)
@click.option(
    "--metric",
    type=click.Choice(
        (
            "Signal2Noise",
            "tTest",
            "Cosine",
            "Euclidean",
            "Manhatten",
            "Pearson",
            "Ratio_of_Classes",
            "Diff_of_Classes",
        )
    ),
    default="Signal2Noise",
    show_default=True,
)
@click.option(
    "--mode",
    type=click.Choice(("Max_probe", "Median_of_probes")),
    default="Max_probe",
    show_default=True,
)
@click.option(
    "--norm",
    type=click.Choice(("meandiv", "None")),
    default="meandiv",
    show_default=True,
)
@click.option("--number-of-permutations", default=1000, show_default=True)
@click.option(
    "--permute",
    type=click.Choice(("phenotype", "gene_set")),
    default="phenotype",
    show_default=True,
)
@click.option(
    "--plot-top-x",
    default=10,
    show_default=True,
    help="number of enrichment plots to generate for each group",
)
@click.option(
    "--rnd-type",
    type=click.Choice(("no_balance", "equalize_and_balance")),
    default="no_balance",
    show_default=True,
)
@click.option("--rnd-seed", default=1234, show_default=True)
@click.option("--seed", help="alias for --rnd-seed", show_default=True, default=None)
@click.option(
    "--scoring-scheme",
    type=click.Choice(("weighted", "classic", "weighted_p2", "weighted_p1.5")),
    default="weighted",
    show_default=True,
)
@click.option("--set-max", default=500, show_default=True)
@click.option("--set-min", default=15, show_default=True)
@click.option(
    "--sort", type=click.Choice(("real", "abs")), default="real", show_default=True
)
@click.option(
    "-n",
    "--number",
    default=9999,
    help="Number of pathways to plot in output",
    show_default=True,
)
@click.option(
    "--group",
    default=None,
    help="""Specify group for GSEA.
                Supersedes main --group option, also allows specifying genes for gene-pathway association""",
)
@click.option(
    "--contrasts",
    default=None,
    help="""""",
)
@click.option(
    "--plot-genes",
    default=False,
    show_default=True,
    is_flag=True,
    help="Plot heatmap of genes in differential gene sets",
)
@click.option(
    "--plot-genes-sig",
    default=False,
    show_default=True,
    is_flag=True,
    help="Plot heatmap of genes in differential gene sets that pass FWER 0.05",
)
@click.option(
    "--annotate",
    type=click.Choice(
        [
            "PSMs",
            "PSMs_u2g",
            "PeptideCount",
            "PeptideCount_S",
            "PeptideCount_S_u2g",
            "PeptideCount_u2g",
            "SRA",
        ]
    ),
    default=None,
    show_default=True,
)
@click.option(
    "--no-homologene-remap",
    is_flag=True,
    default=False,
    help="If not a human species, do not map to human through homologene. Requires pairing with appropriate gene set file",
)
@click.option(
    "--stat-metric",
    type=click.Choice(["FWER p-val", "FDR q-val"]),
    default="FWER p-val",
    show_default=True,
)
@click.pass_context
def gsea(
    ctx,
    show_result,
    collapse,
    gmt,
    only_human,
    geneset,
    metric,
    mode,
    number_of_permutations,
    norm,
    permute,
    plot_top_x,
    rnd_seed,
    seed,
    rnd_type,
    scoring_scheme,
    sort,
    set_max,
    set_min,
    number,
    group,
    contrasts,  # TODO add contrasts  parsing to determine directionality and
    plot_genes,
    plot_genes_sig,
    annotate,
    no_homologene_remap,
    use_cluster2,
    stat_metric,
):
    """
    Run GSEA on specified groups
    """

    gsea_jar = os.path.join(
        os.path.split(os.path.abspath(__file__))[0], "GSEA", "gsea-3.0.jar"
    )

    use_cluster2 = True

    rnd_seed = seed or rnd_seed
    # ===============================================================================================
    stat_metric_cutoff = 0.25 if stat_metric == "FWER p-val" else 0.05
    # ===============================================================================================
    r_source = robjects.r["source"]
    r_file = os.path.join(
        os.path.split(os.path.abspath(__file__))[0], "R", "clusterplot.R"
    )

    r_source(r_file)
    grdevices = importr("grDevices")
    cluster2 = robjects.r["cluster2"]
    # ===============================================================================================
    # TODO move out

    if "--gmt" in sys.argv and "--geneset" not in sys.argv:
        geneset = tuple()

    data_obj = ctx.obj["data_obj"]
    file_fmts = ctx.obj["file_fmts"]

    plt.rc(
        "font",
        **{
            "family": "sans-serif",
            "sans-serif": [
                "DejaVu Sans",
                "Arial",
                "Liberation Sans",
                "Bitstream Vera Sans",
                "sans-serif",
            ],
        },
    )
    if group is None:
        group = data_obj.group  #
        if group is None:
            raise ValueError("Must specify Group")

    # expression = data_obj.areas_log_shifted.copy().fillna(0)
    expression = data_obj.areas_log_shifted.copy().fillna("na")
    expression.index.name = "NAME"
    # if only_human:

    pheno = data_obj.col_metadata.copy()

    # we can groupby one or more groups separated with :
    all_groups = group.split(":")
    logger.info(f"Grouping by {all_groups}")
    gsea_groups = pheno.groupby(all_groups).groups
    for grp, ids in gsea_groups.items():
        pheno.loc[ids, "GSEA_GROUP"] = str.join("", grp)
    # no filenames with dashes allowed for java GSEA
    pheno["GSEA_GROUP"] = pheno["GSEA_GROUP"].str.replace("-", "_")
    _orig_group_value = group
    group = "GSEA_GROUP"

    nsamples = len(pheno.index)

    ngroups = pheno["GSEA_GROUP"].nunique()
    groups = pheno["GSEA_GROUP"].unique()

    # # later
    # pheno_indicator = dict()
    # for ix, grp in enumerate(groups):
    #     pheno_indicator[grp] = ix
    # classes  = list(map(str, [pheno_indicator[grp] for grp in pheno.loc[group]]))

    # java GSEA cannot navigate paths with hyphens
    expression.columns = expression.columns.str.replace("-", "_")
    pheno.index = pheno.index.str.replace("-", "_")

    cls_comparison = ""
    if ngroups == 2:  # reverse it
        cls_comparison = "#{1}_versus_{0}".format(*groups)
    # elif ngroups != 2:
    #     raise ValueError('Must have 2 groups')

    if ngroups < 2:
        raise ValueError("Must have at least 2 groups")

    def maybe_reorder_groups(groups):
        if "cont" in groups[1].lower() and "cont" not in groups[0].lower():
            return groups[::-1]  # switch it
        if "wt" in groups[1].lower() and "wt" not in groups[0].lower():
            return groups[::-1]  # switch it
        return groups

    for groups in itertools.combinations(pheno[group].unique(), 2):
        # # project 731
        # if not (groups[0].startswith("50") or groups[0].startswith("250")):
        #     continue
        # if not (groups[1].startswith("50") or groups[1].startswith("250")):
        #     continue
        # if not groups[0][-3:] == groups[1][-3:]:  # skip different time
        #     continue
        logger.info(groups)
        groups = maybe_reorder_groups(groups)
        logger.info(f"groups after reorder: {groups}")
        # for groups in itertools.combinations(classes, 2):

        # samples0 = classes[groups[0]]
        # samples1 = classes[groups[1]]

        nsamples = pheno[pheno[group].isin(groups)].pipe(len)

        cls_comparison = "#{1}_versus_{0}".format(*groups)

        namegen = partial(
            get_outname,
            name=data_obj.outpath_name,
            taxon=data_obj.taxon,
            non_zeros=data_obj.non_zeros,
            batch=data_obj.batch_applied,
            batch_method="parametric"
            if not data_obj.batch_nonparametric
            else "nonparametric",
            # outpath=os.path.join(data_obj.outpath, "gsea"),
            normtype=data_obj.normtype,
            stat_metric=stat_metric,
            cls=cls_comparison.strip("#"),
        )

        # param_file = os.path.abspath(namegen('gsea_params') + '.txt')

        # # get most recent, sort by name and take last
        # homologene_f = sorted(glob.glob(os.path.join(os.path.split(os.path.abspath(__file__))[0],
        #                                             'data', 'homologene*data')),
        #                     reverse=True)[0]

        # homologene = (pd.read_table(homologene_f, header=None,
        #                             names=('Homologene', 'TaxonID', 'GeneID',
        #                                 'Symbol', 'ProteinGI', 'ProteinAccession'))
        # )
        # # check if we have non-human GeneIDs
        # hgene_query = homologene[ homologene.GeneID.isin(expression.index) ]
        # if hgene_query.TaxonID.nunique() > 1:
        #     raise ValueError('No support for multi-species GSEA')
        # if hgene_query.TaxonID.nunique() == 1 and hgene_query.TaxonID.unique()[0] != 9606:
        #     # remap
        #     print('Remapping {} GeneIDs to human'.format(hgene_query.TaxonID.unique()[0]))
        #     gid_hgene = hgene_query[['GeneID', 'Homologene']].set_index('GeneID')['Homologene'].to_dict()
        #     hgene_hugid = (homologene.query('TaxonID==9606') [['GeneID', 'Homologene']]
        #                 .set_index('Homologene')['GeneID'].to_dict()
        #     )
        #     expression.index = expression.index.map( lambda x: hgene_hugid.get( gid_hgene.get(x) ))

        # _expression = expression.loc[ expression.index.dropna(), pheno[pheno[group].isin(groups)].index
        # ]

        # _expression.index = _expression.index.astype(int)
        # if _expression.index.nunique() < len(_expression.index):
        #     _expression = _expression.groupby(_expression.index).mean()

        # expression = hgene_map(expression)[pheno[pheno[group].isin(groups).index]]

        expression = data_obj.areas_log_shifted.copy()
        expression.index = expression.index.astype(str)
        # java GSEA cannot navigate paths with hyphens
        expression.columns = expression.columns.str.replace("-", "_")
        if not no_homologene_remap:
            expression = hgene_map(expression)
        ix0 = pheno[pheno[group] == (groups[0])].index
        ix1 = pheno[pheno[group] == (groups[1])].index
        ixs = [*ix0, *ix1]  # ixs in pheno are sample names
        expression = expression[ixs]
        # expression = expression.fillna('na')
        expression.index.name = "NAME"

        collapse = "true" if collapse == True else "false"

        if gmt:
            geneset = [*geneset, gmt]

        for gs in geneset:
            # old:
            # outname = namegen("gsea", pathway=os.path.basename(gs))
            # new:
            outname = namegen(
                "gsea",
                pathway=os.path.basename(gs),
                outpath=os.path.join(
                    data_obj.outpath, "gsea", "summaries"
                ),  # outpath for the summary barplot
            )

            # this is outdir for gsea result
            outdir = os.path.abspath(os.path.join(data_obj.outpath, "gsea", "results"))
            report_name = os.path.split(outname)[-1]

            if not os.path.exists(gs):
                f = GENESET_MAPPING.get(gs)
                geneset_file = os.path.join(
                    os.path.split(os.path.abspath(__file__))[0], "GSEA", "genesets", f
                )
            else:  # manually specified geneset file
                geneset_file = gs

            # with NamedTemporaryFile(suffix='.txt') as f, NamedTemporaryFile(mode='w', suffix='.cls') as g:
            # f = namegen("gsea", pathway=os.path.basename(gs)) + ".txt"
            f = (
                namegen(
                    "gsea",
                    outpath=os.path.join(data_obj.outpath, "gsea", "inputs"),
                )
                + ".txt"
            )
            f = Path(f)
            expression.to_csv(f, sep="\t")
            param_file = (
                namegen(
                    "gsea_params",
                    pathway=os.path.basename(gs),
                    outpath=os.path.join(data_obj.outpath, "gsea", "inputs"),
                )
                + ".txt"
            )

            # g_name = namegen("gsea", pathway=os.path.basename(gs)) + ".cls"
            g_name = (
                namegen(
                    "gsea",
                    outpath=os.path.join(data_obj.outpath, "gsea", "inputs"),
                )
                + ".cls"
            )
            with open(g_name, "w") as g:
                # expression.to_csv(f.name, sep="\t")
                # f.close()  # windows compat
                # with open('./results/gsea/kip_dda_kinases_pheno.cls', 'w') as g:

                ## =========================== cls file ================================== #
                # g.write('{} {} 1\n'.format(nsamples, ngroups))

                g.write("{} {} 1\n".format(nsamples, "2"))
                g.write("# {}\n".format(" ".join(groups)))
                # g.write('{}\n'.format(' '.join(classes)))
                # g.write('{}\n'.format(' '.join(pheno[group])))

                _a = pheno[pheno[group] == groups[0]][group].values
                _b = pheno[pheno[group] == groups[1]][group].values
                _out = str.join(" ", [*_a, *_b]) + "\n"
                g.write(
                    _out
                    # "{}\n".format(" ".join(pheno[pheno[group].isin(groups)][group]))
                )

                # g.file.flush()
                # g.close()  # windows compat

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
            rnd_seed\t{rnd_seed}
            gui\tfalse
            """.format(
                collapse=collapse,
                metric=metric,
                mode=mode,
                norm=norm,
                permute=permute,
                rnd_type=rnd_type,
                scoring_scheme=scoring_scheme,
                sort=sort,
                set_max=set_max,
                set_min=set_min,
                nperm=number_of_permutations,
                plot_top_x=plot_top_x,
                cls_comparison=cls_comparison,
                rnd_seed=rnd_seed,
                # rpt_name=report_name,
                res=str(f),
                cls=g.name,  # g is a handle, f is a Path
                gmx=geneset_file,
                outdir=outdir,
                rpt_label=report_name,
            )
            with open(param_file, "w") as param_handle:
                param_handle.write(params)

            try:
                res = subprocess.run(
                    [
                        "java",
                        "-Xmx8192m",
                        "-cp",
                        gsea_jar,
                        "xtools.gsea.Gsea",
                        "-param_file",
                        param_file,
                    ],
                    # stdout=subprocess.PIPE
                )
            except subprocess.CalledProcessError:
                continue

            res.check_returncode()

            folders = [
                os.path.abspath(os.path.join(outdir, f))
                for f in os.listdir(outdir)
                if ("Gsea" in f)
            ]
            folders.sort(key=os.path.getmtime)

            new_folder = folders[-1]
            index = os.path.join(new_folder, "index.html")
            if sys.platform == "darwin":  # check if on OSX
                index = "file://" + index

            if show_result:
                import webbrowser

                webbrowser.open(index)

            # parse result
            # GSEA outputs the summary files of the form:
            # gsea_report_for_[groupname]_[digit_timestamp].xls
            group0 = glob.glob(
                os.path.join(
                    new_folder, "gsea_report_for_{}_[0-9]*.xls".format(groups[0])
                )
            )
            group1 = glob.glob(
                os.path.join(
                    new_folder, "gsea_report_for_{}_[0-9]*.xls".format(groups[1])
                )
            )
            assert len(group0) == len(group1) == 1
            group0_df = pd.read_table(group0[0], index_col="NAME")
            group1_df = pd.read_table(group1[0], index_col="NAME")

            # gsea_sig = group0_df[ group0_df['FWER p-val' ] < .55 ].join(
            #     group1_df[ group1_df['FWER p-val'] < .25 ],
            #     lsuffix='_group0', rsuffix='_group1',
            #     how='outer'
            # )

            cmap = mpl.cm.Reds_r
            bounds = np.linspace(0, 1, 21)
            cnorm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            powernorm = mpl.colors.PowerNorm(0.5, vmin=0, vmax=1)

            gsea_sig = pd.DataFrame()
            cutoff = 0.25
            while len(gsea_sig) < 10 or (
                "NES" in gsea_sig
                and len(gsea_sig) < 5
                and (
                    gsea_sig.query("NES>0").pipe(len) > 1
                    and gsea_sig.query("NES<0").pipe(len) > 1
                )
            ):
                if cutoff > 1:
                    break
                gsea_sig = pd.concat(
                    [
                        group0_df[group0_df[stat_metric] < cutoff],
                        group1_df[group1_df[stat_metric] < cutoff],
                    ]
                )
                cutoff += 0.1

            if gsea_sig.empty:
                print("No gene sets to plot!")
                continue
                # return
            cutoff_val = min(abs(gsea_sig.head(number)["NES"].dropna()))

            tokeep1 = (
                gsea_sig[(gsea_sig["NES"] > 0) & (gsea_sig["NES"] >= cutoff_val)]
                .NES.sort_values(ascending=False)
                .head(number)
                .index
            )
            tokeep2 = (
                gsea_sig[(gsea_sig["NES"] < 0) & (gsea_sig["NES"] <= -cutoff_val)]
                .NES.sort_values(ascending=False)
                .head(number)
                .index
            )
            idx = [x for y in [tokeep1, tokeep2] for x in y]
            gsea_sig = gsea_sig.loc[idx]
            # gesa_sig = gsea_sig.sort_values(by='NES', ascending=False)

            # gsea_sig["color"] = gsea_sig["FWER p-val"].apply(
            gsea_sig["color"] = gsea_sig[stat_metric].apply(
                lambda x: mpl.colors.to_hex(cmap(cnorm(powernorm(x))))
            )
            import textwrap

            gsea_sig.index = [x.replace("HALLMARK_", "") for x in gsea_sig.index]
            gsea_sig.index = gsea_sig.index + [
                "*" if x < stat_metric_cutoff else "" for x in gsea_sig[stat_metric]
            ]
            gsea_sig.index = gsea_sig.index.map(
                lambda x: textwrap.fill(x.replace("_", " "), 24, break_long_words=False)
            )
            # print(gsea_sig[['FWER p-val', 'NES', 'color']])

            # mpl.colors.Normalize(vmin=1.,vmax=1.)

            nes_range = gsea_sig.NES.abs().max() * 2

            # if gsea_sig.NES.max() * gsea_sig.NES.min() < 0:
            #     nes_range = gsea_sig.NES.max() + abs(gsea_sig.NES.min())
            # else:
            #     nes_range = max(abs(gsea_sig.NES.max()), 2.5)

            figwidth = np.round(nes_range * 2.5, decimals=1)
            # edge case when nes_range is zero??
            if np.isnan(figwidth):
                figwidth = 4  # ??
            figheight = 2 + (len(gsea_sig) * 0.8)

            fig = plt.figure(figsize=(figwidth, figheight))
            gs = mpl.gridspec.GridSpec(
                2,
                1,
                width_ratios=[
                    1,
                ],
                height_ratios=[19, 1],
                hspace=0.4,
            )
            # ax0 = plt.subplot(gs[0])
            ax0 = fig.add_subplot(gs[0])
            # ax1 = plt.subplot(gs[1])
            # cax = plt.subplot(gs[1:])
            cax = fig.add_subplot(gs[1:])
            # should always work now
            try:
                gsea_sig["NES"].fillna(0).plot.barh(
                    ax=ax0, color=gsea_sig.color, edgecolor="#222222", linewidth=2
                )
            except Exception as e:
                continue  ## ?? unforseen error
            ax0.axvline(color="#222222")
            # gsea_sig['NES_group1'].fillna(0).plot.barh(colormap=cmap, ax=ax1)
            # ax1.set_yticklabels(())

            gradient = np.apply_along_axis(
                lambda x: cnorm(powernorm(x)), 0, np.linspace(0, 1, 256)
            )
            gradient = np.vstack((gradient, gradient))
            cax.imshow(gradient, aspect="auto", cmap=cmap)
            cax.set_yticklabels(())

            start, end = cax.get_xlim()
            cax.xaxis.set_ticks(np.linspace(start, end, 11))
            cax.set_xticklabels(
                ["{:.2f}".format(x) for x in np.linspace(0, 1, 11)],
                fontdict=dict(size=10),
            )
            # cax.set_xlabel("FWER p-val")
            cax.set_xlabel(stat_metric)
            ax0.grid(axis="x")

            ax0.set_ylabel("")
            ax0.set_xlabel("NES")
            max_nes = gsea_sig.NES.abs().max() + 0.25
            ax0.set_xlim(-max_nes, max_nes)

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
            ax0.text(1, 1.04, groups[1], transform=ax0.transAxes, ha="right")
            ax0.text(
                -0.04,
                -0.12,
                f"* {stat_metric} < {stat_metric_cutoff}",
                fontsize=12,
                transform=ax0.transAxes,
            )
            gs.tight_layout(fig, rect=(0, 0, 1, 0.96))
            # fig.subplots_adjust(left=.4)
            # fig.tight_layout()
            save_multiple(fig, outname, *file_fmts)

            if plot_genes or plot_genes_sig:
                # load genesets
                genesets = dict()
                with open(geneset_file, "r") as f:
                    for line in f:
                        the_geneset = line.strip().split("\t")
                        name = the_geneset[0]
                        # need to do this because GSEA.jar capitalizes all geneset names in output table
                        name = str.upper(name)
                        # genes = the_geneset[2:]

                        thefile = glob.glob(
                            os.path.join(new_folder, "{}*.xls".format(name))
                        )
                        if len(thefile) != 1:
                            continue  # should be OK
                        gsea_res_ = pd.read_table(thefile[0])
                        genesets[name] = (
                            gsea_res_[gsea_res_["CORE ENRICHMENT"] == "Yes"]
                            .PROBE.astype(str)
                            .tolist()
                        )

                # heatmap of individual gene sets

                outpath = os.path.join(
                    # data_obj.outpath, "GSEA_pathways", cls_comparison.strip("#")
                    data_obj.outpath,
                    "gsea",
                    "summaries",
                    cls_comparison.strip("#"),
                )
                outpath = os.path.abspath(outpath)

                missing_values = "masked"
                # outname_func = partial(
                #     get_outname,
                #     name=data_obj.outpath_name,
                #     taxon=data_obj.taxon,
                #     non_zeros=data_obj.non_zeros,
                #     colors_only=data_obj.colors_only,
                #     batch=data_obj.batch_applied,
                #     batch_method="parametric"
                #     if not data_obj.batch_nonparametric
                #     else "nonparametric",
                #     outpath=outpath,
                #     normtype=data_obj.normtype,
                #     missing_values=missing_values,
                #     cls=cls_comparison.strip("#"),
                # )

                if plot_genes_sig:
                    iter_ = gsea_sig[gsea_sig[stat_metric] < 0.25].iterrows()
                else:
                    iter_ = gsea_sig.head(10).iterrows()

                # for ix, row in gsea_sig.iterrows():
                # _expression = data_obj.areas_log_shifted
                for ix, row in iter_:
                    force_plot_genes = False
                    gs = row["GS<br> follow link to MSigDB"]

                    # not this
                    # genes = [maybe_int(x) for x in genesets[gs]]
                    try:
                        genes = genesets[gs]
                    except KeyError:
                        logger.warning(f"Skipping {gs}")
                        continue  # GSEA limits how many results it outputs.
                        # This can happen when there are many differential pathways
                    z_score = "0"
                    row_cluster = True
                    col_cluster = False
                    nclusters = None

                    # col_meta = data_obj.col_metadata.copy()
                    col_meta = pheno
                    # cols in correct order
                    _cols = (
                        col_meta[col_meta[group] == groups[0]].index.tolist()
                        + col_meta[col_meta[group] == groups[1]].index.tolist()
                    )

                    X = expression.replace(to_replace=0, value=np.nan).dropna(
                        axis=0, how="all"
                    )[
                        _cols
                    ]  # drop if all zero
                    _mask = data_obj.mask.copy()
                    _mask.columns = _mask.columns.str.replace("-", "_")
                    if not no_homologene_remap:
                        mask = _mask[expression.columns].pipe(hgene_map, boolean=True)[
                            _cols
                        ]
                    else:
                        mask = _mask[expression.columns][_cols]

                    _expids = ("recno", "runno", "searchno", "label", "rep", "techrep")
                    col_meta = col_meta[
                        [x for x in col_meta.columns if x not in _expids]
                    ]
                    if col_meta.empty:
                        col_meta = robjects.NULL
                    else:
                        col_meta.index.name = "name"
                        col_meta = col_meta.reset_index()
                    main_title = "{}\nNES: {:.2f} FWER pval: {:.2f}".format(
                        gs, row["NES"], row["FWER p-val"]
                    )
                    title_fontsize = 12
                    if len(gs) > 20:
                        title_fontsize -= 1
                    if len(gs) > 30:
                        title_fontsize -= 1
                    if len(gs) > 40:
                        title_fontsize -= 2

                    X.index = X.index.astype(str)
                    X.index.name = "GeneID"
                    genemapper = GeneMapper()
                    symbols = [
                        data_obj.gid_symbol.get(
                            x, genemapper.symbol.get(str(x), str(x))
                        )
                        for x in X.index
                    ]
                    X = X.reset_index().assign(GeneSymbol=symbols)
                    X["GeneID"] = X.GeneID.astype(str)

                    annot_mat = None
                    if annotate:
                        # annot_mat = (data_obj.data.loc[ pd.IndexSlice[X.index.tolist(), annotate], _cols ]
                        # .loc[ pd.IndexSlice[X.index.tolist(), annotate], _cols ]
                        # .reset_index(level=1, drop=True)
                        if not no_homologene_remap:
                            ## TODO cleanup
                            # annot_mat = (data_obj.data.query('Metric == "{}"'.format(annotate))
                            #             .set_index("GeneID")
                            #             .drop('Metric', 1)
                            #             .pipe(hgene_map)
                            #             .round()
                            #             .fillna(0)
                            #             .astype(int)
                            # )

                            _cols = [
                                x for x in data_obj.data.columns if x not in ("Metric")
                            ]
                            annot_mat = (
                                data_obj.data.loc[
                                    (data_obj.data.GeneID.isin(X.GeneID.tolist()))
                                    & (data_obj.data.Metric == annotate)
                                ][_cols]
                                # .reset_index(level=1, drop=True)
                                # .set_index('GeneID')
                                .fillna(0).astype(int)
                            )
                            annot_mat.columns = annot_mat.columns.str.replace("-", "_")
                            # remove columns that aren't in the comparison
                            _toremove = set(annot_mat.columns) - set(X.columns)
                            _keep = [x for x in annot_mat.columns if x not in _toremove]
                            annot_mat = annot_mat[_keep]

                        else:
                            ##TODO implement
                            raise NotImplementedError("")

                            # annot_mat = (data_obj.data.query('Metric == "{}"'.format(annotate))
                            #             .set_index("GeneID")
                            #             .drop('Metric', 1)
                            #             .round()
                            #             .fillna(0)
                            #             .astype(int)
                            # )

                        # annot_mat = (data_obj.data.loc[ data_obj.data.GeneID.isin(X.index) ]
                        #              .query('Metric == "{}"'.format(annotate))[['GeneID', *_cols]]
                        #              .set_index('GeneID')
                        #             .fillna(0)
                        #             .astype(int)
                        # )
                        annot_mat.columns.name = annotate

                    # result = clusterplot(X,
                    if use_cluster2:
                        # rpy2 complient
                        if annot_mat is None:
                            annot_mat = robjects.NULL

                        r_source(r_file)
                        grdevices = importr("grDevices")
                        cluster2 = robjects.r["cluster2"]

                        row_annot_df = None

                        # clean this section up..
                        gene_symbol_fontsize = 7
                        # min_figheight = 6.5
                        # figheight = 12
                        gene_symbols = True
                        # genes is the subselection passed to cluster2
                        if (
                            gene_symbols
                        ):  # make sure there is enough room for the symbols
                            logger.info("gene_symbols is True")
                            figheight = max(
                                ((gene_symbol_fontsize + 2) / 72) * len(genes), 12
                            )
                            logger.info(f"length of matrix is {len(genes)}")
                            logger.info(f"set figheight to {figheight}")
                            if figheight > 218:  # maximum figheight in inches
                                gene_symbol_fontsize = max(218 / figheight, 6)
                                figheight = 218
                        logger.info(f"set figheight to {figheight}")

                        min_figwidth = 4
                        figwidth = (
                            max(min(len(X.columns) / 2, 16), min_figwidth) + 1
                        )  # add an inch for annotations
                        logger.info(f"setting figwidth to {figwidth}")
                        if len(col_meta) > 5:
                            figwidth += 2
                            logger.info(f"setting figwidth to {figwidth}")

                        gr_devices = {
                            ".png": grdevices.png,
                            ".pdf": grdevices.pdf,
                            ".svg": grdevices.svg,
                        }
                        gr_kws = {
                            ".png": dict(
                                width=figwidth, height=figheight, units="in", res=300
                            ),
                            ".pdf": dict(
                                width=figwidth,
                                height=figheight,
                            ),
                            ".svg": dict(
                                width=figwidth,
                                height=figheight,
                            ),
                        }

                        # need to clean up
                        metadata_colorsR = None
                        metadata_color_lists = None
                        metacats = col_meta.dtypes[
                            (col_meta.dtypes == "string")
                            | (col_meta.dtypes == "object")
                        ].index

                        metadata_color_list = list()
                        for metacat in metacats:
                            if (
                                data_obj.metadata_colors is not None
                                and metacat in data_obj.metadata_colors
                            ):
                                # metadata_colorsR = robjects.ListVector([])
                                # for key, x in data_obj.metadata_colors.items():
                                mapping = data_obj.metadata_colors[metacat]
                                entry = robjects.vectors.ListVector(
                                    {metacat: robjects.vectors.ListVector(mapping)}
                                )
                                metadata_color_list.append(entry)
                            else:
                                ncolors = col_meta[metacat].nunique()
                                cmap = sb.color_palette(n_colors=ncolors)
                                color_iter = map(mpl.colors.rgb2hex, cmap)
                                themapping = {
                                    x: c
                                    for x, c in zip(
                                        col_meta[metacat].unique(), color_iter
                                    )
                                }
                                try:
                                    entry = robjects.vectors.ListVector(
                                        {
                                            metacat: robjects.vectors.ListVector(
                                                themapping
                                            )
                                        }
                                    )
                                except Exception as e:  ## ???
                                    logger.error(e)
                                    return
                                metadata_color_list.append(entry)

                        if metadata_color_list:
                            metadata_colorsR = metadata_color_list

                        # outname = outname_func("geneset", pathway=gs)

                        plot_outname = namegen(
                            "gsea",
                            outpath=os.path.join(data_obj.outpath, "gsea", "pathways"),
                        )
                        for file_fmt in ctx.obj["file_fmts"]:
                            grdevice = gr_devices[file_fmt]
                            gr_kw = gr_kws[file_fmt]
                            out = plot_outname + file_fmt

                            grdevice(file=out, **gr_kw)
                            print("Saving", out, "...", end="", flush=True)

                            # heatmap = ret.rx('heatmap')
                            # print(heatmap)

                            # have to put this here to preserve the call to ComplexHeatmap::draw in clusterplot.R

                            result = cluster2(
                                X,
                                main_title=main_title,
                                title_fontsize=title_fontsize,
                                metadata_colors=metadata_colorsR or robjects.NULL,
                                annot_mat=annot_mat,
                                the_annotation=annotate or robjects.NULL,
                                # cmap_name='RdBu_r',
                                # highlight_gids=data_obj.highlight_gids,
                                # highlight_gid_names=data_obj.highlight_gid_names,
                                # force_optimal_ordering=True,
                                # force_plot_genes=force_plot_genes,
                                genes=genes,
                                # gid_symbol=data_obj.gid_symbol,
                                show_gene_symbols=True,
                                z_score=z_score,
                                row_cluster=row_cluster,
                                col_cluster=col_cluster,
                                # metadata=data_obj.config if show_metadata else None,
                                col_data=col_meta,
                                # nclusters=nclusters,
                                show_missing_values=True,
                                # normed=data_obj.normed,
                                gene_symbol_fontsize=gene_symbol_fontsize,
                                seed=1234,
                            )

                            grdevices.dev_off()
                            print("done.", flush=True)

                    else:  # this is not used anymore
                        result = clusterplot(
                            X,
                            main_title=main_title,
                            annot_mat=annot_mat,
                            cmap_name="RdBu_r",
                            highlight_gids=data_obj.highlight_gids,
                            highlight_gid_names=data_obj.highlight_gid_names,
                            force_optimal_ordering=True,
                            force_plot_genes=force_plot_genes,
                            genes=genes,
                            gid_symbol=data_obj.gid_symbol,
                            gene_symbols=True,
                            z_score=z_score,
                            row_cluster=row_cluster,
                            col_cluster=col_cluster,
                            # metadata=data_obj.config if show_metadata else None,
                            col_data=col_meta,
                            nclusters=nclusters,
                            show_missing_values=True,
                            mask=mask,
                            figsize=None,
                            normed=data_obj.normed,
                            gene_symbol_fontsize=12,
                            legend_include=legend_include,
                            legend_exclude=legend_exclude,
                            seed=1234,
                            metadata_colors=data_obj.metadata_colors,
                            circle_col_markers=True,
                            circle_col_marker_size=120,
                            square=False,
                        )

                        g = result["clustermap"]["clustergrid"]

                        outname = outname_func("geneset", pathway=gs)
                        save_multiple(
                            g,
                            outname,
                            *file_fmts,
                        )
                        plt.close(g.fig)


@main.command("ssGSEA")
@click.option(
    "--geneset",
    type=click.Choice(
        (
            "hallmark",
            "go_biological",
            "curated.CP.all",
            "curated.CP.KEGG",
            "curated.CGP",
            "oncogenic.C6",
            "curated.CP.BioCarta",
            "curated.CP.Reactome",
            "curated.CGP",
            "curated.CP.PID",
            "go.All",
            "go.Bio",
            "go.Cell",
            "go.Molecular",
            "motif.gene.sets",
        )
    ),
    default=("hallmark",),
    show_default=True,
    multiple=True,
)
@click.option(
    "-n",
    "--norm",
    type=click.Choice(("rank", "log", "log.rank", "none")),
    default="rank",
    show_default=True,
)
@click.option(
    "--combine_mode",
    type=click.Choice(("combine.off", "combine.replace", "combine.add")),
    default="combine.off",
    show_default=True,
)
@click.option(
    "-w",
    "--weight",
    type=click.FloatRange(0, 100),
    default=0.75,
    show_default=True,
    help="""when weight == 0 all genes have the same weight.
              If weight > 0, the actual values matter and can change the resulting score.
              """,
)
@click.option(
    "-c",
    "--correl",
    type=click.Choice(("rank", "z.score", "symm.rank")),
    default="z.score",
    show_default=True,
)
@click.option(
    "-p",
    "--perm",
    type=int,
    default=1000,
    show_default=True,
)
@click.option(
    "-m",
    "--min-overlap",
    default=10,
    show_default=True,
    help="""minimal overlap between signatures""",
)
@click.option(
    "-r",
    "--rank-plots",
    default=False,
    is_flag=True,
    show_default=True,
    help="Plot rank plots for each gene set",
)
@click.option(
    "-s",
    "--statistic",
    type=click.Choice(("area.under.RES", "Kolmogorov-Smirnov")),
    default="area.under.RES",
    show_default=True,
)
@click.option(
    "-x",
    "--extended-output",
    default=False,
    show_default=True,
    is_flag=True,
    help="If True include additional stats on signature coverage, etc..",
)
@click.option(
    "-g",
    "--globalfdr",
    default=False,
    show_default=True,
    is_flag=True,
    help="Calculate global FDR across all datasets.",
)
@click.option(
    "--z-score/--no-z-score",
    default=True,
    show_default=True,
    is_flag=True,
    help="Calculate global FDR across all datasets.",
)
@click.option(
    "--seed", type=int, default=2009, show_default=True, help="seed for reproducibility"
)
# @click.option('--z-score', default=False, is_flag=True, help='Z-score data before running ssGSEA')
@click.pass_context
def ssGSEA(
    ctx,
    geneset,
    norm,
    combine_mode,
    weight,
    correl,
    perm,
    min_overlap,
    extended_output,
    rank_plots,
    statistic,
    globalfdr,
    seed,
    z_score,
):
    """
    Run ssGSEA with specified geneset
    """
    from .ssGSEA import ssGSEA

    data_obj = ctx.obj["data_obj"]
    file_fmts = ctx.obj["file_fmts"]
    metadata = data_obj.col_metadata

    expression = data_obj.areas_log_shifted
    if z_score:
        expression = expression.apply(sb.matrix.ClusterGrid.z_score, 0)

    namegen = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        outpath=data_obj.outpath,
    )

    # output_prefix = os.path.join(namegen('ssGSEA'), gs)
    for gs in geneset:
        output_dir = os.path.join(data_obj.outpath, "ssGSEA")
        output_prefix = os.path.join(data_obj.outpath, "ssGSEA", gs)
        log_file = os.path.join(output_dir, "ssGSEA_{}.log".format(gs))
        # output_prefix = 'test'
        # log_file = 'test.log'
        if not os.path.exists(output_prefix):
            os.makedirs(output_prefix)

        ssGSEA(
            expression,
            metadata,
            geneset=gs,
            norm=norm,
            combine_mode=combine_mode,
            weight=weight,
            correl=correl,
            nperm=perm,
            min_overlap=min_overlap,
            extended_output=extended_output,
            globalfdr=globalfdr,
            output_prefix=output_prefix,
            log_file=log_file,
            rank_plots=rank_plots,
            seed=seed,
        )

    # default=Path("."),


# @click.option("-p", "--path", help="root path")
# @click.argument("path", type=Path_or_Glob(exists=True, dir_okay=True, file_okay=False))
@main.command("replot_gsea")
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False))
@click.pass_context
def replot_gsea(ctx, path, help=False):
    """ """

    if help is True or (help is False and path is None):
        help_txt = globals()["replot_gsea"].get_help(ctx)
        print(help_txt)
        sys.exit(0)

    from . import plot_gsea_results

    path = Path(path)
    plot_gsea_results.main(path)


@main.command("box")
@click.option(
    "--gene",
    default=None,
    show_default=True,
    multiple=True,
    type=int,
    help="Gene to plot. Multiple allowed",
)
@click.option(
    "--genefile",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.option(
    "--group",
    type=str,
    default=None,
    help="Metadata entry to group color bars. Cannot be used with average",
)
@click.option(
    "--group-order",
    type=str,
    default=None,
    help="""Pipe `|` separated order
to arrange data. For use in conjunction with `--group`
""",
)
@click.option(
    "--retain-order",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Retains original config order""",
)
@click.option(
    "--cmap",
    default=None,
    show_default=True,
    help="""
Any valid, qualitative, colormap? """,
)
@click.option(
    "--linear", default=False, is_flag=True, help="Plot linear values (default log10)"
)
@click.option(
    "--z-score",
    default=False,
    is_flag=True,
    help="Plot zscored values (linear or log10 transformed)",
)
@click.option(
    "--figsize",
    nargs=2,
    type=float,
    default=None,
    show_default=True,
    help="""Optionally specify the figuresize (width, height) in inches
              If not specified, tries to use a reasonable default depending on the number of
              samples.
              """,
)
@click.option("--xtickrotation", default=None, type=int)
@click.option("--xticksize", default=None, type=int)
@click.pass_context
def box(
    ctx,
    group,
    group_order,
    retain_order,
    cmap,
    gene,
    genefile,
    linear,
    z_score,
    figsize,
    xtickrotation,
    xticksize,
):
    if group is None:
        raise ValueError("Must specify group")
    data_obj = ctx.obj["data_obj"]
    metadata = data_obj.col_metadata
    if genefile:
        gene = gene + tuple(parse_gid_file(genefile))
    if len(gene) == 0:
        raise ValueError("Must supply at least 1 gene")

    outpath = os.path.join(data_obj.outpath, "box")
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    outfunc = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        outpath=outpath,
    )
    # if color is not None and average is not None:
    #     raise ValueError('Cannot specify color and average at the same time.')

    data = data_obj.areas_log_shifted.copy()
    data[data_obj.areas == 0] = 0  # fill the zeros back
    data[data_obj.mask] = np.NaN

    try:
        from rpy2.robjects import r
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
    except ModuleNotFoundError:
        print("rpy2 needs to be installed")
        return

    pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], "R", "boxplot.R")

    r_source(r_file)
    grdevices = importr("grDevices")
    Rboxplot = robjects.r["boxplot"]

    plt_size = 6
    gr_devices = {".png": grdevices.png, ".pdf": grdevices.pdf, ".svg": grdevices.svg}
    gr_kws = {
        ".png": dict(width=plt_size, height=plt_size, units="in", res=300),
        ".pdf": dict(
            width=plt_size,
            height=plt_size,
        ),
        ".svg": dict(
            width=plt_size,
            height=plt_size,
        ),
    }

    for g in gene:
        symbol = data_obj.gid_symbol.get(g, "")
        if g not in data.index:
            warn("GeneID {} ({}) not in dataset, skipping..".format(g, symbol))
            continue
        df = data.loc[g].fillna(0).to_frame("Expression").join(metadata).reset_index()

        if linear:
            df = df.apply(lambda x: np.power(10, x))
        elif z_score:
            # TODO: fix this
            df = sb.matrix.ClusterGrid.z_score(df, 0)

        for col in ("runno", "searchno"):
            df[col] = df[col].astype(str)

        title = "{} {}".format(g, symbol)
        outname = outfunc("boxplot_{}_{}".format(symbol, g))

        for file_fmt in ctx.obj["file_fmts"]:
            grdevice = gr_devices[file_fmt]
            gr_kw = gr_kws[file_fmt]

            out = outname + file_fmt

            grdevice(file=out, **gr_kw)
            print("Saving", out, "...", end="", flush=True)

            p = Rboxplot(df, group, title=title)

            grdevices.dev_off()
            print("done.", flush=True)

    # boxplot(data, genes=gene, color=color, cmap=cmap, metadata=col_meta, metadata_colors=data_obj.metadata_colors,
    #         color_order=color_order, linear=linear, z_score=z_score, base_outfunc=outfunc,
    #         file_fmts=ctx.obj['file_fmts'], gid_symbol=data_obj.gid_symbol, figsize=figsize,
    #         xtickrotation=xtickrotation, xticksize=xticksize, retain_order=retain_order
    # )


@main.command("dendrogram")
@click.option("--color", multiple=True)
@click.option("--marker")
@click.option(
    "--linkage",
    type=click.Choice(
        ["single", "complete", "average", "weighted", "centroid", "median", "ward.D2"]
    ),
    default="ward.D2",
    show_default=True,
    help="linkage method for hierarchical clustering",
)
@click.pass_context
def dendrogram(ctx, color, marker, linkage):
    figwidth = 10
    figheight = 10
    try:
        import rpy2.robjects as robjects
    except ModuleNotFoundError:
        print("rpy2 needs to be installed")
        return

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    data_obj = ctx.obj["data_obj"]

    # genes = None
    # if genefile:
    #     genes = parse_gid_file(genefile)

    X = data_obj.areas_log_shifted
    X.index = X.index.astype(str)
    col_meta = data_obj.col_metadata.copy()
    col_meta.index.name = "name"
    col_meta = col_meta.reset_index()

    # ======================================
    outname_func = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        outpath=data_obj.outpath,
    )
    # ======================================

    # genes = None
    # if genefile:
    #     genes = parse_gid_file(genefile)
    #     _tokeep = [x for x in genes if x in X.index]  # preserve order
    #     X = X.loc[_tokeep]

    # now convert to correct orientation
    # variable <color> <shape> gene1 gene2 gene3 ...

    # dfm = df.melt(id_vars=['GeneID', 'GeneSymbol'])
    # dfm = (X.fillna(0).reset_index().melt(id_vars=['GeneID']).merge(col_meta, left_on='variable', right_index=True)
    # )
    # dfm.to_csv(outname_func('pca_input')+'.tsv', sep='\t', index=False)

    for c in color:
        if c in col_meta:
            col_meta[c] = col_meta[c].astype(str)
    if marker in col_meta:
        col_meta[marker] = col_meta[marker].astype(str)

    grdevices = importr("grDevices")
    robjects.pandas2ri.activate()
    r_source = robjects.r["source"]
    r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], "R", "dend.R")

    gr_devices = {".png": grdevices.png, ".pdf": grdevices.pdf, ".svg": grdevices.svg}
    gr_kws = {
        ".png": dict(width=figwidth, height=figheight, units="in", res=300),
        ".pdf": dict(
            width=figwidth,
            height=figheight,
        ),
        ".svg": dict(
            width=figwidth,
            height=figheight,
        ),
    }

    r_source(r_file)
    plotdend = robjects.r["plotdend"]
    outname = outname_func("dend")

    for file_fmt in ctx.obj["file_fmts"]:
        grdevice = gr_devices[file_fmt]
        gr_kw = gr_kws[file_fmt]
        out = outname + file_fmt

        grdevice(file=out, **gr_kw)
        print("Saving", out, "...", end="", flush=True)

        plotdend(
            X.reset_index(),
            col_data=col_meta,
            # color = np.array(color) or robjects.NULL,
            color=np.array(color),
            shape=marker or robjects.NULL,
            linkage=linkage,
        )
        grdevices.dev_off()
        print("done.", flush=True)


@main.command("bar")
@click.option(
    "--gene",
    default=None,
    show_default=True,
    multiple=True,
    type=str,
    help="Gene to plot. Multiple allowed",
)
@click.option(
    "--genefile",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.option(
    "--average",
    is_flag=True,
    default=False,
    show_default=True,
    help="Take average within each subgroup as specified by `--group`",
)
@click.option(
    "--group",
    type=str,
    default=None,
    help="Metadata entry to group color bars. Cannot be used with average",
)
@click.option(
    "--group-order",
    type=str,
    default=None,
    help="""Pipe `|` separated order
to arrange data. For use in conjunction with `--group`
""",
)
@click.option(
    "--retain-order",
    is_flag=True,
    default=True,
    show_default=True,
    help="""Retains original config order. Depreciated (defaults to true)""",
)
@click.option(
    "--cmap",
    default=None,
    show_default=True,
    help="""
Any valid, qualitative, colormap? """,
)
@click.option(
    "--linear/--log",
    default=True,
    is_flag=True,
    help="Toggle between linear and log xform",
)
@click.option(
    "--z-score",
    default=False,
    is_flag=True,
    help="Plot zscored values (linear or log10 transformed)",
)
@click.option(
    "--figsize",
    nargs=2,
    type=float,
    default=None,
    show_default=True,
    help="""Optionally specify the figuresize (width, height) in inches
              If not specified, tries to use a reasonable default depending on the number of
              samples.
              """,
)
@click.option("--xtickrotation", default=None, type=int)
@click.option("--xticksize", default=None, type=int)
@click.pass_context
def bar(
    ctx,
    average,
    group,
    group_order,
    cmap,
    gene,
    genefile,
    linear,
    z_score,
    figsize,
    retain_order,  # depreciated
    xtickrotation,
    xticksize,
):
    log = not linear

    if average == True and group is None:
        raise ValueError("Must specify group with average")

    if group_order is not None:
        group_order = group_order.split("|")

    data_obj = ctx.obj["data_obj"]
    metadata = data_obj.col_metadata

    if genefile:
        gene = gene + tuple(parse_gid_file(genefile))
    if len(gene) == 0:
        raise ValueError("Must supply at least 1 gene")

    outpath = os.path.join(data_obj.outpath, "bar")
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    data = data_obj.areas_log_shifted.copy()
    data[data_obj.areas == 0] = 0  # fill the zeros back
    data[data_obj.mask] = np.NaN
    data.index = data.index.astype(str)

    logger.info(f"batch applied is set to {data_obj.batch_applied}")

    outfunc = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        outpath=outpath,
    )
    # if color is not None and average is not None:
    #     raise ValueError('Cannot specify color and average at the same time.')

    from .rutils import gr_devices, robjects, close_grdevice

    r_source = robjects.r["source"]
    r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], "R", "barplot.R")
    r_source(r_file)
    Rbarplot = robjects.r["barplot"]
    # grdevices = importr("grDevices")
    Rbarplot = robjects.r["barplot"]
    from tackle.containers import GeneMapper

    gm = GeneMapper()

    for g in gene:
        logger.info(f"Plotting {g}")
        if g not in data.index:
            warn("GeneID {} not in dataset, skipping..".format(g))
            continue
        df = data.loc[g].fillna(0).to_frame("Expression").join(metadata).reset_index()
        logger.info(f"{df}")

        if not figsize:
            plt_height = 5
            plt_width = len(df) // 2
            plt_width = max(5, plt_width)
            plt_width = min(24, plt_width)
        else:
            plt_width, plt_height = figsize
        gr_kws = {
            ".png": dict(width=plt_width, height=plt_height, units="in", res=300),
            ".pdf": dict(
                width=plt_width,
                height=plt_height,
            ),
            ".svg": dict(
                width=plt_width,
                height=plt_height,
            ),
        }

        if linear:
            df["Expression"] = df["Expression"].apply(lambda x: np.power(10, x))
            ylab = "Expression (linear)"
        elif log:
            ylab = "log10 Expression"
        #

        if z_score:
            df["Expression"] = sb.matrix.ClusterGrid.z_score(df["Expression"])
            ylab = "z-score " + ylab

        for col in ("runno", "searchno"):
            df[col] = df[col].astype(str)

        symbol = data_obj.gid_symbol.get(
            g,
        )
        if symbol is None:
            symbol = gm.symbol.get(g, "")
        title = "{} {}".format(g, symbol)
        data_xform = (
            "log" if not (linear or z_score) else ("linear" if linear else "zscore")
        )
        outname_str = f"barplot_{g}_{symbol}_{data_xform}"
        if average:
            outname_str += "_average"
        if z_score:
            outname_str += "_zscore"
        outname = outfunc(outname_str)

        for file_fmt in ctx.obj["file_fmts"]:
            grdevice = gr_devices[file_fmt]
            gr_kw = gr_kws[file_fmt]

            out = outname + file_fmt

            grdevice(file=out, **gr_kw)
            print("Saving", out, "...", end="", flush=True)

            df["index"] = pd.Categorical(
                df["index"], df["index"].drop_duplicates(), ordered=True
            )
            p = Rbarplot(
                df,
                average or robjects.NULL,
                group or robjects.NULL,
                group_order=group_order or robjects.NULL,
                title=title,
                ylab=ylab,
            )

            close_grdevice()
            # grdevices.dev_off()
            print("done.", flush=True)


@main.command("genecorr")
@click.option(
    "-g",
    "--gene",
    default=None,
    show_default=True,
    multiple=True,
    type=str,
    help="Gene to plot. Multiple allowed",
)
@click.option(
    "--genefile",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=True,
    multiple=False,
    help="""File of geneids to plot.
              Should have 1 geneid per line. """,
)
@click.option(
    "-m",
    "--method",
    type=click.Choice(["pearson", "spearman"]),
    default="spearman",
    show_default=True,
    help="correlation algorithm to use",
)
@click.option(
    "--fillna/--no-fillna",
    is_flag=True,
    default=True,
    show_default=True,
    help="""whether or not to replace missing observations with zero
              Recommended to keep on with Spearman correlation metric,
              and off when using Pearson correlation metric.
              """,
)
@click.pass_context
def genecorr(
    ctx,
    gene,
    genefile,
    method,
    fillna,
):
    """
    calculate gene-gene expression correlations
    """

    data_obj = ctx.obj["data_obj"]
    metadata = data_obj.col_metadata
    if genefile:
        gene = gene + tuple(parse_gid_file(genefile))
    if len(gene) == 0:
        raise ValueError("Must supply at least 1 gene")

    outpath = os.path.join(data_obj.outpath, "corr")
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    outfunc = partial(
        get_outname,
        name=data_obj.outpath_name,
        taxon=data_obj.taxon,
        non_zeros=data_obj.non_zeros,
        colors_only=data_obj.colors_only,
        batch=data_obj.batch_applied,
        batch_method="parametric"
        if not data_obj.batch_nonparametric
        else "nonparametric",
        normtype=data_obj.normtype,
        outpath=outpath,
    )
    # if color is not None and average is not None:
    #     raise ValueError('Cannot specify color and average at the same time.')

    data = data_obj.areas_log_shifted.copy()
    data[data_obj.areas == 0] = 0  # fill the zeros back
    data[data_obj.mask] = np.NaN
    data.index = data.index.astype(str)

    for g in gene:
        if g not in data.index:
            warn("GeneID {} not in dataset, skipping..".format(g))
            continue

        query = data.loc[g]
        corr_title = f"{method}"
        if fillna:
            data = data.fillna(0)
            corr_title = f"{method}_0fill"

        corr = (
            data.T.corrwith(query, method=method)
            .sort_values(ascending=False)
            .to_frame(corr_title)
        )

        symbol = data_obj.gid_symbol.get(
            g,
        )
        if symbol is None:
            symbol = gm.symbol.get(g, "")

        outname = outfunc(
            "corr_{}_{}_{}{}".format(g, symbol, method, "_0fill" if fillna else "")
        )

        corr.to_csv(outname + ".tsv", sep="\t", index=True)
        print(f"Wrote {outname}.tsv")

        # df = data.loc[g].fillna(0).to_frame("Expression").join(metadata).reset_index()

        # for file_fmt in ctx.obj["file_fmts"]:
        #     grdevice = gr_devices[file_fmt]
        #     gr_kw = gr_kws[file_fmt]

        #     out = outname + file_fmt

        #     grdevice(file=out, **gr_kw)
        #     print("Saving", out, "...", end="", flush=True)

        #     # cannot pass None to r func?
        #     if average is None:
        #         average = np.nan
        #     if group is None:
        #         group = np.nan
        #     if group_order is None:
        #         group_order = np.nan

        #     df["index"] = pd.Categorical(
        #         df["index"], df["index"].drop_duplicates(), ordered=True
        #     )
        #     p = Rbarplot(
        #         df, average, group, group_order=group_order, title=title, ylab=ylab
        #     )

        #     grdevices.dev_off()
        #     print("done.", flush=True)


# @main.command('bar')
# @click.option('--gene',
#               default=None, show_default=True, multiple=True, type=int,
#               help="Gene to plot. Multiple allowed")
# @click.option('--genefile', type=click.Path(exists=True, dir_okay=False),
#               default=None, show_default=True, multiple=False,
#               help="""File of geneids to plot.
#               Should have 1 geneid per line. """)
# @click.option('--average', type=str, default=None, help='Metadata entry to group data and plot average. Cannot be used with color')
# @click.option('--color', type=str, default=None, help='Metadata entry to group color bars. Cannot be used with average')
# @click.option('--color-order', type=str, default=None, help="""Pipe `|` separated order
# to arrange data. For use in conjunction with `--color`
# """)
# @click.option('--retain-order', is_flag=True, default=False, show_default=True,
#               help="""Retains original config order""")
# @click.option('--cmap', default=None, show_default=True, help="""
# Any valid, qualitative, matplotlib colormap. See https://matplotlib.org/examples/color/colormaps_reference.html.
# """)
# @click.option('--linear', default=False, is_flag=True, help='Plot linear values (default log10)')
# @click.option('--z-score', default=False, is_flag=True, help='Plot zscored values (linear or log10 transformed)')
# @click.option('--figsize', nargs=2, type=float, default=None, show_default=True,
#               help='''Optionally specify the figuresize (width, height) in inches
#               If not specified, tries to use a reasonable default depending on the number of
#               samples.
#               ''')
# @click.option('--xtickrotation', default=None, type=int)
# @click.option('--xticksize', default=None, type=int)
# @click.pass_context
# def bar(ctx, average, color, color_order, retain_order, cmap, gene, genefile, linear, z_score,
#         figsize, xtickrotation, xticksize):


#     data_obj = ctx.obj['data_obj']
#     col_meta = data_obj.col_metadata
#     if genefile:
#         gene = gene + tuple(parse_gid_file(genefile))
#     if len(gene) == 0:
#         raise ValueError("Must supply at least 1 gene")

#     outpath = os.path.join(data_obj.outpath, 'bar')
#     if not os.path.exists(outpath):
#         os.makedirs(outpath)

#     outfunc = partial(get_outname, name=data_obj.outpath_name, taxon=data_obj.taxon,
#                       non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
#                       batch=data_obj.batch_applied,
#                       batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
#                       outpath=outpath,
#     )
#     # if color is not None and average is not None:
#     #     raise ValueError('Cannot specify color and average at the same time.')

#     data = data_obj.areas_log_shifted.copy()
#     data[data_obj.areas == 0] = 0 # fill the zeros back
#     data[data_obj.mask] = np.NaN


#     barplot(data, genes=gene, color=color, cmap=cmap, metadata=col_meta, metadata_colors=data_obj.metadata_colors,
#             average=average, color_order=color_order, linear=linear, z_score=z_score, base_outfunc=outfunc,
#             file_fmts=ctx.obj['file_fmts'], gid_symbol=data_obj.gid_symbol, figsize=figsize,
#             xtickrotation=xtickrotation, xticksize=xticksize, retain_order=retain_order
#     )


if __name__ == "__main__":
    main(obj={})
