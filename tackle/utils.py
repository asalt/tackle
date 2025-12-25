"""

"""

import os
import re
import configparser
import glob
import operator as op
import tempfile
import shutil
from pathlib import Path
from collections import OrderedDict, defaultdict, Counter
from functools import lru_cache
from warnings import warn

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



# sb.set_context('notebook', font_scale=1.8)


def clean_categorical(col_data):
    for col in col_data:
        series = col_data[col]
        if isinstance(series.dtype, pd.CategoricalDtype):
            # check to make sure categories are strings
            cats = series.cat.categories
            if not (all(isinstance(x, str) for x in cats)):
                col_data[col] = series.cat.rename_categories(lambda x: str(x))
    return col_data


RESERVED_COLORS = {
    "True": "green",
    "False": "red",
    "NA": "grey",
    "Tbio": "blue",
    "Tchem": "green",
    "Tclin": "pink",
    "Tdark": "black",
}


def get_default_color_mapping(s: pd.Series) -> dict:
    # specials = "label", "plex"
    if s.dtype == "float":
        return
    palette = "bright"
    s_as_numeric = pd.to_numeric(s, errors="coerce")
    has_positive_sum = s_as_numeric.sum() > 0
    limited_unique_values = s_as_numeric.nunique() < 11
    has_no_missing = s_as_numeric.isna().sum() == 0

    if has_positive_sum and limited_unique_values and has_no_missing:  # integer
        s = s_as_numeric
        palette = "light:#4133"

    ncolors = s.nunique()
    cmap = sb.color_palette(palette=palette, n_colors=ncolors)
    color_iter = map(mpl.colors.rgb2hex, cmap)
    themapping = {
        str(x): c for x, c in zip(sorted(s.unique()), color_iter)
    }  # have to make str?
    themapping = {k: v for k, v in themapping.items() if k}

    # update for boolean
    for k in RESERVED_COLORS.keys():
        if k in themapping:
            themapping[k] = RESERVED_COLORS[k]

    # print(s.name, themapping)
    return themapping


STR_DTYPE_COERCION = {
    "TRUE": "True",
    "True": "True",
    "true": "True",
    "FALSE": "False",
    "False": "False",
    "false": "False",
    "NA": "NA",
    "<NA>": "NA",
    "na": "NA",
    "<na>": "NA",
    "nan": "NA",
    "<nan>": "NA",
}


def set_pandas_datatypes(df: pd.DataFrame) -> pd.DataFrame:
    PROTECTED = "plex", "label"

    def decide_dtype(s: pd.Series) -> pd.Series:
        # if s.name == "metDose":
        #     return int(s)
        if isinstance(s, pd.CategoricalDtype):
            cats = s.cat.categories
            if not (all(isinstance(x, str) for x in cats)):
                newcats = [str(x) for x in cats]  # do not need?
                newvalues = [str(x) for x in s]
                return pd.CategoricalDtype(newvalues)
        if s.name in PROTECTED:
            s = s.apply(str)
        elif pd.to_numeric(s, errors="coerce").sum() == 0:
            s = s.apply(str)
            s = s.replace(STR_DTYPE_COERCION)
        else:
            s = pd.to_numeric(s, errors="coerce")
        # if isinstance(s, str):
        #     s = s.replace(old, new)
        return s

    # dfnew = df.copy()
    dfnew = df.apply(decide_dtype)
    # dfnew["aspirinRx"] = dfnew["aspirinRx"].replace(
    #     to_replace={"TRUE": True, "FALSE": False, "NA": np.nan}
    # )

    return dfnew


def fix_name(x):
    return (
        x.replace(":", "_")
        .replace(" ", "_")
        .replace(".", "")
        .replace("/", "dv")
        .replace("+", "")
        .replace("(", "")
        .replace(")", "")
        .replace(")", "")
        .replace("fraction", "frac")
        .replace("treatment", "treat")
        .replace("clustermap", "")
        .replace("normtype", "")
        .replace("zscore_by", "zby")
        .replace("annotate", "annot")
        # .replace("an", "annot")
    )


idx = pd.IndexSlice

N_COLORS = 100
# r_colors = sb.color_palette("coolwarm", n_colors=N_COLORS+1)
r_colors = sb.color_palette("RdBu_r", n_colors=N_COLORS + 1)
STEP = 0.2

from .logger import get_logger as _base_get_logger


def _get_logger(name=__name__):
    return _base_get_logger(name)


logger = _get_logger()


def maybe_int(x):
    try:
        return str(int(x))
    except ValueError:
        # warn('Value {} cannot be converted to int'.format(x))
        return x


def plot_imputed(edata_impute, observed, missing, downshift, scale):
    fig, ax = plt.subplots()
    sb.histplot(edata_impute.stack(), label="All Data", ax=ax, kde=False)
    sb.histplot(observed, label="Observed Values", ax=ax, kde=False)
    sb.histplot(
        edata_impute[missing].stack(), label="Imputed Values", ax=ax, kde=False
    )

    # sb.histplot(
    #     edata_impute[missing].stack(), label="Imputed Values", ax=ax, kde=False
    # )

    ax.legend(fontsize=10, loc="upper right", markerscale=0.4)
    # ax.set_xlim(0, 10)
    title = "downshift : {:.2g} scale : {:.2g}".format(downshift, scale)
    ax.set_title(title)
    # outname = os.path.join('../results/imputation_testing', 'distribution_ds_{:.2g}_scale_{:.2g}'.format(downshift, scale))

    # fig.savefig(outname+'.png', dpi=90)
    # plt.close(fig)


# def impute_missing_old(frame, downshift=2.0, scale=1.0, random_state=1234, make_plot=True):
def impute_missing_old(
    frame, downshift=2.0, scale=1.0, random_state=1234, n_draws=8,
    make_plot=True
):
    """
    frame: is a rectangular expression matrix
    effective_width = scale / sqrt(n_draws)

    """
    # _norm_notna = frame.replace(0, np.NAN).stack()

    observed = frame.replace(0, np.nan).stack().dropna().to_frame()
    #missing = frame.isna()
    missing = observed.isna()

    _norm_notna = observed.stack()
    # _norm_notna += np.abs(_norm_notna.min())
    _mean = _norm_notna.mean()
    _sd = _norm_notna.std()
    _norm = stats.norm(loc=_mean - (_sd * downshift), scale=_sd * scale)
    _number_na = frame.replace(0, np.nan).isna().sum().sum()

    # print(frame.replace(0, np.nan).isna().sum())
    random_value_list = list()
    for i in range(n_draws):
        random_values = _norm.rvs(size=_number_na, random_state=random_state + i)
        random_value_list.append(random_values)

    random_values = np.mean(random_value_list, 0)
    assert random_values.shape == random_value_list[0].shape

    # _areas_log = np.log10(frame.replace(0, np.nan))
    # _areas_log += np.abs(_areas_log.min().min())
    # _areas_log = frame.copy().replace(0, np.nan)
    areas_log = frame.copy()

    start_ix = 0
    for col in areas_log:
        last_ix = areas_log[col].isna().sum()
        # print(_areas_log[col].isna().sum())
        areas_log.loc[areas_log[col].isna(), col] = random_values[
            start_ix : start_ix + last_ix
        ]
        start_ix += last_ix

    if make_plot:
        plot_imputed(areas_log, observed, missing, downshift=downshift, scale=scale)

    return areas_log




def impute_missing_mqish(
    frame,
    downshift: float = 1.8,     # MQ-ish
    effective_width: float = 0.3,
    n_draws: int = 1,
    random_state: int = 1234,
    make_plot: bool = True,
    scale = 0.8
):
    """
    MaxQuant-style MNAR imputation on a rectangular expression matrix (log-intensities).

    effective_width is the final SD as a fraction of the global SD.
    If n_draws > 1, we set scale so that scale / sqrt(n_draws) = effective_width.
    """

    # treat 0 as missing consistently
    frame_na = frame.replace(0, np.nan)

    observed = frame_na.stack().dropna()

    # global mean / sd from observed values only
    _mean = observed.mean()
    _sd = observed.std()

    # choose scale so that the *effective* width after averaging is what we want
    scale = effective_width * np.sqrt(n_draws)

    # downshifted, narrow Gaussian
    norm_dist = stats.norm(loc=_mean - (_sd * downshift), scale=_sd * scale)

    # how many values to impute (Na or 0 originally)
    n_na = frame_na.isna().sum().sum()

    # possibly do multiple draws then average -> same mean, SD / sqrt(n_draws)
    draws = []
    for i in range(n_draws):
        draws.append(norm_dist.rvs(size=n_na, random_state=random_state + i))
    random_values = np.mean(draws, axis=0)

    areas_imputed = frame_na.copy()

    # fill column-by-column
    start_ix = 0
    for col in areas_imputed:
        n_missing_col = areas_imputed[col].isna().sum()
        if n_missing_col == 0:
            continue
        areas_imputed.loc[areas_imputed[col].isna(), col] = random_values[
            start_ix : start_ix + n_missing_col
        ]
        start_ix += n_missing_col

    if make_plot:
        plot_imputed(
            areas_imputed,
            observed,
            frame_na.isna(),  # mask of missing values
            downshift=downshift,
            scale=effective_width,  # report effective width
        )

    return areas_imputed


def impute_missing_lupine(
    frame,
    n_models=None,
    device=None,
    biased=True,
    mode="local",
):
    """
    Impute missing values using the Lupine deep learning model.

    Parameters
    ----------
    frame : pandas.DataFrame
        Rectangular matrix of quantifications with NaN for missing entries.
        Values are expected to be on the same (log) scale as the rest of
        the analysis pipeline (i.e. the same scale used for Gaussian MNAR).
    n_models : int, optional
        Number of ensemble models to fit. If None, uses the value from the
        LUPINE_N_MODELS environment variable, or falls back to 5.
    device : str, optional
        Torch device string, e.g. \"cpu\" or \"cuda\". If None, uses the
        LUPINE_DEVICE environment variable, or \"cpu\".
    biased : bool, optional
        Whether to use Lupine's MNAR-biased batch selection.

    Returns
    -------
    pandas.DataFrame
        Imputed matrix with the same index/columns as the input frame.
    """
    try:
        from lupine.lupine import Lupine
    except ImportError as e:
        raise RuntimeError(
            "Lupine backend selected but the `lupine` package is not importable."
        ) from e

    if mode != "local":
        logger.warning(
            "Lupine joint mode requested but not configured; falling back to local ensemble."
        )

    if n_models is None:
        try:
            n_models = int(os.environ.get("LUPINE_N_MODELS", "5"))
        except ValueError:
            n_models = 5
    if n_models <= 0:
        raise ValueError("n_models must be a positive integer")

    if device is None:
        env_device = os.environ.get("LUPINE_DEVICE")
        if env_device:
            device = env_device
        else:
            # Prefer GPU if available, otherwise fall back to CPU.
            try:
                import torch

                device = "cuda" if torch.cuda.is_available() else "cpu"
            except Exception:
                device = "cpu"

    rows = list(frame.index)
    cols = list(frame.columns)
    mat = frame.values.astype(float)

    gen = np.random.default_rng(seed=18)
    n_layers_hparam_space = [1, 2]
    n_factors_hparam_space = [32, 64, 128, 256]
    n_nodes_hparam_space = [256, 512, 1024, 2048]

    # Create a temporary outpath for Lupine's internal checkpointing.
    out_root = tempfile.mkdtemp(prefix="lupine_", dir=os.getcwd())
    # Ensure trailing separator for the string-concatenation used in LupineBase
    outpath = os.path.join(out_root, "")
    Path(outpath).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(outpath, "tmp")).mkdir(parents=True, exist_ok=True)

    qmats = []
    try:
        for _ in range(n_models):
            n_layers_curr = int(gen.choice(n_layers_hparam_space))
            prot_factors_curr = int(gen.choice(n_factors_hparam_space))
            run_factors_curr = int(gen.choice(n_factors_hparam_space))
            n_nodes_curr = int(gen.choice(n_nodes_hparam_space))
            curr_seed = int(gen.integers(low=1, high=10_000))

            model = Lupine(
                n_prots=mat.shape[0],
                n_runs=mat.shape[1],
                n_prot_factors=prot_factors_curr,
                n_run_factors=run_factors_curr,
                n_layers=n_layers_curr,
                n_nodes=n_nodes_curr,
                rand_seed=curr_seed,
                testing=False,
                biased=biased,
                device=device,
                outpath=outpath,
            )
            model_recon = model.fit_transform(mat)
            qmats.append(model_recon)

        qmats_mean = np.mean(qmats, axis=0)
    finally:
        try:
            shutil.rmtree(out_root)
        except Exception as e:
            logger.warning(f"...{e}")
            # Best-effort cleanup; ignore failures.
            pass

    return pd.DataFrame(qmats_mean, index=rows, columns=cols)


impute_missing = impute_missing_mqish # impute_missing_old
impute_missing = impute_missing_old # impute_missing_old

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
        raise ValueError("Must provide metadata if specifying subgroup")

    if subgroup is None:
        columns = [x for x in df.columns if x not in ("GeneID", "Metric")]
        if isinstance(nonzero_value, float):  # then ratio of total
            threshold = len(columns) * nonzero_value

        # mask = (df.loc[ idx[:, column], :].fillna(0)
        #         .where(lambda x : x != 0)
        #         .count(1)
        #         .where(lambda x: x >= threshold)
        #         .dropna())

        mask = (
            df.loc[df.Metric == column][columns]
            .fillna(0)
            .where(lambda x: x != 0)
            .count(1)
            .where(lambda x: x >= threshold)
            .dropna()
        )

        # gids = mask.index.get_level_values(0)
        gids = df.loc[mask.index].GeneID.values

        # return df.loc[ gids.values ]
        return df[df.GeneID.isin(gids)]

    else:
        all_gids = set()

        # for sample, grp in metadata.T.groupby(subgroup):
        for sample, grp in metadata.groupby(subgroup):
            columns = grp.index

            if isinstance(nonzero_value, float):  # then ratio of total
                threshold = len(columns) * nonzero_value

            mask = (
                df.loc[df.Metric == column][columns]
                .fillna(0)
                .where(lambda x: x != 0)
                .count(1)
                .where(lambda x: x >= threshold)
                .dropna()
            )
            gids = df.loc[mask.index].GeneID

            # mask = (df.loc[ idx[:, column], columns].fillna(0)
            #         .where(lambda x : x != 0)
            #         .count(1)
            #         .where(lambda x: x >= threshold)
            #         .dropna())

            # gids = mask.index.get_level_values(0)
            all_gids |= set(gids)
        return df[df.GeneID.isin(all_gids)]

        # return df.loc[ idx[tuple(all_gids), :], : ]


def filter_sra(df, SRA="S", number_sra=1):
    if SRA is None:
        SRA = "S"

    if SRA == "S":
        sra_list = ("S",)
    elif SRA == "R":
        sra_list = ("S", "R")
    else:
        sra_list = ("S",)

    # mask = ((df.loc[ idx[:, 'SRA'], :].isin(sra_list))
    mask = (
        (df.loc[df.Metric == "SRA"].isin(sra_list))
        .sum(1)
        .where(lambda x: x >= number_sra)
        .dropna()
    )
    # gids = mask.index.get_level_values(0)
    gids = df.loc[mask.index, "GeneID"].values

    # return df.loc[ gids.values ]
    return df.loc[df.GeneID.isin(gids)]


def filter_upept(df, number=1):
    # mask = ((df.loc[ idx[:, 'SRA'], :].isin(sra_list))

    pept_table = (
        df.loc[df.Metric == "PeptideCount_u2g"]
        .drop("Metric", 1)
        .set_index("GeneID")
        .astype(float)
    )
    to_keep = (pept_table > 1).sum(1).where(lambda x: x > 3).dropna().index

    # return df.loc[ gids.values ]
    return df.loc[df.GeneID.isin(to_keep)]


def filter_funcats(df_long, funcats):
    mask = (
        df_long.loc[idx[:, "FunCats"], df_long.columns[0]]
        .str.contains(funcats)
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
    mask = df.loc[df.query('Metric == "TaxonID"').index, df.columns[2]] == taxon
    gids = df.loc[mask.index].GeneID
    # mask = df.loc[ idx[:, ('TaxonID')], : ] == 9606
    # gids = mask[mask].dropna().index.get_level_values(0)
    # return df.loc[ gids.values ]
    return df[df.GeneID.isin(gids)]


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
    if "color" in kwargs:
        color = kwargs.pop("color")
    if "bins" in kwargs:
        kwargs.pop("bins")
    if "edgecolor" in kwargs:
        kwargs.pop("edgecolor")
    X = x[x > 0]
    try:
        nbins = seaborn_bin_calc(X)
    except ZeroDivisionError:
        nbins = 10
    # print(nbins)
    ax.hist(X.values, color="k", bins=nbins, edgecolor="none", **kwargs)

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
    labels = list(reversed(["{:.1f}".format(x) for x in np.arange(1, -1.1, -STEP)]))
    cmap = mpl.colors.ListedColormap(r_colors)
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap)
    cbar.set_ticks(np.arange(0, 2.1, STEP / 2))
    cbar.set_ticklabels(labels)
    ax.set_ylabel("r value")


def make_xaxis(ax, yloc=0, offset=0.05, fmt_str="%1.1f", **props):
    xmin, xmax = ax.get_xlim()
    locs = [loc for loc in ax.xaxis.get_majorticklocs() if loc >= xmin and loc <= xmax]
    (tickline,) = ax.plot(
        locs, [yloc] * len(locs), linestyle="", marker=mpl.lines.TICKDOWN, **props
    )
    (axline,) = ax.plot([xmin, xmax], [yloc, yloc], **props)
    tickline.set_clip_on(False)
    axline.set_clip_on(False)
    for loc in locs:
        ax.text(
            loc,
            yloc - offset,
            fmt_str % loc,
            horizontalalignment="center",
            verticalalignment="top",
        )


def plot_delegator(
    x,
    y,
    stat="pearson",
    filter_zeros=True,
    upper_or_lower="upper",
    colors_only=False,
    shade_correlation=True,
    **kwargs,
):
    if upper_or_lower == "upper":
        func = annotate_stat
    elif upper_or_lower == "lower":
        func = scatter

    # x_nonzero = x[ (~x.isnull()) & ~(x.abs() == np.inf) ].index
    # y_nonzero = y[ (~y.isnull()) & ~(y.abs() == np.inf) ].index
    x_nonzero = x[x > 0].index
    y_nonzero = y[y > 0].index

    nonzeros = list(set(x_nonzero) & set(y_nonzero))
    # nonzeros = list(set(x_nonzero) | set(y_nonzero))

    X, Y = x, y

    if filter_zeros:
        X = x.loc[nonzeros]
        Y = y.loc[nonzeros]

    kwargs["alpha"] = 0.4
    ax = plt.gca()

    if stat == "pearson":
        r = pearson_r(X, Y)
        text = "Pearson"
    elif stat == "spearman":
        r = spearman_r(X, Y)
        text = "Spearman"
    text = "n = {:,}\nr = {:.2f}".format(len(X), r)

    if not shade_correlation:
        pass
    elif np.isnan(r):
        print("Could not calculate r")
    else:
        ax_bg_ix = int(
            round(r + 1, 2) * N_COLORS / 2
        )  # add 1 to shift from -1 - 1 to 0 - 2 for indexing
        ax_bg = r_colors[ax_bg_ix]
        ax.patch.set_facecolor(ax_bg)
        ax.patch.set_alpha(0.5)
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

    size = mpl.rcParams["font.size"]
    ax.annotate(
        text,
        xy=(0.5, 0.5),
        xycoords="axes fraction",
        va="center",
        ha="center",
        size=size,
    )
    # sb.despine(ax=ax, left=True, bottom=True)


def scatter(x, y, ax, xymin, xymax, **kwargs):
    s = 4
    marker = "."
    if "text" in kwargs:
        kwargs.pop("text")
    if "color" in kwargs:
        kwargs.pop("color")
    if "alpha" in kwargs:
        alpha = kwargs.pop("alpha")
    if "s" in kwargs:
        s = kwargs.pop("s")
    if "marker" in kwargs:
        marker = kwargs.pop("marker")

    ax.scatter(x, y, color="#222222", alpha=alpha, s=s, marker=marker, **kwargs)
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
    # rasterized = False
    # if 'rasterized' in save_kwargs:
    #     rasterized = save_kwargs.pop('rasterized')

    for ext in exts:
        out = os.path.abspath(filename + ext)

        if ext == ".pdf":
            rasterized = True
        else:
            rasterized = False

        if verbose:
            print("Saving", out, "...", end="", flush=True)
        fig.savefig(out, dpi=dpi, **save_kwargs)
        if verbose:
            print("done.\n\n", flush=True)


def make_config(path="."):
    config = configparser.ConfigParser()
    config.optionxform = str
    config["Name"] = OrderedDict((("recno", 12345), ("runno", 1), ("searchno", 1)))
    file_ = os.path.join(path, "generic_config.ini")
    with open(file_, "w") as cf:
        print("Creating ", file_)
        config.write(cf)
        cf.write("#runno and searchno are optional, default to 1\n")


@lru_cache()
def read_config(configfile, enforce=True):
    """reads config file and returns the data.
    Needs to have recno at a minimum.
    recno, runno, searchno are used for getting data.
    Other fields are used for PCA plots and (optionally) clustermaps.
    """
    config = configparser.ConfigParser()
    config.optionxform = str
    with open(configfile, "r") as f:
        config.read_file(f)

    sections = config.sections()  # retains order
    FIELDS = ("recno", "runno", "searchno", "label")

    data = defaultdict(
        lambda: dict(runno=1, searchno=1, label="none")
    )  # does not retain order (no guarantee)
    # as of py3.7 I believe it does
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
        if section_key.startswith("__"):
            pass
        elif "recno" not in data[section_key] and enforce:  # record number is required
            print(section_key, "does not have recno defined, skipping")
            data.pop(section_key)

    ordered_data = OrderedDict()
    for key in sections:
        if key not in data.keys():
            continue
        ordered_data[key] = data[key]

    return ordered_data


# genemapper = GeneMapper()


def parse_gid_file(gids, symbol_gid_mapping=None, sheet=0) -> set:
    """
    :gids: collection of files to read and extract geneids
    :symbol_gid_mapping: optional dictionary that maps symbols to geneid
    :sheet: excel sheet to use (when input is excel doc)
    """
    from .containers import get_gene_mapper
    genemapper = get_gene_mapper()

    _df = None
    _dtype = {
        "geneid": str,
        "GeneID": str,
        "genesymbol": str,
        "symbol": str,
        "GeneSymbol": str,
        "Symbol": str,
        "junk": str,
        "GeneID": str,
        "entrez_gene": str,
    }
    if gids.endswith(".csv"):
        _df = pd.read_csv(gids, dtype=_dtype, comment="#")
    elif gids.endswith(".tsv") | gids.endswith(".txt"):
        _df = pd.read_table(gids, dtype=_dtype, comment="#")
    elif gids.endswith(".xlsx"):  # try to parse plain text file
        _df = pd.read_excel(gids, dtype=_dtype, sheet_name=sheet, comment="#")
    else:  # maybe no extension, and is text
        _df = pd.read_table(gids, dtype=_dtype, comment="#")

    _df = _df.rename(columns={x:x.lower() for x in _dtype.keys()})
    #_df = _df.rename(columns={x:"geneid" for x in _dtype.keys()})
    _df.columns = _df.columns.str.lower()


    # If the file has only one column, and its name isn't in the expected set
    if _df.shape[1] == 1 and _df.columns[0].lower() not in [x.lower() for x in _dtype.keys()]:
    # Re-read the file with no header, treat first row as data
        if gids.endswith(".csv"):
            _df = pd.read_csv(gids, dtype=str, header=None, comment="#")
        elif gids.endswith(".tsv") or gids.endswith(".txt"):
            _df = pd.read_table(gids, dtype=str, header=None, comment="#")
        elif gids.endswith(".xlsx"):
            _df = pd.read_excel(gids, dtype=str, sheet_name=sheet, header=None)
        else:
            _df = pd.read_table(gids, dtype=str, header=None, comment="#")

        # Assign column name explicitly
        _df.columns = ["geneid"]

    #

    #
    def get_gid_from_symbol(genesymbol):
        # only works for human

        gid = genemapper.df.query(
            'GeneSymbol == "{}" & TaxonID == "9606"'.format(genesymbol)
        )
        if gid.empty:
            gid = genemapper.df.query(
                'GeneSymbol == "{}" & TaxonID == "10090"'.format(genesymbol)
            )
            if gid.empty:
                warn("Could not find GeneID from genesymbol {}".format(genesymbol))
                return
        # else:
        #     print(genesymbol, gid)
        return gid.index[0]

    #
    #
    if _df is not None and "geneid" in _df:
        return _df.geneid.tolist()
    if _df is not None and "genesymbol" in _df.columns:
        _res = [get_gid_from_symbol(genesymbol) for genesymbol in _df.genesymbol]
        return _res

    # ===================================================================================

    rgx_digit = re.compile(r"\W?(\d+)\W?")
    rgx_word = re.compile(r"([A-Za-z]+\d*)(?=\W)")

    def regex_symbol_xtract(line):
        genesymbol = rgx_word.search(line)
        if genesymbol is None:
            warn("Could not parse GeneID from line {}".format(line))
            return
        # gid = genemapper.df.query('GeneSymbol == "{}" & TaxonID == 9606'.format(line.strip()))
        gid = genemapper.df.query(
            'GeneSymbol == "{}" & TaxonID == "9606"'.format(genesymbol.group())
        )
        if gid.empty:
            warn("Could not parse GeneID from line {}".format(line))
            return
        else:
            return gid.index[0]

    # symbol_gid = {v:k for k, v in genemapper.symbol.items()}

    if symbol_gid_mapping is None:
        symbol_gid_mapping = dict()
    gid_out = list()
    # rgx = re.compile(r'(?<![A-Za-z])(\d+)(?![A-Za-z])')
    # rgx = re.compile(r'(?<=[\w])(\d+)(?![A-Za-z])')
    # rgx_digit = re.compile(r'(?<=\W)(\d+)(?=\W)')
    # rgx_word = re.compile(r"([A-Za-z]+\d*)(?=\W)")
    with open(gids, "r") as iterator:  # try to parse file
        # SPLITCHAR = None
        # if ',' in f:
        #     SPLITCHAR = ','
        # iterator = .split(SPLITCHAR)
        # iterator = f
        if gids.endswith("xlsx") or gids.endswith("xls"):
            _df = pd.read_excel(gids, sheet_name=sheet)
            _valid_cols = [
                x
                for x in [re.search(".*geneid.*", x, flags=re.I) for x in _df.columns]
                if x
            ]
            if _valid_cols and len(_valid_cols) == 1:
                _col = _valid_cols[0].group()
                _df = _df[[_col]]
            elif _valid_cols and len(_valid_cols > 1):
                raise ValueError("Cannot parse {}".format(gids))
            iterator = iter(_df.to_string(index=False).splitlines())
        for line in iterator:
            if line.startswith("#") or not line.strip():
                continue

            linestrip = line.strip()
            try:
                gid = float(linestrip)
                gid = str(int(gid))
                gid_out.append(gid)
            except ValueError:
                # try regex
                try:
                    gid = rgx_digit.search(linestrip).group(1)
                    gid_out.append(str(int(gid)))  ## can make this better
                except AttributeError:
                    # try symbol mapping
                    # TODO: expand from just human
                    genesymbol = rgx_word.search(line)
                    if genesymbol is None:
                        warn("Could not parse GeneID from line {}".format(line))
                        continue
                    # gid = genemapper.df.query('GeneSymbol == "{}" & TaxonID == 9606'.format(line.strip()))
                    gid = genemapper.df.query(
                        'GeneSymbol == "{}" & TaxonID == 9606'.format(
                            genesymbol.group()
                        )
                    )
                    if gid.empty:
                        warn("Could not parse GeneID from line {}".format(line))
                        pass
                    else:
                        gid_out.append(gid.index[0])
                    # warn('Could not parse GeneID from line {}'.format(line))
                    # pass

                    if linestrip.isalnum():
                        gid = regex_symbol_xtract(linestrip)
                        if gid:
                            gid_out.append(gid)
                        continue
            # import ipdb; ipdb.set_trace()
            # try symbol mapping
            # TODO: expand from just human
            # genesymbol = rgx_word.search(line)
            # if genesymbol is None:
            #     warn('Could not parse GeneID from line {}'.format(line))
            #     pass
            # # gid = genemapper.df.query('GeneSymbol == "{}" & TaxonID == 9606'.format(line.strip()))
            # gid = genemapper.df.query('GeneSymbol == "{}" & TaxonID == 9606'.format(genesymbol.group()))
            # if gid.empty:
            #     warn('Could not parse GeneID from line {}'.format(line))
            #     pass
            # else:
            #     gid_out.append(gid.index[0])

    retval = list()
    for gid in gid_out:
        if gid not in retval:
            retval.append(gid)
    # c = Counter()
    # for gid in gid_out:
    #     if c[gid] == 0:
    #         retval.append(gid)
    #         c[gid] += 1

    return set(retval)


def get_file_name(full_file):
    # fname, ext = os.path.splitext(full_file)
    fname, ext = os.path.splitext(os.path.basename(full_file))
    grp = re.search(r"\w+", fname)
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
    if index_col not in df["Metric"].unique():
        return
    # if index_col not in df.index.get_level_values(1).unique():
    #     return
    # selection = df.loc[idx[:, index_col], :]
    selection = df[df.Metric == index_col]
    _cols = [x for x in df.columns if x not in ("GeneID", "Metric")]

    # df.loc[selection.index, _cols] = (
    #     df.loc[selection.index, _cols]
    #     .fillna(method="ffill", axis=1)
    #     .fillna(method="bfill", axis=1)
    # )
    # this should be the new code
    df.loc[selection.index, _cols] = (
        df.loc[selection.index, _cols].ffill(axis=1).bfill(axis=1)
    )

    # df.loc[idx[:, index_col], :] = (selection.fillna(method='ffill', axis=1)
    #                                 .fillna(method='bfill', axis=1)
    # )


def standardize_meta(df):
    # gm = GeneMapper()
    from .containers import get_gene_mapper
    gm = get_gene_mapper()
    # df["GeneSymbol"] = df.index.map(lambda x: gm.symbol.get(str(x), x))
    df["FunCats"] = df.index.map(lambda x: gm.funcat.get(x, ""))
    df["GeneDescription"] = df.index.map(lambda x: gm.description.get(str(x), ""))

    if "TaxonID" not in df:
        logger.warning(f"TaxonID not in df")
        # try to look up by GeneID
        # from mygene import MyGeneInfo

        # mg = MyGeneInfo()

        # taxonids_ = mg.getgenes(df.index.tolist(), fields="taxid", as_dataframe=True)
        tid_lookup = df.GeneID.apply(lambda x: gm.taxon.get(x))
        df["TaxonID"] = tid_lookup

    if "Symbol" in df and "GeneSymbol" in df:
        # logger.warning(f"both Symbol and GeneSymbol are present")
        if (
            df.Symbol.fillna("") != df.GeneSymbol.fillna("")
        ).all():  # then we have a problem and need to rectify
            logger.warning(f"Symbol : {df.Symbol.head()}")
            logger.warning(f"GeneSymbol : {df.GeneSymbol.head()}")
            logger.warning(f"Setting GeneSymbol to Symbol")
        if df.GeneSymbol.isna().sum() > 0 and df.Symbol.isna().sum() == 0:
            df["GeneSymbol"] = df.Symbol
        elif df.Symbol.isna().sum() > 0 and df.GeneSymbol.isna().sum() == 0:
            df["Symbol"] = df.GeneSymbol
        # now replace if that failed
        lookup_dict = df.GeneSymbol.to_dict()
        if df.GeneSymbol.isna().sum() > 0 or df.Symbol.isna().sum() < 0:
            df["GeneSymbol"] = df.index.map(
                lambda x: lookup_dict.get(x, gm.symbol.get(str(x), x))
            )
    # df["FunCats"] = df.index.map(lambda x: gm.funcat.get(x, ""))
    # df["GeneDescription"] = df.index.map(lambda x: gm.description.get(str(x), ""))
    return df


DEFAULT_NAS = [
    "-1.#IND",
    "1.#QNAN",
    "1.#IND",
    "-1.#QNAN",
    "#N/A N/A",
    "#N/A",
    "N/A",
    "n/a",
    "NA",
    "#NA",
    "NULL",
    "null",
    "NaN",
    "-NaN",
    "nan",
    "-nan",
    "",
]


def isna_str(entry):
    return pd.isna(entry) | True if entry in DEFAULT_NAS else False


def parse_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    # expids = ('recno', 'runno', 'searchno')
    expids = tuple()
    metadata_filtered = OrderedDict(
        [(k, v) for k, v in metadata.items() if not k.startswith("__")]
    )
    # col_data = pd.DataFrame.from_dict(metadata, orient='columns').filter(regex='^(?!__)')
    # col_data = pd.DataFrame.from_dict(metadata_filtered, orient='columns')
    col_data = pd.DataFrame(metadata_filtered, columns=metadata_filtered.keys())
    col_data = col_data.loc[[x for x in col_data.index if x not in expids]].T

    for col in col_data.columns:
        if col in ("recno", "runno", "searchno", "label", "plex"):
            col_data[col] = col_data[col].astype(str)
            continue
        # col_data.loc[col_data[col].apply(isna_str), col] = np.NAN
        try:
            col_data[col] = col_data[col].astype(float)
        except ValueError:
            pass
        # try:
        #     col_data[col] = col_data[col].convert_dtypes()
        # except AttributeError:
        #     pass  # for pandas < 1.0

    # do not think this is needed anymore
    # for col in col_data.columns:
    #     if not col_data[col].dtype == np.float:
    #         col_data[col] = col_data[col].fillna('NA')
    return col_data


class iFOT:
    def __init__(self):
        self.file = os.path.join(
            os.path.split(os.path.abspath(__file__))[0], "data", "geneignore.txt"
        )
        self._to_ignore = None

    @property
    def to_ignore(self):
        if self._to_ignore is None:
            self._to_ignore = pd.read_table(
                self.file, comment="#", header=None, names=["GeneID"], dtype=str
            )["GeneID"]
            # a pandas.Series
        return self._to_ignore

    def filter(self, genes):
        return [
            g for g in genes if (self.bool_real(g) and str(g) not in self.to_ignore)
        ]

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


def normalize(
    df,
    name="name",
    ifot=False,
    ifot_ki=False,
    ifot_tf=False,
    median=False,
    quantile75=False,
    quantile90=False,
    genefile_norm=None,
    outcol=None,
    taxon=None,
):
    if ifot:  # normalize by ifot but without keratins
        # nonzero_ix = (df[~df.GeneID.isin(ifot_normalizer.to_ignore)]).index
        nonzero_ix = df[
            (~df.GeneID.isin(ifot_normalizer.to_ignore))
            & (np.isfinite(df.iBAQ_dstrAdj))
        ].index
        norm_ = df.loc[nonzero_ix, "iBAQ_dstrAdj"].sum()
    elif median or quantile75 or quantile90:
        if quantile75:
            q = 0.75
        if quantile90:
            q = 0.90
        if median:
            q = 0.5
        # ifot_normalizer.to_ignore
        nonzero_ix = (
            df[~df.GeneID.isin(ifot_normalizer.to_ignore)]
            .query("iBAQ_dstrAdj > 0")
            .index
        )
        norm_ = df.loc[nonzero_ix, "iBAQ_dstrAdj"].quantile(q=q)

        print(q, norm_)

        if taxon and taxon in UNANNOTATED_TIDS:
            norm_ = 1
    elif ifot_ki:
        norm_ = df.loc[
            df["FunCats"].fillna("").str.contains("KI"), "iBAQ_dstrAdj"
        ].sum()
    elif ifot_tf:
        if taxon and taxon in UNANNOTATED_TIDS:
            norm_ = 1
        else:
            norm_ = df.loc[
                df["FunCats"].fillna("").str.contains("TF"), "iBAQ_dstrAdj"
            ].sum()
    elif genefile_norm:
        gids_for_normalization = parse_gid_file(genefile_norm)
        if not gids_for_normalization:
            warn("No genes found in file: {}".format(genefile_norm))
        overlapping_gids = set(df.index) & set(gids_for_normalization)
        if not overlapping_gids:
            warn("No genes in file {} present in experiment".format(genefile_norm))
        norm_ = df.loc[overlapping_gids, "iBAQ_dstrAdj"].sum()
    else:
        norm_ = 1
    if norm_ == 0:
        # error = '{} has a sum of 0 when trying to normalize, aborting'.format(name)
        error = (
            "{} has a sum of 0 when trying to normalize, skipping normalization".format(
                name
            )
        )
        warn(error)
        # raise click.Abort()
        sum_ = 1
    # print(norm_)
    return df["iBAQ_dstrAdj"] / norm_
    # df[outcol] = df['iBAQ_dstrAdj'] / norm_  # use generic 'area' name for all normalization procedures
    # return df


def genefilter(
    df,
    funcats=None,
    funcats_inverse=None,
    geneid_subset=None,
    ignored_geneid_subset=None,
    fix_histones=True,
):
    if funcats:  # do this after possible normalization
        df = df[df["FunCats"].fillna("").str.contains(funcats, case=False)]
    if funcats_inverse:  # do this after possible normalization
        df = df[~df["FunCats"].fillna("").str.contains(funcats_inverse, case=False)]
    if geneid_subset:  # do this at the end
        df = df.loc[set(df.index) & set(geneid_subset)]
    if ignored_geneid_subset:
        tokeep = list(set(df.index) - set(ignored_geneid_subset))
        df = df.loc[tokeep]
    # valid_ixs = (x for x in df.index if not np.isnan(x))
    # if fix_histones:
    #     tmp = df.query('GeneSymbol.str.contains("HIST")')
    valid_ixs = (x for x in df.index if not pd.isna(x))
    return df.loc[valid_ixs]


def filter_and_assign(
    df,
    name,
    funcats=None,
    funcats_inverse=None,
    geneid_subset=None,
    ignored_geneid_subset=None,
    ifot=False,
    ifot_ki=False,
    ifot_tf=False,
    median=False,
):
    """Filter by funcats and geneid_subset (if given)
    remove NAN GeneIDs"""

    if ifot:  # normalize by ifot but without keratins
        norm_ = df.loc[ifot_normalizer.filter(df.index), "iBAQ_dstrAdj"].sum()
    elif median:
        norm_ = df.loc[ifot_normalizer.filter(df.index), "iBAQ_dstrAdj"].median()
    elif ifot_ki:
        norm_ = df.loc[
            df["FunCats"].fillna("").str.contains("KI"), "iBAQ_dstrAdj"
        ].sum()
    elif ifot_tf:
        norm_ = df.loc[
            df["FunCats"].fillna("").str.contains("TF"), "iBAQ_dstrAdj"
        ].sum()
    else:
        norm_ = 1
    if norm_ == 0:
        error = "{} has a sum of 0 when trying to normalize, aborting".format(name)
        print(error)
        raise click.Abort()
        # sum_ = 1
    df["area"] = (
        df["iBAQ_dstrAdj"] / norm_
    )  # use generic 'area' name for all normalization procedures

    if funcats:  # do this after possible normalization
        df = df[df["FunCats"].fillna("").str.contains(funcats, case=False)]
    if funcats_inverse:  # do this after possible normalization
        df = df[~df["FunCats"].fillna("").str.contains(funcats_inverse, case=False)]
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


def get_outname(
    plottype: str,
    taxon=None,
    tx=None,
    non_zeros=1,
    colors_only=None,
    batch=None,
    batch_method="parametric",
    name=None,  # depreciated
    outpath=".",
    **kwargs,
):
    """
    :colors_only: does nothing, depreciated

    """
    if taxon is None and tx is not None:
        taxon = tx
    print(f"plottype is {plottype}")
    colors = "colors_only" if colors_only else "annotated"
    if "missing_values" in kwargs:
        kwargs.pop("missing_values")

    kwarg_values = list()
    for key, value in filter(str, kwargs.items()):
        _value = str(value).replace(" ", "_").replace("-", "_")
        s = "{}_{}".format(key, _value)
        kwarg_values.append(s)
    kwarg_string = "_".join(kwarg_values) if kwarg_values else ""

    batch_str = ("{}Batch_{}_".format(batch_method, batch) if batch else "").replace(
        "parametricBatch", "paramBatch"
    )
    if (
        bool(plottype) and os.path.split(outpath)[-1] != plottype
    ):  # or dtype, doesn't have to be a plot
        outpath = os.path.join(outpath, plottype)

    if batch_str:
        outpath = os.path.join(outpath, batch_str)
    if taxon != "all":
        outpath = os.path.join(outpath, taxon)

    # "{}".format(kwargs)
    fname = "{}more_{}".format(
        # name,
        # plottype,
        non_zeros,
        kwarg_string,
    ).strip("_")
    fname = (
        fname.replace("noCov_", "")
        .replace("norm_", "")
        .replace("sort_", "")
        .replace("Treatment", "treat")
        .replace("median", "med")
        .replace("__", "_")
    )
    fname = fname.lstrip(r".")

    # here is where we could split the path further
    # this might not be necessary here anymore
    # unsure how often this is hit
    if os.name == "nt" and len(fname) > 260:
        fname = fname[:260]
    elif len(fname) > 333:
        fname = fname[:333]
    # import ipdb; ipdb.set_trace()
    ret = os.path.join(outpath, fname)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    return ret


class TooManyCategories(Exception):
    pass


from collections import OrderedDict

try:
    from collections import Callable
except ImportError:
    from collections.abc import Callable


class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if default_factory is not None and not isinstance(default_factory, Callable):
            raise TypeError("first argument must be callable")
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
            args = (self.default_factory,)
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy

        return type(self)(self.default_factory, copy.deepcopy(self.items()))

    def __repr__(self):
        return "OrderedDefaultDict(%s, %s)" % (
            self.default_factory,
            OrderedDict.__repr__(self),
        )


def hgene_map(expression, keep_orig=False, boolean=False):
    # get most recent, sort by name and take last
    homologene_f = sorted(
        glob.glob(
            os.path.join(
                os.path.split(os.path.abspath(__file__))[0], "data", "homologene*data"
            )
        ),
        reverse=True,
    )[0]

    homologene = pd.read_table(
        homologene_f,
        header=None,
        dtype=str,
        names=(
            "Homologene",
            "TaxonID",
            "GeneID",
            "Symbol",
            "ProteinGI",
            "ProteinAccession",
        ),
    )
    # check if we have non-human GeneIDs
    hgene_query = homologene[homologene.GeneID.isin(expression.index)]
    if hgene_query.TaxonID.nunique() > 1:
        raise ValueError("No support for multi-species GSEA")
    if hgene_query.TaxonID.nunique() == 1 and hgene_query.TaxonID.unique()[0] != "9606":
        # remap
        print("Remapping {} GeneIDs to human".format(hgene_query.TaxonID.unique()[0]))
        gid_hgene = (
            hgene_query[["GeneID", "Homologene"]]
            .set_index("GeneID")["Homologene"]
            .to_dict()
        )
        hgene_hugid = (
            homologene.query('TaxonID=="9606"')[["GeneID", "Homologene"]]
            .set_index("Homologene")["GeneID"]
            .to_dict()
        )
        orig = expression.index
        expression.index = expression.index.map(
            lambda x: hgene_hugid.get(gid_hgene.get(x))
        )
        if keep_orig == True:
            expression["GeneID_orig"] = orig
        else:
            expression = expression.loc[[x for x in expression.index if x]]
    else:
        return expression
    # _expression = expression.loc[ expression.index.dropna(), pheno[pheno[group].isin(groups)].index]
    _expression = expression.loc[expression.index.dropna()]

    _expression.index = _expression.index.astype(int)
    if _expression.index.nunique() < len(_expression.index):
        # take mean of nonzero values
        _expression = _expression.replace(0, np.nan).groupby(_expression.index).mean()
    if boolean:
        _expression = _expression.applymap(bool)
    return _expression


from tempfile import NamedTemporaryFile
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



# from .containers import GeneMapper


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
