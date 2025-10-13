"""
Helpers for running Limma differential expression workflows via rpy2.

The core entry point is :func:`run_limma_pipeline`, which takes preprocessed
expression and phenotype tables and returns a dictionary of pandas DataFrames
indexed by contrast name.  This isolates the heavy rpy2/R interop from the
higher-level container orchestration and keeps the logic testable.
"""

from __future__ import annotations

import itertools
import logging
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Any
import re

import numpy as np
import pandas as pd


_SANITIZE_REPLACEMENTS = (
    (":", "_"),
    (" ", "_"),
    ("-", "_"),
    ("+", "_"),
    ("/", "_"),
    ("?", "qmk"),
)


def _sanitize_name(value: str) -> str:
    """
    Mirror the ad-hoc name munging the CLI previously performed before handing
    identifiers to limma/model.matrix.  Keeping the behaviour identical avoids
    churn in the generated contrast labels and joins.
    """
    for src, dest in _SANITIZE_REPLACEMENTS:
        value = value.replace(src, dest)
    return value


def _sanitize_names(values: Iterable[str]) -> List[str]:
    return [_sanitize_name(v) for v in values]


def normalize_formula_targets(
    formula: Optional[str],
    gene_index: Sequence[Any],
    symbol_lookup: Mapping[str, Sequence[str]],
    logger: Optional[logging.Logger] = None,
) -> Tuple[Optional[str], Optional[List[Any]]]:
    """
    Interpret optional left-hand-side protein references in a limma formula.

    When users pass strings like ``HER2 ~ 0 + model`` we treat ``HER2`` as a
    request to restrict the analysis to that protein (looked up either by
    GeneID or symbol) while the design continues to come from the right-hand
    side.  The function returns the RHS-only formula along with the resolved
    GeneIDs that should be retained in the expression matrix.
    """
    if not formula or "~" not in formula:
        return formula, None

    lhs_raw, rhs_raw = formula.split("~", 1)
    lhs = lhs_raw.strip()
    rhs = rhs_raw.strip()

    if not lhs:
        normalised = rhs if rhs.startswith("~") or not rhs else "~" + rhs
        return (normalised if normalised else None), None

    tokens = [tok.strip() for tok in lhs.split("+") if tok.strip()]
    if not tokens:
        raise ValueError("Formula left-hand side must specify at least one protein.")

    index_map: Dict[str, Any] = {}
    for gid in gene_index:
        gid_str = str(gid)
        index_map.setdefault(gid_str, gid)
        if isinstance(gid, str):
            index_map.setdefault(gid, gid)

    gene_ids: List[Any] = []

    def _lookup_symbol_matches(token: str) -> List[str]:
        if not symbol_lookup:
            return []
        search_tokens = {
            token,
            _sanitize_name(token),
            token.casefold(),
            _sanitize_name(token).casefold(),
        }
        matches: List[str] = []
        for sym, ids in symbol_lookup.items():
            if sym is None:
                continue
            sym_str = str(sym)
            sym_tokens = {
                sym_str,
                _sanitize_name(sym_str),
                sym_str.casefold(),
                _sanitize_name(sym_str).casefold(),
            }
            if sym_tokens & search_tokens:
                matches.extend(ids)
        return matches

    for token in tokens:
        token_str = str(token)
        if token_str in index_map:
            gene_ids.append(index_map[token_str])
            continue
        matches = list(dict.fromkeys(_lookup_symbol_matches(token_str)))
        if not matches:
            raise ValueError(f"Unknown protein '{token}' in formula left-hand side.")
        if len(matches) > 1:
            raise ValueError(
                f"Protein '{token}' maps to multiple GeneIDs: {list(matches)}; "
                "please disambiguate."
            )
        gene_ids.append(matches[0])

    if not rhs:
        raise ValueError(
            "Formula must include a right-hand side when specifying proteins on the left."
        )

    normalised = rhs if rhs.startswith("~") else "~" + rhs
    if logger:
        logger.info("Limma restricting to GeneIDs via formula LHS: %s", gene_ids)
    # Deduplicate while preserving order
    seen = set()
    ordered_ids = []
    for gid in gene_ids:
        if gid in seen:
            continue
        seen.add(gid)
        ordered_ids.append(gid)
    return normalised, ordered_ids


def _resolve_contrasts(
    design_terms: Sequence[str], contrasts_spec: Optional[str]
) -> List[str]:
    """
    Convert either a CSV string of contrasts or an implicit pairwise expansion
    into the list consumed by limma::makeContrasts, matching the historical
    behaviour in containers.stat_model.
    """
    if contrasts_spec:
        if isinstance(contrasts_spec, str):
            raw = [part.strip() for part in contrasts_spec.split(",")]
        else:
            raw = list(contrasts_spec)
        return [c for c in raw if c and not c.startswith("#")]

    candidates = [term for term in design_terms if "Intercept" not in term]
    return ["{} - {}".format(a, b) for a, b in itertools.combinations(candidates, 2)]


def _inject_gene_covariates_in_formula(
    formula: Optional[str],
    edata: pd.DataFrame,
    pheno: pd.DataFrame,
    symbol_lookup: Optional[Mapping[str, Sequence[str]]],
    logger: Optional[logging.Logger] = None,
) -> Tuple[Optional[str], pd.DataFrame]:
    """Create pheno columns from gene IDs/symbols referenced on the RHS.

    Rewrites the formula to use safe variable names. Returns (new_formula, new_pheno).
    """
    if not formula:
        return formula, pheno

    # Tokenize names: numeric or identifier words
    tokens = set(re.findall(r"`[^`]+`|[A-Za-z_.][A-Za-z0-9_.]*|[0-9]+", formula))
    # Skip obvious operators/reserved
    skip = {"I", "C"}
    tokens = [t for t in tokens if t not in skip]

    # Build index lookup supporting str/int ids
    edata_index_map: Dict[str, Any] = {}
    for gid in edata.index:
        s = str(gid)
        edata_index_map.setdefault(s, gid)
        if isinstance(gid, str):
            edata_index_map.setdefault(gid, gid)

    rewritten = formula
    new_pheno = pheno.copy()

    def replace_token(fml: str, old: str, new: str) -> str:
        # Replace whole-token occurrences with safe var
        pattern = re.compile(rf"(?<![A-Za-z0-9_`]){re.escape(old)}(?![A-Za-z0-9_`])")
        return pattern.sub(new, fml)

    # Precompute sanitized sample indices for robust alignment
    pheno_idx_sanitized = [_sanitize_name(str(x)) for x in new_pheno.index]
    replaced: Dict[str, str] = {}
    for tok in sorted(tokens, key=len, reverse=True):
        # Skip if token already a pheno column
        if tok in new_pheno.columns:
            continue

        gene_id_obj: Optional[Any] = None
        # Numeric token interpreted as GeneID
        if tok.isdigit() and tok in edata_index_map:
            gene_id_obj = edata_index_map[tok]
        # Symbol lookup (case-insensitive + sanitized)
        elif symbol_lookup:
            cand_ids = []
            search = {tok, _sanitize_name(tok), tok.casefold(), _sanitize_name(tok).casefold()}
            for sym, ids in symbol_lookup.items():
                if sym is None:
                    continue
                s = str(sym)
                sym_set = {s, _sanitize_name(s), s.casefold(), _sanitize_name(s).casefold()}
                if sym_set & search:
                    cand_ids.extend(ids)
            cand_ids = list(dict.fromkeys(cand_ids))
            if len(cand_ids) == 1 and cand_ids[0] in edata_index_map:
                gene_id_obj = edata_index_map[cand_ids[0]]
            elif len(cand_ids) > 1:
                raise ValueError(
                    f"Ambiguous symbol '{tok}' in formula RHS maps to multiple GeneIDs: {cand_ids}"
                )

        if gene_id_obj is None:
            continue  # Not a gene ref; leave as-is

        safe_var = _sanitize_name(f"GID_{gene_id_obj}")
        # Inject covariate values from edata row for each sample
        values = edata.loc[gene_id_obj]
        # Align to pheno rows (samples) using sanitized names on both sides
        values.index = [_sanitize_name(str(x)) for x in values.index]
        aligned = values.reindex(index=pheno_idx_sanitized)
        new_pheno[safe_var] = pd.to_numeric(aligned.values, errors="coerce")
        rewritten = replace_token(rewritten, tok, safe_var)
        replaced[tok] = safe_var
        if logger:
            logger.info("Injected covariate %s from GeneID %s for formula", safe_var, gene_id_obj)

    # Final sanity: ensure no unresolved numeric tokens remain
    unresolved: List[str] = []
    final_tokens = set(re.findall(r"`[^`]+`|[A-Za-z_.][A-Za-z0-9_.]*|[0-9]+", rewritten))
    for tok in final_tokens:
        if tok in new_pheno.columns:
            continue
        if tok.isdigit():
            unresolved.append(tok)
            continue
        if symbol_lookup:
            s_tok = _sanitize_name(tok)
            # detect if token looks like a symbol we might have intended
            for sym in symbol_lookup.keys():
                if sym is None:
                    continue
                s_sym = _sanitize_name(str(sym))
                if s_tok == s_sym:
                    unresolved.append(tok)
                    break

    if unresolved:
        raise ValueError(
            "Could not resolve formula gene references in RHS: {}. "
            "Ensure those GeneIDs exist in the expression matrix after filtering (non-zeros, taxon, etc.).".format(
                ", ".join(sorted(unresolved))
            )
        )

    return rewritten, new_pheno


def run_limma_pipeline(
    edata: pd.DataFrame,
    pheno: pd.DataFrame,
    *,
    group: Optional[str],
    formula: Optional[str],
    block: Optional[str],
    contrasts: Optional[str],
    logger: Optional[logging.Logger] = None,
    target_gene_ids: Optional[Sequence[str]] = None,
    symbol_lookup: Optional[Mapping[str, Sequence[str]]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Execute a limma linear model using the supplied design specification.

    Parameters
    ----------
    edata:
        Expression matrix (rows = features, columns = samples) already aligned
        with the phenotype table.
    pheno:
        Phenotype / sample metadata table with rows ordered identically to
        ``edata`` columns.
    group, formula:
        Mutually exclusive design specifiers.  When ``group`` is provided,
        ``~0+group`` is used.  Otherwise ``formula`` is passed straight through
        to ``model.matrix``.
    block:
        Optional column name in ``pheno`` used for duplicate correlation.
    contrasts:
        Either ``None`` (auto-generate pairwise contrasts) or a comma-separated
        string describing explicit contrast expressions.
    logger:
        Optional logger instance for debug output.
    target_gene_ids:
        Optional ordered list of GeneIDs; when provided the expression matrix
        is restricted to those rows before fitting.

    Returns
    -------
    Dict[str, DataFrame]
        Mapping of contrast label to result table with log-fold change,
        adjusted p-values, and the original expression columns joined back on.
    """
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    logger = logger or logging.getLogger(__name__)

    importr("limma")
    importr("dplyr", on_conflict="warn")

    # Optionally inject covariates derived from gene expression referenced in formula
    if logger:
        logger.info("Limma: incoming formula: %s", formula)
    formula, pheno = _inject_gene_covariates_in_formula(formula, edata, pheno, symbol_lookup, logger)
    if logger:
        logger.info("Limma: rewritten formula: %s", formula)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        robjects.r.assign("edata", edata)
        robjects.r.assign("pheno", pheno)

    robjects.r("mod0 <- model.matrix(~1, pheno)")

    if target_gene_ids:
        missing = [gid for gid in target_gene_ids if gid not in edata.index]
        if missing:
            raise ValueError(f"Expression matrix missing GeneIDs referenced in formula: {missing}")
        edata = edata.loc[target_gene_ids]
        with localconverter(robjects.default_converter + pandas2ri.converter):
            robjects.r.assign("edata", edata)

    if group and not formula:
        robjects.r(f"mod <- model.matrix(~0+{group}, pheno)")
    elif formula:
        robjects.r(f"mod <- model.matrix({formula}, pheno)")
    else:
        raise ValueError("Must specify one of `group` or `formula` for limma.")

    # Ensure downstream columns are valid variable names.
    edata_cols = list(robjects.r("colnames(edata)"))
    edata_fixed = _sanitize_names(edata_cols)
    if edata_cols != edata_fixed:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            robjects.r.assign("edata_fixed", edata_fixed)
        robjects.r("colnames(edata) <- edata_fixed")

    design_terms = list(robjects.r("colnames(mod)"))
    fixed_terms = _sanitize_names(design_terms)
    if design_terms != fixed_terms:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            robjects.r.assign("fixed_vars", fixed_terms)
        robjects.r("colnames(mod) <- fixed_vars")
    else:
        fixed_terms = design_terms  # already safe

    if block:
        robjects.r(f'block <- as.factor(pheno[["{block}"]])')
        robjects.r("corfit <- duplicateCorrelation(edata, design = mod, block = block)")
        robjects.r("cor <- corfit$consensus")
    else:
        robjects.r("block <- NULL")
        robjects.r("cor <- NULL")

    # Resolve contrasts; if only a single non-intercept term and no contrasts supplied,
    # analyse the direct coefficient for that term instead of building a contrast.
    non_intercept_terms = [t for t in fixed_terms if "Intercept" not in t]
    use_direct_coef = False
    contrast_list = _resolve_contrasts(fixed_terms, contrasts)
    if not contrasts and len(non_intercept_terms) == 1:
        contrast_list = non_intercept_terms
        use_direct_coef = True
    if logger:
        logger.info(
            "Limma: terms=%s, contrasts=%s, direct_coef=%s",
            non_intercept_terms,
            contrast_list,
            use_direct_coef,
        )
    if logger:
        logger.info("limma model matrix terms: %s", fixed_terms)
        logger.info("Contrasts: %s", contrast_list)

    if not contrast_list:
        raise ValueError("No contrasts specified for limma analysis.")

    if use_direct_coef:
        robjects.r("""
fit <- lmFit(as.matrix(edata), mod, block = block, cor = cor)
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
""")
    else:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            robjects.r.assign("contrasts_array", contrast_list)
        robjects.r(
            """
contrasts_matrix <- makeContrasts(contrasts = contrasts_array, levels = mod)
fit <- lmFit(as.matrix(edata), mod, block = block, cor = cor)
fit2 <- contrasts.fit(fit, contrasts_matrix)
fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
"""
        )

    results: Dict[str, pd.DataFrame] = {}
    for idx, contrast_label in enumerate(contrast_list, start=1):
        table = robjects.r(
            f"topTable(fit2, n=Inf, sort.by='none', coef={idx}, confint=TRUE)"
        )
        try:
            result_df = pandas2ri.ri2py(table)
        except AttributeError:
            result_df = pandas2ri.rpy2py(table)

        result_df = result_df.rename(
            columns={
                "logFC": "log2_FC",
                "adj.P.Val": "pAdj",
                "P.Value": "pValue",
            }
        )

        for col in ("log2_FC", "CI.L", "CI.R"):
            if col in result_df.columns:
                result_df[col] = result_df[col].apply(lambda x: x / np.log10(2))

        # Preserve the original order of the expression matrix.
        if len(result_df.index) == len(edata.index):
            result_df.index = pd.Index(edata.index)
        else:
            try:
                result_df.index = result_df.index.astype(edata.index.dtype, copy=False)
            except (TypeError, ValueError):
                pass

        # Join the expression values for downstream export compatibility.
        result_df = result_df.join(edata, how="left")
        results[contrast_label] = result_df

    return results
