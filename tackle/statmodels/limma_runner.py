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

from ..formula_utils import extract_geneid_key as _extract_geneid_key, unquote_backticks as _unquote_backticks


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
        # Support explicit GeneID separators on LHS as well
        gid_key = _extract_geneid_key(token_str)
        if gid_key is not None:
            if gid_key in index_map:
                gene_ids.append(index_map[gid_key])
                continue
            else:
                raise ValueError(f"Unknown GeneID reference '{token_str}' in formula left-hand side.")
        # Support GID tokens (GID_<id> or GID<id>) on LHS as explicit GeneIDs
        m_gid = re.match(r"(?i)^GID_?([A-Za-z0-9_.-]+)$", token_str)
        if m_gid:
            gid_raw = m_gid.group(1)
            if gid_raw in index_map:
                gene_ids.append(index_map[gid_raw])
                continue
            else:
                raise ValueError(f"Unknown GID reference '{token_str}' in formula left-hand side.")
        # Disallow bare numeric tokens on LHS (ambiguous; could be intercept/constant)
        if token_str.isdigit():
            raise ValueError(
                f"Ambiguous numeric token '{token_str}' in formula left-hand side. "
                f"Use explicit 'GeneID_{token_str}' or a gene symbol."
            )
        if token_str in index_map:
            gene_ids.append(index_map[token_str])
            continue
        # Avoid symbol lookup for purely numeric tokens
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
    safe_vars_injected: List[str] = []
    for tok in sorted(tokens, key=len, reverse=True):
        # Skip if token already a pheno column
        if tok in new_pheno.columns:
            continue

        gene_id_obj: Optional[Any] = None
        # Support explicit GeneID separators (GeneID_<id>, GeneID:<id>, GeneID-<id>)
        tok_unquoted = _unquote_backticks(tok)
        gid_key = _extract_geneid_key(tok_unquoted)
        if gid_key is not None:
            if gid_key in edata_index_map:
                gene_id_obj = edata_index_map[gid_key]
            else:
                if logger:
                    logger.warning(
                        "Explicit GeneID reference '%s' did not match any expression row; check filtering and IDs.",
                        _unquote_backticks(tok),
                    )
                else:
                    print(
                        f"Warning: explicit GeneID reference '{_unquote_backticks(tok)}' not found in expression rows."
                    )
                # leave unresolved to be caught below
        # Support already-sanitized variable tokens (GID_<id> or GID<id>) as shorthand
        if gene_id_obj is None:
            m_gid = re.match(r"(?i)^GID_?([A-Za-z0-9_.-]+)$", tok_unquoted)
            if m_gid:
                gid_raw = m_gid.group(1)
                if gid_raw in edata_index_map:
                    gene_id_obj = edata_index_map[gid_raw]
                else:
                    # Don't resolve; will be handled as unresolved below
                    pass
        else:
            tok_unquoted = _unquote_backticks(tok)
            # Do not interpret bare numeric tokens as GeneIDs; only explicit GeneID_ or symbols are allowed.
            # Historically this was allowed; removing to avoid accidental matches with intercepts or df constants.
        # Symbol lookup (case-insensitive + sanitized)
        if gene_id_obj is None and symbol_lookup:
            # Skip symbol lookup for purely numeric tokens (e.g., 0/1/3 in ns(...))
            if tok_unquoted.isdigit():
                continue
            cand_ids = []
            search = {tok_unquoted, _sanitize_name(tok_unquoted), tok_unquoted.casefold(), _sanitize_name(tok_unquoted).casefold()}
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
        safe_vars_injected.append(safe_var)
        if logger:
            logger.info("Injected covariate %s from GeneID %s for formula", safe_var, gene_id_obj)

    # Final sanity: only flag unresolved tokens that look like gene references
    unresolved: List[str] = []
    final_tokens = set(re.findall(r"`[^`]+`|[A-Za-z_.][A-Za-z0-9_.]*|[0-9]+", rewritten))
    sym_keys_sanitized = set()
    if symbol_lookup:
        for sk in symbol_lookup.keys():
            if sk is None:
                continue
            sym_keys_sanitized.add(_sanitize_name(str(sk)))
    for tok in final_tokens:
        if tok in new_pheno.columns:
            continue
        # Candidate gene reference if (a) token matches edata index keys we considered, or
        # (b) matches a known symbol key; otherwise leave numeric constants (e.g., df=3) alone.
        # Ignore pure numeric constants and explicit GeneID_ references that failed lookup
        if tok.isdigit():
            continue
        gid_key = _extract_geneid_key(tok)
        if gid_key is not None:
            # Explicit GeneID reference that failed to resolve should raise
            if gid_key not in edata_index_map:
                unresolved.append(tok)
            # Skip further checks for this token
            continue
        # Also treat bare GID tokens (GID_<id> or GID<id>) as gene-like; if not injected, flag unresolved
        if re.match(r"(?i)^GID_?[A-Za-z0-9_.-]+$", tok):
            gid_raw = _unquote_backticks(tok)
            gid_raw = re.sub(r"(?i)^GID_?", "", gid_raw)
            if gid_raw not in edata_index_map:
                unresolved.append(tok)
            continue
        else:
            candidate_gene = tok in edata_index_map or _sanitize_name(tok) in sym_keys_sanitized
        if candidate_gene:
            unresolved.append(tok)

    if unresolved:
        raise ValueError(
            "Could not resolve formula gene references in RHS: {}. "
            "Ensure those GeneIDs exist in the expression matrix after filtering (non-zeros, taxon, etc.).".format(
                ", ".join(sorted(unresolved))
            )
        )

    # If the formula contains explicit subtraction of injected gene covariates, rewrite those
    # as additive and rely on implicit pairwise contrasts to generate A - B later.
    # Example: '~ 1 + GID_2064 - GID_1956' -> '~ 1 + GID_2064 + GID_1956'
    if safe_vars_injected:
        # Replace only minus directly preceding a GID_ var
        pattern = re.compile(r"-\s*(?=(?:" + "|".join(map(re.escape, set(safe_vars_injected))) + r")\b)")
        new_rewritten = pattern.sub("+ ", rewritten)

        # Also replace minus preceding a function call that contains a GID_ var, e.g. '- scale(GID_123)'
        def _replace_func_negations(s: str) -> str:
            out = []
            i = 0
            n = len(s)
            safe_set = set(safe_vars_injected)
            while i < n:
                m = re.search(r"-\s*([A-Za-z_.][A-Za-z0-9_.]*)\(", s[i:])
                if not m:
                    out.append(s[i:])
                    break
                start = i + m.start()  # index of '-'
                lp = i + m.end() - 1   # index of '('
                # Find matching ')'
                depth = 1
                j = lp + 1
                while j < n and depth > 0:
                    c = s[j]
                    if c == '(':
                        depth += 1
                    elif c == ')':
                        depth -= 1
                    j += 1
                if depth != 0:
                    # Unbalanced; give up and append rest
                    out.append(s[i:])
                    break
                inner = s[lp + 1 : j - 1]
                contains_gid = any(re.search(r"\b" + re.escape(g) + r"\b", inner) for g in safe_set)
                # Prefix before '-'
                out.append(s[i:start])
                call = s[start:j]
                if contains_gid:
                    # flip leading '- ' to '+ '
                    call = re.sub(r"^-\s*", "+ ", call)
                out.append(call)
                i = j
            return ''.join(out)

        newer = _replace_func_negations(new_rewritten)
        if newer != rewritten and logger:
            logger.info(
                "Limma: rewritten formula to include subtracted gene covariates additively for contrasts: %s",
                newer,
            )
        rewritten = newer

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
    joint_export_dir: Optional[str] = None,
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
    # Load splines if using ns() in the formula
    try:
        if formula and "ns(" in str(formula):
            importr("splines")
            if logger:
                logger.info("Limma: loaded 'splines' for ns() in formula")
    except Exception:
        pass
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
    # Ensure R-level syntactic validity (removes parentheses/commas etc.)
    robjects.r("colnames(mod) <- make.names(colnames(mod))")
    fixed_terms = list(robjects.r("colnames(mod)"))

    if block:
        robjects.r(f'block <- as.factor(pheno[["{block}"]])')
        robjects.r("corfit <- duplicateCorrelation(edata, design = mod, block = block)")
        robjects.r("cor <- corfit$consensus")
    else:
        robjects.r("block <- NULL")
        robjects.r("cor <- NULL")

    # Resolve contrasts and direct coefs with mixed models (factors + gene covariates)
    non_intercept_terms = [t for t in fixed_terms if "Intercept" not in t]

    # Identify gene covariate variables present in the formula (e.g., GID_2064)
    gid_vars = sorted(set(re.findall(r"GID_[0-9]+", str(formula) if formula else "")))
    gene_terms = [t for t in non_intercept_terms if any(gid in t for gid in gid_vars)]

    # Tokenize the (rewritten) formula to find pheno variables referenced
    token_candidates = set(re.findall(r"`[^`]+`|[A-Za-z_.][A-Za-z0-9_.]*", str(formula) if formula else ""))
    token_candidates = {_sanitize_name(_unquote_backticks(tok)) for tok in token_candidates}
    token_candidates -= {"ns", "poly", "hinge", "scale", "I", "C"}

    # Map sanitized pheno columns (including injected GID_*)
    pheno_cols_sani = {_sanitize_name(str(c)) for c in pheno.columns}
    used_pheno_tokens = [tok for tok in token_candidates if tok in pheno_cols_sani]

    # Build factor groups from design terms by token prefix; continuous -> single term
    factor_groups: Dict[str, List[str]] = {}
    continuous_terms: List[str] = []
    for tok in used_pheno_tokens:
        if tok.startswith("GID_"):
            continue  # handled as gene_terms
        # Terms that match this variable: exact or prefixed (for 0+factor)
        matches = [t for t in non_intercept_terms if (t == tok or t.startswith(tok))]
        if len(matches) >= 2:
            factor_groups[tok] = matches
        elif len(matches) == 1:
            continuous_terms.append(matches[0])

    # Direct coefficients gather: all gene terms and any single continuous pheno terms
    direct_coef_list: List[str] = []
    direct_coef_list.extend(gene_terms)
    direct_coef_list.extend(continuous_terms)
    # If only one non-intercept term total, prefer direct coef
    if not used_pheno_tokens and len(non_intercept_terms) == 1:
        direct_coef_list = list(non_intercept_terms)

    # Contrasts: explicit or pairwise within each factor group
    contrast_list: List[str] = []
    if contrasts:
        contrast_list = _resolve_contrasts(fixed_terms, contrasts)
    else:
        for _, terms in factor_groups.items():
            for a, b in itertools.combinations([t for t in terms if t in fixed_terms], 2):
                contrast_list.append(f"{a} - {b}")

    # Deduplicate while preserving order
    def _uniq(seq: Iterable[str]) -> List[str]:
        out, seen = [], set()
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out
    direct_coef_list = _uniq(direct_coef_list)
    contrast_list = _uniq(contrast_list)

    if logger:
        logger.info(
            "Limma: terms=%s, contrasts=%s, direct_coef=%s",
            non_intercept_terms,
            contrast_list,
            direct_coef_list,
        )
        logger.info("limma model matrix terms: %s", fixed_terms)
        logger.info("Contrasts: %s", contrast_list)

    if not contrast_list and not direct_coef_list:
        raise ValueError("No contrasts or coefficients selected for limma analysis.")

    robjects.r("print(mod)")
    # Fit once; reuse for direct coefs. Compute separate fit2 for contrasts if any.
    robjects.r(
        """
fit <- lmFit(as.matrix(edata), mod, block = block, cor = cor)
fit2_base <- eBayes(fit, robust = TRUE, trend = TRUE)
"""
    )
    has_contrasts = bool(contrast_list)
    if has_contrasts:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            robjects.r.assign("contrasts_array", contrast_list)
        robjects.r(
            """
contrasts_matrix <- makeContrasts(contrasts = contrasts_array, levels = mod)
fit2_con <- contrasts.fit(fit, contrasts_matrix)
fit2_con <- eBayes(fit2_con, robust = TRUE, trend = TRUE)
"""
        )

    results: Dict[str, pd.DataFrame] = {}

    # Helper to convert an R topTable to our pandas df and join edata
    def _collect_result(table_obj) -> pd.DataFrame:
        try:
            result_df = pandas2ri.ri2py(table_obj)
        except AttributeError:
            result_df = pandas2ri.rpy2py(table_obj)
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
        if len(result_df.index) == len(edata.index):
            result_df.index = pd.Index(edata.index)
        else:
            try:
                result_df.index = result_df.index.astype(edata.index.dtype, copy=False)
            except (TypeError, ValueError):
                pass
        return result_df.join(edata, how="left")

    # Add direct coefficient results (label as 'coef_<term>=<term>' so volcano names improve)
    for coef_label in direct_coef_list:
        table = robjects.r(
            f"topTable(fit2_base, n=Inf, sort.by='none', coef='{coef_label}', confint=TRUE)"
        )
        results[f"coef_{coef_label}={coef_label}"] = _collect_result(table)

    # Add contrast results (label as '<expr>=<expr>' so volcano names include expr and parsing works)
    for contrast_label in contrast_list:
        if not has_contrasts:
            continue
        table = robjects.r(
            f"topTable(fit2_con, n=Inf, sort.by='none', coef='{contrast_label}', confint=TRUE)"
        )
        results[f"{contrast_label}={contrast_label}"] = _collect_result(table)

    # Joint F-tests for any GID variable that produced multiple basis terms
    if bool(formula) and gid_vars and joint_export_dir:
        try:
            import os
            os.makedirs(joint_export_dir, exist_ok=True)
            for cov in gid_vars:
                basis_terms = [t for t in non_intercept_terms if cov in t]
                if len(basis_terms) <= 1:
                    continue
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    robjects.r.assign("joint_coef", basis_terms)
                jtable = robjects.r(
                    "topTable(fit2_base, n=Inf, sort.by='F', coef=joint_coef)"
                )
                try:
                    jdf = pandas2ri.ri2py(jtable)
                except AttributeError:
                    jdf = pandas2ri.rpy2py(jtable)
                jdf = jdf.join(edata, how="left")
                if 'F' in jdf.columns and 't' not in jdf.columns:
                    jdf['t'] = np.sqrt(jdf['F'])
                out_path = os.path.join(joint_export_dir, f"joint_F_{cov}.tsv")
                jdf.to_csv(out_path, sep='\t')
                if logger:
                    logger.info("Limma: wrote joint F-test to %s (terms=%s)", out_path, basis_terms)
        except Exception as e:
            if logger:
                logger.warning("Limma: joint F-test export failed: %s", e)

    return results
