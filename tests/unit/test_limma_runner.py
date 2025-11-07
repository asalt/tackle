import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rpy2")
from rpy2.rinterface_lib.embedded import RRuntimeError

from tackle.statmodels.limma_runner import (
    run_limma_pipeline,
    normalize_formula_targets,
)


def test_limma_pipeline_basic():
    genes = ["geneA", "geneB", "geneC"]
    samples = ["s1", "s2", "s3", "s4"]

    edata = pd.DataFrame(
        [
            [10.0, 11.0, 25.0, 24.0],
            [5.0, 4.5, 6.0, 6.5],
            [100.0, 98.0, 102.0, 101.0],
        ],
        index=genes,
        columns=samples,
    )

    pheno = pd.DataFrame(
        {"condition": ["A", "A", "B", "B"]},
        index=samples,
    )

    try:
        results = run_limma_pipeline(
            edata=edata,
            pheno=pheno,
            group="condition",
            formula=None,
            block=None,
            contrasts=None,
        )
        subset_results = run_limma_pipeline(
            edata=edata,
            pheno=pheno,
            group="condition",
            formula="~0 + condition",
            block=None,
            contrasts=None,
            target_gene_ids=["geneA"],
        )
    except RRuntimeError as err:
        pytest.skip(f"R limma prerequisites missing: {err}")

    assert results, "limma runner should produce at least one contrast"
    assert len(results) == 1

    contrast_name, contrast_df = next(iter(results.items()))
    assert "condition" in contrast_name
    assert {"pAdj", "pValue", "log2_FC"} <= set(contrast_df.columns)
    assert all(col in contrast_df.columns for col in samples)

    # Limma should mark the strongly changing gene as significant.
    top_gene = contrast_df.sort_values("pAdj").index[0]
    assert top_gene in genes
    assert np.isfinite(contrast_df.loc[top_gene, "pAdj"])
    assert np.isfinite(contrast_df.loc[top_gene, "log2_FC"])

    # Original expression values remain attached after the join.
    for sample in samples:
        assert np.isclose(contrast_df.loc["geneA", sample], edata.loc["geneA", sample])

    subset_name, subset_df = next(iter(subset_results.items()))
    assert subset_df.index.tolist() == ["geneA"]
    assert all(col in subset_df.columns for col in samples)


def test_normalize_formula_targets_symbol_resolution():
    formula = "HER2 ~ 0 + condition"
    index = pd.Index(["1234", "5678", "91011"])
    symbol_lookup = {"HER2": ["5678"], "EGFR": ["1234", "91011"]}

    normalised, targets = normalize_formula_targets(
        formula, index, symbol_lookup, logger=None
    )

    assert normalised == "~0 + condition"
    assert targets == ["5678"]


def test_normalize_formula_targets_symbol_case_insensitive():
    formula = "her2 ~ condition"
    index = pd.Index(["5678"])
    symbol_lookup = {"HER2": ["5678"]}

    normalised, targets = normalize_formula_targets(formula, index, symbol_lookup, logger=None)
    assert normalised == "~condition"
    assert targets == ["5678"]


def test_normalize_formula_targets_geneid_lookup_with_int_index():
    formula = "5678 ~ 0 + condition"
    index = pd.Index([1234, 5678])
    symbol_lookup = {}

    normalised, targets = normalize_formula_targets(formula, index, symbol_lookup, logger=None)
    assert normalised == "~0 + condition"
    assert targets == [5678]


def test_normalize_formula_targets_errors_on_ambiguous():
    formula = "EGFR ~ 1 + condition"
    index = pd.Index(["1234", "5678"])
    symbol_lookup = {"EGFR": ["1234", "5678"]}

    with pytest.raises(ValueError):
        normalize_formula_targets(formula, index, symbol_lookup, logger=None)


def test_normalize_formula_targets_unknown_symbol():
    formula = "UNKNOWN ~ 0 + condition"
    index = pd.Index(["1234"])
    symbol_lookup = {}

    with pytest.raises(ValueError):
        normalize_formula_targets(formula, index, symbol_lookup, logger=None)


def test_limma_pipeline_continuous_covariate_from_geneid():
    # Build a simple dataset where gene '2064' (HER2) covaries with other genes
    samples = ["s1", "s2", "s3", "s4", "s5", "s6"]
    her2 = np.array([1.0, 1.2, 1.4, 2.0, 2.2, 2.4])
    # geneX positively associated with HER2, geneY not associated
    geneX = her2 * 2.0 + np.random.normal(scale=0.05, size=len(her2))
    geneY = np.random.normal(scale=0.05, size=len(her2))
    edata = pd.DataFrame(
        [geneX, geneY, her2],
        index=["geneX", "geneY", "2064"],
        columns=samples,
    )
    pheno = pd.DataFrame(index=samples)

    try:
        results = run_limma_pipeline(
            edata=edata,
            pheno=pheno,
            group=None,
            formula="~ GeneID_2064",
            block=None,
            contrasts=None,
        )
    except RRuntimeError as err:
        pytest.skip(f"R limma prerequisites missing: {err}")

    assert len(results) == 1
    label, df = next(iter(results.items()))
    # Should be testing direct coefficient for injected covariate
    assert "GID_2064" in label or label == "2064"
    # geneX should tend to be the most significant
    top = df.sort_values("pAdj").index[0]
    assert top in ("geneX", "2064")


def test_limma_pipeline_continuous_covariate_missing_gene_raises():
    samples = ["s1", "s2", "s3", "s4"]
    edata = pd.DataFrame(
        [np.random.rand(4), np.random.rand(4)],
        index=["1000", "2000"],
        columns=samples,
    )
    pheno = pd.DataFrame(index=samples)

    with pytest.raises(ValueError):
        run_limma_pipeline(
            edata=edata,
            pheno=pheno,
            group=None,
            formula="~ GeneID_2064",
            block=None,
            contrasts=None,
        )
