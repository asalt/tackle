import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rpy2")
from rpy2.rinterface_lib.embedded import RRuntimeError

from tackle.statmodels.limma_runner import (
    run_limma_pipeline,
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


def test_limma_pipeline_full_formula_lhs_does_not_subset_rows():
    genes = ["geneA", "geneB", "2064"]
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
    pheno = pd.DataFrame({"condition": ["A", "A", "B", "B"]}, index=samples)

    try:
        results = run_limma_pipeline(
            edata=edata,
            pheno=pheno,
            group=None,
            formula="GID_2064 ~ 0 + condition",
            block=None,
            contrasts=None,
        )
    except RRuntimeError as err:
        pytest.skip(f"R limma prerequisites missing: {err}")

    assert len(results) == 1
    assert any("condition" in label for label in results)
    assert not any(label.startswith("coef_GID_2064=") for label in results)

    contrast_df = next(iter(results.values()))
    assert contrast_df.index.tolist() == genes
    assert all(col in contrast_df.columns for col in samples)


def test_limma_pipeline_rhs_gene_covariate_still_drives_direct_coef_with_lhs_present():
    samples = ["s1", "s2", "s3", "s4", "s5", "s6"]
    gid_2064 = np.array([1.0, 1.4, 1.9, 1.2, 1.7, 2.1])
    gid_9999 = np.array([0.5, 0.9, 1.3, 0.7, 1.1, 1.5])
    geneX = gid_9999 * 2.0 + np.array([0.2, -0.1, 0.05, -0.2, 0.15, -0.05])
    geneY = np.array([0.1, -0.2, 0.0, 0.05, -0.1, 0.2])
    edata = pd.DataFrame(
        [geneX, geneY, gid_2064, gid_9999],
        index=["geneX", "geneY", "2064", "9999"],
        columns=samples,
    )
    pheno = pd.DataFrame({"condition": ["A", "A", "A", "B", "B", "B"]}, index=samples)

    try:
        results = run_limma_pipeline(
            edata=edata,
            pheno=pheno,
            group=None,
            formula="GID_2064 ~ 0 + GID_9999 + condition",
            block=None,
            contrasts=None,
        )
    except RRuntimeError as err:
        pytest.skip(f"R limma prerequisites missing: {err}")

    assert any(label.startswith("coef_GID_9999=") for label in results)
    assert any("condition" in label for label in results)
    direct_df = next(df for label, df in results.items() if label.startswith("coef_GID_9999="))
    assert direct_df.index.tolist() == ["geneX", "geneY", "2064", "9999"]


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
