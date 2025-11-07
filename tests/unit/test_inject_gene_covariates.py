import pandas as pd
import pytest

pytest.importorskip("rpy2")

from tackle.statmodels.limma_runner import _inject_gene_covariates_in_formula


class _DummyLogger:
    def __init__(self):
        self.warnings = []
        self.infos = []

    def warning(self, msg, *args):
        if args:
            msg = msg % args
        self.warnings.append(msg)

    def info(self, msg, *args):
        if args:
            msg = msg % args
        self.infos.append(msg)


def test_inject_gene_covariates_skips_warning_for_non_gene_tokens():
    edata = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=["101", "202"],
        columns=["s1", "s2"],
    )
    pheno = pd.DataFrame({"condition": ["A", "B"]}, index=["s1", "s2"])
    logger = _DummyLogger()

    rewritten, new_pheno = _inject_gene_covariates_in_formula(
        formula="~ condition",
        edata=edata,
        pheno=pheno,
        symbol_lookup={},
        logger=logger,
    )

    assert rewritten == "~ condition"
    assert new_pheno.equals(pheno)
    assert logger.warnings == []


def test_inject_gene_covariates_warns_for_numeric_gene_tokens():
    edata = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=["101", "202"],
        columns=["s1", "s2"],
    )
    pheno = pd.DataFrame(index=["s1", "s2"])
    logger = _DummyLogger()

    rewritten, new_pheno = _inject_gene_covariates_in_formula(
        formula="~ 101",
        edata=edata,
        pheno=pheno,
        symbol_lookup={},
        logger=logger,
    )

    assert rewritten == "~ GID_101"
    assert "GID_101" in new_pheno.columns
    assert logger.warnings, "expected a warning about deprecated numeric tokens"
    assert "GeneID" in logger.warnings[0]
