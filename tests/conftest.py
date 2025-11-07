import numpy as np
import pandas as pd
import pytest
from types import SimpleNamespace


class _DF(pd.DataFrame):
    @property
    def _constructor(self):
        return _DF

    @property
    def GeneID(self):
        # Provide attribute used by dispatcher to create a 'GeneID' column from the index
        return pd.Series(self.index, index=self.index)


class StubDataObj:
    def __init__(self, outdir, n_genes=30, n_samples=8, seed=42):
        rng = np.random.default_rng(seed)
        genes = [f"g{i+1}" for i in range(n_genes)]
        samples = [f"s{i+1}" for i in range(n_samples)]

        # Two conditions A/B across samples
        conditions = ["A"] * (n_samples // 2) + ["B"] * (n_samples - n_samples // 2)
        self.col_metadata = pd.DataFrame({"condition": conditions}, index=samples)

        expr = rng.normal(loc=0.0, scale=1.0, size=(n_genes, n_samples))
        self.areas_log = _DF(expr, index=genes, columns=samples)
        self.mask = pd.DataFrame(False, index=genes, columns=samples)
        self.gid_symbol = {g: g.upper() for g in genes}

        # Minimal fields consumed by naming and dispatcher
        self.outpath_name = "analysis"
        self.taxon = "test_taxon"
        self.non_zeros = 1
        self.colors_only = False
        self.batch_applied = False
        self.batch_nonparametric = False
        self.normtype = "norm"
        self.outpath = str(outdir)
        self.metadata_colors = None
        self.annotations = None
        self.normed = True
        self.data = pd.DataFrame()


@pytest.fixture
def stub_gene_mapper(monkeypatch):
    import tackle.clusterplot_dispatcher as cdisp
    gm = SimpleNamespace(symbol={})
    monkeypatch.setattr(cdisp, "get_gene_mapper", lambda: gm)
    return gm


@pytest.fixture
def stub_data_obj(tmp_path, stub_gene_mapper):
    return StubDataObj(tmp_path)


@pytest.fixture
def ctx(stub_data_obj):
    return SimpleNamespace(obj={"data_obj": stub_data_obj, "file_fmts": [".png"]})

