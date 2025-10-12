import os
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

from types import SimpleNamespace

import tackle.volcanoplot as volcanoplot


class _DummyConverter:
    def __add__(self, other):
        return self


class _DummyContext:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


class _DummyGrDevices:
    def __init__(self):
        self.open_calls = []

    def png(self, **kwargs):
        self.open_calls.append(("png", kwargs))

    def pdf(self, **kwargs):
        self.open_calls.append(("pdf", kwargs))

    def svg(self, **kwargs):
        self.open_calls.append(("svg", kwargs))

    def dev_off(self):
        return None


class _DummyR(dict):
    def __init__(self):
        super().__init__()
        self["source"] = lambda path: None
        self["volcanoplot"] = lambda *args, **kwargs: None


@pytest.fixture(autouse=True)
def fake_rpy2(monkeypatch):
    converter = _DummyConverter()
    pandas2ri = SimpleNamespace(converter=converter)

    robjects_mod = types.ModuleType("rpy2.robjects")
    robjects_mod.r = _DummyR()
    robjects_mod.NA_Real = float("nan")
    robjects_mod.default_converter = converter
    robjects_mod.pandas2ri = pandas2ri

    packages_mod = types.ModuleType("rpy2.robjects.packages")
    dummy_devices = _DummyGrDevices()
    packages_mod.importr = lambda name: dummy_devices

    conversion_mod = types.ModuleType("rpy2.robjects.conversion")
    conversion_mod.localconverter = lambda _converter: _DummyContext()

    monkeypatch.setitem(sys.modules, "rpy2", types.ModuleType("rpy2"))
    monkeypatch.setitem(sys.modules, "rpy2.robjects", robjects_mod)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.packages", packages_mod)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.conversion", conversion_mod)

    return dummy_devices


@pytest.fixture
def stub_gene_mapper(monkeypatch):
    gm = SimpleNamespace(
        symbol={"geneA": "SYM_A", "geneB": "SYM_B"},
        description={"geneA": "desc A", "geneB": "desc B"},
    )
    monkeypatch.setattr(volcanoplot, "get_gene_mapper", lambda: gm)
    return gm


class _StubDataObj:
    def __init__(self, outdir: Path):
        self.group = "condition"
        self.col_metadata = pd.DataFrame(
            {"condition": ["g0", "g0", "g1", "g1"]},
            index=["s1", "s2", "s3", "s4"],
        )
        self.outpath_name = "analysis"
        self.taxon = "test_taxon"
        self.non_zeros = "N" * 260  # trigger NAME_MAX trimming
        self.colors_only = False
        self.batch_applied = False
        self.normtype = "norm"
        self.batch_nonparametric = False
        self.outpath = str(outdir)
        self.gid_symbol = {"geneA": "SYM_A", "geneB": "SYM_B"}
        self.gid_funcat_mapping = {}
        self.areas_log_shifted = pd.DataFrame()

    def stat_model(self, **kwargs):
        df = pd.DataFrame(
            {
                "log2_FC": [1.5, -1.2],
                "pAdj": [0.01, 0.2],
                "pValue": [0.01, 0.2],
                "t": [2.1, -1.0],
                "FunCats": ["", ""],
            },
            index=["geneA", "geneB"],
        )
        return {"g1 - g0": df}


def test_volcanoplot_writes_safe_csv(tmp_path, stub_gene_mapper):
    data_obj = _StubDataObj(tmp_path)
    ctx = SimpleNamespace(obj={"data_obj": data_obj, "file_fmts": []})

    volcanoplot.volcanoplot(
        ctx=ctx,
        foldchange=1.5,
        expression_data=False,
        number=10,
    )

    csv_files = list((Path(data_obj.outpath) / "volcano").rglob("*.tsv"))
    assert csv_files, "Expected volcanoplot to write a TSV export"

    csv_path = csv_files[0]
    assert csv_path.suffix == ".tsv"
    assert not csv_path.name.endswith(".tsv.tsv")

    limit = os.pathconf(str(csv_path.parent), "PC_NAME_MAX")
    assert len(csv_path.name.encode("utf-8")) <= limit
