import os
from pathlib import Path

import pandas as pd
import pytest

from types import SimpleNamespace

pytest.importorskip("rpy2")
from rpy2.rinterface_lib.embedded import RRuntimeError

import tackle.volcanoplot as volcanoplot


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


def test_volcanoplot_writes_safe_exports(tmp_path, stub_gene_mapper):
    data_obj = _StubDataObj(tmp_path)
    ctx = SimpleNamespace(obj={"data_obj": data_obj, "file_fmts": [".png"]})

    try:
        volcanoplot.volcanoplot(
            ctx=ctx,
            foldchange=1.5,
            expression_data=False,
            number=10,
        )
    except RRuntimeError as err:
        pytest.skip(f"R runtime prerequisites missing: {err}")

    csv_files = list((Path(data_obj.outpath) / "volcano").rglob("*.tsv"))
    assert csv_files, "Expected volcanoplot to write a TSV export"

    csv_path = csv_files[0]
    assert csv_path.suffix == ".tsv"
    assert not csv_path.name.endswith(".tsv.tsv")

    limit = os.pathconf(str(csv_path.parent), "PC_NAME_MAX")
    assert len(csv_path.name.encode("utf-8")) <= limit

    png_files = list((Path(data_obj.outpath) / "volcano").rglob("*.png"))
    assert png_files, "Expected a PNG volcano plot to be generated"

    for plot_path in png_files:
        assert plot_path.exists()
        assert plot_path.suffix == ".png"
        assert not plot_path.name.endswith(".png.png")
        plot_limit = os.pathconf(str(plot_path.parent), "PC_NAME_MAX")
        assert len(plot_path.name.encode("utf-8")) <= plot_limit
