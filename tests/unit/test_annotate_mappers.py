from pathlib import Path

import pandas as pd


def test_annotation_table_get_annotation_filters_membership():
    from tackle.annotate import AnnotationTable, AnnotationTableConfig

    annot = AnnotationTable(AnnotationTableConfig(path=Path("dummy.tsv")))
    annot._df = pd.DataFrame(
        {
            "GeneID": ["1", "2", "3", "4"],
            "GeneSymbol": ["A", "B", "C", "D"],
            "SECRETED": ["SECRETED", "", "0", "False"],
        }
    )

    sub = annot.get_annotation("SECRETED")
    assert sub["GeneID"].tolist() == ["1"]


def test_annotation_table_categories_reads_header_without_loading(tmp_path):
    from tackle.annotate import AnnotationTable, AnnotationTableConfig

    path = tmp_path / "annotations.tsv"
    path.write_text("GeneID\tGeneSymbol\tSECRETED\tNUCLEUS\n1\tA\t\t\n")

    annot = AnnotationTable(AnnotationTableConfig(path=path))
    assert annot._df is None
    assert annot.categories == ["SECRETED", "NUCLEUS"]
    assert annot._df is None


def test_homologene_mapper_map_to_human():
    from tackle.annotate import HomologeneMapper, HomologeneMapperConfig

    hom = HomologeneMapper(HomologeneMapperConfig(path=Path("dummy.data")))
    hom._df = pd.DataFrame(
        {
            "Homologene": ["1", "1", "2", "2"],
            "TaxonID": ["10090", "9606", "10090", "9606"],
            "GeneID": ["m1", "100", "m2", "200"],
            "Symbol": ["M1", "H1", "M2", "H2"],
            "ProteinGI": ["", "", "", ""],
            "ProteinAccession": ["", "", "", ""],
        }
    )

    out = hom.map_to_human(["m1", "m2", "missing"], human_taxon="9606")
    assert out == {"m1": "100", "m2": "200", "missing": None}


def test_annotation_table_map_gene_ids_homologene_fallback():
    from tackle.annotate import (
        AnnotationTable,
        AnnotationTableConfig,
        HomologeneMapper,
        HomologeneMapperConfig,
    )

    annot = AnnotationTable(AnnotationTableConfig(path=Path("dummy.tsv")))
    annot._df = pd.DataFrame(
        {
            "GeneID": ["100", "200"],
            "GeneSymbol": ["H1", "H2"],
            "SECRETED": ["SECRETED", ""],
        }
    )

    hom = HomologeneMapper(HomologeneMapperConfig(path=Path("dummy.data")))
    hom._df = pd.DataFrame(
        {
            "Homologene": ["1", "1", "2", "2"],
            "TaxonID": ["10090", "9606", "10090", "9606"],
            "GeneID": ["m1", "100", "m2", "200"],
            "Symbol": ["M1", "H1", "M2", "H2"],
            "ProteinGI": ["", "", "", ""],
            "ProteinAccession": ["", "", "", ""],
        }
    )

    res = annot.map_gene_ids(
        ["m1", "m2"],
        taxon="10090",
        fallback_to_human=True,
        homologene=hom,
    )

    assert res["GeneID"].tolist() == ["m1", "m2"]
    assert res["GeneSymbol"].tolist() == ["H1", "H2"]


def test_lazy_table_reads_xlsx_path(tmp_path, monkeypatch):
    from tackle.annotate.mappers import LazyTable, LazyTableConfig

    path = tmp_path / "table.xlsx"
    path.write_bytes(b"")

    called = {}

    def fake_read_excel(p, **kwargs):
        called["path"] = Path(p)
        called["kwargs"] = kwargs
        return pd.DataFrame({"x": [1]})

    monkeypatch.setattr(pd, "read_excel", fake_read_excel)

    table = LazyTable(LazyTableConfig(path=path, read_kwargs={"dtype": str}))
    out = table.read()
    assert out["x"].tolist() == [1]
    assert called["path"] == path
    assert called["kwargs"] == {"dtype": str}


def test_gene_mapper_properties_cache_with_slots():
    from tackle.annotate import GeneMapper, GeneMapperConfig

    gm = GeneMapper(GeneMapperConfig(path=Path("dummy.tsv")))
    gm._df = pd.DataFrame(
        {
            "GeneSymbol": pd.Series(["A", "B"], index=["g1", "g2"]),
            "TaxonID": pd.Series(["9606", "9606"], index=["g1", "g2"]),
        }
    )
    first = gm.symbol
    second = gm.symbol
    assert first is second
    assert first == {"g1": "A", "g2": "B"}
