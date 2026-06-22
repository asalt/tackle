import pandas as pd
import pytest


def test_annotations_get_annot_filters_membership(monkeypatch):
    from tackle.containers import Annotations

    annot = Annotations()
    annot._df = pd.DataFrame(
        {
            "GeneID": ["1", "2", "3", "4"],
            "GeneSymbol": ["A", "B", "C", "D"],
            "SECRETED": ["SECRETED", "", "0", "False"],
        }
    )

    sub = annot.get_annot("SECRETED")
    assert sub["GeneID"].tolist() == ["1"]


def test_annotations_get_annot_all_returns_full(monkeypatch):
    from tackle.containers import Annotations

    annot = Annotations()
    annot._df = pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["A", "B"],
            "SECRETED": ["", "SECRETED"],
        }
    )

    assert annot.get_annot("_all").shape[0] == 2


def test_annotations_get_annot_unknown_raises(monkeypatch):
    from tackle.containers import Annotations

    annot = Annotations()
    annot._df = pd.DataFrame({"GeneID": ["1"], "GeneSymbol": ["A"]})

    with pytest.raises(KeyError):
        annot.get_annot("MISSING_COL")


def test_add_annotations_handles_missing_genesymbol(monkeypatch):
    from tackle import containers

    class DummyAnnot:
        def map_gene_ids(self, gids, **kwargs):
            gids = list(gids)
            return pd.DataFrame(
                {
                    "GeneID": gids,
                    "GeneSymbol": ["A"] * len(gids),
                    "TaxonID": ["9606"] * len(gids),
                    "GeneDescription": ["desc"] * len(gids),
                    "FunCats": [""] * len(gids),
                    "SECRETED": ["SECRETED"] * len(gids),
                }
            )

    monkeypatch.setattr(containers, "get_annotation_mapper", lambda: DummyAnnot())

    df = pd.DataFrame({"GeneID": ["1", "2"], "value": [1, 2]})
    out = containers.add_annotations(df, annotations=None)
    assert "GeneSymbol" in out.columns
    assert "SECRETED" in out.columns
    assert "value" in out.columns


def test_synthetic_symbol_from_header_prefers_embedded_symbol():
    from tackle.utils import synthetic_symbol_from_header

    header = "CustomPart|geneid|990000001|taxon|32630|symbol|CustomPart|orig_id|part1|"
    assert synthetic_symbol_from_header(header) == "CustomPart"


def test_synthetic_symbol_from_header_uses_local_alias_map(tmp_path, monkeypatch):
    from tackle import utils

    alias_file = tmp_path / "synthetic-symbols.txt"
    alias_file.write_text(
        "vector_part=Vector_Part\n"
        "sp_vector_part_vector_part=Vector_Part\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("TACKLE_SYNTHETIC_SYMBOL_MAP", str(alias_file))
    monkeypatch.setattr(utils, "_SYNTHETIC_SYMBOL_MAP_CACHE", None)

    header = "meta:p:sp_vector_part_vector_part_geneid_990000001_taxon_32630_"
    assert utils.synthetic_symbol_from_header(header) == "Vector_Part"
