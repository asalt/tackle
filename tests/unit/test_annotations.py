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
