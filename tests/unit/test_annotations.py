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

