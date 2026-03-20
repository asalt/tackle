import pandas as pd

from tackle.clusterplot_dispatcher import (
    _resolve_annotation_filter_gene_ids,
    _resolve_annotation_hits,
    _series_matches_includes,
)


def test_series_matches_includes_numeric_metadata_with_string_values():
    series = pd.Series([0.0, 1.0, 2.0, 1.0], index=["s1", "s2", "s3", "s4"])

    mask = _series_matches_includes(series, ["0", "2"])

    assert mask.to_dict() == {"s1": True, "s2": False, "s3": True, "s4": False}


def test_series_matches_includes_object_numeric_metadata_with_string_values():
    series = pd.Series([0.0, 1.0, 2.0, ""], index=["s1", "s2", "s3", "s4"], dtype=object)

    mask = _series_matches_includes(series, ["1"])

    assert mask.to_dict() == {"s1": False, "s2": True, "s3": False, "s4": False}


def test_resolve_annotation_hits_prefers_direct_hits():
    x_gene_ids = ["101", "102", "103"]
    annot_gene_ids = ["999", "102"]

    hits, strategy = _resolve_annotation_hits(
        x_gene_ids=x_gene_ids,
        annotation_gene_ids=annot_gene_ids,
        hgene_mapper=None,
    )

    assert hits == ["102"]
    assert strategy == "direct"


def test_resolve_annotation_hits_recovers_via_human_remap():
    class DummyMapper:
        def map_to_taxon(self, gids, target_taxon="9606"):
            mapping = {
                ("m1", "9606"): "h1",
                ("m2", "9606"): "h2",
                ("h1", "9606"): "h1",
                ("h2", "9606"): "h2",
                ("h1", "10090"): "m1",
                ("h2", "10090"): "m2",
            }
            return {gid: mapping.get((gid, str(target_taxon))) for gid in gids}

    hits, strategy = _resolve_annotation_hits(
        x_gene_ids=["m1", "m2"],
        annotation_gene_ids=["h1"],
        hgene_mapper=DummyMapper(),
    )

    assert hits == ["m1"]
    assert strategy == "x_to_human"


def test_resolve_annotation_hits_returns_none_when_no_mapping():
    class EmptyMapper:
        def map_to_taxon(self, gids, target_taxon="9606"):
            return {gid: None for gid in gids}

    hits, strategy = _resolve_annotation_hits(
        x_gene_ids=["x1", "x2"],
        annotation_gene_ids=["a1", "a2"],
        hgene_mapper=EmptyMapper(),
    )

    assert hits == []
    assert strategy == "none"


def test_resolve_annotation_filter_gene_ids_unions_hits_and_preserves_order():
    class DummyAnnot:
        def __init__(self, mapping):
            self.mapping = mapping

        def get_annot(self, annotation):
            return pd.DataFrame({"GeneID": self.mapping.get(annotation, [])})

    hits = _resolve_annotation_filter_gene_ids(
        x_gene_ids=["1", "2", "3", "4"],
        annotation_filters=["A", "B"],
        annotator=DummyAnnot({"A": ["2", "4"], "B": ["1"]}),
        hgene_mapper=None,
    )

    assert hits == ["1", "2", "4"]


def test_resolve_annotation_filter_gene_ids_recovers_via_homologene_remap():
    class DummyAnnot:
        def get_annot(self, annotation):
            return pd.DataFrame({"GeneID": ["h1"]})

    class DummyMapper:
        def map_to_taxon(self, gids, target_taxon="9606"):
            mapping = {
                ("m1", "9606"): "h1",
                ("m2", "9606"): "h2",
                ("h1", "9606"): "h1",
                ("h2", "9606"): "h2",
                ("h1", "10090"): "m1",
                ("h2", "10090"): "m2",
            }
            return {gid: mapping.get((gid, str(target_taxon))) for gid in gids}

    hits = _resolve_annotation_filter_gene_ids(
        x_gene_ids=["m1", "m2"],
        annotation_filters=["SECRETED"],
        annotator=DummyAnnot(),
        hgene_mapper=DummyMapper(),
    )

    assert hits == ["m1"]


def test_resolve_annotation_filter_gene_ids_ignores_all_sentinel():
    class DummyAnnot:
        def get_annot(self, annotation):
            return pd.DataFrame({"GeneID": ["2"]})

    hits = _resolve_annotation_filter_gene_ids(
        x_gene_ids=["1", "2"],
        annotation_filters=["_all"],
        annotator=DummyAnnot(),
        hgene_mapper=None,
    )

    assert hits == ["1", "2"]


def test_resolve_annotation_filter_gene_ids_returns_empty_when_no_hits():
    class DummyAnnot:
        def get_annot(self, annotation):
            return pd.DataFrame({"GeneID": ["9"]})

    hits = _resolve_annotation_filter_gene_ids(
        x_gene_ids=["1", "2"],
        annotation_filters=["A"],
        annotator=DummyAnnot(),
        hgene_mapper=None,
    )

    assert hits == []
