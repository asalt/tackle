import numpy as np
import pandas as pd
from types import SimpleNamespace

from tackle.correlationplot import (
    build_correlation_cdesc,
    build_correlation_rdesc,
    correlation_contract,
    compute_pairwise_counts,
    compute_sample_metric,
    exclude_correlation_samples,
    flatten_metadata_fields,
    normalize_linkage,
    normalize_metric,
    prepare_logged_matrix,
    resolve_linkage,
    select_plot_metadata,
    write_correlation_input_gctx,
    write_correlation_sample_gctx,
    zscore_feature_rows,
)
from tackle.gct_io import read_gctx


class _FakeData:
    def __init__(self):
        self.areas_log_shifted = pd.DataFrame(
            {"A": [1.0, 2.0], "B": [3.0, 4.0]},
            index=pd.Index([101, 202], name="GeneID"),
        )
        self.mask = pd.DataFrame(
            {"A": [False, True], "B": [False, False]},
            index=pd.Index([101, 202], name="GeneID"),
        )


def test_prepare_logged_matrix_applies_mask_and_normalizes_identifiers():
    matrix = prepare_logged_matrix(_FakeData())

    assert list(matrix.index) == ["101", "202"]
    assert list(matrix.columns) == ["A", "B"]
    assert np.isnan(matrix.loc["202", "A"])
    assert matrix.loc["202", "B"] == 4.0


def test_exclude_correlation_samples_preserves_order_and_records_unknown_names():
    matrix = pd.DataFrame(
        {"A": [1.0], "B": [2.0], "C": [3.0], "D": [4.0]},
        index=["g1"],
    )

    selected, selection = exclude_correlation_samples(
        matrix,
        ("C", "missing", "C"),
    )

    assert list(selected.columns) == ["A", "B", "D"]
    assert selection == {
        "requested": ["C", "missing"],
        "applied": ["C"],
        "unknown": ["missing"],
    }


def test_zscore_feature_rows_fills_below_minimum_then_restores_missing_values():
    matrix = pd.DataFrame(
        {"A": [1.0, 5.0], "B": [np.nan, 5.0], "C": [3.0, 5.0]},
        index=["variable", "constant"],
    )

    result = zscore_feature_rows(matrix)

    observed_sd = np.std([1.0, 3.0], ddof=1)
    filled = np.array([1.0, 1.0 - observed_sd, 3.0])
    expected = (filled - filled.mean()) / filled.std(ddof=1)
    assert np.isclose(result.loc["variable", "A"], expected[0])
    assert np.isnan(result.loc["variable", "B"])
    assert np.isclose(result.loc["variable", "C"], expected[2])
    assert np.allclose(result.loc["constant"], 0.0)


def test_zscore_feature_rows_keeps_sparse_constant_detections_positive():
    matrix = pd.DataFrame(
        {
            "A": [np.nan, np.nan],
            "B": [10.0, np.nan],
            "C": [10.0, 20.0],
            "D": [np.nan, np.nan],
        },
        index=["two_detections", "one_detection"],
    )

    result = zscore_feature_rows(matrix)

    assert np.isnan(result.loc["two_detections", "A"])
    assert np.isnan(result.loc["two_detections", "D"])
    assert np.allclose(
        result.loc["two_detections", ["B", "C"]],
        [np.sqrt(3) / 2, np.sqrt(3) / 2],
    )
    assert np.isclose(result.loc["one_detection", "C"], 1.5)
    assert result.loc["one_detection", ["A", "B", "D"]].isna().all()


def test_zscore_feature_rows_can_scale_within_multiple_metadata_fields():
    matrix = pd.DataFrame(
        {"A": [1.0], "B": [3.0], "C": [10.0], "D": [14.0]},
        index=["g1"],
    )
    metadata = pd.DataFrame(
        {
            "group": ["x", "x", "y", "y"],
            "batch": [1, 1, 1, 1],
        },
        index=["A", "B", "C", "D"],
    )

    result = zscore_feature_rows(matrix, metadata=metadata, by=("group:batch",))

    assert np.allclose(
        result.loc["g1"].to_numpy(),
        [-1 / np.sqrt(2), 1 / np.sqrt(2), -1 / np.sqrt(2), 1 / np.sqrt(2)],
    )


def test_compute_sample_metric_pearson_and_chord_distance():
    matrix = pd.DataFrame(
        {
            "A": [1.0, 2.0, 3.0],
            "B": [2.0, 4.0, 6.0],
            "C": [3.0, 2.0, 1.0],
        }
    )

    result = compute_sample_metric(matrix, "pearson")

    assert np.isclose(result.metric_matrix.loc["A", "B"], 1.0)
    assert np.isclose(result.metric_matrix.loc["A", "C"], -1.0)
    assert np.isclose(result.clustering_distance_matrix.loc["A", "C"], 2.0)
    assert np.isclose(result.clustering_distance_matrix.loc["A", "B"], 0.0)
    assert result.overlap_counts.loc["A", "C"] == 3


def test_compute_sample_metric_spearman_displays_r_not_distance():
    matrix = pd.DataFrame(
        {
            "A": [1.0, 2.0, 3.0],
            "B": [2.0, 4.0, 6.0],
            "C": [3.0, 2.0, 1.0],
        }
    )

    result = compute_sample_metric(matrix, "spearman")

    assert np.isclose(result.metric_matrix.loc["A", "C"], -1.0)
    assert np.isclose(result.clustering_distance_matrix.loc["A", "C"], 2.0)
    assert not result.metric_matrix.equals(result.clustering_distance_matrix)


def test_correlation_contract_resolves_metric_aliases_and_automatic_linkages():
    contract = correlation_contract()

    assert contract["default_linkage"] == "auto"
    assert contract["linkage_aliases"] == {"weighted": "mcquitty"}
    assert normalize_metric("euclidean") == "l2"
    assert normalize_metric("manhattan") == "l1"
    assert resolve_linkage("l2", "auto") == "ward.D2"
    assert resolve_linkage("pearson", "auto") == "ward.D2"
    assert resolve_linkage("l1", "auto") == "average"
    assert normalize_linkage("ward.D2") == "ward.D2"
    assert "weighted" in contract["linkages"]


def test_pairwise_counts_can_remain_pre_transform_for_all_metrics():
    source = pd.DataFrame(
        {
            "A": [1.0, 5.0, np.nan],
            "B": [1.0, 5.0, 2.0],
            "C": [1.0, 5.0, 3.0],
        },
        index=["constant", "also_constant", "partial"],
    )
    counts = compute_pairwise_counts(source)
    transformed = zscore_feature_rows(source)

    result = compute_sample_metric(
        transformed,
        "pearson",
        overlap_counts=counts,
    )

    assert result.overlap_counts.equals(counts)
    assert counts.loc["A", "B"] == 2
    assert counts.loc["B", "C"] == 3


def test_compute_sample_metric_euclidean_uses_pairwise_complete_rms():
    matrix = pd.DataFrame(
        {
            "A": [1.0, 4.0, np.nan],
            "B": [3.0, 8.0, 1000.0],
            "C": [1.0, np.nan, 7.0],
        },
        index=["g1", "g2", "g3"],
    )

    result = compute_sample_metric(matrix, "euclidean")

    # A/B share two features: sqrt(((1-3)^2 + (4-8)^2) / 2).
    assert np.isclose(result.metric_matrix.loc["A", "B"], np.sqrt(10.0))
    # A/C share only g1; B/C share g1 and g3. No missing-value fill leaks in.
    assert np.isclose(result.metric_matrix.loc["A", "C"], 0.0)
    assert np.isclose(
        result.metric_matrix.loc["B", "C"],
        np.sqrt(((3.0 - 1.0) ** 2 + (1000.0 - 7.0) ** 2) / 2.0),
    )
    assert result.metric_matrix.equals(result.clustering_distance_matrix)
    assert result.overlap_counts.loc["A", "B"] == 2
    assert result.overlap_counts.loc["A", "C"] == 1


def test_compute_sample_metric_l1_uses_pairwise_complete_mean_absolute_distance():
    matrix = pd.DataFrame(
        {
            "A": [1.0, 5.0, np.nan],
            "B": [4.0, 9.0, 1000.0],
            "C": [1.0, np.nan, 7.0],
        },
        index=["g1", "g2", "g3"],
    )

    result = compute_sample_metric(matrix, "manhattan")

    assert result.metric == "l1"
    assert np.isclose(result.metric_matrix.loc["A", "B"], (3.0 + 4.0) / 2.0)
    assert np.isclose(result.metric_matrix.loc["A", "C"], 0.0)
    assert np.isclose(
        result.metric_matrix.loc["B", "C"],
        (3.0 + 993.0) / 2.0,
    )
    assert result.metric_matrix.equals(result.clustering_distance_matrix)


def test_compute_sample_metric_euclidean_leaves_no_overlap_undefined():
    matrix = pd.DataFrame(
        {"A": [1.0, np.nan], "B": [np.nan, 2.0]},
        index=["g1", "g2"],
    )

    result = compute_sample_metric(matrix, "euclidean")

    assert np.isnan(result.metric_matrix.loc["A", "B"])
    assert np.isnan(result.clustering_distance_matrix.loc["A", "B"])
    assert result.overlap_counts.loc["A", "B"] == 0


def test_select_plot_metadata_preserves_repeatable_cut_fields():
    metadata = pd.DataFrame(
        {
            "recno": [1, 2, 3],
            "group": ["A", "A", "B"],
            "batch": ["x", "y", "x"],
            "sample_id": ["s1", "s2", "s3"],
        },
        index=["S1", "S2", "S3"],
    )

    selected, cut_fields = select_plot_metadata(
        metadata,
        ["S1", "S2", "S3"],
        exclude=("batch",),
        cut_by=("group:batch",),
    )

    assert cut_fields == ["group", "batch"]
    assert list(selected.columns) == ["group", "batch"]
    assert flatten_metadata_fields(("group:batch", "group")) == ["group", "batch"]


def test_build_correlation_descriptors_keep_all_sample_and_mapper_metadata(monkeypatch):
    from tackle import containers

    gene_mapper = SimpleNamespace(
        df=pd.DataFrame(
            {
                "GeneSymbol": ["mapped-a", "mapped-b"],
                "GeneDescription": ["desc-a", "desc-b"],
                "TaxonID": [9606, 9606],
            },
            index=pd.Index(["101", "202"], name="GeneID"),
        )
    )
    annotation_mapper = SimpleNamespace(
        df=pd.DataFrame(
            {
                "GeneID": [101, 202],
                "MitoCarta": ["Y", ""],
                "GeneSymbol": ["annotation-a", "annotation-b"],
            }
        )
    )
    monkeypatch.setattr(containers, "get_gene_mapper", lambda: gene_mapper)
    monkeypatch.setattr(containers, "get_annotation_mapper", lambda: annotation_mapper)

    data_obj = SimpleNamespace(
        gid_symbol={101: "carried-a", 202: "carried-b"},
        gid_funcat_mapping={101: "TF", 202: "PM"},
        col_metadata=pd.DataFrame(
            {
                "group": ["x", "y"],
                "recno": [11, 12],
                "continuous": [0.25, 0.75],
            },
            index=["A", "B"],
        ),
    )

    rdesc = build_correlation_rdesc(data_obj, [101, 202])
    cdesc = build_correlation_cdesc(data_obj, ["B", "A"])

    assert rdesc.loc["101", "GeneID"] == "101"
    assert rdesc.loc["101", "GeneSymbol"] == "carried-a"
    assert rdesc.loc["202", "GeneDescription"] == "desc-b"
    assert rdesc.loc["101", "MitoCarta"] == "Y"
    assert rdesc.loc["202", "FunCats"] == "PM"
    assert list(cdesc.index) == ["B", "A"]
    assert list(cdesc.columns) == ["group", "recno", "continuous", "id"]


def test_write_correlation_input_gctx_passes_exact_matrix_and_descriptors(
    monkeypatch, tmp_path
):
    import tackle.correlationplot as correlationplot

    matrix = pd.DataFrame(
        {"A": [1.25, np.nan], "B": [-0.5, 2.75]},
        index=pd.Index(["101", "202"], name="GeneID"),
    )
    data_obj = SimpleNamespace(
        col_metadata=pd.DataFrame({"group": ["x", "y"]}, index=["A", "B"]),
    )
    rdesc = pd.DataFrame(
        {"GeneID": ["101", "202"], "GeneSymbol": ["A1", "B2"], "id": ["101", "202"]},
        index=["101", "202"],
    )
    captured = {}

    monkeypatch.setattr(correlationplot, "build_correlation_rdesc", lambda *_: rdesc)

    def fake_write_gctx(matrix_arg, path_arg, **kwargs):
        captured["matrix"] = matrix_arg
        captured["path"] = path_arg
        captured.update(kwargs)
        return path_arg

    monkeypatch.setattr(correlationplot, "write_gctx", fake_write_gctx)
    output = write_correlation_input_gctx(data_obj, matrix, tmp_path / "input.gctx")

    assert output == tmp_path / "input.gctx"
    assert captured["matrix"].equals(matrix)
    assert list(captured["col_metadata"].columns) == ["group", "id"]
    assert captured["row_metadata"].equals(rdesc)
    assert captured["matrix_dtype"] == "float64"
    assert captured["content_addressed"] is True


def test_write_correlation_sample_gctx_embeds_original_metadata_on_both_axes(
    tmp_path,
):
    samples = ["A", "B"]
    matrix = pd.DataFrame(
        [[1.0, 0.25], [0.25, 1.0]],
        index=samples,
        columns=samples,
    )
    data_obj = SimpleNamespace(
        col_metadata=pd.DataFrame(
            {"group": ["x", "y"], "batch": [1, 2]}, index=samples
        )
    )

    output = write_correlation_sample_gctx(
        data_obj, matrix, tmp_path / "raw_metric.gctx"
    )
    parsed = read_gctx(output)

    assert parsed.row_metadata.equals(parsed.col_metadata)
    assert list(parsed.row_metadata.columns) == ["group", "batch"]
    assert parsed.row_metadata.loc["A", "group"] == "x"


def test_correlation_cli_exposes_metric_linkage_and_filter_contract():
    from tackle.main import main

    command = main.commands["correlation"]
    options = {flag for param in command.params for flag in param.opts}
    params = {param.name: param for param in command.params}

    assert {
        "--metric",
        "--linkage",
        "--sample-exclude",
        "--legend-include",
        "--legend-exclude",
        "--annotate",
        "--file-format",
    } <= options
    assert "--display" not in options
    assert "--method" not in options
    assert "--annotate-cells" not in options
    assert "--correlation-format" not in options
    assert tuple(params["metric"].type.choices) == (
        "l2",
        "l1",
        "pearson",
        "spearman",
        "euclidean",
        "manhattan",
    )
    assert tuple(params["linkage"].type.choices) == (
        "auto",
        "single",
        "complete",
        "average",
        "weighted",
        "centroid",
        "median",
        "ward.D2",
    )
    assert params["metric"].default == "l2"
    assert params["linkage"].default == "auto"
    assert params["sample_exclude"].multiple is True
    assert params["annotate"].default is False
