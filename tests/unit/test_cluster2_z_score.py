import inspect

from click.testing import CliRunner
import pandas as pd
import pytest

from tackle import clusterplot_dispatcher
from tackle.clusterplot_dispatcher import (
    _cluster2_genefile_name,
    _select_cluster2_genes,
    normalize_cluster2_standard_scale,
    normalize_cluster2_z_score,
    standard_scale_cluster2_matrix,
)
from tackle.main import cluster2


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, None),
        ("None", None),
        ("row", "row"),
        ("column", "column"),
        ("col", "column"),
        ("0", "row"),
        ("1", "column"),
    ],
)
def test_normalize_cluster2_z_score(value, expected):
    assert normalize_cluster2_z_score(value) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, None),
        ("None", None),
        ("row", "row"),
        ("column", "column"),
        ("col", "column"),
        ("0", "column"),
        ("1", "row"),
    ],
)
def test_normalize_cluster2_standard_scale(value, expected):
    assert normalize_cluster2_standard_scale(value) == expected


def test_standard_scale_cluster2_matrix_uses_named_axes():
    matrix = pd.DataFrame({"s1": [1.0, 3.0], "s2": [2.0, 6.0]})

    row_scaled = standard_scale_cluster2_matrix(matrix, "row")
    column_scaled = standard_scale_cluster2_matrix(matrix, "column")

    assert row_scaled.sum(axis=1).tolist() == pytest.approx([1.0, 1.0])
    assert column_scaled.sum(axis=0).tolist() == pytest.approx([1.0, 1.0])


def test_select_cluster2_genes_preserves_genefile_order_and_forces_missing_rows():
    matrix = pd.DataFrame(
        {"s1": [1.0, 2.0], "s2": [3.0, 4.0]},
        index=pd.Index(["g1", "g2"], name="GeneID"),
    )
    genes = ["g3", "g1", "g2", "g3"]

    selected = _select_cluster2_genes(matrix, genes, force_plot_genes=True)

    assert selected.index.tolist() == ["g3", "g1", "g2"]
    assert selected.loc["g3"].isna().all()
    pd.testing.assert_series_equal(selected.loc["g1"], matrix.loc["g1"])


def test_select_cluster2_genes_omits_missing_rows_without_force():
    matrix = pd.DataFrame(
        {"s1": [1.0, 2.0]},
        index=pd.Index(["g1", "g2"], name="GeneID"),
    )

    selected = _select_cluster2_genes(
        matrix, ["g3", "g2", "g1"], force_plot_genes=False
    )

    assert selected.index.tolist() == ["g2", "g1"]


def test_force_plot_genes_does_not_restore_histone_dedupe_exclusions():
    matrix = pd.DataFrame(
        {"s1": [1.0]},
        index=pd.Index(["retained_histone"], name="GeneID"),
    )

    selected = _select_cluster2_genes(
        matrix,
        ["collapsed_histone", "unmeasured_gene", "retained_histone"],
        force_plot_genes=True,
        histone_dedupe_excluded_geneids={"collapsed_histone"},
    )

    assert selected.index.tolist() == ["unmeasured_gene", "retained_histone"]
    assert selected.loc["unmeasured_gene"].isna().all()


def test_cluster2_genefile_name_uses_only_the_filename():
    name = _cluster2_genefile_name("/analysis/dv-parent/gene lists/my genes.tsv")

    assert name == "my_genes"
    assert "/" not in name
    assert "dv-parent" not in name


def test_cluster2_geneset_name_removes_path_separators_without_dv_encoding():
    name = _cluster2_genefile_name("GO/pathway\\subset", from_path=False)

    assert name == "GOpathwaysubset"


def test_cluster2_z_score_help_uses_named_default_and_documents_aliases():
    result = CliRunner().invoke(cluster2, ["--help"])
    normalized_help = " ".join(result.output.split())

    assert result.exit_code == 0, result.output
    assert "[None|row|column|col|0|1]" in result.output
    assert "[default: row]" in result.output
    assert "col=column" in normalized_help
    assert "0=row, 1=column" in normalized_help


def test_cluster2_standard_scale_help_retains_opt_in_default_and_aliases():
    result = CliRunner().invoke(cluster2, ["--help"])
    normalized_help = " ".join(result.output.split())

    assert result.exit_code == 0, result.output
    assert "--standard-scale [None|row|column|col|0|1]" in normalized_help
    assert "0=column, 1=row" in normalized_help
    standard_scale_help = normalized_help.split("--standard-scale", 1)[1]
    assert "[default: None]" in standard_scale_help.split("--show-missing-values", 1)[0]


@pytest.mark.parametrize(
    ("option", "expected"),
    [
        ("--keep-description-attributes", True),
        ("--no-keep-description-attributes", False),
    ],
)
def test_cluster2_forwards_keep_description_attributes_by_keyword(
    monkeypatch, option, expected
):
    run_signature = inspect.signature(clusterplot_dispatcher.run)
    captured = {}

    def fake_run(*args, **kwargs):
        bound = run_signature.bind(*args, **kwargs)
        captured.update(bound.arguments)

    monkeypatch.setattr(clusterplot_dispatcher, "run", fake_run)

    result = CliRunner().invoke(cluster2, [option])

    assert result.exit_code == 0, result.output
    assert captured["keep_description_attributes"] is expected
