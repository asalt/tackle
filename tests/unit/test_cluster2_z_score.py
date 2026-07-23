import inspect

from click.testing import CliRunner
import pandas as pd
import pytest

from tackle import clusterplot_dispatcher
from tackle.clusterplot_dispatcher import (
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
