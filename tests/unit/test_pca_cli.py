from click.testing import CliRunner

from tackle.main import pca2


def test_single_pc_option_is_documented_as_repeatable_and_comma_delimited():
    result = CliRunner().invoke(pca2, ["--help"])

    assert result.exit_code == 0
    assert "--test-single-pc TEXT" in result.output
    assert "comma-delimited value such as 1,2" in result.output
    assert "starts at 6 x 7" in result.output
    assert "statistical captions" in " ".join(result.output.split())
    assert "--top-variable-pct" in result.output
    assert "--annotate-fontsize FLOAT RANGE" in result.output
    assert "Font size for sample labels drawn by --annotate" in " ".join(
        result.output.split()
    )
    assert "[default: 2.5; x>0]" in result.output
    assert "after --fillna and before centering/scaling" in " ".join(
        result.output.split()
    )
    assert "100 retains all features" in " ".join(result.output.split())


def test_single_pc_option_requires_a_test_grouping_field():
    result = CliRunner().invoke(pca2, ["--test-single-pc", "1,2"])

    assert result.exit_code == 2
    assert "--test-single-pc require at least one --test-by" in result.output


def test_top_variable_percentage_must_be_greater_than_zero_and_at_most_100():
    zero = CliRunner().invoke(pca2, ["--top-variable-pct", "0"])
    too_large = CliRunner().invoke(pca2, ["--top-variable-percent", "100.1"])

    assert zero.exit_code == 2
    assert "0<x<=100" in zero.output
    assert too_large.exit_code == 2
    assert "0<x<=100" in too_large.output


def test_annotate_fontsize_must_be_positive():
    result = CliRunner().invoke(pca2, ["--annotate-fontsize", "0"])

    assert result.exit_code == 2
    assert "x>0" in result.output
