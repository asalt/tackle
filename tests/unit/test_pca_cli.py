from click.testing import CliRunner

from tackle.main import pca2


def test_single_pc_option_is_documented_as_repeatable_and_comma_delimited():
    result = CliRunner().invoke(pca2, ["--help"])

    assert result.exit_code == 0
    assert "--test-single-pc TEXT" in result.output
    assert "comma-delimited value such as 1,2" in result.output


def test_single_pc_option_requires_a_test_grouping_field():
    result = CliRunner().invoke(pca2, ["--test-single-pc", "1,2"])

    assert result.exit_code == 2
    assert "--test-single-pc require at least one --test-by" in result.output
