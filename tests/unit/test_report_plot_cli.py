from click.testing import CliRunner

from tackle.main import make_deck, make_html


def test_make_html_plot_kind_help_includes_gene_level_plot_families():
    result = CliRunner().invoke(make_html, ["--help"])

    assert result.exit_code == 0
    assert "bar" in result.output
    assert "box" in result.output
    assert "violin" in result.output


def test_make_deck_plot_kind_help_includes_gene_level_plot_families():
    result = CliRunner().invoke(make_deck, ["--help"])

    assert result.exit_code == 0
    assert "bar" in result.output
    assert "box" in result.output
    assert "violin" in result.output
