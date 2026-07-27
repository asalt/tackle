from types import SimpleNamespace

from click.testing import CliRunner

from tackle.main import make_deck, make_html


def test_make_html_plot_kind_help_includes_gene_level_plot_families():
    result = CliRunner().invoke(make_html, ["--help"])

    assert result.exit_code == 0
    assert "bar" in result.output
    assert "box" in result.output
    assert "violin" in result.output
    assert "--not-plot-kind" in result.output


def test_make_html_not_plot_kind_subtracts_from_default_and_explicit_sets(
    monkeypatch, tmp_path
):
    from tackle import html_overview

    captured = []

    def fake_build_html_overview(**kwargs):
        captured.append(kwargs["include_kinds"])
        return SimpleNamespace(
            out_html=tmp_path / "index.html",
            assets_dir=None,
            methods_md=None,
        )

    monkeypatch.setattr(
        html_overview,
        "build_html_overview",
        fake_build_html_overview,
    )
    runner = CliRunner()

    default_result = runner.invoke(
        make_html,
        [
            "--base-dir",
            str(tmp_path),
            "--outdir",
            str(tmp_path / "default-report"),
            "--not-plot-kind",
            "volcano",
        ],
    )
    assert default_result.exit_code == 0, default_result.output
    assert "volcano" not in captured[0]
    assert "pca" in captured[0]
    assert "correlation" in captured[0]

    explicit_result = runner.invoke(
        make_html,
        [
            "--base-dir",
            str(tmp_path),
            "--outdir",
            str(tmp_path / "explicit-report"),
            "--plot-kind",
            "pca",
            "--plot-kind",
            "volcano",
            "--not-plot-kind",
            "volcano",
        ],
    )
    assert explicit_result.exit_code == 0, explicit_result.output
    assert captured[1] == ("pca",)


def test_make_deck_plot_kind_help_includes_gene_level_plot_families():
    result = CliRunner().invoke(make_deck, ["--help"])

    assert result.exit_code == 0
    assert "bar" in result.output
    assert "box" in result.output
    assert "violin" in result.output
