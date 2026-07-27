import inspect

from click.testing import CliRunner

from tackle import clusterplot_dispatcher
from tackle.main import cluster2


def _invoke_with_captured_dispatch(monkeypatch, options, *, global_formats=(".pdf",)):
    run_signature = inspect.signature(clusterplot_dispatcher.run)
    captured = {}

    def fake_run(*args, **kwargs):
        bound = run_signature.bind(*args, **kwargs)
        captured.update(bound.arguments)

    monkeypatch.setattr(clusterplot_dispatcher, "run", fake_run)
    obj = {"file_fmts": tuple(global_formats)}
    result = CliRunner().invoke(cluster2, options, obj=obj)
    return result, captured, obj


def test_cluster2_file_format_help_describes_global_fallback():
    result = CliRunner().invoke(cluster2, ["--help"])
    normalized_help = " ".join(result.output.split())

    assert result.exit_code == 0, result.output
    assert "--file-format [.png|.pdf|.svg]" in normalized_help
    assert "Defaults to the global --file-format setting." in normalized_help


def test_cluster2_inherits_global_file_formats_without_mutating_context(monkeypatch):
    result, captured, obj = _invoke_with_captured_dispatch(
        monkeypatch,
        [],
        global_formats=(".pdf", ".svg"),
    )

    assert result.exit_code == 0, result.output
    assert captured["file_fmts"] == (".pdf", ".svg")
    assert obj["file_fmts"] == (".pdf", ".svg")


def test_cluster2_local_file_formats_override_global_and_are_repeatable(monkeypatch):
    result, captured, obj = _invoke_with_captured_dispatch(
        monkeypatch,
        [
            "--file-format",
            ".png",
            "--file-format",
            ".svg",
        ],
        global_formats=(".pdf",),
    )

    assert result.exit_code == 0, result.output
    assert captured["file_fmts"] == (".png", ".svg")
    assert obj["file_fmts"] == (".pdf",)


def test_cluster2_warns_when_figsize_overrides_individual_dimensions(monkeypatch):
    result, captured, _obj = _invoke_with_captured_dispatch(
        monkeypatch,
        [
            "--figsize",
            "8",
            "9",
            "--figwidth",
            "12",
            "--figheight",
            "13",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (
        "Warning: --figsize was provided with --figwidth or --figheight; "
        "using the --figsize values."
    ) in result.output
    assert captured["figsize"] == (8.0, 9.0)
    assert captured["figwidth"] == 12.0
    assert captured["figheight"] == 13.0


def test_cluster2_help_explains_description_and_dimension_behavior():
    result = CliRunner().invoke(cluster2, ["--help"])
    normalized_help = " ".join(result.output.split())

    assert result.exit_code == 0, result.output
    assert "Keep UniProt-style attributes" in normalized_help
    assert "if present." in normalized_help
    assert "Set figure width in inches and auto-size height." in normalized_help
    assert "Set figure height in inches and auto-size width." in normalized_help
    assert "Overrides --figwidth and --figheight" in normalized_help
