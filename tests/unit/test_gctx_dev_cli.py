from __future__ import annotations

import json

import pandas as pd
from click.testing import CliRunner

import tackle.main as tackle_main
from tackle.gct_io import write_gctx


def _write_demo_gctx(tmp_path):
    matrix = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=["gene_1", "gene_2"],
        columns=["sample_a", "sample_b"],
    )
    return write_gctx(matrix, tmp_path / "demo.gctx")


def test_dev_inspect_gctx_runs_through_standalone_dispatch(tmp_path, monkeypatch):
    target = _write_demo_gctx(tmp_path)
    argv = ["tackle", "dev", "inspect-gctx", str(target)]
    monkeypatch.setattr(tackle_main.sys, "argv", argv)

    result = CliRunner().invoke(tackle_main.main, argv[1:])

    assert result.exit_code == 0, result.output
    assert "GCTX inspection" in result.output
    assert "Logical matrix: 2 rows x 2 columns" in result.output
    assert "Tackle provenance" in result.output
    assert "Content hash algorithm:" in result.output


def test_dev_inspect_gctx_json_and_spec(tmp_path):
    target = _write_demo_gctx(tmp_path)
    command = tackle_main.main.commands["dev"]

    result = CliRunner().invoke(
        command,
        ["inspect-gctx", str(target), "--json", "--spec"],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["inspection"]["logical_shape"] == {"rows": 2, "columns": 2}
    assert payload["inspection"]["provenance"]["available"] is True
    assert payload["spec"]["content_hash"]["fallback_algorithm"] == "blake2b-256"


def test_dev_inspect_gctx_can_print_spec_without_target():
    command = tackle_main.main.commands["dev"]

    result = CliRunner().invoke(command, ["inspect-gctx", "--spec"])

    assert result.exit_code == 0, result.output
    assert "Tackle GCTX provenance specification" in result.output
    assert "fall back to blake2b-256" in result.output
