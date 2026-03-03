import json
from pathlib import Path

from click.testing import CliRunner


def test_cli1_init_writes_analysis_files(tmp_path):
    from tackle.cli1.app import cli

    conf_path = tmp_path / "Example.conf"
    conf_path.write_text("dummy=1\n", encoding="utf-8")

    analysis_dir = tmp_path / "analysis"
    data_dir = tmp_path / "data"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--analysis-dir",
            str(analysis_dir),
            "--conf",
            str(conf_path),
            "--name",
            "Run1",
            "init",
            "--data-dir",
            str(data_dir),
        ],
    )
    assert result.exit_code == 0, result.output

    assert (analysis_dir / "analysis.json").exists()
    assert (analysis_dir / "manifest.sqlite").exists()
    assert (analysis_dir / "manifest.json").exists()
    assert (analysis_dir / "inputs" / "Example.conf").exists()
    assert (analysis_dir / "matrices").exists()

    payload = json.loads((analysis_dir / "analysis.json").read_text(encoding="utf-8"))
    assert payload["analysis_name"] == "Example"
    assert payload["run_name"] == "Run1"
    assert Path(payload["analysis_dir"]) == analysis_dir.resolve()

