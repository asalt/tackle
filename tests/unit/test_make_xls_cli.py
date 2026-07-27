from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from click.testing import CliRunner

from tackle import exporter
from tackle.main import make_xls


def _invoke_make_xls(monkeypatch, args, *, obj=None):
    captured = {}

    def fake_build_export_xlsx(**kwargs):
        captured.update(kwargs)
        return kwargs["out_path"]

    monkeypatch.setattr(exporter, "build_export_xlsx", fake_build_export_xlsx)
    result = CliRunner().invoke(make_xls, args, obj=obj)
    return result, captured


def test_make_xls_defaults_to_summary_under_standalone_export_dir(
    monkeypatch, tmp_path
):
    result, captured = _invoke_make_xls(
        monkeypatch,
        ["--base-dir", str(tmp_path)],
    )

    assert result.exit_code == 0, result.output
    assert captured["base_dir"] == str(tmp_path)
    assert captured["out_path"] == str(tmp_path / "export" / "summary.xlsx")
    assert captured["pheno_df"] is None
    assert captured["meta"]["outpath"] == str(tmp_path)


def test_make_xls_name_and_outpath_are_resolved_independently(
    monkeypatch, tmp_path
):
    destination = tmp_path / "workbooks"
    result, captured = _invoke_make_xls(
        monkeypatch,
        [
            "--base-dir",
            str(tmp_path),
            "--name",
            "review.xlsx",
            "--outpath",
            str(destination),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["base_dir"] == str(tmp_path)
    assert captured["out_path"] == str(destination / "review.xlsx")


def test_make_xls_loaded_context_uses_analysis_export_dir(monkeypatch, tmp_path):
    pheno = pd.DataFrame({"group": ["A", "B"]}, index=["s1", "s2"])
    data_obj = SimpleNamespace(
        outpath=str(tmp_path / "analysis"),
        outpath_name="run1",
        col_metadata=pheno,
        taxon="human",
        non_zeros=1,
        batch_applied=False,
        batch_nonparametric=False,
    )
    result, captured = _invoke_make_xls(
        monkeypatch,
        [],
        obj={"data_obj": data_obj},
    )

    assert result.exit_code == 0, result.output
    assert captured["base_dir"] == data_obj.outpath
    assert captured["out_path"] == str(
        Path(data_obj.outpath) / "export" / "summary.xlsx"
    )
    assert captured["pheno_df"] is pheno
    assert captured["meta"]["analysis_name"] == "run1"


def test_make_xls_name_rejects_paths_and_non_xlsx_extensions(monkeypatch, tmp_path):
    nested_result, nested_capture = _invoke_make_xls(
        monkeypatch,
        [
            "--base-dir",
            str(tmp_path),
            "--name",
            "nested/review.xlsx",
        ],
    )
    extension_result, extension_capture = _invoke_make_xls(
        monkeypatch,
        [
            "--base-dir",
            str(tmp_path),
            "--name",
            "review.xls",
        ],
    )

    assert nested_result.exit_code == 2
    assert "use --outpath to select its directory" in nested_result.output
    assert not nested_capture
    assert extension_result.exit_code == 2
    assert "must end in .xlsx" in extension_result.output
    assert not extension_capture


def test_make_xls_help_describes_name_and_default_output_directory():
    result = CliRunner().invoke(make_xls, ["--help"])
    normalized_help = " ".join(result.output.split())

    assert result.exit_code == 0, result.output
    assert "-n, --name TEXT" in normalized_help
    assert "[default: summary.xlsx]" in normalized_help
    assert "--outpath DIRECTORY" in normalized_help
    assert "Defaults to <analysis outpath>/export" in normalized_help
