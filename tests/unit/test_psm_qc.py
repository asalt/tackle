from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from click.testing import CliRunner

from tackle.main import metrics
from tackle.psm_qc import run_psm_qc


def _write_qual(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _modern_rows() -> list[dict]:
    return [
        {
            "SpectrumFile": "run.raw",
            "Precursor.Id": "AAA2",
            "RTmin": 0.2,
            "Sequence": "AAA",
            "SequenceModi": "AAA",
            "Charge": 2,
            "MissedCleavages": 0,
            "q_value": 0.001,
            "PEP": 0.002,
            "DeltaMassPPM": -0.5,
            "PrecursorArea": 100,
            "PSM_UseFLAG": 1,
            "oriFLAG": 1,
        },
        {
            "SpectrumFile": "run.raw",
            "Precursor.Id": "AAA2",
            "RTmin": 0.2,
            "Sequence": "AAA",
            "SequenceModi": "AAA",
            "Charge": 2,
            "MissedCleavages": 0,
            "q_value": 0.001,
            "PEP": 0.002,
            "DeltaMassPPM": -0.5,
            "PrecursorArea": 100,
            "PSM_UseFLAG": 1,
            "oriFLAG": 0,
        },
        {
            "SpectrumFile": "run.raw",
            "Precursor.Id": "BBB3",
            "RTmin": 1.2,
            "Sequence": "BBB",
            "SequenceModi": "BBB",
            "Charge": 3,
            "MissedCleavages": 1,
            "q_value": 0.02,
            "PEP": 0.04,
            "DeltaMassPPM": 1.0,
            "PrecursorArea": 200,
            "PSM_UseFLAG": 0,
            "oriFLAG": 1,
        },
        {
            "SpectrumFile": "run.raw",
            "Precursor.Id": "CCC3",
            "RTmin": 1.8,
            "Sequence": "CCC",
            "SequenceModi": "CCC[Oxidation]",
            "Charge": 3,
            "MissedCleavages": 1,
            "q_value": 0.003,
            "PEP": 0.006,
            "DeltaMassPPM": 0.5,
            "PrecursorArea": 300,
            "PSM_UseFLAG": 1,
            "oriFLAG": 1,
        },
    ]


def test_psm_qc_scans_each_search_once_and_applies_provenance_flags(
    tmp_path: Path,
) -> None:
    pytest.importorskip("polars")
    source = tmp_path / "raw" / "123_1_4_labelnone_psms_QUAL.tsv"
    _write_qual(source, _modern_rows())
    config = OrderedDict(
        {
            "TMT_A": {"recno": "123", "runno": "1", "searchno": "4"},
            "TMT_B": {"recno": "123.0", "runno": "1", "searchno": "4"},
        }
    )

    artifacts = run_psm_qc(
        config=config,
        data_dir=tmp_path,
        output_stem=tmp_path / "out" / "study_psm_qc",
        file_fmts=(".png",),
        png_res=72,
    )

    summary = pd.read_csv(artifacts.summary_tsv, sep="\t")
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["config_sections"] == "TMT_A;TMT_B"
    assert row["source_file"] == str(source.resolve())
    assert row["psm_use_filter"] == "applied"
    assert row["ori_filter"] == "applied"
    assert row["identity_basis"] == "oriFLAG"
    assert row["rows_raw"] == 4
    assert row["rows_psm_use"] == 3
    assert row["rows_canonical"] == 2
    assert row["unique_sequences"] == 2
    assert row["unique_modified_sequences"] == 2
    assert row["precursor_area_median"] == 200

    ids = pd.read_csv(artifacts.ids_over_time_tsv, sep="\t")
    assert ids[["rt_minute", "ids"]].to_dict("records") == [
        {"rt_minute": 0, "ids": 1},
        {"rt_minute": 1, "ids": 1},
    ]
    charge = pd.read_csv(artifacts.charge_distribution_tsv, sep="\t")
    assert charge[["charge", "count", "fraction"]].to_dict("records") == [
        {"charge": 2, "count": 1, "fraction": 0.5},
        {"charge": 3, "count": 1, "fraction": 0.5},
    ]
    assert {path.name for path in artifacts.plot_paths} == {
        "study_psm_qc_ids_over_time.png",
        "study_psm_qc_composition.png",
    }
    assert all(path.exists() for path in artifacts.plot_paths)


def test_psm_qc_uses_scan_identity_when_oriflag_is_unavailable(
    tmp_path: Path,
) -> None:
    pytest.importorskip("polars")
    rows = [
        {
            "SpectrumFile": "legacy.raw",
            "FirstScan": 10,
            "RTmin": 1.1,
            "Sequence": "AAA",
            "SequenceModi": "AAA",
            "Charge": 2,
        },
        {
            "SpectrumFile": "legacy.raw",
            "FirstScan": 10,
            "RTmin": 1.1,
            "Sequence": "AAA",
            "SequenceModi": "AAA",
            "Charge": 2,
        },
        {
            "SpectrumFile": "legacy.raw",
            "FirstScan": 11,
            "RTmin": 1.8,
            "Sequence": "BBB",
            "SequenceModi": "BBB",
            "Charge": 3,
        },
    ]
    _write_qual(tmp_path / "200_1_1_labelnone_psms_QUAL.tsv", rows)

    artifacts = run_psm_qc(
        config={"Legacy": {"recno": 200, "runno": 1, "searchno": 1}},
        data_dir=tmp_path,
        output_stem=tmp_path / "legacy_psm_qc",
        file_fmts=(),
    )

    summary = pd.read_csv(artifacts.summary_tsv, sep="\t").iloc[0]
    assert summary["psm_use_filter"] == "not_available"
    assert summary["ori_filter"] == "not_available"
    assert summary["identity_basis"] == "SpectrumFile+FirstScan"
    assert summary["rows_raw"] == 3
    assert summary["rows_canonical"] == 2


def test_psm_qc_reports_missing_and_ambiguous_sources_without_guessing(
    tmp_path: Path,
) -> None:
    pytest.importorskip("polars")
    _write_qual(tmp_path / "a" / "301_1_1_labelnone_psms_QUAL.tsv", _modern_rows())
    _write_qual(tmp_path / "b" / "301_1_1_labelnone_psms_QUAL.tsv", _modern_rows())
    config = OrderedDict(
        {
            "Missing": {"recno": 300, "runno": 1, "searchno": 1},
            "Ambiguous": {"recno": 301, "runno": 1, "searchno": 1},
        }
    )

    with pytest.warns(UserWarning):
        artifacts = run_psm_qc(
            config=config,
            data_dir=tmp_path,
            output_stem=tmp_path / "status_psm_qc",
            file_fmts=(),
        )

    summary = pd.read_csv(artifacts.summary_tsv, sep="\t")
    assert summary.set_index("run_id")["status"].to_dict() == {
        "300_1_1": "missing",
        "301_1_1": "ambiguous",
    }
    assert summary["source_file"].isna().all()
    assert not artifacts.plot_paths


def test_metrics_cli_forwards_psm_qc_options(tmp_path: Path, monkeypatch) -> None:
    captured = {}

    def _fake_make_metrics(_data_obj, _file_fmts, **kwargs):
        captured.update(kwargs)
        return {}

    monkeypatch.setattr("tackle.main.make_metrics", _fake_make_metrics)
    data_obj = SimpleNamespace(data_dir=tmp_path)
    result = CliRunner().invoke(
        metrics,
        ["--psm-qc", "--psm-data-dir", str(tmp_path)],
        obj={"data_obj": data_obj, "file_fmts": (".png",), "png_res": 150},
    )

    assert result.exit_code == 0, result.output
    assert captured["psm_qc"] is True
    assert captured["psm_data_dir"] == tmp_path
    assert captured["png_res"] == 150


def test_metrics_full_remains_a_backward_compatible_psm_qc_alias(
    tmp_path: Path, monkeypatch
) -> None:
    captured = {}

    def _fake_make_metrics(_data_obj, _file_fmts, **kwargs):
        captured.update(kwargs)
        return {}

    monkeypatch.setattr("tackle.main.make_metrics", _fake_make_metrics)
    result = CliRunner().invoke(
        metrics,
        ["--full"],
        obj={
            "data_obj": SimpleNamespace(data_dir=tmp_path),
            "file_fmts": (".png",),
            "png_res": 300,
        },
    )

    assert result.exit_code == 0, result.output
    assert captured["psm_qc"] is True
