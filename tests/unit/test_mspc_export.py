from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner
import numpy as np
import pandas as pd
import pytest

from tackle import containers
from tackle.containers import Data
from tackle.main import export as export_command


def test_mspc_export_source_uses_every_loaded_geneid_and_returns_a_copy():
    data = Data.__new__(Data)
    data.data = pd.DataFrame(
        {
            "GeneID": ["g1", "g1", "g2", "g2"],
            "Metric": ["GeneSymbol", "area", "GeneSymbol", "area"],
            "S1": ["A", 10.0, "B", 0.0],
        }
    )
    data.df_filtered = data.data.loc[data.data["GeneID"] == "g1"].set_index(["GeneID", "Metric"])

    source = data._mspc_export_source()

    assert list(source.index.get_level_values("GeneID").unique()) == ["g1", "g2"]
    source.loc[("g1", "area"), "S1"] = 99.0
    original = data.data.loc[
        (data.data["GeneID"] == "g1") & (data.data["Metric"] == "area"),
        "S1",
    ].iloc[0]
    assert original == 10.0


@pytest.fixture
def export_data(tmp_path, monkeypatch):
    data = Data.__new__(Data)
    raw = pd.DataFrame(
        {"S1": [10.0, 0.0, 0.01], "S2": [1000.0, np.nan, 0.02]},
        index=pd.Index(["g1", "g2", "g3"], name="GeneID"),
    )
    frames = []
    for metric, values in (
        ("iBAQ_dstrAdj", raw),
        ("area", raw.div(pd.Series({"S1": 10.0, "S2": 100.0}))),
        (
            "GeneSymbol",
            pd.DataFrame({"S1": ["A", "B", "C"], "S2": ["A", "B", "C"]}, index=raw.index),
        ),
        ("TaxonID", pd.DataFrame("9606", index=raw.index, columns=raw.columns)),
    ):
        frames.append(values.assign(Metric=metric).reset_index())
    data.data = pd.concat(frames, ignore_index=True)
    data.df_filtered = data.data[data.data.GeneID != "g3"].set_index(["GeneID", "Metric"])
    data._areas = data.df_filtered.xs("area", level="Metric").astype(float)
    data._areas_log = pd.DataFrame(42.0, index=data._areas.index, columns=raw.columns)
    data._mask = data._areas.replace(0, np.nan).isna()
    data._gid_symbol = {"g1": "A", "g2": "B", "g3": "C"}
    data.col_metadata = pd.DataFrame(
        {"recno": [1, 2], "runno": [1, 1], "searchno": [1, 1], "label": [0, 0]},
        index=raw.columns,
    )
    for flag in ("ifot", "ifot_ki", "ifot_tf", "trim_mean", "quantile75", "quantile90"):
        setattr(data, flag, False)
    data.median = True
    data.genefile_norm = None
    data.impute_missing_values = True
    data.fill_na_zero = False
    data.annotations = False
    data.outpath = str(tmp_path)
    data.outpath_name = "demo"
    data.taxon = "all"
    data.non_zeros = 1
    data.colors_only = False
    data.batch_applied = "batch"
    data.batch_nonparametric = False
    mapper = SimpleNamespace(
        symbol=data._gid_symbol,
        description={},
        df=pd.DataFrame({"median_isoform_mass": [10, 20, 30]}, index=raw.index),
    )
    monkeypatch.setattr(containers, "get_gene_mapper", lambda: mapper)
    return data


@pytest.mark.parametrize("fill_na_zero", [False, True])
def test_unnormalized_area_uses_raw_filtered_values_without_mutating_analysis(
    export_data, fill_na_zero
):
    data = export_data
    data.fill_na_zero = fill_na_zero
    before = [frame.copy(deep=True) for frame in (data.data, data.df_filtered, data.areas_log)]

    path = Path(data._export_unnormalized_area())
    table = pd.read_table(path, index_col="GeneID")

    assert path.parent == Path(data.outpath) / "export" / "data_area"
    assert path.name == "demo_data_area_nz1_none_unnormalized_linear.tsv"
    assert table.index.tolist() == ["g1", "g2"]
    assert table.columns.tolist() == ["GeneSymbol", "S1", "S2"]
    expected = data.df_filtered.xs("iBAQ_dstrAdj", level="Metric").astype(float)
    pd.testing.assert_frame_equal(table[["S1", "S2"]], expected)
    assert table.GeneSymbol.tolist() == ["A", "B"]
    assert table.loc["g2", "S1"] == 0
    assert pd.isna(table.loc["g2", "S2"])
    for original, current in zip(before, (data.data, data.df_filtered, data.areas_log)):
        pd.testing.assert_frame_equal(original, current)


def test_unnormalized_area_omits_annotations_even_when_globally_enabled(export_data, monkeypatch):
    data = export_data
    data.annotations = True
    calls = []

    def annotate(frame, annotations):
        calls.append(annotations)
        return frame.assign(TaxonID="9606", GeneDescription="extra annotation")

    monkeypatch.setattr(containers, "add_annotations", annotate)

    companion = pd.read_table(data._export_unnormalized_area())
    assert companion.columns.tolist() == ["GeneID", "GeneSymbol", "S1", "S2"]
    assert calls == []

    normal = pd.read_table(data._perform_data_export(level="area"))
    assert "TaxonID" in normal and "GeneDescription" in normal
    assert calls == [True]
    assert data.annotations is True


def test_mspc_creates_and_reports_companion_independently_of_existing_files(export_data):
    data = export_data
    result = CliRunner().invoke(export_command, ["--level", "MSPC"], obj={"data_obj": data})
    assert result.exit_code == 0, result.output
    assert "- area_unnormalized:" in result.output
    main_path, companion_path = map(Path, data._perform_data_export(level="MSPC"))
    assert set(pd.read_table(main_path).GeneID) == {"g1", "g2", "g3"}
    assert set(pd.read_table(companion_path).GeneID) == {"g1", "g2"}
    main_contents = main_path.read_bytes()
    main_mtime = main_path.stat().st_mtime_ns
    companion_path.unlink()

    assert data._perform_data_export(level="MSPC") == (str(main_path), str(companion_path))
    assert main_path.read_bytes() == main_contents
    assert main_path.stat().st_mtime_ns == main_mtime
    companion_mtime = companion_path.stat().st_mtime_ns
    data._perform_data_export(level="MSPC")
    assert companion_path.stat().st_mtime_ns == companion_mtime

    data.df_filtered.loc[("g1", "iBAQ_dstrAdj"), "S2"] = 100.0
    data._perform_data_export(level="MSPC", force=True)
    table = pd.read_table(companion_path, index_col="GeneID")
    assert table.loc["g1", "S2"] == 100.0


@pytest.mark.parametrize("linear", [False, True])
def test_standard_area_export_keeps_analysis_values_and_naming(export_data, linear):
    data = export_data
    path = Path(data._perform_data_export(level="area", linear=linear))
    table = pd.read_table(path, index_col="GeneID")
    assert "med" in path.name
    assert "unnormalized" not in path.name
    assert "Batch" in str(path.parent)
    if linear:
        np.testing.assert_allclose(table.loc["g1", ["S1", "S2"]].astype(float), [1, 10])
    else:
        assert (table[["log10_S1", "log10_S2"]] == 42).all().all()


def test_unnormalized_companion_cannot_overwrite_unnormalized_analysis_export(export_data):
    data = export_data
    data.median = False
    data.batch_applied = None
    analysis_path = Path(data._perform_data_export(level="area"))
    contents = analysis_path.read_bytes()

    companion_path = Path(data._export_unnormalized_area(force=True))

    assert companion_path != analysis_path
    assert analysis_path.read_bytes() == contents


@pytest.mark.parametrize("fill_na_zero", [False, True])
def test_unnormalized_area_preserves_all_zeros(export_data, fill_na_zero):
    data = export_data
    data.fill_na_zero = fill_na_zero
    data.df_filtered.loc[(slice(None), "iBAQ_dstrAdj"), :] = 0

    path = data._export_unnormalized_area()
    values = pd.read_table(path)[["S1", "S2"]]

    assert (values == 0).all().all()


def test_existing_logged_companion_does_not_prevent_linear_export(export_data):
    folder = Path(export_data.outpath) / "export" / "data_area"
    folder.mkdir(parents=True)
    old_path = folder / "demo_data_area_nz1_none_unnormalized.tsv"
    old_contents = "GeneID\tlog10_S1\ng1\t0.045757\n"
    old_path.write_text(old_contents, encoding="utf-8")

    linear_path = Path(export_data._export_unnormalized_area())

    assert linear_path != old_path
    assert pd.read_table(linear_path, index_col="GeneID").loc["g1", "S1"] == 10.0
    assert old_path.read_text(encoding="utf-8") == old_contents


def test_cli1_build_copies_primary_mspc_table_when_export_returns_companion(
    export_data, monkeypatch, tmp_path
):
    from tackle.cli1 import app

    conf = tmp_path / "demo.conf"
    conf.write_text("[S1]\nrecno=1\n", encoding="utf-8")
    analysis_dir = tmp_path / "analysis"
    runner = CliRunner()
    initialized = runner.invoke(
        app.cli,
        ["--analysis-dir", str(analysis_dir), "--conf", str(conf), "--name", "demo", "init"],
    )
    assert initialized.exit_code == 0, initialized.output
    target = analysis_dir / "matrices" / "matrix.mspc.tsv"
    monkeypatch.setattr(app, "_data_from_analysis", lambda *args, **kwargs: export_data)
    monkeypatch.setattr(
        app.AnalysisDir, "canonical_matrix_paths", lambda self: {"mspc_tsv": target}
    )

    built = runner.invoke(app.cli, ["--analysis-dir", str(analysis_dir), "build"])

    assert built.exit_code == 0, built.output
    assert set(pd.read_table(target).GeneID) == {"g1", "g2", "g3"}
