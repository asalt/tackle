from __future__ import annotations

import contextlib
import sys
import types
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import tackle.metrics as metrics_module
import tackle.metrics_html as metrics_html
from tackle.metrics import _miscut_long_frame
from tackle.metrics_html import build_metrics_html_report


def _write_dummy_png(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\nMETRICS")


def _write_js_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("window.Tabulator = function(){ return {}; };", encoding="utf-8")


def _write_css_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(".tabulator{display:block;}", encoding="utf-8")


def test_build_metrics_html_report_inlines_tabulator_and_embeds_plots(
    tmp_path: Path, monkeypatch
) -> None:
    metrics_df = pd.DataFrame(
        {
            "S": [10, 8],
            "R": [2, 1],
            "A": [1, 0],
            "GPGroups": [1200, 980],
            "Total_psms": [25000, 18000],
            "Total_peptides": [4100, 3600],
        },
        index=["SampleA", "SampleB"],
    )
    metadata_df = pd.DataFrame(
        {
            "Genotype": ["WT", "KO"],
            "Age": ["3mo", "24mo"],
        },
        index=["SampleA", "SampleB"],
    )
    miscut_df = pd.DataFrame(
        {
            "miscuts": [0, 1],
            "name": ["SampleA", "SampleB"],
            "Trypsin": [0.82, 0.74],
        }
    )

    js_stub = tmp_path / "vendor" / "tabulator.min.js"
    css_stub = tmp_path / "vendor" / "tabulator.min.css"
    plot_png = tmp_path / "plots" / "demo_metrics_plot.png"
    _write_js_stub(js_stub)
    _write_css_stub(css_stub)
    _write_dummy_png(plot_png)

    monkeypatch.setattr(
        metrics_html,
        "_find_tabulator_bundle_paths",
        lambda: (js_stub, css_stub),
    )

    out_html = tmp_path / "metrics.html"
    result = build_metrics_html_report(
        out_html=str(out_html),
        metrics_df=metrics_df,
        metadata_df=metadata_df,
        miscut_df=miscut_df,
        plot_paths=[plot_png],
        export_stem="demo_metrics",
        title="Metrics HTML",
        analysis_label="DemoAnalysis",
    )

    assert result.exists()
    html_text = result.read_text(encoding="utf-8")
    assert "Metrics HTML" in html_text
    assert "window.Tabulator = function(){ return {}; };" in html_text
    assert ".tabulator{display:block;}" in html_text
    assert "Interactive Metrics Table" in html_text
    assert "Miscut Ratio Detail" in html_text
    assert "Genotype" in html_text
    assert "Total_psms" in html_text
    assert "data:image/png;base64," in html_text


def test_build_metrics_html_report_falls_back_to_static_table_without_tabulator(
    tmp_path: Path, monkeypatch
) -> None:
    metrics_df = pd.DataFrame(
        {
            "GPGroups": [1200, 980],
            "Total_psms": [25000, 18000],
        },
        index=["SampleA", "SampleB"],
    )

    monkeypatch.setattr(metrics_html, "_find_tabulator_bundle_paths", lambda: None)

    out_html = tmp_path / "metrics.html"
    build_metrics_html_report(
        out_html=str(out_html),
        metrics_df=metrics_df,
        metadata_df=None,
        miscut_df=None,
        plot_paths=[],
        export_stem="demo_metrics",
        title="Fallback Metrics HTML",
    )

    html_text = out_html.read_text(encoding="utf-8")
    assert "falling back to a plain HTML table" in html_text
    assert "static-table" in html_text
    assert "window.Tabulator" not in html_text


def test_prepare_metrics_payload_uses_median_cards_not_total_psms() -> None:
    metrics_df = pd.DataFrame(
        {
            "GPGroups": [1200, 980],
            "Total_psms": [25000, 18000],
            "Total_peptides": [4100, 3600],
        },
        index=["SampleA", "SampleB"],
    )

    payload = metrics_html._prepare_metrics_payload(
        metrics_df=metrics_df,
        metadata_df=None,
        miscut_df=None,
    )

    labels = [card["label"] for card in payload["summary_cards"]]
    assert labels == ["Samples", "Median Peptides", "Median GP Groups"]


def test_miscut_long_frame_preserves_real_miscut_categories() -> None:
    frame = pd.DataFrame(
        {
            "SampleA": [8, 2],
            "SampleB": [6, 4],
        },
        index=pd.Index([0, 1], name="TrypsinMiscut"),
    )

    long_frame = _miscut_long_frame(frame, value_name="Trypsin")

    assert list(long_frame.columns) == ["miscuts", "name", "Trypsin"]
    assert set(long_frame["miscuts"].tolist()) == {0, 1}
    counts = long_frame["miscuts"].value_counts().to_dict()
    assert counts == {0: 2, 1: 2}


def test_make_metrics_includes_dist_and_genecount_pngs_in_html(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(metrics_html, "_find_tabulator_bundle_paths", lambda: None)

    def _write_dummy_png(path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"\x89PNG\r\n\x1a\nMETRICS")

    def _fake_get_device(**_kwargs):
        def _device(file):
            _write_dummy_png(Path(file))

        return _device

    def _fake_save_multiple(_fig, outname, *file_fmts, **_kwargs):
        for file_fmt in file_fmts:
            path = Path(str(outname) + str(file_fmt))
            path.parent.mkdir(parents=True, exist_ok=True)
            if str(file_fmt).lower() == ".png":
                _write_dummy_png(path)
            else:
                path.write_text("stub", encoding="utf-8")

    class _DummyConverter:
        def __add__(self, _other):
            return self

    class _DummyR:
        def __init__(self):
            self._handlers = {
                "source": lambda _path: None,
                "metrics": lambda _frame, return_plots=True: OrderedDict(
                    [("overview", object())]
                ),
                "print": lambda _plot: None,
            }

        def __getitem__(self, key):
            return self._handlers[key]

        def assign(self, _name, _value):
            return None

        def __call__(self, _code):
            return None

    dummy_r = _DummyR()
    dummy_converter = _DummyConverter()

    robjects_module = types.ModuleType("rpy2.robjects")
    robjects_module.r = dummy_r
    robjects_module.default_converter = dummy_converter
    robjects_module.pandas2ri = SimpleNamespace(converter=dummy_converter)
    robjects_module.conversion = SimpleNamespace(py2rpy=lambda frame: frame)

    packages_module = types.ModuleType("rpy2.robjects.packages")
    packages_module.importr = lambda _name: SimpleNamespace(dev_off=lambda: None)

    conversion_module = types.ModuleType("rpy2.robjects.conversion")
    conversion_module.localconverter = lambda _converter: contextlib.nullcontext()

    monkeypatch.setitem(sys.modules, "rpy2", types.ModuleType("rpy2"))
    monkeypatch.setitem(sys.modules, "rpy2.robjects", robjects_module)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.packages", packages_module)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.conversion", conversion_module)
    monkeypatch.setattr(metrics_module.grdevice_helper, "get_device", _fake_get_device)
    monkeypatch.setattr(metrics_module, "save_multiple", _fake_save_multiple)

    data_obj = SimpleNamespace(
        metric_values=OrderedDict(
            {
                "SampleA": {
                    "SRA": {"S": 10, "R": 2, "A": 1},
                    "GPGroups": 1200,
                    "PSMs": {"Total": 25000, "u2g": 18000},
                    "Peptides": {
                        "Total": 4100,
                        "u2g": 3200,
                        "Strict": 2200,
                        "Strict_u2g": 1800,
                    },
                    "Area": np.array([1.1, 1.4, 1.8]),
                    "Trypsin": {0: 8, 1: 2},
                    "Trypsin/P": {0: 7, 1: 3},
                },
                "SampleB": {
                    "SRA": {"S": 8, "R": 1, "A": 0},
                    "GPGroups": 980,
                    "PSMs": {"Total": 18000, "u2g": 12000},
                    "Peptides": {
                        "Total": 3600,
                        "u2g": 2800,
                        "Strict": 2000,
                        "Strict_u2g": 1500,
                    },
                    "Area": np.array([0.9, 1.2, 1.6]),
                    "Trypsin": {0: 6, 1: 4},
                    "Trypsin/P": {0: 5, 1: 5},
                },
            }
        ),
        outpath_name="DemoAnalysis",
        outpath=str(tmp_path / "demo"),
        taxon="human",
        non_zeros=1,
        colors_only=False,
        batch_applied=False,
        batch_nonparametric=False,
        normtype="none",
        col_metadata=pd.DataFrame(
            {"Genotype": ["WT", "KO"]},
            index=["SampleA", "SampleB"],
        ),
        data=pd.DataFrame(
            {
                "GeneID": [101, 102, 103],
                "Metric": ["SRA", "SRA", "SRA"],
                "SampleA": ["S", "S", "S"],
                "SampleB": ["S", "R", "S"],
            }
        ),
    )

    outputs = metrics_module.make_metrics(
        data_obj,
        file_fmts=(".png",),
        before_filter=True,
        before_norm=False,
    )

    assert outputs["metrics_html"] is not None
    plot_names = [path.name for path in outputs["plot_pngs"]]
    assert any("overview" in name for name in plot_names)
    assert any("metrics_dist" in name for name in plot_names)
    assert any("metrics_genecounts" in name for name in plot_names)

    html_text = outputs["metrics_html"].read_text(encoding="utf-8")
    assert "Metrics Dist" in html_text
    assert "Metrics Genecounts" in html_text


def test_metrics_plot_dimensions_shrink_for_small_sample_counts() -> None:
    assert metrics_module._metrics_plot_dimensions(3) == (6.8, 7.5)
    assert metrics_module._metrics_plot_dimensions(12) == (9.0, 9.0)
