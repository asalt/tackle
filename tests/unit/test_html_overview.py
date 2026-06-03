from __future__ import annotations

import base64
import gzip
import html as html_lib
import json
import os
import re
from pathlib import Path

import pandas as pd
import pytest

import tackle.html_overview as html_overview
import tackle.telemetry as telemetry
from tackle.html_overview import build_html_overview


def _write_dummy_png(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\nDUMMY")


def _write_js_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("window.Plotly={newPlot:function(){return Promise.resolve();}};", encoding="utf-8")


def _write_css_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("/* stub css */\n.tabulator{display:block;}\n", encoding="utf-8")


def _extract_resource_manifest(html_text: str) -> dict:
    match = re.search(
        r'<script id="tackle-resource-manifest" type="application/json">(.*?)</script>',
        html_text,
        flags=re.S,
    )
    assert match, "resource manifest script not found"
    return json.loads(html_lib.unescape(match.group(1)))


def test_ai_timeout_default_is_sixty_seconds_and_env_can_override(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("TACKLE_HTML_AI_TIMEOUT_SECONDS", raising=False)
    monkeypatch.delenv("TACKLE_MAKE_HTML_AI_TIMEOUT_SECONDS", raising=False)

    assert html_overview._effective_html_ai_timeout_seconds(60.0) == 60.0

    monkeypatch.setenv("TACKLE_HTML_AI_TIMEOUT_SECONDS", "90")
    assert html_overview._effective_html_ai_timeout_seconds(60.0) == 90.0


def test_collect_metrics_ai_context_reads_numeric_values(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    metrics_dir = base / "metrics"
    metrics_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "sample": ["WT-1", "KO-1"],
            "GPGroups": [1200, 1500],
            "Total_psms": [25000, 30000],
            "Total_peptides": [4100, 4600],
        }
    ).to_csv(metrics_dir / "run_metrics_nz1_med.tsv", sep="\t", index=False)

    ctx = html_overview._collect_metrics_ai_context(root=base)

    assert ctx["available"] is True
    assert ctx["table_count"] == 1
    assert ctx["group_summary_available"] is False
    assert ctx["group_summary"] is None
    table = ctx["tables"][0]
    assert table["rows"] == 2
    assert table["numeric_summary"]["GPGroups"]["median"] == 1350.0
    assert table["sample_rows"][0]["sample"] == "WT-1"
    assert table["sample_rows"][1]["Total_psms"] == 30000


def test_ai_summary_prompt_uses_factual_skeleton_and_allows_takeaway() -> None:
    msg = html_overview._build_ai_summary_prompt(
        section_key="metrics",
        section_label="Metrics",
        prompt_obj={
            "section": {"key": "metrics", "label": "Metrics", "count": 3},
            "metrics": {
                "tables": [
                    {
                        "path": "metrics/run_metrics.tsv",
                        "rows": 2,
                        "numeric_summary": {
                            "GPGroups": {"median": 1350.0, "min": 1200.0, "max": 1500.0}
                        },
                    }
                ]
            },
        },
    )

    assert "factual_skeleton" in msg
    assert "GPGroups: median 1350.0, range 1200.0 to 1500.0." in msg
    assert "Do not compare apparent groups such as WT/KO from sample names" in msg
    assert "optionally end with one short takeaway paragraph" in msg
    assert "Use 4-7 short bullets" not in msg


def test_ai_summary_prompt_labels_cluster_heatmap_scaling_without_inference() -> None:
    msg = html_overview._build_ai_summary_prompt(
        section_key="cluster",
        section_label="Clustering",
        prompt_obj={
            "section": {
                "key": "cluster",
                "label": "Clustering",
                "count": 1,
                "plot_paths": [
                    "clustermap/clustermap_nz1_ccT_cut_geno_lward.D2_med_rcT_rdsl_z_0_9401x10.png"
                ],
            }
        },
    )

    assert "row z-scored" in msg
    assert "mean-centered and scaled across samples" in msg
    assert "geno as metadata/coloring/cut variables" in msg
    assert "not proof that the data are genotyping measurements" in msg


def test_summarize_topdiff_ai_context_extracts_topn_and_metric() -> None:
    summary = html_overview._summarize_topdiff_ai_context(
        [
            {
                "key": "genoKO_minus_genoWT",
                "label": "genoKO minus genoWT",
                "count": 2,
                "groups": [
                    {
                        "key": "default",
                        "label": "default",
                        "count": 2,
                        "plot_paths": [
                            "topdiff/genoKO_minus_genoWT/clustermap_dir_b_log2_FC_z_0_100x10.png",
                            "topdiff/genoKO_minus_genoWT/clustermap_dir_b_pValue_z_0_50x10.png",
                        ],
                    }
                ],
            }
        ]
    )

    group = summary["contrasts"][0]["groups"][0]
    assert group["topn_values"] == [50, 100]
    assert group["ranking_metrics"] == ["log2_FC", "pValue"]
    assert "row z-scored" in group["scaling_note"]


@pytest.fixture(autouse=True)
def _disable_tabulator_vendor_download(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(html_overview, "_find_tabulator_bundle_paths", lambda: None)


def test_build_html_overview_bundle(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    pca_dir = base / "pca" / "mouse"
    metrics_dir = base / "metrics"
    topdiff_dir = base / "topdiff" / "cluster"

    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)
    volcano_topdiff_png = (
        volcano_dir
        / "topdiff"
        / "groupA_minus_groupB"
        / "cmap_nz1_dir_b_pValue_z_0_100x32.png"
    )
    _write_dummy_png(volcano_topdiff_png)
    _write_dummy_png(pca_dir / "run1_pca_pc1_vs_pc2.png")
    _write_dummy_png(metrics_dir / "qc_plot.png")
    _write_dummy_png(topdiff_dir / "topdiff_heatmap.png")

    tsv_path = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv"
    df = pd.DataFrame(
        {
            "GeneID": ["1", "2", "3"],
            "GeneSymbol": ["AAA", "BBB", "CCC"],
            "log2_FC": [1.2, -2.3, 0.5],
            "pValue": [0.01, 0.02, 0.5],
        }
    )
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(base_dir=str(base), out_dir=str(out_dir), force=True)

    assert outputs.out_html.exists()
    html_text = outputs.out_html.read_text(encoding="utf-8")
    assert "Volcano Plots" in html_text
    assert "Top Diff Heatmaps" in html_text
    assert "Volcano Plots (1)" in html_text
    assert "Top Diff Heatmaps (2)" in html_text
    assert "assets/plots/volcano/mouse/run1_volcano_nz1_groupA_minus_groupB_pValue.png" in html_text
    assert (
        "assets/plots/volcano/mouse/topdiff/groupA_minus_groupB/cmap_nz1_dir_b_pValue_z_0_100x32.png"
        in html_text
    )
    assert "assets/data/volcano/mouse/run1_volcano_nz1_groupA_minus_groupB.tsv" in html_text
    manifest = _extract_resource_manifest(html_text)
    assert manifest["mode"] == "chunked"
    assert any(entry["kind"] == "plot_image" for entry in manifest["entries"].values())
    assert "data:image/png;base64," not in html_text

    assert (out_dir / "assets" / "plots" / "volcano" / "mouse" / plot_png.name).exists()
    assert (
        out_dir
        / "assets"
        / "plots"
        / "volcano"
        / "mouse"
        / "topdiff"
        / "groupA_minus_groupB"
        / volcano_topdiff_png.name
    ).exists()
    assert (out_dir / "assets" / "data" / "volcano" / "mouse" / tsv_path.name).exists()


def test_build_html_overview_groups_volcano_by_sort_metric(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"

    png_fc = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_sort_log2_FC.png"
    png_pv = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_sort_pValue.png"
    _write_dummy_png(png_fc)
    _write_dummy_png(png_pv)

    tsv_path = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv"
    df = pd.DataFrame(
        {
            "GeneID": ["1", "2", "3"],
            "GeneSymbol": ["AAA", "BBB", "CCC"],
            "log2_FC": [1.2, -2.3, 0.5],
            "pValue": [0.01, 0.02, 0.5],
        }
    )
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(base_dir=str(base), out_dir=str(out_dir), force=True)

    assert outputs.out_html.exists()
    html_text = outputs.out_html.read_text(encoding="utf-8")
    assert "Sort: log2_FC (1)" in html_text
    assert "Sort: pValue (1)" in html_text


def test_build_html_overview_pngquant_missing_falls_back(tmp_path: Path, monkeypatch) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)

    monkeypatch.setattr(html_overview, "_pngquant_available", lambda: False)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        pngquant=True,
    )

    assert outputs.out_html.exists()
    assert (out_dir / "assets" / "plots" / "volcano" / "mouse" / plot_png.name).exists()


def test_build_html_overview_self_contained_uses_optimized_png(tmp_path: Path, monkeypatch) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    plot_png.parent.mkdir(parents=True, exist_ok=True)
    original = b"\x89PNG\r\n\x1a\nORIGINAL"
    optimized = b"\x89PNG\r\n\x1a\nOPTIMIZED"
    plot_png.write_bytes(original)

    monkeypatch.setattr(html_overview, "_pngquant_available", lambda: True)

    def _fake_opt(src, dst, **kwargs):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(optimized)
        return True

    monkeypatch.setattr(html_overview, "_pngquant_optimize", _fake_opt)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        self_contained=True,
        pngquant=True,
    )

    assert outputs.out_html.exists()
    html_text = outputs.out_html.read_text(encoding="utf-8")
    assert base64.b64encode(optimized).decode("ascii") in html_text
    assert base64.b64encode(original).decode("ascii") not in html_text


def test_build_html_overview_can_defer_self_contained_images(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    plot_png.parent.mkdir(parents=True, exist_ok=True)
    plot_png.write_bytes(b"\x89PNG\r\n\x1a\nDEFER")

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        self_contained=True,
        defer_plot_images=True,
    )

    html_text = outputs.out_html.read_text(encoding="utf-8")
    assert 'class="plot-img plot-img-deferred"' in html_text
    assert 'data-resource-id="' in html_text
    manifest = _extract_resource_manifest(html_text)
    assert manifest["mode"] == "inline"
    assert any(
        entry["kind"] == "plot_image" and entry.get("data_url", "").startswith("data:image/png;base64,")
        for entry in manifest["entries"].values()
    )


def test_build_html_overview_pngquant_topdiff_quality_override(tmp_path: Path, monkeypatch) -> None:
    base = tmp_path / "analysis"
    volcano_png = base / "volcano" / "mouse" / "v1_sort_pValue.png"
    topdiff_png = base / "topdiff" / "cluster" / "topdiff_clustermap.png"
    volcano_png.parent.mkdir(parents=True, exist_ok=True)
    topdiff_png.parent.mkdir(parents=True, exist_ok=True)
    volcano_png.write_bytes(b"\x89PNG\r\n\x1a\nVOLC")
    topdiff_png.write_bytes(b"\x89PNG\r\n\x1a\nTOPD")

    monkeypatch.setattr(html_overview, "_pngquant_available", lambda: True)
    seen_quality = {}

    def _fake_opt(src, dst, **kwargs):
        seen_quality[Path(src).name] = kwargs.get("quality")
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(b"\x89PNG\r\n\x1a\nOPT")
        return True

    monkeypatch.setattr(html_overview, "_pngquant_optimize", _fake_opt)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        pngquant=True,
        pngquant_quality="60-70",
        pngquant_topdiff_quality="88-92",
    )

    assert outputs.out_html.exists()
    assert seen_quality[volcano_png.name] == "60-70"
    assert seen_quality[topdiff_png.name] == "88-92"


def test_build_html_overview_embeds_interactive_payload_gz_json(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)

    payload_obj = {
        "name": "demo",
        "version": 1,
        "values": [1, 2, 3],
        "nested": {"a": True, "b": "x"},
    }
    payload_path = tmp_path / "payload.json"
    payload_path.write_text(json.dumps(payload_obj, indent=2), encoding="utf-8")

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        interactive_payload=str(payload_path),
        self_contained=True,
    )

    html_text = outputs.out_html.read_text(encoding="utf-8")
    manifest = _extract_resource_manifest(html_text)
    entry = manifest["entries"]["interactive_payload"]
    m = re.search(r"data:application/gzip;base64,([A-Za-z0-9+/=]+)", entry["data_url"])
    assert m, "embedded gzip data URL not found"
    gz_bytes = base64.b64decode(m.group(1))
    decoded = json.loads(gzip.decompress(gz_bytes).decode("utf-8"))
    assert decoded == payload_obj


def test_build_html_overview_auto_embeds_protein_metadata_payload(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"

    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)

    volcano_tsv = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv"
    volcano_df = pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "FunCats": ["CAT-A", "CAT-B"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
        }
    )
    volcano_tsv.parent.mkdir(parents=True, exist_ok=True)
    volcano_df.to_csv(volcano_tsv, sep="\t", index=False)

    export_df = pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "TaxonID": [10090, 10090],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "FunCats": ["CAT-A", "CAT-B"],
            "IO": ["S", ""],
            "PeptideCount_57515_1_7": [3, 5],
            "PeptideCount_57516_1_7": [4, 6],
            "PeptideCount_u2g_57515_1_7": [2, 3],
            "PeptideCount_u2g_57516_1_7": [1, 2],
            "PSMs_57515_1_7": [8, 10],
            "PSMs_57516_1_7": [6, 12],
        }
    )
    export_dir.mkdir(parents=True, exist_ok=True)
    export_path = export_dir / "demo_data_MSPC.tsv"
    export_df.to_csv(export_path, sep="\t", index=False)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        interactive_protein_metadata=True,
        self_contained=True,
    )

    assert outputs.out_html.exists()
    html_text = outputs.out_html.read_text(encoding="utf-8")
    assert "protein-meta-btn" in html_text
    assert 'data-protein-id="1"' in html_text
    assert 'data-tab="interactive"' not in html_text
    assert 'id="tab-interactive"' not in html_text

    manifest = _extract_resource_manifest(html_text)
    entry = manifest["entries"]["interactive_payload"]
    m = re.search(r"data:application/gzip;base64,([A-Za-z0-9+/=]+)", entry["data_url"])
    assert m, "embedded gzip data URL not found"
    gz_bytes = base64.b64decode(m.group(1))
    decoded = json.loads(gzip.decompress(gz_bytes).decode("utf-8"))

    assert decoded["kind"] == "protein_metadata"
    assert decoded["protein_count"] == 2
    assert decoded["source_table"] == "export/data_MSPC_demo/paramBatch_plex_noCov_/demo_data_MSPC.tsv"
    assert decoded["proteins"]["1"]["gene_symbol"] == "AAA"
    assert decoded["proteins"]["1"]["annotations"]["IO"] == "S"
    assert decoded["proteins"]["1"]["metrics"]["peptide_count_total"] == 4
    assert decoded["proteins"]["1"]["metrics"]["peptide_count_u2g_total"] == 2
    assert decoded["proteins"]["1"]["metrics"]["psms_total"] == 8
    assert decoded["proteins"]["1"]["runs"]["57515_1_7"]["peptide_count"] == 3


def test_interactive_protein_metadata_uses_max_metric_values_instead_of_summing(
    tmp_path: Path,
) -> None:
    base = tmp_path / "analysis"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"
    export_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        {
            "GeneID": ["1", "1"],
            "GeneSymbol": ["AAA", "AAA"],
            "GeneDescription": ["alpha protein", "alpha protein"],
            "PeptideCount_57515_1_7": [3, 5],
            "PeptideCount_57516_1_7": [4, 2],
            "PSMs_57515_1_7": [8, 6],
            "PSMs_57516_1_7": [6, 12],
        }
    ).to_csv(export_dir / "demo_data_MSPC.tsv", sep="\t", index=False)

    payload = html_overview._build_interactive_protein_metadata_payload(
        root=base,
        pandas_low_memory=False,
    )

    assert payload["proteins"]["1"]["metrics"]["peptide_count_total"] == 5
    assert payload["proteins"]["1"]["metrics"]["psms_total"] == 12
    assert payload["proteins"]["1"]["runs"]["57515_1_7"]["peptide_count"] == 5
    assert payload["proteins"]["1"]["runs"]["57516_1_7"]["peptide_count"] == 4
    assert payload["proteins"]["1"]["runs"]["57515_1_7"]["psms"] == 8
    assert payload["proteins"]["1"]["runs"]["57516_1_7"]["psms"] == 12


def test_build_html_overview_chunked_bundles_plotly_tabulator_and_hides_interactive_tab(
    tmp_path: Path, monkeypatch
) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"

    _write_dummy_png(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png")
    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)

    export_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "PeptideCount_57515_1_7": [3, 5],
            "PSMs_57515_1_7": [8, 10],
        }
    ).to_csv(export_dir / "demo_data_MSPC.tsv", sep="\t", index=False)

    plotly_stub = tmp_path / "vendor" / "plotly.min.js"
    tabulator_js = tmp_path / "vendor" / "tabulator.min.js"
    tabulator_css = tmp_path / "vendor" / "tabulator.min.css"
    _write_js_stub(plotly_stub)
    _write_js_stub(tabulator_js)
    _write_css_stub(tabulator_css)
    monkeypatch.setattr(html_overview, "_find_plotly_bundle_path", lambda: plotly_stub)
    monkeypatch.setattr(
        html_overview,
        "_find_tabulator_bundle_paths",
        lambda: (tabulator_js, tabulator_css),
    )

    out_dir = tmp_path / "report" / "html"
    build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        interactive_protein_metadata=True,
    )

    html_text = (out_dir / "index.html").read_text(encoding="utf-8")
    manifest = _extract_resource_manifest(html_text)

    assert 'data-tab="interactive"' not in html_text
    assert 'id="tab-interactive"' not in html_text
    assert 'data-interactive-volcano-trigger="' in html_text
    assert 'data-interactive-table-trigger="' in html_text
    assert "plotly_bundle" in manifest["entries"]
    assert manifest["entries"]["plotly_bundle"]["kind"] == "script_bundle"
    assert manifest["entries"]["plotly_bundle"]["path"] == "assets/interactive/vendor/plotly.min.js"
    assert manifest["entries"]["tabulator_bundle_js"]["kind"] == "script_bundle"
    assert manifest["entries"]["tabulator_bundle_js"]["path"] == "assets/interactive/vendor/tabulator.min.js"
    assert manifest["entries"]["tabulator_bundle_css"]["kind"] == "style_bundle"
    assert manifest["entries"]["tabulator_bundle_css"]["path"] == "assets/interactive/vendor/tabulator.min.css"
    assert (out_dir / "assets" / "interactive" / "vendor" / "plotly.min.js").exists()
    assert (out_dir / "assets" / "interactive" / "vendor" / "tabulator.min.js").exists()
    assert (out_dir / "assets" / "interactive" / "vendor" / "tabulator.min.css").exists()


def test_build_html_overview_chunked_interactive_resources(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"

    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)

    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)

    export_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "PeptideCount_57515_1_7": [3, 5],
            "PSMs_57515_1_7": [8, 10],
        }
    ).to_csv(export_dir / "demo_data_MSPC.tsv", sep="\t", index=False)

    out_dir = tmp_path / "report" / "html"
    build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        interactive_protein_metadata=True,
    )

    html_text = (out_dir / "index.html").read_text(encoding="utf-8")
    manifest = _extract_resource_manifest(html_text)
    assert manifest["mode"] == "chunked"
    assert "data:image/png;base64," not in html_text

    plot_entries = [entry for entry in manifest["entries"].values() if entry["kind"] == "plot_image"]
    assert plot_entries
    assert plot_entries[0]["path"].startswith("assets/plots/")
    assert (out_dir / plot_entries[0]["path"]).exists()

    payload_entry = manifest["entries"]["interactive_payload"]
    assert payload_entry["path"].startswith("assets/interactive/chunks/")
    assert (out_dir / payload_entry["path"]).exists()

    volcano_entries = [entry for entry in manifest["entries"].values() if entry["kind"] == "volcano_dataset"]
    assert volcano_entries
    assert volcano_entries[0]["path"].startswith("assets/interactive/chunks/")
    assert (out_dir / volcano_entries[0]["path"]).exists()


def test_build_html_overview_interactive_volcano_prefers_padj_when_available(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)

    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
            "pAdj": [0.0, 1e-6],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)

    out_dir = tmp_path / "report" / "html"
    build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
    )

    html_text = (out_dir / "index.html").read_text(encoding="utf-8")
    manifest = _extract_resource_manifest(html_text)
    volcano_entry = next(
        entry for entry in manifest["entries"].values() if entry["kind"] == "volcano_dataset"
    )
    payload = json.loads(gzip.decompress((out_dir / volcano_entry["path"]).read_bytes()).decode("utf-8"))
    assert payload["plot_p_col"] == "pAdj"
    assert payload["rows"][1]["pAdj"] == 1e-06


def test_build_html_overview_ai_summary_clamps_timeout_and_survives_timeout(
    tmp_path: Path, monkeypatch
) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)
    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)

    seen_timeouts = []

    def _fake_from_env(*, agent_api, agent_id=None, local_events_path=None):
        return telemetry.AgentTelemetryConfig(
            agent_api=agent_api,
            agent_id="test-agent",
            timeout_seconds=120.0,
            local_events_path=local_events_path,
        )

    def _fake_support_chat(config, **kwargs):
        seen_timeouts.append(config.timeout_seconds)
        raise TimeoutError("timed out")

    monkeypatch.setenv("TACKLE_AGENT_API", "http://127.0.0.1:9999")
    monkeypatch.setattr(telemetry.AgentTelemetryConfig, "from_env", _fake_from_env)
    monkeypatch.setattr(telemetry, "support_chat", _fake_support_chat)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        ai_summary=True,
    )

    assert outputs.out_html.exists()
    html_text = outputs.out_html.read_text(encoding="utf-8")
    assert "AI summary unavailable" in html_text
    assert seen_timeouts
    assert all(timeout == 30.0 for timeout in seen_timeouts)


def test_build_html_overview_ai_summary_prompt_requests_facts_and_cautious_interpretation(
    tmp_path: Path, monkeypatch
) -> None:
    base = tmp_path / "analysis"
    pca_dir = base / "pca" / "mouse"
    _write_dummy_png(pca_dir / "run1_pca_pc1_vs_pc2.png")

    pd.DataFrame(
        {
            "variable": ["S1", "S2"],
            "genotype": ["WT", "CCL2KO"],
            "age": ["3mo", "24mo"],
            "PC1": [1.2, -0.8],
            "PC2": [0.4, -0.1],
        }
    ).to_csv(pca_dir / "run1_scores.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "PC": ["PC1", "PC2"],
            "Variance": [18.29, 12.47],
        }
    ).to_csv(pca_dir / "run1_variance.tsv", sep="\t", index=False)

    seen_messages = []

    def _fake_from_env(*, agent_api, agent_id=None, local_events_path=None):
        return telemetry.AgentTelemetryConfig(
            agent_api=agent_api,
            agent_id="test-agent",
            timeout_seconds=10.0,
            local_events_path=local_events_path,
        )

    def _fake_support_chat(config, **kwargs):
        seen_messages.append(kwargs.get("message") or "")
        return {"message": "- placeholder summary"}

    monkeypatch.setenv("TACKLE_AGENT_API", "http://127.0.0.1:9999")
    monkeypatch.setattr(telemetry.AgentTelemetryConfig, "from_env", _fake_from_env)
    monkeypatch.setattr(telemetry, "support_chat", _fake_support_chat)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        ai_summary=True,
    )

    assert outputs.out_html.exists()
    assert seen_messages
    prompt = seen_messages[0]
    assert "Begin with factual description grounded in the provided JSON" in prompt
    assert "add one clearly hedged interpretation bullet" in prompt
    assert "Do not claim to see or interpret PNG image contents directly" in prompt
    assert "For PCA, prioritize sample counts, metadata groupings, explained variance" in prompt


def test_build_html_overview_ai_summary_emits_telemetry_events(
    tmp_path: Path, monkeypatch
) -> None:
    base = tmp_path / "analysis"
    pca_dir = base / "pca" / "mouse"
    _write_dummy_png(pca_dir / "run1_pca_pc1_vs_pc2.png")

    pd.DataFrame(
        {
            "variable": ["S1", "S2"],
            "genotype": ["WT", "CCL2KO"],
            "age": ["3mo", "24mo"],
            "PC1": [1.2, -0.8],
            "PC2": [0.4, -0.1],
        }
    ).to_csv(pca_dir / "run1_scores.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "PC": ["PC1", "PC2"],
            "Variance": [18.29, 12.47],
        }
    ).to_csv(pca_dir / "run1_variance.tsv", sep="\t", index=False)

    seen_events = []

    def _fake_from_env(*, agent_api, agent_id=None, local_events_path=None):
        return telemetry.AgentTelemetryConfig(
            agent_api=agent_api,
            agent_id="test-agent",
            timeout_seconds=10.0,
            local_events_path=local_events_path,
        )

    def _fake_support_chat(config, **kwargs):
        return {"message": "- placeholder summary"}

    def _fake_emit_event(self, **kwargs):
        seen_events.append(kwargs)
        return True

    monkeypatch.setenv("TACKLE_AGENT_API", "http://127.0.0.1:9999")
    monkeypatch.setattr(telemetry.AgentTelemetryConfig, "from_env", _fake_from_env)
    monkeypatch.setattr(telemetry, "support_chat", _fake_support_chat)
    monkeypatch.setattr(telemetry.AgentTelemetry, "emit_event", _fake_emit_event)

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        ai_summary=True,
    )

    assert outputs.out_html.exists()
    event_types = [str(item.get("type") or "") for item in seen_events]
    assert "tackle.make_html.ai_summary.start" in event_types
    assert "tackle.make_html.ai_summary.section.start" in event_types
    assert "tackle.make_html.ai_summary.section.complete" in event_types
    assert "tackle.make_html.ai_summary.complete" in event_types


def test_build_html_overview_displays_base_dir_relative_to_cwd(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)

    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    _write_dummy_png(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png")

    out_dir = tmp_path / "report" / "html"
    outputs = build_html_overview(base_dir=str(base), out_dir=str(out_dir), force=True)

    html_text = outputs.out_html.read_text(encoding="utf-8")
    expected_rel = os.path.relpath(base.resolve(), Path.cwd().resolve()).replace("\\", "/")
    if expected_rel == ".":
        expected_rel = "./"
    elif not expected_rel.startswith("."):
        expected_rel = f"./{expected_rel}"

    assert f"Base directory: <code>{expected_rel}</code>" in html_text
    assert str(base.resolve()) not in html_text


def test_build_html_overview_ai_summary_uses_cache_until_force_refresh(
    tmp_path: Path, monkeypatch
) -> None:
    base = tmp_path / "analysis"
    pca_dir = base / "pca" / "mouse"
    _write_dummy_png(pca_dir / "run1_pca_pc1_vs_pc2.png")
    pd.DataFrame(
        {
            "variable": ["S1", "S2"],
            "group": ["WT", "KO"],
            "PC1": [1.2, -0.8],
            "PC2": [0.4, -0.1],
        }
    ).to_csv(pca_dir / "run1_scores.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "PC": ["PC1", "PC2"],
            "Variance": [18.29, 12.47],
        }
    ).to_csv(pca_dir / "run1_variance.tsv", sep="\t", index=False)

    seen_messages = []

    def _fake_from_env(*, agent_api, agent_id=None, local_events_path=None):
        return telemetry.AgentTelemetryConfig(
            agent_api=agent_api,
            agent_id="test-agent",
            timeout_seconds=10.0,
            local_events_path=local_events_path,
        )

    def _fake_support_chat(config, **kwargs):
        seen_messages.append(kwargs.get("message") or "")
        return {"message": f"- cached summary #{len(seen_messages)}"}

    monkeypatch.setenv("TACKLE_AGENT_API", "http://127.0.0.1:9999")
    monkeypatch.setattr(telemetry.AgentTelemetryConfig, "from_env", _fake_from_env)
    monkeypatch.setattr(telemetry, "support_chat", _fake_support_chat)

    out_dir = tmp_path / "report" / "html"
    first = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        ai_summary=True,
    )
    first_html = first.out_html.read_text(encoding="utf-8")
    assert "- cached summary #1" in first_html
    assert len(seen_messages) == 1

    monkeypatch.delenv("TACKLE_AGENT_API", raising=False)
    second = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        ai_summary=True,
    )
    second_html = second.out_html.read_text(encoding="utf-8")
    assert "- cached summary #1" in second_html
    assert len(seen_messages) == 1

    cache_files = list((out_dir / ".ai-cache").glob("*.json"))
    assert len(cache_files) == 1
    assert len(cache_files[0].stem) == 64

    monkeypatch.setenv("TACKLE_AGENT_API", "http://127.0.0.1:9999")
    third = build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        ai_summary=True,
        force_ai_summary=True,
    )
    third_html = third.out_html.read_text(encoding="utf-8")
    assert "- cached summary #2" in third_html
    assert len(seen_messages) == 2
