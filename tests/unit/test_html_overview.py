from __future__ import annotations

from pathlib import Path

import pandas as pd

from tackle.html_overview import build_html_overview


def _write_dummy_png(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\nDUMMY")


def test_build_html_overview_bundle(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    pca_dir = base / "pca" / "mouse"
    metrics_dir = base / "metrics"
    topdiff_dir = base / "topdiff" / "cluster"

    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)
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
    assert "assets/plots/volcano/mouse/run1_volcano_nz1_groupA_minus_groupB_pValue.png" in html_text
    assert "assets/data/volcano/mouse/run1_volcano_nz1_groupA_minus_groupB.tsv" in html_text

    assert (out_dir / "assets" / "plots" / "volcano" / "mouse" / plot_png.name).exists()
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
