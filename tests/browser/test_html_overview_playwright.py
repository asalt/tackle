from __future__ import annotations

import contextlib
import http.server
import re
import socketserver
import threading
from functools import partial
from pathlib import Path

import pandas as pd
import pytest

playwright = pytest.importorskip("playwright.sync_api")
expect = playwright.expect
sync_playwright = playwright.sync_playwright

import tackle.html_overview as html_overview
from tackle.html_overview import build_html_overview


def _write_dummy_png(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\nBROWSER")


def _write_plotly_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """
window.__plotlyScriptLoads = (window.__plotlyScriptLoads || 0) + 1;
window.Plotly = {
  newPlot: function(el, data, layout) {
    window.__plotlyNewPlotCalls = (window.__plotlyNewPlotCalls || 0) + 1;
    const totalPoints = (data || []).reduce(function(acc, trace) {
      return acc + (((trace && trace.x) || []).length);
    }, 0);
    const maxMarkerSize = (data || []).reduce(function(acc, trace) {
      const size = (((trace || {}).marker) || {}).size;
      return Math.max(acc, Number(size || 0));
    }, 0);
    el.dataset.stubPoints = String(totalPoints);
    el.dataset.stubTitle = String((((layout || {}).title || {}).text) || "");
    el.dataset.stubYUpper = String((((layout || {}).yaxis || {}).range || [])[1] || "");
    el.dataset.stubHeight = String((layout || {}).height || "");
    el.dataset.stubMarkerSize = String(maxMarkerSize || "");
    el.innerHTML = "<div class='plotly-stub'>Plotly stub</div>";
    el.on = function(name, cb) {
      el.__plotlyHandlers = el.__plotlyHandlers || {};
      el.__plotlyHandlers[name] = cb;
    };
    return Promise.resolve(el);
  }
};
        """.strip(),
        encoding="utf-8",
    )


def _write_tabulator_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """
window.__tabulatorScriptLoads = (window.__tabulatorScriptLoads || 0) + 1;
window.Tabulator = function(el, options) {
  window.__tabulatorTablesCreated = (window.__tabulatorTablesCreated || 0) + 1;
  this.el = el;
  this.options = options || {};
  this.data = Array.isArray((options || {}).data) ? (options || {}).data.slice() : [];
  this.filteredData = this.data.slice();
  this.render = function() {
    const columns = (this.options.columns) || [];
    el.dataset.stubTotalRows = String(this.data.length);
    el.dataset.stubVisibleRows = String(this.filteredData.length);
    el.dataset.stubColumns = String(columns.length);
    el.dataset.stubColumnFields = columns.map(function(col) { return String((col || {}).field || ""); }).join(",");
    el.dataset.stubFrozenFields = columns.filter(function(col) { return !!((col || {}).frozen); }).map(function(col) { return String((col || {}).field || ""); }).join(",");
    el.innerHTML = "<div class='tabulator-stub'>Tabulator stub " + this.filteredData.length + "/" + this.data.length + "</div>";
  };
  this.setFilter = function(fn) {
    this.filteredData = typeof fn === "function" ? this.data.filter(function(row) { return !!fn(row); }) : this.data.slice();
    this.render();
  };
  this.clearFilter = function() {
    this.filteredData = this.data.slice();
    this.render();
  };
  this.getDataCount = function(scope) {
    return scope === "active" ? this.filteredData.length : this.data.length;
  };
  this.download = function(format, fileName) {
    window.__tabulatorDownloads = window.__tabulatorDownloads || [];
    window.__tabulatorDownloads.push({ format: format, fileName: fileName });
  };
  this.redraw = function() {
    this.render();
  };
  this.render();
  return this;
};
        """.strip(),
        encoding="utf-8",
    )


def _write_css_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(".tabulator-stub { display: block; }", encoding="utf-8")


@contextlib.contextmanager
def _serve_directory(root: Path):
    handler = partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    with socketserver.TCPServer(("127.0.0.1", 0), handler) as httpd:
        thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        thread.start()
        try:
            port = httpd.server_address[1]
            yield f"http://127.0.0.1:{port}"
        finally:
            httpd.shutdown()
            thread.join(timeout=5)


@pytest.fixture(autouse=True)
def _disable_tabulator_vendor_download(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(html_overview, "_find_tabulator_bundle_paths", lambda: None)


def test_html_overview_playwright_smoke(tmp_path: Path) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"

    plot_png = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png"
    _write_dummy_png(plot_png)

    volcano_tsv = volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv"
    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "FunCats": ["CAT-A", "CAT-B"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
        }
    ).to_csv(volcano_tsv, sep="\t", index=False)

    export_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "IO": ["S", ""],
            "PeptideCount_57515_1_7": [3, 5],
            "PeptideCount_u2g_57515_1_7": [2, 3],
            "PSMs_57515_1_7": [8, 10],
        }
    ).to_csv(export_dir / "demo_data_MSPC.tsv", sep="\t", index=False)

    out_dir = tmp_path / "report" / "html"
    build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        interactive_protein_metadata=True,
        defer_plot_images=True,
    )

    with _serve_directory(out_dir) as base_url:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch()
            except Exception as exc:  # pragma: no cover
                pytest.skip(f"Playwright browser unavailable: {exc}")

            page = browser.new_page()
            try:
                page.goto(f"{base_url}/index.html")
                page.wait_for_load_state("domcontentloaded")

                assert page.evaluate(
                    "() => typeof window.tackleInteractiveReport?.hydrateDeferredImages"
                ) == "function"
                assert page.evaluate(
                    "() => typeof window.tackleInteractiveReport?.loadPayload"
                ) == "function"
                assert page.evaluate(
                    "() => typeof window.tackleInteractiveReport?.loadVolcanoData"
                ) == "function"

                plot = page.locator("details.plot").first
                img = plot.locator("img.plot-img").first
                assert img.get_attribute("src") in (None, "")

                plot.locator("summary").first.click()
                expect(img).to_have_attribute("src", re.compile(r"^assets/plots/"))

                page.locator("summary", has_text="Top hits (TSV preview)").first.click()
                page.get_by_role("button", name="Details").first.click()
                panel = page.locator(".protein-meta-panel")
                panel.wait_for()
                expect(panel).to_contain_text("alpha protein")

                snapshot = page.evaluate(
                    "() => window.tackleInteractiveReport?.snapshotState?.() || null"
                )
                assert snapshot is not None
                assert "interactive_payload" in snapshot["resourceKeys"]
                assert "interactive_payload" in snapshot["parsedDataKeys"]
                assert any(key.startswith("volcano:") for key in snapshot["parsedDataKeys"])
                assert any(key.startswith("protein:") for key in snapshot["parsedDataKeys"])
                assert any(item.get("loaded") for item in snapshot["plotObjects"])
            finally:
                browser.close()


def test_html_overview_playwright_lazy_loads_plotly_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"

    _write_dummy_png(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png")
    _write_dummy_png(volcano_dir / "run2_volcano_nz1_groupC_minus_groupD_pValue.png")

    pd.DataFrame(
        {
            "GeneID": ["1", "2"],
            "GeneSymbol": ["AAA", "BBB"],
            "GeneDescription": ["alpha protein", "beta protein"],
            "FunCats": ["CAT-A", "CAT-B"],
            "log2_FC": [1.5, -1.1],
            "pValue": [0.01, 0.02],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "GeneID": ["1", "3"],
            "GeneSymbol": ["AAA", "CCC"],
            "GeneDescription": ["alpha protein", "gamma protein"],
            "FunCats": ["CAT-A", "CAT-C"],
            "log2_FC": [0.8, -2.4],
            "pValue": [0.03, 0.005],
        }
    ).to_csv(volcano_dir / "run2_volcano_nz1_groupC_minus_groupD.tsv", sep="\t", index=False)

    export_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "GeneID": ["1", "2", "3"],
            "GeneSymbol": ["AAA", "BBB", "CCC"],
            "GeneDescription": ["alpha protein", "beta protein", "gamma protein"],
            "PeptideCount_57515_1_7": [3, 5, 4],
            "PSMs_57515_1_7": [8, 10, 7],
        }
    ).to_csv(export_dir / "demo_data_MSPC.tsv", sep="\t", index=False)

    plotly_stub = tmp_path / "vendor" / "plotly.min.js"
    _write_plotly_stub(plotly_stub)
    monkeypatch.setattr(html_overview, "_find_plotly_bundle_path", lambda: plotly_stub)

    out_dir = tmp_path / "report" / "html"
    build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        interactive_protein_metadata=True,
        defer_plot_images=True,
    )

    with _serve_directory(out_dir) as base_url:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch()
            except Exception as exc:  # pragma: no cover
                pytest.skip(f"Playwright browser unavailable: {exc}")

            page = browser.new_page()
            try:
                page.goto(f"{base_url}/index.html")
                page.wait_for_load_state("domcontentloaded")

                assert page.locator('[data-tab="interactive"]').count() == 0
                assert page.evaluate("() => window.__plotlyScriptLoads || 0") == 0

                plot_details = page.locator("details.plot").filter(has=page.locator("img.plot-img"))
                first_plot = plot_details.nth(0)

                first_plot.locator("summary").first.click()
                first_button = first_plot.locator("[data-interactive-volcano-trigger]").first
                expect(first_button).to_be_visible()
                first_button.click()

                first_host = first_plot.locator(".interactive-volcano-host").first
                expect(first_host).to_have_attribute("data-rendered", "1")
                expect(first_host.locator(".plotly-stub")).to_have_text("Plotly stub")
                expect(first_host).to_have_attribute("data-point-count", "2")
                expect(first_host).to_have_attribute("data-p-col", "pValue")
                expect(first_host).to_have_attribute("data-stub-height", "560")
                expect(first_host).to_have_attribute("data-stub-marker-size", "11")

                snapshot = page.evaluate(
                    "() => window.tackleInteractiveReport?.snapshotState?.() || null"
                )
                assert snapshot is not None
                assert "plotly_bundle" in snapshot["resourceKeys"]
                assert "plotly_bundle" in snapshot["libraryKeys"]
                assert snapshot["pendingLibraryKeys"] == []
                assert page.evaluate("() => window.__plotlyScriptLoads || 0") == 1
                assert page.evaluate("() => window.__plotlyNewPlotCalls || 0") == 1

                contrast_tabs = page.locator("nav.contrast-tabs .contrast-tab")
                if contrast_tabs.count() > 1:
                    contrast_tabs.nth(1).click()

                active_plot = page.locator(".contrast-panel.active details.plot").filter(
                    has=page.locator("img.plot-img")
                ).first
                active_plot.locator("summary").first.click()
                second_button = active_plot.locator("[data-interactive-volcano-trigger]").first
                expect(second_button).to_be_visible()
                second_button.click()

                second_host = active_plot.locator(".interactive-volcano-host").first
                expect(second_host).to_have_attribute("data-rendered", "1")
                expect(second_host).to_have_attribute("data-point-count", "2")
                assert page.evaluate("() => window.__plotlyScriptLoads || 0") == 1
                assert page.evaluate("() => window.__plotlyNewPlotCalls || 0") == 2
            finally:
                browser.close()


def test_html_overview_playwright_interactive_volcano_sets_dynamic_y_range(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"

    _write_dummy_png(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png")
    pd.DataFrame(
        {
            "GeneID": ["1", "2", "3"],
            "GeneSymbol": ["AAA", "BBB", "CCC"],
            "GeneDescription": ["alpha protein", "beta protein", "gamma protein"],
            "log2_FC": [1.5, -1.1, 0.8],
            "pValue": [0.01, 0.02, 0.03],
            "pAdj": [0.0, 1e-6, 0.02],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)

    plotly_stub = tmp_path / "vendor" / "plotly.min.js"
    _write_plotly_stub(plotly_stub)
    monkeypatch.setattr(html_overview, "_find_plotly_bundle_path", lambda: plotly_stub)

    out_dir = tmp_path / "report" / "html"
    build_html_overview(
        base_dir=str(base),
        out_dir=str(out_dir),
        force=True,
        defer_plot_images=True,
    )

    with _serve_directory(out_dir) as base_url:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch()
            except Exception as exc:  # pragma: no cover
                pytest.skip(f"Playwright browser unavailable: {exc}")

            page = browser.new_page()
            try:
                page.goto(f"{base_url}/index.html")
                page.wait_for_load_state("domcontentloaded")

                plot = page.locator("details.plot").filter(has=page.locator("img.plot-img")).first
                plot.locator("summary").first.click()
                plot.locator("[data-interactive-volcano-trigger]").first.click()

                host = plot.locator(".interactive-volcano-host").first
                expect(host).to_have_attribute("data-rendered", "1")
                expect(host).to_have_attribute("data-p-col", "pAdj")
                y_upper = float(host.get_attribute("data-y-upper") or "0")
                p_stretch = float(host.get_attribute("data-p-stretch") or "0")
                assert p_stretch > 0
                assert y_upper > 6.5
            finally:
                browser.close()


def test_html_overview_playwright_lazy_loads_tabulator_once_and_filters_categories(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base = tmp_path / "analysis"
    volcano_dir = base / "volcano" / "mouse"
    export_dir = base / "export" / "data_MSPC_demo" / "paramBatch_plex_noCov_"

    _write_dummy_png(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB_pValue.png")
    _write_dummy_png(volcano_dir / "run2_volcano_nz1_groupC_minus_groupD_pValue.png")

    pd.DataFrame(
        {
            "GeneID": ["1", "2", "3"],
            "GeneSymbol": ["AAA", "BBB", "CCC"],
            "GeneDescription": ["alpha protein", "beta protein", "gamma protein"],
            "log2_FC": [1.5, -1.1, 0.2],
            "pValue": [0.01, 0.02, 0.3],
        }
    ).to_csv(volcano_dir / "run1_volcano_nz1_groupA_minus_groupB.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "GeneID": ["1", "4"],
            "GeneSymbol": ["AAA", "DDD"],
            "GeneDescription": ["alpha protein", "delta protein"],
            "log2_FC": [0.8, -2.4],
            "pValue": [0.03, 0.005],
        }
    ).to_csv(volcano_dir / "run2_volcano_nz1_groupC_minus_groupD.tsv", sep="\t", index=False)

    export_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "GeneID": ["1", "2", "3", "4"],
            "GeneSymbol": ["AAA", "BBB", "CCC", "DDD"],
            "GeneDescription": ["alpha protein", "beta protein", "gamma protein", "delta protein"],
            "CYTO_NUC": ["CY", "", "NUC", ""],
            "MATRISOME": ["", "Core", "", ""],
            "PeptideCount_57515_1_7": [3, 5, 4, 2],
            "PSMs_57515_1_7": [8, 10, 7, 5],
        }
    ).to_csv(export_dir / "demo_data_MSPC.tsv", sep="\t", index=False)

    tabulator_js = tmp_path / "vendor" / "tabulator.min.js"
    tabulator_css = tmp_path / "vendor" / "tabulator.min.css"
    _write_tabulator_stub(tabulator_js)
    _write_css_stub(tabulator_css)
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
        defer_plot_images=True,
    )

    with _serve_directory(out_dir) as base_url:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch()
            except Exception as exc:  # pragma: no cover
                pytest.skip(f"Playwright browser unavailable: {exc}")

            page = browser.new_page()
            try:
                page.goto(f"{base_url}/index.html")
                page.wait_for_load_state("domcontentloaded")

                assert page.evaluate("() => window.__tabulatorScriptLoads || 0") == 0

                first_plot = page.locator("details.plot").filter(
                    has=page.locator("img.plot-img")
                ).first
                first_plot.locator("summary").first.click()

                first_button = first_plot.locator("[data-interactive-table-trigger]").first
                expect(first_button).to_be_visible()
                first_button.click()

                first_host = first_plot.locator(".interactive-table-host").first
                expect(first_host).to_have_attribute("data-rendered", "1")
                expect(first_host.locator(".tabulator-stub")).to_contain_text("Tabulator stub")
                expect(first_host).to_have_attribute("data-row-count", "3")
                expect(first_host).to_have_attribute("data-visible-count", "3")
                expect(
                    first_host
                ).to_have_attribute(
                    "data-stub-column-fields",
                    "GeneID,GeneSymbol,GeneDescription,CYTO_NUC,MATRISOME,PeptideCountTotal,PSMsTotal,log2_FC,pValue,__protein_details",
                )
                expect(first_host).to_have_attribute("data-stub-frozen-fields", "GeneID,GeneSymbol")
                assert page.evaluate("() => window.__tabulatorScriptLoads || 0") == 1
                assert page.evaluate("() => window.__tabulatorTablesCreated || 0") == 1

                cyto_btn = first_plot.locator(".category-filter-btn", has_text="Cyto/Nuc").first
                matri_btn = first_plot.locator(".category-filter-btn", has_text="MATRISOME").first
                expect(cyto_btn).to_be_visible()
                expect(matri_btn).to_be_visible()

                cyto_btn.click()
                expect(first_host).to_have_attribute("data-visible-count", "2")
                expect(
                    first_plot.locator(".interactive-table-status").first
                ).to_contain_text("2 of 3 rows")

                snapshot = page.evaluate(
                    "() => window.tackleInteractiveReport?.snapshotState?.() || null"
                )
                assert snapshot is not None
                assert "tabulator_bundle_js" in snapshot["resourceKeys"]
                assert "tabulator_bundle_js" in snapshot["libraryKeys"]
                table_states = {
                    item["tableId"]: item
                    for item in snapshot.get("tableObjects", [])
                    if item.get("tableId")
                }
                first_table_id = first_button.get_attribute("data-interactive-table-trigger")
                assert first_table_id in table_states
                assert table_states[first_table_id]["visibleRowCount"] == 2
                assert table_states[first_table_id]["activeCategories"] == ["CYTO_NUC"]

                contrast_tabs = page.locator("nav.contrast-tabs .contrast-tab")
                if contrast_tabs.count() > 1:
                    contrast_tabs.nth(1).click()

                second_plot = page.locator(".contrast-panel.active details.plot").filter(
                    has=page.locator("img.plot-img")
                ).first
                second_plot.locator("summary").first.click()
                second_button = second_plot.locator("[data-interactive-table-trigger]").first
                second_button.click()

                second_host = second_plot.locator(".interactive-table-host").first
                expect(second_host).to_have_attribute("data-rendered", "1")
                expect(second_host).to_have_attribute("data-row-count", "2")
                assert page.evaluate("() => window.__tabulatorScriptLoads || 0") == 1
                assert page.evaluate("() => window.__tabulatorTablesCreated || 0") == 2
            finally:
                browser.close()
