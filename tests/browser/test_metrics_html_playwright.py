from __future__ import annotations

import contextlib
import http.server
import socketserver
import threading
from functools import partial
from pathlib import Path

import pandas as pd
import pytest

playwright = pytest.importorskip("playwright.sync_api")
expect = playwright.expect
sync_playwright = playwright.sync_playwright

import tackle.metrics_html as metrics_html
from tackle.metrics_html import build_metrics_html_report


def _write_tabulator_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """
window.__tabulatorLoads = (window.__tabulatorLoads || 0) + 1;
window.Tabulator = function(target, options) {
  window.__tabulatorTablesCreated = (window.__tabulatorTablesCreated || 0) + 1;
  const el = typeof target === "string" ? document.querySelector(target) : target;
  this.el = el;
  this.options = options || {};
  this.data = Array.isArray((options || {}).data) ? (options || {}).data.slice() : [];
  this.filteredData = this.data.slice();
  this.render = function() {
    const fields = ((this.options.columns) || []).map(function(col) { return String((col || {}).field || ""); }).join(",");
    el.dataset.stubFields = fields;
    el.dataset.stubVisibleRows = String(this.filteredData.length);
    el.innerHTML = "<div class='tabulator-stub'>Rows " + this.filteredData.length + "/" + this.data.length + "</div>";
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
    window.__metricsDownloads = window.__metricsDownloads || [];
    window.__metricsDownloads.push({ format: format, fileName: fileName });
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


def test_metrics_html_playwright_renders_filters_and_downloads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    metrics_df = pd.DataFrame(
        {
            "S": [10, 8],
            "R": [2, 1],
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
    _write_tabulator_stub(js_stub)
    _write_css_stub(css_stub)
    monkeypatch.setattr(
        metrics_html,
        "_find_tabulator_bundle_paths",
        lambda: (js_stub, css_stub),
    )

    out_dir = tmp_path / "report"
    out_html = out_dir / "metrics.html"
    build_metrics_html_report(
        out_html=str(out_html),
        metrics_df=metrics_df,
        metadata_df=metadata_df,
        miscut_df=miscut_df,
        plot_paths=[],
        export_stem="demo_metrics",
        title="Metrics HTML",
        analysis_label="DemoAnalysis",
    )

    with _serve_directory(out_dir) as base_url:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch()
            except Exception as exc:  # pragma: no cover
                pytest.skip(f"Playwright browser unavailable: {exc}")

            page = browser.new_page()
            try:
                page.goto(f"{base_url}/metrics.html")
                page.wait_for_load_state("domcontentloaded")

                host = page.locator("#metrics-table")
                expect(host.locator(".tabulator-stub")).to_have_text("Rows 2/2")
                expect(host).to_have_attribute(
                    "data-stub-fields",
                    "sample_name,meta__0__genotype,meta__1__age,metric__0__s,metric__1__r,metric__2__gpgroups,metric__3__total-psms,metric__4__total-peptides",
                )
                expect(page.locator("#metrics-status")).to_have_text("2 rows")

                page.locator("#metrics-search").fill("KO")
                expect(host.locator(".tabulator-stub")).to_have_text("Rows 1/2")
                expect(page.locator("#metrics-status")).to_have_text("1 of 2 rows")

                page.locator("#metrics-reset").click()
                expect(host.locator(".tabulator-stub")).to_have_text("Rows 2/2")

                page.locator("#metrics-download").click()
                downloads = page.evaluate("() => window.__metricsDownloads || []")
                assert downloads == [{"format": "csv", "fileName": "demo_metrics.csv"}]

                assert page.evaluate("() => window.__tabulatorTablesCreated || 0") == 1
                page.locator("details.miscut-panel summary").click()
                miscut_host = page.locator("#miscut-table")
                expect(miscut_host.locator(".tabulator-stub")).to_have_text("Rows 2/2")
                assert page.evaluate("() => window.__tabulatorTablesCreated || 0") == 2
            finally:
                browser.close()
