from __future__ import annotations

import argparse
import contextlib
import http.server
import json
import re
import socketserver
import threading
import time
from functools import partial
from pathlib import Path
from typing import Any, Dict, Iterable, List

from tackle.html_overview import build_html_overview


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build and/or profile standalone Tackle HTML overview reports in a local browser. "
            "This is meant for sandbox experimentation around eager vs deferred loading."
        )
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--base-dir",
        type=Path,
        help="Analysis directory to build report variants from.",
    )
    group.add_argument(
        "--html",
        action="append",
        type=Path,
        default=None,
        help="Existing HTML file(s) to profile.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path("/tmp/tackle_html_profile"),
        help="Output root when building variants from --base-dir.",
    )
    parser.add_argument(
        "--self-contained",
        action="store_true",
        default=False,
        help="Build self-contained HTML reports.",
    )
    parser.add_argument(
        "--pngquant",
        action="store_true",
        default=False,
        help="Optimize PNGs during build.",
    )
    parser.add_argument(
        "--interactive-protein-metadata",
        action="store_true",
        default=False,
        help="Embed the protein metadata payload during build.",
    )
    parser.add_argument(
        "--compare-deferred",
        action="store_true",
        default=False,
        help="Build both eager and deferred variants from --base-dir for comparison.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        default=False,
        help="Emit machine-readable JSON.",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        default=False,
        help="Only do static HTML analysis; skip Playwright browser profiling.",
    )
    return parser.parse_args()


@contextlib.contextmanager
def _serve_directory(root: Path) -> Iterable[str]:
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


def _build_variant(
    *,
    base_dir: Path,
    out_dir: Path,
    defer_plot_images: bool,
    self_contained: bool,
    pngquant: bool,
    interactive_protein_metadata: bool,
) -> Path:
    build_html_overview(
        base_dir=str(base_dir),
        out_dir=str(out_dir),
        force=True,
        self_contained=bool(self_contained),
        pngquant=bool(pngquant),
        interactive_protein_metadata=bool(interactive_protein_metadata),
        defer_plot_images=bool(defer_plot_images),
    )
    return out_dir / "index.html"


def _analyze_static_html(html_path: Path) -> Dict[str, Any]:
    text = html_path.read_text(encoding="utf-8")
    eager_sources = re.findall(r'<img[^>]+(?<!data-)src="(data:image/png;base64,[^"]+)"', text)
    deferred_sources = re.findall(r'<img[^>]+data-src="(data:image/png;base64,[^"]+)"', text)
    inline_plot_urls = re.findall(r"data:image/png;base64,[A-Za-z0-9+/=]+", text)
    inline_payload_urls = re.findall(r"data:application/gzip;base64,[A-Za-z0-9+/=]+", text)
    manifest_match = re.search(
        r'<script id="tackle-resource-manifest" type="application/json">(.*?)</script>',
        text,
        flags=re.S,
    )
    resource_mode = None
    if manifest_match:
        try:
            manifest = json.loads(manifest_match.group(1))
            resource_mode = manifest.get("mode")
        except Exception:
            resource_mode = None
    return {
        "inline_plot_src_count": int(len(eager_sources)),
        "inline_plot_deferred_count": int(len(deferred_sources)),
        "inline_plot_data_url_chars": int(sum(len(item) for item in inline_plot_urls)),
        "interactive_payload_data_url_chars": int(sum(len(item) for item in inline_payload_urls)),
        "has_state_api": bool("tackleInteractiveReport" in text and "snapshotState" in text),
        "has_resource_cache_api": bool("setCachedResource" in text and "getCachedResource" in text),
        "resource_mode": resource_mode,
    }


def _profile_html(html_path: Path, *, use_browser: bool = True) -> Dict[str, Any]:
    html_path = html_path.resolve()
    result: Dict[str, Any] = {
        "html_path": str(html_path),
        "file_bytes": html_path.stat().st_size,
        "static": _analyze_static_html(html_path),
        "browser": {
            "available": False,
            "error": None,
            "initial_img_src": None,
            "navigation": None,
            "first_plot_open_ms": None,
            "payload": None,
            "protein_details": None,
            "final_state": None,
        },
    }
    if not use_browser:
        return result

    try:
        from playwright.sync_api import expect, sync_playwright
    except Exception as exc:  # pragma: no cover
        result["browser"]["error"] = (
            "Playwright for Python is required to profile browser behavior. "
            f"Import failed: {exc}"
        )
        return result

    report_dir = html_path.parent

    with _serve_directory(report_dir) as base_url:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch()
            except Exception as exc:  # pragma: no cover
                result["browser"]["error"] = f"Unable to launch Chromium via Playwright: {exc}"
                return result

            page = browser.new_page()
            try:
                page.goto(f"{base_url}/{html_path.name}", wait_until="domcontentloaded")
                page.wait_for_load_state("load")

                nav = page.evaluate(
                    """() => {
                      const entry = performance.getEntriesByType("navigation")[0];
                      const snapshot = window.tackleInteractiveReport?.snapshotState?.() || null;
                      return {
                        domContentLoadedMs: entry ? entry.domContentLoadedEventEnd : null,
                        loadMs: entry ? entry.loadEventEnd : null,
                        encodedBodySize: entry ? entry.encodedBodySize : null,
                        transferSize: entry ? entry.transferSize : null,
                        snapshot: snapshot,
                      };
                    }"""
                )

                first_plot = page.locator("details.plot").first
                first_img = first_plot.locator("img.plot-img").first
                initial_img_src = first_img.get_attribute("src")

                image_open_ms = None
                if first_plot.count():
                    t0 = time.perf_counter()
                    first_plot.locator("summary").first.click()
                    expect(first_img).to_have_attribute("src", re.compile(r".+"))
                    image_open_ms = round((time.perf_counter() - t0) * 1000, 3)

                payload_profile = None
                has_payload = page.evaluate(
                    "() => typeof window.tackleInteractiveReport?.loadPayload === 'function'"
                )
                if has_payload:
                    payload_profile = page.evaluate(
                        """async () => {
                          const t0 = performance.now();
                          const payload = await window.tackleInteractiveReport.loadPayload();
                          return {
                            durationMs: Math.round((performance.now() - t0) * 1000) / 1000,
                            kind: payload && payload.kind ? payload.kind : null,
                            proteinCount: payload && payload.protein_count ? payload.protein_count : null,
                          };
                        }"""
                    )

                protein_profile = None
                details_summary = page.locator("summary", has_text="Top hits (TSV preview)").first
                if details_summary.count():
                    details_summary.click()
                    details_button = page.get_by_role("button", name="Details").first
                    if details_button.count():
                        t0 = time.perf_counter()
                        details_button.click()
                        panel = page.locator(".protein-meta-panel").first
                        panel.wait_for()
                        protein_profile = {
                            "durationMs": round((time.perf_counter() - t0) * 1000, 3),
                            "textSample": panel.text_content()[:200] if panel.text_content() else "",
                        }

                final_state = page.evaluate(
                    "() => window.tackleInteractiveReport?.snapshotState?.() || null"
                )

                result["browser"] = {
                    "available": True,
                    "error": None,
                    "initial_img_src": initial_img_src,
                    "navigation": nav,
                    "first_plot_open_ms": image_open_ms,
                    "payload": payload_profile,
                    "protein_details": protein_profile,
                    "final_state": final_state,
                }
                return result
            finally:
                browser.close()


def _format_summary(item: Dict[str, Any]) -> str:
    static = item.get("static") or {}
    browser = item.get("browser") or {}
    nav = browser.get("navigation") or {}
    payload = browser.get("payload") or {}
    final_state = browser.get("final_state") or {}
    return "\n".join(
        [
            f"HTML: {item['html_path']}",
            f"Bytes: {item['file_bytes']}",
            f"Inline plot src count: {static.get('inline_plot_src_count')}",
            f"Inline plot deferred count: {static.get('inline_plot_deferred_count')}",
            f"Inline plot data-url chars: {static.get('inline_plot_data_url_chars')}",
            f"Interactive payload data-url chars: {static.get('interactive_payload_data_url_chars')}",
            f"Resource mode: {static.get('resource_mode')}",
            f"Has resource cache API: {static.get('has_resource_cache_api')}",
            f"Browser available: {browser.get('available')}",
            f"Browser error: {browser.get('error')}",
            f"DOMContentLoaded: {nav.get('domContentLoadedMs')}",
            f"Load: {nav.get('loadMs')}",
            f"First plot open ms: {browser.get('first_plot_open_ms')}",
            f"Payload decode ms: {payload.get('durationMs')}",
            f"Payload kind: {payload.get('kind')}",
            f"Payload proteins: {payload.get('proteinCount')}",
            f"State active tab: {final_state.get('activeTab') if final_state else None}",
            f"State resource keys: {', '.join(final_state.get('resourceKeys', [])) if final_state else ''}",
            f"State parsed keys: {', '.join(final_state.get('parsedDataKeys', [])) if final_state else ''}",
        ]
    )


def main() -> int:
    args = _parse_args()

    html_paths: List[Path] = []
    if args.base_dir is not None:
        base_dir = args.base_dir.expanduser().resolve()
        variants = [("deferred", True)]
        if args.compare_deferred:
            variants.append(("eager", False))
        for label, defer in variants:
            out_dir = (args.out_root / label).expanduser().resolve()
            html_paths.append(
                _build_variant(
                    base_dir=base_dir,
                    out_dir=out_dir,
                    defer_plot_images=defer,
                    self_contained=bool(args.self_contained),
                    pngquant=bool(args.pngquant),
                    interactive_protein_metadata=bool(args.interactive_protein_metadata),
                )
            )
    else:
        html_paths = [Path(p).expanduser().resolve() for p in (args.html or [])]

    results = [_profile_html(path, use_browser=not bool(args.no_browser)) for path in html_paths]

    if args.json:
        print(json.dumps(results, indent=2, sort_keys=True))
    else:
        for item in results:
            print(_format_summary(item))
            print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
