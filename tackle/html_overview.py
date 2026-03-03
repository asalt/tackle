from __future__ import annotations

import base64
import difflib
import html
import re
import shutil
import json
import os
import subprocess
import tempfile
import uuid
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from .deckgen import _summarize_volcano_top_hits
from .exporter import _collect_tsvs, _read_tsv
from .utils import _get_logger

logger = _get_logger(__name__)


@dataclass
class VolcanoTableItem:
    title: str
    source_relpath: str
    source_path: str
    asset_relpath: str
    preview_table_html: str
    meta: Dict[str, Any]
    top_hits: List[Dict[str, Any]] = field(default_factory=list)
    error: Optional[str] = None


@dataclass
class PlotItem:
    title: str
    category: str
    source_relpath: str
    source_path: str
    asset_relpath: str
    volcano_table: Optional[VolcanoTableItem] = None


@dataclass(frozen=True)
class HtmlOverviewOutputs:
    out_dir: Path
    out_html: Path
    assets_dir: Optional[Path]
    self_contained: bool
    total_plots: int
    total_volcano_tables: int


_CATEGORY_LABELS: Dict[str, str] = {
    "metrics": "Metrics",
    "cluster": "Clustering",
    "pca": "PCA",
    "umap": "UMAP",
    "volcano": "Volcano Plots",
    "topdiff-cluster": "Top Diff Heatmaps",
}

_DEFAULT_ORDER: Tuple[str, ...] = ("metrics", "cluster", "pca", "umap", "volcano", "topdiff-cluster")


def _is_under(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except Exception:
        return False


def _classify_plot_png(relpath: str) -> Optional[str]:
    rel_lower = str(relpath).lower().replace("\\", "/")
    if "topdiff" in rel_lower and ("cluster" in rel_lower or "clustermap" in rel_lower):
        return "topdiff-cluster"
    if "umap" in rel_lower:
        return "umap"
    if "volcano" in rel_lower:
        return "volcano"
    if "pca" in rel_lower or "pcaplot" in rel_lower:
        return "pca"
    if "metrics" in rel_lower:
        return "metrics"
    if "cluster" in rel_lower or "clustermap" in rel_lower:
        return "cluster"
    return None


def _readable_title_from_png_rel(relpath: str, *, category: str) -> str:
    p = Path(relpath)
    stem = "_".join(p.stem.split()).replace("__", "_").strip()
    stem = stem.replace("_", " ").replace("-", " ").strip()
    parent = p.parent.name.replace("_", " ").replace("-", " ").strip() if p.parent.name else ""

    category_label = _CATEGORY_LABELS.get(category, category.title())
    if parent and parent.lower() not in {
        "volcano",
        "pca",
        "pcaplot",
        "umap",
        "metrics",
        "cluster",
        "clustermap",
        "topdiff",
    }:
        return f"{category_label}: {parent} / {stem or p.name}"
    return f"{category_label}: {stem or p.name}"


def _readable_title_from_tsv_rel(relpath: str) -> str:
    base = Path(relpath).name
    if base.lower().endswith(".tsv"):
        base = base[:-4]
    # Prefer extracting everything after 'group' if present (matches exporter/make-xls behavior).
    if "group" in base:
        try:
            suffix = base.split("group", 1)[1].lstrip("_")
            label = suffix
        except Exception:
            label = base
    else:
        label = base
    label = label.replace("__", " ").replace("_", " ").strip()
    while "  " in label:
        label = label.replace("  ", " ")
    return label or base


def _slugify_html_id(text: str) -> str:
    raw = "".join(ch.lower() if ch.isalnum() else "-" for ch in str(text))
    while "--" in raw:
        raw = raw.replace("--", "-")
    raw = raw.strip("-")
    return raw or "x"


def _extract_topdiff_contrast(relpath: str) -> Optional[str]:
    parts = [p for p in str(relpath).replace("\\", "/").split("/") if p]
    lowered = [p.lower() for p in parts]
    if "topdiff" not in lowered:
        return None
    idx = lowered.index("topdiff")
    if idx + 1 >= len(parts):
        return None
    return parts[idx + 1] or None


_TOPDIFF_SHAPE_RE = re.compile(r"(?P<rows>\d+)x(?P<cols>\d+)$")


def _extract_topdiff_top_n(relpath: str) -> Optional[int]:
    stem = Path(str(relpath)).stem
    m = _TOPDIFF_SHAPE_RE.search(stem)
    if not m:
        return None
    try:
        return int(m.group("rows"))
    except Exception:
        return None


def _group_topdiff_by_top_n(items: Sequence[PlotItem]) -> List[Dict[str, Any]]:
    by_topn: Dict[int, List[PlotItem]] = defaultdict(list)
    other: List[PlotItem] = []
    for item in items:
        topn = _extract_topdiff_top_n(item.source_relpath)
        if topn is None:
            other.append(item)
        else:
            by_topn[int(topn)].append(item)

    groups: List[Dict[str, Any]] = []
    for topn, plots in sorted(by_topn.items(), key=lambda kv: -kv[0]):
        groups.append(
            {
                "key": _slugify_html_id(f"top-{topn}"),
                "label": f"Top {topn}",
                "count": len(plots),
                "plots": plots,
                "topn": topn,
            }
        )

    if other:
        groups.append(
            {
                "key": "other",
                "label": "Other",
                "count": len(other),
                "plots": other,
                "topn": None,
            }
        )

    return groups


def _group_topdiff_by_contrast(items: Sequence[PlotItem]) -> List[Dict[str, Any]]:
    by_contrast: Dict[str, List[PlotItem]] = defaultdict(list)
    for item in items:
        contrast = _extract_topdiff_contrast(item.source_relpath) or "Other"
        by_contrast[contrast].append(item)

    out: List[Dict[str, Any]] = []
    used_keys = set()
    for contrast, plots in sorted(by_contrast.items(), key=lambda kv: kv[0].lower()):
        key = _slugify_html_id(contrast)
        base = key
        suffix = 2
        while key in used_keys:
            key = f"{base}-{suffix}"
            suffix += 1
        used_keys.add(key)
        out.append(
            {
                "key": key,
                "label": contrast,
                "count": len(plots),
                "groups": _group_topdiff_by_top_n(plots),
            }
        )
    return out


_VOLCANO_CONTRAST_RE = re.compile(
    r"(?:^|_)group_?(?P<contrast>.*?)(?=_(?:imv|fna|dir|lrob|ltrd|sort|pvalue|padj)_|$)",
    flags=re.IGNORECASE,
)


def _extract_volcano_contrast(relpath: str) -> str:
    stem = Path(str(relpath)).stem
    match = _VOLCANO_CONTRAST_RE.search(stem)
    if match:
        contrast = (match.group("contrast") or "").strip("_")
        if contrast:
            return contrast
    return stem or "Other"


def _extract_volcano_sort_metric(relpath: str) -> str:
    stem = Path(str(relpath)).stem
    lower = stem.lower()
    if re.search(r"(?:^|_)sort_?log2(?:_|-)?fc(?:_|$)", lower):
        return "log2_FC"
    if re.search(r"(?:^|_)sort_?pvalue(?:_|$)", lower):
        return "pValue"
    if re.search(r"(?:^|_)log2(?:_|-)?fc(?:_|$)", lower):
        return "log2_FC"
    if re.search(r"(?:^|_)pvalue(?:_|$)", lower):
        return "pValue"
    return "other"


def _group_volcano_by_contrast_and_sort(items: Sequence[PlotItem]) -> List[Dict[str, Any]]:
    by_contrast: Dict[str, Dict[str, List[PlotItem]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for item in items:
        contrast = _extract_volcano_contrast(item.source_relpath)
        sort_metric = _extract_volcano_sort_metric(item.source_relpath)
        by_contrast[contrast][sort_metric].append(item)

    sort_order = {"log2_FC": 0, "pValue": 1, "other": 2}
    sort_labels = {
        "log2_FC": "Sort: log2_FC",
        "pValue": "Sort: pValue",
        "other": "Other",
    }

    out: List[Dict[str, Any]] = []
    used_keys = set()
    for contrast, grouped in sorted(by_contrast.items(), key=lambda kv: kv[0].lower()):
        contrast_key = _slugify_html_id(contrast)
        base = contrast_key
        suffix = 2
        while contrast_key in used_keys:
            contrast_key = f"{base}-{suffix}"
            suffix += 1
        used_keys.add(contrast_key)

        groups: List[Dict[str, Any]] = []
        for metric, plots in sorted(
            grouped.items(), key=lambda kv: (sort_order.get(kv[0], 99), str(kv[0]))
        ):
            metric_key = _slugify_html_id(f"{contrast}-{metric}")
            groups.append(
                {
                    "key": metric_key,
                    "label": sort_labels.get(metric, f"Sort: {metric}"),
                    "count": len(plots),
                    "sort_metric": metric,
                    "plots": sorted(plots, key=lambda p: p.source_relpath),
                }
            )

        out.append(
            {
                "key": contrast_key,
                "label": contrast,
                "count": sum(group["count"] for group in groups),
                "groups": groups,
            }
        )

    return out


def _render_html_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    pieces = ["<table class='data-table'>", "<thead><tr>"]
    for head in headers:
        pieces.append(f"<th>{html.escape(str(head))}</th>")
    pieces.append("</tr></thead><tbody>")
    for row in rows:
        pieces.append("<tr>")
        for value in row:
            pieces.append(f"<td>{html.escape(str(value))}</td>")
        pieces.append("</tr>")
    pieces.append("</tbody></table>")
    return "".join(pieces)


def _df_preview_table_html(df: pd.DataFrame, *, max_rows: int = 30) -> str:
    if df is None or df.empty:
        return "<p class='muted'>No rows.</p>"
    view = df.head(int(max_rows)).copy()
    headers = [str(c) for c in view.columns]
    rows = list(view.itertuples(index=False, name=None))
    return _render_html_table(headers, rows)


def _copy_or_symlink(src: Path, dst: Path, *, copy: bool, force: bool) -> None:
    if (dst.exists() or dst.is_symlink()) and not force:
        raise FileExistsError(f"{dst} already exists (use --force to overwrite)")
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.parent.mkdir(parents=True, exist_ok=True)
    if copy:
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src.resolve())


def _data_uri_for_file(path: Path, *, mime: str) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def _pngquant_available() -> bool:
    return shutil.which("pngquant") is not None


def _pngquant_optimize(
    src: Path,
    dst: Path,
    *,
    quality: Optional[str] = "65-85",
    speed: int = 3,
    strip: bool = True,
) -> bool:
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["pngquant", "--force", "--skip-if-larger", "--output", str(dst)]
    if quality:
        cmd.append(f"--quality={str(quality)}")
    cmd.append(f"--speed={int(speed)}")
    if strip:
        cmd.append("--strip")
    cmd.extend(["--", str(src)])

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode == 0 and dst.exists() and dst.stat().st_size > 0:
        return True
    return False


def _render_j2_template(template_name: str, **context: object) -> str:
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
    except ImportError as e:
        raise RuntimeError(
            "Jinja2 is required for HTML report rendering. Install with: pip install jinja2"
        ) from e

    template_dir = Path(__file__).resolve().parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(enabled_extensions=("html", "xml", "j2")),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = env.get_template(template_name)
    return template.render(**context)


def _collect_plot_items(
    *,
    base_dir: Path,
    filter_contains: Sequence[str] = (),
    include_kinds: Sequence[str] = _DEFAULT_ORDER,
    exclude_dirs: Sequence[Path] = (),
) -> List[PlotItem]:
    root = base_dir.resolve()
    include = {str(x).strip().lower() for x in include_kinds if str(x).strip()}
    needles = [s for s in filter_contains if s]
    excluded = [p.resolve() for p in exclude_dirs if p]

    # Always exclude report outputs to avoid re-ingesting generated assets.
    excluded.extend(
        [
            (root / "report" / "deck").resolve(),
            (root / "report" / "html").resolve(),
        ]
    )

    items: List[PlotItem] = []
    for png in root.rglob("*.png"):
        png_resolved = png.resolve()
        if any(_is_under(png_resolved, ex) for ex in excluded):
            continue

        rel = png_resolved.relative_to(root)
        rel_str = str(rel).replace("\\", "/")
        if needles and not any(n in rel_str for n in needles):
            continue

        category = _classify_plot_png(rel_str)
        if not category or category not in include:
            continue

        items.append(
            PlotItem(
                title=_readable_title_from_png_rel(rel_str, category=category),
                category=category,
                source_relpath=rel_str,
                source_path=str(png_resolved),
                asset_relpath="",
            )
        )

    items.sort(key=lambda x: (x.category, x.source_relpath))
    return items


def _collect_volcano_tables(
    *,
    base_dir: Path,
    assets_dir: Optional[Path],
    write_assets: bool,
    copy_assets: bool,
    force: bool,
    filter_contains: Sequence[str] = (),
    topn: int = 15,
    direction: str = "both",
    p_cutoff: Optional[float] = 0.05,
    fc_cutoff: Optional[float] = None,
    tsv_engine: str = "auto",
    pandas_low_memory: bool = False,
    max_preview_rows: int = 30,
    exclude_dirs: Sequence[Path] = (),
) -> List[VolcanoTableItem]:
    root = base_dir.resolve()
    needles = [s for s in filter_contains if s]
    excluded = [p.resolve() for p in exclude_dirs if p]
    excluded.append((root / "report").resolve())

    tsvs = _collect_tsvs(root, patterns=[("volcano", "*.tsv")])
    out: List[VolcanoTableItem] = []
    for rel, path in tsvs:
        rel_str = str(rel).replace("\\", "/")
        if needles and not any(n in rel_str for n in needles):
            continue
        if any(_is_under(path.resolve(), ex) for ex in excluded):
            continue

        asset_rel = ""
        if write_assets:
            if assets_dir is None:
                raise ValueError("assets_dir is required when write_assets=True")
            dst = assets_dir / "data" / rel_str
            _copy_or_symlink(path, dst, copy=copy_assets, force=force)
            asset_rel = str(dst.relative_to(assets_dir.parent).as_posix())

        title = _readable_title_from_tsv_rel(rel_str)

        error = None
        meta: Dict[str, Any] = {}
        top_hits: List[Dict[str, Any]] = []
        preview_html = ""
        try:
            df = _read_tsv(path, engine=tsv_engine, low_memory=pandas_low_memory)
            try:
                df_for_summary = df
                # For HTML previews, prefer using raw pValue when available so tables are rarely empty
                # and default filtering doesn't hinge on pAdj cutoffs.
                if "pValue" in df_for_summary.columns and "pAdj" in df_for_summary.columns:
                    df_for_summary = df_for_summary.rename(columns={"pAdj": "_pAdj"})
                picked, meta = _summarize_volcano_top_hits(
                    df_for_summary,
                    topn=int(topn),
                    direction=direction,
                    p_cutoff=p_cutoff,
                    fc_cutoff=fc_cutoff,
                )
                if picked is None or picked.empty:
                    fallback_error = None
                    if p_cutoff is not None and float(p_cutoff) > 0:
                        # If filters remove everything (common when pAdj exists but nothing meets cutoff),
                        # fall back to an unfiltered top-N preview so the report still shows something useful.
                        try:
                            picked2, meta2 = _summarize_volcano_top_hits(
                                df_for_summary,
                                topn=int(topn),
                                direction=direction,
                                p_cutoff=None,
                                fc_cutoff=fc_cutoff,
                            )
                            if picked2 is not None and not picked2.empty:
                                p_col = meta.get("p_col") or meta2.get("p_col") or "p"
                                fallback_error = (
                                    f"No rows matched filters ({p_col}≤{p_cutoff}). "
                                    "Showing top hits without p cutoff."
                                )
                                picked, meta = picked2, meta2
                        except Exception:
                            pass

                    if picked is None or picked.empty:
                        error = fallback_error or "No rows matched filters (or no hits in this direction)."
                        preview_html = _df_preview_table_html(df, max_rows=max_preview_rows)
                    else:
                        error = fallback_error
                        try:
                            top_hits = picked.to_dict(orient="records")
                        except Exception:
                            top_hits = []
                        preview_html = _df_preview_table_html(picked, max_rows=max_preview_rows)
                else:
                    try:
                        top_hits = picked.to_dict(orient="records")
                    except Exception:
                        top_hits = []
                    preview_html = _df_preview_table_html(picked, max_rows=max_preview_rows)
            except Exception as e:
                error = str(e)
                preview_html = _df_preview_table_html(df, max_rows=max_preview_rows)
        except Exception as e:
            error = str(e)
            preview_html = "<p class='muted'>Failed to read TSV.</p>"

        out.append(
            VolcanoTableItem(
                title=title,
                source_relpath=rel_str,
                source_path=str(path.resolve()),
                asset_relpath=asset_rel,
                preview_table_html=preview_html,
                meta=meta,
                top_hits=top_hits,
                error=error,
            )
        )

    out.sort(key=lambda x: x.source_relpath)
    return out


def _attach_volcano_tables(volcano_plots: List[PlotItem], volcano_tables: List[VolcanoTableItem]) -> None:
    by_dir: Dict[str, List[VolcanoTableItem]] = defaultdict(list)
    for table in volcano_tables:
        by_dir[str(Path(table.source_relpath).parent).replace("\\", "/")].append(table)

    for plot in volcano_plots:
        plot_dir = str(Path(plot.source_relpath).parent).replace("\\", "/")
        candidates = by_dir.get(plot_dir, [])
        if not candidates:
            continue

        plot_stem = Path(plot.source_relpath).stem
        best: Optional[VolcanoTableItem] = None
        best_ratio = 0.0
        for cand in candidates:
            cand_stem = Path(cand.source_relpath).stem
            ratio = difflib.SequenceMatcher(None, plot_stem, cand_stem).ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best = cand

        if best is not None and best_ratio >= 0.6:
            plot.volcano_table = best


def _collect_pca_ai_context(
    *,
    root: Path,
    tsv_engine: str,
    pandas_low_memory: bool,
    max_files: int = 3,
    max_rows: int = 100,
) -> List[Dict[str, Any]]:
    items: List[Dict[str, Any]] = []
    pca_root = root / "pca"
    if not pca_root.exists():
        return items

    for scores in sorted(pca_root.rglob("*_scores.tsv"))[: int(max_files)]:
        rel_scores = str(scores.relative_to(root)).replace("\\", "/")
        var_path = scores.with_name(scores.name.replace("_scores.tsv", "_variance.tsv"))
        rel_var = str(var_path.relative_to(root)).replace("\\", "/") if var_path.exists() else None
        try:
            df_scores = _read_tsv(scores, engine=tsv_engine, low_memory=pandas_low_memory)
        except Exception:
            continue

        pc_cols = [c for c in df_scores.columns if str(c).upper().startswith("PC")]
        meta_cols = [c for c in df_scores.columns if c not in pc_cols and c != "variable"]
        keep_cols = [c for c in ["variable", *meta_cols, "PC1", "PC2"] if c in df_scores.columns]
        seen_cols = set()
        keep_cols = [c for c in keep_cols if not (c in seen_cols or seen_cols.add(c))]
        scores_rows = (
            df_scores[keep_cols].head(int(max_rows)).to_dict(orient="records") if keep_cols else []
        )

        variance_rows = None
        if var_path.exists():
            try:
                df_var = _read_tsv(var_path, engine=tsv_engine, low_memory=pandas_low_memory)
                variance_rows = df_var.head(10).to_dict(orient="records")
            except Exception:
                variance_rows = None

        items.append(
            {
                "scores_tsv": rel_scores,
                "variance_tsv": rel_var,
                "rows": scores_rows,
                "variance": variance_rows,
            }
        )

    return items


def _truncate_text(text: str, *, max_chars: int) -> str:
    text = str(text)
    if max_chars < 0:
        return text
    if len(text) <= max_chars:
        return text
    return text[:max_chars] + "…"


def build_html_overview(
    *,
    base_dir: str,
    out_dir: str,
    title: Optional[str] = None,
    show_date: bool = True,
    copy_assets: bool = True,
    force: bool = False,
    self_contained: bool = False,
    filter_contains: Sequence[str] = (),
    include_kinds: Sequence[str] = _DEFAULT_ORDER,
    volcano_topn: int = 15,
    volcano_direction: str = "both",
    volcano_p_cutoff: Optional[float] = 0.05,
    volcano_fc_cutoff: Optional[float] = None,
    max_table_rows: int = 30,
    tsv_engine: str = "auto",
    pandas_low_memory: bool = False,
    ai_summary: bool = False,
    ai_model_label: Optional[str] = None,
    pngquant: bool = False,
    pngquant_quality: Optional[str] = "65-85",
    pngquant_speed: int = 3,
    pngquant_strip: bool = True,
) -> HtmlOverviewOutputs:
    root = Path(base_dir).expanduser().resolve()
    if not root.exists() or not root.is_dir():
        raise FileNotFoundError(str(root))

    out_root = Path(out_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    if not out_root.is_dir():
        raise NotADirectoryError(str(out_root))

    out_html = out_root / "index.html"
    if out_html.exists() and not force:
        raise FileExistsError(f"{out_html} already exists (use --force to overwrite)")

    assets_dir: Optional[Path]
    if self_contained:
        assets_dir = None
    else:
        assets_dir = out_root / "assets"
        assets_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "HTML overview: base_dir=%s out_dir=%s copy_assets=%s force=%s self_contained=%s filter=%s pngquant=%s quality=%s speed=%s strip=%s",
        root,
        out_root,
        copy_assets,
        force,
        self_contained,
        list(filter_contains),
        pngquant,
        pngquant_quality,
        pngquant_speed,
        pngquant_strip,
    )

    # Collect images first (exclude our output dir).
    plots = _collect_plot_items(
        base_dir=root,
        filter_contains=filter_contains,
        include_kinds=include_kinds,
        exclude_dirs=(out_root,),
    )

    # Copy/symlink plot assets into the bundle, or embed them as data URIs.
    use_pngquant = bool(pngquant)
    if use_pngquant and not _pngquant_available():
        logger.warning("pngquant requested but not found on PATH; using original PNG files.")
        use_pngquant = False

    temp_png_dir_ctx = (
        tempfile.TemporaryDirectory(prefix="tackle_make_html_pngquant_")
        if (self_contained and use_pngquant)
        else None
    )
    temp_png_dir = Path(temp_png_dir_ctx.name) if temp_png_dir_ctx is not None else None
    try:
        if self_contained:
            for item in plots:
                source_png = Path(item.source_path)
                embed_png = source_png
                if use_pngquant and temp_png_dir is not None:
                    optimized_png = temp_png_dir / item.source_relpath
                    if _pngquant_optimize(
                        source_png,
                        optimized_png,
                        quality=pngquant_quality,
                        speed=int(pngquant_speed),
                        strip=bool(pngquant_strip),
                    ):
                        embed_png = optimized_png
                item.asset_relpath = _data_uri_for_file(embed_png, mime="image/png")
        else:
            if assets_dir is None:
                raise ValueError("assets_dir is required when self_contained=False")
            for item in plots:
                source_png = Path(item.source_path)
                dst = assets_dir / "plots" / item.source_relpath
                wrote_optimized = False
                if use_pngquant:
                    wrote_optimized = _pngquant_optimize(
                        source_png,
                        dst,
                        quality=pngquant_quality,
                        speed=int(pngquant_speed),
                        strip=bool(pngquant_strip),
                    )
                if not wrote_optimized:
                    _copy_or_symlink(source_png, dst, copy=copy_assets, force=force)
                item.asset_relpath = str(dst.relative_to(out_root).as_posix())
    finally:
        if temp_png_dir_ctx is not None:
            temp_png_dir_ctx.cleanup()

    # Volcano TSV summaries (copied into assets/data/...).
    volcano_tables = _collect_volcano_tables(
        base_dir=root,
        assets_dir=assets_dir,
        write_assets=not self_contained,
        copy_assets=copy_assets,
        force=force,
        filter_contains=filter_contains,
        topn=volcano_topn,
        direction=volcano_direction,
        p_cutoff=volcano_p_cutoff,
        fc_cutoff=volcano_fc_cutoff,
        tsv_engine=tsv_engine,
        pandas_low_memory=pandas_low_memory,
        max_preview_rows=max_table_rows,
        exclude_dirs=(out_root,),
    )

    # Attach volcano tables to volcano plot items when possible.
    volcano_plots = [p for p in plots if p.category == "volcano"]
    _attach_volcano_tables(volcano_plots, volcano_tables)

    attached = {p.volcano_table.source_relpath for p in volcano_plots if p.volcano_table is not None}
    unmatched_tables = [t for t in volcano_tables if t.source_relpath not in attached]

    # Group plots by category in a stable order.
    grouped: Dict[str, List[PlotItem]] = defaultdict(list)
    for item in plots:
        grouped[item.category].append(item)

    sections = []
    for key in include_kinds:
        k = str(key).strip().lower()
        if not k:
            continue
        label = _CATEGORY_LABELS.get(k, k.title())
        items = grouped.get(k, [])
        if len(items) == 0:
            continue
        if k == "topdiff-cluster":
            contrasts = _group_topdiff_by_contrast(items)
            sections.append(
                {
                    "key": k,
                    "label": label,
                    "count": len(items),
                    "plots": items,
                    "contrasts": contrasts,
                }
            )
        elif k == "volcano":
            contrasts = _group_volcano_by_contrast_and_sort(items)
            sections.append(
                {
                    "key": k,
                    "label": label,
                    "count": len(items),
                    "plots": items,
                    "contrasts": contrasts,
                }
            )
        else:
            sections.append({"key": k, "label": label, "count": len(items), "plots": items})

    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    generated_at = now if bool(show_date) else ""
    report_title = title or f"Tackle overview: {root.name}"

    ai_section_summaries: Dict[str, str] = {}
    ai_section_errors: Dict[str, str] = {}
    ai_model_label_value = (str(ai_model_label).strip() if ai_model_label else "").strip()
    if not ai_model_label_value:
        ai_model_label_value = (os.environ.get("TACKLE_AGENT_MODEL_LABEL") or "").strip() or "llama3.1 tulu"

    if ai_summary:
        agent_api = (os.environ.get("TACKLE_AGENT_API") or "").strip()
        if not agent_api:
            msg = (
                "AI summary requested but TACKLE_AGENT_API is not configured "
                "(set it in ~/.ispec/tackle-agent.conf or as an env var)."
            )
            for sec in sections:
                ai_section_errors[str(sec["key"])] = msg
            if unmatched_tables:
                ai_section_errors["volcano-tables"] = msg
        else:
            try:
                from .telemetry import AgentTelemetryConfig, support_chat

                cfg = AgentTelemetryConfig.from_env(agent_api=agent_api)
                session_base = f"tackle-make-html:{uuid.uuid4().hex}"

                pca_context = _collect_pca_ai_context(
                    root=root,
                    tsv_engine=tsv_engine,
                    pandas_low_memory=pandas_low_memory,
                )

                volcano_items = [
                    {
                        "tsv": t.source_relpath,
                        "meta": t.meta,
                        "top_hits": t.top_hits[:50],
                        "note": t.error,
                    }
                    for t in volcano_tables[:20]
                ]
                volcano_unmatched_items = [
                    {
                        "tsv": t.source_relpath,
                        "meta": t.meta,
                        "top_hits": t.top_hits[:50],
                        "note": t.error,
                    }
                    for t in unmatched_tables[:20]
                ]

                for sec in sections:
                    key = str(sec.get("key") or "").strip()
                    if not key:
                        continue
                    try:
                        plot_paths = [p.source_relpath for p in sec.get("plots", [])][:30]
                        prompt_obj: Dict[str, Any] = {
                            "analysis_base_dir": str(root),
                            "generated_at": now,
                            "section": {
                                "key": key,
                                "label": sec.get("label"),
                                "count": sec.get("count"),
                                "plot_paths": plot_paths,
                            },
                        }
                        if key == "pca":
                            prompt_obj["pca"] = pca_context
                        if key == "volcano":
                            prompt_obj["volcano"] = volcano_items
                        if key == "topdiff-cluster":
                            contrasts_prompt: List[Dict[str, Any]] = []
                            for c in sec.get("contrasts", []) or []:
                                groups_prompt: List[Dict[str, Any]] = []
                                for g in c.get("groups", []) or []:
                                    plots = g.get("plots", []) or []
                                    groups_prompt.append(
                                        {
                                            "key": g.get("key"),
                                            "label": g.get("label"),
                                            "topn": g.get("topn"),
                                            "count": g.get("count"),
                                            "plot_paths": [
                                                getattr(p, "source_relpath", str(p)) for p in plots
                                            ][:30],
                                        }
                                    )
                                contrasts_prompt.append(
                                    {
                                        "key": c.get("key"),
                                        "label": c.get("label"),
                                        "count": c.get("count"),
                                        "groups": groups_prompt,
                                    }
                                )
                            prompt_obj["topdiff"] = contrasts_prompt

                        prompt_json = json.dumps(prompt_obj, ensure_ascii=False)
                        prompt_json = _truncate_text(prompt_json, max_chars=16000)

                        msg = (
                            "Write a concise section summary for this Tackle HTML report section. "
                            "Use short bullets. "
                            "Only claim what can be supported by the provided JSON (file names + small tables); "
                            "do not hallucinate image contents.\n\n"
                            f"SECTION_JSON: {prompt_json}"
                        )
                        resp = support_chat(
                            cfg,
                            session_id=f"{session_base}:{key}",
                            message=msg,
                            meta={
                                "_queue_force_inline": True,
                                "source": "tackle_make_html",
                                "analysis_outpath": str(root),
                                "report_title": report_title,
                                "section_key": key,
                                "section_label": sec.get("label"),
                                "model_label": ai_model_label_value,
                            },
                        )
                        text = (resp.get("message") or "").strip()
                        if not text:
                            ai_section_errors[key] = "AI summary endpoint returned an empty response."
                        else:
                            ai_section_summaries[key] = text
                    except Exception as e:
                        ai_section_errors[key] = f"{type(e).__name__}: {e}"

                if unmatched_tables:
                    key = "volcano-tables"
                    try:
                        prompt_obj = {
                            "analysis_base_dir": str(root),
                            "generated_at": now,
                            "section": {
                                "key": key,
                                "label": "Volcano Tables (unmatched)",
                                "count": len(unmatched_tables),
                            },
                            "volcano_unmatched": volcano_unmatched_items,
                        }
                        prompt_json = json.dumps(prompt_obj, ensure_ascii=False)
                        prompt_json = _truncate_text(prompt_json, max_chars=16000)
                        msg = (
                            "Write a concise summary for the unmatched volcano TSV tables in this Tackle HTML report. "
                            "Use short bullets and describe what tables are present and any standout top hits.\n\n"
                            f"SECTION_JSON: {prompt_json}"
                        )
                        resp = support_chat(
                            cfg,
                            session_id=f"{session_base}:{key}",
                            message=msg,
                            meta={
                                "_queue_force_inline": True,
                                "source": "tackle_make_html",
                                "analysis_outpath": str(root),
                                "report_title": report_title,
                                "section_key": key,
                                "section_label": "Volcano Tables (unmatched)",
                                "model_label": ai_model_label_value,
                            },
                        )
                        text = (resp.get("message") or "").strip()
                        if not text:
                            ai_section_errors[key] = "AI summary endpoint returned an empty response."
                        else:
                            ai_section_summaries[key] = text
                    except Exception as e:
                        ai_section_errors[key] = f"{type(e).__name__}: {e}"

            except Exception as e:
                msg = f"{type(e).__name__}: {e}"
                for sec in sections:
                    ai_section_errors[str(sec["key"])] = msg
                if unmatched_tables:
                    ai_section_errors["volcano-tables"] = msg

    html_text = _render_j2_template(
        "overview_report.html.j2",
        title=report_title,
        generated_at=generated_at,
        show_date=bool(show_date),
        base_dir=str(root),
        self_contained=bool(self_contained),
        ai_summary_enabled=bool(ai_summary),
        ai_model_label=ai_model_label_value,
        ai_section_summaries=ai_section_summaries,
        ai_section_errors=ai_section_errors,
        sections=sections,
        total_plots=len(plots),
        volcano_tables_total=len(volcano_tables),
        volcano_tables_unmatched=unmatched_tables,
    )

    out_html.write_text(html_text, encoding="utf-8")
    return HtmlOverviewOutputs(
        out_dir=out_root,
        out_html=out_html,
        assets_dir=assets_dir,
        self_contained=bool(self_contained),
        total_plots=len(plots),
        total_volcano_tables=len(volcano_tables),
    )
