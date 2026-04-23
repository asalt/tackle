from __future__ import annotations

import base64
import difflib
import hashlib
import html
import io
import re
import shutil
import json
import os
import subprocess
import tempfile
import uuid
from collections import defaultdict
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.request import Request, urlopen
import gzip

try:
    from tqdm import tqdm
except Exception:
    def tqdm(iterable, *args, **kwargs):  # type: ignore[override]
        return iterable

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
    resource_id: str = ""
    top_hits: List[Dict[str, Any]] = field(default_factory=list)
    error: Optional[str] = None


@dataclass
class PlotItem:
    title: str
    category: str
    source_relpath: str
    source_path: str
    asset_relpath: str
    resource_id: str = ""
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

_PROTEIN_META_BASE_FIELDS: Tuple[str, ...] = (
    "GeneID",
    "TaxonID",
    "GeneSymbol",
    "GeneDescription",
    "Description",
    "FunCats",
    "GeneType",
    "median_isoform_mass",
    "MitoCarta_Pathways",
)

_PROTEIN_META_TAG_FIELDS: Tuple[str, ...] = (
    "IDG",
    "CYTO_NUC",
    "ER_GOLGI",
    "MitoCarta",
    "SurfaceLabel",
    "CellMembrane",
    "MATRISOME",
    "SECRETED",
    "glycomineN",
    "glycomineO",
    "IO",
)

_PROTEIN_META_METRIC_PREFIXES: Tuple[Tuple[str, str], ...] = (
    ("PeptideCount_S_u2g_", "peptide_count_s_u2g"),
    ("PeptideCount_u2g_", "peptide_count_u2g"),
    ("PeptideCount_S_", "peptide_count_s"),
    ("PeptideCount_", "peptide_count"),
    ("PSMs_S_u2g_", "psms_s_u2g"),
    ("PSMs_u2g_", "psms_u2g"),
    ("PSMs_S_", "psms_s"),
    ("PSMs_", "psms"),
)

_PROTEIN_META_METRIC_LABELS: Dict[str, str] = {
    "peptide_count_total": "Max peptide count",
    "peptide_count_u2g_total": "Max unique peptide count",
    "peptide_count_s_total": "Max strict peptide count",
    "peptide_count_s_u2g_total": "Max strict unique peptide count",
    "psms_total": "Max PSMs",
    "psms_u2g_total": "Max unique PSMs",
    "psms_s_total": "Max strict PSMs",
    "psms_s_u2g_total": "Max strict unique PSMs",
}

_TABULATOR_VERSION = "6.3.1"
_TABULATOR_JS_URL = (
    f"https://unpkg.com/tabulator-tables@{_TABULATOR_VERSION}/dist/js/tabulator.min.js"
)
_TABULATOR_CSS_URL = (
    f"https://unpkg.com/tabulator-tables@{_TABULATOR_VERSION}/dist/css/tabulator.min.css"
)


def _is_under(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except Exception:
        return False


def _classify_plot_png(relpath: str) -> Optional[str]:
    rel_lower = str(relpath).lower().replace("\\", "/")
    if rel_lower.startswith("topdiff/") or "/topdiff/" in rel_lower:
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


def _render_html_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    *,
    interactive_protein_lookup: bool = False,
) -> str:
    header_names = [str(head) for head in headers]
    gene_id_idx = header_names.index("GeneID") if "GeneID" in header_names else None
    gene_symbol_idx = header_names.index("GeneSymbol") if "GeneSymbol" in header_names else None
    use_lookup = bool(
        interactive_protein_lookup and (gene_id_idx is not None or gene_symbol_idx is not None)
    )

    table_classes = ["data-table"]
    if use_lookup:
        table_classes.append("protein-data-table")

    pieces = [f"<table class='{' '.join(table_classes)}'>", "<thead><tr>"]
    for head in header_names:
        pieces.append(f"<th>{html.escape(str(head))}</th>")
    if use_lookup:
        pieces.append("<th>Protein</th>")
    pieces.append("</tr></thead><tbody>")
    for row in rows:
        row_values = list(row)
        gene_id = ""
        gene_symbol = ""
        if gene_id_idx is not None and gene_id_idx < len(row_values):
            gene_id = str(row_values[gene_id_idx]).strip()
        if gene_symbol_idx is not None and gene_symbol_idx < len(row_values):
            gene_symbol = str(row_values[gene_symbol_idx]).strip()

        row_attrs: List[str] = []
        if use_lookup and gene_id:
            row_attrs.append(f"data-protein-id=\"{html.escape(gene_id)}\"")
        if use_lookup and gene_symbol:
            row_attrs.append(f"data-protein-symbol=\"{html.escape(gene_symbol)}\"")
        attrs = f" {' '.join(row_attrs)}" if row_attrs else ""

        pieces.append(f"<tr{attrs}>")
        for value in row_values:
            pieces.append(f"<td>{html.escape(str(value))}</td>")
        if use_lookup:
            if gene_id or gene_symbol:
                parts = ["<button class='protein-meta-btn' type='button'"]
                if gene_id:
                    parts.append(f" data-protein-id=\"{html.escape(gene_id)}\"")
                if gene_symbol:
                    parts.append(f" data-protein-symbol=\"{html.escape(gene_symbol)}\"")
                parts.append(">Details</button>")
                pieces.append(f"<td>{''.join(parts)}</td>")
            else:
                pieces.append("<td><span class='muted'>n/a</span></td>")
        pieces.append("</tr>")
    pieces.append("</tbody></table>")
    return "".join(pieces)


def _df_preview_table_html(
    df: pd.DataFrame,
    *,
    max_rows: int = 30,
    interactive_protein_lookup: bool = False,
) -> str:
    if df is None or df.empty:
        return "<p class='muted'>No rows.</p>"
    view = df.head(int(max_rows)).copy()
    headers = [str(c) for c in view.columns]
    rows = list(view.itertuples(index=False, name=None))
    return _render_html_table(
        headers,
        rows,
        interactive_protein_lookup=interactive_protein_lookup,
    )


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


def _data_uri_for_bytes(data: bytes, *, mime: str) -> str:
    encoded = base64.b64encode(data).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_bytes(dst: Path, data: bytes, *, force: bool) -> None:
    if dst.exists() and not force:
        raise FileExistsError(f"{dst} already exists (use --force to overwrite)")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(data)


def _make_resource_id(kind: str, relpath: str) -> str:
    digest = hashlib.sha1(str(relpath).encode("utf-8")).hexdigest()[:12]
    stem = _slugify_html_id(Path(str(relpath)).stem)[:40]
    prefix = _slugify_html_id(kind)[:20]
    return f"{prefix}-{stem}-{digest}"


def _gzip_bytes(data: bytes, *, compresslevel: int = 9, mtime: int = 0) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(
        fileobj=buf,
        mode="wb",
        compresslevel=int(compresslevel),
        mtime=int(mtime),
    ) as handle:
        handle.write(data)
    return buf.getvalue()


def _clean_payload_scalar(value: Any, *, float_digits: Optional[int] = None) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    if hasattr(value, "item"):
        try:
            value = value.item()
        except Exception:
            pass

    if isinstance(value, str):
        cleaned = value.strip()
        return cleaned or None
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        if float_digits is not None:
            value = round(float(value), int(float_digits))
        if float(value).is_integer():
            return int(value)
        return float(value)
    return value


def _clean_json_row(row: Dict[str, Any], *, float_digits: Optional[int] = None) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, value in row.items():
        cleaned = _clean_payload_scalar(value, float_digits=float_digits)
        if cleaned is None:
            continue
        out[str(key)] = cleaned
    return out


def _parse_protein_metric_column(name: str) -> Optional[Tuple[str, str]]:
    text = str(name)
    for prefix, metric_key in _PROTEIN_META_METRIC_PREFIXES:
        if text.startswith(prefix):
            run_id = text[len(prefix) :].strip("_")
            if run_id:
                return metric_key, run_id
    return None


def _merge_protein_metric_value(prior: Any, current: Any) -> Any:
    if prior is None:
        return current
    if current is None:
        return prior
    try:
        merged = max(float(prior), float(current))
    except Exception:
        return current
    return _clean_payload_scalar(merged, float_digits=2)


def _find_protein_metadata_export_tsv(root: Path) -> Optional[Tuple[Path, List[str]]]:
    export_root = root / "export"
    if not export_root.exists() or not export_root.is_dir():
        return None

    ranked: List[Tuple[int, int, str, Path, List[str]]] = []
    for path in sorted(export_root.rglob("*.tsv")):
        try:
            header_df = pd.read_csv(path, sep="\t", nrows=0)
        except Exception:
            continue

        columns = [str(c) for c in header_df.columns]
        if "GeneID" not in columns:
            continue

        score = 0
        lower_path = str(path.relative_to(root)).lower().replace("\\", "/")
        if "data_mspc" in lower_path:
            score += 100
        if "protein" in lower_path:
            score += 10
        if "GeneSymbol" in columns:
            score += 15
        if "GeneDescription" in columns or "Description" in columns:
            score += 15

        metric_cols = [c for c in columns if _parse_protein_metric_column(c) is not None]
        if not metric_cols:
            continue
        score += min(len(metric_cols), 25)

        ranked.append((score, len(columns), lower_path, path, columns))

    if not ranked:
        return None

    ranked.sort(key=lambda item: (-item[0], -item[1], item[2]))
    _, _, _, best_path, best_cols = ranked[0]
    return best_path, best_cols


def _build_interactive_protein_metadata_payload(
    *,
    root: Path,
    pandas_low_memory: bool,
) -> Dict[str, Any]:
    found = _find_protein_metadata_export_tsv(root)
    if found is None:
        raise RuntimeError(
            "Interactive protein metadata needs an export TSV with GeneID and peptide/PSM columns "
            "(for example export/data_MSPC*/.../*.tsv)."
        )

    source_path, source_columns = found
    base_cols = [c for c in _PROTEIN_META_BASE_FIELDS if c in source_columns]
    tag_cols = [c for c in _PROTEIN_META_TAG_FIELDS if c in source_columns]
    metric_cols = [c for c in source_columns if _parse_protein_metric_column(c) is not None]

    usecols: List[str] = []
    for col in [*base_cols, *tag_cols, *metric_cols]:
        if col not in usecols:
            usecols.append(col)

    if "GeneID" not in usecols:
        usecols.insert(0, "GeneID")

    df = pd.read_csv(
        source_path,
        sep="\t",
        usecols=usecols,
        low_memory=pandas_low_memory,
        dtype={"GeneID": str},
    )
    if "GeneID" not in df.columns:
        raise RuntimeError(f"Protein metadata source is missing GeneID: {source_path}")

    proteins: Dict[str, Dict[str, Any]] = {}
    symbol_to_gene_id: Dict[str, str] = {}
    ambiguous_symbols = set()
    pretty_field_names = {
        "TaxonID": "taxon_id",
        "GeneSymbol": "gene_symbol",
        "GeneDescription": "gene_description",
        "Description": "gene_description",
        "FunCats": "fun_cats",
        "GeneType": "gene_type",
        "median_isoform_mass": "median_isoform_mass",
        "MitoCarta_Pathways": "mitocarta_pathways",
    }

    for row in df.to_dict(orient="records"):
        gene_id_raw = _clean_payload_scalar(row.get("GeneID"))
        if gene_id_raw is None:
            continue
        gene_id = str(gene_id_raw)
        protein = proteins.setdefault(gene_id, {"gene_id": gene_id})

        for col in base_cols:
            if col == "GeneID":
                continue
            val = _clean_payload_scalar(row.get(col), float_digits=2)
            if val is None:
                continue
            key = pretty_field_names.get(col, col.lower())
            if key not in protein or protein.get(key) in (None, ""):
                protein[key] = val

        annotations = protein.setdefault("annotations", {})
        for col in tag_cols:
            val = _clean_payload_scalar(row.get(col))
            if val is None:
                continue
            annotations[str(col)] = val
        if not annotations:
            protein.pop("annotations", None)

        metrics = protein.setdefault("metrics", {})
        runs = protein.setdefault("runs", {})
        for col in metric_cols:
            parsed = _parse_protein_metric_column(col)
            if parsed is None:
                continue
            metric_key, run_id = parsed
            val = _clean_payload_scalar(row.get(col), float_digits=2)
            if val is None:
                continue

            total_key = f"{metric_key}_total"
            metrics[total_key] = _merge_protein_metric_value(metrics.get(total_key), val)

            run_metrics = runs.setdefault(run_id, {})
            run_metrics[metric_key] = _merge_protein_metric_value(run_metrics.get(metric_key), val)
        if not metrics:
            protein.pop("metrics", None)
        if not runs:
            protein.pop("runs", None)

        gene_symbol = protein.get("gene_symbol")
        if gene_symbol:
            symbol = str(gene_symbol)
            existing = symbol_to_gene_id.get(symbol)
            if existing is None:
                symbol_to_gene_id[symbol] = gene_id
            elif existing != gene_id:
                ambiguous_symbols.add(symbol)

    for symbol in ambiguous_symbols:
        symbol_to_gene_id.pop(symbol, None)

    return {
        "schema_version": 1,
        "kind": "protein_metadata",
        "analysis_base_dir": str(root),
        "source_table": str(source_path.relative_to(root)).replace("\\", "/"),
        "protein_count": int(len(proteins)),
        "metric_labels": dict(_PROTEIN_META_METRIC_LABELS),
        "proteins": proteins,
        "symbol_to_gene_id": symbol_to_gene_id,
    }


def _pick_interactive_volcano_p_column(columns: Sequence[str]) -> Optional[str]:
    cols = {str(col) for col in columns}
    if "pAdj" in cols:
        return "pAdj"
    if "pValue" in cols:
        return "pValue"
    return None


def _build_volcano_dataset_payload(
    *,
    df: pd.DataFrame,
    source_relpath: str,
    title: str,
    meta: Dict[str, Any],
    top_hits: List[Dict[str, Any]],
    error: Optional[str],
) -> Dict[str, Any]:
    records = [
        _clean_json_row(row, float_digits=10)
        for row in df.to_dict(orient="records")
    ]
    columns = [str(col) for col in df.columns]
    plot_p_col = _pick_interactive_volcano_p_column(columns)
    return {
        "schema_version": 1,
        "kind": "volcano_dataset",
        "source_table": str(source_relpath),
        "title": str(title),
        "row_count": int(len(records)),
        "columns": columns,
        "plot_p_col": str(plot_p_col) if plot_p_col else None,
        "rows": records,
        "meta": dict(meta or {}),
        "top_hits": [
            _clean_json_row(dict(row), float_digits=10) for row in (top_hits or [])
        ],
        "error": str(error) if error else None,
    }


def _encode_interactive_payload_obj_gz_json(payload_obj: Any, *, source_name: str) -> Dict[str, Any]:
    minified = json.dumps(payload_obj, ensure_ascii=False, separators=(",", ":"))
    json_bytes = minified.encode("utf-8")
    gz_bytes = _gzip_bytes(json_bytes, compresslevel=9, mtime=0)
    b64 = base64.b64encode(gz_bytes).decode("ascii")
    sha256 = hashlib.sha256(gz_bytes).hexdigest()
    return {
        "source_name": source_name,
        "mime": "application/gzip",
        "data_url": f"data:application/gzip;base64,{b64}",
        "sha256": sha256,
        "json_bytes": int(len(json_bytes)),
        "gz_bytes": int(len(gz_bytes)),
        "b64_chars": int(len(b64)),
    }


def _effective_resource_mode(*, self_contained: bool, interactive_resource_mode: str) -> str:
    mode = str(interactive_resource_mode or "auto").strip().lower() or "auto"
    if mode not in {"auto", "inline", "chunked"}:
        raise RuntimeError(
            f"Unsupported interactive resource mode: {interactive_resource_mode}"
        )
    if mode == "auto":
        return "inline" if self_contained else "chunked"
    if mode == "chunked" and self_contained:
        raise RuntimeError(
            "Chunked interactive resources are not compatible with --self-contained. "
            "Use --interactive-resource-mode inline or disable --self-contained."
        )
    return mode


def _build_json_resource_entry(
    *,
    resource_id: str,
    payload_obj: Any,
    source_name: str,
    kind: str,
    resource_mode: str,
    out_root: Path,
    force: bool,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    encoded = _encode_interactive_payload_obj_gz_json(payload_obj, source_name=source_name)
    entry = {
        "id": str(resource_id),
        "kind": str(kind),
        "format": "json",
        "compression": "gzip",
        "source_name": str(source_name),
        "sha256": str(encoded["sha256"]),
        "json_bytes": int(encoded["json_bytes"]),
        "bytes": int(encoded["gz_bytes"]),
    }
    if resource_mode == "inline":
        entry["data_url"] = str(encoded["data_url"])
    else:
        relpath = f"assets/interactive/chunks/{resource_id}.json.gz"
        _write_bytes(out_root / relpath, base64.b64decode(encoded["data_url"].split(",", 1)[1]), force=force)
        entry["path"] = relpath
    return entry, encoded


def _build_plot_resource_entry(
    *,
    resource_id: str,
    item: PlotItem,
    resource_mode: str,
    data_url: Optional[str],
    path: Optional[str],
    source_bytes: bytes,
) -> Dict[str, Any]:
    entry = {
        "id": str(resource_id),
        "kind": "plot_image",
        "format": "png",
        "compression": "none",
        "source_name": str(Path(item.source_relpath).name),
        "source_relpath": str(item.source_relpath),
        "title": str(item.title),
        "sha256": _sha256_bytes(source_bytes),
        "bytes": int(len(source_bytes)),
    }
    if resource_mode == "inline":
        if not data_url:
            raise RuntimeError(f"Inline plot resource is missing data URL: {item.source_relpath}")
        entry["data_url"] = str(data_url)
    else:
        if not path:
            raise RuntimeError(f"Chunked plot resource is missing path: {item.source_relpath}")
        entry["path"] = str(path)
    return entry


def _build_file_resource_entry(
    *,
    resource_id: str,
    source_path: Path,
    source_name: str,
    kind: str,
    format: str,
    mime: str,
    resource_mode: str,
    out_root: Path,
    path: Optional[str],
    copy_assets: bool,
    force: bool,
) -> Dict[str, Any]:
    source_bytes = source_path.read_bytes()
    entry = {
        "id": str(resource_id),
        "kind": str(kind),
        "format": str(format),
        "compression": "none",
        "source_name": str(source_name),
        "sha256": _sha256_bytes(source_bytes),
        "bytes": int(len(source_bytes)),
    }
    if resource_mode == "inline":
        entry["data_url"] = _data_uri_for_bytes(source_bytes, mime=mime)
        return entry

    if not path:
        raise RuntimeError(f"Chunked resource is missing a destination path: {source_name}")

    dst = out_root / str(path)
    _copy_or_symlink(source_path, dst, copy=copy_assets, force=force)
    entry["path"] = str(path)
    return entry


def _encode_interactive_payload_gz_json(payload_path: Path) -> Dict[str, Any]:
    text = payload_path.read_text(encoding="utf-8")
    obj = json.loads(text)
    return _encode_interactive_payload_obj_gz_json(obj, source_name=payload_path.name)


def _find_plotly_bundle_path() -> Optional[Path]:
    try:
        import plotly  # type: ignore
    except Exception:
        return None

    candidate = Path(plotly.__file__).resolve().parent / "package_data" / "plotly.min.js"
    if candidate.exists() and candidate.is_file():
        return candidate
    return None


def _download_vendor_asset(*, url: str, dst: Path) -> Path:
    if dst.exists() and dst.is_file():
        return dst

    dst.parent.mkdir(parents=True, exist_ok=True)
    req = Request(
        str(url),
        headers={
            "User-Agent": "tackle-make-html/1.0",
            "Accept": "*/*",
        },
    )
    with urlopen(req, timeout=30) as response:
        data = response.read()
    if not data:
        raise RuntimeError(f"Downloaded empty vendor asset from {url}")

    tmp_dst = dst.with_name(f"{dst.name}.{uuid.uuid4().hex}.tmp")
    tmp_dst.write_bytes(data)
    tmp_dst.replace(dst)
    return dst


def _find_tabulator_bundle_paths() -> Optional[Tuple[Path, Path]]:
    env_js = (os.environ.get("TACKLE_TABULATOR_JS") or "").strip()
    env_css = (os.environ.get("TACKLE_TABULATOR_CSS") or "").strip()
    if env_js or env_css:
        if not (env_js and env_css):
            logger.warning(
                "Tabulator bundle overrides require both TACKLE_TABULATOR_JS and TACKLE_TABULATOR_CSS."
            )
            return None
        js_path = Path(env_js).expanduser().resolve()
        css_path = Path(env_css).expanduser().resolve()
        if not js_path.exists() or not js_path.is_file():
            logger.warning("Tabulator JS override does not exist: %s", js_path)
            return None
        if not css_path.exists() or not css_path.is_file():
            logger.warning("Tabulator CSS override does not exist: %s", css_path)
            return None
        return js_path, css_path

    cache_dir = Path(tempfile.gettempdir()) / "tackle_vendor_assets" / "tabulator" / _TABULATOR_VERSION
    js_path = cache_dir / "tabulator.min.js"
    css_path = cache_dir / "tabulator.min.css"
    try:
        _download_vendor_asset(url=_TABULATOR_JS_URL, dst=js_path)
        _download_vendor_asset(url=_TABULATOR_CSS_URL, dst=css_path)
    except Exception as exc:
        logger.warning("Unable to resolve Tabulator vendor assets: %s", exc)
        return None
    return js_path, css_path


def _effective_html_ai_timeout_seconds(timeout_seconds: float) -> float:
    raw = (
        os.environ.get("TACKLE_HTML_AI_TIMEOUT_SECONDS")
        or os.environ.get("TACKLE_MAKE_HTML_AI_TIMEOUT_SECONDS")
        or ""
    ).strip()
    if raw:
        try:
            return max(1.0, float(raw))
        except ValueError:
            pass
    return max(1.0, min(float(timeout_seconds), 30.0))


def _build_ai_summary_prompt(
    *,
    section_key: str,
    section_label: str,
    prompt_obj: Mapping[str, Any],
    unmatched_volcano_tables: bool = False,
) -> str:
    base_rules = [
        "Write a concise summary for this Tackle HTML report section.",
        "Use 4-7 short bullets.",
        "Begin with factual description grounded in the provided JSON: counts, files, tables, contrasts, sample/group labels, and numeric values when available.",
        "Do not claim to see or interpret PNG image contents directly.",
        "If the JSON supports it, add one clearly hedged interpretation bullet using language like 'suggests', 'is consistent with', or 'may indicate'.",
        "Prefer describing the generated files and exported measurements over repeating long file paths.",
        "If evidence is weak or incomplete, say so plainly instead of guessing.",
    ]

    section_rules: List[str] = []
    key = str(section_key or "").strip().lower()
    if unmatched_volcano_tables or key == "volcano-tables":
        section_rules.extend(
            [
                "For unmatched volcano tables, summarize which TSV outputs are present, what contrasts they appear to represent, and any top-hit patterns visible in the exported rows.",
                "Do not imply that a missing PNG was reviewed.",
            ]
        )
    elif key == "pca":
        section_rules.extend(
            [
                "For PCA, prioritize sample counts, metadata groupings, explained variance, and whether score tables suggest separation or overlap along PC axes.",
                "Any interpretation about separation must be based on scores/variance TSV values and sample annotations, not on the unseen image.",
            ]
        )
    elif key == "volcano":
        section_rules.extend(
            [
                "For volcano sections, describe contrasts, significance metrics/cutoffs, directionality, and notable top-hit genes from the exported tables when present.",
                "Any biological interpretation must stay cautious and be tied to the reported statistics rather than the plot image itself.",
            ]
        )
    elif key == "topdiff-cluster":
        section_rules.extend(
            [
                "For top-diff heatmaps, describe the contrasts, top-N variants, and output organization from the file names and metadata.",
                "Do not claim clustering patterns unless they are supported by accompanying numeric context.",
            ]
        )
    elif key in {"metrics", "cluster", "umap"}:
        section_rules.append(
            f"For {section_label or section_key}, describe what artifacts were produced and any concrete numeric context present in the provided JSON."
        )

    prompt_json = json.dumps(prompt_obj, ensure_ascii=False)
    prompt_json = _truncate_text(prompt_json, max_chars=16000)

    lines = base_rules[:]
    if section_rules:
        lines.append("Section-specific guidance:")
        lines.extend(section_rules)
    lines.append("")
    lines.append(f"SECTION_JSON: {prompt_json}")
    return "\n".join(lines)


def _ai_summary_cache_dir(out_root: Path) -> Path:
    return out_root / ".ai-cache"


def _ai_summary_cache_path(cache_dir: Path, *, message: str) -> Path:
    digest = hashlib.sha256(str(message).encode("utf-8")).hexdigest()
    return cache_dir / f"{digest}.json"


def _load_ai_summary_cache(cache_dir: Path, *, message: str) -> Optional[str]:
    cache_path = _ai_summary_cache_path(cache_dir, message=message)
    if not cache_path.exists() or not cache_path.is_file():
        return None
    try:
        obj = json.loads(cache_path.read_text(encoding="utf-8"))
    except Exception as e:
        logger.warning("AI summary cache read failed: %s (%s: %s)", cache_path, type(e).__name__, e)
        return None

    expected_hash = cache_path.stem
    if str(obj.get("message_sha256") or "") != expected_hash:
        return None
    if str(obj.get("prompt_message") or "") != str(message):
        return None

    text = str(obj.get("response_message") or "").strip()
    return text or None


def _write_ai_summary_cache(
    cache_dir: Path,
    *,
    message: str,
    response: Mapping[str, Any],
    section_key: str,
    section_label: str,
    model_label: str,
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = _ai_summary_cache_path(cache_dir, message=message)
    payload = {
        "schema_version": 1,
        "message_sha256": cache_path.stem,
        "prompt_message": str(message),
        "response": dict(response),
        "response_message": str(response.get("message") or ""),
        "section_key": str(section_key),
        "section_label": str(section_label),
        "model_label": str(model_label),
        "cached_at": datetime.now().astimezone().isoformat(),
    }
    tmp_path = cache_dir / f".{cache_path.stem}.{uuid.uuid4().hex}.tmp"
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(cache_path)
    return cache_path


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


def _pngquant_quality_for_item(
    *,
    item: PlotItem,
    default_quality: Optional[str],
    topdiff_quality: Optional[str],
) -> Optional[str]:
    if item.category == "topdiff-cluster" and topdiff_quality:
        return str(topdiff_quality)
    return str(default_quality) if default_quality else None


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
    out_root: Path,
    assets_dir: Optional[Path],
    interactive_resource_mode: str,
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
    interactive_protein_lookup: bool = False,
    exclude_dirs: Sequence[Path] = (),
) -> Tuple[List[VolcanoTableItem], Dict[str, Dict[str, Any]]]:
    root = base_dir.resolve()
    needles = [s for s in filter_contains if s]
    excluded = [p.resolve() for p in exclude_dirs if p]
    excluded.append((root / "report").resolve())

    tsvs = _collect_tsvs(root, patterns=[("volcano", "*.tsv")])
    out: List[VolcanoTableItem] = []
    resource_entries: Dict[str, Dict[str, Any]] = {}
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
        resource_id = ""
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
                        preview_html = _df_preview_table_html(
                            df,
                            max_rows=max_preview_rows,
                            interactive_protein_lookup=interactive_protein_lookup,
                        )
                    else:
                        error = fallback_error
                        try:
                            top_hits = picked.to_dict(orient="records")
                        except Exception:
                            top_hits = []
                        preview_html = _df_preview_table_html(
                            picked,
                            max_rows=max_preview_rows,
                            interactive_protein_lookup=interactive_protein_lookup,
                        )
                else:
                    try:
                        top_hits = picked.to_dict(orient="records")
                    except Exception:
                        top_hits = []
                    preview_html = _df_preview_table_html(
                        picked,
                        max_rows=max_preview_rows,
                        interactive_protein_lookup=interactive_protein_lookup,
                    )
            except Exception as e:
                error = str(e)
                preview_html = _df_preview_table_html(
                    df,
                    max_rows=max_preview_rows,
                    interactive_protein_lookup=interactive_protein_lookup,
                )

            resource_id = _make_resource_id("volcano-dataset", rel_str)
            payload_obj = _build_volcano_dataset_payload(
                df=df,
                source_relpath=rel_str,
                title=title,
                meta=meta,
                top_hits=top_hits,
                error=error,
            )
            entry, _ = _build_json_resource_entry(
                resource_id=resource_id,
                payload_obj=payload_obj,
                source_name=rel_str,
                kind="volcano_dataset",
                resource_mode=interactive_resource_mode,
                out_root=out_root,
                force=force,
            )
            resource_entries[resource_id] = entry
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
                resource_id=resource_id,
                top_hits=top_hits,
                error=error,
            )
        )

    out.sort(key=lambda x: x.source_relpath)
    return out, resource_entries


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
    interactive_payload: Optional[str] = None,
    interactive_protein_metadata: bool = False,
    interactive_resource_mode: str = "auto",
    defer_plot_images: bool = True,
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
    force_ai_summary: bool = False,
    pngquant: bool = False,
    pngquant_quality: Optional[str] = "65-85",
    pngquant_topdiff_quality: Optional[str] = "85-90",
    pngquant_speed: int = 3,
    pngquant_strip: bool = True,
) -> HtmlOverviewOutputs:
    root = Path(base_dir).expanduser().resolve()
    if not root.exists() or not root.is_dir():
        raise FileNotFoundError(str(root))
    if interactive_payload and interactive_protein_metadata:
        raise RuntimeError(
            "Pass either interactive_payload or interactive_protein_metadata, not both."
        )
    resource_mode = _effective_resource_mode(
        self_contained=bool(self_contained),
        interactive_resource_mode=str(interactive_resource_mode),
    )

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
        "HTML overview: base_dir=%s out_dir=%s copy_assets=%s force=%s self_contained=%s filter=%s pngquant=%s quality=%s topdiff_quality=%s speed=%s strip=%s",
        root,
        out_root,
        copy_assets,
        force,
        self_contained,
        list(filter_contains),
        pngquant,
        pngquant_quality,
        pngquant_topdiff_quality,
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
    optimized_count = 0
    total_before = 0
    total_after = 0
    resource_entries: Dict[str, Dict[str, Any]] = {}
    try:
        if self_contained:
            pbar = tqdm(plots, desc="compressing", disable=not use_pngquant)
            for item in pbar:
                source_png = Path(item.source_path)
                embed_png = source_png
                item.resource_id = _make_resource_id("plot-image", item.source_relpath)
                if use_pngquant and temp_png_dir is not None:
                    optimized_png = temp_png_dir / item.source_relpath
                    quality = _pngquant_quality_for_item(
                        item=item,
                        default_quality=pngquant_quality,
                        topdiff_quality=pngquant_topdiff_quality,
                    )
                    if _pngquant_optimize(
                        source_png,
                        optimized_png,
                        quality=quality,
                        speed=int(pngquant_speed),
                        strip=bool(pngquant_strip),
                    ):
                        before = source_png.stat().st_size
                        after = optimized_png.stat().st_size
                        total_before += before
                        total_after += after
                        optimized_count += 1
                        saved = total_before - total_after
                        ratio = (total_after / total_before) if total_before else 1.0
                        pbar.set_postfix(
                            saved=f"{saved / 1e6:.1f}MB",
                            ratio=f"{ratio:.2f}x",
                        )
                        embed_png = optimized_png
                plot_bytes = embed_png.read_bytes()
                plot_data_url = (
                    _data_uri_for_bytes(plot_bytes, mime="image/png")
                    if resource_mode == "inline"
                    else None
                )
                item.asset_relpath = ""
                resource_entries[item.resource_id] = _build_plot_resource_entry(
                    resource_id=item.resource_id,
                    item=item,
                    resource_mode=resource_mode,
                    data_url=plot_data_url,
                    path=None,
                    source_bytes=plot_bytes,
                )
            if hasattr(pbar, "close"):
                pbar.close()
        else:
            if assets_dir is None:
                raise ValueError("assets_dir is required when self_contained=False")
            pbar = tqdm(plots, desc="compressing", disable=not use_pngquant)
            for item in pbar:
                source_png = Path(item.source_path)
                dst = assets_dir / "plots" / item.source_relpath
                item.resource_id = _make_resource_id("plot-image", item.source_relpath)
                wrote_optimized = False
                if use_pngquant:
                    quality = _pngquant_quality_for_item(
                        item=item,
                        default_quality=pngquant_quality,
                        topdiff_quality=pngquant_topdiff_quality,
                    )
                    wrote_optimized = _pngquant_optimize(
                        source_png,
                        dst,
                        quality=quality,
                        speed=int(pngquant_speed),
                        strip=bool(pngquant_strip),
                    )
                    if wrote_optimized:
                        before = source_png.stat().st_size
                        after = dst.stat().st_size
                        total_before += before
                        total_after += after
                        optimized_count += 1
                        saved = total_before - total_after
                        ratio = (total_after / total_before) if total_before else 1.0
                        pbar.set_postfix(
                            saved=f"{saved / 1e6:.1f}MB",
                            ratio=f"{ratio:.2f}x",
                        )
                if not wrote_optimized:
                    _copy_or_symlink(source_png, dst, copy=copy_assets, force=force)
                item.asset_relpath = str(dst.relative_to(out_root).as_posix())
                plot_bytes = dst.read_bytes()
                plot_data_url = (
                    _data_uri_for_bytes(plot_bytes, mime="image/png")
                    if resource_mode == "inline"
                    else None
                )
                resource_entries[item.resource_id] = _build_plot_resource_entry(
                    resource_id=item.resource_id,
                    item=item,
                    resource_mode=resource_mode,
                    data_url=plot_data_url,
                    path=item.asset_relpath,
                    source_bytes=plot_bytes,
                )
            if hasattr(pbar, "close"):
                pbar.close()
    finally:
        if temp_png_dir_ctx is not None:
            temp_png_dir_ctx.cleanup()

    if use_pngquant:
        if optimized_count > 0 and total_before > 0:
            logger.info(
                "PNG compression summary: %.1f MB -> %.1f MB (%.2fx) [%d/%d optimized]",
                total_before / 1e6,
                total_after / 1e6,
                total_after / total_before,
                optimized_count,
                len(plots),
            )
        else:
            logger.info("PNG compression summary: no PNG files were optimized.")

    interactive_payload_ctx: Optional[Dict[str, Any]] = None
    plotly_resource_id: Optional[str] = None
    tabulator_resource_ids: Optional[Dict[str, str]] = None
    if interactive_payload:
        payload_path = Path(str(interactive_payload)).expanduser()
        if not payload_path.exists() or not payload_path.is_file():
            raise FileNotFoundError(str(payload_path))
        payload_obj = json.loads(payload_path.read_text(encoding="utf-8"))
        payload_resource_id = "interactive_payload"
        payload_entry, payload_encoded = _build_json_resource_entry(
            resource_id=payload_resource_id,
            payload_obj=payload_obj,
            source_name=payload_path.name,
            kind="custom",
            resource_mode=resource_mode,
            out_root=out_root,
            force=force,
        )
        resource_entries[payload_resource_id] = payload_entry
        interactive_payload_ctx = dict(payload_encoded)
        interactive_payload_ctx["kind"] = "custom"
    elif interactive_protein_metadata:
        payload_obj = _build_interactive_protein_metadata_payload(
            root=root,
            pandas_low_memory=bool(pandas_low_memory),
        )
        payload_resource_id = "interactive_payload"
        payload_entry, payload_encoded = _build_json_resource_entry(
            resource_id=payload_resource_id,
            payload_obj=payload_obj,
            source_name=str(payload_obj.get("source_table") or "protein_metadata.auto.json"),
            kind="protein_metadata",
            resource_mode=resource_mode,
            out_root=out_root,
            force=force,
        )
        resource_entries[payload_resource_id] = payload_entry
        interactive_payload_ctx = dict(payload_encoded)
        interactive_payload_ctx["kind"] = "protein_metadata"
        interactive_payload_ctx["protein_count"] = int(payload_obj.get("protein_count") or 0)
        interactive_payload_ctx["source_table"] = str(payload_obj.get("source_table") or "")
    if interactive_payload_ctx is not None:
        interactive_payload_ctx["resource_id"] = "interactive_payload"
        interactive_payload_ctx["resource_mode"] = resource_mode

    interactive_protein_lookup = bool(
        interactive_payload_ctx and interactive_payload_ctx.get("kind") == "protein_metadata"
    )

    # Volcano TSV summaries (copied into assets/data/...).
    volcano_tables, volcano_resource_entries = _collect_volcano_tables(
        base_dir=root,
        out_root=out_root,
        assets_dir=assets_dir,
        interactive_resource_mode=resource_mode,
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
        interactive_protein_lookup=interactive_protein_lookup,
        exclude_dirs=(out_root,),
    )
    resource_entries.update(volcano_resource_entries)

    if resource_mode == "chunked" and volcano_resource_entries:
        plotly_bundle = _find_plotly_bundle_path()
        if plotly_bundle is not None:
            plotly_resource_id = "plotly_bundle"
            resource_entries[plotly_resource_id] = _build_file_resource_entry(
                resource_id=plotly_resource_id,
                source_path=plotly_bundle,
                source_name=plotly_bundle.name,
                kind="script_bundle",
                format="javascript",
                mime="text/javascript",
                resource_mode=resource_mode,
                out_root=out_root,
                path="assets/interactive/vendor/plotly.min.js",
                copy_assets=copy_assets,
                force=force,
            )

        tabulator_bundle = _find_tabulator_bundle_paths()
        if tabulator_bundle is not None:
            tabulator_js_path, tabulator_css_path = tabulator_bundle
            tabulator_resource_ids = {
                "js": "tabulator_bundle_js",
                "css": "tabulator_bundle_css",
            }
            resource_entries[tabulator_resource_ids["js"]] = _build_file_resource_entry(
                resource_id=tabulator_resource_ids["js"],
                source_path=tabulator_js_path,
                source_name=tabulator_js_path.name,
                kind="script_bundle",
                format="javascript",
                mime="text/javascript",
                resource_mode=resource_mode,
                out_root=out_root,
                path="assets/interactive/vendor/tabulator.min.js",
                copy_assets=copy_assets,
                force=force,
            )
            resource_entries[tabulator_resource_ids["css"]] = _build_file_resource_entry(
                resource_id=tabulator_resource_ids["css"],
                source_path=tabulator_css_path,
                source_name=tabulator_css_path.name,
                kind="style_bundle",
                format="css",
                mime="text/css",
                resource_mode=resource_mode,
                out_root=out_root,
                path="assets/interactive/vendor/tabulator.min.css",
                copy_assets=copy_assets,
                force=force,
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
        no_agent_api_msg = (
            "AI summary requested but TACKLE_AGENT_API is not configured "
            "(set it in ~/.ispec/tackle-agent.conf or as an env var)."
        )
        cache_dir = _ai_summary_cache_dir(out_root)
        cfg = None
        support_chat_func = None
        session_base = None
        telemetry_client = None

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

        def _ensure_ai_client():
            nonlocal cfg, support_chat_func, session_base, telemetry_client
            if (
                cfg is not None
                and support_chat_func is not None
                and session_base is not None
                and telemetry_client is not None
            ):
                return cfg, support_chat_func, session_base, telemetry_client
            if not agent_api:
                raise RuntimeError(no_agent_api_msg)

            from .telemetry import (
                AgentTelemetry,
                AgentTelemetryConfig,
                make_local_events_path,
                support_chat,
            )

            cfg = AgentTelemetryConfig.from_env(
                agent_api=agent_api,
                local_events_path=make_local_events_path(root),
            )
            cfg = replace(
                cfg,
                timeout_seconds=_effective_html_ai_timeout_seconds(cfg.timeout_seconds),
            )
            support_chat_func = support_chat
            session_base = f"tackle-make-html:{uuid.uuid4().hex}"
            telemetry_client = AgentTelemetry(cfg, logger=logger)
            logger.info(
                "AI summary enabled for make-html: timeout=%.1fs sections=%d unmatched_tables=%d cache_dir=%s force_refresh=%s",
                cfg.timeout_seconds,
                len(sections),
                len(unmatched_tables),
                cache_dir,
                bool(force_ai_summary),
            )
            return cfg, support_chat_func, session_base, telemetry_client

        def _emit_ai_summary_event(
            *,
            phase: str,
            section_key: Optional[str] = None,
            section_label: Optional[str] = None,
            severity: str = "info",
            value: Optional[Dict[str, Any]] = None,
        ) -> None:
            try:
                _, _, session_base_local, telemetry_client_local = _ensure_ai_client()
            except Exception:
                return

            dims = {
                "analysis_outpath": str(root),
                "report_title": report_title,
                "section_key": section_key or None,
                "section_label": section_label or None,
            }
            dims = {key: val for key, val in dims.items() if val not in {None, ""}}
            try:
                telemetry_client_local.emit_event(
                    type=f"tackle.make_html.ai_summary.{phase}",
                    name=str(section_key or "make-html-ai-summary"),
                    severity=str(severity),
                    correlation_id=session_base_local,
                    dimensions=dims,
                    value=value or {},
                )
            except Exception as e:
                logger.debug(
                    "AI summary telemetry emit failed: phase=%s section=%s error=%r",
                    phase,
                    section_key,
                    e,
                )

        def _resolve_ai_summary(
            *,
            key: str,
            label: str,
            message: str,
            meta: Mapping[str, Any],
        ) -> Tuple[Optional[str], Optional[str]]:
            if not force_ai_summary:
                cached_text = _load_ai_summary_cache(cache_dir, message=message)
                if cached_text is not None:
                    logger.info("AI summary cache hit: section=%s", key)
                    _emit_ai_summary_event(
                        phase="section.complete",
                        section_key=key,
                        section_label=label,
                        value={"cached": True, "prompt_chars": len(message)},
                    )
                    return cached_text, None

            try:
                cfg_local, support_chat_local, session_base_local, _ = _ensure_ai_client()
                logger.info("AI summary request: section=%s", key)
                _emit_ai_summary_event(
                    phase="section.start",
                    section_key=key,
                    section_label=label,
                    value={"cached": False, "prompt_chars": len(message)},
                )
                resp = support_chat_local(
                    cfg_local,
                    session_id=f"{session_base_local}:{key}",
                    message=message,
                    meta=meta,
                )
                text = (resp.get("message") or "").strip()
                if not text:
                    _emit_ai_summary_event(
                        phase="section.failed",
                        section_key=key,
                        section_label=label,
                        severity="error",
                        value={"error": "empty_response"},
                    )
                    return None, "AI summary endpoint returned an empty response."
                _write_ai_summary_cache(
                    cache_dir,
                    message=message,
                    response=resp,
                    section_key=key,
                    section_label=label,
                    model_label=ai_model_label_value,
                )
                logger.info("AI summary complete: section=%s", key)
                _emit_ai_summary_event(
                    phase="section.complete",
                    section_key=key,
                    section_label=label,
                    value={
                        "cached": False,
                        "prompt_chars": len(message),
                        "response_chars": len(text),
                    },
                )
                return text, None
            except Exception as e:
                logger.warning("AI summary failed: section=%s error=%s: %s", key, type(e).__name__, e)
                _emit_ai_summary_event(
                    phase="section.failed",
                    section_key=key,
                    section_label=label,
                    severity="error",
                    value={"error": f"{type(e).__name__}: {e}"},
                )
                return None, f"{type(e).__name__}: {e}"

        _emit_ai_summary_event(
            phase="start",
            value={
                "sections": len(sections),
                "unmatched_volcano_tables": len(unmatched_tables),
                "force_refresh": bool(force_ai_summary),
                "cache_dir": str(cache_dir),
                "ai_model_label": ai_model_label_value,
            },
        )

        for sec in sections:
            key = str(sec.get("key") or "").strip()
            if not key:
                continue
            plot_paths = [p.source_relpath for p in sec.get("plots", [])][:30]
            prompt_obj: Dict[str, Any] = {
                "analysis_base_dir": str(root),
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

            msg = _build_ai_summary_prompt(
                section_key=key,
                section_label=str(sec.get("label") or key),
                prompt_obj=prompt_obj,
            )
            text, err = _resolve_ai_summary(
                key=key,
                label=str(sec.get("label") or key),
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
            if text:
                ai_section_summaries[key] = text
            elif err:
                ai_section_errors[key] = err

        if unmatched_tables:
            key = "volcano-tables"
            prompt_obj = {
                "analysis_base_dir": str(root),
                "section": {
                    "key": key,
                    "label": "Volcano Tables (unmatched)",
                    "count": len(unmatched_tables),
                },
                "volcano_unmatched": volcano_unmatched_items,
            }
            msg = _build_ai_summary_prompt(
                section_key=key,
                section_label="Volcano Tables (unmatched)",
                prompt_obj=prompt_obj,
                unmatched_volcano_tables=True,
            )
            text, err = _resolve_ai_summary(
                key=key,
                label="Volcano Tables (unmatched)",
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
            if text:
                ai_section_summaries[key] = text
            elif err:
                ai_section_errors[key] = err

        _emit_ai_summary_event(
            phase="complete",
            value={
                "sections_total": len(sections) + (1 if unmatched_tables else 0),
                "sections_succeeded": len(ai_section_summaries),
                "sections_failed": len(ai_section_errors),
            },
        )

    html_text = _render_j2_template(
        "overview_report.html.j2",
        title=report_title,
        generated_at=generated_at,
        show_date=bool(show_date),
        base_dir=(lambda _p: ("./" if _p == "." else (_p if _p.startswith(".") else f"./{_p}")))(os.path.relpath(root, Path.cwd().resolve()).replace("\\", "/")),
        self_contained=bool(self_contained),
        interactive_resource_mode=resource_mode,
        defer_plot_images=bool(defer_plot_images),
        interactive_payload=interactive_payload_ctx,
        plotly_resource_id=plotly_resource_id,
        tabulator_resource_ids=tabulator_resource_ids,
        resource_manifest={
            "mode": resource_mode,
            "entries": resource_entries,
            "resource_count": len(resource_entries),
        },
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
