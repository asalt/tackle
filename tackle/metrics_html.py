from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import pandas as pd

from .html_overview import (
    _clean_payload_scalar,
    _data_uri_for_file,
    _find_tabulator_bundle_paths,
    _render_j2_template,
    _slugify_html_id,
)
from .utils import _get_logger

logger = _get_logger(__name__)


def _safe_inline_script_text(text: str) -> str:
    return str(text).replace("</script", "<\\/script")


def _field_key(prefix: str, index: int, label: str) -> str:
    return f"{prefix}__{int(index)}__{_slugify_html_id(label)}"


def _humanize_plot_title(path: Path, *, export_stem: str) -> str:
    stem = path.stem
    export_name = Path(str(export_stem)).name
    if stem.startswith(export_name + "_"):
        stem = stem[len(export_name) + 1 :]
    label = stem.replace("__", " ").replace("_", " ").replace("-", " ").strip()
    while "  " in label:
        label = label.replace("  ", " ")
    return label.title() or path.name


def _clean_scalar(value: Any) -> Any:
    return _clean_payload_scalar(value, float_digits=2)


def _prepare_metadata_frame(
    metadata_df: Optional[pd.DataFrame],
    sample_names: Sequence[str],
) -> pd.DataFrame:
    if metadata_df is None or metadata_df.empty:
        return pd.DataFrame(index=pd.Index(sample_names, name="Sample"))

    meta = metadata_df.copy()
    meta.index = meta.index.map(str)
    meta = meta.reindex(list(sample_names))
    keep_cols: List[str] = []
    for col in meta.columns:
        values = [_clean_scalar(value) for value in meta[col].tolist()]
        if any(value not in (None, "") for value in values):
            keep_cols.append(str(col))
    if not keep_cols:
        return pd.DataFrame(index=pd.Index(sample_names, name="Sample"))
    return meta[keep_cols]


def _prepare_metrics_payload(
    *,
    metrics_df: pd.DataFrame,
    metadata_df: Optional[pd.DataFrame],
    miscut_df: Optional[pd.DataFrame],
) -> Dict[str, Any]:
    metrics_frame = metrics_df.copy()
    metrics_frame.index = metrics_frame.index.map(str)
    sample_names = list(metrics_frame.index)
    metadata_frame = _prepare_metadata_frame(metadata_df, sample_names)

    metrics_schema: List[Dict[str, Any]] = [
        {
            "field": "sample_name",
            "title": str(metrics_frame.index.name or "Sample"),
            "kind": "sample",
            "numeric": False,
            "frozen": True,
        }
    ]
    metadata_field_map: Dict[str, str] = {}
    for idx, col in enumerate(metadata_frame.columns):
        field = _field_key("meta", idx, str(col))
        metadata_field_map[str(col)] = field
        metrics_schema.append(
            {
                "field": field,
                "title": str(col),
                "kind": "metadata",
                "numeric": False,
                "frozen": False,
            }
        )

    metric_field_map: Dict[str, str] = {}
    for idx, col in enumerate(metrics_frame.columns):
        field = _field_key("metric", idx, str(col))
        metric_field_map[str(col)] = field
        metrics_schema.append(
            {
                "field": field,
                "title": str(col),
                "kind": "metric",
                "numeric": bool(pd.api.types.is_numeric_dtype(metrics_frame[col])),
                "frozen": False,
            }
        )

    rows: List[Dict[str, Any]] = []
    for sample in sample_names:
        row: Dict[str, Any] = {"sample_name": str(sample)}
        if sample in metadata_frame.index:
            meta_row = metadata_frame.loc[sample]
            for col in metadata_frame.columns:
                row[metadata_field_map[str(col)]] = _clean_scalar(meta_row.get(col))
        metric_row = metrics_frame.loc[sample]
        for col in metrics_frame.columns:
            row[metric_field_map[str(col)]] = _clean_scalar(metric_row.get(col))
        rows.append(row)

    miscut_payload: Optional[Dict[str, Any]] = None
    if miscut_df is not None and not miscut_df.empty:
        miscut_schema: List[Dict[str, Any]] = []
        miscut_field_map: Dict[str, str] = {}
        for idx, col in enumerate(miscut_df.columns):
            field = _field_key("miscut", idx, str(col))
            miscut_field_map[str(col)] = field
            miscut_schema.append(
                {
                    "field": field,
                    "title": str(col),
                    "kind": "metric",
                    "numeric": bool(pd.api.types.is_numeric_dtype(miscut_df[col])),
                    "frozen": idx == 0,
                }
            )
        miscut_rows: List[Dict[str, Any]] = []
        for row_obj in miscut_df.to_dict(orient="records"):
            row: Dict[str, Any] = {}
            for col in miscut_df.columns:
                row[miscut_field_map[str(col)]] = _clean_scalar(row_obj.get(col))
            miscut_rows.append(row)
        miscut_payload = {
            "rows": miscut_rows,
            "column_schema": miscut_schema,
        }

    def _summary_card(label: str, value: Any, detail: Optional[str] = None) -> Dict[str, str]:
        return {
            "label": str(label),
            "value": str(value),
            "detail": str(detail or ""),
        }

    summary_cards: List[Dict[str, str]] = [_summary_card("Samples", len(sample_names))]

    def _find_metric_column(candidates: Sequence[str]) -> Optional[str]:
        available = {str(col): col for col in metrics_frame.columns}
        for candidate in candidates:
            if candidate in available:
                return str(available[candidate])
        return None

    total_peptides_col = _find_metric_column(("Total_peptides", "PeptideCount", "peptide_count_total"))
    gpgroups_col = _find_metric_column(("GPGroups", "GPGroup", "GeneProducts"))

    if total_peptides_col:
        summary_cards.append(
            _summary_card(
                "Median Peptides",
                int(pd.to_numeric(metrics_frame[total_peptides_col], errors="coerce").fillna(0).median()),
                "per sample",
            )
        )
    if gpgroups_col:
        summary_cards.append(
            _summary_card(
                "Median GP Groups",
                int(pd.to_numeric(metrics_frame[gpgroups_col], errors="coerce").fillna(0).median()),
                "per sample",
            )
        )

    return {
        "metrics": {
            "rows": rows,
            "column_schema": metrics_schema,
        },
        "miscut": miscut_payload,
        "summary_cards": summary_cards,
    }


def _render_static_html_table(column_schema: Sequence[Mapping[str, Any]], rows: Sequence[Mapping[str, Any]]) -> str:
    pieces = ["<table class='static-table'><thead><tr>"]
    for col in column_schema:
        pieces.append(f"<th>{html.escape(str(col.get('title') or ''))}</th>")
    pieces.append("</tr></thead><tbody>")
    for row in rows:
        pieces.append("<tr>")
        for col in column_schema:
            value = row.get(str(col.get("field") or ""), "")
            pieces.append(f"<td>{html.escape('' if value is None else str(value))}</td>")
        pieces.append("</tr>")
    pieces.append("</tbody></table>")
    return "".join(pieces)


def build_metrics_html_report(
    *,
    out_html: str,
    metrics_df: pd.DataFrame,
    metadata_df: Optional[pd.DataFrame] = None,
    miscut_df: Optional[pd.DataFrame] = None,
    plot_paths: Sequence[Path] = (),
    export_stem: str = "metrics",
    title: Optional[str] = None,
    analysis_label: Optional[str] = None,
    before_filter: bool = False,
    before_norm: bool = False,
) -> Path:
    out_path = Path(out_html).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = _prepare_metrics_payload(
        metrics_df=metrics_df,
        metadata_df=metadata_df,
        miscut_df=miscut_df,
    )

    tabulator_css = None
    tabulator_js = None
    tabulator_available = False
    bundle = _find_tabulator_bundle_paths()
    if bundle is not None:
        try:
            tabulator_js_path, tabulator_css_path = bundle
            tabulator_css = tabulator_css_path.read_text(encoding="utf-8")
            tabulator_js = _safe_inline_script_text(
                tabulator_js_path.read_text(encoding="utf-8")
            )
            tabulator_available = True
        except Exception as exc:
            logger.warning("Metrics HTML: unable to inline Tabulator assets: %s", exc)
            tabulator_css = None
            tabulator_js = None
            tabulator_available = False

    plot_items: List[Dict[str, str]] = []
    for raw_path in plot_paths:
        path = Path(raw_path)
        if not path.exists() or not path.is_file() or path.suffix.lower() != ".png":
            continue
        plot_items.append(
            {
                "title": _humanize_plot_title(path, export_stem=export_stem),
                "source_name": path.name,
                "data_url": _data_uri_for_file(path, mime="image/png"),
            }
        )

    report_title = title or f"Tackle Metrics Report: {analysis_label or Path(export_stem).name}"
    context_pills = [
        {
            "label": "Filter",
            "value": "Before filtering" if before_filter else "After filtering",
        },
        {
            "label": "Area",
            "value": "Before normalization" if before_norm else "After normalization",
        },
        {
            "label": "Generated",
            "value": datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
        },
    ]

    html_text = _render_j2_template(
        "metrics_report.html.j2",
        title=report_title,
        analysis_label=analysis_label or Path(export_stem).name,
        export_stem=Path(export_stem).name,
        context_pills=context_pills,
        tabulator_available=tabulator_available,
        tabulator_css=tabulator_css,
        tabulator_js=tabulator_js,
        summary_cards=payload["summary_cards"],
        metrics_payload=payload["metrics"],
        miscut_payload=payload["miscut"],
        plots=plot_items,
        metrics_download_name=f"{Path(export_stem).name}.csv",
        static_metrics_table_html=_render_static_html_table(
            payload["metrics"]["column_schema"],
            payload["metrics"]["rows"],
        ),
        static_miscut_table_html=(
            _render_static_html_table(
                payload["miscut"]["column_schema"],
                payload["miscut"]["rows"],
            )
            if payload["miscut"] is not None
            else ""
        ),
    )
    out_path.write_text(html_text, encoding="utf-8")
    return out_path
