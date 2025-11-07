import os
import re
from pathlib import Path
from typing import Iterable, List, Tuple, Optional, Dict, Any, Sequence

import pandas as pd

from .utils import _get_logger

logger = _get_logger(__name__)

def _sanitize_sheet_name(name: str) -> str:
    # Excel sheet names: max 31 chars; cannot contain: : \ / ? * [ ]
    invalid = r'[:\\/?*\[\]]'
    safe = re.sub(invalid, "_", str(name))
    safe = safe.strip()
    # Avoid empty names
    if not safe:
        safe = "Sheet"
    # Excel also dislikes names that end with a single quote
    if safe.endswith("'"):
        safe = safe[:-1]
    if len(safe) > 31:
        safe = safe[:31]
    return safe


def _sanitize_label(value: str) -> str:
    """Sanitize a label for use as a column prefix without truncation.

    Removes characters that Excel sheet headers may not like but does not enforce
    a 31-character limit (which only applies to sheet names, not headers).
    """
    invalid = r'[:\\/?*\[\]]'
    safe = re.sub(invalid, "_", str(value))
    # Normalize whitespace to underscores
    safe = re.sub(r"\s+", "_", safe)
    return safe.strip("_")

def _collect_tsvs(root: Path, patterns: Iterable[Tuple[str, str]]) -> List[Tuple[str, Path]]:
    logger.info("Exporter: scanning root=%s for patterns=%s", root, list(patterns))
    files: List[Tuple[str, Path]] = []
    for subdir, pat in patterns:
        # Prefer the canonical subdir under the provided root, but fall back to
        # recursively searching for any matching subdir deeper in the tree.
        search_dirs: List[Path] = []
        if subdir:
            d = root / subdir
            if d.exists():
                search_dirs = [d]
            else:
                # Recursively discover candidate subdirs named like 'export'/'volcano'
                search_dirs = [p for p in root.rglob(subdir) if p.is_dir()]
            logger.info(
                "Exporter: searching %d directories for '%s' files: %s",
                len(search_dirs), pat, ", ".join(str(p) for p in search_dirs) or "<none>"
            )
        else:
            search_dirs = [root]

        for d in search_dirs:
            for p in d.rglob(pat):
                try:
                    rel = p.relative_to(root)
                except ValueError:
                    # If p is outside root due to symlinks, fallback to basename
                    rel = Path(subdir) / p.name if subdir else p.name
                files.append((str(rel), p))
    logger.info("Exporter: discovered %d TSV files", len(files))
    # Stable order
    files.sort()
    return files


def build_export_xlsx(
    base_dir: str,
    out_path: str,
    include_export: bool = True,
    include_volcano: bool = True,
    pheno_df: Optional[pd.DataFrame] = None,
    meta: Optional[Dict[str, Any]] = None,
    *,
    engine_preference: Optional[str] = None,
    filter_contains: Optional[Sequence[str]] = None,
    slim_volcano: bool = True,
    volcano_topn: Optional[int] = None,
    merge_volcano: bool = False,
) -> str:
    """
    Create an Excel workbook summarizing common tackle exports under an analysis directory.

    Parameters
    ----------
    base_dir: Root analysis output directory (e.g., data_obj.outpath)
    out_path: Destination .xlsx path
    include_export: Include tables from the 'export/' subdir
    include_volcano: Include volcano result TSVs
    """
    root = Path(base_dir)
    logger.info(
        "Exporter: build_export_xlsx base_dir=%s out_path=%s include_export=%s include_volcano=%s",
        root, out_path, include_export, include_volcano,
    )
    patterns: List[Tuple[str, str]] = []
    if include_export:
        patterns.append(("export", "*.tsv"))
    if include_volcano:
        patterns.append(("volcano", "*.tsv"))

    targets = _collect_tsvs(root, patterns)
    if filter_contains:
        fc = [s for s in filter_contains if s]
        if fc:
            logger.info("Exporter: filtering %d files by substrings: %s", len(targets), ", ".join(fc))
            targets = [(rel, p) for (rel, p) in targets if any(s in str(rel) for s in fc)]
    if not targets:
        logger.error("Exporter: no TSVs found under %s for patterns %s", base_dir, patterns)
        raise FileNotFoundError(f"No TSV exports found under {base_dir} matching \n{patterns}")

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Choose an engine if available
    # Prefer xlsxwriter for faster writes when available
    engine = None
    candidates = ("xlsxwriter", "openpyxl")
    if engine_preference:
        candidates = (engine_preference,) + tuple(c for c in candidates if c != engine_preference)
    for candidate in candidates:
        try:
            __import__(candidate)
            engine = candidate
            break
        except Exception:
            continue
    if engine is None:
        logger.error("Exporter: no Excel writer engine available (need openpyxl or xlsxwriter)")
        raise RuntimeError(
            "No Excel writer engine found. Install one of 'openpyxl' or 'xlsxwriter' to create XLSX exports."
        )
    logger.info("Exporter: using Excel engine '%s'", engine)

    with pd.ExcelWriter(out, engine=engine) as xw:
        sources_rows: List[Dict[str, Any]] = []
        # Manifest sheet removed per request; key metadata remains in logs/context
        # Optional phenotype/metadata sheet first
        if pheno_df is not None and isinstance(pheno_df, pd.DataFrame) and not pheno_df.empty:
            pheno = pheno_df.copy()
            pheno.index.name = pheno.index.name or "Sample"
            pheno_reset = pheno.reset_index()
            pheno_sheet = _sanitize_sheet_name("phenotype")
            logger.info("Exporter: adding sheet '%s' with shape %s", pheno_sheet, pheno_reset.shape)
            pheno_reset.to_excel(xw, sheet_name=pheno_sheet, index=False)
        # Optionally merge all volcano TSVs into a single wide sheet
        volcano_targets = [(rel, p) for (rel, p) in targets if str(rel).startswith("volcano")]
        export_targets = [(rel, p) for (rel, p) in targets if str(rel).startswith("export")]

        def _label_from_rel(relpath: str) -> str:
            base = os.path.basename(relpath)
            name = re.sub(r"\.tsv$", "", base, flags=re.IGNORECASE)
            m = re.search(r"(?<=group)[_]?(.*)$", name)
            lbl = m.group(1) if m else name
            # Use a generous, sanitized label for column prefixes (no 31-char cap needed)
            return _sanitize_label(lbl)

        if merge_volcano and volcano_targets:
            logger.info("Exporter: merging %d volcano tables into a single wide sheet", len(volcano_targets))
            merged: Optional[pd.DataFrame] = None
            seen_labels = set()
            annot_df: Optional[pd.DataFrame] = None  # index = GeneID; columns: FunCats, GeneDescription
            keep_cols_base = [
                "GeneID",
                "GeneSymbol",
                "log2_FC",
                "pAdj",
                "pValue",
                "signedlogP",
                "t",
            ]
            for idx, (rel, path) in enumerate(volcano_targets, start=1):
                try:
                    df = pd.read_csv(path, sep="\t")
                except Exception:
                    df = pd.read_table(path)
                df.columns = [str(c) for c in df.columns]
                if "GeneID" not in df.columns:
                    logger.warning("Exporter: skipping volcano file without GeneID column: %s", rel)
                    continue
                df["GeneID"] = df["GeneID"].astype(str)
                # Capture gene annotations for front matter (once, coalescing across files)
                try:
                    cand = {}
                    if "FunCats" in df.columns:
                        cand["FunCats"] = df["FunCats"].astype(str)
                    if "GeneDescription" in df.columns:
                        cand["GeneDescription"] = df["GeneDescription"].astype(str)
                    elif "Description" in df.columns:
                        cand["GeneDescription"] = df["Description"].astype(str)
                    if cand:
                        ad = pd.DataFrame(cand)
                        ad.insert(0, "GeneID", df["GeneID"].astype(str))
                        ad = ad.drop_duplicates(subset=["GeneID"]).set_index("GeneID")
                        annot_df = ad if annot_df is None else annot_df.combine_first(ad)
                except Exception:
                    logger.debug("Exporter: unable to capture annotations from %s", rel)
                label = _label_from_rel(rel)
                # Ensure uniqueness
                if label in seen_labels:
                    i = 2
                    while f"{label}_{i}" in seen_labels:
                        i += 1
                    label = f"{label}_{i}"
                seen_labels.add(label)
                # Select columns
                cols = [c for c in keep_cols_base if c in df.columns]
                slim_df = df[cols] if slim_volcano and cols else df
                if slim_volcano and cols:
                    logger.info("Exporter: slim-volcano applied to '%s' (columns kept: %s)", rel, cols)
                # Optionally reduce to top-N rows for speed
                if volcano_topn:
                    if "pAdj" in slim_df.columns:
                        slim_df = slim_df.sort_values("pAdj").head(int(volcano_topn))
                    elif "pValue" in slim_df.columns:
                        slim_df = slim_df.sort_values("pValue").head(int(volcano_topn))
                    else:
                        slim_df = slim_df.head(int(volcano_topn))
                # Rename metric columns with label prefix
                rename_map = {}
                for c in slim_df.columns:
                    if c in ("GeneID", "GeneSymbol"):
                        continue
                    rename_map[c] = f"{label}__{c}"
                slim_df = slim_df.rename(columns=rename_map)
                # Merge outer on GeneID; preserve GeneSymbol if present
                if merged is None:
                    merged = slim_df.copy()
                else:
                    merged = pd.merge(merged, slim_df, on=["GeneID", *(["GeneSymbol"] if "GeneSymbol" in merged.columns and "GeneSymbol" in slim_df.columns else [])], how="outer")
                logger.info("Exporter: merged '%s' -> current shape %s", rel, merged.shape if merged is not None else None)
                # Track sources
                try:
                    stat = path.stat(); mtime = stat.st_mtime
                except Exception:
                    mtime = None
                sources_rows.append({
                    "sheet": "volcano_merged",
                    "relative_path": str(rel),
                    "rows": int(df.shape[0]),
                    "cols": int(df.shape[1]),
                    "mtime": mtime,
                    "kind": "volcano",
                })
            if merged is not None:
                # Attach front-matter annotations if available
                if annot_df is not None:
                    try:
                        merged = merged.merge(annot_df.reset_index(), on="GeneID", how="left")
                    except Exception:
                        logger.debug("Exporter: failed to merge annotations onto merged volcano table")
                # Order columns: front matter then metrics grouped by type
                front = [c for c in ("GeneID", "GeneSymbol", "FunCats", "GeneDescription") if c in merged.columns]
                # Discover metric columns with pattern '<label>__<metric>'
                metric_cols = [c for c in merged.columns if c not in front]
                label_metric = {}
                for c in metric_cols:
                    m = re.match(r"(.+?)__([A-Za-z0-9_.]+)$", c)
                    if not m:
                        continue
                    lbl, met = m.group(1), m.group(2)
                    label_metric.setdefault(met, []).append((lbl, c))
                # Desired metric ordering
                metric_order = ["log2_FC", "signedlogP", "pValue", "pAdj", "t"]
                ordered_metric_cols: List[str] = []
                for met in metric_order:
                    pairs = label_metric.get(met, [])
                    # sort by label for deterministic order
                    for _lbl, col in sorted(pairs, key=lambda x: x[0]):
                        ordered_metric_cols.append(col)
                # Append any leftover metric columns not matched above
                leftover = [c for c in metric_cols if c not in ordered_metric_cols]
                ordered_cols = front + ordered_metric_cols + leftover
                merged = merged[ordered_cols]
                sheet_name = _sanitize_sheet_name("volcano_merged")
                logger.info("Exporter: adding merged volcano sheet '%s' with shape %s", sheet_name, merged.shape)
                merged.to_excel(xw, sheet_name=sheet_name, index=False)

        # Write non-volcano (export) tables as individual sheets, and volcano tables
        # individually only if not merging
        remaining = export_targets + ([] if merge_volcano else volcano_targets)
        total = len(remaining)
        for idx, (rel, path) in enumerate(remaining, start=1):
            try:
                df = pd.read_csv(path, sep="\t")
            except Exception:
                # Fallback: try without dtype inference pitfalls
                df = pd.read_table(path)
            # Deduce a concise sheet name from relative path
            base = str(rel)
            base = base.replace(os.sep, "__")
            base = re.sub(r"\.tsv$", "", base, flags=re.IGNORECASE)
            sheet = _sanitize_sheet_name(base)
            # Handle duplicates by appending an index
            suffix = 1
            final = sheet
            while final in xw.sheets:
                suffix += 1
                final = _sanitize_sheet_name(f"{sheet}_{suffix}")
            orig_shape = df.shape
            # Optional volcano slimming: keep only summary columns and/or top-N
            if slim_volcano and str(rel).startswith("volcano"):
                keep_cols = [
                    "GeneID",
                    "GeneSymbol",
                    "log2_FC",
                    "pAdj",
                    "pValue",
                    "signedlogP",
                    "t",
                ]
                cols = [c for c in keep_cols if c in df.columns]
                if cols:
                    df = df[cols]
            if volcano_topn and str(rel).startswith("volcano"):
                # If sorting columns exist, use them; otherwise just head()
                if "pAdj" in df.columns:
                    df = df.sort_values("pAdj", ascending=True).head(int(volcano_topn))
                elif "pValue" in df.columns:
                    df = df.sort_values("pValue", ascending=True).head(int(volcano_topn))
                else:
                    df = df.head(int(volcano_topn))
            logger.info(
                "Exporter: [%d/%d] adding sheet '%s' from '%s' with shape %s (orig %s)",
                idx, total, final, rel, df.shape, orig_shape,
            )
            df.to_excel(xw, sheet_name=final, index=False)

            try:
                stat = path.stat()
                mtime = stat.st_mtime
            except Exception:
                mtime = None
            sources_rows.append({
                "sheet": final,
                "relative_path": str(rel),
                "rows": int(df.shape[0]),
                "cols": int(df.shape[1]),
                "mtime": mtime,
                "kind": "volcano" if str(rel).startswith("volcano") else ("export" if str(rel).startswith("export") else "other"),
            })

        # Add a 'sources' sheet last with a summary of included TSVs
        try:
            if sources_rows:
                sources_df = pd.DataFrame(sources_rows)
                logger.info("Exporter: adding 'sources' sheet with %d entries", len(sources_df))
                sources_df.to_excel(xw, sheet_name=_sanitize_sheet_name("sources"), index=False)
        except Exception:
            logger.exception("Exporter: failed to add sources sheet")

    logger.info("Exporter: wrote workbook %s", out)
    return str(out)
