#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from http.cookiejar import CookieJar
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import HTTPCookieProcessor, Request, build_opener, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
ANNOTATION_DIR = ROOT / "tackle" / "data" / "annotations"
LEGACY_TABLE = ANNOTATION_DIR / "combined_annotations_new.tsv"
SOURCE_CACHE = ANNOTATION_DIR / "source_cache"
GENERATED_DIR = ANNOTATION_DIR / "generated"

HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
SURFACEOME_URL = "https://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"
SURFACEOME_IDS_URL = "https://wlab.ethz.ch/surfaceome/surfaceome_ids.txt"
MATRISOME_URL = "https://matrisomedb.org/"
MATRISOME_SEARCH_URL = "https://matrisomedb.org/pats_search/"

MATRISOME_CORE = ("Collagens", "ECM Glycoproteins", "Proteoglycans")
MATRISOME_ASSOCIATED = (
    "ECM-affiliated Proteins",
    "ECM Regulators",
    "Secreted Factors",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def request_bytes(url: str, *, data: bytes | None = None, headers: dict | None = None) -> bytes:
    req_headers = {
        "User-Agent": "tackle-annotation-builder/0.1",
    }
    if headers:
        req_headers.update(headers)
    req = Request(url, data=data, headers=req_headers)
    with urlopen(req, timeout=120) as response:
        return response.read()


def download(url: str, dest: Path, *, force: bool = False) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    status = "cached"
    if force or not dest.exists():
        status = "downloaded"
        dest.write_bytes(request_bytes(url))
    return {
        "url": url,
        "path": str(dest.relative_to(ROOT)),
        "status": status,
        "sha256": sha256_file(dest),
        "bytes": dest.stat().st_size,
    }


def norm_symbol(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def split_symbols(value) -> list[str]:
    text = norm_symbol(value)
    if not text:
        return []
    tokens = re.split(r"[;,]\s*|\s+or\s+|\s+and\s+", text)
    out = []
    for token in tokens:
        token = token.strip()
        if not token or token.upper() in {"NA", "N/A", "NONE"}:
            continue
        out.append(token)
    return out


def positive_mask(series: pd.Series) -> pd.Series:
    values = series.fillna("").astype(str).str.strip()
    lower = values.str.lower()
    return (values != "") & (~lower.isin({"0", "false", "na", "nan", "none"}))


def add_values_by_symbol(
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    symbols: Iterable[str],
    column: str,
    value: str,
) -> int:
    hits = 0
    for symbol in symbols:
        for geneid in symbol_to_geneids.get(norm_symbol(symbol), ()):
            generated.loc[generated["GeneID"] == geneid, column] = value
            hits += 1
    return hits


def load_legacy() -> pd.DataFrame:
    legacy = pd.read_csv(LEGACY_TABLE, sep="\t", dtype=str).fillna("")
    legacy["GeneID"] = legacy["GeneID"].astype(str)
    legacy["GeneSymbol"] = legacy["GeneSymbol"].astype(str)
    return legacy


def read_hpa(cache: Path, manifest: dict, *, force: bool) -> pd.DataFrame | None:
    try:
        info = download(HPA_URL, cache / "proteinatlas.tsv.zip", force=force)
        with zipfile.ZipFile(cache / "proteinatlas.tsv.zip") as zf:
            members = zf.namelist()
            member = next(x for x in members if x.endswith(".tsv"))
            with zf.open(member) as handle:
                df = pd.read_csv(handle, sep="\t", dtype=str).fillna("")
        info.update({"rows": int(len(df)), "columns": list(df.columns)})
        manifest["sources"]["hpa"] = info
        return df
    except Exception as exc:
        manifest["sources"]["hpa"] = {"url": HPA_URL, "status": "error", "error": str(exc)}
        return None


def first_existing_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lookup = {c.lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def apply_hpa_annotations(
    hpa: pd.DataFrame | None,
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    manifest: dict,
) -> None:
    if hpa is None:
        return
    symbol_col = first_existing_column(hpa, ("Gene name", "Gene", "GeneSymbol"))
    if symbol_col is None:
        manifest["transforms"]["hpa"] = {"status": "error", "error": "no gene symbol column"}
        return

    subcell_cols = [
        c
        for c in hpa.columns
        if "subcellular" in c.lower() or "location" in c.lower()
    ]
    secret_cols = [
        c
        for c in hpa.columns
        if "secretome" in c.lower() or "secreted" in c.lower()
    ]

    stats = {
        "symbol_column": symbol_col,
        "subcellular_columns": subcell_cols,
        "secretome_columns": secret_cols,
        "er_golgi_source_rows": 0,
        "secreted_source_rows": 0,
    }

    for _, row in hpa.iterrows():
        symbols = split_symbols(row.get(symbol_col, ""))
        if not symbols:
            continue
        loc_text = " ; ".join(str(row.get(c, "")) for c in subcell_cols).lower()
        has_er = "endoplasmic reticulum" in loc_text or re.search(r"\ber\b", loc_text)
        has_golgi = "golgi" in loc_text
        if has_er or has_golgi:
            if has_er and has_golgi:
                value = "BOTH"
            elif has_er:
                value = "ENDORETICULUM"
            else:
                value = "GOLGI"
            add_values_by_symbol(generated, symbol_to_geneids, symbols, "ER_GOLGI", value)
            stats["er_golgi_source_rows"] += 1

        secret_text = " ; ".join(str(row.get(c, "")) for c in secret_cols).lower()
        if "secret" in secret_text:
            add_values_by_symbol(generated, symbol_to_geneids, symbols, "SECRETED", "SECRETED")
            add_values_by_symbol(generated, symbol_to_geneids, symbols, "Secreted", "Secreted")
            stats["secreted_source_rows"] += 1

    manifest["transforms"]["hpa"] = stats


def read_surfaceome(cache: Path, manifest: dict, *, force: bool) -> pd.DataFrame | None:
    try:
        info = download(SURFACEOME_URL, cache / "table_S3_surfaceome.xlsx", force=force)
        ids_info = download(SURFACEOME_IDS_URL, cache / "surfaceome_ids.txt", force=force)
        df = pd.read_excel(
            cache / "table_S3_surfaceome.xlsx", dtype=str, header=1
        ).fillna("")
        info.update(
            {
                "rows": int(len(df)),
                "columns": list(df.columns),
                "ids_file": ids_info,
            }
        )
        manifest["sources"]["surfaceome"] = info
        return df
    except Exception as exc:
        manifest["sources"]["surfaceome"] = {
            "url": SURFACEOME_URL,
            "status": "error",
            "error": str(exc),
        }
        return None


def apply_surfaceome_annotations(
    surfaceome: pd.DataFrame | None,
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    manifest: dict,
) -> None:
    if surfaceome is None:
        return
    symbol_col = first_existing_column(
        surfaceome,
        (
            "UniProt gene",
            "Gene",
            "Gene name",
            "Gene names",
            "HGNC symbol",
            "Symbol",
        ),
    )
    label_col = first_existing_column(surfaceome, ("Surfaceome Label", "SURFY", "label"))
    uniprot_subcellular_col = first_existing_column(
        surfaceome, ("UniProt subcellular", "Subcellular location")
    )
    stats = {
        "symbol_column": symbol_col,
        "label_column": label_col,
        "uniprot_subcellular_column": uniprot_subcellular_col,
        "surface_source_rows": 0,
        "nonsurface_source_rows": 0,
        "cell_membrane_source_rows": 0,
    }
    if symbol_col is None:
        stats.update({"status": "error", "error": "no gene symbol column"})
        manifest["transforms"]["surfaceome"] = stats
        return

    for _, row in surfaceome.iterrows():
        symbols = split_symbols(row.get(symbol_col, ""))
        if not symbols:
            continue
        label = str(row.get(label_col, "") if label_col else "").strip().lower()
        if label == "surface" or label in {"1", "true"}:
            add_values_by_symbol(generated, symbol_to_geneids, symbols, "SurfaceLabel", "surface")
            stats["surface_source_rows"] += 1
        elif label == "non-surface" or label == "nonsurface":
            add_values_by_symbol(
                generated, symbol_to_geneids, symbols, "SurfaceLabel", "nonsurface"
            )
            stats["nonsurface_source_rows"] += 1

        subcellular = str(row.get(uniprot_subcellular_col, "") if uniprot_subcellular_col else "")
        if "cell membrane" in subcellular.lower():
            add_values_by_symbol(
                generated, symbol_to_geneids, symbols, "CellMembrane", "Cell Membrane"
            )
            stats["cell_membrane_source_rows"] += 1

    manifest["transforms"]["surfaceome"] = stats


def matrisome_opener():
    cj = CookieJar()
    return build_opener(HTTPCookieProcessor(cj))


def fetch_matrisome_category(opener, csrf: str, category: str) -> list[dict]:
    data = [
        ("csrfmiddlewaretoken", csrf),
        ("seq", ""),
        ("category[]", category),
        ("species[]", "Human"),
    ]
    body = urlencode(data).encode()
    req = Request(
        MATRISOME_SEARCH_URL,
        data=body,
        headers={
            "User-Agent": "tackle-annotation-builder/0.1",
            "Content-Type": "application/x-www-form-urlencoded",
            "X-Requested-With": "XMLHttpRequest",
            "X-CSRFToken": csrf,
            "Referer": MATRISOME_URL,
        },
    )
    with opener.open(req, timeout=120) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return payload.get("data", [])


def read_matrisome(cache: Path, manifest: dict, *, force: bool) -> dict[str, pd.DataFrame]:
    cache.mkdir(parents=True, exist_ok=True)
    out: dict[str, pd.DataFrame] = {}
    source_info = {"url": MATRISOME_URL, "status": "ok", "categories": {}}
    try:
        opener = matrisome_opener()
        html = opener.open(
            Request(MATRISOME_URL, headers={"User-Agent": "tackle-annotation-builder/0.1"}),
            timeout=120,
        ).read().decode("utf-8", errors="ignore")
        match = re.search(r"name='csrfmiddlewaretoken' value='([^']+)'", html)
        if not match:
            raise RuntimeError("CSRF token not found on MatrisomeDB page")
        csrf = match.group(1)
        for category in (*MATRISOME_CORE, *MATRISOME_ASSOCIATED):
            safe = re.sub(r"[^A-Za-z0-9]+", "_", category).strip("_")
            path = cache / f"matrisomedb_{safe}.json"
            status = "cached"
            if force or not path.exists():
                rows = fetch_matrisome_category(opener, csrf, category)
                path.write_text(json.dumps(rows, indent=2, sort_keys=True), encoding="utf-8")
                status = "downloaded"
            rows = json.loads(path.read_text(encoding="utf-8"))
            df = pd.DataFrame(rows).fillna("")
            out[category] = df
            source_info["categories"][category] = {
                "path": str(path.relative_to(ROOT)),
                "status": status,
                "rows": int(len(df)),
                "sha256": sha256_file(path),
            }
    except Exception as exc:
        source_info = {"url": MATRISOME_URL, "status": "error", "error": str(exc)}
    manifest["sources"]["matrisome"] = source_info
    return out


def apply_matrisome_annotations(
    matrisome: dict[str, pd.DataFrame],
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    manifest: dict,
) -> None:
    stats = {"core_source_rows": 0, "associated_source_rows": 0}
    for category, df in matrisome.items():
        if df.empty or "gene" not in df.columns:
            continue
        value = "CORE" if category in MATRISOME_CORE else "ASSOCIATED"
        for symbol in df["gene"].astype(str).dropna().unique():
            add_values_by_symbol(generated, symbol_to_geneids, [symbol], "MATRISOME", value)
        if value == "CORE":
            stats["core_source_rows"] += int(len(df))
        else:
            stats["associated_source_rows"] += int(len(df))
    manifest["transforms"]["matrisome"] = stats


def make_comparison(legacy: pd.DataFrame, generated: pd.DataFrame, columns: list[str]):
    summary_rows = []
    detail_frames = []
    merged = legacy[["GeneID", "GeneSymbol", *[c for c in columns if c in legacy]]].merge(
        generated[["GeneID", *columns]],
        on="GeneID",
        how="outer",
        suffixes=("_legacy", "_generated"),
    ).fillna("")

    if "GeneSymbol" not in merged:
        merged["GeneSymbol"] = ""

    for column in columns:
        lcol = f"{column}_legacy"
        gcol = f"{column}_generated"
        if lcol not in merged:
            merged[lcol] = ""
        if gcol not in merged:
            merged[gcol] = ""
        legacy_pos = positive_mask(merged[lcol])
        generated_pos = positive_mask(merged[gcol])
        overlap = legacy_pos & generated_pos
        legacy_only = legacy_pos & ~generated_pos
        generated_only = generated_pos & ~legacy_pos
        denom = int((legacy_pos | generated_pos).sum())
        summary_rows.append(
            {
                "column": column,
                "legacy_positive": int(legacy_pos.sum()),
                "generated_positive": int(generated_pos.sum()),
                "overlap": int(overlap.sum()),
                "legacy_only": int(legacy_only.sum()),
                "generated_only": int(generated_only.sum()),
                "jaccard": round(float(overlap.sum() / denom), 4) if denom else 1.0,
            }
        )
        details = merged.loc[
            legacy_only | generated_only,
            ["GeneID", "GeneSymbol", lcol, gcol],
        ].copy()
        details.insert(0, "column", column)
        details = details.rename(columns={lcol: "legacy_value", gcol: "generated_value"})
        detail_frames.append(details)

    summary = pd.DataFrame(summary_rows)
    details = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    return summary, details


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild selected annotation columns from public sources and compare to the legacy manual table."
    )
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--outdir", type=Path, default=GENERATED_DIR)
    parser.add_argument("--cache-dir", type=Path, default=SOURCE_CACHE)
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    legacy = load_legacy()
    generated = legacy[["GeneID", "GeneSymbol"]].copy()
    for column in (
        "ER_GOLGI",
        "MATRISOME",
        "SECRETED",
        "Secreted",
        "SurfaceLabel",
        "CellMembrane",
    ):
        generated[column] = ""

    symbol_to_geneids: dict[str, set[str]] = defaultdict(set)
    for _, row in legacy.iterrows():
        symbol = norm_symbol(row.get("GeneSymbol"))
        geneid = norm_symbol(row.get("GeneID"))
        if symbol and geneid:
            symbol_to_geneids[symbol].add(geneid)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "legacy_table": str(LEGACY_TABLE.relative_to(ROOT)),
        "sources": {},
        "transforms": {},
        "citations": {
            "matrisome": {
                "label": "MatrisomeDB",
                "pmid": "31586405",
                "doi": "10.1093/nar/gkz849",
            },
            "matrisome_2_0": {
                "label": "MatrisomeDB 2.0",
                "pmid": "36399478",
                "doi": "10.1093/nar/gkac1009",
            },
            "hpa_subcellular": {
                "label": "Human Protein Atlas subcellular map",
                "pmid": "28495876",
            },
            "hpa_secretome": {
                "label": "Human Protein Atlas human secretome",
                "pmid": "31776213",
            },
            "surfaceome": {
                "label": "The in silico human surfaceome",
                "pmid": "30373828",
                "doi": "10.1073/pnas.1808790115",
            },
        },
    }

    hpa = read_hpa(args.cache_dir, manifest, force=args.force_download)
    apply_hpa_annotations(hpa, generated, symbol_to_geneids, manifest)

    surfaceome = read_surfaceome(args.cache_dir, manifest, force=args.force_download)
    apply_surfaceome_annotations(surfaceome, generated, symbol_to_geneids, manifest)

    matrisome = read_matrisome(args.cache_dir, manifest, force=args.force_download)
    apply_matrisome_annotations(matrisome, generated, symbol_to_geneids, manifest)

    compare_columns = [
        "ER_GOLGI",
        "MATRISOME",
        "SECRETED",
        "Secreted",
        "SurfaceLabel",
        "CellMembrane",
    ]
    summary, details = make_comparison(legacy, generated, compare_columns)

    generated_path = args.outdir / "combined_annotations_rebuilt.tsv"
    summary_path = args.outdir / "annotation_comparison_summary.tsv"
    details_path = args.outdir / "annotation_comparison_by_gene.tsv"
    manifest_path = args.outdir / "annotation_rebuild_manifest.json"

    generated.to_csv(generated_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    details.to_csv(details_path, sep="\t", index=False)
    manifest.update(
        {
            "outputs": {
                "rebuilt_table": str(generated_path.relative_to(ROOT)),
                "comparison_summary": str(summary_path.relative_to(ROOT)),
                "comparison_by_gene": str(details_path.relative_to(ROOT)),
            },
            "comparison_summary": summary.to_dict(orient="records"),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    print(f"Wrote {generated_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {details_path}")
    print(f"Wrote {manifest_path}")
    print("")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
