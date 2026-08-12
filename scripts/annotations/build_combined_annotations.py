#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import time
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from http.cookiejar import CookieJar
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import HTTPCookieProcessor, Request, build_opener, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
ANNOTATION_DIR = ROOT / "tackle" / "data" / "annotations"
LEGACY_TABLE = ANNOTATION_DIR / "combined_annotations_new.tsv"
SOURCE_CACHE = ANNOTATION_DIR / "source_cache"
GENERATED_DIR = ANNOTATION_DIR / "generated"
BACKUP_DIR = ANNOTATION_DIR / "backups"

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
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

DEFAULT_REFRESH_COLUMNS = (
    "ER_GOLGI",
    "SECRETED",
    "SurfaceLabel",
    "CellMembrane",
    "glycomineN",
    "glycomineO",
    "MitoCarta",
    "MitoCarta_Pathways",
    "IDG",
)

SOURCE_BACKED_COLUMNS = (
    *DEFAULT_REFRESH_COLUMNS,
    "MATRISOME",
)

IDENTIFIER_COLUMNS = ("GeneID", "GeneSymbol")
SUSPECT_SURFACEOME_GENEIDS = {
    # Present repeatedly in the Surfaceome master table across unrelated symbols;
    # treat as a source placeholder/anomaly rather than an Entrez GeneID match.
    "101000000",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def serialize_path(path: Path) -> str:
    """Use repository-relative paths when possible, absolute paths otherwise."""
    path = Path(path).resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


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
        "path": serialize_path(dest),
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


def split_geneids(value) -> list[str]:
    text = norm_symbol(value)
    if not text:
        return []
    out = []
    for token in re.split(r"[;,|]\s*|\s+", text):
        token = token.strip()
        if not token or token in SUSPECT_SURFACEOME_GENEIDS:
            continue
        if token.isdigit():
            out.append(token)
    return list(dict.fromkeys(out))


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
    if column not in generated.columns:
        return 0
    hits = 0
    for symbol in symbols:
        for geneid in symbol_to_geneids.get(norm_symbol(symbol), ()):
            generated.loc[generated["GeneID"] == geneid, column] = value
            hits += 1
    return hits


def add_values_by_geneid(
    generated: pd.DataFrame,
    geneids: Iterable[str],
    column: str,
    value: str,
) -> int:
    if column not in generated.columns:
        return 0
    hits = 0
    generated_geneids = set(generated["GeneID"].astype(str))
    for geneid in geneids:
        geneid = str(geneid).strip()
        if not geneid or geneid not in generated_geneids:
            continue
        generated.loc[generated["GeneID"].astype(str) == geneid, column] = value
        hits += 1
    return hits


def add_values_by_geneid_or_symbol(
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    geneids: Iterable[str],
    symbols: Iterable[str],
    column: str,
    value: str,
) -> int:
    geneid_hits = add_values_by_geneid(generated, geneids, column, value)
    symbol_hits = add_values_by_symbol(generated, symbol_to_geneids, symbols, column, value)
    return geneid_hits + symbol_hits


def _merge_er_golgi_value(current: str, value: str) -> str:
    current = str(current or "").strip()
    value = str(value or "").strip()
    if not current:
        return value
    if not value or current == value or current == "BOTH":
        return current
    if value == "BOTH":
        return value
    if {current, value} == {"ENDORETICULUM", "GOLGI"}:
        return "BOTH"
    return value


def add_er_golgi_by_geneid_or_symbol(
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    geneids: Iterable[str],
    symbols: Iterable[str],
    value: str,
) -> int:
    if "ER_GOLGI" not in generated.columns:
        return 0
    target_geneids = set(str(g).strip() for g in geneids if str(g).strip())
    for symbol in symbols:
        target_geneids.update(symbol_to_geneids.get(norm_symbol(symbol), ()))
    hits = 0
    for geneid in target_geneids:
        mask = generated["GeneID"].astype(str) == geneid
        if not bool(mask.any()):
            continue
        current = str(generated.loc[mask, "ER_GOLGI"].iloc[0])
        generated.loc[mask, "ER_GOLGI"] = _merge_er_golgi_value(current, value)
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
    refresh_columns: set[str] | None = None,
) -> None:
    if surfaceome is None:
        return
    refresh_columns = set(refresh_columns or ())
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
    geneid_col = first_existing_column(surfaceome, ("GeneID", "Gene ID", "Entrez GeneID", "Entrez Gene ID"))
    glycomine_n_col = first_existing_column(surfaceome, ("glycomineN sites", "glycomineN"))
    glycomine_o_col = first_existing_column(surfaceome, ("glycomineO sites", "glycomineO"))
    stats = {
        "symbol_column": symbol_col,
        "geneid_column": geneid_col,
        "label_column": label_col,
        "uniprot_subcellular_column": uniprot_subcellular_col,
        "glycomine_n_column": glycomine_n_col,
        "glycomine_o_column": glycomine_o_col,
        "rows_with_direct_geneids": 0,
        "suspect_geneid_rows_skipped": 0,
        "surface_source_rows": 0,
        "nonsurface_source_rows": 0,
        "cell_membrane_source_rows": 0,
        "er_golgi_source_rows": 0,
        "secreted_source_rows": 0,
        "glycomine_n_source_rows": 0,
        "glycomine_o_source_rows": 0,
    }
    if symbol_col is None:
        stats.update({"status": "error", "error": "no gene symbol column"})
        manifest["transforms"]["surfaceome"] = stats
        return

    for _, row in surfaceome.iterrows():
        symbols = split_symbols(row.get(symbol_col, ""))
        raw_geneids = split_geneids(row.get(geneid_col, "") if geneid_col else "")
        if geneid_col and "101000000" in str(row.get(geneid_col, "")):
            stats["suspect_geneid_rows_skipped"] += 1
        if raw_geneids:
            stats["rows_with_direct_geneids"] += 1
        if not symbols and not raw_geneids:
            continue
        label = str(row.get(label_col, "") if label_col else "").strip().lower()
        if "SurfaceLabel" in refresh_columns and (label == "surface" or label in {"1", "true"}):
            add_values_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                "SurfaceLabel",
                "surface",
            )
            stats["surface_source_rows"] += 1
        elif "SurfaceLabel" in refresh_columns and (label == "non-surface" or label == "nonsurface"):
            add_values_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                "SurfaceLabel",
                "nonsurface",
            )
            stats["nonsurface_source_rows"] += 1

        subcellular = str(row.get(uniprot_subcellular_col, "") if uniprot_subcellular_col else "")
        subcellular_lower = subcellular.lower()
        if "CellMembrane" in refresh_columns and "cell membrane" in subcellular.lower():
            add_values_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                "CellMembrane",
                "Cell Membrane",
            )
            stats["cell_membrane_source_rows"] += 1
        has_er = "endoplasmic reticulum" in subcellular_lower or re.search(r"\ber\b", subcellular_lower)
        has_golgi = "golgi" in subcellular_lower
        if "ER_GOLGI" in refresh_columns and (has_er or has_golgi):
            if has_er and has_golgi:
                value = "BOTH"
            elif has_er:
                value = "ENDORETICULUM"
            else:
                value = "GOLGI"
            add_er_golgi_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                value,
            )
            stats["er_golgi_source_rows"] += 1
        if "SECRETED" in refresh_columns and (
            "secreted" in subcellular_lower or "extracellular space" in subcellular_lower
        ):
            add_values_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                "SECRETED",
                "SECRETED",
            )
            stats["secreted_source_rows"] += 1
        glycomine_n = str(row.get(glycomine_n_col, "") if glycomine_n_col else "").strip()
        if (
            "glycomineN" in refresh_columns
            and glycomine_n
            and glycomine_n.lower() not in {"0", "false", "na", "nan", "none"}
        ):
            add_values_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                "glycomineN",
                "N",
            )
            stats["glycomine_n_source_rows"] += 1
        glycomine_o = str(row.get(glycomine_o_col, "") if glycomine_o_col else "").strip()
        if (
            "glycomineO" in refresh_columns
            and glycomine_o
            and glycomine_o.lower() not in {"0", "false", "na", "nan", "none"}
        ):
            add_values_by_geneid_or_symbol(
                generated,
                symbol_to_geneids,
                raw_geneids,
                symbols,
                "glycomineO",
                "O",
            )
            stats["glycomine_o_source_rows"] += 1

    manifest["transforms"]["surfaceome"] = stats


def read_mitocarta(cache: Path, manifest: dict) -> dict[str, pd.DataFrame]:
    files = {
        "human": ("Human.MitoCarta3.0.xls", "A Human MitoCarta3.0", "HumanGeneID"),
        "mouse": ("Mouse.MitoCarta3.0.xls", "A Mouse MitoCarta3.0", "MouseGeneID"),
    }
    out: dict[str, pd.DataFrame] = {}
    source_info: dict[str, dict] = {}
    for species, (filename, sheet_name, geneid_col) in files.items():
        path = cache / filename
        if not path.exists():
            source_info[species] = {
                "path": serialize_path(path),
                "status": "missing",
            }
            continue
        try:
            df = pd.read_excel(path, sheet_name=sheet_name, dtype=str).fillna("")
            out[species] = df
            source_info[species] = {
                "path": serialize_path(path),
                "status": "cached",
                "sheet": sheet_name,
                "geneid_column": geneid_col,
                "rows": int(len(df)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        except Exception as exc:
            source_info[species] = {
                "path": serialize_path(path),
                "status": "error",
                "error": str(exc),
            }
    manifest["sources"]["mitocarta"] = source_info
    return out


def apply_mitocarta_annotations(
    mitocarta: dict[str, pd.DataFrame],
    generated: pd.DataFrame,
    manifest: dict,
) -> None:
    specs = {
        "human": "HumanGeneID",
        "mouse": "MouseGeneID",
    }
    stats: dict[str, dict] = {}
    generated_geneids = set(generated["GeneID"].astype(str))
    for species, geneid_col in specs.items():
        df = mitocarta.get(species)
        species_stats = {
            "geneid_column": geneid_col,
            "source_rows": 0,
            "matched_geneids": 0,
            "mitocarta_values_set": 0,
            "pathway_values_set": 0,
        }
        if df is None or df.empty or geneid_col not in df.columns:
            species_stats["status"] = "missing"
            stats[species] = species_stats
            continue

        for _, row in df.iterrows():
            geneid = str(row.get(geneid_col, "")).strip()
            if not geneid:
                continue
            species_stats["source_rows"] += 1
            if geneid not in generated_geneids:
                continue
            species_stats["matched_geneids"] += 1
            subloc = str(row.get("MitoCarta3.0_SubMitoLocalization", "")).strip() or "unknown"
            pathways = str(row.get("MitoCarta3.0_MitoPathways", "")).strip() or "0"
            if "MitoCarta" in generated.columns:
                generated.loc[generated["GeneID"] == geneid, "MitoCarta"] = subloc
                species_stats["mitocarta_values_set"] += 1
            if "MitoCarta_Pathways" in generated.columns:
                generated.loc[generated["GeneID"] == geneid, "MitoCarta_Pathways"] = pathways
                species_stats["pathway_values_set"] += 1
        stats[species] = species_stats
    manifest["transforms"]["mitocarta"] = stats


def read_idg(cache: Path, manifest: dict) -> pd.DataFrame | None:
    zip_path = cache / "pharos data download.zip"
    csv_path = cache / "IDG_TargetList_Y4.csv"
    source_info: dict[str, dict] = {}

    if zip_path.exists():
        try:
            with zipfile.ZipFile(zip_path) as zf:
                with zf.open("query results.csv") as handle:
                    df = pd.read_csv(handle, dtype=str).fillna("")
                metadata = ""
                if "query metadata.txt" in zf.namelist():
                    with zf.open("query metadata.txt") as handle:
                        metadata = handle.read().decode("utf-8", errors="replace")
            source_info["pharos_zip"] = {
                "path": serialize_path(zip_path),
                "status": "cached",
                "rows": int(len(df)),
                "columns": list(df.columns),
                "sha256": sha256_file(zip_path),
                "bytes": zip_path.stat().st_size,
                "metadata_excerpt": metadata[:2000],
            }
            manifest["sources"]["idg"] = source_info
            return df
        except Exception as exc:
            source_info["pharos_zip"] = {
                "path": serialize_path(zip_path),
                "status": "error",
                "error": str(exc),
            }

    if csv_path.exists():
        try:
            df = pd.read_csv(csv_path, dtype=str).fillna("")
            source_info["target_list"] = {
                "path": serialize_path(csv_path),
                "status": "cached",
                "rows": int(len(df)),
                "columns": list(df.columns),
                "sha256": sha256_file(csv_path),
                "bytes": csv_path.stat().st_size,
            }
            manifest["sources"]["idg"] = source_info
            return df
        except Exception as exc:
            source_info["target_list"] = {
                "path": serialize_path(csv_path),
                "status": "error",
                "error": str(exc),
            }

    if not source_info:
        source_info["status"] = "missing"
    manifest["sources"]["idg"] = source_info
    return None


def apply_idg_annotations(
    idg: pd.DataFrame | None,
    generated: pd.DataFrame,
    symbol_to_geneids: dict[str, set[str]],
    manifest: dict,
) -> None:
    stats = {
        "symbol_column": None,
        "tdl_column": None,
        "source_rows": 0,
        "values_set": 0,
    }
    if idg is None or idg.empty:
        stats["status"] = "missing"
        manifest["transforms"]["idg"] = stats
        return
    symbol_col = first_existing_column(idg, ("Symbol", "Gene", "Gene Symbol", "sym"))
    tdl_col = first_existing_column(idg, ("Target Development Level", "TDL", "IDG"))
    stats["symbol_column"] = symbol_col
    stats["tdl_column"] = tdl_col
    if symbol_col is None or tdl_col is None:
        stats.update({"status": "error", "error": "missing Symbol or Target Development Level column"})
        manifest["transforms"]["idg"] = stats
        return
    for _, row in idg.iterrows():
        symbols = split_symbols(row.get(symbol_col, ""))
        value = str(row.get(tdl_col, "")).strip()
        if not symbols or not value:
            continue
        stats["source_rows"] += 1
        stats["values_set"] += add_values_by_symbol(
            generated,
            symbol_to_geneids,
            symbols,
            "IDG",
            value,
        )
    manifest["transforms"]["idg"] = stats


def _uniprot_cache_path(cache: Path, symbol: str, organism_id: str) -> Path:
    safe_symbol = re.sub(r"[^A-Za-z0-9_.-]+", "_", symbol).strip("_") or "unknown"
    safe_organism = re.sub(r"[^0-9]+", "_", organism_id).strip("_") or "unknown"
    return cache / "uniprot" / f"{safe_symbol}__{safe_organism}.json"


def fetch_uniprot_symbol(
    symbol: str,
    organism_id: str,
    cache: Path,
    *,
    force: bool = False,
    sleep_seconds: float = 0.0,
) -> dict:
    path = _uniprot_cache_path(cache, symbol, organism_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            pass

    query = f"(gene_exact:{symbol}) AND (organism_id:{organism_id})"
    url = UNIPROT_SEARCH_URL + "?" + urlencode({"query": query, "format": "json", "size": "5"})
    payload: dict
    try:
        raw = request_bytes(url)
        payload = json.loads(raw.decode("utf-8"))
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
        payload = {"results": [], "error": str(exc), "query": query, "url": url}
    payload.setdefault("query", query)
    payload.setdefault("url", url)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if sleep_seconds:
        time.sleep(sleep_seconds)
    return payload


def _join_unique_text(values: Iterable[str]) -> str:
    return "; ".join(dict.fromkeys(str(t).strip() for t in values if str(t).strip()))


def _uniprot_entry_evidence_fields(entry: dict) -> dict[str, str]:
    locations: list[str] = []
    location_notes: list[str] = []
    keywords: list[str] = []
    go_cellular_component: list[str] = []
    go_molecular_function: list[str] = []
    go_biological_process: list[str] = []
    for comment in entry.get("comments", []) or []:
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for loc in comment.get("subcellularLocations", []) or []:
            for key in ("location", "topology", "orientation"):
                value = loc.get(key)
                if isinstance(value, dict) and value.get("value"):
                    locations.append(str(value.get("value")))
        for text in comment.get("texts", []) or []:
            if isinstance(text, dict) and text.get("value"):
                location_notes.append(str(text.get("value")))

    for keyword in entry.get("keywords", []) or []:
        value = keyword.get("name") if isinstance(keyword, dict) else keyword
        if value:
            keywords.append(str(value))

    for xref in entry.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "GO":
            continue
        for prop in xref.get("properties", []) or []:
            if prop.get("key") == "GoTerm" and prop.get("value"):
                value = str(prop.get("value"))
                if value.startswith("C:"):
                    go_cellular_component.append(value)
                elif value.startswith("F:"):
                    go_molecular_function.append(value)
                elif value.startswith("P:"):
                    go_biological_process.append(value)

    return {
        "uniprot_locations": _join_unique_text(locations),
        "uniprot_location_notes": _join_unique_text(location_notes),
        "uniprot_keywords": _join_unique_text(keywords),
        "go_cellular_component": _join_unique_text(go_cellular_component),
        "go_molecular_function": _join_unique_text(go_molecular_function),
        "go_biological_process": _join_unique_text(go_biological_process),
    }


def _classify_uniprot_subcellular(text: str) -> dict[str, str]:
    lower = text.lower()
    has_er = "endoplasmic reticulum" in lower or re.search(r"\ber\b", lower)
    has_golgi = "golgi" in lower
    if has_er and has_golgi:
        er_golgi = "BOTH"
    elif has_er:
        er_golgi = "ENDORETICULUM"
    elif has_golgi:
        er_golgi = "GOLGI"
    else:
        er_golgi = ""

    secreted = "SECRETED" if (
        "secreted" in lower
        or "extracellular space" in lower
        or "extracellular region" in lower
    ) else ""

    ecm = "ECM_EVIDENCE" if (
        "extracellular matrix" in lower
        or "basement membrane" in lower
        or "collagen" in lower
        or "matrisome" in lower
    ) else ""

    return {"ER_GOLGI": er_golgi, "SECRETED": secreted, "MATRISOME_EVIDENCE": ecm}


def apply_uniprot_backfill(
    generated: pd.DataFrame,
    manifest: dict,
    evidence_path: Path,
    *,
    refresh_columns: set[str],
    cache: Path,
    organism_ids: tuple[str, ...],
    limit: int,
    force: bool,
    sleep_seconds: float,
    matrisome_evidence: bool,
) -> None:
    target_columns = {"ER_GOLGI", "SECRETED"} & refresh_columns
    collect_matrisome_evidence = bool(matrisome_evidence and "MATRISOME" in generated.columns)
    if not target_columns and not collect_matrisome_evidence:
        return

    candidate = generated.copy()
    masks = []
    if "ER_GOLGI" in target_columns and "ER_GOLGI" in candidate.columns:
        masks.append(~positive_mask(candidate["ER_GOLGI"]))
    if "SECRETED" in target_columns and "SECRETED" in candidate.columns:
        masks.append(~positive_mask(candidate["SECRETED"]))
    if collect_matrisome_evidence and "MATRISOME" in candidate.columns:
        masks.append(~positive_mask(candidate["MATRISOME"]))
    if not masks:
        return

    missing_mask = masks[0]
    for mask in masks[1:]:
        missing_mask = missing_mask | mask
    rows = candidate.loc[missing_mask, ["GeneID", "GeneSymbol"]].drop_duplicates()
    if limit > 0:
        rows = rows.head(limit)

    evidence_rows: list[dict] = []
    stats = {
        "candidate_rows": int(len(rows)),
        "queried_symbols": 0,
        "evidence_rows": 0,
        "secreted_values_set": 0,
        "er_golgi_values_set": 0,
        "matrisome_evidence_rows": 0,
        "limit": int(limit),
        "organism_ids": list(organism_ids),
    }

    for _, row in rows.iterrows():
        geneid = str(row.get("GeneID", "")).strip()
        symbol = str(row.get("GeneSymbol", "")).strip()
        if not symbol:
            continue
        stats["queried_symbols"] += 1
        for organism_id in organism_ids:
            payload = fetch_uniprot_symbol(
                symbol,
                organism_id,
                cache,
                force=force,
                sleep_seconds=sleep_seconds,
            )
            for entry in payload.get("results", []) or []:
                evidence_fields = _uniprot_entry_evidence_fields(entry)
                location_text = evidence_fields["uniprot_locations"]
                if not location_text:
                    continue
                calls = _classify_uniprot_subcellular(location_text)
                if not any(calls.values()):
                    continue
                accession = entry.get("primaryAccession", "")
                entry_genes = []
                for gene in entry.get("genes", []) or []:
                    gene_name = gene.get("geneName", {}) if isinstance(gene, dict) else {}
                    if gene_name.get("value"):
                        entry_genes.append(str(gene_name.get("value")))
                if "SECRETED" in target_columns and calls["SECRETED"]:
                    hits = add_values_by_geneid(generated, [geneid], "SECRETED", calls["SECRETED"])
                    stats["secreted_values_set"] += hits
                if "ER_GOLGI" in target_columns and calls["ER_GOLGI"]:
                    current = str(generated.loc[generated["GeneID"].astype(str) == geneid, "ER_GOLGI"].iloc[0]) if "ER_GOLGI" in generated.columns and bool((generated["GeneID"].astype(str) == geneid).any()) else ""
                    merged = _merge_er_golgi_value(current, calls["ER_GOLGI"])
                    hits = add_values_by_geneid(generated, [geneid], "ER_GOLGI", merged)
                    stats["er_golgi_values_set"] += hits
                if calls["MATRISOME_EVIDENCE"]:
                    stats["matrisome_evidence_rows"] += 1
                evidence_rows.append(
                    {
                        "GeneID": geneid,
                        "GeneSymbol": symbol,
                        "organism_id": organism_id,
                        "accession": accession,
                        "entry_gene_names": "|".join(entry_genes),
                        "ER_GOLGI": calls["ER_GOLGI"],
                        "SECRETED": calls["SECRETED"],
                        "MATRISOME_EVIDENCE": calls["MATRISOME_EVIDENCE"],
                        **evidence_fields,
                    }
                )

    evidence = pd.DataFrame(evidence_rows)
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(evidence_path, sep="\t", index=False)
    stats["evidence_rows"] = int(len(evidence))
    manifest["transforms"]["uniprot_backfill"] = stats


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
                "path": serialize_path(path),
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
        value_changed = overlap & (merged[lcol].astype(str) != merged[gcol].astype(str))
        unchanged = overlap & ~value_changed
        denom = int((legacy_pos | generated_pos).sum())
        summary_rows.append(
            {
                "column": column,
                "legacy_positive": int(legacy_pos.sum()),
                "generated_positive": int(generated_pos.sum()),
                "overlap": int(overlap.sum()),
                "unchanged": int(unchanged.sum()),
                "value_changed": int(value_changed.sum()),
                "legacy_only": int(legacy_only.sum()),
                "generated_only": int(generated_only.sum()),
                "jaccard": round(float(overlap.sum() / denom), 4) if denom else 1.0,
            }
        )
        details = merged.loc[
            legacy_only | generated_only | value_changed,
            ["GeneID", "GeneSymbol", lcol, gcol],
        ].copy()
        details["change_type"] = ""
        details.loc[legacy_only, "change_type"] = "legacy_only"
        details.loc[generated_only, "change_type"] = "generated_only"
        details.loc[value_changed, "change_type"] = "value_changed"
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
    parser.add_argument("--backup-dir", type=Path, default=BACKUP_DIR)
    parser.add_argument(
        "--refresh-column",
        action="append",
        choices=SOURCE_BACKED_COLUMNS,
        dest="refresh_columns",
        help=(
            "Column to refresh from source data. Can be repeated. "
            "Defaults refresh source-backed candidate columns except MATRISOME."
        ),
    )
    parser.add_argument(
        "--include-matrisome",
        action="store_true",
        help="Also refresh MATRISOME from MatrisomeDB. Disabled by default because current source parity needs review.",
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        help="Back up the live combined_annotations_new.tsv and replace it with the candidate table if different.",
    )
    parser.add_argument(
        "--uniprot-backfill",
        action="store_true",
        help="Query cached/remote UniProt records for rows still missing selected SECRETED or ER_GOLGI evidence.",
    )
    parser.add_argument(
        "--uniprot-organism",
        action="append",
        default=None,
        help="UniProt organism_id to query for backfill; can be repeated. Defaults to 9606 and 10090.",
    )
    parser.add_argument(
        "--uniprot-backfill-limit",
        type=int,
        default=500,
        help="Maximum missing annotation rows to query when --uniprot-backfill is enabled; 0 means no limit.",
    )
    parser.add_argument(
        "--uniprot-sleep",
        type=float,
        default=0.0,
        help="Seconds to sleep after each uncached UniProt request.",
    )
    parser.add_argument(
        "--uniprot-matrisome-evidence",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="When --uniprot-backfill is enabled, record ECM-like UniProt evidence for Matrisome review without updating MATRISOME.",
    )
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    refresh_columns = list(args.refresh_columns or DEFAULT_REFRESH_COLUMNS)
    if args.include_matrisome and "MATRISOME" not in refresh_columns:
        refresh_columns.append("MATRISOME")
    refresh_columns = [c for c in dict.fromkeys(refresh_columns) if c in SOURCE_BACKED_COLUMNS]

    legacy = load_legacy()
    generated = legacy.copy()
    for column in refresh_columns:
        if column in generated.columns:
            generated[column] = ""

    symbol_to_geneids: dict[str, set[str]] = defaultdict(set)
    for _, row in legacy.iterrows():
        symbol = norm_symbol(row.get("GeneSymbol"))
        geneid = norm_symbol(row.get("GeneID"))
        if symbol and geneid:
            symbol_to_geneids[symbol].add(geneid)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "legacy_table": serialize_path(LEGACY_TABLE),
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
                "doi": "10.1126/science.aal3321",
            },
            "hpa_secretome": {
                "label": "Human Protein Atlas human secretome",
                "pmid": "31772123",
                "doi": "10.1126/scisignal.aaz0274",
            },
            "surfaceome": {
                "label": "The in silico human surfaceome",
                "pmid": "30373828",
                "doi": "10.1073/pnas.1808790115",
            },
            "mitocarta3": {
                "label": "MitoCarta3.0 mitochondrial protein inventory",
                "pmid": "33174596",
                "doi": "10.1093/nar/gkaa1011",
            },
            "pharos": {
                "label": "Pharos target-development-level annotations",
                "doi": "10.1093/nar/gkac1033",
            },
        },
        "settings": {
            "refresh_columns": refresh_columns,
            "promote": bool(args.promote),
            "include_matrisome": bool(args.include_matrisome),
        },
    }

    if {"ER_GOLGI", "SECRETED"} & set(refresh_columns):
        hpa = read_hpa(args.cache_dir, manifest, force=args.force_download)
        apply_hpa_annotations(hpa, generated, symbol_to_geneids, manifest)

    if {
        "SurfaceLabel",
        "CellMembrane",
        "glycomineN",
        "glycomineO",
        "ER_GOLGI",
        "SECRETED",
    } & set(refresh_columns):
        surfaceome = read_surfaceome(args.cache_dir, manifest, force=args.force_download)
        apply_surfaceome_annotations(
            surfaceome,
            generated,
            symbol_to_geneids,
            manifest,
            refresh_columns=set(refresh_columns),
        )

    if {"MitoCarta", "MitoCarta_Pathways"} & set(refresh_columns):
        mitocarta = read_mitocarta(args.cache_dir, manifest)
        apply_mitocarta_annotations(mitocarta, generated, manifest)

    if "IDG" in refresh_columns:
        idg = read_idg(args.cache_dir, manifest)
        apply_idg_annotations(idg, generated, symbol_to_geneids, manifest)

    if "MATRISOME" in refresh_columns:
        matrisome = read_matrisome(args.cache_dir, manifest, force=args.force_download)
        apply_matrisome_annotations(matrisome, generated, symbol_to_geneids, manifest)

    uniprot_evidence_path = args.outdir / "uniprot_subcellular_evidence.tsv"
    if args.uniprot_backfill:
        apply_uniprot_backfill(
            generated,
            manifest,
            uniprot_evidence_path,
            refresh_columns=set(refresh_columns),
            cache=args.cache_dir,
            organism_ids=tuple(args.uniprot_organism or ("9606", "10090")),
            limit=max(0, args.uniprot_backfill_limit),
            force=args.force_download,
            sleep_seconds=max(0.0, args.uniprot_sleep),
            matrisome_evidence=args.uniprot_matrisome_evidence,
        )

    compare_columns = [
        c for c in generated.columns if c not in IDENTIFIER_COLUMNS
    ]
    summary, details = make_comparison(legacy, generated, compare_columns)

    generated_path = args.outdir / "combined_annotations_rebuilt.tsv"
    summary_path = args.outdir / "annotation_comparison_summary.tsv"
    details_path = args.outdir / "annotation_comparison_by_gene.tsv"
    manifest_path = args.outdir / "annotation_rebuild_manifest.json"

    identical_to_live = legacy.fillna("").astype(str).equals(generated.fillna("").astype(str))

    generated.to_csv(generated_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    details.to_csv(details_path, sep="\t", index=False)
    manifest.update(
        {
            "outputs": {
            "rebuilt_table": serialize_path(generated_path),
            "comparison_summary": serialize_path(summary_path),
            "comparison_by_gene": serialize_path(details_path),
            "uniprot_subcellular_evidence": serialize_path(uniprot_evidence_path)
            if args.uniprot_backfill and uniprot_evidence_path.exists()
            else "",
        },
        "comparison_summary": summary.to_dict(orient="records"),
        "candidate_identical_to_live": bool(identical_to_live),
    }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    promoted = False
    backup_path = None
    if args.promote:
        if identical_to_live:
            print("Candidate is identical to live annotation table; not promoting.")
        else:
            args.backup_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            backup_path = args.backup_dir / f"{LEGACY_TABLE.stem}.{stamp}.tsv"
            shutil.copy2(LEGACY_TABLE, backup_path)
            generated.to_csv(LEGACY_TABLE, sep="\t", index=False)
            promoted = True

    print(f"Wrote {generated_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {details_path}")
    print(f"Wrote {manifest_path}")
    if promoted:
        print(f"Promoted candidate to {LEGACY_TABLE}")
        print(f"Backed up prior live table to {backup_path}")
    elif not args.promote:
        print("Did not promote candidate; rerun with --promote to replace the live table after review.")
    print("")
    print(f"Refresh columns: {', '.join(refresh_columns)}")
    print(f"Candidate identical to live table: {identical_to_live}")
    print("")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
