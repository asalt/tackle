"""Histone exact-duplicate reduction helpers.

The caller passes the long ``Data.data`` frame produced by ``containers.Data``.
That frame has one row per ``GeneID``/``Metric`` pair:

    GeneID | Metric          | sample_a | sample_b | ...
    8336   | TaxonID         | 9606     | 9606     | ...
    8336   | AreaSum_dstrAdj | 123.4    | 0        | ...
    8336   | area            | 12.3     | 0        | ...

Only exact duplicate histone ``AreaSum_dstrAdj`` vectors are reduced, and only
within the same TaxonID. This preserves the existing "shared peptide evidence"
semantics: identical histone rows are represented once in downstream analyses,
but the manifest records all plausible histone identifiers represented by that
selected row.
"""

from __future__ import annotations

import os
from typing import Iterable

import pandas as pd


MANIFEST_COLUMNS = [
    "taxonid",
    "selected_geneid",
    "selected_genesymbol",
    "candidate_geneids",
    "candidate_genesymbols",
]


def _value_columns(data: pd.DataFrame) -> list:
    return [x for x in data.columns if x not in ("GeneID", "Metric")]


def _first_nonempty(row, columns: Iterable) -> str:
    for col in columns:
        value = row[col]
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text and text.lower() != "nan":
            return text
    return ""


def _histone_geneids(histone_info) -> set:
    if not hasattr(histone_info, "df") or "GeneID" not in histone_info.df.columns:
        return set()
    return {
        str(x).strip()
        for x in histone_info.df["GeneID"].dropna().tolist()
        if str(x).strip()
    }


def _histone_symbol_map(histone_info) -> dict:
    if not hasattr(histone_info, "df"):
        return {}
    if not {"GeneID", "GeneSymbol"}.issubset(histone_info.df.columns):
        return {}

    symbol_map = {}
    for _, row in histone_info.df[["GeneID", "GeneSymbol"]].iterrows():
        geneid = str(row["GeneID"]).strip()
        symbol = str(row["GeneSymbol"]).strip()
        if not geneid or geneid.lower() == "nan":
            continue
        if not symbol or symbol.lower() == "nan":
            continue
        symbol_map.setdefault(geneid, [])
        if symbol not in symbol_map[geneid]:
            symbol_map[geneid].append(symbol)
    return symbol_map


def _symbols_for(geneids: Iterable[str], symbol_map: dict) -> str:
    symbols = []
    for geneid in geneids:
        for symbol in symbol_map.get(str(geneid), []):
            if symbol not in symbols:
                symbols.append(symbol)
    return ",".join(symbols)


def _dedupe_key(row, value_cols: list) -> tuple:
    key = []
    for col in value_cols:
        value = row[col]
        if pd.isna(value):
            key.append("<NA>")
        else:
            key.append(value)
    return tuple(key)


def _geneid_to_taxon(data: pd.DataFrame, value_cols: list) -> dict:
    if "TaxonID" not in set(data["Metric"]):
        return {}
    taxon_rows = data.loc[data["Metric"] == "TaxonID", ["GeneID", *value_cols]].copy()
    if taxon_rows.empty:
        return {}
    taxon_rows["_GeneID_str"] = taxon_rows["GeneID"].astype(str)
    taxon_rows["_taxonid"] = taxon_rows.apply(
        lambda row: _first_nonempty(row, value_cols), axis=1
    )
    return dict(zip(taxon_rows["_GeneID_str"], taxon_rows["_taxonid"]))


def _sortable_geneid(geneid):
    text = str(geneid)
    return int(text) if text.isnumeric() else text


def write_histone_dedupe_manifest(records: list, outpath: str | None, logger=None) -> None:
    if not outpath:
        return

    manifest_path = os.path.join(outpath, "context", "histone_dedupe.tsv")
    try:
        os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
        pd.DataFrame(records, columns=MANIFEST_COLUMNS).to_csv(
            manifest_path, sep="\t", index=False
        )
        if records and logger is not None:
            logger.info("Wrote histone dedupe manifest: %s", manifest_path)
    except Exception as exc:
        if logger is not None:
            logger.warning("Failed to write histone dedupe manifest: %s", exc)


def dedupe_histone_rows(
    data: pd.DataFrame,
    histone_info_object,
    outpath: str | None = None,
    logger=None,
) -> pd.DataFrame:
    """Remove exact duplicate histone rows and write a transparent manifest.

    Parameters are intentionally narrow:
    ``data`` is the long ``self.data`` frame, ``histone_info_object`` yields
    histone lookup tables, and ``outpath`` determines where the optional
    ``context/histone_dedupe.tsv`` sidecar is written.
    """

    value_cols = _value_columns(data)
    if not value_cols:
        write_histone_dedupe_manifest([], outpath, logger=logger)
        return data

    deduped = data
    records = []
    taxon_by_geneid = _geneid_to_taxon(deduped, value_cols)

    for histone_info in histone_info_object:
        histone_geneids = _histone_geneids(histone_info)
        if not histone_geneids:
            continue

        symbol_map = _histone_symbol_map(histone_info)
        histone_values = deduped.loc[
            deduped["GeneID"].astype(str).isin(histone_geneids)
            & (deduped["Metric"] == "AreaSum_dstrAdj")
        ].copy()
        if histone_values.empty:
            continue

        histone_values["_GeneID_str"] = histone_values["GeneID"].astype(str)
        histone_values["_taxonid"] = (
            histone_values["_GeneID_str"].map(taxon_by_geneid).fillna("<unknown>")
        )
        histone_values["_geneid_sortable"] = histone_values["GeneID"].apply(
            _sortable_geneid
        )
        histone_values["_dedupe_key"] = histone_values.apply(
            lambda row: _dedupe_key(row, value_cols), axis=1
        )
        histone_values = histone_values.sort_values(
            by=["_taxonid", value_cols[0], "_geneid_sortable"]
        )

        for _, group in histone_values.groupby(["_taxonid", "_dedupe_key"], sort=False):
            if len(group) <= 1:
                continue
            candidate_geneids = [str(x) for x in group["GeneID"].tolist()]
            selected_geneid = candidate_geneids[0]
            records.append(
                {
                    "taxonid": str(group["_taxonid"].iloc[0]),
                    "selected_geneid": selected_geneid,
                    "selected_genesymbol": _symbols_for([selected_geneid], symbol_map),
                    "candidate_geneids": ",".join(candidate_geneids),
                    "candidate_genesymbols": _symbols_for(candidate_geneids, symbol_map),
                }
            )

        keep_rows = histone_values.drop_duplicates(["_taxonid", *value_cols])
        remove_geneids = set(histone_values["_GeneID_str"]) - set(
            keep_rows["_GeneID_str"]
        )
        if remove_geneids:
            deduped = deduped.loc[~deduped["GeneID"].astype(str).isin(remove_geneids)]

    write_histone_dedupe_manifest(records, outpath, logger=logger)
    return deduped
