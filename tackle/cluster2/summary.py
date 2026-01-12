from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


_CLUSTER_COLUMN_ALIASES = {
    "Cluster": "cluster",
    "Kmeans_Cluster": "cluster",
    "kmeans_cluster": "cluster",
    "kmeans_Cluster": "cluster",
    "Kmeans_cluster": "cluster",
    "PAM_cluster": "cluster",
}


@dataclass(frozen=True)
class ClusterSummaryTables:
    runs: pd.DataFrame
    clusters: pd.DataFrame


def _read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path)
    return pd.read_table(path)


def load_cluster_metrics(path: str) -> pd.DataFrame:
    df = _read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    df = df.rename(columns=_CLUSTER_COLUMN_ALIASES)
    if "cluster" not in df.columns:
        raise ValueError(f"{path} is missing a `cluster` column")
    if "sil_width" not in df.columns:
        raise ValueError(f"{path} is missing a `sil_width` column")
    df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").astype("Int64")
    df["sil_width"] = pd.to_numeric(df["sil_width"], errors="coerce")
    return df


def summarize_cluster_metrics(df: pd.DataFrame, *, run_id: str) -> ClusterSummaryTables:
    clean = df.copy()
    clean["run_id"] = run_id

    def _finite_sil(series: pd.Series) -> pd.Series:
        vals = pd.to_numeric(series, errors="coerce")
        return vals[np.isfinite(vals)]

    def _sil_q10(series: pd.Series) -> float:
        vals = _finite_sil(series)
        return float(vals.quantile(0.10)) if len(vals) else float("nan")

    def _sil_neg_frac(series: pd.Series) -> float:
        vals = _finite_sil(series)
        return float((vals < 0).mean()) if len(vals) else float("nan")

    cluster_summary = (
        clean.groupby(["run_id", "cluster"], dropna=False)
        .agg(
            n_genes=("cluster", "size"),
            sil_mean=("sil_width", "mean"),
            sil_median=("sil_width", "median"),
            sil_q10=("sil_width", _sil_q10),
            sil_neg_frac=("sil_width", _sil_neg_frac),
        )
        .reset_index()
        .sort_values(by=["run_id", "cluster"], kind="stable")
    )

    run_summary = (
        clean.groupby("run_id")
        .agg(
            n_genes=("cluster", "size"),
            n_clusters=("cluster", lambda s: int(pd.Series(s).dropna().nunique())),
            sil_mean=("sil_width", "mean"),
            sil_median=("sil_width", "median"),
            sil_q10=("sil_width", _sil_q10),
            sil_neg_frac=("sil_width", _sil_neg_frac),
        )
        .reset_index()
    )

    return ClusterSummaryTables(runs=run_summary, clusters=cluster_summary)


def summarize_cluster_files(
    paths: Iterable[str], *, run_id_prefix: Optional[str] = None
) -> ClusterSummaryTables:
    path_list = [str(p) for p in paths]
    run_tables: list[pd.DataFrame] = []
    cluster_tables: list[pd.DataFrame] = []
    for path in path_list:
        run_id = Path(path).name
        if run_id_prefix:
            run_id = f"{run_id_prefix}{run_id}"
        df = load_cluster_metrics(path)
        tables = summarize_cluster_metrics(df, run_id=run_id)
        run_tables.append(tables.runs.assign(path=path))
        cluster_tables.append(tables.clusters.assign(path=path))

    runs = pd.concat(run_tables, ignore_index=True) if run_tables else pd.DataFrame()
    clusters = (
        pd.concat(cluster_tables, ignore_index=True) if cluster_tables else pd.DataFrame()
    )
    return ClusterSummaryTables(runs=runs, clusters=clusters)

