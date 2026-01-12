from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import pandas as pd


def summarize_silhouette_df(df: pd.DataFrame) -> dict[str, Any]:
    if df.empty:
        return {
            "n_genes": 0,
            "n_clusters": 0,
            "sil_mean": float("nan"),
            "sil_q10": float("nan"),
            "sil_neg_frac": float("nan"),
        }

    clusters = pd.to_numeric(df.get("cluster"), errors="coerce")
    n_clusters = int(pd.Series(clusters).dropna().nunique())

    sil = pd.to_numeric(df.get("sil_width"), errors="coerce")
    sil = sil[np.isfinite(sil)]
    if len(sil) == 0:
        sil_mean = sil_q10 = sil_neg_frac = float("nan")
    else:
        sil_mean = float(sil.mean())
        sil_q10 = float(sil.quantile(0.10))
        sil_neg_frac = float((sil < 0).mean())

    return {
        "n_genes": int(len(df)),
        "n_clusters": n_clusters,
        "sil_mean": sil_mean,
        "sil_q10": sil_q10,
        "sil_neg_frac": sil_neg_frac,
    }


def select_best_k(candidates: Iterable[dict[str, Any]]) -> int:
    rows = list(candidates)
    if not rows:
        raise ValueError("No k candidates provided")

    def finite_or_neg_inf(value: Any) -> float:
        try:
            val = float(value)
        except Exception:
            return float("-inf")
        return val if np.isfinite(val) else float("-inf")

    def finite_or_pos_inf(value: Any) -> float:
        try:
            val = float(value)
        except Exception:
            return float("inf")
        return val if np.isfinite(val) else float("inf")

    def sort_key(row: dict[str, Any]) -> tuple[float, float, float, int]:
        sil_mean = finite_or_neg_inf(row.get("sil_mean"))
        sil_q10 = finite_or_neg_inf(row.get("sil_q10"))
        sil_neg_frac = finite_or_pos_inf(row.get("sil_neg_frac"))
        k = int(row["nclusters"])
        return (-sil_mean, -sil_q10, sil_neg_frac, k)

    best = min(rows, key=sort_key)
    return int(best["nclusters"])

