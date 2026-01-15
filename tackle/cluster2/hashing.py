from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any, Iterable

import pandas as pd


def _stable_json_dumps(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _hash_pandas_object(values: Any) -> bytes:
    return pd.util.hash_pandas_object(values, index=False).values.tobytes()


def hash_cluster2_input_matrix(
    df: pd.DataFrame,
    *,
    geneid_column: str = "GeneID",
    ignore_columns: Iterable[str] = ("GeneID", "GeneSymbol"),
) -> str:
    ignore_set = set(ignore_columns)
    sample_cols = [c for c in df.columns if c not in ignore_set]

    if geneid_column in df.columns:
        gene_ids = df[geneid_column].astype(str)
    else:
        gene_ids = pd.Index(df.index.astype(str))

    mat = df[sample_cols]
    if any(dtype == "object" for dtype in mat.dtypes):
        mat = mat.apply(pd.to_numeric, errors="coerce")

    hasher = hashlib.sha1()
    hasher.update(b"cluster2_input_v1\0")
    hasher.update(str(mat.shape[0]).encode("utf-8"))
    hasher.update(b"\0")
    hasher.update(str(mat.shape[1]).encode("utf-8"))
    hasher.update(b"\0")
    hasher.update(_hash_pandas_object(pd.Index(sample_cols).astype(str)))
    hasher.update(_hash_pandas_object(pd.Index(gene_ids).astype(str)))
    hasher.update(pd.util.hash_pandas_object(mat, index=False).values.tobytes())
    return hasher.hexdigest()


@dataclass(frozen=True)
class Cluster2AutosweepCacheKey:
    key: str
    params: dict[str, Any]

    @classmethod
    def build(cls, *, input_hash: str, **params: Any) -> "Cluster2AutosweepCacheKey":
        payload = {"input_hash": str(input_hash), **params}
        key = hashlib.sha1(_stable_json_dumps(payload).encode("utf-8")).hexdigest()
        return cls(key=key, params=payload)


@dataclass(frozen=True)
class Cluster2RowHclustCacheKey:
    key: str
    params: dict[str, Any]

    @classmethod
    def build(cls, *, input_hash: str, **params: Any) -> "Cluster2RowHclustCacheKey":
        payload = {"input_hash": str(input_hash), **params}
        key = hashlib.sha1(_stable_json_dumps(payload).encode("utf-8")).hexdigest()
        return cls(key=key, params=payload)


def is_autosweep_cache_hit(df_sweep: pd.DataFrame, *, cache_key: str) -> bool:
    if df_sweep is None or not isinstance(df_sweep, pd.DataFrame) or df_sweep.empty:
        return False
    if "sweep_key" not in df_sweep.columns:
        return False
    keys = df_sweep["sweep_key"].astype(str).dropna().unique()
    return len(keys) == 1 and keys[0] == str(cache_key)
