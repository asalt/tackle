from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd


def _json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item"):
        try:
            return obj.item()
        except Exception:
            pass
    return str(obj)


def _string_bytes(values) -> bytes:
    return b"\x1e".join(str(x).encode("utf-8") for x in values)


def compute_imputation_cache_key(
    frame: pd.DataFrame,
    *,
    backend: str,
    gaussian_method: Optional[str] = None,
    lupine_mode: Optional[str] = None,
    params: Optional[Mapping[str, Any]] = None,
) -> str:
    arr = frame.to_numpy(dtype=np.float64, copy=True)
    nan_mask = np.isnan(arr)
    arr_no_nan = np.nan_to_num(arr, nan=0.0)

    meta = {
        "backend": str(backend),
        "gaussian_method": str(gaussian_method) if gaussian_method is not None else None,
        "lupine_mode": str(lupine_mode) if lupine_mode is not None else None,
        "params": dict(params or {}),
    }

    digest = hashlib.sha1()
    digest.update(str(frame.shape).encode("utf-8"))
    digest.update(_string_bytes(frame.index.astype(str)))
    digest.update(b"\x1f")
    digest.update(_string_bytes(frame.columns.astype(str)))
    digest.update(b"\x1f")
    digest.update(nan_mask.tobytes(order="C"))
    digest.update(arr_no_nan.tobytes(order="C"))
    digest.update(
        json.dumps(meta, sort_keys=True, separators=(",", ":"), default=_json_default).encode("utf-8")
    )
    return digest.hexdigest()


def _cache_paths(cache_dir: str | Path, key: str) -> tuple[Path, Path]:
    root = Path(cache_dir).expanduser().resolve()
    return root / f"{key}.pkl.gz", root / f"{key}.json"


def load_imputation_cache(cache_dir: str | Path, key: str) -> Optional[pd.DataFrame]:
    matrix_path, _meta_path = _cache_paths(cache_dir, key)
    if not matrix_path.exists():
        return None
    return pd.read_pickle(matrix_path)


def save_imputation_cache(
    cache_dir: str | Path,
    key: str,
    frame: pd.DataFrame,
    *,
    metadata: Optional[Mapping[str, Any]] = None,
    overwrite: bool = False,
) -> Path:
    matrix_path, meta_path = _cache_paths(cache_dir, key)
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    if matrix_path.exists() and not overwrite:
        return matrix_path

    frame.to_pickle(matrix_path, compression="gzip")

    payload = {
        "cache_key": key,
        "created_at": datetime.now().astimezone().isoformat(),
        "matrix_path": str(matrix_path),
        "shape": [int(frame.shape[0]), int(frame.shape[1])],
        "metadata": dict(metadata or {}),
    }
    meta_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    return matrix_path
