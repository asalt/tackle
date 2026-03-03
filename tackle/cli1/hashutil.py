from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    file_path = Path(path)
    hasher = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def tsv_dims(path: str | Path) -> tuple[int | None, int | None]:
    """Return (rows, cols) for a TSV without loading it into pandas.

    Assumes first row is header, tab-delimited.
    """
    file_path = Path(path)
    try:
        with file_path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            if not header:
                return 0, 0
            cols = len(header.rstrip("\n").split("\t"))
            rows = 0
            for _ in handle:
                rows += 1
        return rows, cols
    except OSError:
        return None, None


def gct_dims(path: str | Path) -> tuple[int | None, int | None]:
    """Return (rows, cols) for GCT v1.2/v1.3.

    GCT files start with a version line, then a line containing: <nrows><tab><ncols>.
    """
    file_path = Path(path)
    try:
        with file_path.open("r", encoding="utf-8") as handle:
            _version = handle.readline()
            dims = handle.readline()
            if not dims:
                return None, None
            parts = dims.strip().split("\t")
            if len(parts) < 2:
                return None, None
            return int(parts[0]), int(parts[1])
    except (OSError, ValueError):
        return None, None

