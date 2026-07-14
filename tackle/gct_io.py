from __future__ import annotations

import json
import hashlib
import os
import re
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd


_HASHED_STEM_RE = re.compile(r"\.[0-9a-f]{20,64}$")


def gct_io_contracts() -> dict[str, Any]:
    """Return tackle's compact GCT/GCTX serialization contract registry.

    Keeping these values behind one function avoids populating every importing
    module with a collection of contract constants.  Callers that persist a
    contract value should read it once at the point where the artifact is made.
    """

    return {
        "gctx": {
            "provenance_schema": 1,
            "hash_contract": "tackle-gct-payload-v2",
            "writer_contract": "tackle-gctx-writer-v4",
            "filename_hash_chars": 20,
            "compression_min_bytes": 4096,
            "attrs": {
                "schema": "tackle_provenance_schema",
                "hash": "tackle_content_hash",
                "hash_algorithm": "tackle_content_hash_algorithm",
                "hash_contract": "tackle_hash_contract",
                "writer_contract": "tackle_writer_contract",
                "matrix_dtype": "tackle_matrix_dtype",
                "row_metadata_order": "tackle_row_metadata_columns",
                "column_metadata_order": "tackle_col_metadata_columns",
                "legacy_sha256": "tackle_content_sha256",
            },
        },
        "gct_text_writer": "tackle-gct-text-writer-v1",
        "tsv_hash": "tackle-tsv-v2",
    }


def _gctx_attr(name: str) -> str:
    return str(gct_io_contracts()["gctx"]["attrs"][name])


def _new_content_hasher(algorithm: str = "auto") -> tuple[str, Any]:
    requested = str(algorithm or "auto").strip().lower()
    if requested in {"auto", "blake3"}:
        try:
            import blake3

            return "blake3", blake3.blake3()
        except ImportError:
            if requested == "blake3":
                raise RuntimeError(
                    "BLAKE3 was requested but the 'blake3' Python package is unavailable."
                )
    if requested in {"auto", "blake2b", "blake2b-256"}:
        return "blake2b-256", hashlib.blake2b(digest_size=32)
    if requested in {"sha256", "sha-256"}:
        return "sha256", hashlib.sha256()
    raise ValueError("hash_algorithm must be auto, blake3, blake2b-256, or sha256")


def selected_content_hash_algorithm() -> str:
    """Return the algorithm selected by the current environment's auto policy."""

    algorithm, _ = _new_content_hasher("auto")
    return algorithm


@dataclass(frozen=True)
class GCTFrames:
    """In-memory representation of a GCT/GCTX dataset."""

    matrix: pd.DataFrame
    row_metadata: pd.DataFrame
    col_metadata: pd.DataFrame
    attributes: Mapping[str, Any]


@dataclass(frozen=True)
class _NormalizedGCT:
    matrix: pd.DataFrame
    row_metadata: pd.DataFrame
    col_metadata: pd.DataFrame


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _validate_identifier(value: str, *, kind: str) -> None:
    if not value:
        raise ValueError(f"{kind} identifiers must not be empty")
    if "\x00" in value:
        raise ValueError(f"{kind} identifier contains a NUL byte: {value!r}")


def _normalize_metadata(
    frame: pd.DataFrame | None,
    ids: pd.Index,
    *,
    dimension: str,
) -> pd.DataFrame:
    if frame is None:
        return pd.DataFrame(index=ids.copy())
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"{dimension} metadata must be a pandas DataFrame")

    result = frame.copy()
    result.index = result.index.map(_text)
    expected = list(ids)
    observed = list(result.index)
    if observed != expected:
        raise ValueError(
            f"{dimension} metadata index must exactly match matrix "
            f"{dimension} identifiers and ordering"
        )
    if result.columns.has_duplicates:
        duplicated = list(result.columns[result.columns.duplicated()].map(_text))
        raise ValueError(
            f"{dimension} metadata contains duplicate field names: {duplicated}"
        )
    result.columns = result.columns.map(_text)
    for field in result.columns:
        if not field:
            raise ValueError(f"{dimension} metadata field names must not be empty")
        if "/" in field or "\x00" in field:
            raise ValueError(
                f"{dimension} metadata field is not GCTX-compatible: {field!r}"
            )

    if "id" in result.columns:
        stored_ids = result["id"].map(
            lambda value: "" if pd.isna(value) else _text(value)
        )
        if list(stored_ids) != expected:
            raise ValueError(
                f"{dimension} metadata has an 'id' field that does not match its index"
            )
        result = result.drop(columns=["id"])
    return result


def normalize_gct_frames(
    matrix: pd.DataFrame,
    *,
    row_metadata: pd.DataFrame | None = None,
    col_metadata: pd.DataFrame | None = None,
) -> GCTFrames:
    """Validate and align a numeric matrix and its two metadata frames."""

    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas DataFrame")
    if matrix.empty or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("GCT matrices must contain at least one row and one column")

    normalized_matrix = matrix.copy()
    normalized_matrix.index = normalized_matrix.index.map(_text)
    normalized_matrix.columns = normalized_matrix.columns.map(_text)
    if normalized_matrix.index.has_duplicates:
        raise ValueError("matrix row identifiers must be unique")
    if normalized_matrix.columns.has_duplicates:
        raise ValueError("matrix column identifiers must be unique")
    for rid in normalized_matrix.index:
        _validate_identifier(rid, kind="row")
    for cid in normalized_matrix.columns:
        _validate_identifier(cid, kind="column")

    try:
        normalized_matrix = normalized_matrix.apply(pd.to_numeric, errors="raise")
        values = normalized_matrix.to_numpy(dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("GCT matrix values must be numeric") from exc
    if np.iscomplexobj(values):
        raise ValueError("GCT matrix values must be real numbers")

    rdesc = _normalize_metadata(row_metadata, normalized_matrix.index, dimension="row")
    cdesc = _normalize_metadata(
        col_metadata, normalized_matrix.columns, dimension="column"
    )
    return GCTFrames(
        matrix=normalized_matrix,
        row_metadata=rdesc,
        col_metadata=cdesc,
        attributes={},
    )


def _update_blob(hasher: Any, value: bytes) -> None:
    hasher.update(len(value).to_bytes(8, byteorder="little", signed=False))
    hasher.update(value)


def _update_string(hasher: Any, value: Any) -> None:
    _update_blob(hasher, _text(value).encode("utf-8"))


def _canonical_float_array(values: Any, dtype: np.dtype) -> np.ndarray:
    array = np.asarray(values, dtype=dtype).copy(order="C")
    array[array == 0] = 0
    if np.issubdtype(array.dtype, np.floating):
        array[np.isnan(array)] = np.nan
    return np.asarray(array, dtype=dtype.newbyteorder("<"), order="C")


def _metadata_storage_vector(series: pd.Series) -> tuple[str, np.ndarray]:
    dtype = series.dtype
    if pd.api.types.is_bool_dtype(dtype):
        values = np.asarray(
            [255 if pd.isna(value) else int(bool(value)) for value in series],
            dtype=np.uint8,
        )
        return "boolean", values
    if pd.api.types.is_integer_dtype(dtype) and not series.isna().any():
        return "integer", np.asarray(series, dtype="<i8")
    if pd.api.types.is_numeric_dtype(dtype):
        return "numeric", _canonical_float_array(series, np.dtype("<f8"))
    if pd.api.types.is_datetime64_any_dtype(dtype):
        values = np.asarray(
            ["-666" if pd.isna(value) else value.isoformat() for value in series],
            dtype=object,
        )
        return "datetime", values

    values = np.asarray(
        ["-666" if pd.isna(value) else _text(value) for value in series],
        dtype=object,
    )
    return "string", values


def _update_metadata_hash(hasher: Any, frame: pd.DataFrame) -> None:
    _update_string(hasher, len(frame.columns))
    for field in frame.columns:
        _update_string(hasher, field)
        storage_type, values = _metadata_storage_vector(frame[field])
        _update_string(hasher, storage_type)
        if storage_type in {"string", "datetime"}:
            for value in values:
                _update_string(hasher, value)
        else:
            _update_blob(hasher, np.ascontiguousarray(values).tobytes(order="C"))


def gct_content_hash(
    matrix: pd.DataFrame,
    *,
    row_metadata: pd.DataFrame | None = None,
    col_metadata: pd.DataFrame | None = None,
    matrix_dtype: str | np.dtype = "float64",
    serialization_contract: str | None = None,
    hash_algorithm: str = "auto",
) -> str:
    """Return tackle's deterministic hash of a complete canonical GCT payload."""

    contract = gct_io_contracts()["gctx"]
    if serialization_contract is None:
        serialization_contract = str(contract["writer_contract"])
    frames = normalize_gct_frames(
        matrix, row_metadata=row_metadata, col_metadata=col_metadata
    )
    dtype = np.dtype(matrix_dtype)
    if dtype not in {np.dtype("float32"), np.dtype("float64")}:
        raise ValueError("matrix_dtype must be float32 or float64")

    actual_algorithm, hasher = _new_content_hasher(hash_algorithm)
    _update_string(hasher, contract["hash_contract"])
    _update_string(hasher, actual_algorithm)
    _update_string(hasher, serialization_contract)
    _update_string(hasher, dtype.name)
    _update_string(hasher, frames.matrix.shape[0])
    _update_string(hasher, frames.matrix.shape[1])
    for rid in frames.matrix.index:
        _update_string(hasher, rid)
    for cid in frames.matrix.columns:
        _update_string(hasher, cid)
    values = _canonical_float_array(frames.matrix.to_numpy(), dtype)
    _update_blob(hasher, values.tobytes(order="C"))
    _update_metadata_hash(hasher, frames.row_metadata)
    _update_metadata_hash(hasher, frames.col_metadata)
    return hasher.hexdigest()


def content_addressed_path(
    path: str | Path,
    content_hash: str,
    *,
    suffix: str | None = None,
) -> Path:
    """Insert a visible hash component immediately before a file suffix."""

    output = Path(path).expanduser()
    if suffix is not None:
        suffix = suffix if suffix.startswith(".") else f".{suffix}"
        if output.suffix.lower() != suffix.lower():
            output = output.with_suffix(suffix)
    clean_stem = _HASHED_STEM_RE.sub("", output.stem)
    hash_chars = int(gct_io_contracts()["gctx"]["filename_hash_chars"])
    return output.with_name(
        f"{clean_stem}.{content_hash[:hash_chars]}{output.suffix}"
    )


def _decode_attribute(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_gctx_content_hash_info(path: str | Path) -> dict[str, str] | None:
    """Read tackle's current or legacy hash and algorithm from the HDF5 header."""

    import h5py

    with h5py.File(Path(path), "r") as handle:
        value = handle.attrs.get(_gctx_attr("hash"))
        algorithm = handle.attrs.get(_gctx_attr("hash_algorithm"))
        if value is None:
            value = handle.attrs.get(_gctx_attr("legacy_sha256"))
            if value is not None:
                algorithm = "sha256"
    if value is None:
        return None
    return {
        "algorithm": _text(_decode_attribute(algorithm or "unknown")),
        "hash": _text(_decode_attribute(value)),
    }


def read_gctx_content_hash(path: str | Path) -> str | None:
    """Read tackle's current or legacy content hash without loading the matrix."""

    info = read_gctx_content_hash_info(path)
    return info["hash"] if info else None


def gctx_provenance_spec() -> dict[str, Any]:
    """Return the machine-readable specification for tackle-authored GCTX files."""

    contract = gct_io_contracts()["gctx"]
    return {
        "name": "tackle GCTX provenance",
        "schema_version": contract["provenance_schema"],
        "gctx_layout": {
            "version": "GCTX1.0",
            "matrix": "/0/DATA/0/matrix",
            "matrix_physical_orientation": "columns_by_rows",
            "row_ids": "/0/META/ROW/id",
            "column_ids": "/0/META/COL/id",
            "row_metadata": "/0/META/ROW/<field>",
            "column_metadata": "/0/META/COL/<field>",
        },
        "content_hash": {
            "policy": "auto",
            "preferred_algorithm": "blake3",
            "fallback_algorithm": "blake2b-256",
            "selected_algorithm": selected_content_hash_algorithm(),
            "supported_algorithms": ["blake3", "blake2b-256", "sha256"],
            "digest_bits": 256,
            "filename_prefix_hex_characters": contract["filename_hash_chars"],
            "purpose": "content addressing and cache invalidation",
            "canonical_payload_order": [
                "hash contract",
                "actual content hash algorithm",
                "serialization/writer contract",
                "matrix storage dtype",
                "logical matrix row count and column count",
                "ordered row identifiers",
                "ordered column identifiers",
                "logical row-major little-endian matrix bytes",
                "ordered row metadata field names, storage types, and values",
                "ordered column metadata field names, storage types, and values",
            ],
            "numeric_normalization": [
                "negative zero is normalized to positive zero",
                "NaN values are normalized before hashing",
            ],
            "not_hashed": [
                "output path and filename",
                "HDF5 chunk geometry",
                "HDF5 compression level",
                "file modification time",
                "physical HDF5 byte layout",
            ],
            "collision_handling": (
                "The filename carries an 80-bit prefix. Reuse additionally requires "
                "the recorded algorithm and complete 256-bit root-attribute hash to "
                "match; a prefix collision raises instead of silently reusing another "
                "payload."
            ),
        },
        "root_attributes": [
            {
                "name": "version",
                "value": "GCTX1.0",
                "description": "GCTX container version.",
            },
            {
                "name": "src",
                "description": "Absolute output path recorded when written.",
            },
            {
                "name": _gctx_attr("schema"),
                "value": contract["provenance_schema"],
                "description": "Version of this tackle provenance attribute schema.",
            },
            {
                "name": _gctx_attr("hash_algorithm"),
                "description": (
                    "Actual algorithm used for this file: blake3, blake2b-256, "
                    "or explicitly requested sha256."
                ),
            },
            {
                "name": _gctx_attr("hash"),
                "description": "Complete 256-bit canonical payload digest in hexadecimal.",
            },
            {
                "name": _gctx_attr("hash_contract"),
                "value": contract["hash_contract"],
                "description": "Canonical logical-payload hashing contract.",
            },
            {
                "name": _gctx_attr("writer_contract"),
                "value": contract["writer_contract"],
                "description": "GCTX serialization contract included in the digest.",
            },
            {
                "name": _gctx_attr("matrix_dtype"),
                "description": "Logical matrix storage dtype, normally float64.",
            },
            {
                "name": _gctx_attr("row_metadata_order"),
                "description": "JSON array preserving row metadata field order.",
            },
            {
                "name": _gctx_attr("column_metadata_order"),
                "description": "JSON array preserving column metadata field order.",
            },
        ],
        "dataset_attributes": {
            "tackle_storage_type": (
                "Logical metadata type: string, datetime, boolean, integer, or numeric."
            ),
            "tackle_string_encoding": "UTF-8 for fixed-width textual datasets.",
        },
        "legacy": {
            "content_hash_attribute": _gctx_attr("legacy_sha256"),
            "algorithm": "sha256",
            "inspection_support": True,
        },
    }


def _inspection_value(value: Any) -> Any:
    value = _decode_attribute(value)
    if isinstance(value, np.ndarray):
        return [_inspection_value(item) for item in value.tolist()]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return _text(value)


def _dataset_inspection(dataset: Any) -> dict[str, Any]:
    return {
        "path": dataset.name,
        "shape": [int(value) for value in dataset.shape],
        "dtype": str(dataset.dtype),
        "compression": dataset.compression,
        "compression_options": _inspection_value(dataset.compression_opts),
        "shuffle": bool(dataset.shuffle),
        "chunks": (
            [int(value) for value in dataset.chunks]
            if dataset.chunks is not None
            else None
        ),
    }


def _metadata_field_order(group: Any, encoded_order: Any) -> list[str]:
    available = [str(field) for field in group.keys() if str(field) != "id"]
    requested: list[str] = []
    if encoded_order is not None and encoded_order != "":
        try:
            parsed = json.loads(_text(_decode_attribute(encoded_order)))
            if isinstance(parsed, list):
                requested = [str(field) for field in parsed]
        except (TypeError, ValueError, json.JSONDecodeError):
            requested = []
    ordered = [field for field in requested if field in available]
    ordered.extend(sorted(field for field in available if field not in ordered))
    return ordered


def inspect_gctx(path: str | Path) -> dict[str, Any]:
    """Inspect a GCTX header and tackle provenance without loading matrix values."""

    import h5py

    source = Path(path).expanduser().resolve()
    if not source.exists() or not source.is_file():
        raise FileNotFoundError(str(source))

    with h5py.File(source, "r") as handle:
        try:
            matrix_dataset = handle["/0/DATA/0/matrix"]
            row_group = handle["/0/META/ROW"]
            col_group = handle["/0/META/COL"]
            row_id_dataset = row_group["id"]
            col_id_dataset = col_group["id"]
        except KeyError as exc:
            raise ValueError(
                f"Not a supported GCTX layout; missing HDF5 object: {exc}"
            ) from exc

        attributes = {
            str(key): _inspection_value(value) for key, value in handle.attrs.items()
        }
        row_fields = _metadata_field_order(
            row_group, attributes.get(_gctx_attr("row_metadata_order"))
        )
        col_fields = _metadata_field_order(
            col_group, attributes.get(_gctx_attr("column_metadata_order"))
        )
        matrix_info = _dataset_inspection(matrix_dataset)
        row_id_info = _dataset_inspection(row_id_dataset)
        col_id_info = _dataset_inspection(col_id_dataset)

    n_rows = int(row_id_info["shape"][0])
    n_columns = int(col_id_info["shape"][0])
    expected_physical_shape = [n_columns, n_rows]
    warnings: list[str] = []
    if matrix_info["shape"] != expected_physical_shape:
        warnings.append(
            "Physical matrix shape does not match column-by-row GCTX orientation: "
            f"observed {matrix_info['shape']}, expected {expected_physical_shape}."
        )

    current_hash = attributes.get(_gctx_attr("hash"))
    legacy_hash = attributes.get(_gctx_attr("legacy_sha256"))
    content_hash = current_hash or legacy_hash
    algorithm = attributes.get(_gctx_attr("hash_algorithm"))
    if not algorithm and legacy_hash:
        algorithm = "sha256"
    writer_contract = attributes.get(_gctx_attr("writer_contract"))
    provenance_available = bool(content_hash or writer_contract)

    filename_match = re.search(r"\.([0-9a-f]{20,64})$", source.stem)
    filename_hash_prefix = filename_match.group(1) if filename_match else None
    filename_hash_matches = (
        str(content_hash).startswith(filename_hash_prefix)
        if content_hash and filename_hash_prefix
        else None
    )

    tackle_attributes = {
        key: value for key, value in attributes.items() if key.startswith("tackle_")
    }
    provenance = {
        "available": provenance_available,
        "message": (
            "Tackle provenance is available."
            if provenance_available
            else "No tackle provenance available; this file may have been written by cmapR or another GCTX writer."
        ),
        "legacy": bool(legacy_hash and not current_hash),
        "schema_version": attributes.get(_gctx_attr("schema")),
        "content_hash_algorithm": algorithm,
        "content_hash": content_hash,
        "hash_contract": attributes.get(_gctx_attr("hash_contract")),
        "writer_contract": writer_contract,
        "matrix_dtype": attributes.get(_gctx_attr("matrix_dtype")),
        "filename_hash_prefix": filename_hash_prefix,
        "filename_hash_matches": filename_hash_matches,
        "attributes": tackle_attributes,
    }

    stat = source.stat()
    return {
        "inspection_schema_version": 1,
        "path": str(source),
        "name": source.name,
        "file_size_bytes": int(stat.st_size),
        "gctx_version": attributes.get("version"),
        "logical_shape": {"rows": n_rows, "columns": n_columns},
        "matrix_dataset": matrix_info,
        "row_ids_dataset": row_id_info,
        "column_ids_dataset": col_id_info,
        "row_metadata_fields": row_fields,
        "column_metadata_fields": col_fields,
        "root_attributes": attributes,
        "provenance": provenance,
        "warnings": warnings,
    }


def _human_bytes(value: int) -> str:
    size = float(value)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{value} B"


def format_gctx_provenance_spec(spec: Mapping[str, Any] | None = None) -> str:
    """Pretty-format the tackle GCTX provenance specification."""

    payload = dict(spec or gctx_provenance_spec())
    hash_spec = payload["content_hash"]
    lines = [
        "Tackle GCTX provenance specification",
        f"  Schema version: {payload['schema_version']}",
        f"  GCTX layout: {payload['gctx_layout']['version']}",
        (
            "  Content hash: "
            f"auto -> {hash_spec['selected_algorithm']} in this environment "
            f"({hash_spec['digest_bits']} bits; "
            f"{hash_spec['filename_prefix_hex_characters']} hex characters in filenames)"
        ),
        (
            "  Policy: prefer "
            f"{hash_spec['preferred_algorithm']}; fall back to "
            f"{hash_spec['fallback_algorithm']}; sha256 is explicitly supported"
        ),
        "  Canonical payload, in order:",
    ]
    lines.extend(
        f"    {index}. {value}"
        for index, value in enumerate(hash_spec["canonical_payload_order"], 1)
    )
    lines.append("  Excluded from the content hash:")
    lines.extend(f"    - {value}" for value in hash_spec["not_hashed"])
    lines.append("  Root attributes:")
    for item in payload["root_attributes"]:
        expected = f" = {item['value']}" if "value" in item else ""
        lines.append(f"    - {item['name']}{expected}: {item['description']}")
    lines.append(f"  Collision handling: {hash_spec['collision_handling']}")
    return "\n".join(lines)


def format_gctx_inspection(
    inspection: Mapping[str, Any], *, include_spec: bool = False
) -> str:
    """Pretty-format an :func:`inspect_gctx` result for terminal display."""

    shape = inspection["logical_shape"]
    matrix = inspection["matrix_dataset"]
    provenance = inspection["provenance"]
    lines = [
        "GCTX inspection",
        f"  Path: {inspection['path']}",
        f"  File size: {_human_bytes(int(inspection['file_size_bytes']))}",
        f"  GCTX version: {inspection.get('gctx_version') or 'unknown'}",
        f"  Logical matrix: {shape['rows']} rows x {shape['columns']} columns",
        f"  Matrix dtype: {matrix['dtype']}",
        f"  Matrix compression: {matrix['compression'] or 'none'}",
        (
            "  Row metadata fields "
            f"({len(inspection['row_metadata_fields'])}): "
            + (", ".join(inspection["row_metadata_fields"]) or "none")
        ),
        (
            "  Column metadata fields "
            f"({len(inspection['column_metadata_fields'])}): "
            + (", ".join(inspection["column_metadata_fields"]) or "none")
        ),
        "",
        "Tackle provenance",
    ]
    if provenance["available"]:
        lines.extend(
            [
                f"  Status: available{' (legacy)' if provenance['legacy'] else ''}",
                f"  Schema version: {provenance.get('schema_version') or 'legacy/unset'}",
                f"  Content hash algorithm: {provenance.get('content_hash_algorithm') or 'unknown'}",
                f"  Content hash: {provenance.get('content_hash') or 'unavailable'}",
                f"  Hash contract: {provenance.get('hash_contract') or 'unavailable'}",
                f"  Writer contract: {provenance.get('writer_contract') or 'unavailable'}",
                (
                    "  Filename hash prefix matches: "
                    + (
                        "yes"
                        if provenance.get("filename_hash_matches") is True
                        else "no"
                        if provenance.get("filename_hash_matches") is False
                        else "not encoded in filename"
                    )
                ),
            ]
        )
    else:
        lines.append(f"  Status: {provenance['message']}")

    for warning in inspection.get("warnings", []):
        lines.append(f"  Warning: {warning}")
    if include_spec:
        lines.extend(["", format_gctx_provenance_spec()])
    return "\n".join(lines)


def _gctx_chunks(
    shape: tuple[int, int], dtype: np.dtype, max_chunk_kb: int
) -> tuple[int, int]:
    physical_rows, physical_cols = shape
    target_elements = max(1, int(max_chunk_kb) * 1024 // dtype.itemsize)
    row_chunk = max(1, min(physical_rows, 1000))
    col_chunk = max(1, min(physical_cols, target_elements // row_chunk))
    return row_chunk, col_chunk


def _write_hdf5_vector(
    group: Any,
    name: str,
    series: pd.Series | Iterable[Any],
    *,
    compression_level: int,
) -> None:
    if not isinstance(series, pd.Series):
        values = np.asarray([_text(value) for value in series], dtype=object)
        storage_type = "string"
    else:
        storage_type, values = _metadata_storage_vector(series)

    if storage_type in {"string", "datetime"}:
        # HDF5 stores variable-length strings in a separate heap whose payload
        # is not compressed by the dataset's gzip filter. Fixed-width UTF-8
        # bytes keep metadata-rich GCTX files compact and remain readable by
        # both cmapR/rhdf5 and cmapPy/h5py.
        encoded = []
        for value in values:
            raw = _text(value).encode("utf-8")
            if b"\x00" in raw:
                raise ValueError(f"GCTX metadata field {name!r} contains a NUL byte")
            encoded.append(raw)
        width = max(1, max(map(len, encoded), default=1))
        array = np.asarray(encoded, dtype=f"S{width}")
    else:
        array = np.asarray(values)

    kwargs: dict[str, Any] = {}
    compression_min_bytes = int(
        gct_io_contracts()["gctx"]["compression_min_bytes"]
    )
    if int(compression_level) > 0 and array.nbytes >= compression_min_bytes:
        kwargs.update(
            compression="gzip",
            compression_opts=int(compression_level),
            shuffle=True,
        )
    dataset = group.create_dataset(name, data=array, **kwargs)
    if storage_type in {"string", "datetime"}:
        dataset.attrs["tackle_string_encoding"] = "utf-8"
    dataset.attrs["tackle_storage_type"] = storage_type


def write_gctx(
    matrix: pd.DataFrame,
    path: str | Path,
    *,
    row_metadata: pd.DataFrame | None = None,
    col_metadata: pd.DataFrame | None = None,
    matrix_dtype: str | np.dtype = "float64",
    hash_algorithm: str = "auto",
    compression_level: int = 6,
    max_chunk_kb: int = 1024,
    content_addressed: bool = True,
    overwrite: bool = False,
) -> Path:
    """Write a cmapR-compatible GCTX and return its exact path.

    Content-addressed files are named ``stem.<20 hex chars>.gctx``. The full
    256-bit digest and actual algorithm are stored as root HDF5 attributes,
    allowing a matching existing file to be reused after reading only its header.
    """

    import h5py

    contract = gct_io_contracts()["gctx"]
    frames = normalize_gct_frames(
        matrix, row_metadata=row_metadata, col_metadata=col_metadata
    )
    dtype = np.dtype(matrix_dtype)
    if dtype not in {np.dtype("float32"), np.dtype("float64")}:
        raise ValueError("matrix_dtype must be float32 or float64")
    if not 0 <= int(compression_level) <= 9:
        raise ValueError("compression_level must be between 0 and 9")

    actual_hash_algorithm, _ = _new_content_hasher(hash_algorithm)
    digest = gct_content_hash(
        frames.matrix,
        row_metadata=frames.row_metadata,
        col_metadata=frames.col_metadata,
        matrix_dtype=dtype,
        serialization_contract=str(contract["writer_contract"]),
        hash_algorithm=actual_hash_algorithm,
    )
    requested = Path(path).expanduser()
    if requested.suffix.lower() != ".gctx":
        requested = requested.with_suffix(".gctx")
    output = (
        content_addressed_path(requested, digest) if content_addressed else requested
    )
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    if output.exists() and not overwrite:
        existing_hash = read_gctx_content_hash_info(output)
        if existing_hash == {
            "algorithm": actual_hash_algorithm,
            "hash": digest,
        }:
            return output
        raise FileExistsError(
            f"Existing GCTX hash does not match its requested payload: {output}"
        )

    temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
    try:
        with h5py.File(temporary, "w") as handle:
            handle.attrs["version"] = np.bytes_("GCTX1.0")
            handle.attrs["src"] = str(output)
            handle.attrs[_gctx_attr("schema")] = contract["provenance_schema"]
            handle.attrs[_gctx_attr("hash")] = digest
            handle.attrs[_gctx_attr("hash_algorithm")] = actual_hash_algorithm
            handle.attrs[_gctx_attr("hash_contract")] = contract["hash_contract"]
            handle.attrs[_gctx_attr("writer_contract")] = contract[
                "writer_contract"
            ]
            handle.attrs[_gctx_attr("matrix_dtype")] = dtype.name
            handle.attrs[_gctx_attr("row_metadata_order")] = json.dumps(
                list(frames.row_metadata.columns), separators=(",", ":")
            )
            handle.attrs[_gctx_attr("column_metadata_order")] = json.dumps(
                list(frames.col_metadata.columns), separators=(",", ":")
            )

            data_group = handle.require_group("/0/DATA/0")
            row_group = handle.require_group("/0/META/ROW")
            col_group = handle.require_group("/0/META/COL")

            physical = _canonical_float_array(frames.matrix.to_numpy(), dtype).T
            matrix_kwargs: dict[str, Any] = {}
            if compression_level and physical.nbytes >= int(
                contract["compression_min_bytes"]
            ):
                matrix_kwargs.update(
                    chunks=_gctx_chunks(physical.shape, dtype, max_chunk_kb),
                    compression="gzip",
                    compression_opts=int(compression_level),
                    shuffle=True,
                )
            data_group.create_dataset(
                "matrix", data=physical, dtype=dtype, **matrix_kwargs
            )

            _write_hdf5_vector(
                row_group,
                "id",
                frames.matrix.index,
                compression_level=int(compression_level),
            )
            _write_hdf5_vector(
                col_group,
                "id",
                frames.matrix.columns,
                compression_level=int(compression_level),
            )
            for field in frames.row_metadata.columns:
                _write_hdf5_vector(
                    row_group,
                    field,
                    frames.row_metadata[field],
                    compression_level=int(compression_level),
                )
            for field in frames.col_metadata.columns:
                _write_hdf5_vector(
                    col_group,
                    field,
                    frames.col_metadata[field],
                    compression_level=int(compression_level),
                )
            handle.flush()
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return output


def _read_hdf5_vector(dataset: Any) -> pd.Series:
    values = dataset[...]
    storage_type = _text(
        _decode_attribute(dataset.attrs.get("tackle_storage_type", ""))
    )
    if values.dtype.kind in {"O", "S", "U"}:
        decoded = [
            value.decode("utf-8") if isinstance(value, bytes) else _text(value)
            for value in values
        ]
        return pd.Series([pd.NA if value == "-666" else value for value in decoded])
    if storage_type == "boolean":
        return pd.Series(
            [pd.NA if int(value) == 255 else bool(value) for value in values],
            dtype="boolean",
        )
    return pd.Series(values)


def read_gctx(path: str | Path) -> GCTFrames:
    """Read a GCTX written by tackle without depending on cmapPy."""

    import h5py

    source = Path(path).expanduser().resolve()
    with h5py.File(source, "r") as handle:
        row_group = handle["/0/META/ROW"]
        col_group = handle["/0/META/COL"]
        rids = list(_read_hdf5_vector(row_group["id"]).astype(str))
        cids = list(_read_hdf5_vector(col_group["id"]).astype(str))
        values = np.asarray(handle["/0/DATA/0/matrix"][...]).T

        row_order = json.loads(
            _text(
                _decode_attribute(
                    handle.attrs.get(_gctx_attr("row_metadata_order"), "[]")
                )
            )
        )
        col_order = json.loads(
            _text(
                _decode_attribute(
                    handle.attrs.get(_gctx_attr("column_metadata_order"), "[]")
                )
            )
        )
        if not row_order:
            row_order = [field for field in row_group.keys() if field != "id"]
        if not col_order:
            col_order = [field for field in col_group.keys() if field != "id"]

        rdesc = pd.DataFrame(
            {
                field: _read_hdf5_vector(row_group[field]).to_numpy()
                for field in row_order
            },
            index=rids,
        )
        cdesc = pd.DataFrame(
            {
                field: _read_hdf5_vector(col_group[field]).to_numpy()
                for field in col_order
            },
            index=cids,
        )
        attributes = {
            _text(key): _decode_attribute(value) for key, value in handle.attrs.items()
        }
    matrix = pd.DataFrame(values, index=rids, columns=cids)
    return GCTFrames(
        matrix=matrix,
        row_metadata=rdesc,
        col_metadata=cdesc,
        attributes=attributes,
    )


def _gct_text(value: Any, *, missing: str = "") -> str:
    if pd.isna(value):
        return missing
    return _text(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_gct(
    matrix: pd.DataFrame,
    path: str | Path,
    *,
    row_metadata: pd.DataFrame | None = None,
    col_metadata: pd.DataFrame | None = None,
    precision: int = 4,
    content_addressed: bool = False,
    overwrite: bool = False,
    chunk_rows: int = 10_000,
) -> Path:
    """Write textual GCT 1.3 through one persistent buffered handle."""

    frames = normalize_gct_frames(
        matrix, row_metadata=row_metadata, col_metadata=col_metadata
    )
    if int(precision) < 0:
        raise ValueError("precision must be non-negative")
    precision = int(precision)
    text_writer_contract = str(gct_io_contracts()["gct_text_writer"])
    digest = gct_content_hash(
        frames.matrix,
        row_metadata=frames.row_metadata,
        col_metadata=frames.col_metadata,
        matrix_dtype="float64",
        serialization_contract=f"{text_writer_contract}:precision={precision}",
    )
    requested = Path(path).expanduser()
    if requested.suffix.lower() != ".gct":
        requested = requested.with_suffix(".gct")
    output = (
        content_addressed_path(requested, digest) if content_addressed else requested
    ).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not overwrite:
        if content_addressed:
            return output
        raise FileExistsError(f"GCT already exists: {output}")

    temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
    float_format = f"%.{precision}f"
    nr, nc = frames.matrix.shape
    nrdesc = frames.row_metadata.shape[1]
    ncdesc = frames.col_metadata.shape[1]
    try:
        with temporary.open(
            "w", encoding="utf-8", newline="", buffering=1024 * 1024
        ) as handle:
            handle.write("#1.3\n")
            handle.write(f"{nr}\t{nc}\t{nrdesc}\t{ncdesc}\n")
            handle.write(
                "\t".join(
                    [
                        "id",
                        *[_gct_text(field) for field in frames.row_metadata.columns],
                        *[_gct_text(cid) for cid in frames.matrix.columns],
                    ]
                )
                + "\n"
            )
            for field in frames.col_metadata.columns:
                values = [
                    _gct_text(value, missing="NA")
                    for value in frames.col_metadata[field]
                ]
                handle.write(
                    "\t".join(
                        [
                            _gct_text(field),
                            *(["na"] * nrdesc),
                            *values,
                        ]
                    )
                    + "\n"
                )

            for start in range(0, nr, max(1, int(chunk_rows))):
                stop = min(nr, start + max(1, int(chunk_rows)))
                ids = pd.Series(frames.matrix.index[start:stop], name="id")
                rdesc = frames.row_metadata.iloc[start:stop].reset_index(drop=True)
                for field in rdesc.columns:
                    if not pd.api.types.is_numeric_dtype(rdesc[field]):
                        rdesc[field] = rdesc[field].map(
                            lambda value: _gct_text(value, missing="")
                        )
                block = pd.concat(
                    [
                        ids.reset_index(drop=True),
                        rdesc,
                        frames.matrix.iloc[start:stop].reset_index(drop=True),
                    ],
                    axis=1,
                )
                block.to_csv(
                    handle,
                    sep="\t",
                    header=False,
                    index=False,
                    na_rep="NA",
                    float_format=float_format,
                    lineterminator="\n",
                )
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return output


def dataframe_content_hash(
    frame: pd.DataFrame,
    *,
    serialization_contract: str | None = None,
    float_format: str = "%.17g",
    hash_algorithm: str = "auto",
) -> str:
    """Hash the exact canonical TSV representation of a DataFrame."""

    if serialization_contract is None:
        serialization_contract = str(gct_io_contracts()["tsv_hash"])
    stable = frame.copy()
    stable.index = stable.index.map(_text)
    stable.columns = stable.columns.map(_text)
    payload = stable.to_csv(
        sep="\t",
        index=True,
        na_rep="NA",
        float_format=float_format,
        lineterminator="\n",
    ).encode("utf-8")
    actual_algorithm, hasher = _new_content_hasher(hash_algorithm)
    _update_string(hasher, serialization_contract)
    _update_string(hasher, actual_algorithm)
    _update_blob(hasher, payload)
    return hasher.hexdigest()


def write_hashed_tsv(
    frame: pd.DataFrame,
    path: str | Path,
    *,
    float_format: str = "%.17g",
    hash_algorithm: str = "auto",
) -> Path:
    """Write a small content-addressed TSV, reusing an identical existing file."""

    tsv_contract = str(gct_io_contracts()["tsv_hash"])
    stable = frame.copy()
    stable.index = stable.index.map(_text)
    stable.columns = stable.columns.map(_text)
    payload = stable.to_csv(
        sep="\t",
        index=True,
        na_rep="NA",
        float_format=float_format,
        lineterminator="\n",
    ).encode("utf-8")
    actual_algorithm, hasher = _new_content_hasher(hash_algorithm)
    _update_string(hasher, tsv_contract)
    _update_string(hasher, actual_algorithm)
    _update_blob(hasher, payload)
    digest = hasher.hexdigest()
    requested = Path(path).expanduser()
    if requested.suffix.lower() != ".tsv":
        requested = requested.with_suffix(".tsv")
    output = content_addressed_path(requested, digest).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        return output
    temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_bytes(payload)
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return output
