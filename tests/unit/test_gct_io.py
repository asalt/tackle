from __future__ import annotations

import builtins
import importlib.util
import re

import numpy as np
import pandas as pd
import pytest

from tackle import gct_io
from tackle.gct_io import (
    format_gctx_inspection,
    format_gctx_provenance_spec,
    gct_content_hash,
    gct_io_contracts,
    gctx_provenance_spec,
    inspect_gctx,
    read_gctx,
    selected_content_hash_algorithm,
    write_gct,
    write_gctx,
    write_hashed_tsv,
)


def _frames():
    rids = [
        "1017",
        "ENSG00000141510",
        "TP53",
        "gene_00047",
        "00123",
        "chr1:100-200",
        "A/B",
    ]
    cids = ["S1", "S2", "sample_03"]
    values = (
        np.arange(1, len(rids) * len(cids) + 1, dtype=float).reshape(
            len(rids), len(cids)
        )
        / 7
    )
    values[1, 2] = np.nan
    matrix = pd.DataFrame(values, index=rids, columns=cids)
    rdesc = pd.DataFrame(
        {
            "GeneSymbol": [
                "GENE1017",
                "TP53",
                "TP53_alias",
                "G47",
                "LEADING",
                "LOC",
                "SLASH",
            ],
            "numeric_meta": np.arange(1, len(rids) + 1),
            "id": rids,
        },
        index=rids,
    )
    cdesc = pd.DataFrame(
        {
            "group": ["A", "A", "B"],
            "batch": [1, 2, 2],
            "id": cids,
        },
        index=cids,
    )
    return matrix, rdesc, cdesc


def test_contract_registry_is_callable_and_keeps_constants_out_of_module_namespace():
    contracts = gct_io_contracts()

    assert contracts["gctx"]["writer_contract"] == "tackle-gctx-writer-v4"
    assert contracts["tsv_hash"] == "tackle-tsv-v2"
    assert not hasattr(gct_io, "GCTX_CONTRACT")
    assert not hasattr(gct_io, "TSV_HASH_CONTRACT")


def test_write_gctx_roundtrips_string_ids_matrix_and_metadata(tmp_path):
    matrix, rdesc, cdesc = _frames()

    output = write_gctx(
        matrix,
        tmp_path / "mixed.gctx",
        row_metadata=rdesc,
        col_metadata=cdesc,
    )
    parsed = read_gctx(output)

    assert re.fullmatch(r"mixed\.[0-9a-f]{20}\.gctx", output.name)
    assert list(parsed.matrix.index) == list(matrix.index)
    assert list(parsed.matrix.columns) == list(matrix.columns)
    assert np.allclose(parsed.matrix, matrix, equal_nan=True)
    assert list(parsed.row_metadata.columns) == ["GeneSymbol", "numeric_meta"]
    assert list(parsed.col_metadata.columns) == ["group", "batch"]
    assert parsed.row_metadata["GeneSymbol"].tolist() == rdesc["GeneSymbol"].tolist()
    assert parsed.col_metadata["group"].tolist() == cdesc["group"].tolist()
    attrs = gct_io_contracts()["gctx"]["attrs"]
    assert parsed.attributes[attrs["matrix_dtype"]] == "float64"
    assert (
        parsed.attributes[attrs["hash_algorithm"]] == selected_content_hash_algorithm()
    )
    assert parsed.attributes[attrs["hash"]] == gct_content_hash(
        matrix,
        row_metadata=rdesc,
        col_metadata=cdesc,
    )


def test_write_gctx_reuses_matching_hash_without_replacing_file(tmp_path, monkeypatch):
    matrix, rdesc, cdesc = _frames()
    output = write_gctx(
        matrix,
        tmp_path / "reused.gctx",
        row_metadata=rdesc,
        col_metadata=cdesc,
    )

    def unexpected_replace(*_args, **_kwargs):
        raise AssertionError("an identical content-addressed GCTX was rewritten")

    monkeypatch.setattr(gct_io.os, "replace", unexpected_replace)
    reused = write_gctx(
        matrix,
        tmp_path / "reused.gctx",
        row_metadata=rdesc,
        col_metadata=cdesc,
    )

    assert reused == output


def test_gctx_hash_changes_with_values_metadata_and_storage_dtype():
    matrix, rdesc, cdesc = _frames()
    baseline = gct_content_hash(matrix, row_metadata=rdesc, col_metadata=cdesc)
    changed_matrix = matrix.copy()
    changed_matrix.iloc[0, 0] += 0.25
    changed_rdesc = rdesc.copy()
    changed_rdesc.loc["TP53", "GeneSymbol"] = "TP53_changed"

    assert (
        gct_content_hash(changed_matrix, row_metadata=rdesc, col_metadata=cdesc)
        != baseline
    )
    assert (
        gct_content_hash(matrix, row_metadata=changed_rdesc, col_metadata=cdesc)
        != baseline
    )
    assert (
        gct_content_hash(
            matrix,
            row_metadata=rdesc,
            col_metadata=cdesc,
            matrix_dtype="float32",
        )
        != baseline
    )


def test_hash_algorithm_is_recorded_and_supports_fallback_and_sha256(tmp_path):
    matrix, rdesc, cdesc = _frames()
    algorithms = ["blake2b-256", "sha256"]
    if importlib.util.find_spec("blake3") is not None:
        algorithms.insert(0, "blake3")
    digests = {
        algorithm: gct_content_hash(
            matrix,
            row_metadata=rdesc,
            col_metadata=cdesc,
            hash_algorithm=algorithm,
        )
        for algorithm in algorithms
    }

    assert all(len(value) == 64 for value in digests.values())
    assert len(set(digests.values())) == len(algorithms)

    fallback = write_gctx(
        matrix,
        tmp_path / "fallback.gctx",
        row_metadata=rdesc,
        col_metadata=cdesc,
        hash_algorithm="blake2b-256",
    )
    parsed = read_gctx(fallback)
    attrs = gct_io_contracts()["gctx"]["attrs"]
    assert parsed.attributes[attrs["hash_algorithm"]] == "blake2b-256"


def test_auto_hash_falls_back_to_hashlib_blake2b(monkeypatch):
    real_import = builtins.__import__

    def import_without_blake3(name, *args, **kwargs):
        if name == "blake3":
            raise ImportError("simulated missing optional blake3")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_blake3)

    assert selected_content_hash_algorithm() == "blake2b-256"


def test_inspect_gctx_reports_layout_provenance_and_spec(tmp_path):
    matrix, rdesc, cdesc = _frames()
    output = write_gctx(
        matrix,
        tmp_path / "inspect.gctx",
        row_metadata=rdesc,
        col_metadata=cdesc,
    )

    inspection = inspect_gctx(output)
    spec = gctx_provenance_spec()
    rendered = format_gctx_inspection(inspection, include_spec=True)

    assert inspection["logical_shape"] == {"rows": 7, "columns": 3}
    assert inspection["row_metadata_fields"] == ["GeneSymbol", "numeric_meta"]
    assert inspection["column_metadata_fields"] == ["group", "batch"]
    assert inspection["provenance"]["available"] is True
    assert (
        inspection["provenance"]["content_hash_algorithm"]
        == selected_content_hash_algorithm()
    )
    assert inspection["provenance"]["filename_hash_matches"] is True
    assert spec["content_hash"]["fallback_algorithm"] == "blake2b-256"
    assert "Canonical payload, in order" in rendered
    assert "fall back to blake2b-256" in format_gctx_provenance_spec(spec)


def test_inspect_external_gctx_reports_no_tackle_provenance(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "external.gctx"
    with h5py.File(path, "w") as handle:
        handle.attrs["version"] = "GCTX1.0"
        handle.create_dataset("/0/DATA/0/matrix", data=np.eye(2))
        handle.create_dataset("/0/META/ROW/id", data=np.asarray([b"R1", b"R2"]))
        handle.create_dataset("/0/META/COL/id", data=np.asarray([b"C1", b"C2"]))

    inspection = inspect_gctx(path)
    rendered = format_gctx_inspection(inspection)

    assert inspection["logical_shape"] == {"rows": 2, "columns": 2}
    assert inspection["provenance"]["available"] is False
    assert "No tackle provenance available" in rendered


def test_gct_validation_rejects_silent_metadata_reindexing(tmp_path):
    matrix, rdesc, cdesc = _frames()

    with pytest.raises(ValueError, match="exactly match"):
        write_gctx(
            matrix,
            tmp_path / "bad.gctx",
            row_metadata=rdesc.iloc[::-1],
            col_metadata=cdesc,
        )


def test_text_gct_is_real_v13_and_cmap_py_can_parse_it(tmp_path):
    parse_gct = pytest.importorskip("cmapPy.pandasGEXpress.parse_gct")
    matrix, rdesc, cdesc = _frames()

    output = write_gct(
        matrix,
        tmp_path / "mixed.gct",
        row_metadata=rdesc,
        col_metadata=cdesc,
        precision=12,
    )
    parsed = parse_gct.parse(str(output))

    assert output.name == "mixed.gct"
    assert list(parsed.data_df.index) == list(matrix.index)
    assert list(parsed.data_df.columns) == list(matrix.columns)
    # cmapPy's text parser materializes the matrix as float32 even when the
    # file contains more precision; this test is about GCT 1.3 compatibility.
    assert np.allclose(parsed.data_df, matrix, equal_nan=True, rtol=0, atol=2e-7)
    assert parsed.row_metadata_df["GeneSymbol"].tolist() == rdesc["GeneSymbol"].tolist()
    assert parsed.col_metadata_df["group"].tolist() == cdesc["group"].tolist()


def test_write_hashed_tsv_reuses_identical_pairwise_counts(tmp_path, monkeypatch):
    counts = pd.DataFrame([[5, 4], [4, 6]], index=["S1", "S2"], columns=["S1", "S2"])
    output = write_hashed_tsv(counts, tmp_path / "pairwise_counts.tsv")

    def unexpected_replace(*_args, **_kwargs):
        raise AssertionError("an identical content-addressed TSV was rewritten")

    monkeypatch.setattr(gct_io.os, "replace", unexpected_replace)
    reused = write_hashed_tsv(counts, tmp_path / "pairwise_counts.tsv")

    assert reused == output
    assert re.fullmatch(r"pairwise_counts\.[0-9a-f]{20}\.tsv", output.name)
