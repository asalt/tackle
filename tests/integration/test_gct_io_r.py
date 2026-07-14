from __future__ import annotations

import json
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from tackle.gct_io import write_gctx


def test_tackle_gctx_with_root_hash_attributes_is_cmapr_compatible(tmp_path):
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")
    check = subprocess.run(
        [rscript, "-e", "stopifnot(requireNamespace('cmapR', quietly=TRUE))"],
        text=True,
        capture_output=True,
    )
    if check.returncode != 0:
        pytest.skip("cmapR is not installed")

    rids = ["00123", "ENSG00000141510", "TP53", "gene_00047", "A/B"]
    cids = ["S1", "sample_02"]
    matrix = pd.DataFrame(
        [[1.25, 2.5], [3.75, np.nan], [5.0, 6.25], [7.5, 8.75], [9.0, 10.5]],
        index=rids,
        columns=cids,
    )
    rdesc = pd.DataFrame(
        {"GeneSymbol": ["LEADING", "TP53", "TP53_alias", "G47", "SLASH"]},
        index=rids,
    )
    cdesc = pd.DataFrame({"group": ["A", "B"], "batch": [1, 2]}, index=cids)
    output = write_gctx(
        matrix,
        tmp_path / "compatibility.gctx",
        row_metadata=rdesc,
        col_metadata=cdesc,
    )

    code = f"""
ds <- cmapR::parse_gctx({json.dumps(str(output))})
stopifnot(identical(ds@rid, c('00123', 'ENSG00000141510', 'TP53', 'gene_00047', 'A/B')))
stopifnot(identical(ds@cid, c('S1', 'sample_02')))
stopifnot(identical(dim(ds@mat), c(5L, 2L)))
stopifnot(isTRUE(all.equal(ds@mat[1, 1], 1.25)))
stopifnot(is.na(ds@mat[2, 2]))
stopifnot(identical(as.character(ds@rdesc$GeneSymbol), c('LEADING', 'TP53', 'TP53_alias', 'G47', 'SLASH')))
stopifnot(identical(as.character(ds@cdesc$group), c('A', 'B')))
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
