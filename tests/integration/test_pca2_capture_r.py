import json
import shutil
import subprocess
from pathlib import Path

import pytest


def test_pca2_returns_exact_pre_prcomp_matrix_without_ordered_column():
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")

    required = ("tidyverse", "ggfortify", "ggrepel")
    package_check = "; ".join(
        f"stopifnot(requireNamespace('{package}', quietly=TRUE))"
        for package in required
    )
    check = subprocess.run(
        [rscript, "-e", package_check],
        text=True,
        capture_output=True,
    )
    if check.returncode != 0:
        pytest.skip("Required pca2 R packages are not installed")

    r_file = Path(__file__).resolve().parents[2] / "tackle" / "R" / "pcaplot.R"
    code = f"""
source({json.dumps(str(r_file))})
d <- data.frame(
  GeneID = rep(c('g1', 'g2', 'g3'), each = 4),
  variable = rep(c('A', 'B', 'C', 'D'), 3),
  value = c(1, 2, 3, 4, 4, 3, 2, 1, 1, 3, 2, 5),
  group = rep(c('x', 'x', 'y', 'y'), 3)
)
result <- pca2(
  d,
  color = 'group',
  center = TRUE,
  scale = FALSE,
  fillna = 'min',
  max_pc = 2,
  export_tables = FALSE,
  return_plots = TRUE,
  return_data = TRUE
)
stopifnot(identical(names(result), c('plots', 'pca_mat', 'scores', 'variance')))
stopifnot('scree' %in% names(result[['plots']]))
stopifnot(identical(dim(result[['pca_mat']]), c(4L, 3L)))
stopifnot(!'ordered' %in% colnames(result[['pca_mat']]))
stopifnot(identical(dim(result[['scores']]), c(4L, 3L)))
stopifnot(identical(rownames(result[['scores']]), rownames(result[['pca_mat']])))
stopifnot(identical(colnames(result[['scores']]), c('PC1', 'PC2', 'PC3')))
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
