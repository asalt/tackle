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
  label = TRUE,
  label_size = 4.25,
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
stopifnot(is.null(result[['plots']][['pc1_vs_pc2']]$labels$tag))
stopifnot(identical(result[['plots']][['pc1_vs_pc2']]$coordinates$clip, 'off'))
note_spec <- attr(
  result[['plots']][['pc1_vs_pc2']],
  'tackle_pca_feature_note',
  exact = TRUE
)
stopifnot(identical(note_spec$labels[[1]], quote(bold(n == 3L))))
stopifnot(identical(note_spec$color, 'grey35'))
stopifnot(identical(note_spec$fontsize, 13))
label_layers <- Filter(
  function(layer) inherits(layer$geom, 'GeomTextRepel'),
  result[['plots']][['pc1_vs_pc2']]$layers
)
stopifnot(length(label_layers) == 1L)
stopifnot(identical(label_layers[[1]]$aes_params$size, 4.25))
stopifnot(identical(label_layers[[1]]$aes_params$fontface, 'bold'))

d_filter <- data.frame(
  GeneID = rep(c('g_high', 'g_missing', 'g_complete', 'g_low'), each = 4),
  variable = rep(c('A', 'B', 'C', 'D'), 4),
  value = c(
    0, 20, 5, 15,
    0, 10, NA, NA,
    0, 8, 0, 8,
    0, 1, 0, 1
  ),
  group = rep(c('x', 'x', 'y', 'y'), 4)
)
result_avg <- pca2(
  d_filter,
  color = 'group',
  center = FALSE,
  scale = FALSE,
  fillna = 'avg',
  top_variable_pct = 50,
  max_pc = 2,
  export_tables = FALSE,
  return_plots = TRUE,
  return_data = TRUE
)
result_min <- pca2(
  d_filter,
  color = 'group',
  center = FALSE,
  scale = FALSE,
  fillna = 'min',
  top_variable_pct = 50,
  max_pc = 2,
  export_tables = FALSE,
  return_plots = TRUE,
  return_data = TRUE
)
stopifnot(identical(dim(result_avg[['pca_mat']]), c(4L, 2L)))
stopifnot(setequal(colnames(result_avg[['pca_mat']]), c('g_high', 'g_complete')))
stopifnot(setequal(colnames(result_min[['pca_mat']]), c('g_high', 'g_missing')))
stopifnot(identical(as.numeric(result_min[['pca_mat']][, 'g_missing']), c(0, 10, 0, 0)))
avg_note_spec <- attr(
  result_avg[['plots']][['pc1_vs_pc2']],
  'tackle_pca_feature_note',
  exact = TRUE
)
stopifnot(identical(avg_note_spec$labels[[1]], quote(bold(n == 2L))))
stopifnot(identical(avg_note_spec$labels[[2]], quote(bold(Var >= P[50]))))
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
