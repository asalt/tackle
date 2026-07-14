import json
import shutil
import subprocess
from pathlib import Path

import pytest


def test_correlation_heatmap_uses_precomputed_dissimilarity_for_both_axes(tmp_path):
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")

    required = ("ComplexHeatmap", "circlize", "jsonlite")
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
        pytest.skip("Required correlation-plot R packages are not installed")

    r_file = Path(__file__).resolve().parents[2] / "tackle" / "R" / "correlationplot.R"
    outname = tmp_path / "correlation"
    code = f"""
source({json.dumps(str(r_file))})
ids <- c('S1', 'S2', 'S3', 'S4')
raw <- matrix(
  c(1, .95, -.7, -.4, .95, 1, -.5, -.2,
    -.7, -.5, 1, .90, -.4, -.2, .90, 1),
  4,
  4,
  dimnames = list(ids, ids)
)
distance <- matrix(
  c(0, 9, 1, 8,
    9, 0, 8, 1,
    1, 8, 0, 7,
    8, 1, 7, 0),
  4,
  4,
  byrow = TRUE,
  dimnames = list(ids, ids)
)
supplied_distance <- distance
distance <- .correlation_distance_for_hclust(distance)
stopifnot(identical(distance, supplied_distance))
clustering <- .correlation_clustering_functions(distance, TRUE, 'ward.D2')

expected <- stats::hclust(stats::as.dist(distance), method = 'ward.D2')
from_display <- stats::hclust(
  stats::as.dist(sqrt(2 * (1 - raw))),
  method = 'ward.D2'
)
row_hclust <- clustering$rows(raw)
# ComplexHeatmap transposes the matrix before invoking column clustering.
column_hclust <- clustering$columns(t(raw))
stopifnot(
  identical(row_hclust$method, 'ward.D2'),
  identical(row_hclust$merge, expected$merge),
  isTRUE(all.equal(row_hclust$height, expected$height)),
  identical(row_hclust$order, expected$order),
  identical(column_hclust$order, row_hclust$order),
  !identical(row_hclust$merge, from_display$merge)
)

# Split callbacks also reuse an identical tree derived only from the supplied
# dissimilarity subset, rather than recomputing a metric from display values.
slice_ids <- c('S1', 'S3', 'S4')
expected_slice <- stats::hclust(
  stats::as.dist(distance[slice_ids, slice_ids, drop = FALSE]),
  method = 'ward.D2'
)
row_slice <- clustering$rows(raw[slice_ids, , drop = FALSE])
column_slice <- clustering$columns(t(raw[, slice_ids, drop = FALSE]))
stopifnot(
  identical(row_slice$merge, expected_slice$merge),
  identical(row_slice$order, expected_slice$order),
  identical(column_slice$order, row_slice$order)
)
stopifnot(identical(.correlation_text_color(c('#ffffff', '#000000')), c('black', 'white')))
stopifnot(
  identical(.correlation_resolve_linkage('l1', 'auto'), 'average'),
  identical(.correlation_resolve_linkage('l2', 'auto'), 'ward.D2'),
  identical(.correlation_resolve_linkage('spearman', 'auto'), 'ward.D2')
)
stopifnot(identical(.correlation_hclust(distance, ids, 'weighted')$method, 'mcquitty'))

counts <- matrix(100L, 4, 4, dimnames = list(ids, ids))
metadata <- data.frame(
  sample = ids,
  group = c('A', 'A', 'B', 'B'),
  batch = c('x', 'y', 'x', 'y'),
  check.names = FALSE
)
written <- correlation_heatmap(
  raw,
  distance_matrix = sqrt(2 * (1 - raw)),
  overlap_counts = counts,
  metadata = metadata,
  metric = 'pearson',
  linkage = 'ward.D2',
  cut_by = 'group',
  cluster = TRUE,
  annotate = TRUE,
  outname = {json.dumps(str(outname))},
  outfiletypes = '.png',
  fig_width = 7,
  fig_height = 7,
  png_res = 72,
  title = 'correlation integration smoke'
)
stopifnot(length(written) == 1L, file.exists(written[[1]]))
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
