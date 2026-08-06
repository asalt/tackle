import json
import shutil
import subprocess
from pathlib import Path

import pytest


def test_cluster2_myzscore_anchors_sparse_detections_below_observed_minimum(tmp_path):
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")

    required = (
        "tidyverse",
        "ComplexHeatmap",
        "circlize",
        "stringr",
        "cluster",
        "dendsort",
    )
    package_check = "; ".join(
        f"stopifnot(requireNamespace('{package}', quietly=TRUE))"
        for package in required
    )
    check = subprocess.run(
        [rscript, "-e", package_check],
        text=True,
        capture_output=True,
        cwd=tmp_path,
    )
    if check.returncode != 0:
        pytest.skip("Required cluster2 R packages are not installed")

    r_file = Path(__file__).resolve().parents[2] / "tackle" / "R" / "clusterplot.R"
    code = f"""
source({json.dumps(str(r_file))})

two_detected <- as.numeric(myzscore(c(NA, 10, 10, NA)))
one_detected <- as.numeric(myzscore(c(NA, NA, 20, NA)))
variable <- as.numeric(myzscore(c(1, NA, 3)))
constant_complete <- as.numeric(myzscore(c(5, 5, 5)))

stopifnot(
  is.na(two_detected[[1]]),
  isTRUE(all.equal(two_detected[2:3], rep(sqrt(3) / 2, 2))),
  is.na(two_detected[[4]]),
  isTRUE(all.equal(one_detected[[3]], 1.5)),
  all(is.na(one_detected[c(1, 2, 4)])),
  variable[[3]] > variable[[1]],
  is.na(variable[[2]]),
  identical(constant_complete, c(0, 0, 0)),
  identical(
    cluster2_legend_title(z_score = "row"),
    "log(iBAQ) zscore"
  ),
  identical(
    cluster2_legend_title(z_score = "row", z_score_by = "group"),
    "log(iBAQ) zscore by group"
  ),
  identical(
    cluster2_legend_title(z_score = "column"),
    "log(iBAQ) sample-wise zscore"
  ),
  identical(
    cluster2_legend_title(z_score = "row", standard_scale = "row"),
    "log(iBAQ) zscore (row-sum scaled)"
  ),
  identical(
    cluster2_legend_title(z_score = "row", standard_scale = "column"),
    "log(iBAQ) zscore (column-sum scaled)"
  ),
  identical(cluster2_legend_title(z_score = NULL), "log(iBAQ)")
)

tiny_data <- data.frame(
  GeneID = c("g1", "g2"),
  GeneSymbol = c("G1", "G2"),
  S1 = c(1, 3),
  S2 = c(2, 4),
  S3 = c(4, 2),
  check.names = FALSE
)
tiny_metadata <- data.frame(
  name = c("S1", "S2", "S3"),
  row.names = c("S1", "S2", "S3")
)
tiny_row_result <- cluster2(
  tiny_data,
  col_data = tiny_metadata,
  z_score = "row",
  row_cluster = FALSE,
  col_cluster = FALSE
)
stopifnot(identical(
  tiny_row_result$heatmap@ht_list[[1]]@matrix_legend_param$title,
  "log(iBAQ) zscore"
))
tiny_row_legend_width <- grid::convertWidth(
  tiny_row_result$heatmap@ht_list[[1]]@matrix_legend_param$legend_width,
  "cm",
  valueOnly = TRUE
)
tiny_row_title_width <- grid::convertWidth(
  1.1 * grid::stringWidth("log(iBAQ) zscore"),
  "cm",
  valueOnly = TRUE
)
stopifnot(
  tiny_row_legend_width >= 3,
  tiny_row_legend_width + 1e-12 >= tiny_row_title_width
)
tiny_row_matrix <- tiny_row_result$heatmap@ht_list[[1]]@matrix
stopifnot(isTRUE(all.equal(
  unname(rowMeans(tiny_row_matrix)),
  rep(0, nrow(tiny_row_matrix)),
  tolerance = 1e-12
)))

tiny_column_result <- cluster2(
  tiny_data,
  col_data = tiny_metadata,
  z_score = "column",
  row_cluster = FALSE,
  col_cluster = FALSE
)
stopifnot(identical(
  tiny_column_result$heatmap@ht_list[[1]]@matrix_legend_param$title,
  "log(iBAQ) sample-wise zscore"
))
tiny_column_legend_width <- grid::convertWidth(
  tiny_column_result$heatmap@ht_list[[1]]@matrix_legend_param$legend_width,
  "cm",
  valueOnly = TRUE
)
tiny_column_title_width <- grid::convertWidth(
  1.1 * grid::stringWidth("log(iBAQ) sample-wise zscore"),
  "cm",
  valueOnly = TRUE
)
stopifnot(
  tiny_column_legend_width >= 3,
  tiny_column_legend_width + 1e-12 >= tiny_column_title_width,
  tiny_column_legend_width > tiny_row_legend_width
)
tiny_column_matrix <- tiny_column_result$heatmap@ht_list[[1]]@matrix
stopifnot(isTRUE(all.equal(
  unname(colMeans(tiny_column_matrix)),
  rep(0, ncol(tiny_column_matrix)),
  tolerance = 1e-12
)))

explicit_missing_data <- data.frame(
  GeneID = "g1",
  GeneSymbol = "G1",
  S1 = 0,
  S2 = 1,
  S3 = NA_real_,
  check.names = FALSE
)
explicit_missing_result <- cluster2(
  explicit_missing_data,
  col_data = tiny_metadata,
  z_score = "row",
  row_cluster = FALSE,
  col_cluster = FALSE
)
explicit_missing_matrix <- explicit_missing_result$heatmap@ht_list[[1]]@matrix
stopifnot(
  is.finite(explicit_missing_matrix[1, "S1"]),
  explicit_missing_matrix[1, "S2"] > explicit_missing_matrix[1, "S1"],
  is.na(explicit_missing_matrix[1, "S3"])
)
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr + result.stdout
