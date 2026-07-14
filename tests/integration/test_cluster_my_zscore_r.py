import json
import shutil
import subprocess
from pathlib import Path

import pytest


def test_cluster2_myzscore_anchors_sparse_detections_below_observed_minimum():
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
  identical(constant_complete, c(0, 0, 0))
)
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
