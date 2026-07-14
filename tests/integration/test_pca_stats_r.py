import json
import shutil
import subprocess
from pathlib import Path

import pytest


def test_r_pca_welch_james_and_r2_match_reference_values():
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")

    r_file = Path(__file__).resolve().parents[2] / "tackle" / "R" / "pca_stats.R"
    code = f"""
source({json.dumps(str(r_file))})
Y <- matrix(c(
  0,1, 1,0, 2,2, 1,3, 3,1,
  2,4, 4,2, 3,5, 5,3, 6,6, 4,7,
  6,8, 7,6, 9,9, 8,7, 10,11, 7,10, 11,12
), ncol = 2, byrow = TRUE)
group <- c(rep('A', 5), rep('B', 6), rep('C', 7))
result <- pca_welch_james(Y, group)
stopifnot(
  identical(result$status, 'ok'),
  isTRUE(all.equal(result$statistic, 18.525774053380193)),
  isTRUE(all.equal(result$numerator_df, 4)),
  isTRUE(all.equal(result$denominator_df, 9.711665712288147)),
  isTRUE(all.equal(result$p_value, 0.0001508437775382289))
)

r2_values <- matrix(c(0,0, 2,0, 8,0, 10,0), ncol = 2, byrow = TRUE)
stopifnot(isTRUE(all.equal(
  pca_euclidean_r2(r2_values, c('A', 'A', 'B', 'B')),
  1 - 4 / 68
)))
geometry <- pca_pairwise_centroid_geometry(
  r2_values, c('A', 'A', 'B', 'B')
)
stopifnot(
  geometry$geometry_group_a == 'A',
  geometry$geometry_group_b == 'B',
  isTRUE(all.equal(geometry$centroid_distance, 8)),
  isTRUE(all.equal(geometry$rms_radius_a, 1)),
  isTRUE(all.equal(geometry$rms_radius_b, 1)),
  isTRUE(all.equal(geometry$pooled_rms_radius, 1)),
  isTRUE(all.equal(geometry$standardized_separation, 8))
)

caption_row <- data.frame(
  group_field = 'group', r2 = 0.625, status = 'ok', numerator_df = 4,
  denominator_df = 7.5, welch_james_f = 3.25, p_adj = 0.0123,
  p_adjust_method = 'holm'
)
stopifnot(grepl('Holm-adjusted p=0.0123', pca_format_test_caption(caption_row), fixed = TRUE))
caption_pairwise <- data.frame(
  group_a = 'A', group_b = 'B', centroid_distance = 0.63,
  rms_radius_a = 0.18, rms_radius_b = 0.14,
  pooled_rms_radius = sqrt((0.18^2 + 0.14^2) / 2),
  standardized_separation = 3.85
)
geometry_caption <- pca_format_test_caption(caption_row, caption_pairwise)
stopifnot(
  grepl('Centroid distance = 0.63', geometry_caption, fixed = TRUE),
  grepl('RMS radii (A, B) = 0.18, 0.14', geometry_caption, fixed = TRUE),
  grepl('Standardized separation = 3.85', geometry_caption, fixed = TRUE)
)

set.seed(20260713)
adaptive_scores <- as.data.frame(matrix(rnorm(12 * 12), nrow = 12))
colnames(adaptive_scores) <- paste0('PC', seq_len(12))
rownames(adaptive_scores) <- sprintf('S%02d', 0:11)
adaptive_metadata <- data.frame(
  group = rep(c('A', 'B', 'C', 'D'), each = 3),
  row.names = rownames(adaptive_scores)
)
adaptive_scopes <- list(
  list(
    name = 'PC1_PC2', pcs = c('PC1', 'PC2'),
    plot_key = 'pc1_vs_pc2', selection = 'displayed'
  ),
  list(
    name = 'leading_estimable_pcs', pcs = colnames(adaptive_scores),
    plot_key = NULL, selection = 'leading_estimable'
  )
)
adaptive <- pca_analyze_separation(
  adaptive_scores, adaptive_metadata, 'group', adaptive_scopes, 'holm'
)
expected_variance <- sum(vapply(
  adaptive_scores[c('PC1', 'PC2')], stats::var, numeric(1)
)) / sum(vapply(adaptive_scores, stats::var, numeric(1))) * 100
stopifnot(
  nrow(adaptive$omnibus) == 1,
  adaptive$omnibus$scope[[1]] == 'PC1_PC2',
  !'is_leading_scope' %in% colnames(adaptive$omnibus),
  adaptive$omnibus$n_pcs[[1]] == 2,
  isTRUE(all.equal(adaptive$omnibus$explained_variance_pct[[1]], expected_variance)),
  adaptive$omnibus$status[[1]] == 'ok',
  nrow(adaptive$pairwise) == 6,
  !'is_leading_scope' %in% colnames(adaptive$pairwise),
  all(adaptive$pairwise$n_pcs == 2),
  all(adaptive$pairwise$status == 'ok')
)

two_group_scores <- adaptive_scores[1:6, 1:2, drop = FALSE]
two_group_metadata <- data.frame(
  group = rep(c('A', 'B'), each = 3),
  row.names = rownames(two_group_scores)
)
two_group <- pca_analyze_separation(
  two_group_scores,
  two_group_metadata,
  'group',
  list(list(
    name = 'PC1_PC2', pcs = c('PC1', 'PC2'),
    plot_key = 'pc1_vs_pc2', selection = 'displayed'
  )),
  'holm'
)
stopifnot(
  nrow(two_group$omnibus) == 1,
  !'centroid_distance' %in% colnames(two_group$omnibus),
  nrow(two_group$pairwise) == 1,
  two_group$pairwise$group_a[[1]] == 'A',
  two_group$pairwise$group_b[[1]] == 'B',
  isTRUE(all.equal(two_group$pairwise$p_adj[[1]], two_group$pairwise$p_value[[1]])),
  is.finite(two_group$pairwise$standardized_separation[[1]])
)
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


def test_r_pairwise_separation_plot_is_renderable(tmp_path: Path):
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")

    package_check = subprocess.run(
        [
            rscript,
            "-e",
            "; ".join(
                f"stopifnot(requireNamespace('{package}', quietly=TRUE))"
                for package in ("ggplot2", "ggrepel", "patchwork")
            ),
        ],
        text=True,
        capture_output=True,
    )
    if package_check.returncode != 0:
        pytest.skip("Required PCA separation plotting packages are not installed")

    r_file = Path(__file__).resolve().parents[2] / "tackle" / "R" / "pca_stats.R"
    output = tmp_path / "pairwise-separation.png"
    code = f"""
source({json.dumps(str(r_file))})
pairwise <- data.frame(
  group_field = rep('group', 3),
  scope = rep('PC1_PC2', 3),
  pcs = rep('PC1,PC2', 3),
  n_pcs = rep(2, 3),
  explained_variance_pct = rep(42.5, 3),
  group_a = c('A', 'A', 'B'),
  group_b = c('B', 'C', 'C'),
  centroid_distance = c(8, 5, 3),
  pooled_rms_radius = c(2, 4, 5),
  standardized_separation = c(4, 1.25, 0.6),
  r2 = c(0.7, 0.3, 0.1),
  p_adjust_method = rep('holm', 3),
  p_adj = c(0.003, 0.08, 0.5),
  status = rep('ok', 3)
)
plot <- pca_plot_pairwise_separation(
  pairwise, group_field = 'group', scope = 'PC1_PC2'
)
stopifnot(inherits(plot, 'patchwork'))
ggplot2::ggsave(
  {json.dumps(str(output))}, plot,
  width = 12.5, height = 7.2, dpi = 72, bg = 'white'
)
stopifnot(file.exists({json.dumps(str(output))}))
"""
    result = subprocess.run(
        [rscript, "-e", code],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
    assert output.stat().st_size > 0
