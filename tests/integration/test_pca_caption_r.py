import json
import shutil
import subprocess
from pathlib import Path

import pytest


def test_r_pca_caption_wraps_to_measured_plot_width_and_expands_height():
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not available")

    r_file = Path(__file__).resolve().parents[2] / "tackle" / "R" / "pca_caption.R"
    code = f"""
source({json.dumps(str(r_file))})
library(ggplot2)
caption <- paste(
  'treatment: R²=0.297; WJ F*(2, 12.71)=10.79; Holm-adjusted p=0.00547',
  'Centroid distance = 5.65',
  'RMS radii (NO_Young, NO_Old) = 5.43, 2.85',
  'Standardized separation = 1.30',
  sep = '\\n'
)
plot <- ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point() +
  labs(caption = caption) +
  theme_classic(base_size = 20) +
  theme(
    plot.caption.position = 'plot',
    plot.caption = element_text(hjust = 0)
  )
prepared <- pca_prepare_plot_for_output(plot, 6, 7, expand_height = TRUE)
stopifnot(
  prepared$original_line_count == 4,
  prepared$wrapped_line_count > prepared$original_line_count,
  prepared$fig_height > 7,
  prepared$caption_width > 5,
  prepared$caption_width < 6,
  prepared$max_line_width <= prepared$caption_width + 1e-8
)
long_token <- paste(rep('long_group_name', 20), collapse = '_')
long_token_plot <- plot + labs(caption = long_token)
long_token_prepared <- pca_prepare_plot_for_output(
  long_token_plot,
  6,
  7,
  expand_height = TRUE
)
stopifnot(
  long_token_prepared$wrapped_line_count > 1,
  long_token_prepared$max_line_width <=
    long_token_prepared$caption_width + 1e-8
)
"""
    result = subprocess.run(
        [rscript, "-e", code],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr or result.stdout
