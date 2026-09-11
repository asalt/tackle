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

long_legend_plot <- ggplot(
  data.frame(
    x = seq_len(6),
    y = seq_len(6),
    group = c(
      'noTPO_noTPO_VCF_UBC',
      'PBS_PBS_C57WT',
      'PBS_PBS_hJAK2WT_mJak2KO',
      'TPO_TPO400_C57WT',
      'TPO_TPO400_hJAK2WT_mJak2KO',
      'TPO_TPO50_VCF_UBC'
    )
  ),
  aes(x, y, colour = group)
) +
  geom_point() +
  coord_fixed() +
  theme_classic(base_size = 20)
long_legend_prepared <- pca_prepare_plot_for_output(
  long_legend_plot,
  6,
  7,
  expand_width = TRUE
)
stopifnot(
  long_legend_prepared$legend_width > 4,
  long_legend_prepared$fixed_width > 5,
  long_legend_prepared$fig_width >= 9
)

count_plot <- ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point() +
  coord_fixed(clip = 'off') +
  scale_x_continuous(sec.axis = sec_axis(~.)) +
  scale_y_continuous(sec.axis = sec_axis(~.)) +
  theme_classic(base_size = 20)
attr(count_plot, 'tackle_pca_feature_note') <- list(
  labels = as.expression(list(bquote(bold(n == 1133L)))),
  color = 'grey35',
  fontsize = 13
)
count_prepared <- pca_prepare_plot_for_output(
  count_plot,
  6.25,
  7,
  expand_height = TRUE,
  expand_width = TRUE
)
stopifnot(
  inherits(count_prepared$plot, 'tackle_pca_gtable'),
  count_prepared$feature_note_line_count == 1L,
  count_prepared$feature_note_height > 0,
  count_prepared$fig_height > 7
)
count_layout <- count_prepared$plot$layout
note_layout <- count_layout[count_layout$name == 'pca-feature-note', , drop = FALSE]
xlab_layout <- count_layout[grepl('^xlab-b', count_layout$name), , drop = FALSE]
panel_layout <- count_layout[grepl('^panel', count_layout$name), , drop = FALSE]
stopifnot(
  nrow(note_layout) == 1L,
  note_layout$t > max(xlab_layout$b),
  note_layout$l == min(panel_layout$l),
  note_layout$r == max(panel_layout$r)
)

filtered_plot <- count_plot
attr(filtered_plot, 'tackle_pca_feature_note') <- list(
  labels = as.expression(list(
    bquote(bold(n == 567L)),
    bquote(bold(Var >= P[50]))
  )),
  color = 'grey35',
  fontsize = 13
)
filtered_prepared <- pca_prepare_plot_for_output(
  filtered_plot,
  6.25,
  7,
  expand_height = TRUE,
  expand_width = TRUE
)
stopifnot(
  filtered_prepared$feature_note_line_count == 2L,
  filtered_prepared$feature_note_height > count_prepared$feature_note_height,
  filtered_prepared$fig_height > count_prepared$fig_height
)

png_path <- tempfile(fileext = '.png')
print_png_path <- tempfile(fileext = '.png')
pdf_path <- tempfile(fileext = '.pdf')
ggsave(
  png_path,
  count_prepared$plot,
  width = count_prepared$fig_width,
  height = count_prepared$fig_height,
  dpi = 100
)
ggsave(
  pdf_path,
  filtered_prepared$plot,
  width = filtered_prepared$fig_width,
  height = filtered_prepared$fig_height
)
png(
  print_png_path,
  width = count_prepared$fig_width,
  height = count_prepared$fig_height,
  units = 'in',
  res = 100
)
print(count_prepared$plot)
dev.off()
stopifnot(
  file.info(png_path)$size > 0,
  file.info(print_png_path)$size > 0,
  file.info(pdf_path)$size > 0
)
"""
    result = subprocess.run(
        [rscript, "-e", code],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr or result.stdout
