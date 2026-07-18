from __future__ import annotations

import json
import textwrap


def render_pca_replay_rmd(
    *,
    title: str,
    gct_relpath: str = "pca_input_pre_svd.gctx",
    context_relpath: str = "pca_replay_context.json",
    include_separation: bool = False,
) -> str:
    """Render a standalone pca2 replay with tackle's score/biplot styling."""

    template = r'''
---
title: __TITLE__
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: show
---

```{r setup, include=FALSE}

required_packages <- c(
  "cmapR",
  "jsonlite",
  "ggplot2",
  "ggfortify",
  "ggrepel",
  "dplyr",
  "readr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Required packages are unavailable: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggfortify)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(tibble)
  library(jsonlite)
  library(cmapR)
})


```

## Replay contract

The GCT matrix is the authoritative input. Its rows are samples and its columns are
features, exactly matching the matrix passed to `prcomp()` by tackle. No missing-value
filling, centering, scaling, normalization, or transposition is performed in this replay.

## Editable replay parameters

This block gathers the captured plotting parameters in one place. Edit `pc_pairs` or
`replot_formats` here before rerunning the PCA and plotting chunks below.

```{r replay-parameters}
ctx <- jsonlite::fromJSON("__CONTEXT__", simplifyVector = FALSE)
plot_params <- ctx$plot_parameters
captured_figsize <- as.numeric(unlist(plot_params$figsize, use.names = FALSE))
fig_width <- if (length(captured_figsize) >= 1 && is.finite(captured_figsize[[1]])) captured_figsize[[1]] else 6
fig_height <- if (length(captured_figsize) >= 2 && is.finite(captured_figsize[[2]])) captured_figsize[[2]] else 7
pc_pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
replot_formats <- as.character(unlist(plot_params$file_formats, use.names = FALSE))
if (length(replot_formats) == 0) replot_formats <- ".png"
out_dir <- "pca_replay_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
knitr::opts_chunk$set(fig.width = fig_width, fig.height = fig_height)
```

## Load the authoritative input

```{r load-authoritative-input}
parse_gct_any <- function(path) {
  exports <- getNamespaceExports("cmapR")
  if ("parse_gctx" %in% exports) return(cmapR::parse_gctx(path))
  if ("parse.gctx" %in% exports) return(cmapR::parse.gctx(path))
  stop("cmapR does not export parse.gctx / parse_gctx")
}

stored_ds <- parse_gct_any("__GCT__")
pca_mat <- as.matrix(stored_ds@mat)
storage.mode(pca_mat) <- "double"

expected_samples <- as.integer(unlist(ctx$matrix$n_samples))
expected_features <- as.integer(unlist(ctx$matrix$n_features))
stopifnot(nrow(pca_mat) == expected_samples)
stopifnot(ncol(pca_mat) == expected_features)
stopifnot(identical(rownames(pca_mat), as.character(unlist(ctx$sample_ids))))
stopifnot(identical(colnames(pca_mat), as.character(unlist(ctx$feature_ids))))
if (any(!is.finite(pca_mat))) {
  stop("The authoritative pre-SVD matrix unexpectedly contains non-finite values.")
}

sample_meta <- as.data.frame(stored_ds@rdesc, check.names = FALSE, stringsAsFactors = FALSE)
rownames(sample_meta) <- as.character(stored_ds@rid)
sample_meta <- sample_meta[rownames(pca_mat), , drop = FALSE]
sample_meta$variable <- rownames(sample_meta)

feature_meta <- as.data.frame(stored_ds@cdesc, check.names = FALSE, stringsAsFactors = FALSE)
rownames(feature_meta) <- as.character(stored_ds@cid)
feature_meta <- feature_meta[colnames(pca_mat), , drop = FALSE]
```

## PCA

These are deliberately the fixed arguments used at tackle's final PCA boundary. The
requested preprocessing settings are recorded in the context and have already been
realized in `pca_mat`.

```{r run-pca}
pca_res <- stats::prcomp(pca_mat, center = FALSE, scale. = FALSE)
variance <- pca_res$sdev^2
variance_ratio <- variance / sum(variance)

scores_df <- as.data.frame(pca_res$x) |>
  tibble::rownames_to_column("variable") |>
  dplyr::left_join(
    sample_meta |> tibble::rownames_to_column("sample_row") |> dplyr::select(-variable),
    by = c("variable" = "sample_row")
  )
variance_df <- data.frame(
  pc = seq_along(variance),
  stdev = pca_res$sdev,
  variance = variance,
  variance_ratio = variance_ratio,
  variance_ratio_cum = cumsum(variance_ratio)
)
readr::write_tsv(scores_df, file.path(out_dir, "pca_scores.tsv"))
readr::write_tsv(variance_df, file.path(out_dir, "pca_variance.tsv"))
knitr::kable(head(scores_df, 10))
```

__SEPARATION_START__
## Heteroscedastic PCA separation tests

When separation testing was requested, these tables recompute the descriptive
Euclidean R2 and Johansen Welch-James approximate-df tests from the exact PCA
score columns. Pairwise rows are produced for every factor having at least two
levels, including the sole comparison for a two-level factor. This is an
independent-samples test; no pairing or blocking is inferred.
The default `leading` scope is reduced to the largest PC1 through PCk block for
which both the omnibus test and every required pairwise test are estimable. The
`explained_variance_pct` column is the percentage of total PCA variance carried
by the exact components in that row. For every two-group comparison, centroid
distance and the two within-group RMS radii are reported along with pooled RMS
radius and standardized separation (centroid distance divided by pooled RMS
radius). Explicitly requested single-PC tests use ordinary Welch one-way ANOVA
and two-sided pairwise Welch t-tests without variance pooling, moderation, or
regularization.

```{r pca-separation-tests, results='asis'}
separation_config <- ctx$separation_testing
separation_enabled <- !is.null(separation_config$enabled) &&
  length(unlist(separation_config$enabled, use.names = FALSE)) > 0 &&
  isTRUE(as.logical(unlist(separation_config$enabled, use.names = FALSE)[[1]]))
separation_omnibus <- data.frame()
separation_pairwise <- data.frame()
single_pc_omnibus <- data.frame()
single_pc_pairwise <- data.frame()
separation_caption_by_plot <- character()
separation_summary_plots <- list()

if (separation_enabled) {
  source("pca_stats.R", local = TRUE)
  adjustment <- as.character(
    unlist(separation_config$p_adjust_method, use.names = FALSE)[[1]]
  )
  separation_result <- pca_analyze_separation(
    scores = as.data.frame(pca_res$x, check.names = FALSE),
    metadata = sample_meta,
    group_fields = separation_config$group_fields,
    scopes = separation_config$resolved_scopes,
    p_adjust_method = adjustment
  )
  separation_omnibus <- separation_result$omnibus
  separation_pairwise <- separation_result$pairwise
  readr::write_tsv(
    separation_omnibus,
    file.path(out_dir, "pca_separation_omnibus.tsv")
  )
  if (nrow(separation_pairwise) > 0) {
    readr::write_tsv(
      separation_pairwise,
      file.path(out_dir, "pca_separation_pairwise.tsv")
    )
  }

  single_pcs <- as.character(
    unlist(separation_config$resolved_single_pcs, use.names = FALSE)
  )
  if (length(single_pcs) > 0) {
    single_pc_result <- pca_analyze_single_pc_separation(
      scores = as.data.frame(pca_res$x, check.names = FALSE),
      metadata = sample_meta,
      group_fields = separation_config$group_fields,
      pcs = single_pcs,
      p_adjust_method = adjustment
    )
    single_pc_omnibus <- single_pc_result$omnibus
    single_pc_pairwise <- single_pc_result$pairwise
    readr::write_tsv(
      single_pc_omnibus,
      file.path(out_dir, "pca_single_pc_omnibus.tsv")
    )
    readr::write_tsv(
      single_pc_pairwise,
      file.path(out_dir, "pca_single_pc_pairwise.tsv")
    )
  }

  for (scope in separation_config$resolved_scopes) {
    plot_key <- as.character(unlist(scope$plot_key, use.names = FALSE))
    scope_name <- as.character(unlist(scope$name, use.names = FALSE)[[1]])
    if (length(plot_key) == 1 && !is.na(plot_key) && nzchar(plot_key)) {
      rows <- separation_omnibus[separation_omnibus$scope == scope_name, , drop = FALSE]
      pair_rows <- separation_pairwise[
        separation_pairwise$scope == scope_name,
        ,
        drop = FALSE
      ]
      separation_caption_by_plot[[plot_key]] <- pca_format_test_caption(
        rows,
        pair_rows
      )
    }
  }

  cat("### Omnibus factor-level tests\n\n")
  print(knitr::kable(separation_omnibus, digits = 4))
  if (nrow(separation_pairwise) > 0) {
    cat("\n\n### Pairwise comparisons\n\n")
    print(knitr::kable(separation_pairwise, digits = 4))
  }
  if (nrow(single_pc_omnibus) > 0) {
    cat("\n\n### Single-PC Welch ANOVA tests\n\n")
    print(knitr::kable(single_pc_omnibus, digits = 4))
    cat("\n\n### Single-PC pairwise Welch t-tests\n\n")
    print(knitr::kable(single_pc_pairwise, digits = 4))
  }

  if (nrow(separation_pairwise) > 0 || nrow(single_pc_pairwise) > 0) {
    cat("\n\n### Pairwise separation summaries\n\n")
    safe_name <- function(value) gsub("[^A-Za-z0-9._-]+", "_", value)
    summary_formats <- as.character(
      unlist(ctx$plot_parameters$file_formats, use.names = FALSE)
    )
    if (length(summary_formats) == 0) summary_formats <- ".png"
    write_summary_plots <- function(pairwise_table, allowed_scopes) {
      if (nrow(pairwise_table) == 0) return(invisible(NULL))
      for (field in unique(pairwise_table$group_field)) {
        for (scope_name in intersect(unique(pairwise_table$scope), allowed_scopes)) {
          pair_rows <- pairwise_table[
            pairwise_table$group_field == field &
              pairwise_table$scope == scope_name,
            ,
            drop = FALSE
          ]
          finite_geometry <- is.finite(pair_rows$centroid_distance) &
            is.finite(pair_rows$pooled_rms_radius) &
            is.finite(pair_rows$standardized_separation)
          pair_rows <- pair_rows[finite_geometry, , drop = FALSE]
          if (nrow(pair_rows) == 0) next
          plot_name <- paste(
            "separation",
            safe_name(field),
            safe_name(scope_name),
            sep = "_"
          )
          summary_plot <- pca_plot_pairwise_separation(
            pair_rows,
            group_field = field,
            scope = scope_name
          )
          separation_summary_plots[[plot_name]] <<- summary_plot
          print(summary_plot)
          for (ext in summary_formats) {
            ggplot2::ggsave(
              file.path(out_dir, paste0("pca_", plot_name, ext)),
              summary_plot,
              width = 12.5,
              height = 7.2,
              dpi = 320,
              bg = "white"
            )
          }
        }
      }
    }
    displayed_scope_names <- vapply(
      Filter(
        function(scope) {
          key <- as.character(unlist(scope$plot_key, use.names = FALSE))
          length(key) == 1 && !is.na(key) && nzchar(key)
        },
        separation_config$resolved_scopes
      ),
      function(scope) as.character(unlist(scope$name, use.names = FALSE)[[1]]),
      character(1)
    )
    write_summary_plots(separation_pairwise, displayed_scope_names)
    write_summary_plots(single_pc_pairwise, single_pcs)
  }
}
```
__SEPARATION_END__

## Tackle-style replot code

The editable defaults below draw PC1/PC2, PC1/PC3, and PC2/PC3 whenever those
components exist. Plot annotations and loading arrows follow the parameters captured
from the original `pca2` invocation.

### Plot settings and metadata

```{r plot-settings}
scalar_or_null <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  value <- as.character(unlist(x, use.names = FALSE)[[1]])
  if (is.na(value) || !nzchar(value) || value == "NULL") return(NULL)
  value
}
flag_value <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  isTRUE(as.logical(unlist(x, use.names = FALSE)[[1]]))
}
integer_value <- function(x, default) {
  if (is.null(x) || length(x) == 0) return(as.integer(default))
  as.integer(unlist(x, use.names = FALSE)[[1]])
}

color_field <- scalar_or_null(plot_params$color)
shape_field <- scalar_or_null(plot_params$marker)
label_points <- flag_value(plot_params$annotate)
show_frame <- flag_value(plot_params$frame)
do_encircle <- flag_value(plot_params$encircle) && !is.null(color_field)
show_loadings <- flag_value(plot_params$show_loadings)
ntop_loadings <- integer_value(plot_params$ntop_loadings, 10)
plot_title <- scalar_or_null(plot_params$title)
normalize_by <- scalar_or_null(ctx$preprocessing$normalize_by)

if (!is.null(color_field) && !color_field %in% colnames(sample_meta)) {
  stop("Captured color field is absent from GCT row metadata: ", color_field)
}
if (!is.null(shape_field) && !shape_field %in% colnames(sample_meta)) {
  stop("Captured marker field is absent from GCT row metadata: ", shape_field)
}
if (!is.null(color_field)) {
  sample_meta[[color_field]] <- factor(
    as.character(sample_meta[[color_field]]),
    levels = sort(unique(as.character(sample_meta[[color_field]]))),
    ordered = TRUE
  )
}

color_values <- NULL
if (!is.null(color_field) && !is.null(ctx$metadata_colors[[color_field]])) {
  color_values <- unlist(ctx$metadata_colors[[color_field]], use.names = TRUE)
}

loading_scores <- as.data.frame(pca_res$rotation)
loading_scores$feature <- rownames(loading_scores)
feature_labels <- data.frame(
  feature = rownames(feature_meta),
  feature_label = rownames(feature_meta),
  stringsAsFactors = FALSE
)
if ("GeneSymbol" %in% colnames(feature_meta)) {
  symbols <- as.character(feature_meta$GeneSymbol)
  keep <- !is.na(symbols) & nzchar(symbols)
  feature_labels$feature_label[keep] <- symbols[keep]
}
```

### PCA plotting function

```{r pca-plot-function}
make_pca_plot <- function(pc_pair) {
  x1 <- pc_pair[[1]]
  x2 <- pc_pair[[2]]
  p <- ggplot2::autoplot(
    pca_res,
    data = sample_meta,
    colour = color_field,
    shape = shape_field,
    label = label_points,
    label.repel = label_points,
    label.label = if (label_points) "variable" else NULL,
    frame = show_frame,
    frame.alpha = 0.5,
    label.size = 2.5,
    size = 4,
    x = x1,
    y = x2,
    max.overlaps = Inf,
    title = plot_title
  ) +
    labs(title = plot_title) +
    ggplot2::theme_classic(base_size = 20) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(override.aes = list(shape = 15, linetype = "solid"))) +
    theme(
      axis.line.x = element_line(linewidth = 1),
      axis.line.y = element_line(linewidth = 1)
    ) +
    scale_x_continuous(sec.axis = sec_axis(~.)) +
    scale_y_continuous(sec.axis = sec_axis(~.)) +
    scale_shape_manual(values = c(16, 17, 15, 7, 9, 12, 13, 14)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.7, show.legend = NA) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.7)

  if (!is.null(color_values)) {
    p <- p +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::scale_fill_manual(values = color_values)
  }

  if (do_encircle) {
    if (!requireNamespace("ggalt", quietly = TRUE)) {
      stop("The replayed pca2 call used --encircle, but R package 'ggalt' is missing.")
    }
    encircle_data <- p$data |>
      dplyr::filter(!is.na(.data[[color_field]]))
    encircle_layer <- ggalt::geom_encircle(
      data = encircle_data,
      aes(
        group = .data[[color_field]],
        fill = .data[[color_field]],
        colour = .data[[color_field]]
      ),
      alpha = 0.25,
      size = 0.25,
      show.legend = FALSE,
      na.rm = TRUE
    )
    p$layers <- c(list(encircle_layer), p$layers)
  }

  if (show_loadings && ntop_loadings > 0) {
    pc_x <- paste0("PC", x1)
    pc_y <- paste0("PC", x2)
    loadings_df <- loading_scores |>
      dplyr::transmute(
        feature = .data$feature,
        loading_x = .data[[pc_x]],
        loading_y = .data[[pc_y]],
        loading_score = abs(loading_x) + abs(loading_y)
      ) |>
      dplyr::arrange(dplyr::desc(loading_score)) |>
      dplyr::slice_head(n = ntop_loadings) |>
      dplyr::left_join(feature_labels, by = "feature")

    n_samples <- nrow(pca_res$x)
    lam <- pca_res$sdev[c(x1, x2)] * sqrt(n_samples)
    score_x_scaled <- pca_res$x[, pc_x] / lam[[1]]
    score_y_scaled <- pca_res$x[, pc_y] / lam[[2]]
    max_score_x <- max(abs(score_x_scaled), na.rm = TRUE)
    max_score_y <- max(abs(score_y_scaled), na.rm = TRUE)
    max_loading_x <- max(abs(loading_scores[[pc_x]]), na.rm = TRUE)
    max_loading_y <- max(abs(loading_scores[[pc_y]]), na.rm = TRUE)
    scaler <- min(max_score_x / max_loading_x, max_score_y / max_loading_y)
    if (!is.finite(scaler) || scaler <= 0) scaler <- 1
    loadings_df <- loadings_df |>
      dplyr::mutate(
        xend = loading_x * scaler * 0.8,
        yend = loading_y * scaler * 0.8
      )

    p <- p +
      geom_segment(
        data = loadings_df,
        aes(x = 0, y = 0, xend = xend, yend = yend),
        inherit.aes = FALSE,
        linewidth = 0.4,
        alpha = 0.7,
        color = "#666666",
        arrow = grid::arrow(length = grid::unit(0.012, "npc"))
      ) +
      ggrepel::geom_text_repel(
        data = loadings_df,
        aes(x = xend, y = yend, label = feature_label),
        inherit.aes = FALSE,
        size = 3.2,
        color = "#4a4a4a",
        max.overlaps = Inf,
        box.padding = 0.25,
        segment.color = "#9a9a9a",
        segment.alpha = 0.55
      )
  }

  if (!is.null(normalize_by)) {
    p <- p + labs(caption = paste("centered within", normalize_by))
  }
__SEPARATION_CAPTION_START__
  plot_key <- paste0("pc", x1, "_vs_pc", x2)
  if (plot_key %in% names(separation_caption_by_plot)) {
    stats_caption <- separation_caption_by_plot[[plot_key]]
    existing_caption <- p$labels$caption
    if (!is.null(existing_caption) && length(existing_caption) > 0 &&
        !is.na(existing_caption[[1]]) && nzchar(existing_caption[[1]])) {
      stats_caption <- paste(existing_caption[[1]], stats_caption, sep = "\n")
    }
    p <- p +
      labs(caption = stats_caption) +
      theme(
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0)
      )
  }
__SEPARATION_CAPTION_END__
  if (!is.null(plot_title)) {
    annotation_label <- paste(strwrap(gsub("_", " ", plot_title), width = 40), collapse = "\n")
    p <- p + annotate(
      "text",
      x = -Inf,
      y = -Inf,
      hjust = 0,
      vjust = 0,
      label = annotation_label
    )
  }
  p
}
```

### Build the requested PCA planes

```{r build-pca-plots}
available_pcs <- ncol(pca_res$x)
pc_pairs <- Filter(function(pair) all(pair <= available_pcs), pc_pairs)
plots <- setNames(
  lapply(pc_pairs, make_pca_plot),
  vapply(pc_pairs, function(pair) paste0("pc", pair[[1]], "_vs_pc", pair[[2]]), character(1))
)
```

## PC1/PC2, PC1/PC3, and PC2/PC3

```{r pca-pairs, results='asis'}
for (plot_name in names(plots)) {
  cat("\n\n### ", toupper(gsub("_", " ", plot_name)), "\n\n", sep = "")
  print(plots[[plot_name]])
  for (ext in replot_formats) {
    if (!startsWith(ext, ".")) ext <- paste0(".", ext)
    ggplot2::ggsave(
      file.path(out_dir, paste0(plot_name, ext)),
      plots[[plot_name]],
      width = fig_width,
      height = fig_height
    )
  }
}
```

## Scree plot

```{r scree}
scree_n <- min(length(variance_ratio), max(10, 3))
scree_idx <- seq_len(scree_n)
scree_df <- tibble::tibble(
  pc_index = scree_idx,
  pc_label = factor(
    paste0("PC", scree_idx),
    levels = paste0("PC", scree_idx),
    ordered = TRUE
  ),
  explained_pct = variance_ratio[scree_idx] * 100,
  cumulative_pct = cumsum(variance_ratio[scree_idx]) * 100
)
scree_plot <- ggplot(scree_df, aes(x = pc_label)) +
  geom_col(aes(y = explained_pct), fill = "#0b7a75", width = 0.72) +
  geom_line(aes(y = cumulative_pct, group = 1), color = "#d67c2c", linewidth = 0.9) +
  geom_point(aes(y = cumulative_pct), color = "#d67c2c", size = 2.2) +
  scale_y_continuous(
    name = "Explained variance (%)",
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.03)),
    sec.axis = sec_axis(~ ., name = "Cumulative variance (%)")
  ) +
  labs(title = "PCA Scree Plot", subtitle = paste("First", scree_n, "principal components"), x = NULL) +
  ggplot2::theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = rel(0.8)))
print(scree_plot)
ggplot2::ggsave(file.path(out_dir, "scree.png"), scree_plot, width = fig_width, height = fig_height)
```

## Captured input parameters

```{r context-table}
parameter_rows <- data.frame(
  parameter = names(ctx$preprocessing),
  value = vapply(ctx$preprocessing, function(x) paste(unlist(x), collapse = ", "), character(1)),
  stringsAsFactors = FALSE
)
knitr::kable(parameter_rows)
```

```{r context-dump, results='asis'}
cat("<details><summary><strong>Full replay context (JSON)</strong></summary><pre>")
cat(jsonlite::toJSON(ctx, pretty = TRUE, auto_unbox = TRUE))
cat("</pre></details>")
```
'''
    rendered = (
        textwrap.dedent(template)
        .lstrip()
        .replace("__TITLE__", json.dumps(str(title)))
        .replace("__GCT__", str(gct_relpath).replace("\\", "/"))
        .replace("__CONTEXT__", str(context_relpath).replace("\\", "/"))
    )
    conditional_blocks = (
        ("__SEPARATION_START__", "__SEPARATION_END__"),
        ("__SEPARATION_CAPTION_START__", "__SEPARATION_CAPTION_END__"),
    )
    for start, end in conditional_blocks:
        if include_separation:
            rendered = rendered.replace(start, "").replace(end, "")
            continue
        before, remainder = rendered.split(start, 1)
        _conditional, after = remainder.split(end, 1)
        rendered = before.rstrip() + "\n" + after.lstrip()
    return rendered


def render_pca_replay_sh(*, rmd_name: str = "replot.Rmd") -> str:
    return textwrap.dedent(
        f'''\
        #!/usr/bin/env bash
        set -euo pipefail
        cd "$(dirname "$0")"
        Rscript -e "rmarkdown::render('{rmd_name}')"
        '''
    )
