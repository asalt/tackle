suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

umap2 <- function(data, outname = "umap", outfiletypes = c(".pdf"),
                  color = NULL, shape = NULL, label = FALSE,
                  showframe = TRUE,
                  color_list = NULL, marker_list = NULL,
                  names_from = "GeneID",
                  fillna = c("min", "avg"),
                  scale = FALSE,
                  center = TRUE,
                  normalize_by = NULL,
                  title = NULL,
                  annot_str = NULL,
                  fig_width = 9,
                  fig_height = 8,
                  export_tables = TRUE,
                  return_plots = FALSE,
                  n_neighbors = 15,
                  min_dist = 0.1,
                  metric = "euclidean",
                  spread = 1,
                  n_components = 2,
                  max_components = 2,
                  n_epochs = 0,
                  random_state = 1234,
                  ...) {
  fillna <- match.arg(fillna)
  n_components <- as.integer(n_components)
  max_components <- as.integer(max_components)
  n_neighbors <- as.integer(n_neighbors)
  n_epochs <- as.integer(n_epochs)
  random_state <- as.integer(random_state)

  if (n_components < 2) {
    stop("n_components must be >= 2")
  }
  if (max_components < 2) {
    stop("max_components must be >= 2")
  }
  if (n_neighbors < 2) {
    stop("n_neighbors must be >= 2")
  }

  fillna_func <- switch(fillna,
    avg = function(x) {
      vals <- x[is.finite(x)]
      if (length(vals) == 0) return(0)
      mean(vals)
    },
    min = function(x) {
      vals <- x[is.finite(x)]
      if (length(vals) == 0) return(0)
      min(vals)
    },
    stop("Unsupported fillna type")
  )

  safe_scale <- function(x, center = TRUE, scale = TRUE) {
    x <- as.numeric(x)
    if (length(x) == 0) return(x)
    if (all(is.na(x))) return(rep(0, length(x)))

    if (isTRUE(scale)) {
      sdx <- stats::sd(x, na.rm = TRUE)
      if (is.na(sdx) || sdx == 0) {
        if (isTRUE(center)) {
          return(x - mean(x, na.rm = TRUE))
        }
        return(x)
      }
    }

    as.numeric(base::scale(x, center = center, scale = scale))
  }

  if (is.null(normalize_by)) {
    forumap <- data %>%
      group_by(GeneID) %>%
      mutate(fill_value = fillna_func(value)) %>%
      mutate(value = if_else(is.na(value), fill_value, value)) %>%
      mutate(value = safe_scale(value, scale = !!scale, center = !!center)) %>%
      ungroup() %>%
      mutate(value = if_else(is.na(value), 0, value))
  } else {
    forumap <- data %>%
      group_by(GeneID, !!sym(normalize_by)) %>%
      mutate(fill_value = fillna_func(value)) %>%
      mutate(value = if_else(is.na(value), fill_value, value)) %>%
      mutate(value = safe_scale(value, scale = !!scale, center = !!center)) %>%
      ungroup() %>%
      mutate(value = if_else(is.na(value), 0, value))
  }
  .forumap <- forumap %>% pivot_wider(id_cols = c(variable, !!color, !!shape), names_from = !!names_from, values_from = value)

  if (!is.null(color)) {
    .forumap <- .forumap %>%
      mutate(!!color := factor(!!sym(color), levels = sort(unique(!!sym(color)))), ordered = TRUE)
  }

  drop_cols <- c("variable")
  if (!is.null(color) && is.character(color) && color != "") drop_cols <- c(drop_cols, color)
  if (!is.null(shape) && is.character(shape) && shape != "") drop_cols <- c(drop_cols, shape)
  drop_cols <- unique(drop_cols)

  umap_mat <- .forumap %>% select(-any_of(drop_cols)) %>% as.matrix()
  rownames(umap_mat) <- as.character(.forumap$variable)

  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("R package 'uwot' is required for UMAP. Install with install.packages('uwot').")
  }

  set.seed(random_state)
  umap_args <- list(
    X = umap_mat,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    spread = spread,
    n_components = n_components,
    init = "spectral",
    ret_model = FALSE,
    verbose = TRUE
  )
  if (n_epochs > 0) {
    umap_args$n_epochs <- n_epochs
  }
  umap_res <- do.call(uwot::umap, umap_args)
  umap_res <- as.data.frame(umap_res)
  colnames(umap_res) <- paste0("UMAP", seq_len(ncol(umap_res)))
  umap_res <- umap_res %>% tibble::rownames_to_column("variable")

  meta_cols <- c("variable")
  if (!is.null(color) && is.character(color) && color != "") meta_cols <- c(meta_cols, color)
  if (!is.null(shape) && is.character(shape) && shape != "") meta_cols <- c(meta_cols, shape)
  if (!is.null(normalize_by) && is.character(normalize_by) && normalize_by != "") meta_cols <- c(meta_cols, normalize_by)
  meta_cols <- unique(meta_cols)
  meta_df <- forumap %>% select(any_of(meta_cols)) %>% distinct()
  umap_plot_df <- umap_res %>% left_join(meta_df, by = "variable")

  if (!is.null(export_tables) && export_tables == TRUE &&
      !is.null(outname) && !is.na(outname) && outname != "") {
    embedding_out <- paste0(outname, "_embedding.tsv")
    params_out <- paste0(outname, "_params.tsv")

    params_df <- tibble::tibble(
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      metric = metric,
      spread = spread,
      n_components = n_components,
      max_components = max_components,
      n_epochs = n_epochs,
      random_state = random_state,
      center = center,
      scale = scale,
      fillna = fillna,
      normalize_by = ifelse(is.null(normalize_by), "", normalize_by)
    )

    readr::write_tsv(umap_plot_df %>% select(any_of(meta_cols), starts_with("UMAP")), embedding_out)
    readr::write_tsv(params_df, params_out)
  }

  .color_list <- color_list
  if (!is.null(.color_list) && length(.color_list) > 0) {
    .color_list <- unlist(.color_list)
  } else {
    .color_list <- NULL
  }

  max_plot_components <- min(max_components, n_components)
  if (max_plot_components < 2) {
    stop("Need at least 2 components to plot UMAP.")
  }
  dim_combos <- combn(seq_len(max_plot_components), 2)

  plots <- list()
  for (i in seq_len(ncol(dim_combos))) {
    x_idx <- dim_combos[[1, i]]
    y_idx <- dim_combos[[2, i]]
    x_col <- paste0("UMAP", x_idx)
    y_col <- paste0("UMAP", y_idx)

    p <- ggplot(umap_plot_df, aes(x = .data[[x_col]], y = .data[[y_col]]))
    if (!is.null(color) && !is.null(shape)) {
      p <- p + geom_point(
        aes(color = .data[[color]], fill = .data[[color]], shape = .data[[shape]]),
        size = 4, alpha = 0.9, stroke = 0.8
      )
    } else if (!is.null(color)) {
      p <- p + geom_point(
        aes(color = .data[[color]], fill = .data[[color]]),
        size = 4, alpha = 0.9, stroke = 0.8, shape = 21
      )
    } else if (!is.null(shape)) {
      p <- p + geom_point(
        aes(shape = .data[[shape]]),
        size = 4, alpha = 0.9, stroke = 0.8
      )
    } else {
      p <- p + geom_point(size = 4, alpha = 0.9, stroke = 0.8, color = "#2d2d2d")
    }

    if (!is.null(color) && isTRUE(showframe)) {
      p <- p + stat_ellipse(
        aes(group = .data[[color]], fill = .data[[color]]),
        geom = "polygon",
        alpha = 0.14,
        color = NA,
        show.legend = FALSE,
        type = "norm"
      )
    }

    if (isTRUE(label)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = variable),
        size = 3.4,
        max.overlaps = Inf,
        box.padding = 0.22,
        segment.alpha = 0.55
      )
    }

    p <- p +
      labs(
        title = title,
        x = x_col,
        y = y_col
      ) +
      ggplot2::theme_classic(base_size = 20) +
      coord_fixed(ratio = 1) +
      geom_hline(yintercept = 0, color = "grey50", show.legend = NA) +
      geom_vline(xintercept = 0, color = "grey50")

    if (!is.null(shape)) {
      p <- p + scale_shape_manual(values = c(16, 17, 15, 7, 9, 12, 13, 14))
    }
    if (!is.null(.color_list)) {
      p <- p +
        ggplot2::scale_color_manual(values = .color_list) +
        ggplot2::scale_fill_manual(values = .color_list)
    }
    if (!is.null(annot_str)) {
      p <- p + annotate("text",
        x = -Inf, y = -Inf, hjust = 0, vjust = 0,
        label = annot_str %>% stringr::str_replace_all("_", " ") %>% str_wrap(width = 40)
      )
    }
    if (!is.null(normalize_by)) {
      p <- p + labs(caption = paste("centered within ", normalize_by))
    }

    key <- paste0("umap", x_idx, "_vs_", y_idx)
    if (isTRUE(return_plots)) {
      plots[[key]] <- p
    } else if (!is.null(outname) && !is.na(outname) && outname != "") {
      for (ext in outfiletypes) {
        out <- paste0(outname, key, ext)
        print(paste("Saving", out))
        device <- NULL
        if (ext == ".pdf" || ext == "pdf") {
          device <- grDevices::cairo_pdf
        }
        ggplot2::ggsave(out, p, device = device, width = fig_width, height = fig_height)
      }
    }
  }

  if (isTRUE(return_plots)) {
    return(plots)
  }

  invisible(NULL)
}
