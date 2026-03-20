suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggrepel))
# library(PCAtools)

## library(cluster)

## ; variable <color> <shape> gene1 gene2 gene3 ...
pca2 <- function(data, outname = "pca", outfiletypes = c(".pdf"),
                 color = NULL, shape = NULL, label = FALSE,
                 encircle = FALSE,
                 showframe = TRUE, frame.type = "t",
                 max_pc = 2, color_list = NULL, marker_list = NULL,
                 names_from = "GeneID",
                 fillna = c("min", "avg"),
                 scale = FALSE,
                 center = TRUE,
                 normalize_by = NULL,
                 show_loadings = FALSE,
                 ntop_loadings = 10,
                 title = NULL,
                 annot_str = NULL,
                 fig_width = 9,
                 fig_height = 8,
                 export_tables = TRUE,
                 return_plots = FALSE,
                 ...) {
  fillna <- match.arg(fillna)


  ## exprs_long <- data %>% pivot_longer(c(-GeneSymbol, -GeneID)) %>%
  ##   left_join(col_data, by='name')
  # TODO: fix bug when `color` or `shape` is set to an r function (e.g. `repeat`)


  # rownames(data) <- data[["GeneID"]]
  # .scalefunc = FALSE
  # if (scale == TRUE){
  #  .scalefunc <- ~ apply(., 2, sd, na.rm = TRUE)
  # }

  # forpca <- data %>%
  #   pivot_wider(id_cols = c(variable, !!normalize_by, !!color, !!shape), names_from = !!names_from, values_from = value)


  #   fillna_func <- switch(fillna,
  #     avg = function(x) mean(x, na.rm = TRUE),
  #     min = function(x) min(x, na.rm = TRUE),
  #     stop("Unsupported fillna type")
  #   )
  #
  #   if (is.null(normalize_by)) {
  #     forpca <- data %>%
  #       group_by(GeneID) %>%
  #       mutate(value = if_else(is.na(value), fillna_func(value), value)) %>%
  #       mutate(value = scale(value, scale = !!scale, center = !!center)) %>%
  #       ungroup() %>%
  #       mutate(value = if_else(is.na(value), 0, value))
  #
  #   } else {
  #     forpca <- data %>%
  #       group_by(GeneID, !!sym(normalize_by)) %>%
  #       mutate(value = if_else(is.na(value), fillna_func(value), value)) %>%
  #       mutate(value = scale(value, scale = !!scale, center = !!center)) %>%
  #       ungroup() %>%
  #       mutate(value = if_else(is.na(value), 0, value))
  #   }

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
    forpca <- data %>%
      group_by(GeneID) %>%
      mutate(fill_value = fillna_func(value)) %>%
      mutate(value = if_else(is.na(value), fill_value, value)) %>%
      # select(-fill_value) %>%
      mutate(value = safe_scale(value, scale = !!scale, center = !!center)) %>%
      ungroup() %>%
      mutate(value = if_else(is.na(value), 0, value))
  } else {
    forpca <- data %>%
      group_by(GeneID, !!sym(normalize_by)) %>%
      mutate(fill_value = fillna_func(value)) %>%
      mutate(value = if_else(is.na(value), fill_value, value)) %>%
      # select(-fill_value) %>%
      mutate(value = safe_scale(value, scale = !!scale, center = !!center)) %>%
      ungroup() %>%
      mutate(value = if_else(is.na(value), 0, value))
  }
  print("Scaling done")

  feature_map <- NULL
  if ("GeneSymbol" %in% colnames(forpca)) {
    feature_map <- forpca %>%
      select(GeneID, GeneSymbol) %>%
      distinct() %>%
      mutate(
        GeneID = as.character(.data$GeneID),
        GeneSymbol = as.character(.data$GeneSymbol)
      )
  }

  .forpca <- forpca %>% pivot_wider(id_cols = c(variable, !!color, !!shape), names_from = !!names_from, values_from = value)

  # if (!is.null(color) && any(names(forpca) == color)) {
  if (!is.null(color)) {
    .forpca <- .forpca %>%
      mutate(!!color := factor(!!sym(color), levels = sort(unique(!!sym(color)))), ordered = TRUE)
  }

  # .forpca <- .forpca %>% select(-variable, -!!color, -!!shape)
  # .forpca <- forpca %>% mutate(!!color := factor(.data[[color]], levels = sort(unique(.data[[color]])), ordered = TRUE))

  # browser()
  drop_cols <- c("variable")
  if (!is.null(color) && is.character(color) && color != "") drop_cols <- c(drop_cols, color)
  if (!is.null(shape) && is.character(shape) && shape != "") drop_cols <- c(drop_cols, shape)
  drop_cols <- unique(drop_cols)

  pca_mat <- .forpca %>% select(-any_of(drop_cols)) %>% as.matrix()
  rownames(pca_mat) <- as.character(.forpca$variable)

  pca_res <- prcomp(pca_mat, scale. = FALSE, center = FALSE)
  variance <- pca_res$sdev^2
  variance_ratio <- variance / sum(variance)
  print("SVD done")

  loading_scores <- as.data.frame(pca_res$rotation)
  loading_scores$feature <- rownames(loading_scores)

  if (!is.null(export_tables) && export_tables == TRUE &&
      !is.null(outname) && !is.na(outname) && outname != "") {
    meta_cols <- c("variable")
    if (!is.null(color) && is.character(color) && color != "") meta_cols <- c(meta_cols, color)
    if (!is.null(shape) && is.character(shape) && shape != "") meta_cols <- c(meta_cols, shape)
    if (!is.null(normalize_by) && is.character(normalize_by) && normalize_by != "") meta_cols <- c(meta_cols, normalize_by)
    meta_cols <- unique(meta_cols)

    meta_df <- forpca %>% select(any_of(meta_cols)) %>% distinct()

    scores_df <- as.data.frame(pca_res$x) %>%
      tibble::rownames_to_column("variable") %>%
      left_join(meta_df, by = "variable") %>%
      select(any_of(meta_cols), everything())

    variance_df <- data.frame(
      pc = seq_along(variance),
      stdev = pca_res$sdev,
      variance = variance,
      variance_ratio = variance_ratio,
      variance_ratio_cum = cumsum(variance_ratio)
    )

    scores_out <- paste0(outname, "_scores.tsv")
    variance_out <- paste0(outname, "_variance.tsv")

    readr::write_tsv(scores_df, scores_out)
    readr::write_tsv(variance_df, variance_out)
  }

  scree_n <- min(length(variance_ratio), max(10, max_pc))
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
    geom_col(
      aes(y = explained_pct),
      fill = "#0b7a75",
      width = 0.72
    ) +
    geom_line(
      aes(y = cumulative_pct, group = 1),
      color = "#d67c2c",
      linewidth = 0.9
    ) +
    geom_point(
      aes(y = cumulative_pct),
      color = "#d67c2c",
      size = 2.2
    ) +
    scale_y_continuous(
      name = "Explained variance (%)",
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.03)),
      sec.axis = sec_axis(~ ., name = "Cumulative variance (%)")
    ) +
    labs(
      title = if (!is.null(title) && !is.na(title) && title != "") {
        paste(title, "Scree Plot")
      } else {
        "PCA Scree Plot"
      },
      subtitle = paste("First", scree_n, "principal components"),
      x = NULL
    ) +
    ggplot2::theme_classic(base_size = 20) +
    theme(
      axis.line.x = element_line(linewidth = 1),
      axis.line.y = element_line(linewidth = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.8)),
      axis.text.x = element_text(angle = 0, vjust = 0.8)
    )

  #   .x <- forpca %>%
  # 	  group_by(!!normalize_by) %>%
  # 	  select(-variable, -!!color, -!!shape)

  # .x %>% scale(center = !!center, scale=.scalefunc)

  ## mutate(!!color = as.character(color), !!shape = as.character(!!shape))
  label_column <- NULL
  label_repel <- FALSE
  if (label) {
    label_column <- "variable"
    label_repel <- TRUE
  }
  #   if (!is.null(normalize_by))
  # 	  {
  #        	forpca %<>% scale(x, center = !!center, scale = apply(x, 2, sd, na.rm = TRUE))
  # 		#forpca %<>% group_by(!!normalize_by) %>% rowwise(scale( center = !!center, scale = !!scale) )
  # 	  }


  ## pca_res <- prcomp(forpca%>%select(-name, -model, -subtype), scale. = FALSE, center=TRUE)
  # pca_res <- prcomp(forpca %>% select(-variable, -!!color, -!!shape), scale. = !!scale, center = !!center)


  # PCAtools::pca(pca_res, )
  # pca(forpca%>%select(-variable, -!!color))
  # pca_res <- PCAtools::pca(forpca %>% select(-variable, -!!color, -!!shape),
  #   metadata = forpca,
  #   scale = scale, center = center,
  # )


  pc_combos <- combn(1:max_pc, 2)
  # browser()

  do_encircle <- isTRUE(encircle) &&
    !is.null(color) && is.character(color) && color != ""
  if (isTRUE(encircle) && !do_encircle) {
    warning("encircle=TRUE requested but `color` is NULL/empty; skipping encircle.")
  }
  if (do_encircle) {
    suppressPackageStartupMessages(suppressMessages({
      if (!requireNamespace("ggalt", quietly = TRUE)) {
        stop("R package 'ggalt' is required for encircle. Install with install.packages('ggalt').")
      }
      library(ggalt)
    }))
  }

  plots <- list()
  for (i in 1:dim(pc_combos)[2]) {
    x1 <- pc_combos[[1, i]]
    x2 <- pc_combos[[2, i]]

    # Ensure a non-empty color palette is available for manual scales.
    # `color_list` may arrive as a (named) list via rpy2; unlist to an atomic vector.
    .color_list <- color_list
    if (!is.null(.color_list) && length(.color_list) > 0) {
      .color_list <- unlist(.color_list)
    } else {
      .color_list <- NULL
    }

    # .....
    # TODO redo
    # .color <- expr(color)
    # color <- "repeat"
      p <- autoplot(pca_res,
        data = .forpca,
        colour = color,
        shape = shape,
        label = label,
      label.repel = label_repel,
      label.label = label_column,
      frame = showframe,
      ## frame.colour = color,
      ## frame.type = 'convex',
      frame.alpha = .5,
      scale_color_manual = .color_list,
      label.size = 2.5,
      size = 4,
      x = x1, y = x2,
      max.overlaps = Inf,
      title = title, # why doesn't this work?
    ) +
      labs(title = title) +
      ggplot2::theme_classic(base_size = 20) +
      coord_fixed(ratio = 1) +
      guides(color = guide_legend(override.aes = list(shape = 15, linetype = "solid"))) +
      theme(
        axis.line.x = element_line(linewidth = 1),
        ## axis.line.x.bottom = element_line(size=1),
        axis.line.y = element_line(linewidth = 1),
        ## axis.line.y.right = element_line(size=1),
      ) +
      scale_x_continuous(sec.axis = sec_axis(~.)) + # hack: https://stackoverflow.com/questions/63055640/how-to-keep-top-and-bottom-axes-in-ggplot-while-removing-the-panel-border
      scale_y_continuous(sec.axis = sec_axis(~.)) +
      scale_shape_manual(values = c(16, 17, 15, 7, 9, 12, 13, 14)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.7, show.legend = NA) +
      geom_vline(xintercept = 0, color = "grey50", linewidth = 0.7)
    if (!is.null(.color_list)) {
      p <- p +
        ggplot2::scale_color_manual(values = .color_list) +
        ggplot2::scale_fill_manual(values = .color_list)
    }
    if (do_encircle && !is.null(p$data) && is.data.frame(p$data) && color %in% colnames(p$data)) {
      encircle_data <- p$data %>%
        filter(!is.na(.data[[color]]))
      encircle_layer <- ggalt::geom_encircle(
        data = encircle_data,
        aes(
          group = .data[[color]],
          fill = .data[[color]],
          colour = .data[[color]]
        ),
        alpha = 0.25,
        size = 0.25,
        show.legend = FALSE,
        na.rm = TRUE
      )
      # Put encircles behind points/labels.
      p$layers <- c(list(encircle_layer), p$layers)
    }
    if (isTRUE(show_loadings) && is.numeric(ntop_loadings) && ntop_loadings > 0) {
      pc_x <- paste0("PC", x1)
      pc_y <- paste0("PC", x2)
      if (all(c(pc_x, pc_y) %in% colnames(loading_scores))) {
        loadings_df <- loading_scores %>%
          transmute(
            feature = .data$feature,
            loading_x = .data[[pc_x]],
            loading_y = .data[[pc_y]],
            loading_score = abs(loading_x) + abs(loading_y)
          ) %>%
          arrange(desc(loading_score)) %>%
          slice_head(n = as.integer(ntop_loadings))

        if (!is.null(feature_map)) {
          loadings_df <- loadings_df %>%
            mutate(feature = as.character(.data$feature)) %>%
            left_join(feature_map, by = c("feature" = "GeneID")) %>%
            mutate(
              feature_label = if_else(
                !is.na(.data$GeneSymbol) & .data$GeneSymbol != "",
                .data$GeneSymbol,
                .data$feature
              )
            )
        } else {
          loadings_df <- loadings_df %>%
            mutate(feature_label = as.character(.data$feature))
        }

        # Scale loading arrows to match ggfortify::autoplot.prcomp default (scale=1),
        # which rescales scores by sdev*sqrt(n_samples). This keeps biplot vectors
        # from distorting the score axis ranges.
        n_samples <- nrow(pca_res$x)
        lam <- pca_res$sdev[c(x1, x2)] * sqrt(n_samples)
        score_x_scaled <- pca_res$x[, pc_x] / lam[1]
        score_y_scaled <- pca_res$x[, pc_y] / lam[2]
        max_score_x <- max(abs(score_x_scaled), na.rm = TRUE)
        max_score_y <- max(abs(score_y_scaled), na.rm = TRUE)

        max_loading_x <- max(abs(loading_scores[[pc_x]]), na.rm = TRUE)
        max_loading_y <- max(abs(loading_scores[[pc_y]]), na.rm = TRUE)
        scaler <- min(max_score_x / max_loading_x, max_score_y / max_loading_y)
        if (!is.finite(scaler) || scaler <= 0) {
          scaler <- 1
        }

        loadings_df <- loadings_df %>%
          mutate(
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
    }
    if (!is.null(annot_str)) {
      p <- p + annotate("text",
        x = -Inf, y = -Inf, hjust = 0, vjust = 0,
        label = annot_str %>% stringr::str_replace_all("_", " ") %>% str_wrap(width = 40)
      )
    }

    if (!is.null(normalize_by)) p <- p + labs(caption = paste("centered within ", normalize_by))


    # THIS doesn't work
    # probably cannot do this (easily) with autoplot
    # if (!is.null(color) & (encircle == TRUE) ) {
    #     p <- p +
    #         ggalt::geom_encircle(
    #                            mapping = aes(#group = !!sym(color),
    #                                          fill = !!sym(color),
    #                                          #colour = !!sym(color)),
    #                            alpha = .2,
    #                            size = 1,
    #                            show.legend = F,
    #                            na.rm = T
    #                            ) +
    #     scale_color_manual(values = color_list) +
    #     scale_fill_manual(values = color_list)
    # }

    key <- paste0(paste0("pc", x1), "_vs_", paste0("pc", x2))

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
      } # end of outfiletypes ext
    }
  } # end of pc combo loop

  if (isTRUE(return_plots)) {
    plots[["scree"]] <- scree_plot
  } else if (!is.null(outname) && !is.na(outname) && outname != "") {
    for (ext in outfiletypes) {
      out <- paste0(outname, "scree", ext)
      print(paste("Saving", out))
      device <- NULL
      if (ext == ".pdf" || ext == "pdf") {
        device <- grDevices::cairo_pdf
      }
      ggplot2::ggsave(out, scree_plot, device = device, width = fig_width, height = fig_height)
    }
  }

  if (isTRUE(return_plots)) {
    return(plots)
  }

  invisible(NULL)
}

## annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1) +
##     annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1)
