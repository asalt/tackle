suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggfortify))
# library(PCAtools)

## library(cluster)

## ; variable <color> <shape> gene1 gene2 gene3 ...
pca2 <- function(data, outname = "pca", outfiletypes = c(".pdf"),
                 color = NULL, shape = NULL, label = FALSE,
                 encircle = TRUE,
                 showframe = TRUE, frame.type = "t",
                 max_pc = 2, color_list = NULL, marker_list = NULL,
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
    avg = function(x) mean(x, na.rm = TRUE),
    min = function(x) min(x, na.rm = TRUE),
    stop("Unsupported fillna type")
  )

  safe_scale <- function(x, center = TRUE, scale = TRUE) {
    as.numeric(scale(x, center = center, scale = scale))
  }

  if (is.null(normalize_by)) {
    forpca <- data %>%
      group_by(GeneID) %>%
      mutate(fill_value = fillna_func(value)) %>%
      mutate(value = if_else(is.na(value), fill_value, value)) %>%
      # select(-fill_value) %>%
      mutate(value = safe_scale(value, scale = !!scale, center = !!center)) %>%
      ungroup() # %>%
    # mutate(value = if_else(is.na(value), 0, value))
  } else {
    forpca <- data %>%
      group_by(GeneID, !!sym(normalize_by)) %>%
      mutate(fill_value = fillna_func(value)) %>%
      mutate(value = if_else(is.na(value), fill_value, value)) %>%
      # select(-fill_value) %>%
      mutate(value = safe_scale(value, scale = !!scale, center = !!center)) %>%
      ungroup() # %>%
    # mutate(value = if_else(is.na(value), 0, value))
  }
  print("Scaling done")



  # extra guard for nas. shouldn't be necessary now that we fill nas with zeros or minval or avg val
  nas <- forpca[is.na(forpca$value), ]
  names_to_remove <- unique(nas$GeneID)
  forpca <- forpca %>% filter(!GeneID %in% names_to_remove)
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
  print("SVD done")

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

    variance <- pca_res$sdev^2
    variance_ratio <- variance / sum(variance)
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

  plots <- list()
  for (i in 1:dim(pc_combos)[2]) {
    x1 <- pc_combos[[1, i]]
    x2 <- pc_combos[[2, i]]

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
      scale_color_manual = color_list,
      label.size = 2.5,
      size = 4,
      x = x1, y = x2,
      max.overlaps = Inf,
      title = title, # why doesn't this work?
    ) +
      labs(title = title) +
      ggplot2::scale_color_manual(values = color_list) +
      ggplot2::scale_fill_manual(values = color_list) +
      ggplot2::theme_classic(base_size = 20) +
      coord_fixed(ratio = 1) +
      guides(color = guide_legend(override.aes = list(shape = 15, linetype = "solid"))) +
      theme(
        axis.line.x = element_line(size = 1),
        ## axis.line.x.bottom = element_line(size=1),
        axis.line.y = element_line(size = 1),
        ## axis.line.y.right = element_line(size=1),
      ) +
      scale_x_continuous(sec.axis = sec_axis(~.)) + # hack: https://stackoverflow.com/questions/63055640/how-to-keep-top-and-bottom-axes-in-ggplot-while-removing-the-panel-border
      scale_y_continuous(sec.axis = sec_axis(~.)) +
      scale_shape_manual(values = c(16, 17, 15, 7, 9, 12, 13, 14)) +
      geom_hline(yintercept = 0, color = "grey50", show.legend = NA) +
      geom_vline(xintercept = 0, color = "grey50")
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
    return(plots)
  }

  invisible(NULL)
}

## annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1) +
##     annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1)
