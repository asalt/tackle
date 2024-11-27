library(ggplot2)
library(rlang)
library(purrr)
library(magrittr)
library(tidyverse)
library(ggfortify)
# library(PCAtools)

## library(cluster)

## ; variable <color> <shape> gene1 gene2 gene3 ...
pca2 <- function(data, outname = "pca", outfiletypes = c(".pdf"),
                 color = NULL, shape = NULL, label = FALSE,
                 encircle = TRUE,
                 showframe = TRUE, frame.type = "t",
                 max_pc = 2, color_list = NULL, marker_list = NULL,
                 names_from = "GeneID",
                 scale = FALSE,
                 center = TRUE,
                 normalize_by = NULL,
                 title = NULL,
                 annot_str = NULL,
                 fig_width=9,
                 fig_height=8,
                 ...) {

  ## exprs_long <- data %>% pivot_longer(c(-GeneSymbol, -GeneID)) %>%
  ##   left_join(col_data, by='name')
  # TODO: fix bug when `color` or `shape` is set to an r function (e.g. `repeat`)


  # rownames(data) <- data[["GeneID"]]
  # .scalefunc = FALSE
  # if (scale == TRUE){
  #  .scalefunc <- ~ apply(., 2, sd, na.rm = TRUE)
  # }
  # browser()

  # forpca <- data %>%
  #   pivot_wider(id_cols = c(variable, !!normalize_by, !!color, !!shape), names_from = !!names_from, values_from = value)


  if (is.null(normalize_by)) {
    forpca <- data %>%
      group_by(GeneID) %>%
      mutate(value = scale(value, scale = !!scale, center = !!center)) %>%
      ungroup()


    # pivot_wider(id_cols = c(variable, !!color, !!shape), names_from = !!names_from, values_from = value)
    # select(-variable, -!!color, -!!shape) %>% scale( scale = !!scale, center = !!center)

    # .forpca <- forpca %>% pivot_wider(id_cols = c(variable, !!color, !!shape), names_from = !!names_from, values_from = value)
    # .forpca <- forpca %>% pivot_wider(id_cols=variable, names_from = !!names_from, values_from = value) %>% select(-variable)
    # .forpca <- forpca %>% select(-variable, -!!color, -!!shape) %>% scale( scale = !!scale, center = !!center)
  } else {
    forpca <- data %>%
      group_by(GeneID, !!sym(normalize_by)) %>%
      mutate(value = scale(value, scale = !!scale, center = !!center)) %>%
      ungroup()
    # .forpca <- forpca %>% pivot_wider(names_from = !!names_from, values_from = value) %>% select(-variable)
    # .forpca <- forpca %>% select(-variable, -!!color, -!!shape) %>% scale( scale = !!scale, center = !!center)
    # forpca <- forpca %>% group_by(normalize_by) %>% select(-variable, -!!color, -!!shape) %>% scale( scale = .scalefunc, center = !!center)
    # .forpca <-
  }
  nas <- forpca[ is.na(forpca$value), ]
  names_to_remove <- unique(nas$GeneID)
  forpca %<>% filter(!GeneID %in% names_to_remove)
  forpca %<>% pivot_wider(id_cols = c(variable, !!color, !!shape), names_from = !!names_from, values_from = value)

  if (!is.null(color) && any(names(forpca) == color)) {
      forpca <- forpca %>%
        mutate(!!color := as.factor(!!sym(color)))
  }

  .forpca <- forpca %>% select(-variable, -!!color, -!!shape)
  pca_res <- prcomp(.forpca, scale. = F, center = F)

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

  for (i in 1:dim(pc_combos)[2]) {
    x1 <- pc_combos[[1, i]]
    x2 <- pc_combos[[2, i]]


    for (ext in outfiletypes) {
      # .....
      # TODO redo
      # .color <- expr(color)
      # color <- "repeat"
      p <- autoplot(pca_res,
        data = forpca,
        colour = color,
        shape = shape,
        label = label, label.repel = label_repel,
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
      if (!is.null(annot_str)) p <- p + annotate("text", x = -Inf, y = -Inf, hjust = 0, vjust = 0,
                                                 label = annot_str %>% stringr::str_replace_all("_", " ") %>% str_wrap(width = 40))

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

      # print(p)


      out <- paste0(outname, paste0("pc", x1), "_vs_", paste0("pc", x2), ext)
      print(paste("Saving", out))
      device <- NULL
      if (ext == ".pdf" || ext == "pdf") {
        device <- grDevices::cairo_pdf
      }
      ggplot2::ggsave(out, p, device = device, width = fig_width, height = fig_height)
    }
  }
}

## annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1) +
##     annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1)
