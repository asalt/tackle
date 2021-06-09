library(ggplot2)
library(tidyverse)
library(ggfortify)

## library(cluster)

## ; variable <color> <shape> gene1 gene2 gene3 ...
pca2 <- function(data, outname = 'pca', outfiletypes = c('.pdf'),
                 color=NULL, shape=NULL, label=FALSE,
                 showframe = TRUE, frame.type = 't',
                 max_pc = 2, color_list = NULL, marker_list = NULL,
                 names_from = "GeneID",
                 ...
                 ) {

  ## exprs_long <- data %>% pivot_longer(c(-GeneSymbol, -GeneID)) %>%
  ##   left_join(col_data, by='name')
  #TODO: fix bug when `color` or `shape` is set to an r function (e.g. `repeat`)

  forpca <- data %>%
    pivot_wider(id_cols=c(variable, !!color, !!shape), names_from= !!names_from, values_from=value)


    ## mutate(!!color = as.character(color), !!shape = as.character(!!shape))
  label_column <- NULL
  label_repel <- FALSE
  if (label) {
    label_column <- 'variable'
    label_repel <- TRUE
  }


  ## pca_res <- prcomp(forpca%>%select(-name, -model, -subtype), scale. = FALSE, center=TRUE)
  pca_res <- prcomp(forpca %>% select(-variable, -!!color, -!!shape), scale. = FALSE, center=TRUE)


  pc_combos <- combn(1:max_pc, 2)

  for (i in 1:dim(pc_combos)[2]){

    x1 <- pc_combos[[1, i]]
    x2 <- pc_combos[[2, i]]


    for (ext in outfiletypes) {


      #.....
      # TODO redo
      p <- autoplot(pca_res,
                    data = forpca, colour = color, shape = shape,
                    label = label, label.repel = label_repel,
                    label.label = label_column,
                    frame = showframe,
                    ## frame.colour = color,
                    ## frame.type = 'convex',
                    frame.alpha = .5,
                    scale_color_manual=color_list,
                    label.size = 2.5,
                    size = 4,
                    x = x1, y = x2,
                    max.overlaps=Inf
                    ) +
        ggplot2::scale_color_manual(values = color_list) +
        ggplot2::scale_fill_manual(values = color_list) +
        ggplot2::theme_classic(base_size = 20) +
        coord_fixed(ratio = 1) +
        guides(color = guide_legend(override.aes = list(shape=15, linetype = 'solid')))+
        theme(axis.line.x = element_line(size=1),
              ## axis.line.x.bottom = element_line(size=1),
              axis.line.y = element_line(size=1),
              ## axis.line.y.right = element_line(size=1),
        ) +
        scale_x_continuous(sec.axis = sec_axis(~.)) + # hack: https://stackoverflow.com/questions/63055640/how-to-keep-top-and-bottom-axes-in-ggplot-while-removing-the-panel-border
        scale_y_continuous(sec.axis = sec_axis(~.)) +
        geom_hline(yintercept = 0, color = "grey50", show.legend = NA) +
        geom_vline(xintercept = 0, color = "grey50" )

      print(p)


      out <- paste0(outname, paste0('pc', x1), '_vs_', paste0('pc', x2), ext)
      print(paste("Saving", out))
      ggplot2::ggsave(out, p, width = 9, height = 8)
    }
  }

}

        ## annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1) +
        ##     annotate(geom = "segment", y = Inf, yend = Inf, x = -Inf, xend = Inf, size = 1)
