library(ggplot2)
library(tidyverse)
library(ggfortify)

## library(cluster)

## ; variable <color> <shape> gene1 gene2 gene3 ...
pca2 <- function(data, outname = 'pca', outfiletypes = c('.pdf'),
                 color=NULL, shape=NULL, label=FALSE,
                 showframe = TRUE, frame.type = 't',
                 max_pc = 2, color_list = NULL, marker_list = NULL
                 ) {

  ## exprs_long <- data %>% pivot_longer(c(-GeneSymbol, -GeneID)) %>%
  ##   left_join(col_data, by='name')

  forpca <- data %>%
    pivot_wider(id_cols=c(variable, !!color, !!shape), names_from=GeneID, values_from=value)


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


      p <- autoplot(pca_res,
                    data = forpca, colour = color, shape = shape,
                    label = label, label.repel = label_repel,
                    label.label = label_column,
                    frame = showframe,
                    ## frame.colour = color,
                    ## frame.type = 'convex',
                    frame.alpha = .5,
                    scale_color_manual=color_list,
                    label.size = 3,
                    size = 4,
                    x = x1, y = x2
                    ) +
        ggplot2::scale_color_manual(values = color_list) +
        ggplot2::scale_fill_manual(values = color_list) +
        ggplot2::theme_classic(base_size = 16)
      print(p)


      out <- paste0(outname, paste0('pc', x1), '_vs_', paste0('pc', x2), ext)
      print(paste("Saving", out))
      ggplot2::ggsave(out, p, width = 7, height = 7)
    }
  }

}
