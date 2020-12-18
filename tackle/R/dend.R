library(dplyr)
library(ggdendro)
library(circlize)
library(dendextend)
library(ComplexHeatmap)



plotdend <- function(data, col_data, color=NULL, shape=NULL, linkage='average'){


  ## print(linkage)
  ## dend <- dist( select(data, -GeneID, -GeneSymbol) %>% t, method = 'euclidean') %>%
  ##   hclust(method = linkage)
  ## col_data <- read_tsv('./tmp_meta.tsv')


  dist_clust <- dist( select(data, -GeneID) %>% t, method = 'euclidean') %>%
    hclust(method = linkage)
  dend <- as.dendrogram(dist_clust)

  ddata_x <- ggdendro::dendro_data(dend)
  ## labs <- label(ddata_x) %>% left_join(dfm, by=c("label" = "rowname"))
  ## browser()
  labs <- ddata_x$label %>% left_join(col_data, by=c("label" = "name"))

  ## dendextend::labels_colors(dend) <- as.factor(labs[['model']]) %>% as.numeric

  # this is for coloring the text
  dendextend::labels_colors(dend) <- as.factor(labs[[color[[1]]]]) %>% as.numeric()

  ## palette(brewer.pal(8, 'Dark2'))     # six color rainbow


  ## pdf(paste0('./results/dendrogram', '.pdf'),
  ##     width=12, height=12)
  ## ## dend1 <- color_branches(dend, h = 3)
  ## circlize_dendrogram(dend, labels_track_height = NA,  dend_track_height = 0.5,)
  ## dev.off()

  i <- 1
  color_numeric_factors = list()
  for (thecolor in color){
    color_numeric_factor <- as.factor(labs[[thecolor]]) %>% as.numeric()
    color_numeric_factors[[i]] <- color_numeric_factor
    i <- i + 1
  }


  ## browser()

  ## ====================================================================================
  ## pdf('dendrogram.pdf', width=10, height=10)

  n <- length(labels(dend))
  ## circos.par(track.margin = c(.01, .01), track.height = 0.1)

  circos.initialize("foo", xlim = c(0, n))

  thelabels <- labels(dend)
  the_cexs <- list()
  for (i in 1:length(thelabels)) {
    thelabel <- thelabels[i]
    thelength <- nchar(thelabel)

    if (thelength < 13) {
      the_cexs[[i]] <- 1.4

    } else if (thelength >= 13 & thelength < 17 ) {
      the_cexs[[i]] <- 1.0
    } else {

      the_cexs[[i]] <- 0.7
    }
  }
  the_cexs <- unlist(the_cexs)

  ## browser()
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(1:n-0.5, rep(0, n), thelabels, col = labels_colors(dend),
                cex=the_cexs,
                facing = "clockwise", niceFacing = TRUE, adj = c(0.46, 0.5))
  }, bg.border = NA, track.height = 0.1)


  for (i in 1:length(color_numeric_factors)) {
    color_numeric_factor <- color_numeric_factors[[i]]
    ## print(i)
      circos.track(ylim = c(0, 1), cell.padding = c(.02, 1, .02, 1), panel.fun = function(x, y) {
          circos.rect(1:n - 0.8, rep(.5, n), 1:n - 0.2, .3, col = color_numeric_factor, border = NA)
      }, bg.border = NA)
  }

  max_height <- attr(dend, "height")

  circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
    circos.dendrogram(dend, max_height = max_height)
  }, track.height = 0.4, bg.border = NA)


  lgds <- list()
  xs <- c(10, 20, 30)

  ## for (i in 1:length(color)) {

  lgd_subtype <- Legend(
      at = unique(labs[[color[1]]]),
        type = "points",
      legend_gp = gpar(col = unique(color_numeric_factor)),
      title_position = "topleft", title = ""
  )
  lgd_list_vertical <- packLegend(lgd_subtype)
    ## lgds[[i]] <- lgd_list_vertical

    draw(lgd_list_vertical, x = unit(xs[i], "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

  ## browser()
    ## draw(lgds[[1]], lgds[[2]], x = unit(xs[1], "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
    ## draw(lgds[2], x = unit(xs[2], "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

  ## circlize_plot()
  ## draw(lgd_list_vertical, x = unit(10, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
  ## draw(lgds, x = unit(10, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

  ## circos.par(track.margin = c(.01, .01), track.height = 0.3)
  ## dev.off()

  ## ====================================================================================
}
