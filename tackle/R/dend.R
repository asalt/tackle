library(ComplexHeatmap)



dend <- dist( select(data, -GeneID, -GeneSymbol) %>% t, method = 'euclidean') %>%
  hclust(method = 'ward.D2')

col_data <- read_tsv('./tmp_meta.tsv')


dist_clust <- dist( select(data, -GeneID, -GeneSymbol) %>% t, method = 'euclidean') %>%
  hclust(method = 'ward.D2')
dend <- as.dendrogram(dist_clust)

ddata_x <- ggdendro::dendro_data(dend)
## labs <- label(ddata_x) %>% left_join(dfm, by=c("label" = "rowname"))
labs <- ddata_x$label %>% left_join(col_data, by=c("label" = "name"))

dendextend::labels_colors(dend) <- as.factor(labs[['model']]) %>% as.numeric

## palette(brewer.pal(8, 'Dark2'))     # six color rainbow


## pdf(paste0('./results/dendrogram', '.pdf'),
##     width=12, height=12)
## ## dend1 <- color_branches(dend, h = 3)
## circlize_dendrogram(dend, labels_track_height = NA,  dend_track_height = 0.5,)
## dev.off()

subtype_col <- as.factor(labs[['subtype']]) %>% as.numeric



## ====================================================================================
pdf('dendrogram.pdf', width=10, height=10)

n <- length(labels(dend))
circos.initialize("foo", xlim = c(0, n))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(1:n-0.8, rep(.5, n), 1:n-0.2, .6, col = subtype_col, border = NA)
}, bg.border = NA)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:n-0.5, rep(0, n), labels(dend), col = labels_colors(dend),
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA, track.height = 0.1)
max_height = attr(dend, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend, max_height = max_height)
}, track.height = 0.5, bg.border = NA)

lgd_subtype = Legend(at = unique(labs[['subtype']]),
                     type='points',
                     legend_gp=gpar(col=unique(subtype_col)),
                     title_position = "topleft", title = "Subtype")
lgd_list_vertical <- packLegend(lgd_subtype)

## circlize_plot()
draw(lgd_list_vertical, x = unit(10, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))

dev.off()

## ====================================================================================
