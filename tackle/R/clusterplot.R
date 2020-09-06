library(tidyverse)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(cluster)

myzscore <- function(value, minval=NA, remask=TRUE){
  mask <- is.na(value)
  if (is.na(minval)) minval <- min(value, na.rm=TRUE)
  value[is.na(value)] <- minval
  out = scale(value)
  if (remask==TRUE)
    out[mask] <- NA
  return(out)
}

dist_no_na <- function(mat){
  mat[is.na(mat)] <- min(mat, na.rm=TRUE)
  edist <- dist(mat)
  return(edist)
}

order_colmeta <-function(annot, the_order, name='X'){
  ## see https://stackoverflow.com/a/26003971 for !!name:=something magic
  res <- annot %>% filter(get(name) %in% the_order) %>%
    mutate(!!name:=factor(get(name), levels=the_order, ordered=TRUE)) %>%
    arrange(get(name))
  res
}

cluster2 <- function(data, annot_mat=NA, cmap_name=NA,
                     the_annotation=NA,
                     genes=NA, highlight_gids=NA,
                     highlight_gid_names=NA,
                     gids_to_annotate=NA,
                     nclusters=NA,
                     show_gene_symbols=FALSE,
                     z_score=NA, z_score_by=NA, standard_scale=NA,
                     mask=NA, show_missing_values=TRUE, max_autoclusters=30,
                     row_cluster=TRUE, col_cluster=TRUE, seed=NA,
                     metadata=NA, col_data=NA, figsize=NA, normed=NA,
                     linkage='average', gene_symbol_fontsize=8,
                     legend_include=NA, legend_exclude=NA,
                     metadata_colors=NA, circle_col_markers=FALSE,
                     circle_col_marker_size=12,
                     force_plot_genes=FALSE, main_title='',
                     order_by_abundance=FALSE){

  if (!is.na(genes)){
    data <- data %>% filter(GeneID %in% genes)
  }

  ## print(data %>% pivot_longer(c(-GeneSymbol, -GeneID)))
  ## print(col_data)

  exprs_long <- data %>% pivot_longer(c(-GeneSymbol, -GeneID)) %>%
    left_join(col_data, by='name', copy=TRUE)
  if (z_score == '0' & is.na(z_score_by)){
    exprs_long <- exprs_long %>% mutate(value = na_if(value,0)) %>% group_by(GeneID) %>%
      mutate(zscore = myzscore(value), zscore_impute=myzscore(value, remask=FALSE)) %>%
      ungroup
  }
  else if (z_score == '0' & !is.na(z_score_by)){
    exprs_long <- exprs_long %>% mutate(value = na_if(value,0)) %>% group_by(GeneID, !!as.name(z_score_by)) %>%
      mutate(zscore = myzscore(value), zscore_impute=myzscore(value, remask=FALSE)) %>%
      ungroup
  }

  toplot <- exprs_long %>% pivot_wider(id_cols=c(GeneID, GeneSymbol), values_from=zscore, names_from=name)
  col_order <- toplot %>% select(-GeneID, -GeneSymbol) %>% colnames
  col_data <- col_data %>% mutate(name = factor(name, levels=col_order, ordered=TRUE)) %>% arrange(name)


  ## if ('HS_ratio' %in% colnames(col_data)){
  ##   col_data[['HS_ratio']] <- cut(col_data[['HS_ratio']],
  ##                                 breaks=0:8/8, labels = paste0('<', 1:8/8),
  ##                                 ordered_result=TRUE
  ##                                 )
  ##   HS_ratio=colorRamp2(.breaks, c=c('white', 'black'))
  ##   .numeric_vals <- as.numeric(as.factor(levels(col_data$HS_ratio)))
  ##   .breaks <- c(.numeric_vals[1], .numeric_vals[length(.numeric_vals)])
  ## }

  ## col_anno_df <- data.frame(select(col_data, -name))
  ## rownames(col_anno_df) <- col_data$name
  ## print(col_anno_df)
  ## ## not sure why this does not work
  ## ## col_annot <- HeatmapAnnotation(df=col_anno_df,
  ## col_annot <- HeatmapAnnotation(df=col_anno_df,
  ##                                which='column',
  ##                                na_col='white'
  ##                                )






  ## Add more args here
  col_data_args <- as.list(col_data%>%select(-name))
  ## col_data_args <- as.list(col_data%>%select(response))
  col_data_args[['na_col']] = 'white'

  ## col_data_args[['col']] = list(
  ##   HS_ratio=colorRamp2(c(0:1), c('white', 'grey30')),
  ##   PAM50.Mar2020 = PAM50_colors ,
  ##   PAM50 = PAM50_colors ,
  ##   response = response_colors,
  ##   TNBCtype = Lehmann_TNBC_colors,
  ##   cohort = c(
  ##     HER2='purple',
  ##     Luminal='blue',
  ##     TNBC='red'
  ##   )
  ## )



  ## print(metadata_colors)
  if (!is.na(metadata_colors)){
    col_data_args[['col']] = list()
    ## print(names(metadata_colors[[1]]))
    for (i in 1:length(metadata_colors)){
      for (entry_name in names(metadata_colors[[i]])){
        entry_values <- names(metadata_colors[[i]][[entry_name]])
        ## print(entry_name)
        ## print(entry_values)
        ## print(paste(entry_name, entry_values))
        col_data_args[['col']][[entry_name]] <- list()
        for (key in entry_values){
            final_val <- metadata_colors[[i]][[entry_name]][[key]]
            ## print(paste(entry_name, key, final_val))
            ## print('****************')
            col_data_args[['col']][[entry_name]][[key]] <- final_val[1]
            ## need to make this atomic
            ## line 456  if(is.atomic(col)) {
        }
        col_data_args[['col']][[entry_name]] <- unlist(col_data_args[['col']][[entry_name]])
        ## print(is.atomic(col_data_args[['col']][[entry_name]]))
        ## THIS IS FALSE< needs to be TRUE
      }
    }
  }
  ## print(col_data_args)

  col_annot <- do.call(ComplexHeatmap::HeatmapAnnotation, col_data_args)

  ## print(colorRamp2(c(1,length(levels(col_data$HS_ratio))), c=c('white', 'black'))(1:8))


  quantiles <- exprs_long %>% select(zscore) %>% quantile(na.rm=TRUE, probs=seq(0,1,.025))
  minval <- quantiles[['2.5%']]
  maxval <- quantiles[['97.5%']]
  ## minval <- exprs_long %>% select(zscore) %>% min(na.rm=TRUE) *.95
  ## maxval <- exprs_long %>% select(zscore) %>% max(na.rm=TRUE) *.95
  col = colorRamp2(c(minval,0,maxval), c('blue', "white", 'red'))
  ## print(minval)
  ## print(maxval)
  ## print(col)

  cell_fun <- NA
  if (!is.na(annot_mat)){
    annot_mat <- annot_mat %>% mutate(GeneID = factor(GeneID, levels=toplot$GeneID, ordered=TRUE)) %>%
      arrange(GeneID)
    cell_fun <- function(j, i, x, y, width, height, fill) {
      row = i
      col = j+1
      grid.text(sprintf("%.0f", annot_mat[row, col]), x, y, gp = gpar(fontsize = 6.5))
  }
  }

  ## right gene symbol annotations
  gene_annot <- NULL
  if (!is.na(gids_to_annotate)) {
    boolean_ixs <- toplot$GeneID %in% gids_to_annotate
    ixs <- which(boolean_ixs) # note that dplyr pipe %>% to `which` does not work!!
    thelabels <- toplot %>% filter(GeneID %in% gids_to_annotate) %>% pull(GeneSymbol)

    gene_annot <- ComplexHeatmap::rowAnnotation(
                                    genes = ComplexHeatmap::anno_mark(
                                                              at = ixs,
                                                              labels = thelabels
                                                            )
                                  )
  }

  ht <- Heatmap(toplot %>% dplyr::select(-GeneID, -GeneSymbol),
                name='mat',
                ## row_split = cbind(kout_genes$cluster),
                ## column_split = cbind(kout_samples$cluster),
                col = col,
                cell_fun = cell_fun,
                ## top_annotation=hc,
                ## right_annotation=gene_annot,
                ## left_annotation = row_annot,
                top_annotation=col_annot,
                column_title_rot = 0,
                cluster_columns=col_cluster,
                cluster_rows=row_cluster,
                cluster_column_slices = TRUE,
                show_row_names=show_gene_symbols,
                clustering_method_rows=linkage,
                clustering_method_columns=linkage,
                clustering_distance_rows=dist_no_na,
                clustering_distance_columns=dist_no_na,
                row_labels=toplot$GeneSymbol,
                row_names_gp = gpar(fontsize = gene_symbol_fontsize),
                column_names_gp = gpar(fontsize = 10),
                border=FALSE,
                column_names_rot=90,
                column_title = main_title,
                column_title_side='top',
                column_names_side='top',
                show_parent_dend_line=TRUE,
                row_dend_width = unit(.8, "in"),
                heatmap_legend_param=list(title='zscore'),
                right_annotation=gene_annot,
                )
  draw(ht, padding = unit(c(10, 2, 2, 2), "mm"))
  if (!is.na(annot_mat)){
    decorate_heatmap_body("mat", {
      grid.text(paste('Annotation:', the_annotation), unit(1, "cm"), unit(-5, "mm"))
    })
  }

  ret <- list(heatmap=ht)
  ret

}


## ========================= COLORS ================================
## remove this later!
PAM50_colors <- c(
  "Basal"= "firebrick",
  "Claudin-low"= "yellow",
  "HER2-E"= "plum",
  "HER2"= "plum",
  "LumB"= "powderblue",
  "LumA"= "darkblue",
  "ERROR"= "grey",
  "Normal"= "white",
  "Her2"= "plum")
response_colors <- c(
  "no"= "#1D835A",
  "ypCR"= "#B84822",
  'LowerProlif' = '#1768a6',
  "NA"= "white"
)
Lehmann_TNBC_colors <- c(
  "IM"  = 'black',
  "BL1" = 'red',
  "BL2" = 'blue',
  "M"   = 'grey',
  "MSL" = 'purple',
  "LAR" = 'green',
  'UNS' = 'white'
)
## =================================================================
