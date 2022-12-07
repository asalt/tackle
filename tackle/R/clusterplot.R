# clusterplot.R
# vim: ts=2 expandtab
#
#
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(cluster)

ht_opt$message <- FALSE

myzscore <- function(value, minval = NA, remask = TRUE) {
  mask <- is.na(value)
  if (is.na(minval)) minval <- min(value, na.rm = TRUE)

  if (minval == Inf) {
    minval <- 0
  }

  value[is.na(value)] <- minval
  # todo make smaller than min val
  out <- scale(value)

  # if all NA:
  if (sum(!is.finite(out)) == length(out)) {
    out[, 1] <- 0
  }

  ## if (sum(!is.finite(value)) == length(value)) {
  ##   value[!is.finite(value)] <- 0
  ## }

  if (remask == TRUE) {
    out[mask] <- NA
  }
  return(out)
}



dist_no_na <- function(mat) {
  mat[is.na(mat)] <- min(mat, na.rm = TRUE) * .9
  edist <- dist(mat)
  return(edist)
}

order_colmeta <- function(annot, the_order, name = "X") {
  ## see https://stackoverflow.com/a/26003971 for !!name:=something magic
  res <- annot %>%
    filter(get(name) %in% the_order) %>%
    mutate(!!name := factor(get(name),
      levels = the_order,
      ordered = TRUE
    )) %>%
    arrange(get(name))
  res
}

cluster2 <- function(data, annot_mat = NULL, cmap_name = NULL,
                     the_annotation = NULL,
                     genes = NULL,
                     row_annot_df = NULL,
                     gids_to_annotate = NULL,
                     nclusters = NULL,
                     cluster_func = NULL,
                     show_gene_symbols = FALSE,
                     z_score = NULL, z_score_by = NULL,
                     standard_scale = NULL,
                     row_dend_side = "left",
                     row_names_side = "right",
                     mask = NULL, show_missing_values = TRUE, max_autoclusters = 30,
                     row_cluster = TRUE, col_cluster = TRUE, seed = NA,
                     metadata = NULL, col_data = NULL, figsize = NULL, normed = NULL,
                     linkage = "average", gene_symbol_fontsize = 8,
                     metadata_colors = NULL, circle_col_markers = FALSE,
                     circle_col_marker_size = 12,
                     force_plot_genes = FALSE, main_title = "",
                     title_fontsize = 9,
                     cluster_row_slices = TRUE,
                     cluster_col_slices = TRUE,
                     cut_by = NULL,
                     order_by_abundance = FALSE) {
  ht_opt$message <- FALSE
  # preserve column order if col_cluster is disabled
  col_data[["name"]] <- factor(col_data[["name"]], ordered = TRUE, levels = col_data$name)

  if (!is.null(cluster_func)) {
    if (cluster_func == "PAM") {
      cluster_func <- cluster::pam
    } else if (tolower(cluster_func) == "kmeans") cluster_func <- kmeans
    ## ?
  }

  if (!is.null(genes)) {
    data <- data %>% filter(GeneID %in% genes)
  }

  ## print(data %>% pivot_longer(c(-GeneSymbol, -GeneID)))
  ## print(col_data)

  # in case a row is duplicated
  exprs_long <- data %>%
    distinct(GeneID, .keep_all = TRUE) %>%
    pivot_longer(c(-GeneSymbol, -GeneID))

  if (!is.null(col_data)) {
    exprs_long <- exprs_long %>%
      left_join(col_data, by = "name", copy = TRUE) %>%
      mutate(name = factor(name, levels = col_data$name, ordered = TRUE))
  }


  if (is.null(z_score)) {
    # do nothing
  } else if (is.null(z_score_by) & z_score == "0") {
    exprs_long <- exprs_long %>%
      mutate(value = na_if(value, 0)) %>%
      group_by(GeneID) %>%
      mutate(zscore = myzscore(value), zscore_impute = myzscore(value, remask = FALSE)) %>%
      ungroup()
  } else if (!is.null(z_score_by) & z_score == "0") {
    exprs_long <- exprs_long %>%
      mutate(value = na_if(value, 0)) %>%
      group_by(GeneID, !!as.name(z_score_by)) %>%
      mutate(zscore = myzscore(value), zscore_impute = myzscore(value, remask = FALSE)) %>%
      ungroup()
  }
  ## else if (is.null(z_score)) {
  ## }

  if (is.null(z_score)) {
    toplot <- exprs_long %>% pivot_wider(id_cols = c(GeneID, GeneSymbol), values_from = value, names_from = name)
    ## deal with this later
    tocluster <- exprs_long %>% pivot_wider(
      id_cols = c(GeneID, GeneSymbol),
      values_from = value, names_from = name
    )
  } else if (!is.null(z_score)) {
    toplot <- exprs_long %>% pivot_wider(id_cols = c(GeneID, GeneSymbol), values_from = zscore, names_from = name)
    tocluster <- exprs_long %>% pivot_wider(
      id_cols = c(GeneID, GeneSymbol),
      values_from = zscore_impute, names_from = name
    )
  }


  ##   exprs_long %>%
  ## dplyr::group_by(GeneID, GeneSymbol, name) %>%
  ## dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  ## dplyr::filter(n > 1L)


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


  ## ===============  ROW ANNOTATION ============================================
  row_data_args <- NULL
  if (!is.null(row_annot_df)) {
    row_annot_df[["GeneID"]] <- rownames(row_annot_df)
    row_annot_df <- row_annot_df %>%
      mutate(GeneID = factor(row_annot_df[["GeneID"]], levels = data$GeneID, ordered = TRUE)) %>%
      arrange(GeneID)

    row_data_args <- as.list(select(row_annot_df, -GeneID))
    # row_data_args[["na_col"]] <- "white"
    row_data_args[["border"]] <- FALSE
    row_data_args[["which"]] <- "row"
    row_data_args[["annotation_legend_param"]] <- list()
    ## rotate all legends to horizontal (place on bottom via draw, below)

    for (thename in names(select(row_annot_df, -GeneID))) {
      ## row_data_args[["annotation_legend_param"]][[thename]] <- list(direction = "horizontal", nrow = 2)
      row_data_args[["annotation_legend_param"]][[thename]] <- list(direction = "horizontal")
    }
    row_data_args[["annotation_name_side"]] <- "top"
    row_data_args[["gp"]] <- gpar(fontsize = 8, col = NA)
    row_data_args[["annotation_name_gp"]] <- gpar(fontsize = 8)
    row_data_args[["annotation_width"]] <- unit(.15, "in")
  }
  ## ===============  COLUMN ANNOTATION ============================================

  # browser()
  col_data_args <- NULL # will get defined if not is.na col_data
  if (!is.null(col_data)) {
    col_order <- toplot %>%
      select(-GeneID, -GeneSymbol) %>%
      colnames()
    col_data <- col_data %>%
      filter(name %in% colnames(data)) %>%
      mutate(name = factor(name, levels = col_order, ordered = TRUE)) %>%
      arrange(name)


    # col_data_args[["na_col"]] <- "white"


    # col_annot_df <- NULL
    # col_data_args[["annotation_legend_param"]][["at"]][["chr"]] <- as.matrix(colnames(col_data))
    # if (!is.null(col_data) && !is.null(metadata_colors)) {
    #   browser()

    #   # col_annot <- do.call(ComplexHeatmap::HeatmapAnnotation, col_data_args)
    #   col_annot <- columnAnnotation(
    #     df = as.data.frame(col_data),
    #     col = metadata_colors
    #     #annotation_legend_param = col_data_args$annotation_legend_param
    #   )
    # }

    ## Add more args here
    col_data_args <- as.list(col_data %>% select(-name))
    col_data_args <- col_data_args[order(names(col_data_args))]
    # $metadata_colors <- NULL
    ## Custom colers
    if (!is.null(metadata_colors)) {
      col_data_args[["col"]] <- list()
      row_data_args[["col"]] <- list()
      ## print(names(metadata_colors[[1]]))
      for (i in 1:length(metadata_colors)) {
        for (entry_name in names(metadata_colors[[i]])) {
          entry_values <- names(metadata_colors[[i]][[entry_name]])
          ## print(entry_name)
          ## print(entry_values)
          ## print(paste(entry_name, entry_values))
          col_data_args[["col"]][[entry_name]] <- list()
          row_data_args[["col"]][[entry_name]] <- list()
          for (key in entry_values) {
            final_val <- metadata_colors[[i]][[entry_name]][[key]]
            # if (final_val %in% c("NA", "<NA>", "NaN", "nan"))
            ## print(paste(entry_name, key, final_val))
            ## print('****************')
            if (entry_name %in% names(col_data_args)) {
              col_data_args[["col"]][[entry_name]][[key]] <- final_val[1]
            }
            ## if (entry_name == 'IDG') browser()
            if (entry_name %in% names(row_data_args)) { # cannot put color mappings that do not exist
              if (key %in% row_data_args[[entry_name]]) {
                row_data_args[["col"]][[entry_name]][[key]] <- final_val[1]
              }
            }
            ## need to make this atomic
            ## line 456  if(is.atomic(col)) {
          }
          col_data_args[["col"]][[entry_name]] <- unlist(col_data_args[["col"]][[entry_name]])
          row_data_args[["col"]][[entry_name]] <- unlist(row_data_args[["col"]][[entry_name]])
          ## print(is.atomic(col_data_args[['col']][[entry_name]]))
          ## THIS NEEDS TO BE TRUE
        }
      }
    }
  } # if (!is.null(col_data))

  # now handle continuous
  ##
  col_data_simple_annot_names <- intersect(col_data_args %>% names(), col_data %>% names())

  # col_data_args %>% sapply(is.atomic)

  col_data_color_names <- col_data_args$col %>% names()
  # for (thename in col_data_color_names) {
  #   col_data_args[["annotation_legend_param"]][[thename]] <- list(fontsize = 8)
  # }
  .missing_color_definition <- setdiff(col_data_simple_annot_names, col_data_color_names)
  # browser()

  # .PROTECTED_ENTRIES <- c("label", "replicate")
  .PROTECTED_ENTRIES <- c()
  for (.entry in .missing_color_definition) {
    if (.entry %in% .PROTECTED_ENTRIES) {
      next
    }
    .values <- col_data_args[[.entry]]
    if (typeof(.values) == "character") {
      # browser()
      next
    }
    .minval <- min(.values, na.rm = T) * 0.975
    .maxval <- max(.values, na.rm = T) * 1.025
    col_data_args[["col"]][[.entry]] <- circlize::colorRamp2(
      breaks = c(.minval, .maxval),
      colors = c("#829ad3", "#2c303f")
    )
  }
  # missing_in_col_
  # browser()

  # l_simple_anno <- sapply(col_data_args, is.atomic)

  # names(col_data_args)
  # col_data_args$annotation_legend_param

  # print(col_data_args)
  # print(row_data_args)

  ## ========================================================================
  ## Now make the annotations, with all arguments populated
  row_annot <- NULL
  if (!is.null(row_annot_df)) { # only if we have row data to plot
    row_annot <- do.call(ComplexHeatmap::HeatmapAnnotation, row_data_args)
    # TODO continuous variables
  }


  # col_data_args[["annotation_legend_param"]][["at"]][["chr"]] <- as.matrix(colnames(col_data))
  if (!is.null(col_data)) {
    col_annot <- do.call(ComplexHeatmap::HeatmapAnnotation, col_data_args)
    # col_annot <- columnAnnotation(
    #   df = as.data.frame(col_data),
    #   col = list(category = color_mapping_function)
    # annotation_legend_param = col_data_args$annotation_legend_param
    # )
    # all column annotations
    # .thenames <- names(col_annot@annot_list)
    # for (.name in .thenames){
    #   .theannot <- col_annot@anno_list[.name]
    #   if .theannot$
    # }
  }
  ## ========================================================================

  ## print(colorRamp2(c(1,length(levels(col_data$HS_ratio))), c=c('white', 'black'))(1:8))

  if (is.null(z_score)) {
    quantiles <- exprs_long %>%
      select(value) %>%
      quantile(na.rm = TRUE, probs = seq(0, 1, .025))
    ## minval <- quantiles[["2.5%"]]
    minval <- 0
    ## maxval <- quantiles[["97.5%"]]
    maxval <- quantiles[["95%"]]
    col <- colorRamp2(c(minval, 0, maxval), c("#FFFFCC", "orange", "red"))
  } else {
    quantiles <- exprs_long %>%
      select(zscore) %>%
      quantile(na.rm = TRUE, probs = seq(0, 1, .025))
    minval <- quantiles[["2.5%"]]
    maxval <- quantiles[["97.5%"]]
    col <- colorRamp2(c(minval, 0, maxval), c("blue", "white", "red"))
  }

  ## quantiles <- exprs_long %>% select(zscore) %>% quantile(na.rm=TRUE, probs=seq(0,1,.025))
  ## minval <- exprs_long %>% select(zscore) %>% min(na.rm=TRUE) *.95
  ## maxval <- exprs_long %>% select(zscore) %>% max(na.rm=TRUE) *.95
  ## print(minval)
  ## print(maxval)
  ## print(col)

  cell_fun <- NULL
  if (!is.null(annot_mat)) {
    annot_mat <- annot_mat %>%
      mutate(GeneID = factor(GeneID, levels = toplot$GeneID, ordered = TRUE)) %>%
      arrange(GeneID)
    cell_fun <- function(j, i, x, y, width, height, fill) {
      row <- i
      col <- j + 1
      grid.text(sprintf("%.0f", annot_mat[row, col]), x, y, gp = gpar(fontsize = 6.5))
    }
  }

  ## right gene symbol annotations
  gene_annot <- NULL
  if (!is.null(gids_to_annotate)) {
    boolean_ixs <- toplot$GeneID %in% gids_to_annotate
    ixs <- which(boolean_ixs) # note that dplyr pipe %>% to `which` does not work!!
    thelabels <- toplot %>%
      filter(GeneID %in% gids_to_annotate) %>%
      pull(GeneSymbol)

    gene_annot <- ComplexHeatmap::rowAnnotation(
      genes = ComplexHeatmap::anno_mark(
        at = ixs,
        labels = thelabels,
        labels_gp = gpar(fontsize = 9)
      )
    )
  }


  ## now do kmeans / pam clustering if specified
  discrete_clusters <- NULL
  sil_df <- NULL
  if (!is.null(cluster_func)) {
    set.seed(1234)
    ## X <- dplyr::select(tocluster, -GeneID, -GeneSymbol) %>% as.matrix
    X <- tocluster[, 3:length(tocluster)]
    clusters <- cluster_func(X, nclusters)
    discrete_clusters <- cbind(clusters$cluster)
    ## discrete_clusters <- cluster_func(tocluster %>% dplyr::select(-GeneID, -GeneSymbol), nclusters)
    ## TODO fix this!
    ## this is fixed??
    ## row_cluster <- FALSE
    dis <- dist_no_na(X)^2
    # dis <- dist_no_na(X)
    # browser()
    sil <- cluster::silhouette(clusters$cluster, dis)
    dev.new()
    pdf("silhouette.pdf", width = 10, height = 20)
    print(plot(sil, col = "grey"))
    dev.off()
    sil_df <- cbind(toplot$GeneID, toplot$GeneSymbol, sil) %>%
      as.data.frame() %>%
      rename(GeneID = V1, GeneSymbol = V2)
    # arrange(c(cluster, sil_width))
  }



  column_split <- NULL
  # test:
  ## cut_by <- 'response'
  if (!is.null(cut_by)) {
    col_data[[cut_by]] <- col_data[[cut_by]] %>% factor(ordered = TRUE)
    column_split <- col_data[[cut_by]] %>% factor(ordered = TRUE, levels = unique(col_data[[cut_by]]))
    # column_split <- col_data[, c(cut_by, 'name')]
  }


  ## print(head(toplot[col_data$name]))
  print("##################################")
  # print(cluster_col_slices)

  mat <- toplot[col_data$name]
  # library(seriation)
  # library(cluster)
  # o1 <- agnes(dist_no_na(mat), method = "ward")
  # o2 <- agnes(dist_no_na(t(mat)), method = "ward")
  # # o1 <- seriate(dist_no_na(mat), method = "GW_ward")
  # # o2 <- seriate(dist_no_na(t(mat)), method = "GW_ward")
  # browser()
  # o1 <- seriate(o1$diss, method = "GW_ward")
  # o2 <- seriate(o2$diss, method = "GW_ward")
  # row_cluster <- as.dendrogram(o1[[1]])
  # col_cluster <- as.dendrogram(o2[[1]])

  ## ht <- Heatmap(toplot %>% dplyr::select(-GeneID, -GeneSymbol),
  ht <- Heatmap(toplot[col_data$name],
    name = "mat",
    row_split = discrete_clusters,
    ## column_split = cbind(kout_samples$cluster),
    column_split = column_split,
    col = col,
    border = TRUE,
    cell_fun = cell_fun,
    row_dend_side = row_dend_side, #' right',
    row_names_side = row_names_side, #' left',
    ## top_annotation=hc,
    ## right_annotation=gene_annot,
    ## left_annotation = row_annot,
    top_annotation = col_annot,
    column_title_rot = 0,
    cluster_columns = col_cluster,
    cluster_rows = row_cluster,
    cluster_row_slices = cluster_row_slices,
    cluster_column_slices = cluster_col_slices,
    show_row_names = show_gene_symbols,
    clustering_method_rows = linkage,
    clustering_method_columns = linkage,
    clustering_distance_rows = dist_no_na,
    clustering_distance_columns = dist_no_na,
    column_dend_reorder = TRUE,
    row_labels = toplot$GeneSymbol,
    row_names_gp = gpar(fontsize = gene_symbol_fontsize),
    column_names_gp = gpar(fontsize = 9),
    ## border = FALSE,
    column_names_rot = 90,
    column_title = main_title,
    column_title_gp = gpar(fontsize = title_fontsize),
    column_title_side = "top",
    column_names_side = "top",
    show_parent_dend_line = TRUE,
    row_dend_width = unit(1.2, "in"),
    heatmap_legend_param = list(title = ifelse(is.null(z_score), "log(iBAQ)", "zscore")),
    right_annotation = gene_annot,
    left_annotation = row_annot
  )

  ht <- ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = unit(c(10, 2, 2, 2), "mm"))
  ht_row_order <- row_order(ht)

  if (!is.null(annot_mat)) {
    xunit <- ifelse(row_cluster == TRUE, 1, 2.4)
    print(xunit)
    decorate_heatmap_body("mat", {
      grid.text(paste("Annotation:", the_annotation), unit(xunit, "cm"), unit(-5, "mm"))
    })
  }


  ret <- list(heatmap = ht, sil_df, ht_row_order)
  ret
}


## ========================= COLORS ================================
## remove this later!
PAM50_colors <- c(
  "Basal" = "firebrick",
  "Claudin-low" = "yellow",
  "HER2-E" = "plum",
  "HER2" = "plum",
  "LumB" = "powderblue",
  "LumA" = "darkblue",
  "ERROR" = "grey",
  "Normal" = "white",
  "Her2" = "plum"
)
response_colors <- c(
  "no" = "#1D835A",
  "ypCR" = "#B84822",
  "LowerProlif" = "#1768a6",
  "NA" = "white"
)
Lehmann_TNBC_colors <- c(
  "IM"  = "black",
  "BL1" = "red",
  "BL2" = "blue",
  "M"   = "grey",
  "MSL" = "purple",
  "LAR" = "green",
  "UNS" = "white"
)
## =================================================================
