suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))

.correlation_annotation_colors <- function(metadata, configured = list()) {
  colors <- list()
  base_palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )

  for (field in colnames(metadata)) {
    values <- metadata[[field]]
    if (is.numeric(values)) {
      finite_values <- values[is.finite(values)]
      if (length(finite_values) == 0) next
      value_range <- range(finite_values)
      if (value_range[[1]] == value_range[[2]]) {
        value_range <- value_range + c(-0.5, 0.5)
      }
      colors[[field]] <- circlize::colorRamp2(
        value_range,
        c("#dbe9f6", "#18375b")
      )
      next
    }

    values <- as.character(values)
    values[is.na(values) | values == ""] <- "NA"
    metadata[[field]] <- values
    levels <- sort(unique(values))
    mapping <- character()
    if (!is.null(configured[[field]])) {
      mapping <- unlist(configured[[field]], use.names = TRUE)
    }
    missing_levels <- setdiff(levels, names(mapping))
    if (length(missing_levels) > 0) {
      generated <- if (length(missing_levels) <= length(base_palette)) {
        base_palette[seq_along(missing_levels)]
      } else {
        grDevices::hcl.colors(length(missing_levels), palette = "Dark 3")
      }
      names(generated) <- missing_levels
      mapping <- c(mapping, generated)
    }
    if ("NA" %in% names(mapping)) mapping[["NA"]] <- "#bdbdbd"
    colors[[field]] <- mapping[levels]
  }
  list(metadata = metadata, colors = colors)
}

.correlation_color_function <- function(mat, metric) {
  if (metric %in% c("pearson", "spearman")) {
    return(circlize::colorRamp2(
      c(-1, 0, 1),
      c("#2166ac", "#f7f7f7", "#b2182b")
    ))
  }
  max_value <- max(mat[is.finite(mat)], na.rm = TRUE)
  if (!is.finite(max_value) || max_value <= 0) max_value <- 1
  circlize::colorRamp2(
    c(0, max_value / 2, max_value),
    c("#ffffcc", "#fd8d3c", "#800026")
  )
}

.correlation_text_color <- function(fill) {
  rgb <- grDevices::col2rgb(fill, alpha = FALSE) / 255
  linear <- ifelse(
    rgb <= 0.04045,
    rgb / 12.92,
    ((rgb + 0.055) / 1.055)^2.4
  )
  luminance <- 0.2126 * linear[1, ] + 0.7152 * linear[2, ] + 0.0722 * linear[3, ]
  ifelse(luminance < 0.42, "white", "black")
}

.correlation_value_label <- function(value) {
  ifelse(is.finite(value), sprintf("%.2f", value), "NA")
}

.correlation_layer_fun <- function(annotate, mat, overlap_counts, fontsize = 9) {
  if (!isTRUE(annotate)) return(NULL)
  force(mat)
  force(overlap_counts)
  force(fontsize)
  function(j, i, x, y, width, height, fill, slice_r, slice_c) {
    if (slice_r > slice_c) {
      lower_triangle <- rep(TRUE, length(i))
    } else if (slice_r < slice_c) {
      lower_triangle <- rep(FALSE, length(i))
    } else {
      # i/j are original matrix indices. Derive ranks from the rendered cell
      # coordinates so the triangle remains correct after clustering.
      x_value <- round(grid::convertX(x, "npc", valueOnly = TRUE), 10)
      y_value <- round(grid::convertY(y, "npc", valueOnly = TRUE), 10)
      column_position <- match(x_value, sort(unique(x_value)))
      row_position <- match(y_value, sort(unique(y_value), decreasing = TRUE))
      lower_triangle <- row_position > column_position
    }

    values <- mat[cbind(i, j)]
    counts <- overlap_counts[cbind(i, j)]
    labels <- .correlation_value_label(values)
    labels[lower_triangle] <- paste0("n = ", counts[lower_triangle])
    grid::grid.text(
      labels,
      x,
      y,
      gp = grid::gpar(
        fontsize = fontsize,
        fontface = "bold",
        col = .correlation_text_color(fill)
      )
    )
  }
}

.correlation_distance_for_hclust <- function(distance_matrix) {
  distance_matrix <- as.matrix(distance_matrix)
  storage.mode(distance_matrix) <- "double"
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be square")
  }
  if (is.null(rownames(distance_matrix)) || is.null(colnames(distance_matrix)) ||
      !identical(rownames(distance_matrix), colnames(distance_matrix))) {
    stop("distance_matrix rows and columns must have identical sample ordering")
  }
  if (!isTRUE(all.equal(
    distance_matrix,
    t(distance_matrix),
    tolerance = sqrt(.Machine$double.eps),
    check.attributes = FALSE
  ))) {
    stop("distance_matrix must be symmetric")
  }

  # The metric/dissimilarity has already been computed. Only make undefined
  # pairs usable by hclust and remove negligible negative roundoff; do not
  # derive or transform a distance from the displayed matrix.
  finite_off_diagonal <- distance_matrix[
    upper.tri(distance_matrix) & is.finite(distance_matrix)
  ]
  replacement <- if (length(finite_off_diagonal) > 0) {
    max(finite_off_diagonal)
  } else {
    1
  }
  if (!is.finite(replacement) || replacement <= 0) replacement <- 1
  distance_matrix[!is.finite(distance_matrix)] <- replacement
  tolerance <- sqrt(.Machine$double.eps)
  if (any(distance_matrix < -tolerance)) {
    stop("distance_matrix cannot contain negative dissimilarities")
  }
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0
  distance_matrix
}

.correlation_single_leaf <- function(label) {
  structure(
    1L,
    members = 1L,
    midpoint = 0,
    height = 0,
    label = as.character(label),
    leaf = TRUE,
    class = "dendrogram"
  )
}

.correlation_resolve_linkage <- function(metric, linkage) {
  if (identical(linkage, "auto")) {
    return(if (identical(metric, "l1")) "average" else "ward.D2")
  }
  linkage
}

.correlation_hclust <- function(distance_matrix, ids, linkage) {
  ids <- as.character(ids)
  if (length(ids) == 0) stop("No samples were supplied for clustering")
  if (length(ids) == 1) return(.correlation_single_leaf(ids[[1]]))
  subset <- distance_matrix[ids, ids, drop = FALSE]
  choices <- c(
    "single", "complete", "average", "weighted", "centroid", "median", "ward.D2"
  )
  linkage <- match.arg(linkage, choices)
  # cluster2 calls this option "weighted"; base R names the same WPGMA
  # implementation "mcquitty".
  effective_linkage <- if (identical(linkage, "weighted")) "mcquitty" else linkage
  # The selected dissimilarity is supplied by Python and is passed directly to
  # hclust; no metric is inferred from the values shown in the heatmap.
  stats::hclust(stats::as.dist(subset), method = effective_linkage)
}

.correlation_clustering_functions <- function(distance_matrix, cluster, linkage) {
  if (!isTRUE(cluster)) {
    return(list(rows = FALSE, columns = FALSE))
  }

  # ComplexHeatmap invokes clustering separately for rows, columns, and split
  # slices. Cache each sample subset so both axes reuse the same dendrogram and
  # therefore the same ordering from the supplied dissimilarity matrix.
  cache <- new.env(parent = emptyenv())
  cluster_samples <- function(x) {
    ids <- as.character(rownames(x))
    key <- paste(c(length(ids), ids), collapse = "\r")
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(
        key,
        .correlation_hclust(distance_matrix, ids, linkage),
        envir = cache
      )
    }
    get(key, envir = cache, inherits = FALSE)
  }
  list(
    rows = cluster_samples,
    # ComplexHeatmap passes the transposed matrix to a column-clustering
    # function, so the sample identifiers are row names here as well.
    columns = cluster_samples
  )
}

correlation_heatmap <- function(
  metric_matrix,
  distance_matrix = NULL,
  overlap_counts = NULL,
  metadata = NULL,
  metric = c("l2", "l1", "pearson", "spearman"),
  linkage = c("auto", "ward.D2", "complete", "average", "single", "weighted", "centroid", "median"),
  cut_by = NULL,
  cluster = FALSE,
  annotate = FALSE,
  outname = "correlation",
  outfiletypes = c(".pdf"),
  metadata_colors_json = "{}",
  title = NULL,
  fig_width = NULL,
  fig_height = NULL,
  png_res = 300
) {
  metric <- match.arg(metric)
  linkage <- .correlation_resolve_linkage(metric, match.arg(linkage))
  mat <- as.matrix(metric_matrix)
  storage.mode(mat) <- "double"
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("metric_matrix must have sample row and column names")
  }
  if (!setequal(rownames(mat), colnames(mat))) {
    stop("metric_matrix row and column sample names differ")
  }
  mat <- mat[colnames(mat), colnames(mat), drop = FALSE]

  if (is.null(distance_matrix)) {
    stop("distance_matrix is required; correlation_heatmap never recomputes dissimilarity from metric_matrix")
  }
  distance_matrix <- as.matrix(distance_matrix)
  if (!setequal(rownames(distance_matrix), colnames(mat)) ||
      !setequal(colnames(distance_matrix), colnames(mat))) {
    stop("distance_matrix sample names differ from metric_matrix")
  }
  distance_matrix <- distance_matrix[colnames(mat), colnames(mat), drop = FALSE]
  distance_matrix <- .correlation_distance_for_hclust(distance_matrix)

  if (is.null(overlap_counts)) {
    overlap_counts <- matrix(NA_integer_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  } else {
    overlap_counts <- as.matrix(overlap_counts)
    overlap_counts <- overlap_counts[rownames(mat), colnames(mat), drop = FALSE]
  }

  metadata_df <- NULL
  configured_colors <- list()
  if (!is.null(metadata) && nrow(metadata) > 0) {
    metadata_df <- as.data.frame(metadata, check.names = FALSE, stringsAsFactors = FALSE)
    if ("sample" %in% colnames(metadata_df)) {
      rownames(metadata_df) <- as.character(metadata_df$sample)
      metadata_df$sample <- NULL
    }
    metadata_df <- metadata_df[colnames(mat), , drop = FALSE]
    if (ncol(metadata_df) == 0) metadata_df <- NULL
  }

  if (!is.null(metadata_colors_json) && nzchar(metadata_colors_json) &&
      requireNamespace("jsonlite", quietly = TRUE)) {
    configured_colors <- jsonlite::fromJSON(
      metadata_colors_json,
      simplifyVector = FALSE
    )
  }

  annotation_colors <- list()
  top_annotation <- NULL
  left_annotation <- NULL
  split_df <- NULL
  if (!is.null(metadata_df)) {
    annotation_spec <- .correlation_annotation_colors(metadata_df, configured_colors)
    metadata_df <- annotation_spec$metadata
    annotation_colors <- annotation_spec$colors
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df = metadata_df,
      col = annotation_colors,
      show_annotation_name = FALSE,
      na_col = "#bdbdbd"
    )
    left_annotation <- ComplexHeatmap::rowAnnotation(
      df = metadata_df,
      col = annotation_colors,
      show_annotation_name = FALSE,
      na_col = "#bdbdbd"
    )

    if (!is.null(cut_by) && length(cut_by) > 0) {
      missing_cut <- setdiff(cut_by, colnames(metadata_df))
      if (length(missing_cut) > 0) {
        stop("cut_by field(s) absent from metadata: ", paste(missing_cut, collapse = ", "))
      }
      split_df <- metadata_df[, cut_by, drop = FALSE]
      for (field in colnames(split_df)) {
        values <- as.character(split_df[[field]])
        split_df[[field]] <- factor(values, levels = unique(values), ordered = TRUE)
      }
    }
  } else if (!is.null(cut_by) && length(cut_by) > 0) {
    stop("cut_by requires metadata")
  }

  n_samples <- ncol(mat)
  n_annotations <- if (is.null(metadata_df)) 0 else ncol(metadata_df)
  cell_inches <- if (isTRUE(annotate)) 0.30 else 0.22
  if (is.null(fig_width)) {
    fig_width <- max(7, 3.5 + cell_inches * n_samples + 0.32 * n_annotations)
  }
  if (is.null(fig_height)) {
    fig_height <- max(7, 3.0 + cell_inches * n_samples + 0.32 * n_annotations)
  }

  heatmap_colors <- .correlation_color_function(mat, metric)
  legend_title <- if (identical(metric, "l2")) {
    "L2 RMS\ndistance"
  } else if (identical(metric, "l1")) {
    "L1 mean absolute\ndistance"
  } else {
    paste0(tools::toTitleCase(metric), "\ncorrelation")
  }

  clustering_functions <- .correlation_clustering_functions(
    distance_matrix,
    cluster,
    linkage
  )
  cell_fontsize <- if (n_samples > 30) {
    4.2
  } else if (n_samples > 20) {
    4.5
  } else if (n_samples > 10) {
    5.25
  } else {
    6.5
  }

  make_heatmap <- function() {
    ComplexHeatmap::Heatmap(
      mat,
      name = legend_title,
      col = heatmap_colors,
      na_col = "#d9d9d9",
      cluster_rows = clustering_functions$rows,
      cluster_columns = clustering_functions$columns,
      row_dend_reorder = FALSE,
      column_dend_reorder = FALSE,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      row_split = split_df,
      column_split = split_df,
      top_annotation = top_annotation,
      left_annotation = left_annotation,
      show_row_names = n_samples <= 100,
      show_column_names = n_samples <= 100,
      row_names_gp = grid::gpar(fontsize = if (n_samples > 50) 6 else 8),
      column_names_gp = grid::gpar(fontsize = if (n_samples > 50) 6 else 8),
      column_names_rot = 45,
      rect_gp = grid::gpar(col = "white", lwd = 0.35),
      layer_fun = .correlation_layer_fun(
        annotate,
        mat,
        overlap_counts,
        fontsize = cell_fontsize
      ),
      heatmap_legend_param = list(title_position = "topcenter")
    )
  }

  outfiletypes <- as.character(outfiletypes)
  written <- character()
  for (ext in outfiletypes) {
    if (!startsWith(ext, ".")) ext <- paste0(".", ext)
    outfile <- paste0(outname, ext)
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    if (identical(ext, ".pdf")) {
      grDevices::cairo_pdf(outfile, width = fig_width, height = fig_height)
    } else if (identical(ext, ".png")) {
      grDevices::png(
        outfile,
        width = fig_width,
        height = fig_height,
        units = "in",
        res = png_res,
        type = "cairo"
      )
    } else if (identical(ext, ".svg")) {
      grDevices::svg(outfile, width = fig_width, height = fig_height)
    } else {
      stop("Unsupported output type: ", ext)
    }
    tryCatch(
      ComplexHeatmap::draw(
        make_heatmap(),
        column_title = title,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        merge_legends = TRUE,
        padding = grid::unit(c(8, 8, 5, 8), "mm")
      ),
      finally = grDevices::dev.off()
    )
    written <- c(written, normalizePath(outfile, mustWork = FALSE))
  }
  written
}
