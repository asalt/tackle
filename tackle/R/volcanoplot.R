# Load packages
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(graphics))
suppressMessages(library(ggrepel))


yaxis.choices <- c("pValue", "pAdj")
number_by.choices <- c("abs_log2_FC", "log2_FC", "pValue", "hybrid")
direction.choices <- c("both", "up", "down")
label_size_by.choices <- c("fixed", "density")

.volcano_density_label_sizes <- function(x,
                                         y,
                                         min_size = 2.4,
                                         max_size = 4.0) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_size <- as.numeric(min_size[[1]])
  max_size <- as.numeric(max_size[[1]])
  if (!is.finite(min_size) || min_size <= 0) {
    stop("min_size must be a positive finite number")
  }
  if (!is.finite(max_size) || max_size < min_size) {
    stop("max_size must be finite and greater than or equal to min_size")
  }
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (length(x) == 0) return(numeric())

  midpoint <- (min_size + max_size) / 2
  sizes <- rep(midpoint, length(x))
  usable <- is.finite(x) & is.finite(y)
  n_usable <- sum(usable)
  if (n_usable == 0) return(sizes)
  if (n_usable == 1) {
    sizes[usable] <- max_size
    return(sizes)
  }

  coordinates <- cbind(x[usable], y[usable])
  distances <- as.matrix(stats::dist(coordinates))
  diag(distances) <- Inf
  neighbors <- min(
    3L,
    max(1L, floor(sqrt(n_usable))),
    n_usable - 1L
  )
  local_spacing <- apply(
    distances,
    1,
    function(values) sort(values, partial = neighbors)[[neighbors]]
  )
  spacing_range <- range(local_spacing)
  if (diff(spacing_range) <= sqrt(.Machine$double.eps)) {
    spacing_score <- rep(0.5, n_usable)
  } else {
    spacing_score <- (local_spacing - spacing_range[[1]]) / diff(spacing_range)
  }
  sizes[usable] <- min_size + spacing_score * (max_size - min_size)
  sizes
}

.select_volcano_labels <- function(label_candidates,
                                   max_labels = 35,
                                   number_by = number_by.choices,
                                   direction = direction.choices,
                                   max_labels_left = NULL,
                                   max_labels_right = NULL) {
  number_by <- match.arg(number_by, number_by.choices)
  direction <- match.arg(direction, direction.choices)
  max_labels <- as.integer(max_labels)

  normalize_side_count <- function(value, option_name) {
    if (is.null(value)) return(NULL)
    value <- as.integer(value[[1]])
    if (is.na(value) || value < 0) {
      stop(option_name, " must be a non-negative integer")
    }
    value
  }
  max_labels_left <- normalize_side_count(max_labels_left, "max_labels_left")
  max_labels_right <- normalize_side_count(max_labels_right, "max_labels_right")

  # One explicit side consumes part of the overall --number budget and the
  # opposite side receives the remainder. With both sides explicit, the
  # overall budget is intentionally ignored.
  total_budget <- if (is.na(max_labels)) 0L else max(0L, max_labels)
  if (!is.null(max_labels_left) && is.null(max_labels_right)) {
    max_labels_right <- max(0L, total_budget - max_labels_left)
  } else if (is.null(max_labels_left) && !is.null(max_labels_right)) {
    max_labels_left <- max(0L, total_budget - max_labels_right)
  }

  if (direction == "up") {
    label_candidates <- label_candidates %>% filter(log2_FC > 0)
  } else if (direction == "down") {
    label_candidates <- label_candidates %>% filter(log2_FC < 0)
  }
  if (nrow(label_candidates) == 0) return(character())

  rank_side <- function(candidates, side, count) {
    if (count <= 0 || nrow(candidates) == 0) return(character())
    candidates <- if (side == "left") {
      candidates %>% filter(log2_FC < 0)
    } else {
      candidates %>% filter(log2_FC > 0)
    }
    if (nrow(candidates) == 0) return(character())

    ranked <- if (number_by == "hybrid") {
      candidates %>%
        mutate(
          .fc_rank = min_rank(desc(abs(log2_FC))),
          .p_rank = min_rank(pValue),
          .best_rank = pmin(.fc_rank, .p_rank),
          .combined_rank = .fc_rank + .p_rank
        ) %>%
        arrange(
          .best_rank,
          .combined_rank,
          pValue,
          desc(abs(log2_FC))
        )
    } else if (number_by == "abs_log2_FC") {
      candidates %>% arrange(desc(abs(log2_FC)))
    } else if (number_by == "log2_FC" && side == "left") {
      candidates %>% arrange(log2_FC)
    } else if (number_by == "log2_FC") {
      candidates %>% arrange(desc(log2_FC))
    } else {
      candidates %>% arrange(pValue)
    }
    rownames(head(ranked, count))
  }

  # Preserve the historical --number selection exactly unless a side-specific
  # override is explicitly supplied.
  legacy_labels <- character()
  if (!is.na(max_labels) && max_labels > 0) {
    if (number_by == "abs_log2_FC") {
      legacy_labels <- label_candidates %>%
        arrange(desc(abs(log2_FC))) %>%
        head(max_labels) %>%
        rownames()
    } else if (number_by == "log2_FC" && direction == "both") {
      legacy_labels <- c(
        rank_side(label_candidates, "left", floor(max_labels / 2)),
        rank_side(label_candidates, "right", ceiling(max_labels / 2))
      )
    } else if (number_by == "log2_FC") {
      side <- if (direction == "down") "left" else "right"
      legacy_labels <- rank_side(label_candidates, side, max_labels)
    } else if (number_by == "pValue" && direction == "both") {
      legacy_labels <- c(
        rank_side(label_candidates, "left", floor(max_labels / 2)),
        rank_side(label_candidates, "right", ceiling(max_labels / 2))
      )
    } else if (number_by == "pValue") {
      legacy_labels <- label_candidates %>%
        arrange(pValue) %>%
        head(max_labels) %>%
        rownames()
    } else if (number_by == "hybrid" && direction == "both") {
      legacy_labels <- c(
        rank_side(label_candidates, "left", floor(max_labels / 2)),
        rank_side(label_candidates, "right", ceiling(max_labels / 2))
      )
    } else if (number_by == "hybrid") {
      side <- if (direction == "down") "left" else "right"
      legacy_labels <- rank_side(label_candidates, side, max_labels)
    }
  }

  if (is.null(max_labels_left) && is.null(max_labels_right)) {
    return(legacy_labels)
  }

  legacy_fc <- label_candidates[legacy_labels, "log2_FC", drop = TRUE]
  left_labels <- legacy_labels[legacy_fc < 0]
  right_labels <- legacy_labels[legacy_fc > 0]
  neutral_labels <- legacy_labels[legacy_fc == 0]
  if (!is.null(max_labels_left)) {
    left_labels <- rank_side(label_candidates, "left", max_labels_left)
  }
  if (!is.null(max_labels_right)) {
    right_labels <- rank_side(label_candidates, "right", max_labels_right)
  }
  unique(c(left_labels, neutral_labels, right_labels))
}

volcanoplot <- function(X, max_labels = 35,
                        max_labels_left = NULL,
                        max_labels_right = NULL,
                        pch = 21, cex = 1.0,
                        alpha = 1.,
                        fc_cutoff = 4, sig = 0.05, label_cex = 1,
                        show_all = FALSE, yaxis = yaxis.choices,
                        direction = direction.choices,
                        group0 = "", group1 = "",
                        sig_metric = "pAdj",
                        number_by = "abs_log2_FC",
                        max_fc = NULL,
                        bg_marker_color = "#222222",
                        force_highlight_geneids = FALSE,
                        annot_cex = 1.,
                        marker_cex = 1.0,
                        point_size = NULL,
                        color_down = "blue",
                        color_up = "red",
                        comparison_wrap_width = NULL,
                        global_xmax = NULL,
                        global_ymax = NULL,
                        x_label_override = NULL,
                        y_label_override = NULL,
                        label_size_by = label_size_by.choices,
                        label_size_min = 2.4,
                        label_size_max = 4.0,
                        semantic_svg = FALSE,
                        verbose = TRUE,
                        ...) {
  if (!is.null(point_size)) {
    cex <- point_size
  }
  POINT_SIZE <- marker_cex

  ploty <- match.arg(yaxis, yaxis.choices)
  number_by <- match.arg(number_by, number_by.choices)
  direction <- match.arg(direction, direction.choices)
  label_size_by <- match.arg(label_size_by, label_size_by.choices)
  linear_fc_cutoff <- fc_cutoff

  if (sig_metric == "pAdj") {
    sig_filter_str <- paste0("FDR<", sig)
  } else {
    sig_filter_str <- paste0("p<", sig)
  }

  X$FC <- 2^abs(X[, "log2_FC"])
  # drop na
  X <- X[!is.na(X$t), ]

  Sig <- ifelse(X[, sig_metric] < sig & abs(X[, "FC"]) > fc_cutoff,
    sig_filter_str, "N.S."
  )
  X[, "Sig"] <- Sig

  X[, "usd"] <- bg_marker_color
  X[(X[, sig_metric] < sig & X$FC > fc_cutoff & X$log2_FC < 0), "usd"] <- color_down
  X[(X[, sig_metric] < sig & X$FC > fc_cutoff & X$log2_FC > 0), "usd"] <- color_up
  X[, "alpha"] <- .20 # new column
  if (!"highlight" %in% colnames(X)) {
    X[, "highlight"] <- FALSE
  }
  highlight_mask <- !is.na(X$highlight) & X$highlight
  X[highlight_mask, "alpha"] <- alpha

  X[, "label"] <- FALSE # new column

  pAdj <- X[, "pAdj"][!is.na(X[, "pAdj"])]
  stretch <- min(pAdj[pAdj > 0]) / 2
  X[, "pAdj"] <- X[, "pAdj"] + stretch
  pvalues <- X[, "pValue"][!is.na(X[, "pValue"])]
  stretch <- min(pvalues[pvalues > 0]) / 2
  X[, "pValue"] <- X[, "pValue"] + stretch

  ## ======================================================================
  ## = calculations for deciding on which dots to show ===============
  ## ======================================================================
  label_candidates <- X %>% filter(Sig != "N.S.")
  to_label <- .select_volcano_labels(
    label_candidates,
    max_labels = max_labels,
    number_by = number_by,
    direction = direction,
    max_labels_left = max_labels_left,
    max_labels_right = max_labels_right
  )
  ## ======================================================================

  X[(highlight_mask & X$log2_FC > 0 & X[, sig_metric] < sig), "label"] <- TRUE
  X[(highlight_mask & X$log2_FC < 0 & X[, sig_metric] < sig), "label"] <- TRUE
  ## ======================================================================

  if (isTRUE(force_highlight_geneids)) {
    X[X$highlight == TRUE, "label"] <- TRUE # label these specifically requested genes to be highlighted
  }

  X[to_label, "label"] <- TRUE # label these from FC and pval thresholds

  X[to_label, "alpha"] <- alpha #
  if (show_all == FALSE) {
    X[(X[, "Sig"] == "N.S.") & (X[, "highlight"] == FALSE), "label"] <- FALSE
  }

  ## ymax <- max(-log10(X[, 'pValue'])) * 1.05
  if (is.null(global_xmax) || isTRUE(is.na(global_xmax))) {
    xmax <- X[, "log2_FC"] %>%
      abs() %>%
      max()
  } else {
    xmax <- global_xmax
  }


  ymax <- max(-log10(X[, ploty])) * 1.08
  if (!is.null(global_ymax) && !isTRUE(is.na(global_ymax))) {
    ymax <- global_ymax
  }
  if (isTRUE(verbose)) {
    print(paste0("ymax: ", ymax))
  }

  ## ratio_sig <- paste0( dim( filter(X, Sig == sig_filter_str) )[1], '/', dim(X)[1] )
  if (direction == "both") {
    value_sig <- dim(X[X$Sig == sig_filter_str, ])[1]
  } else if (direction == "up") {
    value_sig <- dim(X[(X$Sig == sig_filter_str) & (X$log2_FC > 0), ])[1]
  } else if (direction == "down") {
    value_sig <- dim(X[(X$Sig == sig_filter_str) & (X$log2_FC < 0), ])[1]
  }

  ratio_sig <- paste0(value_sig, "/", dim(X)[1])

  if (direction == "both") {
    spacer <- ""
  } else if (direction == "up") {
    spacer <- "up"
  } else if (direction == "down") {
    spacer <- "down"
  }
  footnote <- paste(ratio_sig, "sig.", spacer, "at", sig_filter_str)
  if (fc_cutoff != 0) {
    footnote <- paste(footnote, "and", linear_fc_cutoff, "F.C.")
  }

  ylabel_full <- eval(expression(substitute(paste("-log"[10], " ", ploty), list(ploty = ploty))))
  if (!is.null(y_label_override)) {
    ylabel_full <- y_label_override
  }
  x_label_default <- expression(paste("log"[2], " Fold Change"))
  xlabel_full <- if (!is.null(x_label_override)) x_label_override else x_label_default

  annot_size <- 4.0
  max_nchar <- max(nchar(group0), nchar(group1))
  format_group_label <- function(group, width = 30) {
    group <- str_replace_all(group, "_+", " ")
    group <- str_replace_all(group, "\\s*\\n\\s*", " ")
    group <- str_squish(group)
    group <- str_replace_all(group, "\\s*\\+\\s*", " + ")
    group <- str_replace_all(group, "\\s*:\\s*", " : ")

    if (str_count(group, "(?<=\\S)\\s*-\\s*(?=\\S)|(?<=\\S)-(?=\\S)") == 1) {
      parts <- str_split(group, "(?<=\\S)\\s*-\\s*(?=\\S)|(?<=\\S)-(?=\\S)", n = 2, simplify = TRUE)
      left <- str_squish(parts[1])
      right <- str_squish(parts[2])

      if (nzchar(left) && nzchar(right) && str_detect(left, "\\s") && str_detect(right, "\\s")) {
        one_line <- paste(left, right, sep = " - ")
        if (nchar(one_line) <= width) {
          return(one_line)
        }
        if (max(nchar(left), nchar(right)) <= width) {
          return(paste0(left, " -\n", right))
        }
        left <- stringr::str_wrap(left, width = width, whitespace_only = TRUE)
        right <- stringr::str_wrap(right, width = width, whitespace_only = TRUE)
        return(paste0(left, " -\n", right))
      }
    }

    stringr::str_wrap(group, width = width, whitespace_only = TRUE)
  }
  if (!is.null(comparison_wrap_width) && !is.na(comparison_wrap_width) && comparison_wrap_width > 0) {
    wrap_width <- as.integer(comparison_wrap_width)
  } else {
    wrap_width <- 30
    if ((max_nchar) > 45) wrap_width <- 26
    if ((max_nchar) > 70) wrap_width <- 22
  }
  group0 <- format_group_label(group0, width = wrap_width)
  group1 <- format_group_label(group1, width = wrap_width)
  if ((max_nchar) > 15) annot_size <- annot_size - .3
  if ((max_nchar) > 25) annot_size <- annot_size - .5
  if ((max_nchar) > 35) annot_size <- annot_size - .3
  if ((max_nchar) > 45) annot_size <- annot_size - .3
  if ((max_nchar) > 60) annot_size <- annot_size - .3
  # if ((max_nchar) > 80) annot_size <- annot_size - .4
  annot_size <- annot_size * annot_cex
  side_annot_y <- max(ymax * 0.02, 0.15) # 

  X[, "label_size_plot"] <- 3.2 * label_cex
  label_mask <- !is.na(X$label) & X$label
  if (label_size_by == "density" && any(label_mask)) {
    x_denominator <- max(2 * xmax, .Machine$double.eps)
    y_denominator <- max(ymax, .Machine$double.eps)
    label_x <- (X[label_mask, "log2_FC"] + xmax) / x_denominator
    label_y <- -log10(X[label_mask, ploty]) / y_denominator
    X[label_mask, "label_size_plot"] <- .volcano_density_label_sizes(
      label_x,
      label_y,
      min_size = label_size_min,
      max_size = label_size_max
    ) * label_cex
  }

  outline_color <- "#444444"
  highlight_outline_color <- "purple"
  highlight_stroke <- 0.8
  fillable_pch <- c(21, 22, 23, 24, 25)
  use_fill <- pch %in% fillable_pch

  if (isTRUE(semantic_svg) && !requireNamespace("ggiraph", quietly = TRUE)) {
    stop("ggiraph is required for semantic volcano SVG output")
  }

  base_points <- if (isTRUE(semantic_svg) && use_fill) {
    ggiraph::geom_point_interactive(
      data = X[!highlight_mask, ],
      mapping = aes(
        fill = usd,
        data_id = GeneID,
        tooltip = GeneSymbol
      ),
      color = outline_color,
      size = POINT_SIZE * cex,
      show.legend = FALSE,
      pch = pch
    )
  } else if (isTRUE(semantic_svg)) {
    ggiraph::geom_point_interactive(
      data = X[!highlight_mask, ],
      mapping = aes(
        color = usd,
        data_id = GeneID,
        tooltip = GeneSymbol
      ),
      size = POINT_SIZE * cex,
      show.legend = FALSE,
      pch = pch
    )
  } else if (use_fill) {
    geom_point(
      mapping = aes(fill = usd),
      color = outline_color,
      size = POINT_SIZE * cex,
      show.legend = FALSE,
      pch = pch
    )
  } else {
    geom_point(
      mapping = aes(color = usd),
      size = POINT_SIZE * cex,
      show.legend = FALSE,
      pch = pch
    )
  }

  highlight_points <- NULL
  if (any(highlight_mask)) {
    highlight_points <- if (isTRUE(semantic_svg)) {
      ggiraph::geom_point_interactive(
        data = X[highlight_mask, ],
        mapping = aes(
          fill = usd,
          data_id = GeneID,
          tooltip = GeneSymbol
        ),
        color = highlight_outline_color,
        stroke = highlight_stroke,
        size = POINT_SIZE * cex,
        show.legend = FALSE,
        pch = if (use_fill) pch else 21
      )
    } else {
      geom_point(
        data = X[highlight_mask, ],
        mapping = aes(fill = usd),
        color = highlight_outline_color,
        stroke = highlight_stroke,
        size = POINT_SIZE * cex,
        show.legend = FALSE,
        pch = if (use_fill) pch else 21
      )
    }
  }

  label_points <- if (isTRUE(semantic_svg)) {
    ggiraph::geom_text_repel_interactive(
      data = X[X$label == TRUE, ],
      aes(
        label = GeneSymbol,
        alpha = alpha,
        size = label_size_plot,
        data_id = GeneID,
        tooltip = GeneSymbol
      ),
      color = "black",
      min.segment.length = .15,
      point.padding = 1e-3,
      box.padding = .1,
      fontface = "bold",
      segment.size = .35,
      segment.alpha = .4,
      max.overlaps = Inf,
      seed = 1234,
      show.legend = FALSE
    )
  } else {
    geom_text_repel(
      data = X[X$label == TRUE, ],
      aes(
        label = GeneSymbol,
        alpha = alpha,
        size = label_size_plot
      ),
      color = "black",
      min.segment.length = .15,
      point.padding = 1e-3,
      box.padding = .1,
      fontface = "bold",
      segment.size = .35,
      segment.alpha = .4,
      max.overlaps = Inf,
      seed = 1234,
      show.legend = FALSE
    )
  }

  p <- ggplot(X, aes(log2_FC, -log10(get(ploty)), alpha = alpha)) +
    base_points +
    highlight_points +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_fill_identity() +
    scale_color_identity() +
    coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax), clip = "off") +
    label_points +
    annotate("text", c(-xmax * 0.98, xmax * 0.98), c(side_annot_y, side_annot_y),
      label = c(group0, group1),
      color = c(color_down, color_up),
      fontface = "bold",
      lineheight = 0.95,
      size = annot_size,
      hjust = c(0, 1), vjust = c(0, 0)
    ) +
    labs(
      x = xlabel_full,
      y = ylabel_full,
      caption = footnote
    ) +
    theme_classic(base_size=18) +
    theme(
      plot.caption = element_text(color = grey(.5), size = 12),
      plot.margin = margin(5.5, 18, 10, 18)
    )

  print(p)
}
