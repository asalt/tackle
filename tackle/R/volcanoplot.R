# Load packages
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(graphics))
suppressMessages(library(ggrepel))


yaxis.choices <- c("pValue", "pAdj")
number_by.choices <- c("abs_log2_FC", "log2_FC", "pValue")
direction.choices <- c("both", "up", "down")

volcanoplot <- function(X, max_labels = 35,
                        pch = 21, cex = 0.35,
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
                        marker_cex = 1.4,
                        color_down = "blue",
                        color_up = "red",
                        global_xmax = NULL,
                        global_ymax = NULL,
                        x_label_override = NULL,
                        y_label_override = NULL,
                        verbose = TRUE,
                        ...) {
  POINT_SIZE <- marker_cex

  ploty <- match.arg(yaxis, yaxis.choices)
  number_by <- match.arg(number_by, number_by.choices)
  direction <- match.arg(direction, direction.choices)
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
  if (direction != "up") {
    X[(X[, sig_metric] < sig & X$FC > fc_cutoff & X$log2_FC < 0), "usd"] <- color_down
  }
  if (direction != "down") {
    X[(X[, sig_metric] < sig & X$FC > fc_cutoff & X$log2_FC > 0), "usd"] <- color_up
  }
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
  to_label <- c()
  max_labels <- as.integer(max_labels)
  if (max_labels > 0) {
    label_candidates <- X %>% filter(Sig != "N.S.")
    if (direction == "up") {
      label_candidates <- label_candidates %>% filter(log2_FC > 0)
    } else if (direction == "down") {
      label_candidates <- label_candidates %>% filter(log2_FC < 0)
    }

    if (nrow(label_candidates) > 0) {
      if (number_by == "abs_log2_FC") {
        to_label <- label_candidates %>%
          arrange(desc(abs(log2_FC))) %>%
          head(max_labels) %>%
          rownames()
      } else if (number_by == "log2_FC" && direction == "both") {
        up_n <- ceiling(max_labels / 2)
        down_n <- floor(max_labels / 2)
        to_label_up <- label_candidates %>%
          filter(log2_FC > 0) %>%
          arrange(desc(log2_FC)) %>%
          head(up_n) %>%
          rownames()
        to_label_down <- label_candidates %>%
          filter(log2_FC < 0) %>%
          arrange(log2_FC) %>%
          head(down_n) %>%
          rownames()
        to_label <- c(to_label_down, to_label_up)
      } else if (number_by == "log2_FC") {
        to_label <- label_candidates %>%
          arrange(desc(log2_FC)) %>%
          head(max_labels) %>%
          rownames()
        if (direction == "down") {
          to_label <- label_candidates %>%
            arrange(log2_FC) %>%
            head(max_labels) %>%
            rownames()
        }
      } else if (number_by == "pValue" && direction == "both") {
        up_n <- ceiling(max_labels / 2)
        down_n <- floor(max_labels / 2)
        to_label_up <- label_candidates %>%
          filter(log2_FC > 0) %>%
          arrange(pValue) %>%
          head(up_n) %>%
          rownames()
        to_label_down <- label_candidates %>%
          filter(log2_FC < 0) %>%
          arrange(pValue) %>%
          head(down_n) %>%
          rownames()
        to_label <- c(to_label_down, to_label_up)
      } else if (number_by == "pValue") {
        to_label <- label_candidates %>%
          arrange(pValue) %>%
          head(max_labels) %>%
          rownames()
      }
    }
  }
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
  add_newlines <- function(group) {
    str_replace_all(group, "(\\+|\\-|:)(?=\\S)", " \\1 ")
  }
  group0 <- add_newlines(group0) %>% stringr::str_wrap(18, whitespace_only = TRUE)
  group1 <- add_newlines(group1) %>% stringr::str_wrap(18, whitespace_only = TRUE)
  if ((max_nchar) > 15) annot_size <- annot_size - .5
  if ((max_nchar) > 25) annot_size <- annot_size - .75
  if ((max_nchar) > 35) annot_size <- annot_size - .5
  if ((max_nchar) > 45) annot_size <- annot_size - .5
  if ((max_nchar) > 60) annot_size <- annot_size - .5
  annot_size <- annot_size * annot_cex
  .annot_space <- 0
  if (annot_cex >= 1.8) .annot_space <- .2

  outline_color <- "#444444"
  highlight_outline_color <- "purple"
  highlight_stroke <- 0.8
  fillable_pch <- c(21, 22, 23, 24, 25)
  use_fill <- pch %in% fillable_pch

  base_points <- if (use_fill) {
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
    highlight_points <- geom_point(
      data = X[highlight_mask, ],
      mapping = aes(fill = usd),
      color = highlight_outline_color,
      stroke = highlight_stroke,
      size = POINT_SIZE * cex,
      show.legend = FALSE,
      pch = if (use_fill) pch else 21
    )
  }

  p <- ggplot(X, aes(log2_FC, -log10(get(ploty)), alpha = alpha)) +
    base_points +
    highlight_points +
    scale_alpha_identity() +
    scale_fill_identity() +
    scale_color_identity() +
    ylim(0, ymax) +
    xlim(-xmax, xmax) +
    geom_text_repel(
      data = X[X$label == TRUE, ],
      aes(label = GeneSymbol, alpha = alpha), color = "black", min.segment.length = .15,
      point.padding = 1e-3,
      box.padding = .1, cex = label_cex,
      segment.size = .35, segment.alpha = .4,
      max.overlaps = Inf,
      seed = 1234,
      show.legend = FALSE
    ) +
    annotate("label", c(-xmax - .annot_space, xmax + .annot_space), c(0, 0),
      label = c(group0, group1),
      color = "white",
      fill = c(color_down, color_up),
      linewidth = 0,
      size = annot_size,
      hjust = c(0, 1), vjust = c(0, 0)
    ) +
    labs(
      x = xlabel_full,
      y = ylabel_full,
      caption = footnote
    ) +
    theme_classic() +
    theme(plot.caption = element_text(color = grey(.5), size = 10))

  print(p)
}
