# Load packages
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))

## library(ggplot2)
library(graphics)

# Install ggrepel package if needed
# install.packages("ggrepel")

makeFootnote <- function(footnoteText = format(Sys.time(), "%d %b %Y"),
                         size = 1.0, color = grey(.5)) {
  suppressMessages(
    require(grid)
  )
  pushViewport(viewport())
  grid.text(
    label = footnoteText,
    x = unit(1, "npc") - unit(2, "mm"),
    y = unit(2, "mm"),
    just = c("right", "bottom"),
    gp = gpar(cex = size, col = color)
  )
  popViewport()
}
yaxis.choices <- c("pValue", "pAdj")
number_by.choices <- c("abs_log2_FC", "log2_FC", "pValue")

# :fc_cutoff: cutoff for absolute fold change cutoff
volcanoplot <- function(X, max_labels = 35,
                        pch = 21, cex = 0.35,
                        alpha = 1.,
                        fc_cutoff = 4, sig = 0.05, label_cex = 1,
                        show_all = FALSE, yaxis = yaxis.choices,
                        group0 = "", group1 = "",
                        sig_metric = "pAdj",
                        number_by = "abs_log2_FC",
                        max_fc = NULL,
                        # bg_marker_color = "#22222288",
                        bg_marker_color = "#222222",
                        force_highlight_geneids = FALSE,
                        annot_cex = 1.,
                        marker_cex = 1.4,
                        ...) {
  POINT_SIZE <- marker_cex
  # pch <- 21 # if we want a fillable shape
  # pch <- 16 # keep it as fillable shape

  ploty <- match.arg(yaxis, yaxis.choices)
  number_by <- match.arg(number_by, number_by.choices)
  linear_fc_cutoff <- fc_cutoff
  ## if (fc_cutoff > 0) {
  ##   fc_cutoff <- abs(log2(fc_cutoff))
  ## }

  if (sig_metric == "pAdj") {
    sig_filter_str <- paste0("FDR<", sig)
  } else {
    sig_filter_str <- paste0("p<", sig)
  }

  ## X <- mutate(X, Sig = ifelse(X$pAdj < sig& abs(X[, 'log2_Fold_Change']) > fc_cutoff,
  ##                             sig_filter_str, "N.S."))
  ## library(tibble)
  ## browser()
  ## X <- as.tibble(X)
  ## X <- mutate(X, Sig = ifelse(X[,sig_metric] < sig& abs(X[, 'log2_Fold_Change']) > fc_cutoff,
  ##                             sig_filter_str, "N.S."))

  X$FC <- 2^abs(X[, "log2_FC"])
  # drop na
  X <- X[!is.na(X$t), ]

  Sig <- ifelse(X[, sig_metric] < sig & abs(X[, "FC"]) > fc_cutoff,
    sig_filter_str, "N.S."
  )
  X[, "Sig"] <- Sig


  ## X[ , 'usd' ] = '#222222bb'
  ## X[ , 'usd' ] = '#22222222'
  # X[, "usd"] <- "#22222288"
  X[, "usd"] <- bg_marker_color
  ## X[ , 'usd' ] = '#88888888'
  ## X[ (X$pAdj < sig & X$log2_Fold_Change > fc_cutoff), 'usd' ] = 'red'
  ## X[ (X$pAdj < sig & X$log2_Fold_Change < -fc_cutoff), 'usd' ] = 'blue'
  X[(X[, sig_metric] < sig & X$FC > fc_cutoff & X$log2_FC < 0), "usd"] <- "blue"
  X[(X[, sig_metric] < sig & X$FC > fc_cutoff & X$log2_FC > 0), "usd"] <- "red"
  X[, "alpha"] <- .20 # new column
  ## X[ X$highlight == TRUE, 'usd' ] = "#67ff3d"
  ## X[ X$highlight == TRUE, 'usd' ] = "purple"
  ## X[ X$highlight == TRUE, 'usd' ] = "#00ab25"
  # .sig <- ifelse(force_highlight_geneids == TRUE, 1.1, sig)
  # if (force_highlight_geneids == TRUE) .sig <- 1.1

  # X[(X$highlight == TRUE & X$log2_FC > 0 & X[, sig_metric] < sig), "usd"] <- "red"
  # X[(X$highlight == TRUE & X$log2_FC < 0 & X[, sig_metric] < sig), "usd"] <- "blue"

  X[(X$highlight == TRUE & X$log2_FC > 0 & X[, sig_metric] < sig), "usd"] <- "purple"
  X[(X$highlight == TRUE & X$log2_FC < 0 & X[, sig_metric] < sig), "usd"] <- "purple"


  X[(X$highlight == TRUE & X$log2_FC > 0 & X[, sig_metric] > sig), "usd"] <- "purple"
  X[(X$highlight == TRUE & X$log2_FC < 0 & X[, sig_metric] > sig), "usd"] <- "purple"
  X[, "usd"] <- as.factor(X[, "usd"])

  # print(X[X$GeneID=="659985970",])
  # wait we don't need this here?
  # do we?
  X[(X$highlight == TRUE & X$log2_FC > 0 & X[, sig_metric] < sig), "alpha"] <- alpha
  X[(X$highlight == TRUE & X$log2_FC < 0 & X[, sig_metric] < sig), "alpha"] <- alpha
  X[(X$highlight == TRUE & X[, sig_metric] > sig), "alpha"] <- alpha # force highlight set alpha to 1

  X$outline_width <- 0.2
  X$outline <- "black"
  X[, "label"] <- FALSE # new column



  ## X <- mutate(X, label = ifelse(X$qValue < 0.05, "FDR<0.05", "N.S."))

  pAdj <- X[, "pAdj"][!is.na(X[, "pAdj"])]
  stretch <- min(pAdj[pAdj > 0]) / 2
  X[, "pAdj"] <- X[, "pAdj"] + stretch
  pvalues <- X[, "pValue"][!is.na(X[, "pValue"])]
  stretch <- min(pvalues[pvalues > 0]) / 2
  X[, "pValue"] <- X[, "pValue"] + stretch

  ## ======================================================================
  ## = calculations for deciding on which dots to show ===============
  ## ======================================================================
  if (number_by == "abs_log2_FC") {
    to_label <- head(
      order(abs(X[X$highlight == TRUE, "log2_FC"]), X[X$highlight == TRUE, "pValue"],
        decreasing = c(TRUE, FALSE)
      ),
      max_labels
    )
    max_labels <- max_labels - length(to_label)

    to_label <- head(
      order(abs(X[, "log2_FC"]), X[, "pValue"], decreasing = c(TRUE, FALSE)),
      max_labels
    )
  } else if (number_by == "log2_FC") {
    to_label1 <- X %>%
      filter(log2_FC < 0 & FC > fc_cutoff & Sig != "N.S.") %>%
      arrange(log2_FC) %>%
      head(round(max_labels, 2)) %>%
      rownames()
    to_label2 <- X %>%
      filter(log2_FC > 0 & FC > fc_cutoff & Sig != "N.S.") %>%
      arrange(-log2_FC) %>%
      head(round(max_labels, 2)) %>%
      rownames()


    to_label <- c(to_label1, to_label2)
  } else {
    to_label1 <- X %>%
      filter(log2_FC < 0 & FC > fc_cutoff & Sig != "N.S.") %>%
      arrange(pValue) %>%
      head(round(max_labels, 2)) %>%
      rownames()
    to_label2 <- X %>%
      filter(log2_FC > 0 & FC > fc_cutoff & Sig != "N.S.") %>%
      arrange(pValue) %>%
      head(round(max_labels, 2)) %>%
      rownames()
    to_label <- c(to_label1, to_label2)
  }

  X[(X$highlight == TRUE & X$log2_FC > 0 & X[, sig_metric] < sig), "label"] <- TRUE
  X[(X$highlight == TRUE & X$log2_FC < 0 & X[, sig_metric] < sig), "label"] <- TRUE
  ## ======================================================================

  # X[(X$highlight == TRUE & X[, sig_metric] < .sig), "label"] <- TRUE
  # X[(X$highlight == TRUE & X[, sig_metric] < .sig), "outline"] <- "000000"
  # X[(X$highlight == TRUE & X[, sig_metric] < .sig), "outline_width"] <- .7
  # X[(X$highlight == TRUE & X[, sig_metric] < .sig), ] %>% print()

  if (force_highlight_geneids == TRUE) {
    X[X$highlight == TRUE, "label"] <- TRUE # label these specifically requested genes to be highlighted
  }

  X[to_label, "label"] <- TRUE # label these from FC and pval thresholds

  # print(X[, X$label == TRUE])

  X[to_label, "alpha"] <- alpha #
  # if we want to force highlight and circle these
  # X[to_label, "highlight"] <- 1. # label these
  if (show_all == FALSE) {
    # X[X[, "Sig"] == "N.S.", "label"] <- FALSE
    X[(X[, "Sig"] == "N.S.") & (X[, "highlight"] == FALSE), "label"] <- FALSE
  }
  # except these? # done abocve
  # X[X$highlight == TRUE, "label"] <- TRUE

  ## ymax <- max(-log10(X[, 'pValue'])) * 1.05
  ymax <- max(-log10(X[, ploty])) * 1.08
  xmax <- X[, "log2_FC"] %>%
    abs() %>%
    max()

  # if (!is.null(max_fc)) {
  #   xmax <- max_fc + .2
  # }

  ## ratio_sig <- paste0( dim( filter(X, Sig == sig_filter_str) )[1], '/', dim(X)[1] )
  ratio_sig <- paste0(dim(X[X$Sig == sig_filter_str, ])[1], "/", dim(X)[1])
  footnote <- paste(ratio_sig, "sig. at", sig_filter_str)
  if (fc_cutoff != 0) {
    footnote <- paste(footnote, "and", linear_fc_cutoff, "F.C.")
  }

  ## footnote <- ''
  ylabel_full <- eval(expression(substitute(paste("-log"[10], " ", ploty), list(ploty = ploty))))

  annot_size <- 4.0
  max_nchar <- max(nchar(group0), nchar(group1))
  add_newlines <- function(group) {
    str_replace(group, "\\+", " \\+\n") %>%
      str_replace(":", "\n") %>%
      str_replace("\\-", " \\-\n")
  }
  group0 <- add_newlines(group0)
  group1 <- add_newlines(group1)
  if ((max_nchar) > 15) annot_size <- annot_size - .5
  if ((max_nchar) > 25) annot_size <- annot_size - .75
  if ((max_nchar) > 35) annot_size <- annot_size - .5
  if ((max_nchar) > 45) annot_size <- annot_size - .5
  if ((max_nchar) > 60) annot_size <- annot_size - .5
  annot_size <- annot_size * annot_cex
  .annot_space <- 0
  if (annot_cex >= 1.8) .annot_space <- .2 # why this?
  # print(annot_size)

  .add_border <- function() {
    geom_point(
      data = X[X$highlight == TRUE, ],
      mapping = aes(color = outline, stroke = outline_width, fill = usd),
      pch = 21, size = POINT_SIZE, cex = cex, show.legend = FALSE
    )
  }
  .add_borderless <- function() {
    geom_point(
      data = X[X$highlight == TRUE, ],
      mapping = aes(color = usd, fill = usd),
      pch = 16, size = POINT_SIZE, cex = cex, show.legend = FALSE
    )
  }
  print(pch)

  if (pch == 21) .geom_point_highlight_TRUE <- .add_border
  if (pch == 16) .geom_point_highlight_TRUE <- .add_borderless

  #     geom_point(mapping=aes(color=usd), size = POINT_SIZE, cex = cex, show.legend = FALSE, pch = 16) +
  #     ## g eom_point(data = X[X$highlight == TRUE, ], mapping=aes(color=outline, stroke=outline_width, fill=usd), pch=21, size = POINT_SIZE, cex = cex, show.legend = FALSE, pch = pch) +
  #     geom_point(data = X[X$highlight == TRUE, ], size = POINT_SIZE, cex = cex, show.legend = FALSE, pch = pch, alpha=1.) +



  #  geom_point(mapping=aes(color=usd), size = POINT_SIZE, cex = cex, show.legend = FALSE, pch = 16) +

  # geom_point(data = X[X$highlight == TRUE, ],
  #            mapping=aes(color=outline, stroke=outline_width , fill=  usd ),
  #            pch=21, size = POINT_SIZE, cex = cex, show.legend = FALSE)+


  p <- ggplot(X, aes(log2_FC, -log10(get(ploty)), alpha = alpha, color = usd)) +
    geom_point(mapping = aes(color = usd), size = POINT_SIZE, cex = cex, show.legend = FALSE, pch = 16) +
    .geom_point_highlight_TRUE() +
    scale_alpha_identity() +
    scale_fill_identity() +
    scale_color_identity() +
    ylim(0, ymax) +
    xlim(-xmax, xmax) +
    geom_text_repel(
      data = X[X$label == TRUE, ],
      aes(label = GeneSymbol, alpha = alpha, color = usd), min.segment.length = .15,
      point.padding = 1e-3,
      box.padding = .1, cex = label_cex,
      segment.size = .35, segment.alpha = .4,
      max.overlaps = Inf,
      seed = 1234,
      show_legend = FALSE,
    ) +
    # annotate("text",  c(-xmax, xmax), c(ymax*.98, ymax*.98), label = c(group0, group1),
    annotate("text", c(-xmax - .annot_space, xmax + .annot_space), c(0, 0),
      label = c(group0, group1),
      size = annot_size,
      hjust = c(0, 1), vjust = c(0, 0), color = c("blue", "red")
    ) +
    labs(
      x = expression(paste("log"[2], " Fold Change")),
      y = ylabel_full,
      caption = footnote
    ) +
    theme_classic() +
    theme(plot.caption = element_text(color = grey(.5), size = 10))

  # theme_base() +


  print(p)

  ## print(X %>% filter(Sig=='FDR<0.05') %>% dim)
  ## print(X %>% dim)

  ## footnote <- paste( ratio_sig, 'sig. with', 2**fc_cutoff, 'F.C.' )
  ## makeFootnote( footnote, size = .5 )
}

## geom_text_repel(data = filter( X, label == TRUE ),
