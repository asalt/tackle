# Load packages
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))

## library(ggplot2)
library(graphics)

# Install ggrepel package if needed
# install.packages("ggrepel")

makeFootnote <- function(footnoteText = format(Sys.time(), "%d %b %Y"),
                        size = 1.0, color = grey(.5))
{
  suppressMessages(
    require(grid)
  )
  pushViewport(viewport())
  grid.text(label = footnoteText ,
            x = unit(1,"npc") - unit(2, "mm"),
            y = unit(2, "mm"),
            just = c("right", "bottom"),
            gp = gpar(cex = size, col = color))
  popViewport()
}
yaxis.choices <- c('pValue', 'pAdj')
number_by.choices <- c('abs_log2_Fold_Change', 'log2_Fold_Change', 'pValue')

# :fc_cutoff: cutoff for absolute fold change cutoff
volcanoplot <- function(X, max_labels = 35,
                        pch = 16, cex = 0.35,
                        fc_cutoff = 4, sig = 0.05, label_cex = 1,
                        show_all = FALSE, yaxis = yaxis.choices,
                        group0 = '', group1 = '',
                        sig_metric = 'pAdj',
                        number_by = 'abs_log2_Fold_Change',
                        ...){

  ploty <- match.arg(yaxis, yaxis.choices)
  number_by <- match.arg(number_by, number_by.choices)
  linear_fc_cutoff <- fc_cutoff
  fc_cutoff <- log2(fc_cutoff)

  if (sig_metric == 'pAdj') {
    sig_filter_str <- paste0('FDR<', sig)
  }
  else {
    sig_filter_str <- paste0('p<', sig)
  }

  ## X <- mutate(X, Sig = ifelse(X$pAdj < sig& abs(X[, 'log2_Fold_Change']) > fc_cutoff,
  ##                             sig_filter_str, "N.S."))
  ## library(tibble)
  ## browser()
  ## X <- as.tibble(X)
  ## X <- mutate(X, Sig = ifelse(X[,sig_metric] < sig& abs(X[, 'log2_Fold_Change']) > fc_cutoff,
  ##                             sig_filter_str, "N.S."))

  Sig <- ifelse(X[,sig_metric] < sig& abs(X[, 'log2_Fold_Change']) > fc_cutoff,
               sig_filter_str, "N.S.")
  X[, 'Sig'] <- Sig


  ## X[ , 'usd' ] = '#222222bb'
  ## X[ , 'usd' ] = '#22222222'
  X[ , 'usd' ] = '#22222288'
  ## X[ (X$pAdj < sig & X$log2_Fold_Change > fc_cutoff), 'usd' ] = 'red'
  ## X[ (X$pAdj < sig & X$log2_Fold_Change < -fc_cutoff), 'usd' ] = 'blue'
  X[ (X[,sig_metric] < sig & X$log2_Fold_Change < -fc_cutoff), 'usd' ] = 'blue'
  X[ (X[,sig_metric] < sig & X$log2_Fold_Change > fc_cutoff), 'usd' ] = 'red'
  X[ X$highlight == TRUE, 'usd' ] = 'purple'
  X[, 'usd'] <- as.factor(X[, 'usd'])

  ## X <- mutate(X, label = ifelse(X$qValue < 0.05, "FDR<0.05", "N.S."))

  pAdj <- X[, 'pAdj'][ !is.na( X[, 'pAdj'] ) ]
  stretch <- min( pAdj[ pAdj > 0 ] ) / 2
  X[, 'pAdj'] <- X[, 'pAdj'] + stretch
  pvalues <- X[, 'pValue'][ !is.na( X[, 'pValue'] ) ]
  stretch <- min( pvalues[ pvalues > 0 ] ) / 2
  X[, 'pValue'] <- X[, 'pValue'] + stretch

  ## to_label <- head(order( abs(X[, 'log2_Fold_Change']), X[, 'pAdj'], decreasing = c(TRUE, FALSE) ),
  ##                  max_labels
  ##                  )

  if (number_by == 'abs_log2_Fold_Change'){
    to_label <- head(order( abs(X[, 'log2_Fold_Change']), X[, 'pValue'], decreasing = c(TRUE, FALSE) ),
                     max_labels
                     )
  }
  else if (number_by == 'log2_Fold_Change'){

    to_label1 <- head(order( (X[, 'log2_Fold_Change']), X[, 'pValue'], decreasing = c(TRUE, FALSE) ),
                      round(max_labels/2)
                      )
    to_label2 <- head(order( (X[, 'log2_Fold_Change']), X[, 'pValue'], decreasing = c(FALSE, FALSE) ),
                      round(max_labels/2)
                      )
    to_label <- c(to_label1, to_label2)

  }
  else{
    to_label <- head(order( X[, 'pValue'], abs(X[, 'log2_Fold_Change']), decreasing = c(FALSE, TRUE) ),
                     max_labels
                     )
  }


  X[, 'label'] <- FALSE  # new column
  X[to_label, 'label'] <- TRUE  # label these

  if (show_all == FALSE){
    X[ X[, 'Sig'] == 'N.S.', 'label'] <- FALSE
  }

  X[ X$highlight == TRUE, 'label' ] <- TRUE # also label these no matter what

  ## ymax <- max(-log10(X[, 'pValue'])) * 1.05
  ymax <- max(-log10(X[, ploty])) * 1.05
  xmax <- X[, 'log2_Fold_Change'] %>% abs %>% max

  ## ratio_sig <- paste0( dim( filter(X, Sig == sig_filter_str) )[1], '/', dim(X)[1] )
  ratio_sig <- paste0( dim( X[X$Sig == sig_filter_str,] )[1], '/', dim(X)[1] )
  footnote <- paste( ratio_sig, 'sig. at', sig_filter_str, 'and',  linear_fc_cutoff, 'F.C.' )
  ylabel_full <- eval(expression(substitute(paste('-log'[10],' ', ploty), list(ploty=ploty))))

  p = ggplot(X, aes(log2_Fold_Change, -log10(get(ploty)), col=usd)) +
    theme_base() +
    geom_point(size = 1, cex = cex, show.legend = FALSE, pch=pch) +
    scale_colour_identity() +
    geom_text_repel(data = X[X$label==TRUE, ],
                    aes(label = GeneSymbol),  min.segment.length = .05,
                    point.padding = 1e-6,
                    box.padding = .1, cex = label_cex,
                    segment.size = .35, segment.alpha = .4
                    ) +
    annotate("text",  c(-xmax, xmax), c(ymax*.98, ymax*.98), label = c(group0, group1),
             hjust = c(0, 1), vjust = c(0,0), color = c('blue', 'red')) +
    labs(x = expression(paste('log'[2], ' Fold Change')),
         y=ylabel_full,
         caption = footnote) +
    theme(plot.caption = element_text(color = grey(.5), size=10)) +
    ylim(0, ymax)

  print(p)

  ## print(X %>% filter(Sig=='FDR<0.05') %>% dim)
  ## print(X %>% dim)

  ## footnote <- paste( ratio_sig, 'sig. with', 2**fc_cutoff, 'F.C.' )
  ## makeFootnote( footnote, size = .5 )


}

    ## geom_text_repel(data = filter( X, label == TRUE ),
