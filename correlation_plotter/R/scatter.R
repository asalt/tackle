library(graphics)

#' adapted from PerformanceAnalytics::chart.Correlation
scattermat <- function(x, histogram = FALSE, interactive=TRUE,
                       method = c("pearson", "kendall",
                                  "spearman"), ...)
{

  nrows <- dim(x)[2]
  plt_size <- nrows * .75
  if (interactive){
    dev.new(width = plt_size, height = plt_size)
  }

  ii <- seq(-1, 1, .10)
  colors = colorRampPalette(c('blue', 'white', 'red'))(length(ii))


  ## x = checkData(R, method = "matrix")
  if (missing(method))
      method = method[1]
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs",
      method, cex.cor, ...) {
      usr <- par("usr")
      on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use = use, method = method)
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste(prefix, txt, sep = "")

      overlap_length <- dim( na.omit(matrix(c(x,y), ncol=2)) )[1]

      txt <- paste(txt, paste('n=', overlap_length, sep=''), sep='\n')

      if (missing(cex.cor))
          cex <- 0.7/strwidth(txt)
      ## test <- cor.test(x, y, method = method)
      ## Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
      ##     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***",
      ##         "**", "*", ".", " "))
      ## text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
      bgcolor <- colors[findInterval(r, ii)]
      u <- par('usr')
      names(u) <- c("xleft", "xright", "ybottom", "ytop")
      do.call(rect, c(col = bgcolor, as.list(u)))
      text( 0.5, 0.5, txt, cex = cex )
      text( 0.5, 0.5, txt, cex = cex )

      }
  f <- function(t) {
      dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  hist.panel = function(x, ...) {
      par(new = TRUE)
      hist(x, col = "light gray", probability = TRUE, axes = FALSE,
          main = "", breaks = "FD")
      lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
      ## rug(x)
  }
  if (histogram)
      pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
            diag.panel = hist.panel, method = method, pch=16, col='#22222288',
            ## cex.labels = colnames(x),
            ...)
  ## else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
  else pairs(x, gap = 0, smooth = FALSE, upper.panel = panel.cor,
              method = method, pch=16, col='#22222288',
              ## cex.labels = colnames(x),
              ...)
}

## scattermat(edata)
