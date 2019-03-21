#' Same calculation as sva::f.pvalue but also accounting for
#' batch variance removal (such as from ComBat)
#' This reduces the degrees of freedom by the number of covariates.
#' see https://support.bioconductor.org/p/63007/#63115
pvalue.batch <- function (dat, mod, mod0, ncov=0)
{
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
                    t(mod))
  rss1 <- rowSums(resid * resid)
  rm(resid)
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*%
                     t(mod0))
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)
  fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1 - ncov))
  p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1 - ncov))
  return(p)
}
