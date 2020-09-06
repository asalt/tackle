library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ComplexHeatmap)
library(circlize) # colorRamp2

## This is copied directly from the end of ssgsea-gui.R
## then slightly modified to load and save files from/in appropriate directorys


if (!require("pacman")) install.packages ("pacman")
library(cmapR)
ssgsea_f <- file.path(dirname(sys.frame(1)$ofile),'..', 'ssGSEA2.0', 'src', 'ssGSEA2.0.R')
source(ssgsea_f)

library(gplots)

make_heatmap <- function(input_ds, geneset='hallmark', basedir='.'){

  input_ds <- '/mnt/e/projects/ELLIS-POL_TNBC/scripts/results/KIP_breastcancer_baseline/mednorm/ssGSEA/hallmark/input_dataset.gct'

  input_gct <- parse_gctx(input_ds)

  ## gene/site ids
  gn_input <- input_gct@rid
  ## if(dups)
  ##   gn.input <-  sub('_[0-9]{1,4}$', '', gn.input)

  ## sample names
  all_samp <- input_gct@cid

  ## expression data only
  input <- input_gct@mat
  ## load data
  ## ============================================================

  name <- dir(basedir, pattern=paste( geneset, '-scores.gct', sep=''))
  score_f <- file.path(basedir, name)

  gsea_score_gct <- parse_gctx(score_f)
  gsea_score <- gsea_score_gct@mat

  name <- dir(basedir, pattern=paste( geneset, '-pvalues.gct', sep=''))
  rawpval_f <- file.path(basedir, name)
  gsea_rawpval_gct <- parse_gctx(rawpval_f)
  gsea_rawpval <- gsea_rawpval_gct@mat

  name <- dir(basedir, pattern=paste( geneset, '-fdr-pvalues.gct', sep=''))
  pval_f <- file.path(basedir, name)
  gsea_pval_gct <- parse_gctx(pval_f)
  gsea_pval <- gsea_pval_gct@mat

  ## ============================================================

  gsea_pvall <- gsea_pval %>% as.data.frame %>% rownames_to_column('GeneSet') %>%
      pivot_longer(-GeneSet, names_to='id', values_to='FDR_pval', )

  gsea_rawpvall <- gsea_rawpval %>% as.data.frame %>% rownames_to_column('GeneSet') %>%
      pivot_longer(-GeneSet, names_to='id', values_to='pval', )

  gsea_long <- gsea_score %>% as.data.frame %>% rownames_to_column('GeneSet') %>%
      pivot_longer(-GeneSet, names_to='id', values_to='GSEA_Score', )  %>%
      left_join(gsea_pvall, by=c('GeneSet', 'id')) %>%
      left_join(gsea_rawpvall, by=c('GeneSet', 'id')) %>%
      mutate(log_FDRpVal = -log10(FDR_pval), log_pVal = -log10(pval))

  metrics <- gsea_long %>% select(-id) %>% group_by(GeneSet) %>% summarize_each(funs(mean, sd)) %>%
      arrange(log_FDRpVal_sd)

  ## filter by greater than 1st quartile of logFDRpVal standard deviation
  ## cutoffs <- summary(metrics$log_FDRpVal_sd)
  cutoffs <- summary(metrics$GSEA_Score_sd)
  ## mostdiff <- metrics %>% filter(GSEA_Score_sd > cutoffs[3]) %>% select(GeneSet)
  mostdiff <- metrics %>% select(GeneSet)

  ## ## gene set names
  ## all_gs <- rownames(gsea_score)

  toplot_l <- gsea_long
  toplot_l <- gsea_long %>% filter(GeneSet %in% mostdiff$GeneSet)
  toplot <- toplot_l %>% pivot_wider(id_cols='GeneSet', names_from='id', values_from='GSEA_Score')
  sigs_mat <- toplot_l %>% pivot_wider(id_cols='GeneSet', names_from='id', values_from='FDR_pval')


  annot <- apply(sigs_mat%>%select(-GeneSet), 1:2, function(x) ifelse( x < .05, '*', ''))

  col_meta <- gsea_score_gct@cdesc

  minval <- toplot %>% select(-GeneSet) %>% min(na.rm=TRUE) *.95
  maxval <- toplot %>% select(-GeneSet) %>% max(na.rm=TRUE) *.95
  ## col <- colorRamp2(c(minval, 0, maxval), c("#377EB8", "white", "#E41A1C"))
  col <- colorRamp2(c(minval,0,maxval), c('blue', "white", 'red'))
  print(minval)
  print(maxval)

  ht <- Heatmap(toplot %>% select(-GeneSet),
                ## row_split = cbind(kout_genes$cluster),
                ## column_split = cbind(col_meta$SampleCluster),
                col = col,
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                    if (annot[i, j] == '*') grid.text(annot[i, j], x, y)
                },
                ## top_annotation=hc,
                ## right_annotation=gene_annot,
                ## left_annotation = row_annot,
                ## top_annotation=col_annot, ##TODO: add this back
                column_title_rot = 0,
                cluster_columns=TRUE,
                cluster_rows=TRUE,
                cluster_column_slices = TRUE,
                show_row_names=TRUE,
                clustering_method_rows='ward.D2',
                clustering_method_columns='ward.D2',
                ## clustering_distance_rows=dist_no_na,
                ## clustering_distance_columns=dist_no_na,
                row_labels=toplot$GeneSet,
                row_names_gp = gpar(fontsize = 7.5),
                column_names_gp = gpar(fontsize = 7.5),
                border=FALSE,
                column_names_rot=90,
                name='hm',
                ## column_title = 'KIP PRM',
                column_title_side='top',
                column_names_side='top',
                show_parent_dend_line=TRUE,
                heatmap_legend_param=list(title='GSEA score')
                )
  width <- max(dim(toplot)[2] * .2, 6)
  height <- max(dim(toplot)[1] * .2, 8)
  pdf(file.path(basedir,  paste0(geneset, '.pdf')), width=width, height=height)
  print(ht)
  dev.off()

  ## ============================================================


}

make_plots <- function(input.ds, gene.set.database, basedir='.', geneset='hallmark', rank_plots = TRUE){
    ## #########################################################
    ##              rank plots
    ## #########################################################
    ## flag to indicate presence of duplicated ids in gct file
    ## e.g. in case of gene-centric-redundant signature analysis

    ## dups=f
    ## if(file.exists( sub('\\.gct', '_unique.gct', input.ds)))
    ##   dups <- t

    ## ## input dataset
    ## if(dups)
    ##    input.ds <-sub('\\.gct', '_unique.gct', input.ds)
    input.gct <- parse_gctx(input.ds)

    ## gene/site ids
    gn.input <- input.gct@rid
    ## if(dups)
    ##   gn.input <-  sub('_[0-9]{1,4}$', '', gn.input)

    ## sample names
    all.samp <- input.gct@cid

    ## expression data only
    input <- input.gct@mat


    ## import enrichment scores and p-values
    ## gsea.score.gct <- parse.gctx(dir('.', pattern=paste( i, '-scores(_[0-9]*x[0-9*]|)', '.gct', sep='')))
    ## gsea.score.gct <- parse.gctx(file.path(basedir, dir(basedir, pattern=paste( geneset, '-scores(_[0-9]*x[0-9*]|)',
    ##                                                                            '.gct', sep='')
    gsea.score.gct <- parse.gctx(file.path(basedir, dir(basedir, pattern=paste( geneset, '-scores.gct', sep=''))
                                           )
                                 )
    gsea.score <- gsea.score.gct@mat

    ## gsea.pval.gct <- parse.gctx(file.path(basedir, dir(basedir, pattern=paste( geneset,  '-fdr-pvalues(_[0-9]*x[0-9*]|)',
    ##                                                                           sep='')
    gsea.pval.gct <- parse.gctx(file.path(basedir, dir(basedir, pattern=paste( geneset,  '-fdr-pvalues.gct', sep=''))
                                          )
                                )
    gsea.pval <- gsea.pval.gct@mat

    ## gene set names
    all.gs <- rownames(gsea.score)

    ## plot heatmap of significant gene sets
    ## my_palette <- colorRampPalette(c('blue', 'white', 'red'))(n=299)
    annot <- apply(gsea.pval, 1:2, function(x) ifelse( x < .05, '*', ''))


    col_distance = dist(t(gsea.score), method = 'euclidean')
    col_cluster = hclust(col_distance, method = 'ward.D2')

    row_distance = dist(gsea.score, method = 'euclidean')
    row_cluster = hclust(row_distance, method = 'ward.D2')


    ## print(file.path(basedir,  paste0(geneset, '.pdf')))
    pdf(file.path(basedir,  paste0(geneset, '.pdf')), 8.5, 12)

    heatmap.2(gsea.score, trace='none', margins=c(9,17), col=bluered(99),
              cellnote = annot,
              density.info = 'none',
              Rowv = as.dendrogram(row_cluster),
              Colv = as.dendrogram(col_cluster),
              key.title = '', key.xlab = 'NES',
              na.color = 'grey', notecol='black'
              )
    text(0.08, y = 0.04, labels = '* p.adj < .05', adj = NULL,
         pos = NULL, offset = 0, vfont = NULL,
         cex = 1, col = NULL, font = NULL)

    dev.off()

    ## #############################################
    ## import signature database
    signat.all <- unlist(lapply(gene.set.database, readLines))
    signat.all <- strsplit(signat.all, '\t')
    names(signat.all) <- sapply(signat.all, function(x)x[1])
    signat.all <- lapply(signat.all, function(x) x[-c(1,2)])

    ## keep only scored signatures
    signat <- signat.all[all.gs]

    ## create sub-folder
    outdir <- file.path(basedir, 'rank-plots')
    dir.create(outdir)

    ## loop over gene sets
    for(gs in 1:length(all.gs)){

        gs.name <- all.gs[gs]

        ## pdf(paste('rank-plots/', make.names( chopstring( gsub('\\:|\\/\\\t', ' ', gs.name), nchar=20, add.dots=f)) ,'_2.pdf', sep=''), 9.5, 9.5)
        pdf(file.path(outdir, paste(make.names( chopString( gsub('\\:|\\/\\\t', ' ', gs.name), nChar=20, add.dots=FALSE)),
                                    '_2.pdf', sep='')), 9.5, 9.5)
        par(mfrow=c(3, 3))
        for(samp in 1:length(all.samp)){

            ## extract results
            samp.name <- all.samp[samp]

            ## gsea results
            score <- gsea.score[gs.name, samp.name]
            pval <- gsea.pval[gs.name, samp.name]

            ## extract data
            data.expr <- input[, samp.name ]

            valid.idx <- which( !(is.na( data.expr ) | is.infinite(data.expr)) )

            data.expr <- data.expr[ valid.idx ]
            gn <- gn.input[ valid.idx ]

            ## order
            ord.idx <- order(data.expr, decreasing=TRUE)

            ##gn <- row.names(input)[ord.idx]
            gn <- gn[ ord.idx ]
            data.expr <- data.expr[ ord.idx ]


            plot( data.expr, pch=20, col='darkgrey', lwd=4, type='l', xlab='rank', ylab='expression', main=paste(gs.name, samp.name, sep='\n'), ylim=range(data.expr), yaxs='i')
						abline(h=0, lty='dashed', lwd=2, col='grey70')

            ## #########################################################
            ##  ptm signatures?
            if(length(grep(';u$|;d$', signat[[gs.name]], value=t)) > 0){

                ## locations
                gsea.tmp.u <- sub(';u$','',grep(';u$', signat[[gs.name]], value=t))
                loc.u <- na.omit(match(gsea.tmp.u, gn))

                gsea.tmp.d <- sub(';d$','',grep(';d$',  signat[[gs.name]], value=t))
                loc.d <- na.omit(match(gsea.tmp.d, gn))

                if(!is.null(loc.u)){

                    rug(loc.u, col='darkred', side=3, lwd=3, ticksize=0.02)
                    points(loc.u, data.expr[loc.u], col=my.col2rgb('darkred',  150), pch=16, cex=2)
                }
                if(!is.null(loc.d)){
                    rug(loc.d, col='darkblue', side=1, lwd=3, ticksize=0.02)
                    points(loc.d, data.expr[loc.d], col=my.col2rgb('darkblue',  150), pch=16, cex=2)

                }
                ## some info
                legend('bottom', legend=c(paste('no. down-regulated in signature:', length(grep(';d$', signat[[gs.name]]))),
                                          paste('no. found in data set:', length(loc.d))
                                          ), inset=.05, bty='n', text.col='darkblue')

                legend('top', legend=c(paste('no. up-regulated in signature:', length(grep(';u$', signat[[gs.name]]))),
                                       paste('no. found in data set:', length(loc.u))
                                       ), inset=.05, bty='n', text.col='darkred')
            } else {## end if signature

                ## ####################################################
                ## regular gene set
                loc <- which(gn %in% signat[[gs.name]])
                rug(loc, col=my.col2rgb('darkred',  50), side=3, lwd=2, ticksize=0.02)
								points(loc, data.expr[loc], col=my.col2rgb('darkred',  150), pch=16, cex=2)

                ## box plot
                loc.quart <- quantile(loc)
                rug(loc.quart, col='darkblue', side=3, lwd=2, ticksize=0.03)
                rect( loc.quart[2], max(data.expr)-0.04*max(data.expr-min(data.expr)), loc.quart[4], max(data.expr), border='darkblue', lwd=2, col=NA )
                rect( loc.quart[1], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[2], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )
                rect( loc.quart[4], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[5], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )

                ## some info
                legend('bottom', legend=c(paste('no. in signature:', length( signat[[gs.name]])),
                                          paste('no. found in data set (non-redund.):', sum(signat[[gs.name]] %in% gn)),
                                          paste('no. found in data set (redundant):', length(loc))
                                       ), inset=.05, bty='n', text.col='darkred')

            }

          legend('right',
                 legend=paste('nes=', round(score, 3), ' (p.adj=', round(pval, 5), ')', sep=''),
                 bty='n', inset=.22, cex=1.25)

	}
	par(mfrow=c(1, 1))
	dev.off()
    } ## end loop over gene sets
}
