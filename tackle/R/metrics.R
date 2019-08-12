library(dplyr)
library(ggplot2)
library(ggpubr)
## library(stringr)
library(tidyr)


metrics <- function(df, savename=NULL, exts=NULL, ...){


  df$Sample <- factor(rownames(df))


  df_sra<- df[,c('Sample', 'S', 'R', 'A')] %>%
    gather(S, R, A, key='SRA_metric', value='SRA')
  df_sra$SRA_metric <- factor(df_sra$SRA_metric, levels=c('A', 'R', 'S'))

  df_psms<- df[, c('Sample', 'Total_psms', 'u2g_psms')] %>%
    gather(Total_psms, u2g_psms, key='PSM_metric', value='PSMs')
  df_psms$PSM_metric <- factor(df_psms$PSM_metric, levels=c('Total_psms', 'u2g_psms'))


  df_pepts <- df[, c('Sample', 'Total_peptides', 'u2g_peptides', 'Strict', 'Strict_u2g')] %>%
    gather(Total_peptides, u2g_peptides, Strict, Strict_u2g, key='Peptide_metric', value='Peptides')
  df_pepts$Peptide_metric <- factor(df_pepts$Peptide_metric,
                                    levels=c('Total_peptides', 'u2g_peptides',
                                             'Strict', 'Strict_u2g'))


  ## green = 'darkgreen'; yellow = 'gold'; red ='firebrick'

  p0 <- ggplot(data=df_sra, aes(x=factor(Sample, levels=df$Sample), y=SRA, fill=SRA_metric)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c('firebrick', 'gold', 'darkgreen')) +
    theme(text = element_text(size=14),
          legend.position = 'top',
          axis.text.x=element_blank(),
          legend.title=element_text(size=8),
          legend.text = element_text(size=8),
          ) +
  xlab(NULL)


  p1 <- ggplot(data=df, aes(x=factor(Sample, levels=df$Sample), y=GPGroups)) +
    geom_bar(stat="identity") +
    theme(text = element_text(size=12),
          axis.text.x=element_blank(),
          ) +
    xlab(NULL)


  p2 <- ggplot(data=df_psms, aes(x=factor(Sample, levels=df$Sample), y=PSMs, fill=PSM_metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(text = element_text(size=14),
          legend.position = 'top',
          legend.title=element_text(size=8),
          legend.text = element_text(size=8),
          axis.text.x = element_text(angle=90, hjust=1, vjust=.4)
          ) +
    xlab(NULL)


  p3 <- ggplot(data=df_pepts, aes(x=factor(Sample, levels=df$Sample), y=Peptides, fill=Peptide_metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(text = element_text(size=14),
          legend.position = 'top',
          legend.title=element_text(size=8),
          legend.text = element_text(size=8),
          axis.text.x = element_text(angle=90, hjust=1, vjust=.4)
          ) +
    xlab(NULL)


  pfinal <- ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend=FALSE, align='v', heights=c(1,1.40))

  if (!is.null(savename)){

    if (is.null(exts)){
      exts = c('png')
    }

    for (ext in exts){

      outname <- paste(savename, ext, sep='.')

      thewidth <- dim(df)[1] %/% 3 %>% max(5) %>% min(16)
      print(thewidth)
      ggsave(outname, height=9, width=thewidth, units='in')

    }

  }

}
