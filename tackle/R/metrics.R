library(dplyr)
library(ggplot2)
library(ggpubr)
## library(stringr)
library(tidyr)


metrics <- function(df, savename=NULL, exts=NULL, ...){


  options(bitmapType='cairo')
  df$Sample <- factor(rownames(df))


  df_sra<- df[,c('Sample', 'S', 'R', 'A')] %>%
    gather(S, R, A, key='SRA_metric', value='SRA')
  df_sra$SRA_metric <- factor(df_sra$SRA_metric, levels=c('S', 'R', 'A'),
                              ordered=TRUE)
  df_sra <- df_sra %>% group_by(Sample) %>% mutate(  cumsum = cumsum(SRA) )

  df_psms<- df[, c('Sample', 'Total_psms', 'u2g_psms')] %>%
    gather(Total_psms, u2g_psms, key='PSM_metric', value='PSMs')
  df_psms$PSM_metric <- factor(df_psms$PSM_metric, levels=c('Total_psms', 'u2g_psms'))


  df_pepts <- df[, c('Sample', 'Total_peptides', 'u2g_peptides', 'Strict', 'Strict_u2g')] %>%
    gather(Total_peptides, u2g_peptides, Strict, Strict_u2g, key='Peptide_metric', value='Peptides')
  df_pepts$Peptide_metric <- factor(df_pepts$Peptide_metric,
                                    levels=c('Total_peptides', 'u2g_peptides',
                                             'Strict', 'Strict_u2g'))

  ## thewidth <- dim(df)[1] %/% 2 %>% max(9)
  thewidth <- dim(df)[1] %/% 2 %>% max(9) %>% min(24)
  annot_scale = (12/thewidth)
  print(thewidth, annot_scale)
  overflow_width <- (24 / (dim(df)[1] %/% 2)) %>% min(1)
  print(dim(df)[1])
  print(overflow_width)

  ## green = 'darkgreen'; yellow = 'gold'; red ='firebrick'

  p0 <- ggplot(data=df_sra, aes(x=factor(Sample, levels=df$Sample),
                                y=SRA,
                                fill=factor(SRA_metric, levels=c('A', 'R', 'S')))) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c('firebrick', 'gold', 'darkgreen')) +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y=cumsum, label=SRA, color=SRA_metric), size=1.44*annot_scale,
              vjust=1.6, show.legend=FALSE) +
    scale_colour_manual(values=c("white", "black", 'white')) +
    theme_classic() +
    theme(text = element_text(size=14),
          legend.position = 'top',
          axis.text.x=element_blank(),
          legend.title=element_text(size=8),
          legend.text = element_text(size=8),
          ) +
    labs(fill='SRA\nMetric') +
  xlab(NULL)


  p1 <- ggplot(data=df, aes(x=factor(Sample, levels=df$Sample), y=GPGroups)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=GPGroups), vjust=1.6, color="white", size=1.24*annot_scale,
              show.legend=FALSE) +
    theme_classic() +
    theme(text = element_text(size=12),
          axis.text.x=element_blank(),
          ) +
    xlab(NULL)


  p2 <- ggplot(data=df_psms, aes(x=factor(Sample, levels=df$Sample), y=PSMs, fill=PSM_metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    theme(text = element_text(size=14),
          legend.position = 'top',
          legend.title=element_text(size=8),
          legend.text = element_text(size=8),
          axis.text.x = element_text(angle=90, hjust=1, vjust=.4, size = 12*overflow_width)
          ) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    labs(fill='PSM\nMetric') +
    xlab(NULL)


  p3 <- ggplot(data=df_pepts, aes(x=factor(Sample, levels=df$Sample), y=Peptides, fill=Peptide_metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    theme(text = element_text(size=14),
          legend.position = 'top',
          legend.title=element_text(size=8),
          legend.text = element_text(size=8),
          axis.text.x = element_text(angle=90, hjust=1, vjust=.4, size = 12*overflow_width)
          ) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    labs(fill='Peptide\nMetric') +
    xlab(NULL)


  pfinal <- ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend=FALSE, align='v', heights=c(1,1.40))

  if (!is.null(savename)){

    if (is.null(exts)){
      exts = c('png')
    }

    for (ext in exts){

      outname <- paste(savename, ext, sep='.')
      print(paste('Saving', outname))

      ggsave(outname, height=9, width=thewidth, units='in', limitsize = FALSE, dpi = 300)

    }

  }

}
