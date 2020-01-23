## library(plyr)
library(dplyr)
library(forcats)
library(data.table)
library(moments)
library(ggplot2)


barplot <- function(df, color=NA, group=NA, group_order=NA, title = '', ylab='log10 Expression'){


  ## dfl <- df %>% filter(GeneSymbol == gene) %>%
  ##   select(-GeneID, -GeneSymbol) %>%
  ##   t
  ## colnames(dfl) <- 'Expression'
  ## dfl <- dfl %>% merge(select(meta, PDX, E2.dependence, Tumor, treat), by=0)

  df[is.na(df$Expression), 'Expression'] <- 0

  ## geom_dotplot(binaxis='y', stackdir='center',  dotsize=.5, position=position_dodge(0.8), alpha=.5) +


  ## if (!is.na(group)){
  ##   df <- df %>% group_by(get(group)) %>% summarize(mean=mean(Expression), sd=sd(Expression))
  ## } else{
  ##   group <- ''
  ## }

  if (!is.na(group_order)){
    df[[group]] <- factor(df[[group]], levels=group_order)
    df <- df %>% mutate(index=fct_reorder(index, as.numeric(get(group)))
                        )
    df$index <- factor(df$index, levels=df$index, ordered=TRUE)
  }


  ## browser()

  ## p <- ggplot(data = df, aes_string(x=factor(group), y="Expression", fill=group)) +
  ## p <- ggplot(data = df, aes_string(x=fct_reorder('index', as.numeric(group_order)), y="Expression", fill=group)) +
  ## p <- ggplot(data = df, aes_string(x='index', y="Expression", fill=group)) +
  p <- ggplot(data = df, aes(x=factor(index, level=df$index), y=Expression, fill=get(group))) +
    geom_bar(stat='identity') +
    xlab(group) + ylab(ylab) + ggtitle(title) +
    guides(fill=guide_legend(title='')) +
    theme_light()+
    theme(
      text = element_text(size=16),
      legend.position='bottom',
      plot.title = element_text(size=16, hjust = 0.5),
      axis.title.x = element_text(size=12,hjust=0.5),
      axis.title.y = element_text(size=14,vjust=1),
      axis.text.x = element_text(size=12,color='black'),
      axis.text.y = element_text(size=14, color='black'),
      panel.grid.major.y = element_line(size = 0.3, linetype = "dotted", color="darkgrey"),
      panel.grid.minor.y = element_line(size = 0.3, linetype = "dotted", color="darkgrey"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank()
      )
  print(p)
  p

}
