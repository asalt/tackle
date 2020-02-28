## library(plyr)
library(dplyr)
library(forcats)
library(data.table)
library(moments)
library(ggplot2)


barplot <- function(df, average=FALSE, group=NA, group_order=NA, title = '', ylab='log10 Expression'){


  ## dfl <- df %>% filter(GeneSymbol == gene) %>%
  ##   select(-GeneID, -GeneSymbol) %>%
  ##   t
  ## colnames(dfl) <- 'Expression'
  ## dfl <- dfl %>% merge(select(meta, PDX, E2.dependence, Tumor, treat), by=0)


  df[is.na(df$Expression), 'Expression'] <- 0

  ## geom_dotplot(binaxis='y', stackdir='center',  dotsize=.5, position=position_dodge(0.8), alpha=.5) +


  if (!is.na(group) & average==TRUE){
    # the order is important here, sd before mean, because we're renaming Expression
    df <- df %>% group_by(get(group)) %>% summarize(sd=sd(Expression), Expression=mean(Expression))
    colnames(df)[1] <- 'index'
    group <- 'index' # we rename it here since that's now our averaged column
  }


  ## I'm bad at reordering things
  if (!is.na(group_order) & average==FALSE){
    df[[group]] <- factor(df[[group]], levels=group_order, ordered=TRUE)
    df <- df %>% mutate(index=fct_reorder(index, as.numeric(get(group)))
                        )
    df$index <- factor(df$index, levels=df$index, ordered=TRUE)
  } else if (!is.na(group_order) & average == TRUE){
    df[[group]] <- factor(df[[group]], levels=group_order, ordered=TRUE)
    df <- arrange(df, as.numeric(df$index))
  }


  ## p <- ggplot(data = df, aes_string(x=factor(group), y="Expression", fill=group)) +
  ## p <- ggplot(data = df, aes_string(x=fct_reorder('index', as.numeric(group_order)), y="Expression", fill=group)) +
  ## p <- ggplot(data = df, aes_string(x='index', y="Expression", fill=group)) +

  ## guides(fill=guide_legend(title=''), ) +


  ## p <- ggplot(data = df, aes(x=factor(index, level=df$index), y=Expression, fill=get(group))) +
  ## p <- ggplot(data = df, aes(x=index, y=Expression, fill=get(group))) +
  p <- ggplot(data = df, aes(x=index, y=Expression, fill=group)) +
    geom_bar(stat='identity') +
    xlab(NULL) + ylab(ylab) + ggtitle(title) +
    theme_light()+
    guides(fill=guide_legend(title='')) +
    theme(
      text = element_text(size=16),
      legend.position='bottom',
      plot.title = element_text(size=16, hjust = 0.5),
      axis.title.x = element_text(size=12,hjust=0.5),
      axis.title.y = element_text(size=14,vjust=1),
      axis.text.x = element_text(size=12,color='black', angle=90),
      axis.text.y = element_text(size=14, color='black'),
      panel.grid.major.y = element_line(size = 0.3, linetype = "dotted", color="darkgrey"),
      panel.grid.minor.y = element_line(size = 0.3, linetype = "dotted", color="darkgrey"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank()
      )


  if (!is.na(group_order)){
    p <- p + scale_fill_discrete(name='', breaks=group_order,
                                 labels=group_order)
      }

  if (average==TRUE){
    p <- p + geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                    position=position_dodge(.9))
  }

  print(p)
  p

}
