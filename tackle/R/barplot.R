## library(plyr)
library(dplyr)
library(forcats)
library(data.table)
library(moments)
library(ggplot2)
library(ggpubr)
library(paletteer)

barplot <- function(df, average = FALSE, group = NULL, group_order = NULL, title = "", ylab = "log10 Expression") {


  ## dfl <- df %>% filter(GeneSymbol == gene) %>%
  ##   select(-GeneID, -GeneSymbol) %>%
  ##   t
  ## colnames(dfl) <- 'Expression'
  ## dfl <- dfl %>% merge(select(meta, PDX, E2.dependence, Tumor, treat), by=0)

  ## browser()

  df[is.na(df$Expression), "Expression"] <- 0

  ## geom_dotplot(binaxis='y', stackdir='center',  dotsize=.5, position=position_dodge(0.8), alpha=.5) +

  ## if (!is.na(group)) df[[group]] <- make.names(df[[group]]) %>% gsub('\\.', '', . )
  ## if (!is.na(group_order)) group_order <- make.names(group_order) %>% gsub("\\.", "", .)
  # browser()

  if (!is.null(group) & average == TRUE) {
    # the order is important here, sd before mean, because we're renaming Expression
    df <- df %>%
      group_by(get(group)) %>%
      summarize(sd = sd(Expression), Expression = mean(Expression))
    colnames(df)[1] <- "index"
    group <- "index" # we rename it here since that's now our averaged column
  }
  ## TODO fix when group is null


  ## I'm bad at reordering things
  if (!is.null(group_order) & average == FALSE) {
    df[[group]] <- factor(df[[group]], levels = group_order, ordered = TRUE)
    ## df <- df %>% mutate(index=fct_reorder(index, as.numeric(get(group))))
    df <- arrange(df, get(group))
    df$index <- factor(df$index, levels = df$index, ordered = TRUE)
  } else if (!is.null(group_order) & average == TRUE) {
    df[[group]] <- factor(df[[group]], levels = group_order, ordered = TRUE)
    df <- arrange(df, as.numeric(df$index))
  } else if (is.null(group_order) & average == FALSE) {
    # automatic ordering
    df[[group]] <- factor(df[[group]], ordered = TRUE)
    df <- df %>% arrange(!!!group)
    ## TODO: fix this
    ## df[['index']] <- factor(df[['index']], ordered = TRUE)
    ## df <- arrange(df, as.numeric(df$index))
    ## ??
    df <- df %>% mutate(index = factor(index, ordered = TRUE))
  }


  ## p <- ggplot(data = df, aes_string(x=factor(group), y="Expression", fill=group)) +
  ## p <- ggplot(data = df, aes_string(x=fct_reorder('index', as.numeric(group_order)), y="Expression", fill=group)) +
  ## p <- ggplot(data = df, aes_string(x='index', y="Expression", fill=group)) +

  ## guides(fill=guide_legend(title=''), ) +


  ## p <- ggplot(data = df, aes(x=factor(index, level=df$index), y=Expression, fill=get(group))) +
  ## p <- ggplot(data = df, aes(x=index, y=Expression, fill=get(group))) +

  ## scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +

  p <- ggplot(data = df, aes_string(x = "index", y = "Expression", fill = group)) +
    geom_bar(stat = "identity") +
    xlab(NULL) +
    ylab(ylab) +
    ggtitle(title) +
    guides(fill = guide_legend(title = "")) +
    labs(title = title, reverse = TRUE) +
    ylab(ylab) +
    scale_fill_paletteer_d("ggthemes::calc") +
    theme_light() +
    theme(
      text = element_text(size = 16),
      legend.position = "bottom",
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.x = element_text(size = 12, hjust = 0.5),
      axis.title.y = element_text(size = 14, vjust = 1),
      axis.text.x = element_text(size = 12, color = "black", angle = 90),
      axis.text.y = element_text(size = 14, color = "black"),
      panel.grid.major.y = element_line(size = 0.3, linetype = "dotted", color = "darkgrey"),
      panel.grid.minor.y = element_line(size = 0.3, linetype = "dotted", color = "darkgrey"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank()
    )
  # theme_light() +
  #+
  # ggpubr::theme_pubr() +
  # paletteer_d("ggthemes::Classic_10")


  if (!is.null(group_order)) {
    p2 <- p + scale_fill_discrete(
      name = "", breaks = group_order,
      labels = group_order
    )
  }

  if (average == TRUE) {
    p <- p + geom_errorbar(aes(ymin = Expression - sd, ymax = Expression + sd),
      width = .2,
      position = position_dodge(.9)
    )
  }

  print(p)
  p
}
