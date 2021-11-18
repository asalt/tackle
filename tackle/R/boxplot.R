
## library(plyr)
## library(dplyr)
library(data.table)
library(moments)
library(ggplot2)


boxplot <- function(df, group, title = "") {


  ## dfl <- df %>% filter(GeneSymbol == gene) %>%
  ##   select(-GeneID, -GeneSymbol) %>%
  ##   t
  ## colnames(dfl) <- 'Expression'
  ## dfl <- dfl %>% merge(select(meta, PDX, E2.dependence, Tumor, treat), by=0)

  df[is.na(df$Expression), "Expression"] <- 0



  p <- ggplot(data = df, aes_string(x = factor(group), y = "Expression", fill = group)) +
    geom_boxplot() +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, position = position_dodge(0.8), alpha = .5) +
    xlab(group) +
    ylab("log10 Expression") +
    ggtitle(title) +
    guides(fill = guide_legend(title = "")) +
    theme_light() +
    theme(
      text = element_text(size = 16),
      legend.position = "bottom",
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.x = element_text(size = 12, hjust = 0.5),
      axis.title.y = element_text(size = 14, vjust = 1),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      panel.grid.major.y = element_line(size = 0.3, linetype = "dotted", color = "darkgrey"),
      panel.grid.minor.y = element_line(size = 0.3, linetype = "dotted", color = "darkgrey"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank()
    )
  print(p)
  p
}
