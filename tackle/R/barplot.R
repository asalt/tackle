## library(plyr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(moments))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(colorspace))
# suppressPackageStartupMessages(library(paletteer))

barplot <- function(df, average = NULL, group = NULL, group_order = NULL, title = "", ylab = "log10 Expression", cmap=NULL) {
  # Check if 'Expression' column exists
  if(!"Expression" %in% colnames(df)) {
    stop("The dataframe must contain an 'Expression' column.")
  }

  # Replace NA in Expression with 0
  df <- df %>%
    mutate(Expression = ifelse(is.na(Expression), 0, Expression))

# Handle grouping and averaging
  if (!is.null(group)) {
    
    # Check if the grouping column exists
    if(!group %in% colnames(df)) {
      stop(paste("The grouping column", group, "does not exist in the dataframe."))
    }
    
    # If averaging is requested
    if (!is.null(average)) {
      df <- df %>%
        group_by({{ group }}) %>%
        summarize(
          Expression = mean(Expression, na.rm = TRUE),
          sd = sd(Expression, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        rename(index = {{ group }})
      group <- "index"  # Update group to the new summarized column
      
    } else {
      # If group_order is provided, set factor levels accordingly
      if (!is.null(group_order)) {
        df <- df %>%
          mutate({{ group }} := factor(!!sym(group), levels = group_order , ordered = TRUE))
      } else {
        # If no group_order, set the group as an ordered factor based on appearance
        df <- df %>%
          mutate({{ group }} := as.factor(!!sym(group)))
      }
      
      # Optional: Create an 'index' for plotting if needed
      # df <- df %>%
      #   arrange({{ group }}) %>%
      #   mutate(index = factor({{ group }}, levels = levels({{ group }}), ordered = TRUE))
    }
  } else {
    # If no grouping is specified, ensure 'index' exists for plotting
    # if(!"index" %in% colnames(df)) {
    #   df <- df %>%
    #     mutate(index = factor(1:n()))
    # }
  }


  ## geom_dotplot(binaxis='y', stackdir='center',  dotsize=.5, position=position_dodge(0.8), alpha=.5) +

  ## TODO fix when group is null
  ## I'm bad at reordering things
  # if (!is.null(group) && !is.null(average)) {
  #   # the order is important here, sd before mean, because we're renaming Expression
  #   df <- df %>%
  #     group_by(get(group)) %>%
  #     summarize(sd = sd(Expression), Expression = mean(Expression))
  #   colnames(df)[1] <- "index"
  #   group <- "index" # we rename it here since that's now our averaged column
  # }
  # else if (!is.null(group) && !is.null(group_order) && is.null(average)) {
  #   df[[group]] <- factor(df[[group]], levels = group_order, ordered = TRUE)
  #   ## df <- df %>% mutate(index=fct_reorder(index, as.numeric(get(group))))
  #   df <- arrange(df, get(group))
  #   df$index <- factor(df$index, levels = df$index, ordered = TRUE)
  # }
  #  else if (!is.null(group_order) && !is.null(average)) {
  #   df[[group]] <- factor(df[[group]], levels = group_order, ordered = TRUE)
  #   df <- arrange(df, as.numeric(df$index))
  # } else if (!is.null(group) && is.null(group_order) && is.null(average)) {
  #   # automatic ordering
  #   df[[group]] <- factor(df[[group]], ordered = TRUE)
  #   df <- df %>% arrange(!!!group)
  #   ## TODO: fix this
  #   ## df[['index']] <- factor(df[['index']], ordered = TRUE)
  #   ## df <- arrange(df, as.numeric(df$index))
  #   ## ??
  #   df <- df %>% mutate(index = factor(index, ordered = TRUE))
  # }

    # scale_fill_paletteer_d("ggthemes::calc") +
  ## scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +

  p <- ggplot(data = df, aes_string(x = "index", y = "Expression", fill = group)) +
    geom_bar(stat = "identity") +
    xlab(NULL) +
    ylab(ylab) +
    ggtitle(title) +
    guides(fill = guide_legend(title = "")) +
    labs(title = title, reverse = TRUE) +
    ylab(ylab) +
    theme_light() +
    colorspace::scale_fill_discrete_qualitative(palette = "Dynamic") +
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

  if (!is.null(average)) {
    p <- p + geom_errorbar(aes(ymin = Expression - sd, ymax = Expression + sd),
      width = .2,
      position = position_dodge(.9)
    )
  }

  print(p)
  p
}
