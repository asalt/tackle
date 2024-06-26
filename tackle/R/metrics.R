suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))



metrics <- function(df, savename = NULL, exts = NULL, ...) {
  options(bitmapType = "cairo")
  df[["Sample"]] <- factor(rownames(df))


  df_sra <- df[, c("Sample", "S", "R", "A")] %>%
    gather(S, R, A, key = "SRA_metric", value = "SRA")
  df_sra$SRA_metric <- factor(df_sra$SRA_metric,
    levels = c("S", "R", "A"),
    ordered = TRUE
  )
  df_sra <- df_sra %>%
    group_by(Sample) %>%
    mutate(cumsum = cumsum(SRA))

  df_psms <- df[, c("Sample", "Total_psms", "u2g_psms")] %>%
    gather(Total_psms, u2g_psms, key = "PSM_metric", value = "PSMs")
  df_psms$PSM_metric <- factor(df_psms$PSM_metric, levels = c("Total_psms", "u2g_psms"))


  df_pepts <- df[, c("Sample", "Total_peptides", "u2g_peptides", "Strict", "Strict_u2g")] %>%
    gather(Total_peptides, u2g_peptides, Strict, Strict_u2g, key = "Peptide_metric", value = "Peptides")
  df_pepts$Peptide_metric <- factor(df_pepts$Peptide_metric,
    levels = c(
      "Total_peptides", "u2g_peptides",
      "Strict", "Strict_u2g"
    )
  )

  ## thewidth <- dim(df)[1] %/% 2 %>% max(9)
  thewidth <- dim(df)[1] %/% 2 %>%
    max(9) %>%
    min(24)
  annot_scale <- (15 / thewidth)
  # print(thewidth, annot_scale)
  overflow_width <- (24 / (dim(df)[1] %/% 2)) %>% min(1)
  legend_text_size <- 10
  axis_text_size <- 18
  # print(dim(df)[1])
  # print(overflow_width)

  ## green = 'darkgreen'; yellow = 'gold'; red ='firebrick'

  p0 <- ggplot(data = df_sra, aes(
    x = factor(Sample, levels = df$Sample),
    y = SRA,
    fill = factor(SRA_metric, levels = c("A", "R", "S"))
  )) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("firebrick", "gold", "darkgreen")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = cumsum, label = SRA, color = SRA_metric),
      size = 1.44 * annot_scale,
      vjust = 1.3, show.legend = FALSE
    ) +
    scale_colour_manual(values = c("white", "black", "white")) +
    theme_classic() +
    theme(
      text = element_text(size = axis_text_size + 2),
      legend.position = "top",
      # axis.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4, size = axis_text_size * overflow_width),
      axis.text.y = element_text(size = axis_text_size),
      legend.title = element_text(size = legend_text_size),
      legend.text = element_text(size = legend_text_size),
    ) +
    labs(fill = "SRA\nMetric", y = "Gene Products") +
    xlab(NULL)
  # coord_flip()


  p1 <- ggplot(data = df, aes(x = factor(Sample, levels = df$Sample), y = GPGroups)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = GPGroups),
      vjust = 1.4, color = "white", size = 1.44 * annot_scale,
      show.legend = FALSE
    ) +
    theme_classic() +
    theme(
      text = element_text(size = axis_text_size + 2),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4, size = axis_text_size * overflow_width),
      axis.text.y = element_text(size = axis_text_size),
      # axis.text.x = element_blank(),
    ) +
    xlab(NULL)
  # coord_flip()


  p2 <- ggplot(data = df_psms, aes(x = factor(Sample, levels = df$Sample), y = PSMs, fill = PSM_metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_classic() +
    theme(
      text = element_text(size = axis_text_size + 2),
      legend.position = "top",
      legend.title = element_text(size = legend_text_size + 2),
      legend.text = element_text(size = legend_text_size),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4, size = axis_text_size * overflow_width),
      axis.text.y = element_text(size = axis_text_size),
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(fill = "PSM\nMetric") +
    xlab(NULL)
  # coord_flip()


  p3 <- ggplot(data = df_pepts, aes(x = factor(Sample, levels = df$Sample), y = Peptides, fill = Peptide_metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_classic() +
    theme(
      text = element_text(size = axis_text_size + 2),
      legend.position = "top",
      legend.title = element_text(size = legend_text_size),
      legend.text = element_text(size = legend_text_size),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4, size = axis_text_size * overflow_width)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(fill = "Peptide\nMetric") +
    xlab(NULL)
  # coord_flip()


  # pfinal <- ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend=FALSE, align='v', heights=c(1,1.40))

  if (!is.null(savename)) {
    if (is.null(exts)) {
      exts <- c("png")
    }

    # plots <- c(p0, p1, p2, p3)
    # plot_names <- c('SRA', 'GPGroups', 'PSMs', 'Peptides')
    # purrr::map2(plots, plot_names, ~ save)
    make_outname <- function(x, ext) {
      paste(savename, x, sep = "_") %>%
        paste(ext, sep = ".")
    }

    for (ext in exts) {
      ggsave(make_outname("SRA", ext),
        plot = print(p0), height = 9, width = thewidth, units = "in", limitsize = FALSE, dpi = 300
      )

      ggsave(make_outname("GPGroups", ext), plot = print(p1), height = 9, width = thewidth, units = "in", limitsize = FALSE, dpi = 300)

      ggsave(make_outname("PSMs", ext),
        plot = p2, height = 9, width = thewidth, units = "in", limitsize = FALSE, dpi = 300
      )

      ggsave(make_outname("Peptides", ext),
        plot = p3, height = 9, width = thewidth, units = "in", limitsize = FALSE, dpi = 300
      )

      # ggsave(outname, height=9, width=thewidth, units='in', limitsize = FALSE, dpi = 300)
    }
  }
}
