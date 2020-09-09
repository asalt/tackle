install.packages("statmod")
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("ggridges")
install.packages(c("ComplexHeatmap", "circlize"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }

BiocManager::install("limma")
BiocManager::install("sva")
