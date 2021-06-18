# these need to be installed within the R environment
# note for osx, xcode needs to be installed from the appstore AND
# xcode-select --install must  be run afterwards

if (!require("pacman")) install.packages("pacman")

library(pacman)
pacman::p_load("devtools", "tidyverse", "ggplot2", "ggrepel", "ggthemes", "statmod", "verification", "ggridges")


# ggplot2
# dplyr
# ggrepel
# ggthemes
# bioconductor
# statmod
# verification
# sva # via bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("sva", "limma"))

# believe these are requirements for ssGSEA
## "rhdf5", "prada")
## library(devtools)
## devtools::install_github("cmap/cmapR")
