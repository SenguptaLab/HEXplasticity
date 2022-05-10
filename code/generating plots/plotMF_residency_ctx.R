#install repo
devtools::install_github("SenguptaLab/MF.matR")

#open repo
library(MF.matR)
library(tidyverse)

#will prompt to select any csv file in individual folder 
#this script will produce a heatmap, a residency plot with overlaid luminance, a chemotaxis index value csv file
plotResidency_stripes(FileFilter = "", heatmap_limits = c(0,2.5,5), heatmap_palette = "Blues")
