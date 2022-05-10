#install repo
devtools::install_github("SenguptaLab/MF.matR")

#open repo
library(MF.matR)
library(tidyverse)

#this script runs on a group of folders to obtain average residency plots for each condition/genotype
#place a text placeholder file in the file path directly above the subfolders that will be analyzed by the script. this tells the script to grab all the folders that are in the same folder as the placeholder file
plotAvgResidency(fillcolor = 'steelblue3', bordercolor = 'black')

#The above script will produce an average residency plot with an expanded y axis (starting at 0). The heatmaps axis starts from 1. Due to noisy data at the edges of the chip, only data from 1-15mm y-position is used. If you want the same axis for the heatmap and the average residency plot, adjust the y axis of the residency plot by inputting the following:
object <- plotAvgResidency()
object[[2]] + scale_y_continuous(limits = c(1,15), expand = c(0,0), breaks = c(1,5,10,15))