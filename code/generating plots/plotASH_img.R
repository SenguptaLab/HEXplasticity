#install repo
devtools::install_github("SenguptaLab/MF.matR")

#open repo
library(MF.matR)
library(tidyverse)

#will prompt to select folder with imaging data
#will generate raw data csv file
#plotting traces for same genotype, different condition
ash_hex<-plotGCaMP_multi(cue = HEX, 
                         genotype = WT,
                         neuron = ASH, 
                         center_on_pulse = 'ON',
                         show.plots = FALSE,
                         backsub = TRUE,
                         heatmap_limits = c(-0.1,0,1.2))

ash_hex$data 
ash_hex$plot

#will prompt to select folder with imaging data
#will generate raw data csv file
#plotting traces for same genotype, different condition
ash_hex_iaa<-plotGCaMP_multi(cue = HEX-IAA, 
                         genotype = HEX-IAA,
                         neuron = AWC, 
                         center_on_pulse = 'ON',
                         show.plots = FALSE,
                         backsub = TRUE,
                         heatmap_limits = c(-0.1,0,1.2))

ash_hex_iaa$data 
ash_hex_iaa$plot

#rbind by condition
alldata <- rbind(ash_hex$data, ash_hex_iaa$data)

#unnest data to analyze for mean & SEM
unnested <- unnest(alldata)

# get mean and SEM for all genotypes at once
alldata.mean <- unnested %>%
  unnest() %>%
  group_by(time, genotype) %>%
  summarise(meanDelF = mean(delF),
            n = n(),
            sem.low = meanDelF - sd(delF) / n()^0.5, #this is for standard error
            sem.high = meanDelF + sd(delF) / n()^0.5)

#check mean +/- SEM 
alldata.mean %>%
  ggplot(aes(x = time, y = meanDelF)) + #plot a graph
  geom_line(aes(colour = genotype), alpha = 1) + #plot the mean
  geom_ribbon(aes(ymin = sem.low, ymax = sem.high, fill = genotype), alpha = 0.2) + #plot the standard error
  coord_cartesian(ylim = c(-.3, .30)) + #y axis limits
  geom_segment(aes(x = 30, xend = 60, y = 1.5, yend = 1.5), colour = "black", size = .5) + #indicates when the odor was on
  coord_cartesian(xlim = c(0, 90)) + #plot the x axis limits
  theme_classic() 

