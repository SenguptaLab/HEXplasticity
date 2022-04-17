
awc_hex<-plotGCaMP_multi(cue = HEX, 
                         genotype = WT,
                         neuron = AWC, 
                         center_on_pulse = 'ON',
                         show.plots = FALSE,
                         backsub = TRUE,
                         heatmap_limits = c(-0.5,0,0.5))

awc_hex$data 
awc_hex$plot

awc_hex_iaa<-plotGCaMP_multi(cue = HEX-IAA, 
                                genotype = WT,
                                neuron = AWC, 
                                center_on_pulse = 'ON',
                                show.plots = FALSE,
                                backsub = TRUE,
                                heatmap_limits = c(-0.5,0,0.5))
awc_hex_iaa$data 
awc_hex_iaa$plot 

alldata <- rbind(awc_hex$data, awc_hex_iaa$data)
#for two genotypes: alldata <- rbind(wt$data, HEX$data)

unnested <- unnest(alldata)
#unnest data to analyze for mean & SEM

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
  #geom_line(linetype = "dashed")+
  geom_ribbon(aes(ymin = sem.low, ymax = sem.high, fill = genotype), alpha = 0.2) + #plot the standard error
  coord_cartesian(ylim = c(-.3, 2)) + #y axis limits
  geom_segment(aes(x = 30, xend = 60, y = 2, yend = 2), colour = "black", size = .5) + #indicates when the odor was on
  coord_cartesian(xlim = c(0, 90)) + #plot the x axis limits
  theme_classic() 
