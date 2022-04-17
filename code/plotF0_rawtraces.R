library(tidyverse)
theme_set(theme_classic())
#choose file in folder to analyze
setwd(dirname(file.choose()))

#######WT hex and hex-IAA fig 2A
WT_hex_2A <- vroom::vroom(file.choose())
WT_hexIAA_2A <- vroom::vroom(file.choose())

WT_hex_2A_F0 <- WT_hex_2A %>%
  group_by(animal, animal_num) %>%
  filter(time > 0, time < 5.25) %>%
  summarise(F0 = mean(MeanGCaMP), cue = "hex")

WT_hexIAA_2A_F0 <- WT_hexIAA_2A %>%
  group_by(animal, animal_num) %>%
  filter(time > 0, time < 5.25) %>%
  summarise(F0 = mean(MeanGCaMP), cue = "hex-IAA")

#### join both  for F0 plot####
(p2A <- rbind(WT_hex_2A_F0, WT_hexIAA_2A_F0) %>%
  ggplot(aes(x = cue, y = F0)) +
  ggbeeswarm::geom_quasirandom(aes(color = cue), 
                               width= 0.1,
                               size = 0.85) +
  scale_color_manual(values = c("grey", "lightblue")) +
  coord_cartesian(ylim = c(0,1500))) +
  stat_summary(aes(color = cue),
               geom = "errorbar", 
               fun.data = "mean_se", 
               width = 0.125) +
  stat_summary(aes(color = cue),
               geom = "crossbar", 
               fun = mean,
               width = .3,
               size = .25)

ggsave("./2A_F0plot.pdf")


rbind(WT_hex_2A_F0, WT_hexIAA_2A_F0) %>%
  write_csv(., file = "./2A_WT_F0.csv")

rbind(WT_hex_2A, WT_hexIAA_2A) %>%
  ggplot(aes(x = time, y = MeanGCaMP)) +
  geom_line(aes(group = animal, color = cue))

#######

###### make a function to make this more concise ####
Extract_F0s <- function(data,
                        cue,
                        genotype,
                        ...){
  data <- vroom::vroom(file.choose())
  data %>%
    group_by(animal, animal_num) %>%
    filter(time > 0, time < 5.25) %>%
    summarise(F0 = mean(MeanGCaMP), 
              cue = cue,
              genotype = genotype) %>%
    distinct()
}

##### for WT, odr-3, odr-1
odr3_hex <- Extract_F0s()
odr3_hexIAA <- Extract_F0s()
WT_comp_hex <- Extract_F0s()


### plot raw traces: ####
#WT from fig 2A:
# hex = 8/10, 8/31, hex-iaa = 8/11, 8/28
rbind(vroom::vroom(file.choose()), #hex
                  vroom::vroom(file.choose())) %>%
  ggplot(aes(x = time, y = MeanGCaMP)) +
  geom_line(aes(color = cue, group = animal), alpha = 0.75) +
  scale_color_manual(values = c("black", "steelblue1"))

#WT corressponding to odr-1, odr-3
# hex = 8/16, 8/28, 8/31, hex-iaa = 7/26, 7/29
WTmerged <- rbind(vroom::vroom(file.choose()), #hex
      vroom::vroom(file.choose()), #hex
      vroom::vroom(file.choose())) # hex- IAA
WTmerged %>%
  mutate(cue = case_when(
    cue == "hex" ~ "HEX",
    TRUE ~ cue)) %>%
  ggplot(aes(x = time, y = MeanGCaMP)) +
  geom_line(aes(color = cue, group = animal), alpha = 0.75) +
  scale_color_manual(values = c("black", "steelblue1"))

#odr-3
rbind(vroom::vroom(file.choose()), vroom::vroom(file.choose())) %>%
  ggplot(aes(x = time, y = MeanGCaMP)) +
  geom_line(aes(color = cue, group = animal), alpha = 0.75) +
  scale_color_manual(values = c("seagreen", "palegreen"))

#odr-1
rbind(vroom::vroom(file.choose()), vroom::vroom(file.choose())) %>%
  ggplot(aes(x = time, y = MeanGCaMP)) +
  geom_line(aes(color = cue, group = animal), alpha = 0.75) +
  scale_color_manual(values = c("mediumorchid4", "mediumorchid1"))

