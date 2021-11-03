#' plotGCaMP_multi
#'
#' Function is a wrapper for exp.fit.all.log.lin which outputs corrected GCaMP signals as well as plots showing
#' original and corrected signals. Also plots average trace for GCaMP signal. Inputs are matfiles with a "signal" field,
#' or .csv files from imageJ analysis.
#' The script searches recursively for matfiles (*.mat) or csv files (\*.csv) from a file (can be any placeholder file). Need to exclude
#' matfiles that are not GCaMP files either using FileFilter, or by putting only relevant files in the folder. Default
#' startpulse assumes 400ms delay for camera recording.
#' Requires max_delta helper function
#' time etc...
#' @param FileFilter string to search/subset filenames
#' @param matlab whether to use matfiles or data from ImageJ quantification. Defaults to FALSE (matfiles)
#' @param genotype label the genotype for these data
#' @param cue label the stimulus cue.
#' @param food label to food cue
#' @param startPulse begin of stimulus
#' @param endPulse endPulse time of stimulus
#' @param center_on_pulse optional parameter to center delF values by the mean of the stimulus duration
#' 'OFF' = Bring values to mean delF of 2nd half of pulse duration, order by OFF responses
#' 'ON' = Bring values to mean delF of 2nd half pre-pulse duration, order by ON responses
#' @param show.plots render plots for baseline correction - defaults to TRUE
#' @param use.Fmax normalize all amplitudes to within 0-1
#' @param neuron neuron being analyzed
#' @param linear optional argument piped into exp.fit.all.log.lin include a linear term in the fit?
#' @param heatmap_limits optional 3-value vector defining the color scale and y axis limits, ie c(-1,0,2)
#' @param folderpath user supplied path for recursive file search. If missing, it will prompt for a fiel selection.
#' @param backsub do you want to subtract background ie subtract ROI of background from the mean cell ROI? Defaults to TRUE
#' @param exp.fit do you want to use raw deltaF/F values - ie exponentially corrected? Defaults to TRUE, which means it will be corrected
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @importFrom magrittr "%$%"
#' @export
#' @examples data <- plotGCaMP_multi(N2, genotype = N2, cue = octanol)
#'
plotGCaMP_multi <- function(FileFilter,
                            matlab = FALSE,
                            genotype = genotype,
                            cue = cue,
                            food = OP50,
                            startPulse = 29.5,
                            endPulse = 59.5,
                            center_on_pulse = "none",
                            show.plots = TRUE,
                            use.Fmax = FALSE,
                            neuron = GCAMP,
                            linear = FALSE,
                            nls = TRUE,
                            backsub = FALSE,
                            exp.fit = FALSE,
                            heatmap_limits = "auto",
                            folderPath,
                            ...) {
  
  message("this function is for analyzing calcium imaging data, for details,
          go to https://github.com/mikeod38/MF.matR")
  library(tidyverse)
  library(magrittr)
  library(patchwork)
  library(scales)
  FileFilter <- quo_name(enquo(FileFilter)) # make Filter usable inside other functions
  genotype <- quo_name(enquo(genotype))
  cue <- quo_name(enquo(cue))
  food <- quo_name(enquo(food))
  neuron <- quo_name(enquo(neuron))
  center_on_pulse <- quo_name(enquo(center_on_pulse))
  #heatmap_limits <- enquo(heatmap_limits)
  
  if (genotype == "genotype")
  {stop("error: You left out the 'genotype' argument, was this intentional?")}
  
  if (cue == "cue")
  {stop("error: You left out the 'cue' argument, was this intentional?")}
  
  
  if (!center_on_pulse %in% c("none", "ON", "OFF")) {
    stop("error: 'center_on_pulse' must be either 'none', 'ON'or 'OFF'")}
  
  
  if(missing(folderPath)) {
    folderPath <- dirname(file.choose())
  }
  message(paste("using files at or below the folder:", basename(folderPath)))
  
  #### import, format and correct for photobleaching ####
  if(matlab == TRUE) {
    files <- list.files(file.path(folderPath), pattern = "*.mat", recursive = TRUE)
    files <- files[stringr::str_detect(files, pattern = paste0(FileFilter))]
    filenames <- files
    files <- file.path(folderPath, files)
    #df <- data.frame(x = 1, genotype = genotype, cue = cue)
    
    data <- purrr::map(files, ~ exp.fit.all.log.lin(filename = .,
                                                    skip.time = 10,
                                                    show.plots = show.plots,
                                                    nls = nls,
                                                    startPulse = startPulse,
                                                    endPulse = ))
  } else { #for matlab = FALSE
    neuronfiles <- list.files(file.path(folderPath), pattern = "*neuron_results.csv", recursive = TRUE)
    neuronfiles <- neuronfiles[stringr::str_detect(neuronfiles, pattern = paste0(FileFilter))]
    neuronfilenames <- neuronfiles
    
    neuronfiles <- file.path(folderPath, neuronfiles)
    
    message("merging imageJ files")
    
    purrr::map(
      neuronfiles,
      ~ merge_FIJI_data(neuronfile = ., show.plots = show.plots, backsub = backsub))
    
    
    files <- list.files(file.path(folderPath), pattern = "*ImageJ_data.csv", recursive = TRUE)
    files <- files[stringr::str_detect(files, pattern = paste0(FileFilter))]
    filenames <- files
    files <- file.path(folderPath, files)
    message("correcting photobleaching")
    data <- purrr::map(files, ~ exp.fit.all.log.lin(
      filename = .,
      skip.time = 10,
      show.plots = show.plots,
      matlab = FALSE,
      linear = linear,
      nls = nls
    ))
    
  }
  
  
  data %<>% tibble(data = .,
                   animal = filenames,
                   animal_num = factor(seq(from = 1, to = length(filenames))),
                   genotype = genotype,
                   cue = cue,
                   food = food,
                   neuron = neuron)
  
  if (data %>%
      unnest(cols = c(data)) %>%
      group_by(animal) %>%
      tally() %$%
      unique(n) %>%
      length() != 1) stop("error: One or more of your video files are of different length, check these files")
  data %>%
    unnest(cols = c(data)) %>%
    group_by(animal) %>%
    tally() %>%
    filter(n < max(n)) %>%
    select(animal)
  
  
  
  #### recenter mean values ####
  if(center_on_pulse == "OFF") {
    
    means <- data %>% unnest(cols = c(data)) %>%
      group_by(animal) %>%
      filter(time > (startPulse+endPulse)/2 & time < endPulse) %>%
      summarise(mean_pulse_delF = mean(delF))
    
    data <- left_join(data, means, by = "animal")
    
    data %<>% unnest(cols = c(data)) %>%
      group_by(animal, animal_num) %>%
      mutate(delF = delF - mean_pulse_delF) %>%
      nest()
  }
  
  if(center_on_pulse == "ON") {
    
    means <- data %>% unnest(cols = c(data)) %>%
      group_by(animal) %>%
      filter(time > startPulse/2 & time < startPulse) %>%
      summarise(mean_pulse_delF = mean(delF))
    
    data <- left_join(data, means, by = "animal")
    
    data %<>% unnest(cols = c(data)) %>%
      group_by(animal,animal_num) %>%
      mutate(delF = delF - mean_pulse_delF) %>%
      nest()
  }
  
  #### arrange heat map settings ####
  
  if(center_on_pulse == "OFF") {
    message("centering on OFF response")
    plot_order <- data %>%
      unnest(cols = c(data)) %>%
      group_by(animal, animal_num, .drop = TRUE) %>%
      summarise(maxD = max_delta(delF, end = endPulse)) %>%
      arrange(maxD)
    
    # breaks = c(-.5,0,1.5)
    # labels = c("-.5", "0", "1.5")
    # limits = c(-0.5,1.5)
  }
  
  if(center_on_pulse %in% c("ON", "none")) {
    message("centering on ON response")
    plot_order <- data %>% unnest(cols = c(data)) %>%
      group_by(animal, animal_num, .drop = TRUE) %>%
      summarise(maxD = max_delta(delF, end = startPulse)) %>%
      arrange(maxD)
    
    # breaks = c(-.5,0,1.5)
    # labels = c("-.5", "0", "1.5")
    # limits = c(-0.5,1.5)
  }
  
  if(use.Fmax == TRUE) {
    range <- data %>% unnest(cols = c(data)) %>%
      group_by(animal) %>%
      summarise(max_delF = max(delF), Fo = quantile(delF, 0.05))
    
    data <- full_join(data, range)
    
    data %<>% unnest(cols = c(data)) %>%
      group_by(animal,animal_num) %>%
      mutate(delF = (delF - Fo) / (max_delF - Fo)) %>%
      nest()
    
    # breaks = c(0,0.5,1)
    # labels = c("0", "0.5", "1")
    # limits = c(0,1)
    
  }
  
  if(!is.numeric(heatmap_limits)) { # using auto calc unless a numeric vector input
    breaks <- round(
      data %>% unnest(cols = c(data)) %$% quantile(delF, c(0.05, 0.5, 0.99)),
      2
    )
    labels <- as.character(breaks)
    limits <- breaks[c(1,3)]
  } else {
    breaks <- heatmap_limits
    labels <- as.character(breaks)
    limits <- breaks[c(1,3)]
  }
  
  if(exp.fit == FALSE) {
    plot1 <- data %>% unnest(cols = c(data)) %>%
      ggplot(aes(x = time, y = delF)) +
      geom_line(aes(group = animal), alpha = 0.1) +
      geom_smooth(method = "loess", span = 0.05) +
      theme_classic() +
      geom_segment(aes(x = !!startPulse,
                       y = limits[2],
                       xend = !!endPulse,
                       yend = limits[2])) +
      annotate(geom = "text",
               label = cue,
               y = 1.2*limits[2],
               x = (startPulse + endPulse)/2,
               size = 6) +
      annotate(geom = "text",
               label = genotype,
               y = 1.2*limits[2],
               x = 10,
               fontface = "italic",
               size = 6) +
      coord_cartesian(ylim = c(limits[1], 1.5*limits[2])) +
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.title.x = element_blank())
    
    plot2 <-  full_join(data, plot_order) %>%
      unnest(cols = c(data)) %>%
      group_by(animal_num) %>%
      ggplot(aes(x = time, y = fct_reorder(animal_num, maxD))) +
      geom_tile(aes(fill = delF)) +
      scale_fill_viridis_c(option = "magma",
                           breaks = breaks,
                           labels = labels,
                           limits = limits,
                           oob =squish) +
      theme_classic() +
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.text.y = element_blank()) +
      labs(y = "Animal number")
  } else {
    plot1 <- data %>% unnest(cols = c(data)) %>%
      ggplot(aes(x = time, y = signal)) +
      geom_line(aes(group = animal), alpha = 0.1) +
      geom_smooth(method = "loess", span = 0.05) +
      theme_classic() +
      geom_segment(aes(x = !!startPulse,
                       y = limits[2],
                       xend = !!endPulse,
                       yend = limits[2])) +
      annotate(geom = "text",
               label = cue,
               y = 1.2*limits[2],
               x = (startPulse + endPulse)/2,
               size = 6) +
      annotate(geom = "text",
               label = genotype,
               y = 1.2*limits[2],
               x = 10,
               fontface = "italic",
               size = 6) +
      coord_cartesian(ylim = c(limits[1], 1.5*limits[2])) +
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.title.x = element_blank())
    
    plot2 <-  full_join(data, plot_order) %>%
      unnest(cols = c(data)) %>%
      group_by(animal_num) %>%
      ggplot(aes(x = time, y = fct_reorder(animal_num, maxD))) +
      geom_tile(aes(fill = signal)) +
      scale_fill_viridis_c(option = "magma",
                           breaks = breaks,
                           labels = labels,
                           limits = limits,
                           oob =squish) +
      theme_classic() +
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.text.y = element_blank()) +
      labs(y = "Animal number")
  }
  
  
  plots <- plot1 + plot2 + plot_layout(ncol = 1, heights = c(2,1))
  
  if(use.Fmax == TRUE) {
    ggsave(plots, filename = file.path(folderPath,paste0(genotype,"_",cue,neuron,"_delFmaxplots.png")),
           width = 11, height = 8.5, units = "in")
  } else {
    ggsave(plots, filename = file.path(folderPath,paste0(genotype,"_",cue,"_",neuron,"_",food,"_plots.png")),
           width = 11, height = 8.5, units = "in")
  }
  write_csv(unnest(data, cols = c(data)), path = file.path(folderPath,paste0(genotype,"_",cue,"_",neuron,"_",food,"_rawdata.csv")))
  
  
  return(list(data = dplyr::full_join(data, plot_order), plot = plots))
}