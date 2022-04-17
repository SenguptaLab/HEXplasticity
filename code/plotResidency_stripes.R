#' plotResidency_stripes
#'
#' Plots residency data by inside, outside stripe. Generates a heatmap by experiment
#' @param FileFilter string to search/subset filenames
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @importFrom magrittr "%$%"
#'
#' @export
#' @examples data <- plotResidency_stripes()
#'
plotResidency_stripes <- function(FileFilter,
                                  folderPath,
                                  arena_size = 16.1,
                                  frame_rate = 2,
                                  vid.length = 20,
                                  y_bins = 50,
                                  y_max = 5,
                                  heatmap_limits,
                                  heatmap_palette = "Greys",
                                  plot.direction = 1,
                                  ...) {
  
  message("select a file in the folder you want to analyze")
  library(tidyverse)
  library(scales)
  ggplot2::theme_set(theme_classic())
  #### making an interactive option for the base folder:
  if(missing(folderPath)) {
    folderPath <- dirname(file.choose())
  }
  message(paste("using files at or below the folder:", basename(folderPath)))
  plot.direction <- plot.direction
  heatmap_palette <- quo_name(enquo(heatmap_palette))
  breaks <- heatmap_limits
  labels <- as.character(breaks)
  limits <- breaks[c(1,3)]
  
  #### select csv files that match the pattern, with an optional FileFilter
  files <- list.files(file.path(folderPath), pattern = "*.csv", recursive = TRUE)
  files <- files[stringr::str_detect(files, pattern = paste0(FileFilter))]
  filenames <- files
  files <- file.path(folderPath, files)
  files
  
  #get pixelsize in pixels / mm
  pixelsize <- files %>%
    stringr::str_subset(., pattern = "preprocess.csv", negate = FALSE) %>%
    read_csv(.) %>% select(pixelSize) %>% as.numeric()
  
  #### correct and determine cue boundaries by luminance:
  luminance <- files %>%
    stringr::str_subset(., pattern = "luminance.csv", negate = FALSE) %>%
    data.table::fread(., select = c(1,frame_rate*vid.length*60)) %>%
    tibble()
  
  luminance %<>%
    mutate(ypos = row_number(),
           y_mm = ypos / pixelsize) %>%
    rename(frame1 = V1, frameLast = V2400) %>%
    pivot_longer(cols = starts_with("frame"),
                 names_to = "frame",
                 values_to = "luminance")
  
  ybinwidth <- arena_size / max(luminance$ypos)
  
  # linear regression to compensate for luminance gradient across devices
  luminance <- luminance %>%
    group_by(frame) %>%
    nest() %>%
    # map a linreg to each frame values, then get just the slope and intercept values
    mutate(lmobj = map(data, function(.x) {
      lm(luminance ~ ypos, data = .x) %>% broom::tidy()
    })) %>%
    unnest(lmobj) %>%
    select(-(5:7)) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    rename(Intercept = `(Intercept)`, slope = ypos) %>%
    #now correct the slope and binarize:
    unnest(cols = c(data)) %>%
    mutate(norm_lum = luminance - (ypos*slope + Intercept)) %>%
    # the data are still wobbly, so I will use a Loess fit to smooth these
    # out and use these to binarize the luminance data
    nest() %>%
    mutate(smoothed = map(data, function(.x) {
      loess(norm_lum ~ ypos, data = .x, span = 0.1) %>% broom::augment()
    })) %>%
    unnest(smoothed) %>%
    mutate(lum_bin = case_when(
      .fitted >= 0 ~ "buffer",
      TRUE ~ "dye"
    ))
  
  #check if boundaries are greater than one worm length:
  first_bound <- luminance %>%
    filter(lum_bin == "dye") %>%
    slice(1) %$% (ypos[1] - ypos[2]) * ybinwidth > 0.5
  
  try(if(first_bound) stop("boundaries unstable"))
  
  second_bound <- luminance %>%
    filter(ypos > max(ypos)/2, lum_bin == "buffer") %>%
    slice(1) %$% (ypos[1] - ypos[2]) * ybinwidth > 0.5
  
  try(if(second_bound) stop("boundaries unstable"))
  #
  #
  luminance <- luminance %>%
    filter(frame == "frame1") %>%
    ungroup() %>%
    select(ypos, lum_bin)
  
  ypos_all <- luminance %>%
    mutate(y_mm = ypos / pixelsize) %>%
    select(y_mm, lum_bin)
  
  #### use y position data to generate relative residence
  ymat <- files %>%
    stringr::str_subset(., pattern = "ymat.csv", negate = FALSE)  %>%
    data.table::fread() %>%
    tibble() %>%
    mutate(worm = row_number()) %>%
    pivot_longer(cols = -worm, names_to = "time", values_to = "ypos") %>%
    separate("time", into = c("spacer", "time"), sep = "V") %>%
    mutate(time = as.numeric(time),
           ypos = round(ypos,0)) %>%
    select(worm, time, ypos)
  
  #merge luminance and yposition:
  raw_residence <- full_join(luminance, ymat) %>%
    filter(!is.na(worm), !is.na(ypos)) %>%
    #censor outermost 1mm due to edge effects
    filter(!ypos < pixelsize, !ypos > max(ypos)-pixelsize) %>%
    mutate(y_mm = ypos / pixelsize)
  
  # calculate relative che index:
  #
  left_bound <-
    ypos_all %>%
    # raw_residence %>%
    filter(lum_bin == "dye") %>% slice(1) %>%
    select(y_mm) %>% as.numeric()
  
  right_bound <-
    ypos_all %>%
    # raw_residence %>%
    filter(y_mm > 6, lum_bin == "buffer") %>% slice(1) %>%
    select(y_mm) %>% as.numeric()
  
  length_buffer <- (left_bound - 1) + (max(ypos_all$y_mm) - right_bound)
  length_dye <- right_bound - left_bound
  total_length <- max(ypos_all$y_mm) - 1
  
  track_counts <- raw_residence %>%
    count(y_mm) %>%
    full_join(ypos_all) %>%
    mutate(n = ifelse(is.na(n),0,n)) %>%
    group_by(lum_bin) %>%
    summarize(n = sum(n)) %>%
    pivot_wider(names_from = lum_bin, values_from = n) %>%
    mutate(n_tracks = buffer + dye,
           length_buffer = length_buffer,
           length_dye = length_dye,
           length_total = total_length,
           norm_tot_buf = (buffer / (length_buffer/total_length)),
           norm_tot_dye = dye / (length_dye/total_length),
           index = (norm_tot_dye - norm_tot_buf ) / (norm_tot_dye + norm_tot_buf))
  
  #convert missing y positions to 0
  rel_residence <- raw_residence %>%
    count(y_mm) %>%
    full_join(ypos_all) %>%
    mutate(n = ifelse(is.na(n),0,n)) %>%
    mutate(ybin = cut(y_mm, y_bins),
           ybin_numeric = as.numeric(as.factor(ybin))) %>%
    group_by(ybin,ybin_numeric,lum_bin) %>%
    summarize(count = sum(n),) %>%
    ungroup() %>%
    mutate(y_mm = seq(0,16.2, length.out = nrow(.)))
  
  mean_res <- mean(rel_residence$count)
  
  rel_residence <- rel_residence %>%
    mutate(relres = count / mean_res)
  #to plot histogram
  
  p.histogram <- raw_residence %>%
    ggplot(aes(y = y_mm)) +
    geom_histogram(aes(x = stat(count) / mean(count), fill = lum_bin), bins = y_bins) +
    labs(y = "position (mm)",
         x = "relative residence") +
    scale_fill_manual(values = c("grey", "lightblue")) +
    coord_cartesian(xlim=c(0,y_max), ylim = c(1,15), expand = FALSE) +
    theme(legend.position = "bottom")
  
  # to plot heatmap
  p.heatmap <- rel_residence %>%
    ggplot(aes(y = y_mm)) +
    geom_raster(aes(fill = relres, x = 1)) +
    #     fill = stat(count) / mean(count),
    #     color = stat(count) / mean(count)),
    # na.rm = FALSE,
    # geom = "tile",
    # position = "identity",
    # bins = y_bins) +
    coord_cartesian(xlim=c(0.5,1.5), ylim = c(0,16.2), expand = FALSE) +
    labs(fill = "relative residence",
         y = "position (mm)") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    scale_y_continuous(breaks = c(0,5,10,15)) +
    scale_fill_distiller(
      breaks = breaks,
      labels = labels,
      limits = limits,
      oob=squish,
      palette = heatmap_palette,
      direction = plot.direction)
  
  
  write_csv(track_counts, file = file.path(folderPath,paste0(basename(folderPath),"index.csv")))
  write_csv(raw_residence, file = file.path(folderPath,paste0(basename(folderPath),"raw_residence.csv")))
  write_csv(rel_residence, file = file.path(folderPath,paste0(basename(folderPath),"rel_residence.csv")))
  
  
  ggsave(plot = p.heatmap,
         filename = file.path(folderPath,paste0(basename(folderPath),"_heatmap.pdf")),
         width = 1,
         height = 4,
         units = "in")
  
  ggsave(plot = p.histogram,
         filename = file.path(folderPath,paste0(basename(folderPath),"_histogram.pdf")),
         width = 4,
         height = 4,
         units = "in")
  
  return(list(raw_residence, rel_residence, p.histogram, p.heatmap))
  
}