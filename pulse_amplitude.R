#' pulse_amplitudes
#'
#' Using output from plotGCaMP_multi() raw data - calculate the amplitude of stimulus addition and
#' removal calcium responses
#' @param FileFilter string to search/subset filenames
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @importFrom magrittr "%$%"
#'
#' @export
#' @examples data <- pulse_amplitudes()
#'

pulse_amplitude <- function(stim_on = 30,
                            stim_off = 60,
                            method = "max") {
  
  message("Select a file in the folder you want to analyze")
  filename <- file.choose()
  folderPath <- dirname(filename)
  message(paste("opening raw data in folder", folderPath))
  
  data <- readr::read_csv(filename)
  
  if(method == "average") {
    
    (prePulse <- data %>%
       filter(time > 27 & time < 29) %>%
       group_by(animal, animal_num) %>%
       summarize(prePulse = mean(delF)))
    
    (stim_ON <- data %>%
        filter(time > 30.5 & time < 35) %>%
        group_by(animal, animal_num) %>%
        summarize(stim_ON = mean(delF)))
    
    (preOFF <- data %>%
        filter(time > 57 & time < 59) %>%
        group_by(animal, animal_num) %>%
        summarize(preOFF = mean(delF)))
    
    (stim_OFF <- data %>%
        filter(time > 60.5 & time < 65) %>%
        group_by(animal, animal_num) %>%
        summarize(stim_OFF = mean(delF)))
    
    data <- cbind(prePulse,
                  stim_ON = stim_ON$stim_ON,
                  preOFF = preOFF$preOFF,
                  stim_OFF = stim_OFF$stim_OFF) %>%
      mutate(stim_ON_delF = stim_ON - prePulse,
             stim_OFF_delF = stim_OFF - preOFF)
  } else {
    
    (prePulse <- data %>%
       filter(time > (stim_on - 3) & time < (stim_on - 1)) %>%
       group_by(animal, animal_num) %>%
       summarize(prePulse = mean(delF)))
    
    (preOFF <- data %>%
        filter(time > (stim_off - 3) & time < (stim_off - 1)) %>%
        group_by(animal, animal_num) %>%
        summarize(preOFF = mean(delF)))
    
    stim_ON_data <- data %>%
      filter(time > (stim_on + 0.5) & time < (stim_on + 10.5)) %>%
      # add back a column of pre-pulse mean values
      full_join(prePulse) %>%
      mutate(difference = delF - prePulse) %>%
      group_by(animal, animal_num) %>%
      summarize(time = time, difference = difference, maxAbs = max(abs(difference))) %>%
      filter(abs(difference) == maxAbs) %>%
      mutate(stim_ON_delF = difference) %>%
      select(animal, animal_num, time_at_maxON = time, stim_ON_delF)
    
    stim_OFF_data <- data %>%
      filter(time > (stim_off + 0.5) & time < (stim_off + 10.5)) %>%
      # add back a column of pre-pulse mean values
      full_join(prePulse) %>%
      mutate(difference = delF - prePulse) %>%
      group_by(animal, animal_num) %>%
      summarize(time = time, difference = difference, maxAbs = max(abs(difference))) %>%
      filter(abs(difference) == maxAbs) %>%
      mutate(stim_OFF_delF = difference) %>%
      select(animal, animal_num, time_at_maxOFF = time, stim_OFF_delF)
    
    data <- full_join(stim_ON_data, stim_OFF_data)
  }
  
  write_csv(data, file.path(folderPath, paste0(basename(filename),"_",method,"_amplitudes.csv")))
  
}