#get F0 for each dataset - vroom is basically fast read_csv
#the group_by makes sure I operate on each animal separately
#the filter is to get the first 5s
#the summarize is to get the mean of the raw F values (MeanGCaMP)
#'distinct' just gets rid of the repeated lines from summarize output

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
#can call as view (Extract_F0s()) to see the data as a table