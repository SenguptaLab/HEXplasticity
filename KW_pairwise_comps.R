#KW then posthoc pairwise Wilcoxon test followed by manual adjustment for p value by Benjamini & Hochberg method
# KW test (nonparametric equivalent for ANOVA), pairwise wilcox test (nonparametric equivalent of t-test)
#Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends 
#the two-samples Wilcoxon test in the situation where there are more than two groups
#posthoc tests for pairwise comarisons
#he "BH" (aka "fdr") and "BY" method of Benjamini, Hochberg, and Yekutieli control the false discovery rate, 
#the expected proportion of false discoveries amongst the rejected hypotheses. The false discovery rate is a less 
#stringent condition than the family-wise error rate, so these methods are more powerful than the others (e.g. Bonferroni)

KW_pairwise_comps <- function(dataset) {
  library(tidyverse)
  data <- read_csv(file.choose())
  data <- data %>%
    #janitor::clean_names() %>%
    pivot_longer(names_to = "group", values_to = "amplitude", cols = everything())
  KWout <- kruskal.test(data$amplitude ~ data$group)
  comparisons <- pairwise.wilcox.test(data$amplitude, data$group,
                                      p.adjust.method = "BH", pool.sd = FALSE)
  return(list(KWout, comparisons))
}
KW_pairwise_comps()
