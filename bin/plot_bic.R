library(dplyr)
library(ggplot2)

options(digits = 10)

t <-
  read.table("~/Dropbox/treemix_rc/04_testing/python/concat/bic.txt", fill = T, sep = "\t") %>%
  rename(Run = V1, BIC = V2) %>%
  mutate(nTrees = gsub("_.*", "", Run)) %>%
  mutate(Run = gsub(".*_mast_", "", Run))

t %>%
  ggplot(aes(x = nTrees, y=BIC)) +
  geom_jitter(width = 0.25, size = 1.5) +
  theme_light()
