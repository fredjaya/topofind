#!/usr/bin/Rscript

args <-  commandArgs(trailingOnly=TRUE)
path <- args[1]

library(dplyr)

models <- 
  read.table(path, header = T, check.names = F) %>%
  select(-`No.`) %>% 
  mutate(SubModel = gsub("\\+.+$", "", Model)) %>%
  mutate(Freq = if_else(grepl("\\+F", Model), T, F)) %>%
  mutate(Invariant = if_else(grepl("\\+I", Model), T, F)) %>%
  mutate(FreeRates = if_else(grepl("\\+R", Model), gsub("^.+\\+R", "", Model), "1")) %>%
  mutate(FreeRates = as.numeric(FreeRates))

write.table(models, "models_parsed.tsv", sep = "\t", quote = F)

bic_rk <- 
  models %>%
  group_by(FreeRates) %>%
  summarise(Model = Model[which.min(BIC)], BIC = min(BIC))

write.table(bic_rk, "best_bic_per_rk.tsv", sep = "\t", quote = F)

# Is the BIC of the best +R2 model better than R1?
bic_rk_only <- bic_rk %>% pull(BIC)
r2_lt_r1 <- bic_rk_only[1] > bic_rk_only[2]

write.table(r2_lt_r1, "r2_lt_r1", quote = F, row.names = F, col.names = F)
