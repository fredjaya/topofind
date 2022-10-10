#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=TRUE)
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

write.table(models, "modelfinder_parsed.tsv", sep = "\t", quote = F, row.names = F)

bic_rk <- 
  models %>%
  group_by(FreeRates) %>%
  summarise(Model = Model[which.min(BIC)], BIC = min(BIC)) %>%
  select(Model, BIC)

write.table(bic_rk, "models_summary.tsv", sep = "\t", quote = F, row.names = F)
