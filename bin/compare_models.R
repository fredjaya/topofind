library(dplyr)

# args1 <- params.out/aln.simpleName

models <- 
  read.table("~/Dropbox/treemix_rc/04_testing/mast_sim/data1a/models3", header = T, check.names = F) %>%
  select(-`No.`) %>% 
  mutate(SubModel = gsub("\\+.+$", "", Model)) %>%
  mutate(Freq = if_else(grepl("\\+F", Model), T, F)) %>%
  mutate(Invariant = if_else(grepl("\\+I", Model), T, F)) %>%
  mutate(FreeRates = if_else(grepl("\\+R", Model), gsub("^.+\\+R", "", Model), "1")) %>%
  mutate(FreeRates = as.numeric(FreeRates))

write.table(models, "models_parsed.tsv", sep = "\t")

bic_rk <- 
  models %>%
  group_by(FreeRates) %>%
  summarise(Model = Model[which.min(BIC)], BIC = min(BIC))

write.table(bic_rk, "best_bic_per_rk.tsv", sep = "\t")