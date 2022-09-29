#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=TRUE)
sitelh <- args[1]
alninfo <- args[2]
prefix <- args[3]

library(ggplot2)
library(MixtureModelHMM)

# Run HMM
hmm_sitelh <- run_HMM(site_info = sitelh, aln_info = alninfo)

# Plot site assignment
ggsave(filename = paste(prefix, "_site_assignment.png", sep = ""), device = "png",
       plot = hmm_sitelh$alignment_plot)

# Output partition file
save_partitioning_scheme(hmm_result = hmm_sitelh, output_filename = paste(prefix, ".partition", sep = ""))
