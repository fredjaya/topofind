#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
sitelh <- args[1]
alninfo <- args[2]

library(ggplot2)
library(MixtureModelHMM)

# Run HMM
hmm_sitelh <- run_HMM(site_info = sitelh, aln_info = alninfo)

# Plot site assignment
ggsave(filename = "site_assignment.png", device = "png",
       plot = hmm_sitelh$alignment_plot)

# Output partition file
save_partitioning_scheme(hmm_result = hmm_sitelh, output_filename = "r2.partition")
