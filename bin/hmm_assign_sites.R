#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=TRUE)
sitelh <- args[1]
alninfo <- args[2]

library(ggplot2)
library(MixtureModelHMM)

# Run HMM
hmm_sitelh <- run_HMM(
  site_info = "~/Dropbox/treemix_rc/04_testing/mast_sim/data1a/t1_r2.sitelh",
  aln_info = "~/Dropbox/treemix_rc/04_testing/mast_sim/data1a/t1_r2.alninfo"
  )

# Plot site assignment
hmm_sitelh$alignment_plot
ggsave(filename = "_site_assignment.png", device = "png")

# Output partition file
save_partitioning_scheme(hmm_result = hmm_sitelh, output_filename = ".partition")
