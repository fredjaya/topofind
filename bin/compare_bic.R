#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
step_0 <- args[1]
bic_0 <- args[2]
step_1 <- args[3]
bic_1 <- args[4]

if (bic_0 < bic_1) {
       cat(step_0)
} else if (bic_1 < bic_0) {
	cat(step_1)
} else {
	cat("uh oh")
}
