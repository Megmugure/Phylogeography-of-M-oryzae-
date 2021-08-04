#! /usr/bin/env Rscript

setwd("/home/mwanjiku/R_analysis/R_output")
library(adegenet)
library(parallel)


dat2 <- fasta2genlight("/home/mwanjiku/merged_vcf/fasta_file3.fasta", chunkSize = 20, parallel = TRUE, n.cores = 20)
save(dat2,file = "output2")


