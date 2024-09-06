#!/usr/bin/Rscript

library(LEA)
library(lfmm)
library(vcfR)
library(qvalue)
setwd("./")
getwd()

.libPaths("~/R/library/")

##IMPUTATION FOR LFMM#####
ped2lfmm("./myplink.ped")
ped2geno("./myplink.ped")
snmf_result <- snmf("./myplink.geno", K= 1:4, project = "continue", repetitions = 5,
                    CPU = 2, alpha = 10,entropy = T )
b <- cross.entropy(snmf_result, K=2)
min_run <- which.min(b)
impute(snmf_result, "./myplink.lfmm", method = "mode", K= 2, run = min_run)

print("Everything is OK")
