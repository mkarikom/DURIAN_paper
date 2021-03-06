# Code adapted from https://github.com/favilaco/deconv_benchmark:
# 
# Copyright 2019 Francisco Avila Cobos, 2021 Matt Karikomi
library(doParallel)
library(DURIAN,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")

### arguments
datapath = Sys.getenv("DATAPATH")
etclib = Sys.getenv("ETCLIB")

nsimulations = as.integer(Sys.getenv("NSIMINTERNAL")) # the number of times to regenerate a set of npseudobulk pseudobulk samples from the given sc source data

### external code
source(etclib)

analyze_ursm_tmp_files(pbdir=datapath)
