library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(MuSiC)
library(Biobase)
library(xbioc)
library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(reshape2, include.only = c("melt"))
library(gridExtra)

savepath = Sys.getenv("SAVEPATH")
outputmaster = Sys.getenv("RUNMASTER")
sparsityparam = Sys.getenv("SPARSITY_PARAM")

dir.create(savepath,recursive=TRUE)

### external code
metricslib="slurm/scrabble_helper_functions/library_cluster_metrics.R"
source(metricslib)

run_clusterMetrics_final_nested_sim(
    savepath_master=savepath,
    outputmaster=outputmaster,
    subdir="scrabble_fam",
    excludemodels=c("DrImpute","mtSCRABBLE"),
    maxwallclock=20,
    sparsityparam=sparsityparam)