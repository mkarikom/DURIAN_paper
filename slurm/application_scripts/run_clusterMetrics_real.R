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
library(ggpubr)
library(ggh4x)

outputmaster = Sys.getenv("OUTPUTMASTER")
savepath = Sys.getenv("SAVEPATH")

dir.create(savepath,recursive=TRUE)

### external code
metricslib="slurm/application_scripts/test/library_cluster_metrics.R"
source(metricslib)

metrics_df = run_clusterMetrics_final_nested_real(
    savepath_master=savepath,
    outputmaster=outputmaster,
    subdir="scrabble_fam",
    excludemodels=c(),
    maxwallclock=500,
    exclude.sa=c(1e-2))