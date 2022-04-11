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

### arguments
sourcepath = Sys.getenv("SOURCEPATH")
prefixdelim = sapply(strsplit(Sys.getenv("PREFIXDELIM"),":")[[1]],strsplit,",")
imethod = Sys.getenv("IMETHOD")
savepath = Sys.getenv("SAVEPATH")
outputmaster = Sys.getenv("OUTPUTMASTER")
CellChatDB = get(Sys.getenv("CELLCHATDB"))
pairedpriority = Sys.getenv("PAIREDPRIORITY")
commonsuff = Sys.getenv("COMMONSUFF")

dir.create(savepath,recursive=TRUE)

### external code
source(Sys.getenv("MULTISETLIB"))

run_signaling_multiset(
    savepath_master=savepath,
    outputmaster=outputmaster,
    sourcepath=sourcepath,
    CellChatDB=CellChatDB,
    prefixes=prefixdelim,
    paired_priority=pairedpriority,
    common=commonsuff)