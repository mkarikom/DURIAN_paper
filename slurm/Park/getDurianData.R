# load the libraries

library(dplyr)
library(Seurat)
library(scater)
library(edgeR)
library(gridExtra)
library(cowplot)
library(data.table)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(RColorBrewer)
library(GEOquery)
library(biomaRt)
library(org.Hs.eg.db)
library(ArrayExpress)
library(readr)
library(AnnotationHub)
library(tximport)
library(DURIAN)

project_dir = "/share/crsp/lab/cellfate/mkarikom/DURIAN_paper_clean"

data_dir = file.path(project_dir,"slurm/Park/durian_data")

#############################################################################
# Beckerman 2017 (https://doi.org/10.1038/nm.4287)
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81492)
#############################################################################
print("acquiring bulk data")

bulkdata <- readRDS(url("https://xuranw.github.io/MuSiC/data/Mousebulkeset.rds","rb"))

Beckerman_T = exprs(bulkdata)[,1:3]
write.csv(t,file = file.path(data_dir,"Beckerman_T.csv"))

#############################################################################
# Park 2018 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107585) 
# single-cell (https://doi.org/10.1126/science.aar2131)
#############################################################################
# this data was downloaded separately

gse=getGEO(filename="/share/crsp/lab/cellfate/mkarikom/DURIAN_paper_clean/slurm/Park/raw_data/GSE107585_series_matrix.txt")

scdata <- readRDS(url("https://xuranw.github.io/MuSiC/data/Mousesubeset.rds","rb"))
pdata = pData(scdata@phenoData)
pdata = data.frame(cellID=row.names(pdata),cellType=pdata$cellType,sampleID=pdata$sampleID)
pdata = pdata[!(pdata$cellType %in% c("Novel1","Novel2","Novel3")),]

rownames(pdata) = pdata$cellID
c = exprs(scdata)
colnames(c) = stringr::str_replace(colnames(c), "\\.", "-")

ParkWT_pDataC = pdata[pdata$sampleID %in% which(gse@phenoData@data[,"genotype:ch1"]=="WT"),]
ParkWT_C = data.frame(c[,colnames(c) %in% rownames(ParkWT_pDataC)])
colnames(ParkWT_C) = stringr::str_replace(colnames(ParkWT_C), "\\.", "-")

rownames(ParkWT_pDataC) = stringr::str_replace(rownames(ParkWT_pDataC), "-", "\\.")
ParkWT_pDataC$cellID = rownames(ParkWT_pDataC)

write.csv(ParkWT_pDataC,file = file.path(data_dir,"ParkWT_pDataC.csv"))
write.csv(ParkWT_C,file = file.path(data_dir,"ParkWT_C.csv"))

# # use seurat vst thresholding for genes
seur = CreateSeuratObject(counts=ParkWT_C)
seur = NormalizeData(seur)
genes = VariableFeatures(FindVariableFeatures(seur, selection.method = "vst", nfeatures = 5000))

commgenes = intersect(genes,rownames(Beckerman_T))
ParkWT.VST_C = subsetsc(ParkWT_C,geneids=commgenes,return_obj=TRUE,nsd=3)
ParkWT.VST_pDataC = ParkWT_pDataC[ParkWT_pDataC$cellID %in% colnames(ParkWT.VST_C),]

write.csv(ParkWT.VST_pDataC,file = file.path(data_dir,"ParkWT.VST_pDataC.csv"))
write.csv(ParkWT.VST_C,file = file.path(data_dir,"ParkWT.VST_C.csv"))
write.csv(Beckerman_T[commgenes,],file = file.path(data_dir,"Beckerman.VST_T.csv"))

set.seed(1)
target = 1000
sampn = min(nrow(ParkWT_C),target)

ParkWT.VST1K05_C = as.data.frame(t(as.data.frame(t(ParkWT_C)) %>% sample_n(sampn)))
ParkWT.VST1K05_C = subsetsc(scremoutlier(ParkWT.VST1K05_C),generate=0.05,return_obj=TRUE,nsd=3)
ParkWT.VST1K05_pDataC = ParkWT.VST_pDataC[colnames(ParkWT.VST1K05_C),]

write.csv(ParkWT.VST1K05_C,file=file.path(data_dir,"ParkWT.VST1K05_C.csv"))
write.csv(ParkWT.VST1K05_pDataC,file=file.path(data_dir,"ParkWT.VST1K05_pDataC.csv"))
write.csv(Beckerman_T[rownames(Beckerman_T) %in% rownames(ParkWT.VST1K05_C),],file=file.path(data_dir,"Beckerman.VST1K05_T.csv"))