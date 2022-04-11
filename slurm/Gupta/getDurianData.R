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
# library(DropletUtils)

# temp_data_dir = "/tmp/mkarikom/Gupta"
# dir.create(file.path(temp_data_dir,"durian_data"),recursive=TRUE)
# dir.create(file.path(temp_data_dir,"raw_data"),recursive=TRUE)
data_dir = "/dfs5/bio/mkarikom/temp/DURIAN/slurm/Gupta"
dir.create(file.path(data_dir,"durian_data"),recursive=TRUE)
dir.create(file.path(data_dir,"raw_data"),recursive=TRUE)
source("/dfs5/bio/mkarikom/temp/DURIAN/slurm/scrabble_helper_functions/library_scrabble.R")
celltypes = ("alpha-beta-gamma-delta-epsilon")

#############################################################################
# Biggs 2018 (https://doi.org/10.7554/eLife.36468.001) 
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110459) at e13.5 
#############################################################################
print("acquiring bulk data")

# destdir=file.path(temp_data_dir,"raw_data","GSE110459")
destdir=file.path(data_dir,"raw_data","GSE110459")
dir.create(destdir,recursive = TRUE)
getGEOSuppFiles(GEO="GSE110459", makeDirectory = FALSE, baseDir = destdir,fetch_files = TRUE, filter_regex = NULL)
df = read.csv(file.path(data_dir,"raw_data","GSE110459","GSE110459_raw_counts.txt.gz"), row.names=1, sep="")

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl",host="uswest.ensembl.org")
mgi_symbols = getBM(attributes=c('mgi_symbol','ensembl_gene_id'), 
      filters = 'ensembl_gene_id', 
      values = rownames(df), 
      mart = ensembl)

nondup = mgi_symbols[!duplicated(mgi_symbols[,1]),]
df_nondup = df[nondup$ensembl_gene_id,]
rownames(df_nondup) = nondup$mgi_symbol
print("writing bulk data")
# write.csv(df_nondup,file = gzfile(file.path(temp_data_dir,"durian_data","BiggsBulk_T.csv.gz")))
# write.csv(df_nondup,file = file.path(temp_data_dir,"durian_data","BiggsBulk_T.csv"))
write.csv(df_nondup,file = file.path(data_dir,"durian_data","BiggsBulk_T.csv"))

#############################################################################
# Gupta 2019 (https://doi.org/10.1016/j.devcel.2018.11.032) 
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122043) 
#############################################################################
print("acquiring sc data")
# load the labels from suoqin
e13_fn = "/dfs5/bio/mkarikom/code/DURIAN/slurm/Gupta/sc_data/GSE122043_Suoqin/data_embryonic_E13allCells_clustering_CellChat.rda"
e14_fn = "/dfs5/bio/mkarikom/code/DURIAN/slurm/Gupta/sc_data/GSE122043_Suoqin/data_embryonic_E14_subclusteringDC_BASAL_MELA_CellChat.rda"
load(e13_fn, e13 <- new.env())
load(e14_fn, e14 <- new.env())

# get pDataC info for suoqin data
e13_donors = sapply(strsplit(names(e13$data_embryonic$labels),"_"),function(x) x[1])
e14_donors = sapply(strsplit(names(e14$data_embryonic$labels),"_"),function(x) x[1])

e13_ids = names(e13$data_embryonic$labels)
e14_ids = names(e14$data_embryonic$labels)

e13_labels = unname(e13$data_embryonic$labels)
e14_labels = unname(e14$data_embryonic$labels)

e13_pDataC = data.frame(cellID=e13_ids,cellType=e13_labels,sampleID=e13_donors)
e14_pDataC = data.frame(cellID=e14_ids,cellType=e14_labels,sampleID=e14_donors)

rownames(e13_pDataC)=e13_pDataC$cellID
rownames(e14_pDataC)=e14_pDataC$cellID

print("writing sc metadata")
write.csv(e13_pDataC,file = file.path(data_dir,"durian_data","GuptaE13SC_pDataC.csv"))
write.csv(e14_pDataC,file = file.path(data_dir,"durian_data","GuptaE14SC_pDataC.csv"))

print("writing sc data(dfs)")
write.csv(e13$data_embryonic$data,file = file.path(data_dir,"durian_data","GuptaE13SC_C.csv"))
write.csv(e14$data_embryonic$data,file = file.path(data_dir,"durian_data","GuptaE14SC_C.csv"))

# #############################################################################
# # write cpm to csv
# #############################################################################

write.csv(e13_pDataC,file = file.path(data_dir,"durian_data","GuptaE13SC.cpm_pDataC.csv"))
write.csv(e14_pDataC,file = file.path(data_dir,"durian_data","GuptaE14SC.cpm_pDataC.csv"))

write.csv(edgeR::cpm(e13$data_embryonic$data),file = file.path(data_dir,"durian_data","GuptaE13SC.cpm_C.csv"))
write.csv(edgeR::cpm(e14$data_embryonic$data),file = file.path(data_dir,"durian_data","GuptaE14SC.cpm_C.csv"))

write.csv(edgeR::cpm(df_nondup),file = file.path(data_dir,"durian_data","BiggsBulk.cpm_T.csv"))

# #############################################################################
# # thresholding
# #############################################################################

# # use simple thresholding for genes
# BaronSC.DM.islet05_C = subsetsc(scremoutlier(BaronSC.DM.islet_C),generate=0.05,return_obj=TRUE,nsd=3)
# BaronSC.H.islet05_C = subsetsc(scremoutlier(BaronSC.H.islet_C),generate=0.05,return_obj=TRUE,nsd=3)
# BaronSC.DM.islet05_pDataC = BaronSC.DM.islet_pDataC[colnames(BaronSC.DM.islet05_C),]
# BaronSC.H.islet05_pDataC = BaronSC.H.islet_pDataC[colnames(BaronSC.H.islet05_C),]
# comgenes = intersect(rownames(BaronSC.H.islet05_C),rownames(BaronSC.DM.islet05_C))
# # save
# write.csv(BaronSC.DM.islet05_pDataC,file = file.path(data_dir,"BaronSC.DM.islet05_pDataC.csv"))
# write.csv(BaronSC.H.islet05_pDataC,file = file.path(data_dir,"BaronSC.H.islet05_pDataC.csv"))
# write.csv(BaronSC.H.islet05_C[comgenes,],file = file.path(data_dir,"BaronSC.DM.islet05_C.csv"))
# write.csv(BaronSC.DM.islet05_C[comgenes,],file = file.path(data_dir,"BaronSC.DM.islet05_C.csv"))

# # use seurat vst thresholding for genes
seur13 = CreateSeuratObject(counts=e13$data_embryonic$data)
seur14 = CreateSeuratObject(counts=e14$data_embryonic$data)
seur13 = NormalizeData(seur13)
seur14 = NormalizeData(seur14)
genes13 = VariableFeatures(FindVariableFeatures(seur13, selection.method = "vst", nfeatures = 5000))
genes14 = VariableFeatures(FindVariableFeatures(seur14, selection.method = "vst", nfeatures = 5000))
comgenes = intersect(genes13,genes14)
comgenes = intersect(comgenes,rownames(df_nondup))
e13.VST_C = subsetsc(e13$data_embryonic$data,geneids=comgenes,return_obj=TRUE,nsd=3)
e14.VST_C = subsetsc(e14$data_embryonic$data,geneids=comgenes,return_obj=TRUE,nsd=3)
e13.VST_pDataC = e13_pDataC[colnames(e13.VST_C),]
e14.VST_pDataC = e14_pDataC[colnames(e14.VST_C),]
# # save
write.csv(e13.VST_pDataC,file = file.path(data_dir,"durian_data","GuptaE13SC.VST_pDataC.csv"))
write.csv(e14.VST_pDataC,file = file.path(data_dir,"durian_data","GuptaE14SC.VST_pDataC.csv"))
write.csv(edgeR::cpm(e13.VST_C[comgenes,]),file = file.path(data_dir,"durian_data","GuptaE13SC.VST_C.csv"))
write.csv(edgeR::cpm(e14.VST_C[comgenes,]),file = file.path(data_dir,"durian_data","GuptaE14SC.VST_C.csv"))
write.csv(edgeR::cpm(df_nondup[comgenes,]),file = file.path(data_dir,"durian_data","BiggsBulk.VST_T.csv"))

# subsample 1000 cells
target = 1000
set.seed(1)
scdata = as.data.frame(t(read.csv("slurm/Gupta/durian_data/GuptaE13SC.VST_C.csv",row.names=1)))
sampn = min(nrow(scdata),target)
GuptaE13SC.VST1K05_C = as.data.frame(t(scdata %>% sample_n(sampn)))
GuptaE13SC.VST1K05_C = subsetsc(scremoutlier(GuptaE13SC.VST1K05_C),generate=0.05,return_obj=TRUE,nsd=3)
GuptaE13SC.VST1K05_pDataC = read.csv("slurm/Gupta/durian_data/GuptaE13SC.VST_pDataC.csv",row.names=1)[colnames(GuptaE13SC.VST1K05_C),]

write.csv(GuptaE13SC.VST1K05_C,file="slurm/Gupta/durian_data/GuptaE13SC.VST1K05_C.csv")
write.csv(GuptaE13SC.VST1K05_pDataC,file="slurm/Gupta/durian_data/GuptaE13SC.VST1K05_pDataC.csv")