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
suoqin_data = file.path(project_dir,"slurm/Gupta/suoqin_data")
data_dir = file.path(project_dir,"slurm/Gupta/durian_data_hpc")
sc_data_dir = file.path(project_dir,"slurm/Gupta/sc_data")
sapply(list(data_dir,
            sc_data_dir),dir.create,recursive = TRUE);

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
write.csv(df_nondup,file = file.path(data_dir,"BiggsBulk_T.csv"))

#############################################################################
# Gupta 2019 (https://doi.org/10.1016/j.devcel.2018.11.032) 
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122043)
#############################################################################
# this data was provided by Suoqin Jin (private communication)

e13_C = read.csv(file.path(suoqin_data,"GuptaE13SC_C.csv.gz"),row.names=1)
e13_pDataC = read.csv(file.path(suoqin_data,"GuptaE13SC_pDataC.csv.gz"),row.names=1)

# #############################################################################
# # thresholding
# #############################################################################

# # use seurat vst thresholding for genes
seur13 = CreateSeuratObject(counts=e13_C)
seur13 = NormalizeData(seur13)
genes13 = VariableFeatures(FindVariableFeatures(seur13, selection.method = "vst", nfeatures = 5000))

commgenes = intersect(genes13,rownames(df_nondup))
e13.VST_C = subsetsc(e13_C,geneids=commgenes,return_obj=TRUE,nsd=3)
e13.VST_pDataC = e13_pDataC[colnames(e13.VST_C),]

write.csv(e13.VST_pDataC,file = file.path(data_dir,"GuptaE13SC.VST_pDataC.csv"))
write.csv(e13.VST_C,file = file.path(data_dir,"GuptaE13SC.VST_C.csv"))
write.csv(df_nondup[commgenes,],file = file.path(data_dir,"BiggsBulk.VST_T.csv"))