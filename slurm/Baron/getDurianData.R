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

data_dir = "/dfs5/bio/mkarikom/temp/DURIAN/slurm/Baron/raw_data"
durian_data_dir = "/dfs5/bio/mkarikom/temp/DURIAN/slurm/Baron/durian_data"
dir.create(data_dir,recursive=TRUE)
dir.create(durian_data_dir,recursive=TRUE)
source("/dfs5/bio/mkarikom/temp/DURIAN/slurm/scrabble_helper_functions/library_scrabble.R")
celltypes = ("alpha-beta-gamma-delta-epsilon")
#############################################################################
# Segerstolpe 2016 (https://doi.org/10.1016/j.cmet.2016.08.020) 
# bulk (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5060/) 
#############################################################################
destdir=file.path(data_dir,"E-MTAB-5060")
dir.create(destdir,recursive = TRUE);oldwd=getwd()
E_MTAB_5060 = getAE("E-MTAB-5060", type = "full",path=destdir)
SegerstolpeBulk_fn = file.path(data_dir,"E-MTAB-5060","pancreas_refseq_rpkms_wholeislet.txt")
SegerstolpeBulk_metadata_fn = file.path(data_dir,"E-MTAB-5060","E-MTAB-5060.sdrf.txt")
SegerstolpeBulk = read.csv(SegerstolpeBulk_fn,sep='\t',header=FALSE)
Segerstolpe_cnames = as.character(unname(SegerstolpeBulk[1,2:8]))
Segerstolpe_hgnc = SegerstolpeBulk[2:dim(SegerstolpeBulk)[1],1]
Segerstolpe_refseq = SegerstolpeBulk[2:dim(SegerstolpeBulk)[1],2]
Segerstolpe_tpkbm = SegerstolpeBulk[2:dim(SegerstolpeBulk)[1],3:9]
Segerstolpe_reads = SegerstolpeBulk[2:dim(SegerstolpeBulk)[1],10:dim(SegerstolpeBulk)[2]]
SegerstolpeBulk = data.frame(matrix(0, length(unique(Segerstolpe_hgnc)), length(Segerstolpe_cnames)))
colnames(SegerstolpeBulk) = Segerstolpe_cnames
rownames(SegerstolpeBulk) = sort(unique(Segerstolpe_hgnc))

# find the mean of all duplicated hgnc symbols
for(i in rownames(SegerstolpeBulk)){
  # rowind = match(Segerstolpe_hgnc,i)[!is.na(match(Segerstolpe_hgnc,i))]
  rowind = which(!is.na(match(Segerstolpe_hgnc,i)))
  SegerstolpeBulk[i,] <- ceiling(colMeans(Segerstolpe_reads[rowind,]))
}


SegerstolpeBulk_metadata_fn = file.path(data_dir,"E-MTAB-5060","E-MTAB-5060.sdrf.txt")
SegerstolpeBulkMeta = read.csv(SegerstolpeBulk_metadata_fn,sep='\t',header=TRUE)
metacnames = c("Source.Name",
                "Characteristics.individual.",
                "Characteristics.disease.",
                "Characteristics.sex.",
                "Characteristics.age.",
                "Characteristics.body.mass.index.")
SegerstolpeBulkMeta = distinct(SegerstolpeBulkMeta[,metacnames])
colnames(SegerstolpeBulkMeta) = c("sampleID","patientID","disease","sex","age","bmi")


#############################################################################
# Baron 2016 (https://doi.org/10.1016/j.cels.2016.08.011) 
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133) 
#############################################################################
# get SC metadata by hand
BaronSCMeta <- data.frame(species = character(0),
                          sampleID = character(0),
                          age = integer(0),
                          sex = character(0),
                          bmi = double(0),
                          t2diabetes = integer(0))
BaronSCMeta[1,] <- c("human","GSM2230757",17,"M",21.5,0)
BaronSCMeta[2,] <- c("human","GSM2230758",51,"F",21.1,0)
BaronSCMeta[3,] <- c("human","GSM2230759",38,"M",27.5,0)
BaronSCMeta[4,] <- c("human","GSM2230760",59,"F",29.9,1)


destdir=file.path(data_dir,"GSE84133")
dir.create(destdir,recursive = TRUE)
getGEOSuppFiles(GEO="GSE84133", makeDirectory = FALSE, baseDir = destdir,
                fetch_files = TRUE, filter_regex = NULL)
untar(tarfile=file.path(data_dir,"GSE84133","GSE84133_RAW.tar"),exdir=file.path(data_dir,"GSE84133","RAW"))

rawfiles = sort(list.files(file.path(data_dir,"GSE84133","RAW"))) # first 4 are human
df = read.csv(file.path(data_dir,"GSE84133","RAW",rawfiles[1]))
gsm = strsplit(rawfiles[1],"_")[[1]][1]
namesdf = c("sampleID","cellID",colnames(df)[2:dim(df)[2]])  
Baron_raw = cbind(gsm,df)
colnames(Baron_raw) = namesdf
for(i in 2:4){
  cat("\n getting raw file ",i,"\n")
  df = read.csv(file.path(data_dir,"GSE84133","RAW",rawfiles[i]))
  gsm = strsplit(rawfiles[i],"_")[[1]][1]
  namesdf = c("sampleID","cellID",colnames(df)[2:dim(df)[2]])  
  new_raw = cbind(gsm,df)
  colnames(new_raw) = namesdf
  Baron_raw = rbind(Baron_raw,new_raw)
}

BaronSC = Baron_raw[,5:dim(Baron_raw)[2]]
rownames(BaronSC) = Baron_raw[,"cellID"]
BaronSC_C = t(BaronSC)
BaronSC_pDataC = data.frame(cellID=Baron_raw$cellID,cellType=Baron_raw$assigned_cluster,sampleID=Baron_raw$sampleID)
rownames(BaronSC_pDataC) = BaronSC_pDataC$cellID


BaronSCMeta = merge(x = BaronSCMeta, y = BaronSC_pDataC, on = "sampleID")
rownames(BaronSCMeta) = BaronSCMeta$cellID
BaronSCMeta = BaronSCMeta[,which(names(BaronSCMeta)!="cellID")]


#############################################################################
# separate disease and healthy
#############################################################################
SegerstolpeBulk.H_T = SegerstolpeBulk[SegerstolpeBulkMeta$sampleID[SegerstolpeBulkMeta$disease=="normal"]]
SegerstolpeBulk.DM_T = SegerstolpeBulk[SegerstolpeBulkMeta$sampleID[SegerstolpeBulkMeta$disease!="normal"]]

BaronSC.DM_pDataC = BaronSC_pDataC[rownames(BaronSCMeta[BaronSCMeta$t2diabetes=="1",]),]
BaronSC.H_pDataC = BaronSC_pDataC[rownames(BaronSCMeta[BaronSCMeta$t2diabetes=="0",]),]

BaronSC.DM_C = BaronSC_C[,rownames(BaronSCMeta[BaronSCMeta$t2diabetes=="1",])]
BaronSC.H_C = BaronSC_C[,rownames(BaronSCMeta[BaronSCMeta$t2diabetes=="0",])]

#############################################################################
# write all data to csv
#############################################################################
write.csv(SegerstolpeBulkMeta,file = file.path(data_dir,"SegerstolpeBulkMeta.csv"))
write.csv(SegerstolpeBulk,file = file.path(data_dir,"SegerstolpeBulk.csv"))
write.csv(BaronSCMeta,file = file.path(data_dir,"BaronSCMeta.csv"))
write.csv(BaronSC,file = file.path(data_dir,"BaronSC.csv"))

write.csv(SegerstolpeBulk.H_T,file = file.path(durian_data_dir,"SegerstolpeBulk.H_T.csv"))
write.csv(SegerstolpeBulk.DM_T,file = file.path(durian_data_dir,"SegerstolpeBulk.DM_T.csv"))

write.csv(BaronSC.DM_pDataC,file = file.path(durian_data_dir,"BaronSC.DM_pDataC.csv"))
write.csv(BaronSC.H_pDataC,file = file.path(durian_data_dir,"BaronSC.H_pDataC.csv"))

write.csv(BaronSC.DM_C,file = file.path(durian_data_dir,"BaronSC.DM_C.csv"))
write.csv(BaronSC.H_C,file = file.path(durian_data_dir,"BaronSC.H_C.csv"))

write.csv(edgeR::cpm(SegerstolpeBulk.H_T),file = file.path(durian_data_dir,"SegerstolpeBulk.H.cpm_T.csv"))
write.csv(edgeR::cpm(SegerstolpeBulk.DM_T),file = file.path(durian_data_dir,"SegerstolpeBulk.DM.cpm_T.csv"))

write.csv(BaronSC.DM_pDataC,file = file.path(durian_data_dir,"BaronSC.DM.cpm_pDataC.csv"))
write.csv(BaronSC.H_pDataC,file = file.path(durian_data_dir,"BaronSC.H.cpm_pDataC.csv"))

write.csv(edgeR::cpm(BaronSC.DM_C),file = file.path(durian_data_dir,"BaronSC.DM.cpm_C.csv"))
write.csv(edgeR::cpm(BaronSC.H_C),file = file.path(durian_data_dir,"BaronSC.H.cpm_C.csv"))


#############################################################################
# islet cells
#############################################################################

BaronSC.DM.islet_pDataC = BaronSC.DM_pDataC %>% dplyr::filter(cellType %in% strsplit(celltypes,"-")[[1]])
BaronSC.H.islet_pDataC = BaronSC.H_pDataC %>% dplyr::filter(cellType %in% strsplit(celltypes,"-")[[1]])
BaronSC.DM.islet_C = BaronSC.DM_C[,BaronSC.DM.islet_pDataC$cellID]
BaronSC.H.islet_C = BaronSC.H_C[,BaronSC.H.islet_pDataC$cellID]

# use simple thresholding for genes
BaronSC.DM.islet05_C = subsetsc(scremoutlier(BaronSC.DM.islet_C),generate=0.05,return_obj=TRUE,nsd=3)
BaronSC.H.islet05_C = subsetsc(scremoutlier(BaronSC.H.islet_C),generate=0.05,return_obj=TRUE,nsd=3)
BaronSC.DM.islet05_pDataC = BaronSC.DM.islet_pDataC[colnames(BaronSC.DM.islet05_C),]
BaronSC.H.islet05_pDataC = BaronSC.H.islet_pDataC[colnames(BaronSC.H.islet05_C),]
comgenes = intersect(rownames(BaronSC.H.islet05_C),rownames(BaronSC.DM.islet05_C))
# save
write.csv(BaronSC.DM.islet05_pDataC,file = file.path(durian_data_dir,"BaronSC.DM.islet05_pDataC.csv"))
write.csv(BaronSC.H.islet05_pDataC,file = file.path(durian_data_dir,"BaronSC.H.islet05_pDataC.csv"))
write.csv(BaronSC.H.islet05_C[comgenes,],file = file.path(durian_data_dir,"BaronSC.DM.islet05_C.csv"))
write.csv(BaronSC.DM.islet05_C[comgenes,],file = file.path(durian_data_dir,"BaronSC.DM.islet05_C.csv"))

# use seurat vst thresholding for genes
seurDM = CreateSeuratObject(counts=BaronSC.DM.islet_C)
seurH = CreateSeuratObject(counts=BaronSC.H.islet_C)
seurDM = NormalizeData(seurDM)
seurH = NormalizeData(seurH)
genesDM = VariableFeatures(FindVariableFeatures(seurDM, selection.method = "vst", nfeatures = 5000))
genesH = VariableFeatures(FindVariableFeatures(seurH, selection.method = "vst", nfeatures = 5000))
comgenes = intersect(genesDM,genesH)
BaronSC.DM.isletVST_C = subsetsc(scremoutlier(BaronSC.DM.islet_C),geneids=comgenes,return_obj=TRUE,nsd=3)
BaronSC.H.isletVST_C = subsetsc(scremoutlier(BaronSC.H.islet_C),geneids=comgenes,return_obj=TRUE,nsd=3)
BaronSC.DM.isletVST_pDataC = BaronSC.DM.islet_pDataC[colnames(BaronSC.DM.isletVST_C),]
BaronSC.H.isletVST_pDataC = BaronSC.H.islet_pDataC[colnames(BaronSC.H.isletVST_C),]
# save
write.csv(BaronSC.DM.isletVST_pDataC,file = file.path(durian_data_dir,"BaronSC.DM.isletVST_pDataC.csv"))
write.csv(BaronSC.H.isletVST_pDataC,file = file.path(durian_data_dir,"BaronSC.H.isletVST_pDataC.csv"))
write.csv(BaronSC.H.isletVST_C[comgenes,],file = file.path(durian_data_dir,"BaronSC.H.isletVST_C.csv"))
write.csv(BaronSC.DM.isletVST_C[comgenes,],file = file.path(durian_data_dir,"BaronSC.DM.isletVST_C.csv"))

# subsample 1000 cells
target = 1000
set.seed(1)
scdata = as.data.frame(t(read.csv("slurm/Baron/durian_data/BaronSC.DM.isletVST_C.csv",row.names=1)))
sampn = min(nrow(scdata),target)
BaronSC.DM.isletVST1K05_C = as.data.frame(t(scdata %>% sample_n(sampn)))
BaronSC.DM.isletVST1K05_C = subsetsc(scremoutlier(BaronSC.DM.isletVST1K05_C),generate=0.05,return_obj=TRUE,nsd=3)
BaronSC.DM.isletVST1K05_pDataC = read.csv("slurm/Baron/durian_data/BaronSC.DM.isletVST_pDataC.csv",row.names=1)[colnames(BaronSC.DM.isletVST1K05_C),]

write.csv(BaronSC.DM.isletVST1K05_C,file="slurm/Baron/durian_data/BaronSC.DM.isletVST1K05_C.csv")
write.csv(BaronSC.DM.isletVST1K05_pDataC,file="slurm/Baron/durian_data/BaronSC.DM.isletVST1K05_pDataC.csv")

set.seed(1)
scdata = as.data.frame(t(read.csv("slurm/Baron/durian_data/BaronSC.H.isletVST_C.csv",row.names=1)))
sampn = min(nrow(scdata),target)
BaronSC.H.isletVST1K05_C = as.data.frame(t(scdata %>% sample_n(sampn)))
BaronSC.H.isletVST1K05_C = subsetsc(scremoutlier(BaronSC.H.isletVST1K05_C),generate=0.05,return_obj=TRUE,nsd=3)
BaronSC.H.isletVST1K05_pDataC = read.csv("slurm/Baron/durian_data/BaronSC.H.isletVST_pDataC.csv",row.names=1)[colnames(BaronSC.H.isletVST1K05_C),]

write.csv(BaronSC.H.isletVST1K05_C,file="slurm/Baron/durian_data/BaronSC.H.isletVST1K05_C.csv")
write.csv(BaronSC.H.isletVST1K05_pDataC,file="slurm/Baron/durian_data/BaronSC.H.isletVST1K05_pDataC.csv")
