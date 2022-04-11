# Code adapted from https://github.com/favilaco/deconv_benchmark:
# 
# Copyright 2019 Francisco Avila Cobos, 2021 Matt Karikomi
library(doParallel)
library(dplyr)
library(Matrix)
library(purrr)

################################################################
# arguments
################################################################

durianlib = Sys.getenv("DURIANLIB")
scdatadir = Sys.getenv("SCDATADIR")
scprefix = Sys.getenv("SCPREFIX")
datapath = Sys.getenv("DATAPATH")
pseudocount = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ngene = as.integer(Sys.getenv("NGENE"))
ncell = as.integer(Sys.getenv("NCELL"))
sim_dropout_rate = as.numeric(Sys.getenv("DROPOUTRATE"))
lambda = as.numeric(Sys.getenv("LAMBDA"))
nzthresh = as.numeric(Sys.getenv("NZTHRESH"))
trainrate = as.numeric(Sys.getenv("PBTRAINRATE"))

### Read data and metadata
print("reading data")
data = as.matrix(read.csv(file.path(scdatadir,paste0(scprefix,"C.csv")),row.names = 1))
full_phenoData = read.csv(file.path(scdatadir,paste0(scprefix,"pDataC.csv")),row.names = 1)
rownames(full_phenoData) = full_phenoData$cellID

if(dim(data)[1] > 0){
  print("data is read")
}else{
  print("problem reading data")
}
if(subsetcelltypes!=""){
  full_phenoData = full_phenoData %>% dplyr::filter(cellType %in% strsplit(subsetcelltypes,"-")[[1]])
  data = data[,full_phenoData$cellID]
}
################################################################
# QC
################################################################

# # First: cells with library size, mitochondrial or ribosomal content further than three MAD away were discarded
# filterCells <- function(filterParam){
#   cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
#   cellsToRemove
# }

# libSizes <- colSums(data)
# gene_names <- rownames(data)

# mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
# rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

# mtPercent <- colSums(data[mtID, ])/libSizes
# rbPercent <- colSums(data[rbID, ])/libSizes

# lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
#   unlist() %>% 
#   unique() -> cellsToRemove

# if(length(cellsToRemove) != 0){
#   data <- data[,-cellsToRemove]
#   full_phenoData <- full_phenoData[-cellsToRemove,]
# }
# Keep only "detectable" genes: at least 5% of cells (regardless of the group) have a read/UMI count different from 0
keep <- which(Matrix::rowSums(data > 0) >= round(nzthresh * ncol(data)))
data = data[keep,]

################################################################
# scale the data
################################################################
data = edgeR::cpm(data)

################################################################
# take a subset of genes
################################################################
set.seed(pseudocount)
if(dim(data)[1] > ngene){
  geneinds = sample(x = 1:dim(data)[1], size = ngene, replace = FALSE)
  data = data[geneinds,]
}
if(dim(data)[2] > ncell){
  cellinds = sample(x = 1:dim(data)[2], size = ncell, replace = FALSE)
  data = data[,cellinds]
  full_phenoData = full_phenoData[cellinds,]
}

################################################################
# train test split
################################################################
ncells = dim(data)[2]
ngenes = dim(data)[1]
splitsize = as.integer(floor(ncells*trainrate))

set.seed(pseudocount)
traininds = sample(x = 1:ncells,size = splitsize,replace = FALSE)
testinds = setdiff(1:ncells,traininds)
testing_phenoData = full_phenoData[testinds,]

training = data[,traininds]
testing = data[,testinds]
testing.orig = testing

################################################################
# use the sampleID of the training to make pseudobulk 
################################################################
train_phenoData = full_phenoData[traininds,]
bulkids = sort(unique(train_phenoData$sampleID))
T = matrix(0,nrow=ngenes,ncol=length(bulkids))
P = matrix(0,nrow=length(unique(full_phenoData$cellType)),ncol=length(bulkids))
rownames(T) = rownames(training)
colnames(T) = bulkids
rownames(P) = sort(unique(full_phenoData$cellType))
colnames(P) = bulkids
for(i in 1:length(bulkids)){
  T[,i] = rowMeans(training[,train_phenoData$sampleID == bulkids[i]])
  bulktable = table(train_phenoData$cellType[train_phenoData$sampleID == bulkids[i]])
  P[names(bulktable),i] = unname(bulktable)/length(which(train_phenoData$sampleID == bulkids[i]))
}
################################################################
# get celltype-specific dropout
################################################################
celltypes = sort(unique(full_phenoData$cellType))
for(i in 1:length(celltypes)){
  print(paste0("downsampling cell type ",i))
  ids = full_phenoData$cellID[which(full_phenoData$cellType == celltypes[i])]
  p_celltype = exp(-lambda*rowMeans(data[,ids])^2)
  testids = testing_phenoData$cellID[which(testing_phenoData$cellType == celltypes[i])]
  dropgenes = names(which(rbernoulli(n=1,p=p_celltype)))
  for(j in testids){
    testing[dropgenes,j] = 0
  }
}

################################################################
# save the data
################################################################

write.csv(T,file.path(datapath,"T.csv"))
write.csv(P,file.path(datapath,"trueP.csv"))
write.csv(testing.orig,file.path(datapath,"trueC.csv"))
write.csv(testing_phenoData,file.path(datapath,"pDataC.csv"))
write.csv(testing,file.path(datapath,"C.csv"))