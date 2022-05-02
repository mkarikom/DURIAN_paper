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
data_dir = file.path(project_dir,"slurm/He/durian_data_hpc")
sc_data_dir = file.path(project_dir,"slurm/He/sc_data")
quant_bulk_data_dir = file.path(project_dir,"slurm/He/salmon_output_sbatch")

sapply(list(data_dir,
            sc_data_dir,
            quant_bulk_data_dir),dir.create,recursive = TRUE);

#############################################################################
# Suarez-Farinas 2015 (https://doi.org/10.1016/j.jaci.2015.03.003) 
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65832) 
#############################################################################
# # run the non-sbatch script to download the data with Aspera
# note this was run locally on fuse-mounted dfs5 due to UDP error on hpc
# $ ./downloadSRA.sh

# # run the sbatch script to quantitate the single-ended, single-run fastq with Salmon
# $ ./runsalmon.sl

# construct the transcript ID -> gene ID map
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v95") # will output "retrieve record with 'object[["AH67950"]]'"
edb = ah[["AH67950"]]
mykeys <- head(keys(org.Hs.eg.db, keytype="ENSEMBL"))
ensid <- keys(edb)
mapping_ids <- AnnotationDbi::select(edb, 
            keys=ensid,
            columns=c("TXIDVERSION","GENEID"));

mapping_names <- AnnotationDbi::select(edb, 
            keys=ensid,
            columns=c("TXIDVERSION","GENENAME")) %>% dplyr::select(c("TXIDVERSION","GENENAME"));

# create a single matrix of transcript counts from the salmon output
bulk_samples = list.dirs(quant_bulk_data_dir,recursive=FALSE)

# create a single matrix of transcript counts from the salmon output
# use the gene-level quant from AnnotationDbi -> tximport
bulk_list_geneid = list()
bulk_list_genename = list()
bulk_list_tx = list()
counter = 0
for(i in bulk_samples){
      counter = counter + 1
      cat(paste0("\n processing bulk sample: ",counter,"\n"))
      print("\n")

      # gene level quant (gene ids)
      fn = file.path(i,"quant.sf")
      txi <- tximport(fn, type = "salmon", tx2gene = mapping_ids)
      geneids=rownames(txi$counts)
      txcounts=unname(txi$counts)
      tempquant = data.frame(geneId=geneids,counts=txcounts)
      pathlen = length(stringr::str_split(i,"/")[[1]])
      samplename = substr(stringr::str_split(i,"/")[[1]][pathlen],12,47)
      colnames(tempquant) = c("geneId",samplename)
      bulk_list_geneid[[counter]] = tempquant

      # gene level quant (gene names)
      fn = file.path(i,"quant.sf")
      txi <- tximport(fn, type = "salmon", tx2gene = mapping_names)
      genenames=rownames(txi$counts)
      txcounts=unname(txi$counts)
      tempquant = data.frame(geneName=genenames,counts=txcounts)
      colnames(tempquant) = c("geneName",samplename)
      bulk_list_genename[[counter]] = tempquant

      # tx level quant
      rawquant = read.table(fn, header = TRUE, sep = "", dec = ".")[,c("Name","NumReads")]
      colnames(rawquant) = c("txId",samplename)
      bulk_list_tx[[counter]] = rawquant
}
bulk_concat_geneid = bulk_list_geneid %>% purrr::reduce(inner_join,by='geneId')
bulk_concat_genename = bulk_list_genename %>% purrr::reduce(inner_join,by='geneName')
bulk_concat_tx = bulk_list_tx %>% purrr::reduce(inner_join,by='txId')

# round the gene symbol quant
bulk_concat_rounded = bulk_concat_genename
for(i in 2:dim(bulk_concat_genename)[2]){
bulk_concat_rounded[,i] = ceiling(bulk_concat_rounded[,i])
}
rownames(bulk_concat_rounded) = bulk_concat_rounded$geneName
bulk_concat_rounded = subset(bulk_concat_rounded,select = -c(geneName))
bulk_status = unname(sapply(names(bulk_concat_rounded),function(x){
            strsplit(x,"_")[[1]][3]
            }))

#############################################################################
# He 2020 (https://doi.org/10.1016/j.jaci.2020.01.042) 
# bulk (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147424) 
#############################################################################
# get SC data
LSfigshare = 25956518
NLfigshare =  25954199
system2("wget",args=c(paste0("-P ",sc_data_dir),paste0("https://ndownloader.figshare.com/files/",LSfigshare)))
system2("wget",args=c(paste0("-P ",sc_data_dir),paste0("https://ndownloader.figshare.com/files/",NLfigshare)))
sc_nl = readRDS(file.path(sc_data_dir,NLfigshare))
sc_ls = readRDS(file.path(sc_data_dir,LSfigshare))

barcodes_nl = colnames(sc_nl@data)
labels_nl = sc_nl@idents
donors_nl = unname(sapply(barcodes_nl,function(x){
stringr::str_split(x,"_")[[1]][1]
}))
pDataC_nl = data.frame(cellID=barcodes_nl,cellType=labels_nl,sampleID=donors_nl)
rownames(pDataC_nl) = pDataC_nl$cellID

barcodes_ls = colnames(sc_ls@data)
labels_ls = sc_ls@idents
donors_ls = unname(sapply(barcodes_ls,function(x){
stringr::str_split(x,"_")[[1]][1]
}))
pDataC_ls = data.frame(cellID=barcodes_ls,cellType=labels_ls,sampleID=donors_ls)
rownames(pDataC_ls) = pDataC_ls$cellID

write.csv(sc_ls@data,file = file.path(data_dir,"HeLS_C.csv"))
write.csv(sc_nl@data,file = file.path(data_dir,"HeNL_C.csv"))

write.csv(pDataC_ls,file = file.path(data_dir,"HeLS_pDataC.csv"))
write.csv(pDataC_nl,file = file.path(data_dir,"HeNL_pDataC.csv"))

# save to durian dir
write.csv(bulk_concat_rounded,file = file.path(data_dir,"Suarez_T.csv"))
write.csv(bulk_concat_rounded[,which(bulk_status=="NL")],file = file.path(data_dir,"SuarezNL_T.csv"))
write.csv(bulk_concat_rounded[,which(bulk_status=="LS")],file = file.path(data_dir,"SuarezLS_T.csv"))

# # use seurat vst thresholding for genes
seurLS = CreateSeuratObject(counts=sc_ls@data)
seurNL = CreateSeuratObject(counts=sc_nl@data)

seurLS = NormalizeData(seurLS)
seurNL = NormalizeData(seurNL)

genesLS = VariableFeatures(FindVariableFeatures(seurLS, selection.method = "vst", nfeatures = 5000))
genesNL = VariableFeatures(FindVariableFeatures(seurNL, selection.method = "vst", nfeatures = 5000))

commgenes = intersect(intersect(genesLS,rownames(bulk_concat_rounded)),genesNL)
ls.VST_C = subsetsc(as.matrix(sc_ls@data),geneids=commgenes,return_obj=TRUE,nsd=3)
nl.VST_C = subsetsc(as.matrix(sc_nl@data),geneids=commgenes,return_obj=TRUE,nsd=3)

ls.VST_pDataC = pDataC_ls[colnames(ls.VST_C),]
nl.VST_pDataC = pDataC_nl[colnames(nl.VST_C),]

write.csv(ls.VST_pDataC,file = file.path(data_dir,"HeLS.VST_pDataC.csv"))
write.csv(nl.VST_pDataC,file = file.path(data_dir,"HeNL.VST_pDataC.csv"))

write.csv(ls.VST_C,file = file.path(data_dir,"HeLS.VST_C.csv"))
write.csv(nl.VST_C,file = file.path(data_dir,"HeNL.VST_C.csv"))

write.csv(bulk_concat_rounded[commgenes,which(bulk_status=="LS")],file = file.path(data_dir,"SuarezLS.VST_T.csv"))
write.csv(bulk_concat_rounded[commgenes,which(bulk_status=="NL")],file = file.path(data_dir,"SuarezNL.VST_T.csv"))

write.csv(bulk_concat_rounded[commgenes,which(bulk_status=="LS")[1:5]],file = file.path(data_dir,"SuarezLS1to5.VST_T.csv"))
write.csv(bulk_concat_rounded[commgenes,which(bulk_status=="NL")[1:5]],file = file.path(data_dir,"SuarezNL1to5.VST_T.csv"))