# compare NL single cell to NL bulk (treatment = bulk->single-cell)
# form pseudobulk from each single-cell donor
# 1) use edgeR for DEG
# 2) compare within vs between groups using bigPint::plotSM,bigPint::plotPCP:
#     1) sc NL vs sc LS
#     2) bulk NL vs sc NL
#     3) bulk LS vs sc LS

library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(GGally)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(RColorBrewer)
library(GEOquery)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(scran)
library(Seurat)
library(apeglm)
library(pheatmap)
library(ggh4x)
library(DURIAN)

print("setup environment")
projdir = "/share/crsp/lab/cellfate/mkarikom/DURIAN_paper_clean"
durian_data_dir = file.path(projdir,"slurm/He/durian_data")
basedir = file.path("slurm/He/durian_data")
plotdir = file.path(basedir,"plots")
backupdir = file.path(basedir,"backup")
raw_data_dir = file.path(projdir,"slurm/He/sc_data_raw")
dir.create(raw_data_dir,recursive = TRUE)

dir.create(plotdir,recursive = TRUE)
dir.create(backupdir, recursive = TRUE)
mincells = 500 # the minimum number of cells to make pseudobulk from He

sessionInfo() %>% capture.output(file=file.path(backupdir,"session_info.txt"));

seur_min_cells = 10
seur_min_features = 200
deseq_min_counts_gene = 10
durian_max_cell_opts = c(2500,5000) # the number of cells to subsample for signaling analysis
############################################################################################
# Check the QC data to see if the high number of dropouts in one sample
##########################################################################################

# LS
print("load suarez ls")
Suarez_ls = read.csv(file.path(durian_data_dir,"SuarezLS_T.csv"),row.names=1)
# colnames(Suarez_ls) = paste0("Bulk_ls.",1:(dim(Suarez_ls)[2]))

print("load jin ls")
He_ls_pDataC = read.csv(file.path(durian_data_dir,"HeLS_pDataC.csv"),row.names=1)
#
He_ls = read.csv(file.path(durian_data_dir,"HeLS_C.csv"),row.names=1)

reads = He_ls

# create pseudobulk from ls sc data (one pseudobulk per sample)
He_ls_ids = sort(unique(He_ls_pDataC$sampleID))
pseudo_reads = as.data.frame(matrix(0, dim(reads)[1], length(He_ls_ids)))
rownames(pseudo_reads) = rownames(reads)
#
for(i in 1:length(He_ls_ids)){
  count = 0
  metarows = which(He_ls_pDataC$sampleID==He_ls_ids[i])
  cellids = intersect(He_ls_pDataC$cellID[metarows],colnames(reads))
  pseudo_reads[,i] = rowSums(reads[,cellids])
}
colnames(pseudo_reads) = paste0("PBulk_ls.",1:length(He_ls_ids))
He_ls_pseudo = pseudo_reads
#
# NL
print("load suarez nl")
Suarez_nl = read.csv(file.path(durian_data_dir,"SuarezNL_T.csv"),row.names=1)
# colnames(Suarez_nl) = paste0("Bulk_nl.",1:(dim(Suarez_nl)[2]))

print("load jin nl")
He_nl_pDataC = read.csv(file.path(durian_data_dir,"HeNL_pDataC.csv"),row.names=1)
#
He_nl = read.csv(file.path(durian_data_dir,"HeNL_C.csv"),row.names=1)

reads = He_nl

# create pseudobulk from nl sc data (one pseudobulk per sample)
He_nl_ids = sort(unique(He_nl_pDataC$sampleID))
pseudo_reads = as.data.frame(matrix(0, dim(reads)[1], length(He_nl_ids)))
rownames(pseudo_reads) = rownames(reads)
#
for(i in 1:length(He_nl_ids)){
  count = 0
  metarows = which(He_nl_pDataC$sampleID==He_nl_ids[i])
  cellids = intersect(He_nl_pDataC$cellID[metarows],colnames(reads))
  pseudo_reads[,i] = rowSums(reads[,cellids])
}
colnames(pseudo_reads) = paste0("PBulk_nl.",1:length(He_nl_ids))
He_nl_pseudo = pseudo_reads
#
#
##########################################################################################
# Check the RAW GEO data to see if the high number of dropouts in one sample
##########################################################################################

# destdir=file.path(raw_data_dir,"GSE147424")
# dir.create(destdir,recursive = TRUE)
getGEOSuppFiles(GEO="GSE147424", makeDirectory = FALSE, baseDir = raw_data_dir,
                fetch_files = TRUE, filter_regex = NULL)
untar(tarfile=file.path(raw_data_dir,"GSE147424_RAW.tar"),exdir=file.path(raw_data_dir,"gsm"))

rawfiles = sort(list.files(file.path(raw_data_dir,"gsm"))) # first 4 are human
df = read.csv(file.path(raw_data_dir,"gsm",rawfiles[1]))
gsm = strsplit(rawfiles[1],"_")[[1]][1]
namesdf = c("sampleID","cellID",colnames(df)[2:dim(df)[2]])  
Baron_raw = cbind(gsm,df)
colnames(Baron_raw) = namesdf

print("load he raw")
rawdata_dir = file.path(raw_data_dir,"gsm")
list.files(rawdata_dir)
donors_orig = NULL
dnames = c()
files = list.files(rawdata_dir)
seur_list = list()
count = 0
for(i in 1:length(files)){
  cat("\n getting donor ",i,"\n")
  donor = read.csv(file.path(rawdata_dir,files[i]))
  genes = donor$X
  donor = donor[,-1]
  rownames(donor) = genes
  if(dim(donor)[2] > mincells){
    count = count+1
    # create the seurat object and do qc
    seur_list[[count]] = CreateSeuratObject(counts = donor, project = "He", min.cells = seur_min_cells, min.features = seur_min_features)
    seur_list[[count]][["percent.mt"]] = PercentageFeatureSet(seur_list[[count]], pattern = "^MT-")
    seur_list[[count]] = subset(seur_list[[count]], subset = nFeature_RNA > 0 & nFeature_RNA < 2500 & percent.mt < 5)
        
    seur.sce = as.SingleCellExperiment(seur_list[[count]])
    seur.sce = computeSumFactors(seur.sce)
    seur_list[[count]] = seur.sce
    scaling = t(replicate(dim(seur.sce@assays@data$counts)[1], seur.sce$sizeFactor))
    donor_filt = as.data.frame(seur.sce@assays@data$counts * scaling)
    genes = rownames(donor_filt)
    # concatenate the pseudobulk vectors
    dnames[count] = strsplit(colnames(donor_filt)[2],"_")[[1]][1]
    if(is.null(donors_orig)){
      donors_orig = data.frame(gene=genes,x=rowSums(donor_filt))
    }else{
      donors_orig = merge(donors_orig,data.frame(gene=genes,x=rowSums(donor_filt)),by="gene")
    }
    
  }
}
colnames(donors_orig) = c("gene",dnames)

# get the metadata
filtered_donors = colnames(donors_orig)[-1]


gse147424 <- getGEO('GSE147424',GSEMatrix=TRUE)
gse=gse147424$GSE147424_series_matrix.txt.gz@phenoData@data
samp = unname(sapply(gse$title,function(x) strsplit(x,"_")[[1]][1]))
status = unname(sapply(gse$title,function(x) strsplit(x,"_")[[1]][2]))[match(filtered_donors,samp)]
status_type = sort(unique(status))


donor_list = list()
for(i in 1:length(status_type)){
  inds = which(status==status_type[i])+1
  donor_list[[status_type[i]]] = donors_orig[,c(1,inds)]
  colnames(donor_list[[status_type[i]]]) = c("gene",paste0("PBulk_raw_",tolower(status_type[i]),".",1:length(inds)))
}

Suarez_ls$gene = rownames(Suarez_ls)
Suarez_nl$gene = rownames(Suarez_nl)
He_ls_pseudo$gene = rownames(He_ls_pseudo)
He_nl_pseudo$gene = rownames(He_nl_pseudo)

bulkdata_orig=merge(Suarez_ls, Suarez_nl,by="gene")
bulkdata_orig=merge(bulkdata_orig, He_ls_pseudo,by="gene")
bulkdata_orig=merge(bulkdata_orig, He_nl_pseudo,by="gene")
bulkdata_orig=merge(bulkdata_orig, donor_list$H,by="gene")
bulkdata_orig=merge(bulkdata_orig, donor_list$LS,by="gene")
bulkdata_orig=merge(bulkdata_orig, donor_list$NL,by="gene")

# set up the differential expression
sample_group = factor(c(
  rep("Suarez_ls",(dim(Suarez_ls)[2]-1)),
  rep("Suarez_nl",(dim(Suarez_nl)[2]-1)),
  rep("He_ls",(dim(He_ls_pseudo)[2]-1)),
  rep("He_nl",(dim(He_nl_pseudo)[2]-1)),
  rep("He_h_raw",(dim(donor_list[["H"]])[2]-1)),
  rep("He_ls_raw",(dim(donor_list[["LS"]])[2]-1)),
  rep("He_nl_raw",(dim(donor_list[["NL"]])[2]-1))))
study_group = factor(c(
  rep("Suarez",(dim(Suarez_ls)[2]-1)),
  rep("Suarez",(dim(Suarez_nl)[2]-1)),
  rep("Jin",(dim(He_ls_pseudo)[2]-1)),
  rep("Jin",(dim(He_nl_pseudo)[2]-1)),
  rep("He",(dim(donor_list[["H"]])[2]-1)),
  rep("He",(dim(donor_list[["LS"]])[2]-1)),
  rep("He",(dim(donor_list[["NL"]])[2]-1))))
modality_group = factor(c(
  rep("Bulk",(dim(Suarez_ls)[2]-1)),
  rep("Bulk",(dim(Suarez_nl)[2]-1)),
  rep("Pseudobulk",(dim(He_ls_pseudo)[2]-1)),
  rep("Pseudobulk",(dim(He_nl_pseudo)[2]-1)),
  rep("Pseudobulk",(dim(donor_list[["H"]])[2]-1)),
  rep("Pseudobulk",(dim(donor_list[["LS"]])[2]-1)),
  rep("Pseudobulk",(dim(donor_list[["NL"]])[2]-1))))
disease_group = factor(c(
  rep("Lesion",(dim(Suarez_ls)[2]-1)),
  rep("Non_Lesion",(dim(Suarez_nl)[2]-1)),
  rep("Lesion",(dim(He_ls_pseudo)[2]-1)),
  rep("Non_Lesion",(dim(He_nl_pseudo)[2]-1)),
  rep("Healthy",(dim(donor_list[["H"]])[2]-1)),
  rep("Lesion",(dim(donor_list[["LS"]])[2]-1)),
  rep("Non_Lesion",(dim(donor_list[["NL"]])[2]-1))))

genes = bulkdata_orig$gene
sampleid = as.factor(colnames(bulkdata_orig[,-1]))
bulkdata_orig = t(bulkdata_orig[,-1])
colnames(bulkdata_orig) = genes

ids = data.frame(ids=sampleid)
rownames(ids) = sampleid
bulkdata_orig = merge(ids,bulkdata_orig,by=0)
bulkdata_orig = bulkdata_orig[,-1]

sample_group = data.frame(ids=sampleid,sample_group=sample_group)
study_group = data.frame(ids=sampleid,study_group=study_group)
modality_group = data.frame(ids=sampleid,modality_group=modality_group)
disease_group = data.frame(ids=sampleid,disease_group=disease_group)

bulkdata_orig = merge(sample_group,bulkdata_orig,by="ids")
bulkdata_orig = merge(study_group,bulkdata_orig,by="ids")
bulkdata_orig = merge(modality_group,bulkdata_orig,by="ids")
bulkdata_orig = merge(disease_group,bulkdata_orig,by="ids")

##########################################################################
# remove any genes expressed in less than 10 bulk/pbulk
##########################################################################
minbulk = max(colSums(bulkdata_orig[,genes]>0))
metacol = setdiff(colnames(bulkdata_orig),genes)
bulkdata_orig[,genes] = round(bulkdata_orig[,genes])
indbulk = which(colSums(bulkdata_orig[,6:dim(bulkdata_orig)[2]]>0) < minbulk)
bulkdata_orig = bulkdata_orig[,which(colnames(bulkdata_orig)!=genes[indbulk])]

##########################################################################
# reshape and plot
##########################################################################

melted = reshape2::melt(bulkdata_orig,ids=c(ids,sample_group,study_group,modality_group,disease_group))

# log all the couunts
melted$value = log(melted$value+1)

bp = ggplot(melted, aes(x=ids, y=value, fill=study_group)) +
  geom_boxplot(aes(fill=study_group)) +
  facet_nested(. ~ modality_group + disease_group, drop=TRUE, scales = "free", space = "free") +
  theme_bw() +
  theme(
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank()) +
  labs(title="Un-Normalized Log Pseudobulk and Bulk Expression", x = "bulk/pseudobulk donor", y= "log expression")
bp
ggsave(plot=bp,file=file.path(plotdir,"logdata.pdf"))

modality_submetrics = list()
disease_submetrics = list()

##########################################################################
# disease
##########################################################################
vartitle = "disease"
submetrics = modality_submetrics
bulkdata = bulkdata_orig[which(bulkdata_orig$study_group!="Jin"),] # remove the cellchat subset of the He data
bulkdata = bulkdata[which(bulkdata$disease_group!="Healthy"),] # only use tissue match unaffected for controls
group = droplevels(bulkdata$disease)

# save to the durian data folder
lsbulk = bulkdata[intersect(which(bulkdata$modality_group=="Bulk"), which(bulkdata$disease_group=="Lesion")),]
lsbulkids = lsbulk$ids
lsbulk = t(lsbulk[,-c(1:5)])
colnames(lsbulk) = lsbulkids

nlbulk = bulkdata[intersect(which(bulkdata$modality_group=="Bulk"), which(bulkdata$disease_group=="Non_Lesion")),]
nlbulkids = nlbulk$ids
nlbulk = t(nlbulk[,-c(1:5)])
colnames(nlbulk) = nlbulkids

labels = c()
count0 = 0
for(i in 1:length(levels(group))){
  count = 0
  for(j in 1:length(group)){
    if(group[j] == levels(group)[i]){
      count = count + 1
      count0 = count0 + 1
      labels[count0] = paste0(levels(group)[i],".",count)
    }
  }
}
#
# input data for deseq2
flipdata = ceiling(t(bulkdata[,-c(1:5)]))
colnames(flipdata) = labels
flipsamples = data.frame(group=group,modality=bulkdata$modality_group)
rownames(flipsamples) = colnames(flipdata)
#
dds = DESeqDataSetFromMatrix(countData = flipdata,
                            colData = flipsamples,
                            design = ~ modality + group)
# filter low count genes
keep = rowSums(counts(dds)) >= deseq_min_counts_gene
dds = dds[keep,]

# differential expression analysis
dds = DESeq(dds)

res = results(dds,name=resultsNames(dds)[3]) # get the genes for lesion vs non-lesion (the first term is intercept, the second is the modality)

# calculate shrinkage to rank genes
resLFC = lfcShrink(dds, coef=resultsNames(dds)[3], type="apeglm")

degs = as.data.frame(resLFC) %>% dplyr::filter(padj <= 1e-2)

# log fold change plot
pdf(file.path(plotdir,"lfc_shrink.pdf"))
plotMA(resLFC, ylim=c(-2,2),alpha=1e-2)
dev.off()

##########################################################################
# reshape and plot after rescaling
##########################################################################
`%notin%` = Negate(`%in%`)
dds_genes = which(colnames(bulkdata) %in% rownames(degs))

bulkdata_rescaled = bulkdata[,c(which(colnames(bulkdata) %notin% genes),dds_genes)]
normcounts = t(counts(dds, normalized=TRUE))

dds_gene_name = colnames(bulkdata)[dds_genes]

bulkdata_rescaled[,dds_gene_name] = normcounts[,dds_gene_name]
melted_resc = melt(bulkdata_rescaled,ids=c(ids,sample_group,study_group,modality_group,disease_group))

# log all the couunts
melted_resc$value = log(melted_resc$value+1)

bp_resc = ggplot(melted_resc, aes(x=ids, y=value, fill=study_group)) +
  geom_boxplot(aes(fill=study_group)) +
  facet_nested(. ~ modality_group + disease_group, drop=TRUE, scales = "free", space = "free") +
  theme_bw() +
  theme(
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank()) +
  labs(title="Sequential SC and PBulk Scaling -> Log Pseudobulk and Bulk Expression", x = "bulk/pseudobulk donor", y= "log expression")
bp_resc
ggsave(plot=bp_resc,file=file.path(plotdir,"logdata_resc.pdf"))

##########################################################################
# plot real-valued transformation of the data
##########################################################################
vsd = vst(dds, blind=FALSE)
head(assay(vsd), 3)
select = order(rowMeans(counts(dds,normalized=TRUE)),
              decreasing=TRUE)[1:100]
df = as.data.frame(colData(dds)[,c("group","modality")])
pdf(file.path(plotdir,"recluster_allgenes_hm.pdf"),width = 10,height = 5)
p=pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
        cluster_cols=FALSE, annotation_col=df)
dev.off()

# pca scatterplot
pdf(file.path(plotdir,"recluster_allgenes_scatter.pdf"),width = 10,height = 5)
plotPCA(vsd, intgroup=c("group", "modality"))
dev.off()





# get the common genes
seur_no_H = seur_list[which(status != "H")] # Correa only has NL so we cant use H
status_no_H = status[which(status != "H")]
genes_no_H = rownames(seur_no_H[[1]])
seur_obj = list()
for(i in 1:length(seur_no_H)){
  seur_obj[[i]] = as.Seurat(seur_no_H[[i]])
  genes_no_H = intersect(genes_no_H,rownames(seur_no_H[[i]]))
}
# filter on common genes
for(i in 1:length(seur_no_H)){
  seur_no_H[[i]] = seur_no_H[[i]][genes_no_H,]
}

# independently normalize and find variable features
for(i in 1:length(seur_obj)){
  seur_obj[[i]] = NormalizeData(seur_obj[[i]])
  seur_obj[[i]] = FindVariableFeatures(seur_obj[[i]], selection.method = "vst", nfeatures = 2000)
}

# add status and donor index
for(i in 1:length(seur_obj)){
  seur_obj[[i]]@meta.data$status = rep(status_no_H[[i]],dim(seur_obj[[i]]@meta.data)[1])
}

# select features that are repeatedly variable across datasets for integration
features = SelectIntegrationFeatures(object.list = seur_obj)

skin.anchors = FindIntegrationAnchors(object.list = seur_obj, anchor.features = features)

skin.combined = IntegrateData(anchorset = skin.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(skin.combined) = "integrated"



#############################################################################################################
#
# create seurat object using bulk/pseudo DEG
#
#############################################################################################################

skin.combined <- SetIdent(skin.combined, value = skin.combined@meta.data$status)

seur_no_H_degs = seur_list[which(status != "H")] # Correa only has NL so we cant use H
status_no_H_degs = status[which(status != "H")]
genes_no_H_degs = intersect(rownames(seur_no_H_degs[[1]]),rownames(degs))
seur_obj_degs = list()
for(i in 1:length(seur_no_H_degs)){
  seur_obj_degs[[i]] = as.Seurat(seur_no_H_degs[[i]])
  genes_no_H_degs = intersect(genes_no_H_degs,rownames(seur_no_H_degs[[i]]))
}
# filter on common genes
for(i in 1:length(seur_no_H_degs)){
  seur_no_H_degs[[i]] = seur_no_H_degs[[i]][genes_no_H_degs,]
}

# independently normalize and find variable features
for(i in 1:length(seur_obj_degs)){
  seur_obj_degs[[i]] = NormalizeData(seur_obj_degs[[i]])
  seur_obj_degs[[i]] = FindVariableFeatures(seur_obj_degs[[i]], selection.method = "vst", nfeatures = 2000)
}

# add status and donor index
for(i in 1:length(seur_obj_degs)){
  seur_obj_degs[[i]]@meta.data$status = rep(status_no_H_degs[[i]],dim(seur_obj_degs[[i]]@meta.data)[1])
}

# select features that are repeatedly variable across datasets for integration
features.degs = SelectIntegrationFeatures(object.list = seur_obj_degs)

skin.anchors.degs = FindIntegrationAnchors(object.list = seur_obj_degs, anchor.features = features.degs)

skin.combined.degs.orig = IntegrateData(anchorset = skin.anchors.degs)


##########################################################################
# save/load the output
##########################################################################

saveRDS(bulkdata_orig,file.path(backupdir,"bulkdata_orig.RDS"))
saveRDS(status,file.path(backupdir,"status.RDS"))
saveRDS(genes,file.path(backupdir,"genes.RDS"))
saveRDS(skin.combined,file.path(backupdir,"skin.combined.noH.RDS"))
saveRDS(degs,file.path(backupdir,"degs.RDS"))
saveRDS(dds,file.path(backupdir,"dds.RDS"))
saveRDS(seur_list,file.path(backupdir,"seur_list.RDS"))
saveRDS(skin.combined.degs.orig,file.path(backupdir,"skin.combined.degs.orig.RDS"))


bulkdata_orig = readRDS(file.path(backupdir,"bulkdata_orig.RDS"))
skin.combined = readRDS(file.path(backupdir,"skin.combined.noH.RDS"))
degs = readRDS(file.path(backupdir,"degs.RDS"))
dds = readRDS(file.path(backupdir,"dds.RDS"))
seur_list = readRDS(file.path(backupdir,"seur_list.RDS"))
status = readRDS(file.path(backupdir,"status.RDS"))
genes = readRDS(file.path(backupdir,"genes.RDS"))
skin.combined.degs.orig = readRDS(file.path(backupdir,"skin.combined.degs.orig.RDS"))
rts = readRDS(file.path(backupdir,"skin.combined.degs.orig.RDS"))


#############################################################################################################
# identify the clusters present in the bulk/pseudo DEG clustering using the markers from He et al
#############################################################################################################
# we use table E2 from the supplement of He et al
skin.combined.degs = skin.combined.degs.orig
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(skin.combined.degs) = "integrated"

library("readxl")
library(httr)
tf=tempfile()
GET("https://ars.els-cdn.com/content/image/1-s2.0-S0091674920301822-mmc9.xlsx", write_disk(tf <- tempfile(fileext = ".xlsx")))
he_markers = read_excel(tf)

# use the scSorter package (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02281-7)
topgenes <- head(VariableFeatures(skin.combined.degs), 2000)
expr = GetAssayData(skin.combined.degs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

anno = he_markers[,c("cluster","gene","p_val_adj")]
colnames(anno) = c("Type","Marker","Weight")
# anno$Weight = anno$Weight + .Machine$double.xmin # make sure we can invert this
anno$Weight = 1 # just set this to a constant
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

library(scSorter)
rts <- scSorter(as.matrix(expr), as.data.frame(anno))

# add rts labels to seurat object
skin.combined.degs@meta.data$pred_type = rts$Pred_Type

# save the labels
saveRDS(rts,file.path(backupdir,"rts.RDS"))
saveRDS(skin.combined.degs,file.path(backupdir,"skin.combined.degs.RDS"))


#############################################################################################################
#
# run clustering on bulk/pseudo DEG, seurat sensitivity 0.01
#
#############################################################################################################
myres = .01
plotdir = file.path(basedir,paste0("plots-sense",myres))
dir.create(plotdir)

# bar plots of cell abundance:
# > table(skin.combined.degs@meta.data$pred_type[which(skin.combined.degs@meta.data$seurat_clusters==4)])/length(which(skin.combined.degs@meta.data$seurat_clusters==4))

# Run the standard workflow for visualization and clustering
skin.combined.degs = ScaleData(skin.combined.degs, verbose = FALSE)
skin.combined.degs = RunPCA(skin.combined.degs, npcs = 30, verbose = FALSE)
skin.combined.degs = RunUMAP(skin.combined.degs, reduction = "pca", dims = 1:30)
skin.combined.degs = FindNeighbors(skin.combined.degs, reduction = "pca", dims = 1:30)
skin.combined.degs = FindClusters(skin.combined.degs, resolution = myres)
skin.combined.degs.markers = FindAllMarkers(skin.combined.degs)

clusters = skin.combined.degs@meta.data$seurat_clusters

levels(clusters)[match("0",levels(clusters))] <- "VEC"
levels(clusters)[match("1",levels(clusters))] <- "FB"
levels(clusters)[match("2",levels(clusters))] <- "TC"
levels(clusters)[match("3",levels(clusters))] <- "KC-MEL_NC_SC"
levels(clusters)[match("4",levels(clusters))] <- "MAC-DC"
levels(clusters)[match("5",levels(clusters))] <- "PC-vSMC"

skin.combined.degs@meta.data$mapped_clusters = clusters

# Visualization label status
p1.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE,group.by = "status")
p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, repel = TRUE,group.by="mapped_clusters")

pdf(file.path(plotdir,"recluster_disease_degs.pdf"),width = 10,height = 5)
p1.degs + p2.degs
dev.off()

# Visualization label donors
skin.combined.degs@meta.data$ident_status = paste0(skin.combined.degs@meta.data$orig.ident,"_",skin.combined.degs@meta.data$status)
p1.degs = DimPlot(skin.combined.degs, reduction = "umap", group.by = "ident_status")
p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, repel = TRUE,group.by="mapped_clusters")

pdf(file.path(plotdir,"recluster_donor_degs.pdf"),width = 10,height = 5)
p1.degs + p2.degs
dev.off()

# Visualization predicted cell type
p1.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, group.by = "pred_type")
p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, repel = TRUE,group.by="mapped_clusters")

pdf(file.path(plotdir,"recluster_celltype_degs.pdf"),width = 10,height = 5)
p1.degs + p2.degs
dev.off()

pdf(file.path(plotdir,"recluster_fb1_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("fib1"=which(skin.combined.degs@meta.data$pred_type == "FB.1")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_fb2_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("fib2"=which(skin.combined.degs@meta.data$pred_type == "FB.2")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_fb3_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("fib3"=which(skin.combined.degs@meta.data$pred_type == "FB.3")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_tc_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("tc"=which(skin.combined.degs@meta.data$pred_type == "TC")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_macdc_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("mac-dc"=which(skin.combined.degs@meta.data$pred_type == "MAC-DC")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_kc1_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("kc1"=which(skin.combined.degs@meta.data$pred_type == "KC.1")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_kc2_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("kc2"=which(skin.combined.degs@meta.data$pred_type == "KC.2")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_kc3_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("kc3"=which(skin.combined.degs@meta.data$pred_type == "KC.3")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_PCvSMC_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("PCvSMC"=which(skin.combined.degs@meta.data$pred_type == "PC-vSMC")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"inter_cluster_types.pdf"),width = 10,height = 5)
    ggplot(data=skin.combined.degs@meta.data, aes(seurat_clusters))+
    geom_bar(aes(fill=as.factor(pred_type)), position="fill")
dev.off()

pdf(file.path(plotdir,"inter_donor_types.pdf"),width = 10,height = 5)
    ggplot(data=skin.combined.degs@meta.data, aes(orig.ident))+
    geom_bar(aes(fill=as.factor(pred_type)), position="fill")
dev.off()

pdf(file.path(plotdir,"inter_disease_types.pdf"),width = 10,height = 5)
    ggplot(data=skin.combined.degs@meta.data, aes(status))+
    geom_bar(aes(fill=as.factor(pred_type)), position="fill")
dev.off()


# save the object
saveRDS(skin.combined.degs,file.path(backupdir,paste0("skin.combined.degs",myres,".RDS")))
#############################################################################################################
# output the single-cell expression data/metadata for DURIAN imputation
#############################################################################################################

# get the NL data
seur_nl = skin.combined.degs[,skin.combined.degs@meta.data$status == "NL"]
nl_pDataC = data.frame(colnames(seur_nl),seur_nl@meta.data$mapped_clusters,seur_nl@meta.data$orig.ident)
colnames(nl_pDataC) = c("cellID","cellType","sampleID")
rownames(nl_pDataC) = nl_pDataC$cellID

comgenes = intersect(rownames(seur_nl@assays$originalexp@data),rownames(nlbulk))
inds = subsetsc(x=as.matrix(seur_nl@assays$originalexp@data),geneids=comgenes,nsd=3)

write.csv(seur_nl@assays$originalexp@data[inds$gene,inds$cell],file.path(durian_data_dir,paste0("HeNL.sense",myres,"_C.csv")))
write.csv(nl_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeNL.sense",myres,"_pDataC.csv")))
write.csv(nlbulk[inds$gene,],file.path(durian_data_dir,paste0("SuarezNL.sense",myres,"_T.csv")))

write.csv(edgeR::cpm(seur_nl@assays$originalexp@data[inds$gene,inds$cell]),file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm_C.csv")))
write.csv(nl_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm_pDataC.csv")))
write.csv(edgeR::cpm(nlbulk[inds$gene,]),file.path(durian_data_dir,paste0("SuarezNL.sense",myres,".cpm_T.csv")))

for(durian_max_cells in durian_max_cell_opts){
  set.seed(42)
  maxcells = min(length(inds$cell),durian_max_cells)
  cells = sample(1:length(inds$cell), maxcells, replace=F)
  write.csv(edgeR::cpm(seur_nl@assays$originalexp@data[inds$gene,inds$cell][,cells]),file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm.sub",maxcells,"_C.csv")))
  write.csv(nl_pDataC[inds$cell,][cells,],file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm.sub",maxcells,"_pDataC.csv")))
}
# get the LS data
seur_ls = skin.combined.degs[,skin.combined.degs@meta.data$status == "LS"]
ls_pDataC = data.frame(colnames(seur_ls),seur_ls@meta.data$mapped_clusters,seur_ls@meta.data$orig.ident)
colnames(ls_pDataC) = c("cellID","cellType","sampleID")
rownames(ls_pDataC) = ls_pDataC$cellID

comgenes = intersect(rownames(seur_ls@assays$originalexp@data),rownames(lsbulk))
inds = subsetsc(x=as.matrix(seur_ls@assays$originalexp@data),geneids=comgenes,nsd=3)

write.csv(seur_ls@assays$originalexp@data[inds$gene,inds$cell],file.path(durian_data_dir,paste0("HeLS.sense",myres,"_C.csv")))
write.csv(ls_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeLS.sense",myres,"_pDataC.csv")))
write.csv(lsbulk[inds$gene,],file.path(durian_data_dir,paste0("SuarezLS.sense",myres,"_T.csv")))

write.csv(edgeR::cpm(seur_ls@assays$originalexp@data[inds$gene,inds$cell]),file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm_C.csv")))
write.csv(ls_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm_pDataC.csv")))
write.csv(edgeR::cpm(lsbulk[inds$gene,]),file.path(durian_data_dir,paste0("SuarezLS.sense",myres,".cpm_T.csv")))

for(durian_max_cells in durian_max_cell_opts){
  set.seed(42)
  maxcells = min(length(inds$cell),durian_max_cells)
  cells = sample(1:length(inds$cell), maxcells, replace=F)
  write.csv(edgeR::cpm(seur_ls@assays$originalexp@data[inds$gene,inds$cell][,cells]),file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm.sub",maxcells,"_C.csv")))
  write.csv(ls_pDataC[inds$cell,][cells,],file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm.sub",maxcells,"_pDataC.csv")))
}

#############################################################################################################
#
# run clustering on bulk/pseudo DEG, seurat sensitivity 0.1
#
#############################################################################################################
myres = .1
plotdir = file.path(basedir,paste0("plots-sense",myres))
dir.create(plotdir)

# bar plots of cell abundance:
# > table(skin.combined.degs@meta.data$pred_type[which(skin.combined.degs@meta.data$seurat_clusters==4)])/length(which(skin.combined.degs@meta.data$seurat_clusters==4))

# Run the standard workflow for visualization and clustering
skin.combined.degs = ScaleData(skin.combined.degs, verbose = FALSE)
skin.combined.degs = RunPCA(skin.combined.degs, npcs = 30, verbose = FALSE)
skin.combined.degs = RunUMAP(skin.combined.degs, reduction = "pca", dims = 1:30)
skin.combined.degs = FindNeighbors(skin.combined.degs, reduction = "pca", dims = 1:30)
skin.combined.degs = FindClusters(skin.combined.degs, resolution = myres)
skin.combined.degs.markers = FindAllMarkers(skin.combined.degs)


clusters = skin.combined.degs@meta.data$seurat_clusters

levels(clusters)[match("0",levels(clusters))] <- "FB"
levels(clusters)[match("1",levels(clusters))] <- "VEC"
levels(clusters)[match("2",levels(clusters))] <- "TC"
levels(clusters)[match("3",levels(clusters))] <- "MAC-DC"
levels(clusters)[match("4",levels(clusters))] <- "KC-SGC"
levels(clusters)[match("5",levels(clusters))] <- "PC-vSMC"
levels(clusters)[match("6",levels(clusters))] <- "LEC"
levels(clusters)[match("7",levels(clusters))] <- "MEL_NC_SC"

skin.combined.degs@meta.data$mapped_clusters = clusters

# Visualization label status
p1.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE,group.by = "status")
p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, repel = TRUE,group.by="mapped_clusters")

pdf(file.path(plotdir,"recluster_disease_degs.pdf"),width = 10,height = 5)
p1.degs + p2.degs
dev.off()

# Visualization label donors
skin.combined.degs@meta.data$ident_status = paste0(skin.combined.degs@meta.data$orig.ident,"_",skin.combined.degs@meta.data$status)
p1.degs = DimPlot(skin.combined.degs, reduction = "umap", group.by = "ident_status")
p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, repel = TRUE,group.by="mapped_clusters")

pdf(file.path(plotdir,"recluster_donor_degs.pdf"),width = 10,height = 5)
p1.degs + p2.degs
dev.off()

# Visualization predicted cell type
p1.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, group.by = "pred_type")
# p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE)
p2.degs = DimPlot(skin.combined.degs, reduction = "umap", label = TRUE, repel = TRUE,group.by="mapped_clusters")

pdf(file.path(plotdir,"recluster_celltype_degs.pdf"),width = 10,height = 5)
p1.degs + p2.degs
dev.off()

pdf(file.path(plotdir,"recluster_fb1_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("fib1"=which(skin.combined.degs@meta.data$pred_type == "FB.1")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_fb2_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("fib2"=which(skin.combined.degs@meta.data$pred_type == "FB.2")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_fb3_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("fib3"=which(skin.combined.degs@meta.data$pred_type == "FB.3")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_tc_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("tc"=which(skin.combined.degs@meta.data$pred_type == "TC")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_macdc_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("mac-dc"=which(skin.combined.degs@meta.data$pred_type == "MAC-DC")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_kc1_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("kc1"=which(skin.combined.degs@meta.data$pred_type == "KC.1")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_kc2_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("kc2"=which(skin.combined.degs@meta.data$pred_type == "KC.2")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_kc3_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("kc3"=which(skin.combined.degs@meta.data$pred_type == "KC.3")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"recluster_PCvSMC_degs.pdf"),width = 10,height = 5)
DimPlot(skin.combined.degs, reduction = "umap", group.by = "pred_type",
        cells.highlight = list("PCvSMC"=which(skin.combined.degs@meta.data$pred_type == "PC-vSMC")),
        split.by = "status")
dev.off()

pdf(file.path(plotdir,"inter_cluster_types.pdf"),width = 10,height = 5)
    ggplot(data=skin.combined.degs@meta.data, aes(seurat_clusters))+
    geom_bar(aes(fill=as.factor(pred_type)), position="fill")
dev.off()

pdf(file.path(plotdir,"inter_donor_types.pdf"),width = 10,height = 5)
    ggplot(data=skin.combined.degs@meta.data, aes(orig.ident))+
    geom_bar(aes(fill=as.factor(pred_type)), position="fill")
dev.off()

pdf(file.path(plotdir,"inter_disease_types.pdf"),width = 10,height = 5)
    ggplot(data=skin.combined.degs@meta.data, aes(status))+
    geom_bar(aes(fill=as.factor(pred_type)), position="fill")
dev.off()


# save the object
saveRDS(skin.combined.degs,file.path(backupdir,paste0("skin.combined.degs",myres,".RDS")))
#############################################################################################################
# output the single-cell expression data/metadata for DURIAN imputation
#############################################################################################################

# get the NL data
seur_nl = skin.combined.degs[,skin.combined.degs@meta.data$status == "NL"]
nl_pDataC = data.frame(colnames(seur_nl),seur_nl@meta.data$mapped_clusters,seur_nl@meta.data$orig.ident)
colnames(nl_pDataC) = c("cellID","cellType","sampleID")
rownames(nl_pDataC) = nl_pDataC$cellID

comgenes = intersect(rownames(seur_nl@assays$originalexp@data),rownames(nlbulk))
inds = subsetsc(x=as.matrix(seur_nl@assays$originalexp@data),geneids=comgenes,nsd=3)

write.csv(seur_nl@assays$originalexp@data[inds$gene,inds$cell],file.path(durian_data_dir,paste0("HeNL.sense",myres,"_C.csv")))
write.csv(nl_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeNL.sense",myres,"_pDataC.csv")))
write.csv(nlbulk[inds$gene,],file.path(durian_data_dir,paste0("SuarezNL.sense",myres,"_T.csv")))

write.csv(edgeR::cpm(seur_nl@assays$originalexp@data[inds$gene,inds$cell]),file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm_C.csv")))
write.csv(nl_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm_pDataC.csv")))
write.csv(edgeR::cpm(nlbulk[inds$gene,]),file.path(durian_data_dir,paste0("SuarezNL.sense",myres,".cpm_T.csv")))

for(durian_max_cells in durian_max_cell_opts){
  set.seed(42)
  maxcells = min(length(inds$cell),durian_max_cells)
  cells = sample(1:length(inds$cell), maxcells, replace=F)
  write.csv(edgeR::cpm(seur_nl@assays$originalexp@data[inds$gene,inds$cell][,cells]),file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm.sub",maxcells,"_C.csv")))
  write.csv(nl_pDataC[inds$cell,][cells,],file.path(durian_data_dir,paste0("HeNL.sense",myres,".cpm.sub",maxcells,"_pDataC.csv")))
}
# get the LS data
seur_ls = skin.combined.degs[,skin.combined.degs@meta.data$status == "LS"]
ls_pDataC = data.frame(colnames(seur_ls),seur_ls@meta.data$mapped_clusters,seur_ls@meta.data$orig.ident)
colnames(ls_pDataC) = c("cellID","cellType","sampleID")
rownames(ls_pDataC) = ls_pDataC$cellID

comgenes = intersect(rownames(seur_ls@assays$originalexp@data),rownames(lsbulk))
inds = subsetsc(x=as.matrix(seur_ls@assays$originalexp@data),geneids=comgenes,nsd=3)

write.csv(seur_ls@assays$originalexp@data[inds$gene,inds$cell],file.path(durian_data_dir,paste0("HeLS.sense",myres,"_C.csv")))
write.csv(ls_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeLS.sense",myres,"_pDataC.csv")))
write.csv(lsbulk[inds$gene,],file.path(durian_data_dir,paste0("SuarezLS.sense",myres,"_T.csv")))

write.csv(edgeR::cpm(seur_ls@assays$originalexp@data[inds$gene,inds$cell]),file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm_C.csv")))
write.csv(ls_pDataC[inds$cell,],file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm_pDataC.csv")))
write.csv(edgeR::cpm(lsbulk[inds$gene,]),file.path(durian_data_dir,paste0("SuarezLS.sense",myres,".cpm_T.csv")))

for(durian_max_cells in durian_max_cell_opts){
  set.seed(42)
  maxcells = min(length(inds$cell),durian_max_cells)
  cells = sample(1:length(inds$cell), maxcells, replace=F)
  write.csv(edgeR::cpm(seur_ls@assays$originalexp@data[inds$gene,inds$cell][,cells]),file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm.sub",maxcells,"_C.csv")))
  write.csv(ls_pDataC[inds$cell,][cells,],file.path(durian_data_dir,paste0("HeLS.sense",myres,".cpm.sub",maxcells,"_pDataC.csv")))
}
