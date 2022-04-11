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

source("slurm/scrabble_helper_functions/library_cluster_metrics.R")
source("slurm/scrabble_helper_functions/library_extra_plots.R")
source("slurm/scrabble_helper_functions/library_scrabble_clusterMetrics_clValid.R")

metafiles = c(
  "slurm/Baron/durian_data/BaronSC.H_pDataC.csv",
  "slurm/Baron/durian_data/BaronSC.DM_pDataC.csv",
  "slurm/He/durian_data/HeNL_pDataC.csv",
  "slurm/He/durian_data/HeLS_pDataC.csv",
  "slurm/Gupta/durian_data/GuptaE13SC_pDataC.csv")

allmetadata = NULL
for(mf in metafiles){
  meta = read.csv(mf,row.names=1)
  row.names(meta) = meta$cellID
  allmetadata = rbind(allmetadata,meta)
}

merged_list = list()
pairs_list = list()
celltypes_list = list()
pathways_list = list()

#### Baron durian vs drimpute ####
outputdir = "slurm/fig_panel_scripts/fig05"
dir.create(outputdir,recursive=TRUE)
datadir = "slurm/Baron/output.clusterMetrics.free.OuterMetrics1K05/pref_BaronSC.DM.isletVST1K05;SegerstolpeBulk.DM.cpm,sub_,suff_OuterMetrics1K05/output_fit"
CellChatDB = CellChatDB.human
# consider drimpute and durian with alpha=1 
modeldirs = list.files(datadir,full.names = TRUE)[regexpr("((DURIAN).*(sA_1))|(DrImpute)",list.files(datadir)) > 0]
cellchats_list = list()
for(i in 1:length(modeldirs)){
    params = get_model_params(tail(strsplit(modeldirs[i],"/")[[1]],n=1))
    imputedC = read.csv(file.path(modeldirs[i],"imputed_C.csv"),row.names=1)
    pDataC = allmetadata[colnames(imputedC),]

    # meta$sampleID = paste(meta$cellID,params$imputemodel,sep=".")
 
    # rownames(meta) = meta$cellID
    # rownames(imputedC) = meta$cellID
    #############################################################################
    # select disease subsets and normalize
    #############################################################################
    data = normalizeData(as.matrix(imputedC))
    meta = pDataC[,c("cellType","sampleID")]

    print("creating cellchat object")
    cellchat <- createCellChat(object = data, meta = meta, group.by = "cellType")

    #############################################################################
    # set ligrec database
    #############################################################################

    cellchat@DB <- CellChatDB

    #############################################################################
    # preprocess the data
    #############################################################################
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    #############################################################################
    # compute probs without preprocessing
    #############################################################################
    cellchat <- computeCommunProb(cellchat,raw.use = TRUE)


    #############################################################################
    # compute pathway probs
    #############################################################################
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    cellchats_list[[params$imputemodel]] = cellchat
}
merged = mergeCellChat(cellchats_list, add.names = names(cellchats_list),cell.prefix=TRUE)
merged_list[["baron"]] = merged
pairs_list[["baron"]] = cellchats_list
celltypes_list[["baron"]] = levels(pairs_list[["baron"]][[2]]@idents)
pathways_list[["baron"]] = pairs_list[["baron"]][[2]]@netP$pathways

#############################################################################
# differential circle plots
#############################################################################
pdf(file.path(outputdir,"baron_circle.pdf"))
  par(mfrow = c(1,2), xpd=TRUE)
  p = netVisual_diffInteraction(merged_list[["baron"]], weight.scale = T)
  netVisual_diffInteraction(merged_list[["baron"]], weight.scale = T, measure = "weight")
dev.off()

#############################################################################
# differential heatmaps
#############################################################################
pdf(file.path(outputdir,"baron_heat.pdf"))
  gg1 <- netVisual_heatmap(merged_list[["baron"]])
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(merged_list[["baron"]], measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
dev.off()

#############################################################################
# differential bubble plots
#############################################################################

pdf(file.path(outputdir,"baron_bubble.pdf"),width=15.5,height=10.0)
  netVisual_bubble(merged_list[["baron"]], remove.isolate = FALSE,comparison=c(1,2),angle.x=45)
dev.off()

pdf(file.path(outputdir,"baron_bubble_wnt.pdf"),width=7.5,height=3.0)
  netVisual_bubble(merged_list[["baron"]], signaling = c("WNT"), remove.isolate = FALSE,comparison=c(1,2),angle.x=45)
dev.off()

#############################################################################
# differential violin plots
#############################################################################
pdf(file.path(outputdir,"baron_violin_wnt.pdf"))
  plotGeneExpression(merged_list[["baron"]], signaling = "WNT", split.by = "datasets", colors.ggplot = T)
dev.off()

#### He durian vs drimpute ####
outputdir = "slurm/fig_panel_scripts/fig05"
dir.create(outputdir,recursive=TRUE)
datadir = "slurm/He/output.clusterMetrics.free.OuterMetrics1K05/pref_HeSC.LS.cpm1K05;SuarezLS.cpm,sub_,suff_OuterMetrics1K05/output_fit"
CellChatDB = CellChatDB.human
# consider drimpute and durian with alpha=1 
modeldirs = list.files(datadir,full.names = TRUE)[regexpr("((DURIAN).*(sA_1))|(DrImpute)",list.files(datadir)) > 0]
cellchats_list = list()
for(i in 1:length(modeldirs)){
    params = get_model_params(tail(strsplit(modeldirs[i],"/")[[1]],n=1))
    imputedC = read.csv(file.path(modeldirs[i],"imputed_C.csv"),row.names=1)
    pDataC = allmetadata[colnames(imputedC),]

    # meta$sampleID = paste(meta$cellID,params$imputemodel,sep=".")
 
    # rownames(meta) = meta$cellID
    # rownames(imputedC) = meta$cellID
    #############################################################################
    # select disease subsets and normalize
    #############################################################################
    data = normalizeData(as.matrix(imputedC))
    meta = pDataC[,c("cellType","sampleID")]

    print("creating cellchat object")
    cellchat <- createCellChat(object = data, meta = meta, group.by = "cellType")

    #############################################################################
    # set ligrec database
    #############################################################################

    cellchat@DB <- CellChatDB

    #############################################################################
    # preprocess the data
    #############################################################################
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    #############################################################################
    # compute probs without preprocessing
    #############################################################################
    cellchat <- computeCommunProb(cellchat,raw.use = TRUE)


    #############################################################################
    # compute pathway probs
    #############################################################################
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    cellchats_list[[params$imputemodel]] = cellchat
}
merged = mergeCellChat(cellchats_list, add.names = names(cellchats_list),cell.prefix=TRUE)
merged_list[["he"]] = merged
pairs_list[["he"]] = cellchats_list
celltypes_list[["he"]] = levels(pairs_list[["he"]][[2]]@idents)
pathways_list[["he"]] = pairs_list[["he"]][[2]]@netP$pathways


#############################################################################
# differential circle plots
#############################################################################
pdf(file.path(outputdir,"he_circle.pdf"))
  par(mfrow = c(1,2), xpd=TRUE)
  p = netVisual_diffInteraction(merged_list[["he"]], weight.scale = T)
  netVisual_diffInteraction(merged_list[["he"]], weight.scale = T, measure = "weight")
dev.off()

#############################################################################
# differential heatmaps
#############################################################################
pdf(file.path(outputdir,"he_heat.pdf"))
  gg1 <- netVisual_heatmap(merged_list[["he"]])
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(merged_list[["he"]], measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
dev.off()

#############################################################################
# differential bubble plots
#############################################################################

pdf(file.path(outputdir,"he_bubble.pdf"),width=50.5,height=40.0)
  netVisual_bubble(merged_list[["he"]], remove.isolate = FALSE,comparison=c(1,2),angle.x=45)
dev.off()

pdf(file.path(outputdir,"he_bubble_sub.pdf"),width=6.5,height=18.0)
  netVisual_bubble(merged_list[["he"]], remove.isolate = FALSE,comparison=c(1,2),angle.x=45,sources.use=c("Inflam. DC","Inflam. FIB"),targets.use=c("Inflam. DC","Inflam. FIB"))
dev.off()



#############################################################################
# differential violin plots
#############################################################################
pdf(file.path(outputdir,"he_violin.pdf"))
  plotGeneExpression(merged_list[["he"]], split.by = "datasets", colors.ggplot = T)
dev.off()

#### gupta durian vs drimpute ####
outputdir = "slurm/fig_panel_scripts/fig05"
dir.create(outputdir,recursive=TRUE)
datadir = "slurm/Gupta/output.clusterMetrics.free.OuterMetrics1K05/pref_GuptaE13SC.VST1K05;BiggsBulk.VST,sub_,suff_OuterMetrics1K05/output_fit"
CellChatDB = CellChatDB.mouse
# consider drimpute and durian with alpha=1
modeldirs = list.files(datadir,full.names = TRUE)[regexpr("((DURIAN).*(sA_1))|(DrImpute)",list.files(datadir)) > 0]
cellchats_list = list()
for(i in 1:length(modeldirs)){
    params = get_model_params(tail(strsplit(modeldirs[i],"/")[[1]],n=1))
    imputedC = read.csv(file.path(modeldirs[i],"imputed_C.csv"),row.names=1)
    pDataC = allmetadata[colnames(imputedC),]

    # meta$sampleID = paste(meta$cellID,params$imputemodel,sep=".")
 
    # rownames(meta) = meta$cellID
    # rownames(imputedC) = meta$cellID
    #############################################################################
    # select disease subsets and normalize
    #############################################################################
    data = normalizeData(as.matrix(imputedC))
    meta = pDataC[,c("cellType","sampleID")]

    print("creating cellchat object")
    cellchat <- createCellChat(object = data, meta = meta, group.by = "cellType")

    #############################################################################
    # set ligrec database
    #############################################################################

    cellchat@DB <- CellChatDB

    #############################################################################
    # preprocess the data
    #############################################################################
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    #############################################################################
    # compute probs without preprocessing
    #############################################################################
    cellchat <- computeCommunProb(cellchat,raw.use = TRUE)


    #############################################################################
    # compute pathway probs
    #############################################################################
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    cellchats_list[[params$imputemodel]] = cellchat
}
merged = mergeCellChat(cellchats_list, add.names = names(cellchats_list),cell.prefix=TRUE)
merged_list[["gupta"]] = merged
pairs_list[["gupta"]] = cellchats_list
celltypes_list[["gupta"]] = levels(pairs_list[["gupta"]][[2]]@idents)
pathways_list[["gupta"]] = pairs_list[["gupta"]][[2]]@netP$pathways


#############################################################################
# differential circle plots
#############################################################################
pdf(file.path(outputdir,"gupta_circle.pdf"))
par(mfrow = c(1,2), xpd=TRUE)
p = netVisual_diffInteraction(merged_list[["gupta"]], weight.scale = T)
netVisual_diffInteraction(merged_list[["gupta"]], weight.scale = T, measure = "weight")
dev.off()

#############################################################################
# differential heatmaps
#############################################################################
pdf(file.path(outputdir,"gupta_heatmap.pdf"))
  gg1 <- netVisual_heatmap(merged_list[["gupta"]])
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(merged_list[["gupta"]], measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
dev.off()

#############################################################################
# differential bubble plots
#############################################################################

pdf(file.path(outputdir,"gupta_bubble_ncwnt_sub.pdf"),width=6.5,height=3.0)
  netVisual_bubble(merged_list[["gupta"]], signaling = c("ncWNT"),sources.use = c(4,5,6), targets.use = c(1,2,4,5,6), remove.isolate = FALSE,comparison=c(1,2),angle.x=45)
dev.off()

pdf(file.path(outputdir,"gupta_bubble_col.pdf"),width=6.5,height=6.5)
  netVisual_bubble(merged_list[["gupta"]], signaling = c("COLLAGEN"),sources.use = c(4,5,6), targets.use = c(1,2,4,5,6), remove.isolate = FALSE,comparison=c(1,2),angle.x=45)
dev.off()

#############################################################################
# differential violin plots
#############################################################################

palette1 = readRDS("slurm/fig_panel_scripts/palette1.RDS")
pdf(file.path(outputdir,"gupta_violin_ncwnt_sub.pdf"),width=6.5,height=5.5)
  plotGeneExpression(merged_list[["gupta"]], signaling = c("ncWNT"), split.by = "datasets", colors.ggplot = T)
dev.off()

pdf(file.path(outputdir,"gupta_violin_ncwnt_sub_thin.pdf"),width=3.5,height=5.5)
  plotGeneExpression(merged_list[["gupta"]], signaling = c("ncWNT"), split.by = "datasets", colors.ggplot = T)
dev.off()

#############################################################################
# differential information flow
#############################################################################
merged_list[["gupta"]] <- computeNetSimilarityPairwise(merged_list[["gupta"]], type = "structural")
#> Compute signaling network similarity for datasets 1 2
merged_list[["gupta"]] <- netEmbedding(merged_list[["gupta"]], type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
merged_list[["gupta"]] <- netClustering(merged_list[["gupta"]], type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
gg1 = rankNet(merged_list[["gupta"]], mode = "comparison", stacked = T, do.stat = TRUE)
ggsave(plot=gg1,file=file.path(outputdir,"gupta_inf_flow.pdf"),width=4.5,height=5)
