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
library(umap)
library(ggpubr)



outputmaster = Sys.getenv("RUNMASTER")

### external code
source("slurm/scrabble_helper_functions/library_extra_plots.R")
source("slurm/scrabble_helper_functions/library_cluster_metrics.R")


paramsets = list.files(outputmaster,full.names=TRUE,include.dirs=TRUE)
for(paramind in 1:length(paramsets)){
    dropsets = list.files(paramsets[paramind],full.names=TRUE,include.dirs=TRUE)
    dropsets_short = list.files(paramsets[paramind],include.dirs=TRUE)
    
    fit_inds = grep("output_fit",dropsets)

    for(fit_ind in fit_inds){
        dataparams = get_model_params(dropsets_short[fit_ind])
        simdirs = list.files(dropsets[fit_ind],full.names=TRUE)
        for(simdir in simdirs[1:2]){
            modeldirs = list.files(simdir,full.names=TRUE)
            modeldirs = modeldirs[grep("imputemodel_",modeldirs)]
            imethods = sapply(modeldirs,strsplit,"/") %>% sapply(last) %>% unname
            ind_dropmodel = grep("imputemodel_dropout",modeldirs)
            pDataC = read.csv(file.path(simdir,"pDataC.csv"),row.names=1)
            trueC = read.csv(file.path(simdir,"trueC.csv"),row.names=1)
            if(file.exists(file.path(modeldirs[ind_dropmodel],"imputation_loss.csv"))){
                droprate_backup = read.csv(file.path(modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$Dropout
                for(methodind in 1:length(imethods)){
                    myparams = get_model_params(imethods[methodind])
                    imputedir = modeldirs[methodind]
                    imethod = myparams[["imputemodel"]]
                    ldfname = file.path(imputedir,paste0(imethod,"_logdf.csv"))
                    scldfname = file.path(imputedir,"outerStats",paste0(imethod,"_logdf.csv")) # if we are running scrabble and things did not finish
                    impname = file.path(imputedir,"imputed_C.csv")
                    if(file.exists(impname)){
                        dir.create(file.path(imputedir,"cluster_plots"))
                        imputedC = read.csv(impname,row.names=1)
                        mergeC = trueC
                        mergeC[rownames(imputedC),colnames(imputedC)] = imputedC
                        p_lr=logratio_plot(imputedC=mergeC,pDataC=pDataC,trueC=trueC,plottitle="MA Plot")
                        p_umap=umap_plot(dat=mergeC,meta=pDataC,plottitle="UMAP Plot",ptsize=2)
                        p_tsne=tsne_plot(dat=mergeC,meta=pDataC,plottitle="tSNE Plot",col="sampleID",ptsize=2)

                        p_lr_control=logratio_plot(imputedC=trueC,pDataC=pDataC,trueC=trueC,plottitle="MA Plot (True)")
                        p_umap_control=umap_plot(dat=trueC,meta=pDataC,plottitle="UMAP Plot (True)",ptsize=2)
                        p_tsne_control=tsne_plot(dat=trueC,meta=pDataC,plottitle="tSNE Plot (True)",col="sampleID",ptsize=2)

                        ggsave(plot=p_lr,filename=file.path(imputedir,"cluster_plots","ma.pdf"),width=5,height=5)
                        ggsave(plot=p_umap,filename=file.path(imputedir,"cluster_plots","umap.pdf"),width=5,height=5)
                        ggsave(plot=p_tsne,filename=file.path(imputedir,"cluster_plots","tsne.pdf"),width=5,height=5)
                        
                        ggsave(plot=p_lr_control,filename=file.path(imputedir,"cluster_plots","ma_control.pdf"),width=5,height=5)
                        ggsave(plot=p_umap_control,filename=file.path(imputedir,"cluster_plots","umap_control.pdf"),width=5,height=5)
                        ggsave(plot=p_tsne_control,filename=file.path(imputedir,"cluster_plots","tsne_control.pdf"),width=5,height=5)
                    }
                }
            }

        }
    }
}

