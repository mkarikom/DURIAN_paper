###############################################
# delete this
###############################################
# Sys.setenv(PBULKDIRCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/output.final.baron")
# Sys.setenv(SUMMARYCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/realdata_scripts/test/summaryfinalcombine")
# Sys.setenv(SPARSITY_PARAM="lambda")
###############################################

source("slurm/scrabble_helper_functions/library_cluster_metrics.R")
source("slurm/scrabble_helper_functions/library_extra_plots.R")

library(dplyr)
library(rmcorr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(scales)
library(umap)
library(viridis)
library(jcolors)

excludemodels = c()
quicksims = c(5)
quickdparams = c("0:0:0:0:6.5:6.5:6.5:6.5",1e-6)
umap_min_dist = 0.9
outputdir = file.path("slurm","fig_panel_scripts","fig02_scatter")
dir.create(outputdir,recursive=TRUE)

# saved benchmark directory
# backupdircombine = "slurm/durian_pseudobulk_OuterStats_clValid/output.final.splatter,n_5"
# outputdir = "slurm/durian_pseudobulk_OuterStats_clValid/output.final.splatter,n_5,output_summaryCombineParam"
# sparsityparam = "dmid"

backupdircombine = "slurm/durian_pseudobulk_BaronOuterStatsAllNested_clValidInternal/output.final.baron,n_5"
sparsityparam = "lambda"
sparsityvar = "λ"

dirnamescombine = list.files(backupdircombine,full.names=TRUE)
df.orig = NULL # the dataframe of all output, to be filled
idx = 1
for(ii in 1:length(dirnamescombine)){
  if(length(grep("_fit",dirnamescombine[ii]))>0){
    backupdir = dirnamescombine[ii]
    bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
    dirnames = list.files(backupdir)
    for(i in 1:length(dirnames)){
      params1 = get_model_params(str_split(dirnames[i],pattern=",")[[1]])
      params1b = get_model_params(bdname)
      params = bind_rows(c(params1b,params1))

      if(params$sim %in% quicksims & params[,sparsityparam] %in% quickdparams){
        modeldirs = list.dirs(file.path(backupdir,dirnames[i]),recursive=FALSE,full.names=FALSE)
        ind_dropmodel = grep("imputemodel_dropout",modeldirs)
        if(file.exists(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"))){
          droprate_backup = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$Dropout
          dropoutC = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputed_C.csv"),row.names=1)
          trueC = read.csv(file.path(backupdir,dirnames[i],"trueC.csv"),row.names=1)
          umap_res <- umap(t(trueC),min_dist=umap_min_dist)
          tsne <- Rtsne::Rtsne(t(trueC),check_duplicates=FALSE)

          for(j in 1:length(modeldirs)){
            if(length(grep("imputemodel",modeldirs[j]))>0){
              print(paste0("found impute folder ",modeldirs[j]))
              logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
              imputedfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputed_C.csv")
              if(file.exists(imputedfn)){
                imputedC = read.csv(imputedfn,row.names=1)

                diffC = imputedC-trueC

                normDiffC = apply(diffC, 2, function(x){sqrt(sum(x^2))})

                pDataC = read.csv(file.path(backupdir,dirnames[i],"pDataC.csv"),row.names=1)
                pDataC = cbind(pDataC,as.data.frame(normDiffC))

                print("creating umap")
                df_umap = cbind(umap_res$layout,pDataC)
                colnames(df_umap)=c("U1","U2",colnames(pDataC))
                df_umap$cellid = rownames(df_umap)
                print("umap created")

                df_tsne = cbind(tsne$Y,df_umap)
                colnames(df_tsne)=c("T1","T2",colnames(df_umap))


                df_umap = df_tsne

                print(paste0("found impute model ",logfn))
                params2 = as.data.frame(get_model_params(modeldirs[j]))
                dataparams = cbind(params,params2)
                runlog = read.csv(logfn,row.names=1)

                runlog = cbind(runlog,dataparams)
                runlog$dropout_rate = droprate_backup

                if(length(grep("DURIAN",runlog$modelname))>0){
                  runlog$modelfam=c("DURIAN")
                  runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
                }else if(length(grep("mtSCRABBLE",runlog$modelname))>0){
                  runlog$modelfam=c("Control")
                  runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
                }else if(length(grep("SCRABBLE",runlog$modelname))>0){
                  runlog$modelfam=c("Existing Methods")
                  runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
                }else if(length(grep("dropout",runlog$modelname))>0){
                  runlog$modelfam=c("Control")
                }else{
                  runlog$modelfam=c("Existing Methods")
                }

                print(paste0("updating df (",class(runlog) ,") with runlog: \n"))
                print(runlog)
                print("rownames:")
                print(rownames(runlog))

                stretchrl = NULL
                for(rowi in 1:nrow(df_umap)){
                  stretchrl = rbind(stretchrl,runlog)
                }
                # runlog_umap = runlog %>% dplyr::slice(rep(1:n(), each = nrow(df_umap)))
                runlog_umap = cbind(stretchrl,df_umap)

                print("calling rbind")
                df.orig = plyr::rbind.fill(df.orig,runlog_umap)
              }
            }
          }
        }
      }
    }
  }
}
df.orig$logRMSE = log(df.orig$RMSE)
df.orig$logENORM = log(df.orig$ENORM)


df.orig$dropoutlevel = paste0(sparsityvar,"=",df.orig[,sparsityparam])
df.filter = df.orig %>% filter(!modelname %in% excludemodels)
df_downsamp = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
df_downsamp$strategy = "Down-Sampling"

write.csv(df_downsamp,file.path(outputdir,"df_downsamp.csv"))


backupdircombine = "slurm/durian_pseudobulk_OuterStats_clValid/output.final.splatter,n_5"
sparsityparam = "dmid"
sparsityvar = "x0"

# backupdircombine = "slurm/durian_pseudobulk_BaronOuterStatsAllNested_clValidInternal/output.final.baron,n_5"
# outputdir = "slurm/durian_pseudobulk_BaronOuterStatsAllNested_clValidInternal/output.final.baron,n_5,output_summaryCombineParam"
# sparsityparam = "lambda"
# sparsityvar = "λ"

dirnamescombine = list.files(backupdircombine,full.names=TRUE)
df.orig = NULL # the dataframe of all output, to be filled
idx = 1
for(ii in 1:length(dirnamescombine)){
  if(length(grep("_fit",dirnamescombine[ii]))>0){
    backupdir = dirnamescombine[ii]
    bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
    dirnames = list.files(backupdir)
    for(i in 1:length(dirnames)){
      params1 = get_model_params(str_split(dirnames[i],pattern=",")[[1]])
      params1b = get_model_params(bdname)
      params = bind_rows(c(params1b,params1))

      if(params$sim %in% quicksims & params[,sparsityparam] %in% quickdparams){
        modeldirs = list.dirs(file.path(backupdir,dirnames[i]),recursive=FALSE,full.names=FALSE)
        ind_dropmodel = grep("imputemodel_dropout",modeldirs)
        if(file.exists(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"))){
          droprate_backup = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$Dropout
          dropoutC = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputed_C.csv"),row.names=1)
          trueC = read.csv(file.path(backupdir,dirnames[i],"trueC.csv"),row.names=1)
          umap_res <- umap(t(trueC),min_dist=umap_min_dist)
          tsne <- Rtsne::Rtsne(t(trueC),check_duplicates=FALSE)
          for(j in 1:length(modeldirs)){
            if(length(grep("imputemodel",modeldirs[j]))>0){
              print(paste0("found impute folder ",modeldirs[j]))
              logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
              imputedfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputed_C.csv")
              if(file.exists(imputedfn)){
                imputedC = read.csv(imputedfn,row.names=1)

                diffC = imputedC-trueC

                normDiffC = apply(diffC, 2, function(x){sqrt(sum(x^2))})

                pDataC = read.csv(file.path(backupdir,dirnames[i],"pDataC.csv"),row.names=1)
                pDataC = cbind(pDataC,as.data.frame(normDiffC))

                print("creating umap")
                df_umap = cbind(umap_res$layout,pDataC)
                colnames(df_umap)=c("U1","U2",colnames(pDataC))
                df_umap$cellid = rownames(df_umap)
                print("umap created")

                df_tsne = cbind(tsne$Y,df_umap)
                colnames(df_tsne)=c("T1","T2",colnames(df_umap))
                
                df_umap = df_tsne

                print(paste0("found impute model ",logfn))
                params2 = as.data.frame(get_model_params(modeldirs[j]))
                dataparams = cbind(params,params2)
                runlog = read.csv(logfn,row.names=1)

                runlog = cbind(runlog,dataparams)
                runlog$dropout_rate = droprate_backup

                if(length(grep("DURIAN",runlog$modelname))>0){
                  runlog$modelfam=c("DURIAN")
                  runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
                }else if(length(grep("mtSCRABBLE",runlog$modelname))>0){
                  runlog$modelfam=c("Control")
                  runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
                }else if(length(grep("SCRABBLE",runlog$modelname))>0){
                  runlog$modelfam=c("Existing Methods")
                  runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
                }else if(length(grep("dropout",runlog$modelname))>0){
                  runlog$modelfam=c("Control")
                }else{
                  runlog$modelfam=c("Existing Methods")
                }

                print(paste0("updating df (",class(runlog) ,") with runlog: \n"))
                print(runlog)
                print("rownames:")
                print(rownames(runlog))

                stretchrl = NULL
                for(rowi in 1:nrow(df_umap)){
                  stretchrl = rbind(stretchrl,runlog)
                }
                # runlog_umap = runlog %>% dplyr::slice(rep(1:n(), each = nrow(df_umap)))
                runlog_umap = cbind(stretchrl,df_umap)

                print("calling rbind")
                df.orig = plyr::rbind.fill(df.orig,runlog_umap)
              }
            }
          }
        }
      }
    }
  }
}
df.orig$logRMSE = log(df.orig$RMSE)
df.orig$logENORM = log(df.orig$ENORM)



df.orig$dropoutlevel = paste0(sparsityvar,"=",unlist(lapply(strsplit(df.orig[,sparsityparam],":"),tail,n=1)))
df.filter = df.orig %>% filter(!modelname %in% excludemodels)
df_sim = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
df_sim$strategy = "Simulation"

write.csv(df_downsamp,file.path(outputdir,"df_sim.csv"))


df = plyr::rbind.fill(df_sim,df_downsamp)

minalpha = min(1-1/log(df$normDiffC))
maxalpha = max(1-1/log(df$normDiffC))

df$normalpha = (1-1/log(df$normDiffC) - minalpha) / (maxalpha - minalpha)

df$mean_dropout = round(df$mean_dropout,digits=3)
nmodels = length(unique(df$modelname))
palette1 = scales::hue_pal()(nmodels)
names(palette1) = sort(unique(df$modelname))

write.csv(df,file.path(outputdir,"df.csv"))

df_tmp = df %>% filter(strategy=="Down-Sampling",mean_dropout=="0.914",sA %in% c(1,NA))
p = ggplot(df_tmp,aes(x=U1, y=U2,color=1-1/log(normDiffC),alpha=normalpha)) + 
scale_color_distiller(palette="PuBuGn",direction=1)+
      geom_point(size=0.25)+
      facet_nested(~modelfam + modelname)
p = p +
     theme_bw() +
     labs(alpha="Error Threshold",color="Cell-wise\nError") +
     rotate_x_text(90) +
     theme(
       axis.text=element_text(size=3,face="bold"),
       axis.title=element_text(size=5,face="bold"),
       # legend.position = "none",
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
       strip.text.x = element_text(size = 5, face = "bold"),
       strip.text.y = element_text(size = 5, face = "bold"),
        legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.spacing=unit(0,"lines"))
ggsave(plot=p,file=file.path(outputdir,"facet_scatter_alpha_ds.pdf"),width=7,height=2)

df_tmp = df %>% filter(strategy=="Simulation",mean_dropout=="0.928",sA %in% c(1,NA))
p = ggplot(df_tmp %>% slice_sample(n = 5000),aes(x=U1, y=U2,color=1-1/log(normDiffC),alpha=normalpha)) + 
scale_color_distiller(palette="PuBuGn",direction=1)+
      geom_point(size=0.25)+
      facet_nested(~modelfam + modelname)
p = p +
     theme_bw() +
     labs(alpha="Error Threshold",color="Cell-wise\nError") +
     rotate_x_text(90) +
     theme(
       axis.text=element_text(size=3,face="bold"),
       axis.title=element_text(size=5,face="bold"),
       # legend.position = "none",
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
       strip.text.x = element_text(size = 5, face = "bold"),
       strip.text.y = element_text(size = 5, face = "bold"),
        legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.spacing=unit(0,"lines"))
ggsave(plot=p,file=file.path(outputdir,"facet_scatter_alpha_sim.pdf"),width=7,height=2)


logistic_trans <- function(x,l=1,k=1,x0=0){
  l/(1+exp(-k*(x-x0)))
}

normalize_trans <- function(x){
  minx = min(x)
  maxx = max(x)
  (x-minx)/(maxx-minx)
}

df_tmp = df %>% filter(strategy=="Simulation",mean_dropout=="0.928",sA %in% c(1,NA))
df_tmp$normlogistic = normalize_trans(logistic_trans(log(df_tmp$normDiffC),k=1,x0=0,l=1))
# df_tmp$normlogistic = normalize_trans(logistic_trans(log(df_tmp$normDiffC),k=1,x0=10,l=1))
df_tmp$normexp = normalize_trans(logistic_trans(log(df_tmp$normDiffC),k=1,x0=10,l=1))
p = ggplot(df_tmp,aes(x=U1, y=U2,color=normlogistic,alpha=normexp)) + 
      # geom_point(size=0.25)+
      geom_point(aes(size=log(normDiffC)))+
      facet_nested(strategy~modelfam + imputemodel)
p = p +
     theme_bw() +
     rotate_x_text(90) +
     labs(alpha="Error Threshold",color="Cell-wise\nError") +
     scale_color_jcolors_contin(palette="pal12")+
     theme(
       axis.text=element_text(size=5,face="bold"),
       axis.title=element_text(size=10,face="bold"),
       # legend.position = "none",
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
       strip.text.x = element_text(size = 5, face = "bold"),
       strip.text.y = element_text(size = 5, face = "bold"),
        legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.spacing=unit(0,"lines"))+
         scale_size(range = c(0.01,1))+
         scale_alpha(range = c(0.01,.8))

ggsave(plot=p,file=file.path(outputdir,"facet_scatter_alpha_sim_normlogistic.pdf"),width=6.5,height=1.8)

df_tmp = df %>% filter(strategy=="Down-Sampling",mean_dropout=="0.914",sA %in% c(1,NA))
df_tmp$normlogistic = normalize_trans(logistic_trans(log(df_tmp$normDiffC),k=1,x0=0,l=1))
# df_tmp$normlogistic = normalize_trans(logistic_trans(log(df_tmp$normDiffC),k=1,x0=10,l=1))
df_tmp$normexp = normalize_trans(logistic_trans(log(df_tmp$normDiffC),k=1,x0=10,l=1))
p = ggplot(df_tmp,aes(x=U1, y=U2,color=normlogistic,alpha=normexp)) + 
      # geom_point(size=0.25)+
      geom_point(aes(size=log(normDiffC)))+
      facet_nested(strategy~modelfam + imputemodel)
p = p +
     theme_bw() +
     rotate_x_text(90) +
     labs(alpha="Error Threshold",color="Cell-wise\nError") +
     scale_color_jcolors_contin(palette="pal12")+
     theme(
       axis.text=element_text(size=5,face="bold"),
       axis.title=element_text(size=10,face="bold"),
       # legend.position = "none",
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
       strip.text.x = element_text(size = 5, face = "bold"),
       strip.text.y = element_text(size = 5, face = "bold"),
        legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.spacing=unit(0,"lines"))+
         scale_size(range = c(0.01,1))+
         scale_alpha(range = c(0.01,.8))
ggsave(plot=p,file=file.path(outputdir,"facet_scatter_alpha_ds_normlogistic.pdf"),width=6.5,height=1.8)