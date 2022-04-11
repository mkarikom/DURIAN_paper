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

excludemodels = c()
quicksims = c(1:5)
quickdparams = c("0:0:0:0:6.5:6.5:6.5:6.5",1e-6)
umap_min_dist = 0.7

outputdir = file.path("slurm","fig_panel_scripts","fig03_ma")
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

          for(j in 1:length(modeldirs)){
            if(length(grep("imputemodel",modeldirs[j]))>0){
              print(paste0("found impute folder ",modeldirs[j]))
              logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
              imputedfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputed_C.csv")
              if(file.exists(imputedfn) & file.exists(logfn)){
                imputedC = read.csv(imputedfn,row.names=1)

                diffC = imputedC-trueC

                normDiffC = apply(diffC, 2, function(x){sqrt(sum(x^2))})

                pDataC = read.csv(file.path(backupdir,dirnames[i],"pDataC.csv"),row.names=1)
                pDataC = cbind(pDataC,as.data.frame(normDiffC))

                ylims=c(-7,7)
                mthresh=c(-2,2)

                print("filling MA table")
                A = 0.5*rowMeans(log2(trueC+1) + log2(imputedC+1))
                # M = rowMeans(log2(trueC+1) - log2(imputedC+1))
                M = rowMeans(log2(imputedC+1) - log2(trueC+1))
                mrate = abs(M/ylims[1])
                plotdf = data.frame(A=A,M=M)
                inbounds=c()
                for(iii in 1:length(M)){
                    if(M[iii] <= mthresh[1]){
                        inbounds = c(inbounds,"Negative")
                    }else if(M[iii] >= mthresh[2]){
                        inbounds = c(inbounds,"Positive")
                    }else{
                        inbounds = c(inbounds,"Below threshold")
                    }
                }
                plotdf$inbounds = as.factor(inbounds)
                plotdf$mrate = mrate
                
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

                runlog_ma = cbind(runlog,plotdf)

                print("calling rbind")
                df.orig = plyr::rbind.fill(df.orig,runlog_ma)
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
          for(j in 1:length(modeldirs)){
            print(modeldirs[j])
            if(length(grep("imputemodel",modeldirs[j]))>0){
              print(paste0("found impute folder ",modeldirs[j]))
              logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
              imputedfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputed_C.csv")
              if(file.exists(imputedfn) & file.exists(logfn)){
                imputedC = read.csv(imputedfn,row.names=1)

                diffC = imputedC-trueC

                normDiffC = apply(diffC, 2, function(x){sqrt(sum(x^2))})

                pDataC = read.csv(file.path(backupdir,dirnames[i],"pDataC.csv"),row.names=1)
                pDataC = cbind(pDataC,as.data.frame(normDiffC))

                ylims=c(-7,7)
                mthresh=c(-2,2)

                print("filling MA table")
                A = 0.5*rowMeans(log2(trueC+1) + log2(imputedC+1))
                # M = rowMeans(log2(trueC+1) - log2(imputedC+1))
                M = rowMeans(log2(imputedC+1) - log2(trueC+1))
                mrate = abs(M/ylims[1])
                plotdf = data.frame(A=A,M=M)
                inbounds=c()
                for(iii in 1:length(M)){
                    if(M[iii] <= mthresh[1]){
                        inbounds = c(inbounds,"Negative")
                    }else if(M[iii] >= mthresh[2]){
                        inbounds = c(inbounds,"Positive")
                    }else{
                        inbounds = c(inbounds,"Below threshold")
                    }
                }
                plotdf$inbounds = as.factor(inbounds)
                plotdf$mrate = mrate

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

                runlog_ma = cbind(runlog,plotdf)

                print("calling rbind")
                df.orig = plyr::rbind.fill(df.orig,runlog_ma)
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

df$mean_dropout = round(df$mean_dropout,digits=3)
nmodels = length(unique(df$modelname))
palette1 = scales::hue_pal()(nmodels)
names(palette1) = sort(unique(df$modelname))

write.csv(df,file.path(outputdir,"df.csv"))

df_tmp = df %>% filter(strategy=="Down-Sampling",sA %in% c(1e-2,NA),sim==1)
p=ggplot(df_tmp %>% sample_n(6000),aes(x=A,y=M,color=inbounds,alpha=mrate))+
# p=ggplot(df_tmp,aes(x=A,y=M,color=inbounds,alpha=mrate))+
      geom_point(size=0.25)+
        geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5)+
        geom_hline(yintercept=2, linetype="dashed",color = "red", size=1)+
        geom_hline(yintercept=-2, linetype="dashed",color = "blue", size=1)+
        ylim(ylims)+
        scale_colour_manual(values = c("grey30","blue", "red"))+
        theme_bw()+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
        xlab("Average counts")+
        ylab("Log ratio")+
      facet_nested(strategy~modelfam + imputemodel)
p = p +
     theme_bw() +
     rotate_x_text(90) +
     labs(color="Error Sign") +
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
       panel.spacing=unit(0,"lines")) +
      guides(color = guide_legend(override.aes = list(size=3)))

ggsave(plot=p,file=file.path(outputdir,"facet_ma_alpha_ds.pdf"),width=6.5,height=1.8)

df_tmp = df %>% filter(strategy=="Simulation",sA %in% c(1e-2,NA),sim==1)
p=ggplot(df_tmp,aes(x=A,y=M,color=inbounds,alpha=mrate))+
      geom_point(size=0.25)+
        geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5)+
        geom_hline(yintercept=2, linetype="dashed",color = "red", size=1)+
        geom_hline(yintercept=-2, linetype="dashed",color = "blue", size=1)+
        ylim(ylims)+
        scale_colour_manual(values = c("grey30","blue", "red"))+
        theme_bw()+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
        xlab("Average counts")+
        ylab("Log ratio")+
      facet_nested(strategy~modelfam + imputemodel)
p = p +
     theme_bw() +
     rotate_x_text(90) +
     labs(color="Error Sign") +
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
       panel.spacing=unit(0,"lines")) + 
       guides(color = guide_legend(override.aes = list(size=3)))
ggsave(plot=p,file=file.path(outputdir,"facet_ma_alpha_sim.pdf"),width=6.5,height=1.8)

