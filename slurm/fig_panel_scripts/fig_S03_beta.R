project_dir = "/share/crsp/lab/cellfate/mkarikom/DURIAN_paper_clean"
source("slurm/scrabble_helper_functions/library_cluster_metrics.R")

library(dplyr)
library(rmcorr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(scales)
library(xtable)

excludemodels = c()


# backupdircombine = "slurm/durian_pseudobulk/output/durian_pseudobulk_DownsampPB/output.final.baron,n_50"
# sparsityparam = "lambda"
# sparsityvar = "λ"

# dirnamescombine = list.files(backupdircombine,full.names=TRUE)
# df.orig = NULL # the dataframe of all output, to be filled
# idx = 1
# for(ii in 1:length(dirnamescombine)){
#   if(length(grep("_fit",dirnamescombine[ii]))>0){
#     backupdir = dirnamescombine[ii]
#     bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
#     dirnames = list.files(backupdir)
#     for(i in 1:length(dirnames)){
#       params1 = get_model_params(str_split(dirnames[i],pattern=",")[[1]])
#       params1b = get_model_params(bdname)
#       params = bind_rows(c(params1b,params1))

#       modeldirs = list.dirs(file.path(backupdir,dirnames[i]),recursive=FALSE,full.names=FALSE)
#       ind_dropmodel = grep("imputemodel_dropout",modeldirs)
#       if(file.exists(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"))){
#         droprate_backup = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$Dropout
#         for(j in 1:length(modeldirs)){
#           if(length(grep("imputemodel",modeldirs[j]))>0){
#             print(paste0("found impute folder ",modeldirs[j]))
#             logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
#             if(file.exists(logfn)){
#               print(paste0("found impute model ",logfn))
#               params2 = as.data.frame(get_model_params(modeldirs[j]))
#               dataparams = cbind(params,params2)
#               runlog = read.csv(logfn,row.names=1)

#               runlog = cbind(runlog,dataparams)
#               runlog$dropout_rate = droprate_backup

#               logfile = list.files(path=file.path(backupdir,dirnames[i],modeldirs[j]),pattern = "\\logdf.csv$", ignore.case = TRUE)
#               if(length(logfile)==0){
#                 ursmlog = list.files(path=file.path(backupdir,dirnames[i],modeldirs[j]),pattern = "^run.log$", ignore.case = TRUE)
#                 runtimeline = grep("Gibbs-EM finished in", readLines(file.path(backupdir,dirnames[i],modeldirs[j],ursmlog)), value = TRUE)
#                 runtime = as.numeric(stringr::str_extract(runtimeline, "\\d+\\.*\\d*"))
#               }else{
#                 logdata = read.csv(file.path(backupdir,dirnames[i],modeldirs[j],logfile),row.names=1)
#                 runtime = logdata$wallclock[nrow(logdata)]
#               }
#               runlog$runtime = rep(runtime,nrow(runlog))

#               if(length(grep("DURIAN",runlog$modelname))>0){
#                 runlog$modelfam=c("DURIAN")
#                 runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
#               }else if(length(grep("mtSCRABBLE",runlog$modelname))>0){
#                 runlog$modelfam=c("Control")
#                 runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
#               }else if(length(grep("SCRABBLE",runlog$modelname))>0){
#                 runlog$modelfam=c("Existing Methods")
#                 runlog$modelname = paste0(paste0(runlog$modelname,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
#               }else if(length(grep("dropout",runlog$modelname))>0){
#                 runlog$modelfam=c("Control")
#               }else{
#                 runlog$modelfam=c("Existing Methods")
#               }
#               df.orig = plyr::rbind.fill(df.orig,runlog)
#             }
#           }
#         }
#       }
#     }
#   }
# }
# df.orig$logRMSE = log(df.orig$RMSE)
# df.orig$logENORM = log(df.orig$ENORM)
# 
# df.orig$dropoutlevel = paste0(sparsityvar,"=",df.orig[,sparsityparam])
# df.filter = df.orig %>% dplyr::filter(!modelname %in% excludemodels)
# df_downsamp = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
# df_downsamp$strategy = "Down-Sampling"

backupdircombine = "/share/crsp/lab/cellfate/mkarikom/DURIAN_paper_clean/slurm/durian_pseudobulk/output/pseudobulk_SimPBbetas,gProb_0.1-0.1-0.5-0.3,deProb_0.3-0.3-0.3-0.3,bLoc_0.1,bScale_0.2,dLoc_0.5,dScale_0.5/output.final.splatter,n_5"
sparsityparam = "dmid"
sparsityvar = "x0"

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

      modeldirs = list.dirs(file.path(backupdir,dirnames[i]),recursive=FALSE,full.names=FALSE)
      ind_dropmodel = grep("imputemodel_dropout",modeldirs)
      if(file.exists(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"))){
        droprate_backup = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$Dropout
        for(j in 1:length(modeldirs)){
          if(length(grep("imputemodel",modeldirs[j]))>0){
            print(paste0("found impute folder ",modeldirs[j]))
            logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
            if(file.exists(logfn)){
              print(paste0("found impute model ",logfn))
              params2 = as.data.frame(get_model_params(modeldirs[j]))
              dataparams = cbind(params,params2)
              runlog = read.csv(logfn,row.names=1)

              runlog = cbind(runlog,dataparams)
              runlog$dropout_rate = droprate_backup

              logfile = list.files(path=file.path(backupdir,dirnames[i],modeldirs[j]),pattern = "\\logdf.csv$", ignore.case = TRUE)
              if(length(logfile)==0){
                ursmlog = list.files(path=file.path(backupdir,dirnames[i],modeldirs[j]),pattern = "^run.log$", ignore.case = TRUE)
                runtimeline = grep("Gibbs-EM finished in", readLines(file.path(backupdir,dirnames[i],modeldirs[j],ursmlog)), value = TRUE)
                runtime = as.numeric(stringr::str_extract(runtimeline, "\\d+\\.*\\d*"))
              }else{
                logdata = read.csv(file.path(backupdir,dirnames[i],modeldirs[j],logfile),row.names=1)
                runtime = logdata$wallclock[nrow(logdata)]
              }
              runlog$runtime = rep(runtime,nrow(runlog))

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
              df.orig = plyr::rbind.fill(df.orig,runlog)
            }
          }
        }
      }
    }
  }
}
df.orig$logRMSE = log(df.orig$RMSE)
df.orig$logENORM = log(df.orig$ENORM)


print("sleeping")
Sys.sleep(1000)


##################################################
#
# plot DURIAN BETA
#
##################################################

outputdir = file.path("slurm","fig_panel_scripts","fig_S03_beta")
dir.create(outputdir,recursive=TRUE)

df.orig$dropoutlevel = paste0(sparsityvar,"=",unlist(lapply(strsplit(df.orig[,sparsityparam],":"),tail,n=1)))
df.filter = df.orig %>% dplyr::filter(sA==1e-2 & sB %in% c(0,1e-7,1e-6,1e-5,1e-4) & modelfam=="DURIAN")
df_sim = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
df_sim$strategy = "Simulation"

df = df_sim
write.csv(df,file.path(outputdir,"df.csv"))

df$mean_dropout = round(df$mean_dropout,digits=2)

melt_tmp = reshape2::melt(df,measure.vars=c("logRMSE","logENORM"))

p = ggplot(melt_tmp,
        aes(
          x=as.factor(modelname), 
          y=value,
          fill=as.factor(sB),
          color=as.factor(sB),
          alpha=0.5)) + 
      geom_boxplot(
        aes(
          x=as.factor(modelname),
          y=value,
          fill=as.factor(sB))) +
      facet_nested(variable ~ imputemodel,scale="free") 

set_palette(p, palette1)
p = p +
    theme_bw() +
    rotate_x_text(30) +
    labs(color="β",y="Error",alpha="del",fill="del") +
    theme(axis.text=element_text(size=10,face="bold"),
      axis.title.y=element_text(size=13,face="bold"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
      # legend.position = "none",
      panel.border=element_rect(colour="black",size=0.5,fill=NA),
        strip.text.x = element_text(size = 8, face = "bold", angle = 0),
        strip.text.y = element_text(size = 10, face = "bold", angle = 270),      
        legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=8,face="bold"), #change legend title font size
      legend.text = element_text(size=8,face="bold"),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing=unit(0,"lines"))

ggsave(plot=p,file=file.path(outputdir,"betas.pdf"),width=5,height=5)