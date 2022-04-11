###############################################
# delete this
###############################################
# Sys.setenv(PBULKDIRCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/output.final.baron")
# Sys.setenv(SUMMARYCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/realdata_scripts/test/summaryfinalcombine")
# Sys.setenv(SPARSITY_PARAM="lambda")
###############################################


source("slurm/scrabble_helper_functions/library_cluster_metrics.R")
source("slurm/scrabble_helper_functions/library_extra_plots.R")
source("slurm/scrabble_helper_functions/library_scrabble_clusterMetrics_clValid.R")

# saved benchmark directory
backupdirs = c("slurm/pseudobulk_OuterStats_clValidStabGprob,gProb_0.1-0.1-0.5-0.3,deProb_0.3-0.3-0.3-0.3,bLoc_0.1,bScale_0.2,dLoc_0.5,dScale_0.5/output.final.splatter,n_50","slurm/durian_pseudobulk_BaronOuterStatsAllNested/output.final.baron,n_50")
# output directory
outputdir = "slurm/fig_panel_scripts/figS01n50"
# the model parameter controlling the sparsity
sparsityparams = c("lambda","dmid")
excludemodels = c()
quicksims = 1:30 # for some reason, r plotting will not work if the dataframe is too large
quickdparams = c("0:0:0:0:6.5:6.5:6.5:6.5","0:0:0:0:5.5:5.5:5.5:5.5",1e-6,5e-6)

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggh4x)

plotdir = file.path(outputdir,"iter_heatmaps")

dir.create(plotdir,recursive=TRUE)
dir.exists(plotdir)
parseargs <- function(input){
  if(!is.na(strtoi(input))){
    output=strtoi(input)
  }else if(!is.na(as.numeric(input))){
    output=as.numeric(input)
  }else{
    output=input
  }
  output
}

foundfolders=list()
missingfolders=list()
countfound = 0
countmissing = 0


df = NULL # the dataframe of all output, to be filled
for(backupdir in backupdirs){
  dirnamescombine = list.files(backupdir,full.names=TRUE)
  idx = 1
  for(ii in 1:length(dirnamescombine)){
    if(length(grep("_fit",dirnamescombine[ii]))>0){
      backupdir = dirnamescombine[ii]
      bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
      paramvec.orig = str_split(bdname,pattern=",")[[1]]
      dirnames = list.files(backupdir)
      for(i in 1:length(dirnames)){
        print(paste0("dirname ",dirnames[i]))
        paramvec = str_split(dirnames[i],pattern=",")[[1]]
        paramvec = c(paramvec,paramvec.orig)
        params = list()
        for(j in 1:length(paramvec)){
          newparam = str_split(paramvec[j],pattern="_")[[1]]
          params[[eval(newparam[1])]] = parseargs(newparam[2])
        }
        if(is.null(params$dmid)){
          sparsityparam = "lambda"
        }else if(is.null(params$lambda)){
          sparsityparam = "dmid"
        }
        if(params$sim %in% quicksims & params[sparsityparam] %in% quickdparams){

          modeldirs = list.dirs(file.path(backupdir,dirnames[i]),recursive=FALSE,full.names=FALSE)
          ind_dropmodel = grep("imputemodel_dropout",modeldirs)
          if(file.exists(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"))){
            droprate_backup = read.csv(file.path(backupdir,dirnames[i],modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$Dropout
            for(j in 1:length(modeldirs)){
              if(length(grep("imputemodel",modeldirs[j]))>0){
                if(length(grep("DURIAN",modeldirs[j]))>0){

                  modelparamvec = str_split(modeldirs[j],pattern=",")[[1]]
                  modelparams = list()
                  for(jj in 1:length(modelparamvec)){
                    newparam = str_split(modelparamvec[jj],pattern="_")[[1]]
                    modelparams[[eval(newparam[1])]] = parseargs(newparam[2])
                  }
                  deconv_method = strsplit(modelparams$imputemodel,"\\.")[[1]][2]
                  model_fn = file.path(backupdir,dirnames[i],modeldirs[j])
                  logfn = file.path(model_fn,paste0(modelparams$imputemodel,"_logdf.csv"))
                  if(file.exists(logfn)){
                    countfound=countfound+1
                    foundfolders[[countfound]] = model_fn
                    logdf = read.csv(logfn,row.names=1)
                    logdf$dropout_rate = droprate_backup
                    # add model info to the table
                    # # paramdf2 = unlist(params) %>% slice(rep(1:n(), each = dim(logdf)[1])) %>% data.frame()    
                    paramdf = as.data.frame(matrix(rep(params,dim(logdf)[1]),nrow=dim(logdf)[1],byrow=TRUE))
                    names(paramdf) = names(params)
                    for(icol in 1:ncol(paramdf)){
                      paramdf[,icol] = unlist(paramdf[,icol])
                    }
                    paramdf$durian_rmse_final = rep(logdf$durian_rmse[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$durian_errnorm_final = rep(logdf$errnorm[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$durian_cellcor_final = rep(logdf$durian_cellcor[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$durian_genecor_final = rep(logdf$durian_genecor[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$durian_mean_cellcor_final = rep(logdf$durian_mean_cellcor[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$durian_mean_genecor_final = rep(logdf$durian_mean_genecor[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$deconv_rmse_final = rep(logdf$deconv_rmse[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$deconv_cor_celltype_final = rep(logdf$deconv_cor_celltype[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$deconv_cor_mean_celltype_final = rep(logdf$deconv_cor_mean_celltype[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$deconv_cor_bulksample_final = rep(logdf$deconv_cor_bulksample[dim(logdf)[1]],dim(paramdf)[1])
                    paramdf$deconv_cor_mean_bulksample_final = rep(logdf$deconv_cor_mean_bulksample[dim(logdf)[1]],dim(paramdf)[1])

                    paramdf$durian_rmse_init = rep(logdf$durian_rmse[1],dim(paramdf)[1])
                    paramdf$durian_errnorm_init = rep(logdf$errnorm[1],dim(paramdf)[1])
                    paramdf$durian_cellcor_init = rep(logdf$durian_cellcor[1],dim(paramdf)[1])
                    paramdf$durian_genecor_init = rep(logdf$durian_genecor[1],dim(paramdf)[1])
                    paramdf$durian_mean_cellcor_init = rep(logdf$durian_mean_cellcor[1],dim(paramdf)[1])
                    paramdf$durian_mean_genecor_init = rep(logdf$durian_mean_genecor[1],dim(paramdf)[1])
                    paramdf$deconv_rmse_init = rep(logdf$deconv_rmse[2],dim(paramdf)[1])
                    paramdf$deconv_cor_celltype_init = rep(logdf$deconv_cor_celltype[2],dim(paramdf)[1])
                    paramdf$deconv_cor_mean_celltype_init = rep(logdf$deconv_cor_mean_celltype[2],dim(paramdf)[1])
                    paramdf$deconv_cor_bulksample_init = rep(logdf$deconv_cor_bulksample[2],dim(paramdf)[1])
                    paramdf$deconv_cor_mean_bulksample_init = rep(logdf$deconv_cor_mean_bulksample[2],dim(paramdf)[1])

                    # add delta info to the table
                    paramdf$durian_rmse_delta = rep(paramdf$durian_rmse_final[1] - paramdf$durian_rmse_init[1],dim(paramdf)[1])
                    paramdf$durian_logrmse_delta = rep(log(paramdf$durian_rmse_final[1]) - log(paramdf$durian_rmse_init[1]),dim(paramdf)[1])
                    paramdf$durian_errnorm_delta = rep(paramdf$durian_errnorm_final[1] - paramdf$durian_errnorm_init[1],dim(paramdf)[1])
                    paramdf$durian_logerrnorm_delta = rep(log(paramdf$durian_errnorm_final[1]) - log(paramdf$durian_errnorm_init[1]),dim(paramdf)[1])
                    paramdf$durian_cellcor_delta = rep(paramdf$durian_cellcor_final[1] - paramdf$durian_cellcor_init[1],dim(paramdf)[1])
                    paramdf$durian_genecor_delta = rep(paramdf$durian_genecor_final[1] - paramdf$durian_genecor_init[1],dim(paramdf)[1])
                    paramdf$durian_mean_cellcor_delta = rep(paramdf$durian_mean_cellcor_final[1] - paramdf$durian_mean_cellcor_init[1],dim(paramdf)[1])
                    paramdf$durian_mean_genecor_delta = rep(paramdf$durian_mean_genecor_final[1] - paramdf$durian_mean_genecor_init[1],dim(paramdf)[1])
                    paramdf$deconv_rmse_delta = rep(paramdf$deconv_rmse_final[1] - paramdf$deconv_rmse_init[1],dim(paramdf)[1])
                    paramdf$deconv_cor_celltype_delta = rep(paramdf$deconv_cor_celltype_final[1] - paramdf$deconv_cor_celltype_init[1],dim(paramdf)[1])
                    paramdf$deconv_cor_mean_celltype_delta = rep(paramdf$deconv_cor_mean_celltype_final[1] - paramdf$deconv_cor_mean_celltype_init[1],dim(paramdf)[1])
                    paramdf$deconv_cor_bulksample_delta = rep(paramdf$deconv_cor_bulksample_final[1] - paramdf$deconv_cor_bulksample_init[1],dim(paramdf)[1])
                    paramdf$deconv_cor_mean_bulksample_delta = rep(paramdf$deconv_cor_mean_bulksample_final[1] - paramdf$deconv_cor_mean_bulksample_init[1],dim(paramdf)[1])


                    paramdf$durian_idx = rep(idx,dim(paramdf)[1])
                    paramdf$deconv_method = rep(deconv_method,dim(paramdf)[1])
                    dfnew = cbind(logdf,paramdf)

                    # normalized metrics
                    dfnew$durian_logrmse_diffRatioInit = abs((log(dfnew$durian_rmse) - log(dfnew$durian_rmse_final[1]))/dfnew$durian_logrmse_delta[1])
                    dfnew$durian_rmse_diffRatioInit = abs((dfnew$durian_rmse - dfnew$durian_rmse_final[1])/dfnew$durian_rmse_delta[1])
                    dfnew$durian_logerrnorm_diffRatioInit = abs((log(dfnew$errnorm) - log(dfnew$durian_errnorm_final[1]))/dfnew$durian_logerrnorm_delta[1])
                    dfnew$durian_errnorm_diffRatioInit = abs((dfnew$errnorm - dfnew$durian_errnorm_final[1])/dfnew$durian_errnorm_delta[1])
                    dfnew$durian_rmse_ratioInit = dfnew$durian_rmse / dfnew$durian_rmse_init[1]
                    dfnew$durian_errnorm_ratioInit = dfnew$errnorm / dfnew$durian_errnorm_init[1]
                    # dfnew$durian_cellcor_norm = dfnew$durian_cellcor / dfnew$durian_cellcor_final[1]
                    # dfnew$durian_genecor_norm = dfnew$durian_genecor / dfnew$durian_genecor_final[1]
                    # dfnew$durian_mean_cellcor_norm = dfnew$durian_mean_cellcor / dfnew$durian_mean_cellcor_final[1]
                    # dfnew$durian_mean_genecor_norm = dfnew$durian_mean_genecor / dfnew$durian_mean_genecor_final[1]
                    dfnew$deconv_rmse_diffRatioInit = abs((dfnew$deconv_rmse - dfnew$deconv_rmse_final[1])/dfnew$deconv_rmse_delta[1])
                    dfnew$deconv_rmse_ratioInit = dfnew$deconv_rmse /dfnew$deconv_rmse_init[1]
                    # dfnew$deconv_cor_celltype_norm = dfnew$deconv_cor_celltype / dfnew$deconv_cor_celltype_final[1]
                    # dfnew$deconv_cor_mean_celltype_norm = dfnew$deconv_cor_mean_celltype / dfnew$deconv_cor_mean_celltype_final[1]
                    # dfnew$deconv_cor_bulksample_norm = dfnew$deconv_cor_bulksample / dfnew$deconv_cor_bulksample_final[1]
                    # dfnew$deconv_cor_mean_bulksample_norm = dfnew$deconv_cor_mean_bulksample / dfnew$deconv_cor_mean_bulksample_final[1]
                    dfnew$scrabbleLoss_ratioInit = dfnew$scrabbleLoss / dfnew$scrabbleLoss[2]

                    df = plyr::rbind.fill(df,dfnew)
                    idx = idx+1
                  }else{
                    countmissing=countmissing+1
                    missingfolders[[countmissing]] = model_fn
                    print(paste0("missing output for: ",model_fn))
                  }
                }
              }
            }
          }
        }

      }
    }
  }
}

# add sparsity
df.ds = df %>% filter(!is.na(lambda)) %>% group_by(lambda) %>% mutate(mean_dropout = mean(dropout_rate))
df.sim = df %>% filter(!is.na(dmid)) %>% group_by(dmid) %>% mutate(mean_dropout = mean(dropout_rate))
df.sim$strategy = "Simulation"
df.ds$strategy = "Down-Sampling"
df = plyr::rbind.fill(df.ds,df.sim)

df$iter = as.integer(df$iter)
df$mean_dropout = round(df$mean_dropout,digits=3)
df$mean_dropout = as.factor(df$mean_dropout)
df$deconv_method = as.factor(df$deconv_method)
df$durian_logrmse = log(df$durian_rmse)
df$durian_logerrnorm = log(df$errnorm)
write.csv(df,file=file.path(plotdir,"df.csv"))
df.orig = df

# # df[,sparsityparam] = paste0(sparsityparam,"=",df[,sparsityparam])
# # df[,sparsityparam] = as.factor(df[,sparsityparam])
# df[,"mean_dropout"] = paste0("dropout=",df[,"mean_dropout"])
# df[,"deconv_method"] = paste0("D-step=",df[,"deconv_method"])

# paletteLvs = unique(as.vector(as.matrix(df[,c("mean_dropout","deconv_method")])))
# df[,sparsityparam] = factor(df[,sparsityparam],levels=paletteLvs)
# df[,"deconv_method"] = factor(df[,"deconv_method"],levels=paletteLvs)
# if(length(grep("lambda",names(df)))>0){
#   myPal = brewer.pal(n = length(paletteLvs), name = "Set1")
# }else{
#   myPal = brewer.pal(n = length(paletteLvs), name = "Dark2")
# }
# names(myPal) = paletteLvs
# # myPal = c("purple","orange","yellow","brown")
# # names(myPal) = paletteLvs


yvars = c("Gene","Cell","logGene","logCell","RMSE","logRMSE","ENORM", "logENORM","MeanGene","MeanCell","logMeanGene","logMeanCell")

nmodels = length(unique(df$modelname))
ncompare = choose(nmodels,2)

m2durian = c("durian_logrmse","durian_logerrnorm","durian_rmse_diffRatioInit","durian_errnorm_diffRatioInit","durian_rmse_ratioInit","durian_errnorm_ratioInit","durian_mean_cellcor")
m2deconv = c("deconv_cor_mean_celltype","deconv_rmse_diffRatioInit","deconv_rmse_ratioInit")

plotdir1 = file.path(plotdir,"nolineMulti")
dir.create(plotdir1)
plotdir2 = file.path(plotdir,"nolineMultiFlip")
dir.create(plotdir2)

write.csv(df,file.path(plotdir,"backup_df.csv"))

df.filt = df %>% filter(iter>0,!is.na(durian_rmse),!is.na(deconv_rmse),!is.na(durian_rmse_diffRatioInit),!is.na(deconv_rmse_diffRatioInit))
df.melt = reshape2::melt(df.filt,measure.vars=c("durian_rmse_diffRatioInit","deconv_rmse_diffRatioInit"))
df.melt.imp = df.melt %>% filter(variable=="durian_rmse_diffRatioInit") %>% mutate(aspect = "Imputation")
df.melt.deconv = df.melt %>% filter(variable=="deconv_rmse_diffRatioInit") %>% mutate(aspect = "Deconvolution")
df.melt = rbind(df.melt.imp,df.melt.deconv)

df.ds = df.melt %>% filter(!is.na(lambda)) %>% mutate(strategy = "Down-Sampling")
df.sim = df.melt %>% filter(!is.na(dmid)) %>% mutate(strategy = "Simulation")
df.melt = rbind(df.ds,df.sim)


df.melt.filt = df.melt %>% filter(aspect=="Imputation")
# filter outliers
upper_bound <- quantile(df.melt.filt$value, 0.9750)
use_ind <- which(df.melt.filt$value < upper_bound)

p = ggplot(
  df.melt.filt[use_ind,],
  aes(x = iter, y=scrabbleLoss,z = value, color=deconv_method)) +
        stat_summary_hex(
          fun=mean,bins=5,aes(color=deconv_method),
          size=1) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=1,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Greens",guide = guide_colourbar(title="Impute\nErr")) +
        facet_nested(aspect ~ strategy+deconv_method,scales="free")+
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Deconv\nMethod") +
        scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1))+
        theme(
        axis.title=element_text(size = 8,face="bold"),
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
        axis.text = element_text(size = 8,face="bold"),
        strip.text.x = element_text(size = 8, face = "bold", angle = 0),
        strip.text.y = element_text(size = 7, face = "bold", angle = 270),
      legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0,"lines")) 
ggsave(plot=p+theme(legend.position = "none"),file=file.path(plotdir,"summary_heatmap_impute.pdf"),width=6.5,height=2)
ggsave(plot=p,file=file.path(plotdir,"summary_heatmap_impute_leg.pdf"),width=6.5,height=2)

saveRDS(p,file=file.path(plotdir,"summary_heatmap_impute.RDS"))

df.melt.filt = df.melt %>% filter(aspect=="Deconvolution")
# filter outliers
upper_bound <- quantile(df.melt.filt$value, 0.95)
use_ind <- which(df.melt.filt$value < upper_bound)

p = ggplot(
  df.melt.filt[use_ind,],
  aes(x = iter, y=scrabbleLoss,z = value, color=deconv_method)) +
        stat_summary_hex(
          fun=mean,bins=5,aes(color=deconv_method),
          size=1) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=1,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Blues",guide = guide_colourbar(title="Deconv\nErr")) +
        facet_nested(aspect ~ strategy+deconv_method,scales="free")+
        scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1))+
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Deconv\nMethod") +
        theme(
        axis.title=element_text(size = 8,face="bold"),
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
        axis.text = element_text(size = 8,face="bold"),
        strip.text.x = element_text(size = 8, face = "bold", angle = 0),
        strip.text.y = element_text(size = 7, face = "bold", angle = 270),
      legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0,"lines")) 
ggsave(plot=p+theme(legend.position = "none"),file=file.path(plotdir,"summary_heatmap_deconv.pdf"),width=6.5,height=2)
ggsave(plot=p,file=file.path(plotdir,"summary_heatmap_deconv_leg.pdf"),width=6.5,height=2)
saveRDS(p,file=file.path(plotdir,"summary_heatmap_deconv.RDS"))

df.melt.filt = df.melt %>% filter(aspect=="Imputation")
# filter outliers
upper_bound <- quantile(df.melt.filt$value, 0.95)
use_ind <- which(df.melt.filt$value < upper_bound)

p = ggplot(
  df.melt.filt[use_ind,],
  aes(x = iter, y=scrabbleLoss,z = value, color=mean_dropout)) +
        stat_summary_hex(
          bins=5,fun=mean,aes(color=mean_dropout),
          size=1) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=1,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Greens",guide = guide_colourbar(title="Impute\nErr")) +
        facet_nested(aspect ~ strategy+mean_dropout,scales="free")+
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Deconv\nMethod") +
        theme(
        axis.title=element_text(size = 8,face="bold"),
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
        axis.text = element_text(size = 8,face="bold"),
        strip.text.x = element_text(size = 8, face = "bold", angle = 0),
        strip.text.y = element_text(size = 7, face = "bold", angle = 270),
      legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0,"lines")) 
ggsave(plot=p+theme(legend.position = "none"),file=file.path(plotdir,"summary_heatmap_impute_dropfacet.pdf"),width=6.5,height=2)
ggsave(plot=p,file=file.path(plotdir,"summary_heatmap_impute_leg_dropfacet.pdf"),width=6.5,height=2)

df.melt.filt = df.melt %>% filter(aspect=="Deconvolution")
# filter outliers
upper_bound <- quantile(df.melt.filt$value, 0.95)
use_ind <- which(df.melt.filt$value < upper_bound)

p = ggplot(
  df.melt.filt[use_ind,],
  aes(x = iter, y=scrabbleLoss,z = value, color=mean_dropout)) +
        stat_summary_hex(
          bins=5,fun=mean,aes(color=mean_dropout),
          size=1) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=1,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Blues",guide = guide_colourbar(title="Deconv\nErr")) +
        facet_nested(aspect ~ strategy+mean_dropout,scales="free")+
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Deconv\nMethod") +
        theme(
        axis.title=element_text(size = 8,face="bold"),
       panel.border=element_rect(colour="black",size=0.5,fill=NA),
        axis.text = element_text(size = 8,face="bold"),
        strip.text.x = element_text(size = 8, face = "bold", angle = 0),
        strip.text.y = element_text(size = 7, face = "bold", angle = 270),
      legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=7,face="bold"), #change legend title font size
      legend.text = element_text(size=5,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0,"lines")) 
ggsave(plot=p+theme(legend.position = "none"),file=file.path(plotdir,"summary_heatmap_deconv_dropfacet.pdf"),width=6.5,height=2)
ggsave(plot=p,file=file.path(plotdir,"summary_heatmap_deconv_leg_dropfacet.pdf"),width=6.5,height=2)

