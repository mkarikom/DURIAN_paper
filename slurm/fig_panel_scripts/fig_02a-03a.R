###############################################
# delete this
###############################################
# Sys.setenv(PBULKDIRCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/output.final.baron")
# Sys.setenv(SUMMARYCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/realdata_scripts/test/summaryfinalcombine")
# Sys.setenv(SPARSITY_PARAM="lambda")
###############################################

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

backupdircombine = "slurm/durian_pseudobulk_BaronOuterStatsAllNested/output.final.baron,n_50"
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

df.orig$dropoutlevel = paste0(sparsityvar,"=",df.orig[,sparsityparam])
df.filter = df.orig %>% filter(!modelname %in% excludemodels)
df_downsamp = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
df_downsamp$strategy = "Down-Sampling"

backupdircombine = "slurm/pseudobulk_OuterStats_clValidStabGprob,gProb_0.1-0.1-0.5-0.3,deProb_0.3-0.3-0.3-0.3,bLoc_0.1,bScale_0.2,dLoc_0.5,dScale_0.5/output.final.splatter,n_50"
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



df.orig$dropoutlevel = paste0(sparsityvar,"=",unlist(lapply(strsplit(df.orig[,sparsityparam],":"),tail,n=1)))
df.filter = df.orig %>% filter(!modelname %in% excludemodels)
df_sim = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
df_sim$strategy = "Simulation"


# ##################################################
# #
# # plot 1e-3,1e-5
# #
# #
# ##################################################

# df = plyr::rbind.fill(df_sim,df_downsamp)
# df$mean_dropout = round(df$mean_dropout,digits=3)
# df =  filter(df,(sB==1e-3 & sG==1e-5) | (sB==1e-6 & sG==1e-4) | (is.na(sB) & is.na(sG)))

# nmodels = length(unique(df$modelname))
# palette1 = scales::hue_pal()(nmodels)
# names(palette1) = sort(unique(df$modelname))


# saveRDS(palette1,file.path("slurm","fig_panel_scripts","palette1.RDS"))

# outputdir = file.path("slurm","fig_panel_scripts","fig02n50.1e-3.1e-5")
# dir.create(outputdir,recursive=TRUE)
# melt_tmp = reshape2::melt(df,measure.vars=c("logRMSE"))

# p = ggplot(melt_tmp,
#         aes(
#           x=as.factor(modelname), 
#           y=value,
#           fill=as.factor(imputemodel),
#           color=as.factor(imputemodel),
#           alpha=0.5)) + 
#       geom_boxplot(
#         aes(
#           x=as.factor(modelname),
#           y=value,
#           fill=as.factor(imputemodel))) +
#       # facet_grid(rows=vars(mean_dropout),cols=vars(modelfam),scale="free") 
#       facet_nested(strategy + mean_dropout ~ modelfam,scale="free") 
#       # stat_summary(fun=mean, geom="point", shape=18, size=1,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
#       #                 position = position_dodge(0.75),stroke=1.5)+
#       # stat_summary(fun=mean, geom="text", size=1,vjust=-0.5,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
#       #                 guides(fill = FALSE)
# set_palette(p, palette1)
# p = p +
#     theme_bw() +
#     rotate_x_text(30) +
#     labs(color="Model Family",y="Mean Error") +
#     theme(axis.text=element_text(size=5,face="bold"),
#       axis.title.y=element_text(size=10,face="bold"),
#       axis.title.x=element_blank(),
#       # legend.position = "none",
#       panel.border=element_rect(colour="black",size=0.5,fill=NA),
#         strip.text.x = element_text(size = 8, face = "bold", angle = 0),
#         strip.text.y = element_text(size = 7, face = "bold", angle = 270),      
#         legend.key.size = unit(5, 'mm'), #change legend key size
#       legend.key.height = unit(4, 'mm'), #change legend key height
#       legend.key.width = unit(4, 'mm'), #change legend key width
#       legend.title = element_text(size=7,face="bold"), #change legend title font size
#       legend.text = element_text(size=5,face="bold"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.spacing=unit(0,"lines"))

# ggsave(plot=p,file=file.path(outputdir,"logRMSE.pdf"),width=6.5,height=6)


# outputdir = file.path("slurm","fig_panel_scripts","fig03n50.1e-3.1e-5")
# dir.create(outputdir,recursive=TRUE)
# melt_tmp = reshape2::melt(df,measure.vars=c("logENORM"))

# p = ggplot(melt_tmp,
#         aes(
#           x=as.factor(modelname), 
#           y=value,
#           fill=as.factor(imputemodel),
#           color=as.factor(imputemodel),
#           alpha=0.5)) + 
#       geom_boxplot(
#         aes(
#           x=as.factor(modelname),
#           y=value,
#           fill=as.factor(imputemodel))) +
#       # facet_grid(rows=vars(mean_dropout),cols=vars(modelfam),scale="free") 
#       facet_nested(strategy + mean_dropout ~ modelfam,scale="free") 
#       # stat_summary(fun=mean, geom="point", shape=18, size=1,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
#       #                 position = position_dodge(0.75),stroke=1.5)+
#       # stat_summary(fun=mean, geom="text", size=1,vjust=-0.5,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
#       #                 guides(fill = FALSE)
# set_palette(p, palette1)
# p = p +
#     theme_bw() +
#     rotate_x_text(30) +
#     labs(color="Model Family",y="Max Error") +
#     theme(axis.text=element_text(size=5,face="bold"),
#       axis.title.y=element_text(size=10,face="bold"),
#       axis.title.x=element_blank(),
#       # legend.position = "none",
#       panel.border=element_rect(colour="black",size=0.5,fill=NA),
#         strip.text.x = element_text(size = 8, face = "bold", angle = 0),
#         strip.text.y = element_text(size = 7, face = "bold", angle = 270),
#         legend.key.size = unit(5, 'mm'), #change legend key size
#       legend.key.height = unit(4, 'mm'), #change legend key height
#       legend.key.width = unit(4, 'mm'), #change legend key width
#       legend.title = element_text(size=7,face="bold"), #change legend title font size
#       legend.text = element_text(size=5,face="bold"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.spacing=unit(0,"lines"))

# ggsave(plot=p,file=file.path(outputdir,"logENORM.pdf"),width=6.5,height=6)

# ################################################
# #
# # stats
# #
# ################################################



# df.tab.rmse = group_by(df, strategy, modelname) %>%
#   summarise(
#     count = n(),
#     mean = mean(logRMSE, na.rm = TRUE),
#     sd = sd(logRMSE, na.rm = TRUE)
#   ) %>% arrange(mean,.by_group = TRUE)
# print(df.tab.rmse,n=100)
# wilcox.tab.rmse = pairwise.wilcox.test(df$logRMSE, df$modelname,p.adjust.method = "BH")

# df.tab.enorm = group_by(df, strategy, modelname) %>%
#   summarise(
#     count = n(),
#     mean = mean(logENORM, na.rm = TRUE),
#     sd = sd(logENORM, na.rm = TRUE)
#   ) %>% arrange(mean,.by_group = TRUE)
# print(df.tab.enorm,n=100)
# wilcox.tab.enorm = pairwise.wilcox.test(df$logENORM, df$modelname,p.adjust.method = "BH")

# df.tab.rmse$modelname = gsub("\r?\n|\r", ",", df.tab.rmse$modelname)
# write.table(df.tab.rmse,file=file.path(outputdir,"summary_rmse.tsv"),row.names=FALSE,sep="\t")
# df.tab.rmse$modelname = gsub("\r?\n|\r", ",", df.tab.rmse$modelname)
# rownames(wilcox.tab.rmse$p.value) = gsub("\r?\n|\r", ",", rownames(wilcox.tab.rmse$p.value))
# colnames(wilcox.tab.rmse$p.value) = gsub("\r?\n|\r", ",", colnames(wilcox.tab.rmse$p.value))
# write.table(as.data.frame(wilcox.tab.rmse$p.value) %>% tibble::rownames_to_column(" "),file=file.path(outputdir,"wilcox_rmse.tsv"),sep="\t",row.names=FALSE)

# df.tab.enorm$modelname = gsub("\r?\n|\r", ",", df.tab.enorm$modelname)
# write.table(df.tab.enorm,file=file.path(outputdir,"summary_enorm.tsv"),row.names=FALSE,sep="\t")
# df.tab.enorm$modelname = gsub("\r?\n|\r", ",", df.tab.enorm$modelname)
# rownames(wilcox.tab.enorm$p.value) = gsub("\r?\n|\r", ",", rownames(wilcox.tab.enorm$p.value))
# colnames(wilcox.tab.enorm$p.value) = gsub("\r?\n|\r", ",", colnames(wilcox.tab.enorm$p.value))
# write.table(as.data.frame(wilcox.tab.enorm$p.value) %>% tibble::rownames_to_column(" "),file=file.path(outputdir,"wilcox_enorm.tsv"),sep="\t",row.names=FALSE)


##################################################
#
# plot mean err 1e-2,1e-5,1e-5
#
##################################################

df = plyr::rbind.fill(df_sim,df_downsamp)
df =  filter(df,(sA==1e-2 & sB==1e-5 & sG==1e-5) | (sA==1 & sB==1e-6 & sG==1e-4) | (is.na(sB) & is.na(sG)))
df$mean_dropout = round(df$mean_dropout,digits=3)

nmodels = length(unique(df$modelname))
palette1 = scales::hue_pal()(nmodels)
names(palette1) = sort(unique(df$modelname))


saveRDS(palette1,file.path("slurm","fig_panel_scripts","palette1.RDS"))

outputdir = file.path("slurm","fig_panel_scripts","fig02n50.1e-5.1e-5")
dir.create(outputdir,recursive=TRUE)
melt_tmp = reshape2::melt(df,measure.vars=c("logRMSE"))

p = ggplot(melt_tmp,
        aes(
          x=as.factor(modelname), 
          y=value,
          fill=as.factor(imputemodel),
          color=as.factor(imputemodel),
          alpha=0.5)) + 
      geom_boxplot(
        aes(
          x=as.factor(modelname),
          y=value,
          fill=as.factor(imputemodel))) +
      # facet_grid(rows=vars(mean_dropout),cols=vars(modelfam),scale="free") 
      facet_nested(strategy + mean_dropout ~ modelfam,scale="free") 
      # stat_summary(fun=mean, geom="point", shape=18, size=1,
      #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
      #                 position = position_dodge(0.75),stroke=1.5)+
      # stat_summary(fun=mean, geom="text", size=1,vjust=-0.5,
      #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
      #                 guides(fill = FALSE)
set_palette(p, palette1)
p = p +
    theme_bw() +
    rotate_x_text(30) +
    labs(color="Model\nFamily",y="Mean Error",alpha="del",fill="del") +
    theme(axis.text=element_text(size=5,face="bold"),
      axis.title.y=element_text(size=13,face="bold"),
      axis.title.x=element_blank(),
      # legend.position = "none",
      panel.border=element_rect(colour="black",size=0.5,fill=NA),
        strip.text.x = element_text(size = 11, face = "bold", angle = 0),
        strip.text.y = element_text(size = 10, face = "bold", angle = 270),      
        legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=13,face="bold"), #change legend title font size
      legend.text = element_text(size=8,face="bold"),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing=unit(0,"lines"))

ggsave(plot=p,file=file.path(outputdir,"logRMSE.pdf"),width=6.5,height=6)

################################################
#
# stats mean err 1e-2,1e-5,1e-5
#
################################################
df$strategy = as.factor(df$strategy)
df$mean_dropout = as.factor(df$mean_dropout)

df.tab.rmse = group_by(df, modelname) %>%
  summarise(
    count = n(),
    mean = mean(logRMSE, na.rm = TRUE),
    sd = sd(logRMSE, na.rm = TRUE)
  ) %>% arrange(mean,.by_group = TRUE)
xtab = xtable(df.tab.rmse, type = "latex",digits=-3) 
caption(xtab) = "Stats for benchmark mean error (RMSE)"
label(xtab) = "tab:rmsesummary"
print(xtab, file = file.path(outputdir,"rmse_stats.tex"))

wilcox.tab.rmse = pairwise.wilcox.test(df$logRMSE, df$modelname,p.adjust.method = "BH")
xtab = xtable(as.data.frame(wilcox.tab.rmse$p.value), type = "latex",digits=-3) 
caption(xtab) = "Adjusted p-values for paired Wilcox tests on benchmark mean error (RMSE)"
label(xtab) = "tab:rmsewilcox"
print(xtab, file = file.path(outputdir,"rmse_wilcox.tex"),rotate.colnames=TRUE)

# individual plots
groupdf = df %>% group_by(strategy,mean_dropout) %>% group_split()
for(i in 1:length(groupdf)){
  strategy = unique(groupdf[[i]]$strategy)
  mean_dropout = unique(groupdf[[i]]$mean_dropout)
  df.tab.rmse = group_by(groupdf[[i]], strategy, mean_dropout,modelname) %>%
  summarise(
    count = n(),
    mean = mean(logRMSE, na.rm = TRUE),
    sd = sd(logRMSE, na.rm = TRUE)
  ) %>% arrange(mean,.by_group = TRUE)

  xtab = xtable(df.tab.rmse, type = "latex",digits=-3)
  caption(xtab) = paste0("Stats for benchmark mean error (RMSE) - Strategy:",strategy,", Mean Dropout:",mean_dropout)
  label(xtab) = paste0("tab:rmsesummary_",strategy,"_",mean_dropout)
  print(xtab, file = file.path(outputdir,paste0("rmse_stats_",strategy,"_",mean_dropout,".tex")))

  wilcox.tab.rmse = pairwise.wilcox.test(groupdf[[i]]$logRMSE, groupdf[[i]]$modelname,p.adjust.method = "BH")
  xtab = xtable(as.data.frame(wilcox.tab.rmse$p.value), type = "latex",digits=-3) 
  caption(xtab) = paste0("Adjusted p-values for paired Wilcox tests on benchmark mean error (RMSE) - Strategy:",strategy,", Mean Dropout:",mean_dropout)
  label(xtab) = paste0("tab:rmsewilcox_",strategy,"_",mean_dropout)
  print(xtab, file = file.path(outputdir,paste0("rmse_wilcox_",strategy,"_",mean_dropout,".tex")),rotate.colnames=TRUE)
}

groupdf = df %>% group_by(strategy) %>% group_split()
for(i in 1:length(groupdf)){
  strategy = unique(groupdf[[i]]$strategy)
  mean_dropout = unique(groupdf[[i]]$mean_dropout)
  df.tab.rmse = group_by(groupdf[[i]], strategy, modelname) %>%
  summarise(
    count = n(),
    mean = mean(logRMSE, na.rm = TRUE),
    sd = sd(logRMSE, na.rm = TRUE)
  ) %>% arrange(mean,.by_group = TRUE)

  xtab = xtable(df.tab.rmse, type = "latex",digits=-3) 
  caption(xtab) = paste0("Stats for benchmark mean error (RMSE) - Strategy:",strategy)
  label(xtab) = paste0("tab:rmsesummary_",strategy)
  print(xtab, file = file.path(outputdir,paste0("rmse_stats_",strategy,".tex")))

  wilcox.tab.rmse = pairwise.wilcox.test(groupdf[[i]]$logRMSE, groupdf[[i]]$modelname,p.adjust.method = "BH")
  xtab = xtable(as.data.frame(wilcox.tab.rmse$p.value), type = "latex",digits=-3) 
  caption(xtab) = paste0("Adjusted p-values for paired Wilcox tests on benchmark mean error (RMSE) - Strategy:",strategy)
  label(xtab) = paste0("tab:rmsewilcox_",strategy)
  print(xtab, file = file.path(outputdir,paste0("rmse_wilcox_",strategy,".tex")),rotate.colnames=TRUE)
}

################################################
#
# plot max err 1e-2,1e-5,1e-5
#
################################################
df = plyr::rbind.fill(df_sim,df_downsamp)
df =  filter(df,(sA==1e-2 & sB==1e-5 & sG==1e-5) | (sA==1 & sB==1e-6 & sG==1e-4) | (is.na(sB) & is.na(sG)))
df$mean_dropout = round(df$mean_dropout,digits=3)

nmodels = length(unique(df$modelname))
palette1 = scales::hue_pal()(nmodels)
names(palette1) = sort(unique(df$modelname))

saveRDS(palette1,file.path("slurm","fig_panel_scripts","palette1.RDS"))

outputdir = file.path("slurm","fig_panel_scripts","fig03n50.1e-5.1e-5")
dir.create(outputdir,recursive=TRUE)
melt_tmp = reshape2::melt(df,measure.vars=c("logENORM"))

p = ggplot(melt_tmp,
        aes(
          x=as.factor(modelname), 
          y=value,
          fill=as.factor(imputemodel),
          color=as.factor(imputemodel),
          alpha=0.5)) + 
      geom_boxplot(
        aes(
          x=as.factor(modelname),
          y=value,
          fill=as.factor(imputemodel))) +
      # facet_grid(rows=vars(mean_dropout),cols=vars(modelfam),scale="free") 
      facet_nested(strategy + mean_dropout ~ modelfam,scale="free") 
      # stat_summary(fun=mean, geom="point", shape=18, size=1,
      #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
      #                 position = position_dodge(0.75),stroke=1.5)+
      # stat_summary(fun=mean, geom="text", size=1,vjust=-0.5,
      #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
      #                 guides(fill = FALSE)
set_palette(p, palette1)
p = p +
    theme_bw() +
    rotate_x_text(30) +
    labs(color="Model\nFamily",y="Max Error",alpha="del",fill="del") +
    theme(axis.text=element_text(size=5,face="bold"),
      axis.title.y=element_text(size=13,face="bold"),
      axis.title.x=element_blank(),
      # legend.position = "none",
      panel.border=element_rect(colour="black",size=0.5,fill=NA),
        strip.text.x = element_text(size = 11, face = "bold", angle = 0),
        strip.text.y = element_text(size = 10, face = "bold", angle = 270),      
      legend.key.size = unit(5, 'mm'), #change legend key size
      legend.key.height = unit(4, 'mm'), #change legend key height
      legend.key.width = unit(4, 'mm'), #change legend key width
      legend.title = element_text(size=13,face="bold"), #change legend title font size
      legend.text = element_text(size=8,face="bold"),
      legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing=unit(0,"lines"))

ggsave(plot=p,file=file.path(outputdir,"logENORM.pdf"),width=6.5,height=6)

################################################
#
# stats max err 1e-2,1e-5,1e-5
#
################################################
df$strategy = as.factor(df$strategy)
df$mean_dropout = as.factor(df$mean_dropout)

df.tab.enorm = group_by(df, modelname) %>%
  summarise(
    count = n(),
    mean = mean(logENORM, na.rm = TRUE),
    sd = sd(logENORM, na.rm = TRUE)
  ) %>% arrange(mean,.by_group = TRUE)
xtab = xtable(df.tab.enorm, type = "latex",digits=-3) 
caption(xtab) = "Stats for benchmark max error (L2 norm)"
label(xtab) = "tab:enormsummary"
print(xtab, file = file.path(outputdir,"enorm_stats.tex"))

wilcox.tab.enorm = pairwise.wilcox.test(df$logENORM, df$modelname,p.adjust.method = "BH")
xtab = xtable(as.data.frame(wilcox.tab.enorm$p.value), type = "latex",digits=-3) 
caption(xtab) = "Adjusted p-values for paired Wilcox tests on benchmark max error (L2 norm)"
label(xtab) = "tab:enormwilcox"
print(xtab, file = file.path(outputdir,"enorm_wilcox.tex"),rotate.colnames=TRUE)

# individual plots
groupdf = df %>% group_by(strategy,mean_dropout) %>% group_split()
for(i in 1:length(groupdf)){
  strategy = unique(groupdf[[i]]$strategy)
  mean_dropout = unique(groupdf[[i]]$mean_dropout)
  df.tab.enorm = group_by(groupdf[[i]], strategy, modelname) %>%
  summarise(
    count = n(),
    mean = mean(logENORM, na.rm = TRUE),
    sd = sd(logENORM, na.rm = TRUE)
  ) %>% arrange(mean,.by_group = TRUE)

  xtab = xtable(df.tab.enorm, type = "latex",digits=-3) 
  caption(xtab) = paste0("Stats for benchmark max error (L2 norm) - Strategy:",strategy,", Mean Dropout:",mean_dropout)
  label(xtab) = paste0("tab:enormsummary_",strategy,"_",mean_dropout)
  print(xtab, file = file.path(outputdir,paste0("enorm_stats_",strategy,"_",mean_dropout,".tex")))

  wilcox.tab.enorm = pairwise.wilcox.test(df$logENORM, df$modelname,p.adjust.method = "BH")
  xtab = xtable(as.data.frame(wilcox.tab.enorm$p.value), type = "latex",digits=-3) 
  caption(xtab) = paste0("Adjusted p-values for paired Wilcox tests on benchmark max error (L2 norm) - Strategy:",strategy,", Mean Dropout:",mean_dropout)
  label(xtab) = paste0("tab:enormwilcox_",strategy,"_",mean_dropout)
  print(xtab, file = file.path(outputdir,paste0("enorm_wilcox_",strategy,"_",mean_dropout,".tex")),rotate.colnames=TRUE)
}

groupdf = df %>% group_by(strategy) %>% group_split()
for(i in 1:length(groupdf)){
  strategy = unique(groupdf[[i]]$strategy)
  mean_dropout = unique(groupdf[[i]]$mean_dropout)
  df.tab.enorm = group_by(groupdf[[i]], strategy, modelname) %>%
  summarise(
    count = n(),
    mean = mean(logENORM, na.rm = TRUE),
    sd = sd(logENORM, na.rm = TRUE)
  ) %>% arrange(mean,.by_group = TRUE)

  xtab = xtable(df.tab.enorm, type = "latex",digits=-3) 
  caption(xtab) = paste0("Stats for benchmark max error (L2 norm) - Strategy:",strategy)
  label(xtab) = paste0("tab:enormsummary_",strategy)
  print(xtab, file = file.path(outputdir,paste0("enorm_stats_",strategy,".tex")))

  wilcox.tab.enorm = pairwise.wilcox.test(df$logENORM, df$modelname,p.adjust.method = "BH")
  xtab = xtable(as.data.frame(wilcox.tab.enorm$p.value), type = "latex",digits=-3) 
  caption(xtab) = paste0("Adjusted p-values for paired Wilcox tests on benchmark max error (L2 norm) - Strategy:",strategy)
  label(xtab) = paste0("tab:enormwilcox_",strategy)
  print(xtab, file = file.path(outputdir,paste0("enorm_wilcox_",strategy,".tex")),rotate.colnames=TRUE)
}

# ################################
# #
# # plot all models
# #
# ################################

# df = plyr::rbind.fill(df_sim,df_downsamp)
# df$mean_dropout = round(df$mean_dropout,digits=3)
# commonmods = intersect(df_sim$modelname,df_downsamp$modelname)
# df = filter(df,modelname %in% commonmods)

# nmodels = length(unique(df$modelname))
# palette1 = scales::hue_pal()(nmodels)
# names(palette1) = sort(unique(df$modelname))


# saveRDS(palette1,file.path("slurm","fig_panel_scripts","palette1.RDS"))

# outputdir = file.path("slurm","fig_panel_scripts","fig02n50.all")
# dir.create(outputdir,recursive=TRUE)
# melt_tmp = reshape2::melt(df,measure.vars=c("logRMSE"))

# p = ggplot(melt_tmp,
#         aes(
#           x=as.factor(modelname), 
#           y=value,
#           fill=as.factor(imputemodel),
#           color=as.factor(imputemodel),
#           alpha=0.5)) + 
#       geom_boxplot(
#         aes(
#           x=as.factor(modelname),
#           y=value,
#           fill=as.factor(imputemodel))) +
#       # facet_grid(rows=vars(mean_dropout),cols=vars(modelfam),scale="free") 
#       facet_nested(strategy + mean_dropout ~ modelfam,scale="free") 
#       # stat_summary(fun=mean, geom="point", shape=18, size=1,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
#       #                 position = position_dodge(0.75),stroke=1.5)+
#       # stat_summary(fun=mean, geom="text", size=1,vjust=-0.5,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
#       #                 guides(fill = FALSE)
# set_palette(p, palette1)
# p = p +
#     theme_bw() +
#     rotate_x_text(30) +
#     labs(color="Model Family",y="Mean Error") +
#     theme(axis.text=element_text(size=5,face="bold"),
#       axis.title.y=element_text(size=10,face="bold"),
#       axis.title.x=element_blank(),
#       # legend.position = "none",
#       panel.border=element_rect(colour="black",size=0.5,fill=NA),
#         strip.text.x = element_text(size = 8, face = "bold", angle = 0),
#         strip.text.y = element_text(size = 7, face = "bold", angle = 270),      
#         legend.key.size = unit(5, 'mm'), #change legend key size
#       legend.key.height = unit(4, 'mm'), #change legend key height
#       legend.key.width = unit(4, 'mm'), #change legend key width
#       legend.title = element_text(size=7,face="bold"), #change legend title font size
#       legend.text = element_text(size=5,face="bold"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.spacing=unit(0,"lines"))

# ggsave(plot=p,file=file.path(outputdir,"logRMSE.pdf"),width=15.5,height=6)


# outputdir = file.path("slurm","fig_panel_scripts","fig03n50.all")
# dir.create(outputdir,recursive=TRUE)
# melt_tmp = reshape2::melt(df,measure.vars=c("logENORM"))

# p = ggplot(melt_tmp,
#         aes(
#           x=as.factor(modelname), 
#           y=value,
#           fill=as.factor(imputemodel),
#           color=as.factor(imputemodel),
#           alpha=0.5)) + 
#       geom_boxplot(
#         aes(
#           x=as.factor(modelname),
#           y=value,
#           fill=as.factor(imputemodel))) +
#       # facet_grid(rows=vars(mean_dropout),cols=vars(modelfam),scale="free") 
#       facet_nested(strategy + mean_dropout ~ modelfam,scale="free") 
#       # stat_summary(fun=mean, geom="point", shape=18, size=1,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
#       #                 position = position_dodge(0.75),stroke=1.5)+
#       # stat_summary(fun=mean, geom="text", size=1,vjust=-0.5,
#       #                 aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
#       #                 guides(fill = FALSE)
# set_palette(p, palette1)
# p = p +
#     theme_bw() +
#     rotate_x_text(30) +
#     labs(color="Model Family",y="Max Error") +
#     theme(axis.text=element_text(size=5,face="bold"),
#       axis.title.y=element_text(size=10,face="bold"),
#       axis.title.x=element_blank(),
#       # legend.position = "none",
#       panel.border=element_rect(colour="black",size=0.5,fill=NA),
#         strip.text.x = element_text(size = 8, face = "bold", angle = 0),
#         strip.text.y = element_text(size = 7, face = "bold", angle = 270),
#         legend.key.size = unit(5, 'mm'), #change legend key size
#       legend.key.height = unit(4, 'mm'), #change legend key height
#       legend.key.width = unit(4, 'mm'), #change legend key width
#       legend.title = element_text(size=7,face="bold"), #change legend title font size
#       legend.text = element_text(size=5,face="bold"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.spacing=unit(0,"lines"))

# ggsave(plot=p,file=file.path(outputdir,"logENORM.pdf"),width=14.5,height=6)

# ################################################
# #
# # stats
# #
# ################################################

# df.tab.rmse = group_by(df, strategy, modelname) %>%
#   summarise(
#     count = n(),
#     mean = mean(logRMSE, na.rm = TRUE),
#     sd = sd(logRMSE, na.rm = TRUE)
#   ) %>% arrange(mean,.by_group = TRUE)
# print(df.tab.rmse,n=100)
# wilcox.tab.rmse = pairwise.wilcox.test(df$logRMSE, df$modelname,p.adjust.method = "BH")

# df.tab.enorm = group_by(df, strategy, modelname) %>%
#   summarise(
#     count = n(),
#     mean = mean(logENORM, na.rm = TRUE),
#     sd = sd(logENORM, na.rm = TRUE)
#   ) %>% arrange(mean,.by_group = TRUE)
# print(df.tab.enorm,n=100)
# wilcox.tab.enorm = pairwise.wilcox.test(df$logENORM, df$modelname,p.adjust.method = "BH")

# df.tab.rmse$modelname = gsub("\r?\n|\r", ",", df.tab.rmse$modelname)
# write.table(df.tab.rmse,file=file.path(outputdir,"summary_rmse.tsv"),row.names=FALSE,sep="\t")
# df.tab.rmse$modelname = gsub("\r?\n|\r", ",", df.tab.rmse$modelname)
# rownames(wilcox.tab.rmse$p.value) = gsub("\r?\n|\r", ",", rownames(wilcox.tab.rmse$p.value))
# colnames(wilcox.tab.rmse$p.value) = gsub("\r?\n|\r", ",", colnames(wilcox.tab.rmse$p.value))
# write.table(as.data.frame(wilcox.tab.rmse$p.value) %>% tibble::rownames_to_column(" "),file=file.path(outputdir,"wilcox_rmse.tsv"),sep="\t",row.names=FALSE)

# df.tab.enorm$modelname = gsub("\r?\n|\r", ",", df.tab.enorm$modelname)
# write.table(df.tab.enorm,file=file.path(outputdir,"summary_enorm.tsv"),row.names=FALSE,sep="\t")
# df.tab.enorm$modelname = gsub("\r?\n|\r", ",", df.tab.enorm$modelname)
# rownames(wilcox.tab.enorm$p.value) = gsub("\r?\n|\r", ",", rownames(wilcox.tab.enorm$p.value))
# colnames(wilcox.tab.enorm$p.value) = gsub("\r?\n|\r", ",", colnames(wilcox.tab.enorm$p.value))
# write.table(as.data.frame(wilcox.tab.enorm$p.value) %>% tibble::rownames_to_column(" "),file=file.path(outputdir,"wilcox_enorm.tsv"),sep="\t",row.names=FALSE)