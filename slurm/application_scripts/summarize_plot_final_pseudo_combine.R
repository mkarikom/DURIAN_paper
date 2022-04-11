###############################################
# delete this
###############################################
# Sys.setenv(PBULKDIRCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/output.final.baron")
# Sys.setenv(SUMMARYCOMBINE="/dfs5/bio/mkarikom/temp/DURIAN/slurm/durian_pseudobulk/realdata_scripts/test/summaryfinalcombine")
# Sys.setenv(SPARSITY_PARAM="lambda")
###############################################

# saved benchmark directory
backupdircombine = Sys.getenv("OUTPUTMASTER")
# output directory
outputdir = Sys.getenv("SUMMARYCOMBINE")
# the model parameter controlling the sparsity
sparsityparam = Sys.getenv("SPARSITY_PARAM")

library(dplyr)
library(rmcorr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggh4x)

plotdir = file.path(outputdir,"plots")

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


dirnamescombine = list.files(backupdircombine,full.names=TRUE)
df.orig = NULL # the dataframe of all output, to be filled
idx = 1
for(ii in 1:length(dirnamescombine)){
  if(length(grep("_fit",dirnamescombine[ii]))>0){
    # browser()
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
      params = bind_rows(params)
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
              runlog = read.csv(logfn,row.names=1)
              runlog = cbind(runlog,params)
              runlog$dropout_rate = droprate_backup
              if(length(grep("DURIAN",runlog$modelname))>0){
                runlog$modelfam=c("DURIAN")
              }else if(length(grep("dropout",runlog$modelname))>0 || length(grep("mtSCRABBLE",runlog$modelname))>0){
                runlog$modelfam=c("Control")
              }else{
                runlog$modelfam=c("Existing Methods")
              }
              df.orig = rbind(df.orig,runlog)
            }
          }
        }
      }
    }
  }
}
df.orig$logRMSE = log(df.orig$RMSE)
df.orig$logMeanGene = -log(df.orig$MeanGene)
df.orig$logMeanCell = -log(df.orig$MeanCell)
df.orig$logGene = -log(df.orig$Gene)
df.orig$logCell = -log(df.orig$Cell)
df.orig$logENORM = log(df.orig$ENORM)

df.orig[,sparsityparam] = as.factor(df.orig[,sparsityparam])
write.csv(df.orig,file=file.path(plotdir,"df.orig.csv"))

excludemodels = c()

yvars = c("Gene","Cell","logGene","logCell","RMSE","logRMSE","ENORM", "logENORM","MeanGene","MeanCell","logMeanGene","logMeanCell")

nmodels = length(unique(df.orig$modelname))
ncompare = choose(nmodels,2)

palette1 = scales::hue_pal()(nmodels)
names(palette1) = sort(unique(df.orig$modelname))

df.filter = df.orig %>% filter(!modelname %in% excludemodels)

df = df.filter %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))

# browser()
for(j in 1:length(yvars)){

  p = ggplot(df,
          aes(
            x=as.factor(modelname), 
            y=get(yvars[j]),
            fill=as.factor(modelname))) + 
        geom_boxplot(
          aes(
            x=as.factor(modelname),
            y=get(yvars[j]),
            fill=as.factor(modelname),
            alpha=as.factor(mean_dropout))) +
        scale_alpha_discrete(
          range = c(0.1, 0.8),
          guide = guide_legend(override.aes = list(fill = "black")),name="dropout rate") +
        facet_grid(~modelfam,scale="free") +
        stat_summary(fun=mean, geom="point", shape=18, size=8,
                        aes(group=as.factor(mean_dropout),x=as.factor(modelname), y=get(yvars[j])),
                        position = position_dodge(.75),stroke=1.5)+
        stat_summary(fun=mean, geom="text", size=8,vjust=-0.5,
                        aes(group=as.factor(mean_dropout),x=as.factor(modelname), label=round(..y..,2)),position = position_dodge(.75)) +
                        guides(fill = FALSE)  
  set_palette(p, palette1)
  p = p +
      theme_bw() +
        rotate_x_text(45) +
      theme(axis.text=element_text(size=25,face="bold"),
        # legend.position = "none",
        axis.title=element_blank(),
        panel.border=element_rect(colour="black",size=1,fill=NA),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=20,face="bold"), #change legend title font size
        legend.text = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        ) 
  ggsave(plot=p,file=file.path(plotdir,paste0(yvars[j],"_combineShort.pdf")),width=15,height=8)
  ggsave(plot=p,file=file.path(plotdir,paste0(yvars[j],"_combine.pdf")),width=15,height=13)
}