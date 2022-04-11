# saved benchmark directory
backupdir = Sys.getenv("PBULKDIR")
# output directory
outputdir = Sys.getenv("SUMMARYFINAL")
# the model parameter controlling the sparsity
sparsityparam = Sys.getenv("SPARSITY_PARAM")

library(dplyr)
library(rmcorr)
library(stringr)
library(ggplot2)
library(ggpubr)

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

bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
paramvec.orig = str_split(bdname,pattern=",")[[1]]

dirnames = list.files(backupdir)
df = NULL # the dataframe of all output, to be filled
idx = 1
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
  for(j in 1:length(modeldirs)){
    if(length(grep("imputemodel",modeldirs[j]))>0){
      print(paste0("found impute folder ",modeldirs[j]))
      logfn = file.path(backupdir,dirnames[i],modeldirs[j],"imputation_loss.csv")
      if(file.exists(logfn)){
        print(paste0("found impute model ",logfn))
        runlog = read.csv(logfn,row.names=1)
        runlog = cbind(runlog,params)
        df = rbind(df,runlog)
      }
    }
  }
}
df$logRMSE = log(df$RMSE)

write.csv(df,file=file.path(plotdir,"df.csv"))

excludemodels = c()

yvars = c("Gene","Cell","RMSE","logRMSE","ENORM","MeanGene","MeanCell")
sparsitygroups = unique(df[,eval(sparsityparam)])

sparsitygroup_list = df %>%
            filter(!modelname %in% excludemodels) %>%
            group_by(get(sparsityparam)) %>%
            group_split()

plotpvals <- function(plotdir=NULL,df=NULL,yvar=NULL,prefix=NULL,labscale=1.0,labplus=0.1,pvals=FALSE,plotwidth=11,plotheight=11){
  nmodels = length(unique(df$modelname))
  ncompare = choose(nmodels,2)

  palette1 = scales::hue_pal()(nmodels)
  names(palette1) = sort(unique(df$modelname))
  
  p = ggplot(df, aes(x=modelname, y=get(yvar),color=modelname)) + 
        geom_boxplot() +
        stat_summary(fun.y=mean, geom="point", shape=18, size=7) +
        rotate_x_text(45) +
        theme(axis.text=element_text(size=25,face="bold"),
            legend.position = "none",
            axis.title=element_blank())
  set_palette(p, palette1)

  ymin = min(df[,eval(yvar)])
  ymax = max(df[,eval(yvar)])
  yrange = ymax-ymin
  ylabmin = ymax + labplus*yrange
  ylabmax = ylabmin + labscale*yrange
  plabpos = seq(ylabmin,ylabmax,length=ncompare)

  if(pvals){
    stat.test = compare_means(as.formula(paste(yvar, "~ modelname")), data = df,method = "t.test")
    stat.test <- stat.test %>%
    mutate(y.position = plabpos)
    p = p + 
      stat_pvalue_manual(stat.test, label = "p = {p.adj}",label.size=2,bracket.size = 0.1,tip.length=0.005)
    ggsave(plot=p,file=file.path(plotdir,paste0(prefix,yvar,"_labp.pdf")),width=plotwidth,height=plotheight)
    write.csv(stat.test,file=file.path(plotdir,paste0(prefix,yvar,"_stat.csv")))
  }else{
    # browser()
    ggsave(plot=p,file=file.path(plotdir,paste0(prefix,yvar,"_labp.pdf")),width=plotwidth,height=plotheight)
  }
}

for(i in 1:length(sparsitygroup_list)){
  sparsityval = sparsitygroups[i]
  for(j in 1:length(yvars)){
    plotpvals(plotdir=plotdir,df=sparsitygroup_list[[i]],yvar=yvars[j],prefix=paste0("pv.",sparsityparam,".",sparsityval,"_"),pvals=TRUE)
    plotpvals(plotdir=plotdir,df=sparsitygroup_list[[i]],yvar=yvars[j],prefix=paste0(sparsityparam,".",sparsityval,"_"))
  }
}

summarystats = df %>% 
     filter(!modelname %in% excludemodels) %>%
     group_by(eval(sparsityparam),modelname) %>% 
     mutate(
       mean.gene = mean(Gene),
       mean.cell = mean(Cell),
       mean.meangene = mean(MeanGene),
       mean.meancell = mean(MeanCell),
       mean.rmse = mean(RMSE),
       mean.logrmse = mean(logRMSE),
       mean.enorm = mean(ENORM)) %>% 
    ungroup() %>%
    group_by(eval(sparsityparam),modelname) %>% 
    distinct(modelname,mean.gene,mean.cell,mean.meangene,mean.meancell,mean.rmse,mean.logrmse,mean.enorm)
     
write.csv(summarystats,file=file.path(plotdir,"summary.csv"))
