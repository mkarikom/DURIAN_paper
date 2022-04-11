# saved benchmark directory
backupdir = Sys.getenv("PBULKDIR")
# output directory
outputdir = Sys.getenv("SUMMARYITER")
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
# define the variables for grouping
# groupvars = c("scAlpha","scBeta","scGamma","scrgthresh","ldagthresh") # what variables define the groups?

# groupvars = c("scAlpha","scBeta","scGamma","scrgthresh","ldagthresh","am","asd","dm","dsd") # what variables define the groups?
groupvars = c(sparsityparam,"deconv_method") # what variables define the groups?


# define the variables for coloring of histograms
# m2vec = c("durian_cellcor","durian_genecor","durian_rmse","deconv_cor","deconv_rmse")
m2vec = c("durian_cellcor","durian_rmse","durian_mean_cellcor","errnorm","deconv_cor_celltype","deconv_cor_bulksample","deconv_rmse")

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
  paramvec = str_split(dirnames[i],pattern=",")[[1]]
  paramvec = c(paramvec,paramvec.orig)
  params = c()
  for(j in 1:length(paramvec)){
    newparam = str_split(paramvec[j],pattern="_")[[1]]
    params[eval(newparam[1])] = parseargs(newparam[2])
  }
  modeldirs = list.dirs(file.path(backupdir,dirnames[i]),recursive=FALSE,full.names=FALSE)
  for(j in 1:length(modeldirs)){
    if(length(grep("imputemodel",modeldirs[j]))>0){
      if(length(grep("DURIAN",modeldirs[j]))>0){
        deconv_method = strsplit(modeldirs[j],"\\.")[[1]][2]
        model_fn = file.path(backupdir,dirnames[i],modeldirs[j])
        logfn = file.path(model_fn,"durian_logdf.csv")
        if(file.exists(logfn)){
          print(paste0("found impute model ",logfn))
          countfound=countfound+1
          foundfolders[[countfound]] = model_fn
          logdf = read.csv(logfn,row.names=1)
          rmse_delta = logdf$durian_rmse[dim(logdf)[1]] - logdf$durian_rmse[1]
          gene_delta = logdf$durian_genecor[dim(logdf)[1]] - logdf$durian_genecor[1]
          cell_delta = logdf$durian_cellcor[dim(logdf)[1]] - logdf$durian_cellcor[1]
          deconv_rmse_delta = logdf$deconv_rmse[dim(logdf)[1]] - logdf$deconv_rmse[2]
          deconv_celltype_delta = logdf$deconv_cor_celltype[dim(logdf)[1]] - logdf$deconv_cor_celltype[2]
          deconv_bulk_delta = logdf$deconv_cor_bulksample[dim(logdf)[1]] - logdf$deconv_cor_bulksample[2]
          params = params[sort(names(params))]
          paramdf = as.data.frame(matrix(rep(params,dim(logdf)[1]),nrow=dim(logdf)[1],byrow=TRUE))
          names(paramdf) = sort(names(params))
          paramdf$durian_idx = rep(idx,dim(paramdf)[1])
          paramdf$rmse_delta = rep(rmse_delta,dim(logdf)[1])
          paramdf$gene_delta = rep(gene_delta,dim(logdf)[1])
          paramdf$cell_delta = rep(cell_delta,dim(logdf)[1])
          paramdf$deconv_method = rep(deconv_method,dim(logdf)[1])
          paramdf$deconv_rmse_delta = rep(deconv_rmse_delta,dim(logdf)[1])
          paramdf$deconv_celltype_delta = rep(deconv_celltype_delta,dim(logdf)[1])
          paramdf$deconv_bulk_delta = rep(deconv_bulk_delta,dim(logdf)[1])
          df = rbind(df,cbind(logdf,paramdf))
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

write.csv(df,file=file.path(plotdir,"df.csv"))

# vs_key = df %>% group_by(deconv_method,durian_idx) %>% group_split()
vs_key = df %>% group_by(deconv_method) %>% group_split()

getrmcor <- function(df,participant="durian_idx",meas1="iter",meas2="deconv_cor"){
    print(length(unique(df[,eval(participant)])))
    print(df[,c(participant,meas1,meas2)])
    check = df %>% group_by(eval(participant)) %>% group_split()
    # browser()
    length(unique(df[,eval(meas2)]))
    rmc = rmcorr(
      participant=participant,
      measure1=meas1,
      measure2=meas2,
      dataset=df)
    tib = tibble::tibble(rmc$r,rmc$p)
    names(tib) = c(paste0(meas2,"_r"),paste0(meas2,"_p"))
    return(list(tib=tib,rmc=rmc))
}
dfrows = list()
counter0=0
for(i in 1:length(vs_key)){
  counter0=counter0+1
  dftmp = list()
  keyvals = paste(names(vs_key[[i]][groupvars]),unique(vs_key[[i]][groupvars]),sep="_",collapse=",")
  dir.create(file.path(plotdir,keyvals))
  counter=0
  for(j in 1:length(m2vec)){
    counter = counter+1
    min_iterlength_idx = min(unlist(lapply(df %>% group_by(get("durian_idx")) %>% group_split(),FUN=function(y){
      nrow(y %>% filter(if_all(m2vec[j], ~ !is.na(.x))))
    })))
    if(min_iterlength_idx > 3){
      tmplist = getrmcor(df=as.data.frame(vs_key[i][[1]]),meas2=m2vec[j])
      tmp = tmplist[["tib"]]
      dftmp[[counter]] = tmp
      p = ggplot(vs_key[[i]], aes(x=iter, y=scrabbleLoss,z=get(m2vec[j]))) +
        stat_summary_hex(fun=mean,bins=5) +
        guides(fill=guide_colourbar(title=m2vec[j]))
      ggsave(plot = p,file.path(plotdir,keyvals,paste0(m2vec[j],"_heatmap.pdf")),width = 4, height = 3)

      pdf(file.path(plotdir,keyvals,paste0(m2vec[j],"_rmc.pdf")),width = 4, height = 3)
        plot(tmplist[["rmc"]])
      dev.off()
    }
  }
  dfrows[[counter0]] = bind_cols(dftmp)
  dfrows[[counter0]] = bind_cols(dfrows[[counter0]],unique(vs_key[i][[1]][,groupvars]))
}
trend_df = bind_rows(dfrows)

write.csv(trend_df,file.path(outputdir,paste0(paste0(strsplit(outputdir,"_")[[1]][-c(1:6)],collapse="_"),"_ancova.csv")))