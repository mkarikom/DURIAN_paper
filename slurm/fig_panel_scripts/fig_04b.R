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

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(umap)

# myclust = parallel::makeCluster(14,type="FORK")
# doParallel::registerDoParallel(cl = myclust)

getimputeratio <- function(obs,orig){
  orig.zeros = which(orig == 0)
  orig.nonzero = which(orig > 0)
  obs.nonzero = which(obs > 0)
  obs.zimputed = intersect(orig.zeros,obs.nonzero)
  obs.nzimputed = intersect(orig.nonzero,obs.nonzero)

  zrate = sum(as.matrix(obs)[obs.zimputed]) / sum(as.matrix(obs))
  nzrate = sum(as.matrix(obs)[obs.nzimputed]) / sum(as.matrix(obs))
  list(zrate=zrate,nzrate=nzrate)
}

excludemodels = c()
quicksims = c(1,42)
outputdir = file.path("slurm","fig_panel_scripts","fig04","scatter")
dir.create(outputdir,recursive=TRUE)

backupdirs = c(
  "slurm/He/output.clusterMetrics.free.OuterMetrics1K05",
  "slurm/Baron/output.clusterMetrics.free.OuterMetrics1K05")
metafiles = c(
  "slurm/Baron/durian_data/BaronSC.H_pDataC.csv",
  "slurm/Baron/durian_data/BaronSC.DM_pDataC.csv",
  "slurm/He/durian_data/HeNL_pDataC.csv",
  "slurm/He/durian_data/HeLS_pDataC.csv")


allmetadata = NULL
for(mf in metafiles){
  meta = read.csv(mf,row.names=1)
  row.names(meta) = meta$cellID
  allmetadata = rbind(allmetadata,meta)
}

# tp<-foreach(i=1:nrow(df.orig),  .inorder = FALSE, .export = c(),
#           .packages = c())%dopar%{
# }

df.orig = NULL # the dataframe of all output, to be filled
df.orig.final = NULL
df.orig.converged = NULL
cstats_count = 0
cstats_list = list()
for(bb in 1:length(backupdirs)){
  print("reading backupdir ")
  print(backupdirs[bb])
  dirnamescombine = list.files(backupdirs[bb],full.names=TRUE)
  for(ii in 1:length(dirnamescombine)){
    print("reading dirname ")
    print(dirnamescombine[ii])
    if(length(grep("pref_",dirnamescombine[ii]))>0){
      backupdir = dirnamescombine[ii]
      bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
      modeldirs = list.files(file.path(backupdir,"output_fit"))
      dropindex = grep("dropout",modeldirs)
      dropparams = sapply(strsplit(file.path(backupdir,"output_fit",modeldirs[dropindex]),"/"),tail,n=1)
      dropfn = file.path(backupdir,"output_fit",modeldirs[dropindex],"imputed_C.csv")
      dropfnmeta = file.path(backupdir,"output_fit",modeldirs[dropindex],"pDataC.csv")
      maxsimrep = max(unlist(lapply(lapply(lapply(lapply(strsplit(dropfn,"/"),tail,n=2),"[",1),get_model_params),"[","simrep")))
      
      params1b = get_model_params(bdname)
      prefix = strsplit(params1b$pref,";")[[1]]
      for(j in 1:length(modeldirs)){
        print("reading modeldir ")
        print(modeldirs[j])
        params1 = get_model_params(modeldirs[j])
        params = bind_rows(c(params1b,params1))
        if(all(params$simrep %in% quicksims)){
          if(length(grep("imputemodel",modeldirs[j]))>0){
            # print(paste0("found impute folder ",file.path(backupdir,"output_fit",modeldirs[j])))
            logfn = file.path(backupdir,"output_fit",modeldirs[j],paste0(params$imputemodel,"_logdf.csv"))
            imputedfn = file.path(backupdir,"output_fit",modeldirs[j],"imputed_C.csv")
            metafn = file.path(backupdir,"output_fit",modeldirs[j],"pDataC.csv")
            # print(paste0("looking for ",params$simrep, ", maxsimrep=",maxsimrep))
            if(file.exists(imputedfn) && file.exists(logfn)){
              cstats_count = cstats_count + 1
              runlog = read.csv(logfn,row.names=1)
              params2 = as.data.frame(get_model_params(modeldirs[j]))
              dataparams = cbind(params,params2)
              runlog = read.csv(logfn,row.names=1)
              runlog = cbind(runlog,dataparams)
              imputedC = read.csv(imputedfn,row.names=1)
              metadata = allmetadata[colnames(imputedC),]
              umap_res = umap(t(imputedC))
              df_umap = cbind(umap_res$layout,metadata)
              colnames(df_umap) = c("U1","U2",colnames(metadata))


              runlog = cbind(tail(runlog,n=1),df_umap)

              runlog$sparsity = getsparsity(imputedC)
              runlog$ngene = nrow(imputedC)
              runlog$ncell = ncol(imputedC)
              runlog$cstats_count = cstats_count
              if(runlog$pref=="HeSC.LS.cpm1K05;SuarezLS.cpm"){
                runlog$tissuetype = "Skin"
                runlog$disease.status = "Eczema"
                runlog$labeling = "Published"
              }else if(runlog$pref=="HeSC.LS.sense01.cpm1K05;SuarezLS.sense0.1.cpm"){
                runlog$tissuetype = "Skin"
                runlog$disease.status = "Eczema"
                runlog$labeling = "Re-Labeled"
              }else if(runlog$pref=="HeSC.NL.cpm1K05;SuarezNL.cpm"){
                runlog$tissuetype = "Skin"
                runlog$disease.status = "Healthy"
                runlog$labeling = "Published"
              }else if(runlog$pref=="HeSC.NL.sense01.cpm1K05;SuarezNL.sense0.1.cpm"){
                runlog$tissuetype = "Skin"
                runlog$disease.status = "Healthy"
                runlog$labeling = "Re-Labeled"
              }else if(runlog$pref=="BaronSC.DM.isletVST1K05;SegerstolpeBulk.DM.cpm"){
                runlog$tissuetype = "Pancreas"
                runlog$disease.status = "Diabetes"
                runlog$labeling = "Published"
              }else if(runlog$pref=="BaronSC.H.isletVST1K05;SegerstolpeBulk.H.cpm"){
                runlog$tissuetype = "Pancreas"
                runlog$disease.status = "Healthy"
                runlog$labeling = "Published"
              }else{
                runlog$tissuetype = NA
                runlog$disease.status = NA
                runlog$labeling = NA
              }

              if(length(grep("DURIAN",runlog$imputemodel))>0){
                runlog$modelfam=c("DURIAN")
                runlog$modelname = paste0(paste0(runlog$imputemodel,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
              }else if(length(grep("mtSCRABBLE",runlog$imputemodel))>0){
                runlog$modelfam=c("Control")
                runlog$modelname = paste0(paste0(runlog$imputemodel,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
              }else if(length(grep("SCRABBLE",runlog$imputemodel))>0){
                runlog$modelfam=c("Existing Methods")
                runlog$modelname = paste0(paste0(runlog$imputemodel,"\n"),paste(runlog$sA,runlog$sB,runlog$sG,sep=","))
              }else if(length(grep("dropout",runlog$imputemodel))>0){
                runlog$modelfam=c("Control")
                runlog$modelname=runlog$imputemodel
              }else{
                runlog$modelfam=c("Existing Methods")
                runlog$modelname=runlog$imputemodel
              }
              print(paste0("dim runlog=",dim(runlog)))
              stretchrl = NULL
              for(rowi in 1:nrow(df_umap)){
                stretchrl = rbind(stretchrl,runlog)
              }
              df.orig = plyr::rbind.fill(df.orig,stretchrl)
            }
          }
        }else{
          print(paste0("simrep ",params$simrep," not in quicksims"))
        }
      }
    }
  }
}

print("finished compiling results")
write.csv(df.orig,file.path(outputdir,"df.orig.csv"))


df.temp = df.orig  %>% filter(labeling=="Published",sA %in% c(1,NA)) %>% distinct(.keep_all = TRUE) 
df.temp$modelname = as.factor(df.temp$modelname)
df.temp$modelname = factor(df.temp$modelname,levels=c("dropout","DURIAN.MuSiC\n1,1e-06,1e-04","DrImpute","SCRABBLE\n1,1e-06,1e-04"))
write.csv(df.temp,gzfile(file.path(outputdir,"df.tmp.csv.gz")))
p = ggplot(df.temp,aes(x=U1, y=U2,color=cellType)) + 
      geom_point(size=0.1)+
      facet_nested(tissuetype+disease.status~modelname,scales="free")
p = p +
     theme_bw() +
     rotate_x_text(90) +
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
         scale_alpha(range = c(0.01,.8))+
  guides(color = guide_legend(override.aes = list(size=3)))+
labs(color="Cell Type") 

ggsave(plot=p,file=file.path(outputdir,"facet_scatter.pdf"),width=6.5,height=3.0)
