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
library(rmcorr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(scales)
library(umap)
library(viridis)
library(jcolors)

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
quicksims = 1:100
outputdir = file.path("slurm","fig_panel_scripts","fig04","alldata")
dir.create(outputdir,recursive=TRUE)
cmetrics = c("within.cluster.ss","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch","widestgap","sindex","sparsity")

backupdirscombine = c(
  "slurm/Baron/output.clusterMetrics.free.OuterMetricsResampleShort",
  "slurm/Baron/output.clusterMetrics.free.OuterMetricsResample",
  "slurm/Gupta/output.clusterMetrics.free.OuterMetricsResample",
  "slurm/Gupta/output.clusterMetrics.free.OuterMetricsResampleShort",
  "slurm/He/output.clusterMetrics.free.OuterMetricsResample",
  "slurm/He/output.clusterMetrics.free.OuterMetricsResampleShort")

sourcepaths = c("slurm/Baron/durian_data","slurm/Gupta/durian_data","slurm/He/durian_data")

df.orig = NULL # the dataframe of all output, to be filled
df.orig.final = NULL
df.orig.converged = NULL
cstats_count = 0
cstats_list = list()
for(bb in 1:length(backupdirscombine)){
  dirnamescombine = list.files(backupdirscombine[bb],full.names=TRUE)
  for(ii in 1:length(dirnamescombine)){
    if(length(grep("pref_",dirnamescombine[ii]))>0){
      backupdir = dirnamescombine[ii]
      bdname = strsplit(backupdir,"/")[[1]][length(strsplit(backupdir,"/")[[1]])]
      modeldirs = list.files(file.path(backupdir,"output_fit"))
      dropindex = grep("dropout",modeldirs)
      dropparams = sapply(strsplit(file.path(backupdir,"output_fit",modeldirs[dropindex]),"/"),tail,n=1)
      dropfn = file.path(backupdir,"output_fit",modeldirs[dropindex],"imputed_C.csv")
      dropfnmeta = file.path(backupdir,"output_fit",modeldirs[dropindex],"pDataC.csv")
      maxsimrep = max(unlist(lapply(lapply(lapply(lapply(strsplit(dropfn,"/"),tail,n=2),"[",1),get_model_params),"[","simrep")))
      drop_list = vector("list",maxsimrep)
      umap_dropout_list = vector("list",maxsimrep)
      dropout_meta_list = vector("list",maxsimrep)
      cstats_dropout_list = vector("list",maxsimrep)
      sparsity_init_list = vector("list",maxsimrep)
      for(dfni in 1:length(dropfn)){
        print(paste0("found dropout folder ",dfni))
        dparam = get_model_params(dropparams[dfni])
        if(file.exists(dropfn[dfni]) & file.exists(dropfnmeta[dfni])){
          meta = read.csv(dropfnmeta[dfni],row.names=1)
          dropout_meta_list[[dparam$simrep]] = meta
          drop_list[[dparam$simrep]] = read.csv(dropfn[dfni],row.names=1)
          umap_dropout_list[[dparam$simrep]] = umap(t(drop_list[[dparam$simrep]]))
          cstats_dropout_list[[dparam$simrep]] = fpc::cluster.stats(dist(umap_dropout_list[[dparam$simrep]]$layout),as.integer(as.factor(meta$cellType)))
          cstats_dropout_list[[dparam$simrep]]$sparsity = getsparsity(drop_list[[dparam$simrep]])
          sparsity_init_list[[dparam$simrep]] = getsparsity(drop_list[[dparam$simrep]])
        }else{
          print(paste0("missing dropout info ",dfni))
        }
      }

      params1b = get_model_params(bdname)
      prefix = strsplit(params1b$pref,";")[[1]]
      for(j in 1:length(modeldirs)){
        params1 = get_model_params(modeldirs[j])
        params = bind_rows(c(params1b,params1))
        if(all(params$simrep %in% quicksims)){
          if(length(grep("imputemodel",modeldirs[j]))>0){
            print(paste0("found impute folder ",file.path(backupdir,"output_fit",modeldirs[j])))
            logfn = file.path(backupdir,"output_fit",modeldirs[j],paste0(params$imputemodel,"_logdf.csv"))
            imputedfn = file.path(backupdir,"output_fit",modeldirs[j],"imputed_C.csv")
            metafn = file.path(backupdir,"output_fit",modeldirs[j],"pDataC.csv")
            print(paste0("looking for ",params$simrep, ", maxsimrep=",maxsimrep))
            if(params$simrep <= maxsimrep && file.exists(imputedfn) && file.exists(logfn) && !is.null(cstats_dropout_list[[params$simrep]])){
              cstats_count = cstats_count + 1
              runlog = read.csv(logfn,row.names=1)

              imputedC = read.csv(imputedfn,row.names=1)
              metadata = dropout_meta_list[[params$simrep[1]]]


              umap_res = umap(t(imputedC))
              cstats = fpc::cluster.stats(dist(umap_res$layout),as.integer(as.factor(metadata$cellType)))
              cstats_list[[cstats_count]] = cstats
              runlog = cbind(runlog,params,t(unlist(cstats[cmetrics])),as.data.frame(getimputeratio(imputedC,drop_list[[dparam$simrep]])))
              runlog[1,c("zrate","nzrate")] = c(NA,NA)

              runlog$sparsity = getsparsity(imputedC)
              runlog$ngene = nrow(imputedC)
              runlog$ncell = ncol(imputedC)
              runlog$cstats_count = cstats_count
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

              # if(runlog$imputemodel[1] == "DURIAN.MuSiC"){
              #   runlog[1,cmetrics] = cstats_dropout_list[[runlog$simrep[1]]][cmetrics]
              #   runlog[1,c("zrate","nzrate")] = c(NA,NA)

              #   emdirs.full = list.files(file.path(backupdir,"output_fit",modeldirs[j]),full.names=TRUE)[grep("emIter_",list.files(file.path(backupdir,"output_fit",modeldirs[j])))]
              #   emdirs = list.files(file.path(backupdir,"output_fit",modeldirs[j]))[grep("emIter_",list.files(file.path(backupdir,"output_fit",modeldirs[j])))]
              #   for(emdir in 1:(length(emdirs)-1)){
              #     emiter = get_model_params(emdirs[emdir])$emIter
              #     if(emiter > 0){
              #       print(paste0("recording iteration ",emiter," stats"))
              #       imputeIter = read.csv(file.path(emdirs.full[emdir],"imputed_C.csv"),row.names=1)
              #       umapiter = umap(t(imputeIter))
              #       cstatsiter = fpc::cluster.stats(dist(umapiter$layout),as.integer(as.factor(metadata$cellType)))
              #       cstatsiter$sparsity = getsparsity(imputeIter)
              #       runlog[emiter,cmetrics] = cstatsiter[cmetrics]
              #       runlog[emiter,c("zrate","nzrate")] = getimputeratio(imputeIter,drop_list[[dparam$simrep]])[c("zrate","nzrate")]
              #     }
              #   }
              # }

              df.orig = plyr::rbind.fill(df.orig,runlog)
              rlstart = runlog[1,]
              rlstart = plyr::rbind.fill(
                reshape2::melt(cbind(rlstart,diameter = cstats_dropout_list[[runlog$simrep[1]]]$diameter[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$diameter)]),measure.vars="diameter"),
                reshape2::melt(cbind(rlstart,average.distance = cstats_dropout_list[[runlog$simrep[1]]]$average.distance[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$average.distance)]),measure.vars="average.distance"),
                reshape2::melt(cbind(rlstart,median.distance = cstats_dropout_list[[runlog$simrep[1]]]$median.distance[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$median.distance)]),measure.vars="median.distance"),
                reshape2::melt(cbind(rlstart,separation = cstats_dropout_list[[runlog$simrep[1]]]$separation[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$separation)]),measure.vars="separation"),
                reshape2::melt(cbind(rlstart,average.toother = cstats_dropout_list[[runlog$simrep[1]]]$average.toother[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$average.toother)]),measure.vars="average.toother"),
                reshape2::melt(cbind(rlstart,clus.avg.silwidths = cstats_dropout_list[[runlog$simrep[1]]]$clus.avg.silwidths[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$clus.avg.silwidths)]),measure.vars="clus.avg.silwidths"),
                reshape2::melt(cbind(rlstart,cwidegap = cstats_dropout_list[[runlog$simrep[1]]]$cwidegap[!is.na(cstats_dropout_list[[runlog$simrep[1]]]$cwidegap)]),measure.vars="cwidegap"),
                reshape2::melt(cbind(rlstart,separation.matrix = cstats_dropout_list[[runlog$simrep[1]]]$separation.matrix[upper.tri(cstats_dropout_list[[runlog$simrep[1]]]$separation.matrix)][!is.na(cstats_dropout_list[[runlog$simrep[1]]]$separation.matrix[upper.tri(cstats_dropout_list[[runlog$simrep[1]]]$separation.matrix)])]),measure.vars="separation.matrix"),
                reshape2::melt(cbind(rlstart,ave.between.matrix = cstats_dropout_list[[runlog$simrep[1]]]$ave.between.matrix[upper.tri(cstats_dropout_list[[runlog$simrep[1]]]$ave.between.matrix)][!is.na(cstats_dropout_list[[runlog$simrep[1]]]$ave.between.matrix[upper.tri(cstats_dropout_list[[runlog$simrep[1]]]$ave.between.matrix)])]),measure.vars="ave.between.matrix")
              )
              rltail = tail(runlog,n=1)
              rltail = plyr::rbind.fill(
                reshape2::melt(cbind(rltail,diameter = cstats$diameter[!is.na(cstats$diameter)]),measure.vars="diameter"),
                reshape2::melt(cbind(rltail,average.distance = cstats$average.distance[!is.na(cstats$average.distance)]),measure.vars="average.distance"),
                reshape2::melt(cbind(rltail,median.distance = cstats$median.distance[!is.na(cstats$median.distance)]),measure.vars="median.distance"),
                reshape2::melt(cbind(rltail,separation = cstats$separation[!is.na(cstats$separation)]),measure.vars="separation"),
                reshape2::melt(cbind(rltail,average.toother = cstats$average.toother[!is.na(cstats$average.toother)]),measure.vars="average.toother"),
                reshape2::melt(cbind(rltail,clus.avg.silwidths = cstats$clus.avg.silwidths[!is.na(cstats$clus.avg.silwidths)]),measure.vars="clus.avg.silwidths"),
                reshape2::melt(cbind(rltail,cwidegap = cstats$cwidegap[!is.na(cstats$cwidegap)]),measure.vars="cwidegap"),
                reshape2::melt(cbind(rltail,separation.matrix = cstats$separation.matrix[upper.tri(cstats$separation.matrix)][!is.na(cstats$separation.matrix[upper.tri(cstats$separation.matrix)])]),measure.vars="separation.matrix"),
                reshape2::melt(cbind(rltail,ave.between.matrix = cstats$ave.between.matrix[upper.tri(cstats$ave.between.matrix)][!is.na(cstats$ave.between.matrix[upper.tri(cstats$ave.between.matrix)])]),measure.vars="ave.between.matrix")
              )
              df.orig.converged = plyr::rbind.fill(df.orig.converged,rltail)
              if(runlog$imputemodel != "dropout"){
                print("calling rbind")
                runlog_ends = rbind(rlstart,rltail)
                runlog_ends[1,cmetrics] = cstats_dropout_list[[runlog$simrep[1]]][cmetrics]
                runlog_ends$sparsity[1] = sparsity_init_list[[runlog$simrep[1]]]
                df.orig.final = plyr::rbind.fill(df.orig.final,runlog_ends)
              }else{
                df.orig.final = plyr::rbind.fill(df.orig.final,runlog)
              }
            }
          }
        }else{
          print(paste0("simrep ",params$simrep," not in quicksims"))
        }
      }
    }
  }
}

ds1 = unname(sapply(sapply(df.orig$pref,strsplit,"\\."),"[[",1))
ds2 = unname(sapply(sapply(df.orig$pref,strsplit,"\\."),"[[",2))
ds3 = unname(sapply(sapply(df.orig$pref,strsplit,"\\."),"[[",3))
df.orig$dataset = paste0(ds1,".",ds2,".",ds3)

ds1 = unname(sapply(sapply(df.orig.final$pref,strsplit,"\\."),"[[",1))
ds2 = unname(sapply(sapply(df.orig.final$pref,strsplit,"\\."),"[[",2))
ds3 = unname(sapply(sapply(df.orig.final$pref,strsplit,"\\."),"[[",3))
df.orig.final$dataset = paste0(ds1,".",ds2,".",ds3)

ds1 = unname(sapply(sapply(df.orig.converged$pref,strsplit,"\\."),"[[",1))
ds2 = unname(sapply(sapply(df.orig.converged$pref,strsplit,"\\."),"[[",2))
ds3 = unname(sapply(sapply(df.orig.converged$pref,strsplit,"\\."),"[[",3))
df.orig.converged$dataset = paste0(ds1,".",ds2,".",ds3)

print("finished compiling results")
write.csv(df.orig,file.path(outputdir,"df.orig.csv"))
write.csv(df.orig.final,file.path(outputdir,"df.orig.final.csv"))
write.csv(df.orig.converged,file.path(outputdir,"df.orig.converged.csv"))

browser()


finalstats = c("diameter","average.distance","median.distance","separation","average.toother","clus.avg.silwidths","cwidegap","separation.matrix","ave.between.matrix")
cmetrics = c("within.cluster.ss","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch","widestgap","sindex","sparsity")


df.orig.converged = read.csv(file.path(outputdir,"df.orig.converged.csv"))

df.orig.converged.unmelt = df.orig.converged[!names(df.orig.converged) %in% c("value", "variable")] %>% distinct(.keep_all = TRUE)
df.orig.converged.remelt = reshape2::melt(df.orig.converged.unmelt,measure.vars=cmetrics)
renames = colnames(df.orig.converged.remelt[!names(df.orig.converged.remelt) %in% c("value", "variable")])
df.converged = rbind(df.orig.converged[c(renames,"value","variable")],df.orig.converged.remelt)%>% distinct(.keep_all = TRUE)




tissuetype = c()
status = c()
for(i in 1:nrow(df.converged)){
  if(df.converged$dataset[i]=="BaronSC.DM.isletVST;SegerstolpeBulk"){
    tissuetype[i] = "Pancreas"
    status[i] = "Diabetic"
  }else if(df.converged$dataset[i]=="BaronSC.H.isletVST;SegerstolpeBulk"){
    tissuetype[i] = "Pancreas"
    status[i] = "Healthy"
  }else if(df.converged$dataset[i]=="GuptaE13SC.VST;BiggsBulk.VST"){
    tissuetype[i] = "Skin"
    status[i] = "Embryo"
  }else if(df.converged$dataset[i]=="HeLS.sense0.1"){
    tissuetype[i] = "Skin"
    status[i] = "Eczema"
  }else if(df.converged$dataset[i]=="HeNL.sense0.1"){
    tissuetype[i] = "Skin"
    status[i] = "Healthy"
  }
}

df.converged$tissuetype = tissuetype
df.converged$status = status

# write.csv(df.converged,file=gzfile(file.path(outputdir,"df.converged.csv.gz")))
outputdir = file.path("slurm","fig_panel_scripts","fig04","alldata")
df.converged = read.csv(file=gzfile(file.path(outputdir,"df.converged.csv.gz")))

library(ggplot2)
library(ggpubr)
library(ggh4x)
library(dplyr)

stats.use = c("cwidegap","diameter","dunn","sindex","within.cluster.ss")
dir.create(file.path(outputdir,"box_allvars"))
fn = file.path(outputdir,"box_allvars","summary_cluster.pdf")
df.filt = df.converged%>%filter(variable %in% stats.use) %>% distinct()
p = ggplot(
  data=df.converged%>%filter(variable %in% stats.use),
  aes(modelname,log(value),color=imputemodel))+
  geom_boxplot()+
facet_nested(variable ~ tissuetype + status,scales="free")+
    theme_bw() +
        rotate_x_text(30) +
            labs(color="Model Family",y="log(metric)") +
    theme(axis.text=element_text(size=5,face="bold"),
      axis.title.y=element_text(size=10,face="bold"),
      axis.title.x=element_blank(),
      # legend.position = "none",
      panel.border=element_rect(colour="black",size=0.5,fill=NA),
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
      ggsave(plot=p,file=fn,width=6.5,height=6.5)