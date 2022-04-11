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
outputdir = Sys.getenv("SUMMARYCOMBINE_ITER")
# the model parameter controlling the sparsity
sparsityparam = Sys.getenv("SPARSITY_PARAM")

library(dplyr)
library(rmcorr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
library(grid)

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
df = NULL # the dataframe of all output, to be filled
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
                paramdf$durian_errnorm_delta = rep(paramdf$durian_errnorm_final[1] - paramdf$durian_errnorm_init[1],dim(paramdf)[1])
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
                dfnew$durian_rmse_diffRatioInit = abs((dfnew$durian_rmse - dfnew$durian_rmse_final[1])/dfnew$durian_rmse_delta[1])
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

                df = rbind(df,dfnew)
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

df = df %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout = mean(dropout_rate))
df$mean_dropout = round(df$mean_dropout,digits=3)
df$mean_dropout = as.factor(df$mean_dropout)
df$deconv_method = as.factor(df$deconv_method)
df$durian_logrmse = log(df$durian_rmse)
df$durian_logerrnorm = log(df$errnorm)
write.csv(df,file=file.path(plotdir,"df.csv"))
df = as.data.frame(df)
df.orig = df

# df[,sparsityparam] = paste0(sparsityparam,"=",df[,sparsityparam])
# df[,sparsityparam] = as.factor(df[,sparsityparam])
df[,"mean_dropout"] = paste0("dropout=",df[,"mean_dropout"])
df[,"deconv_method"] = paste0("D-step=",df[,"deconv_method"])

paletteLvs = unique(as.vector(as.matrix(df[,c("mean_dropout","deconv_method")])))
df[,sparsityparam] = factor(df[,sparsityparam],levels=paletteLvs)
df[,"deconv_method"] = factor(df[,"deconv_method"],levels=paletteLvs)
if(length(grep("lambda",names(df)))>0){
  myPal = brewer.pal(n = length(paletteLvs), name = "Set1")
}else{
  myPal = brewer.pal(n = length(paletteLvs), name = "Dark2")
}
names(myPal) = paletteLvs
# myPal = c("purple","orange","yellow","brown")
# names(myPal) = paletteLvs

excludemodels = c()

yvars = c("Gene","Cell","logGene","logCell","RMSE","logRMSE","ENORM", "logENORM","MeanGene","MeanCell","logMeanGene","logMeanCell")

nmodels = length(unique(df$modelname))
ncompare = choose(nmodels,2)

m2durian = c("durian_logrmse","durian_logerrnorm","durian_rmse_diffRatioInit","durian_errnorm_diffRatioInit","durian_rmse_ratioInit","durian_errnorm_ratioInit","durian_mean_cellcor")
m2deconv = c("deconv_cor_mean_celltype","deconv_rmse_diffRatioInit","deconv_rmse_ratioInit")

plotdir1 = file.path(plotdir,"nolineMulti")
dir.create(plotdir1)
plotdir2 = file.path(plotdir,"nolineMultiFlip")
dir.create(plotdir2)

for(j in 1:length(m2durian)){
  for(jj in 1:length(m2deconv)){
    p1 = filter(df,iter>0,!is.na(deconv_rmse)) %>%
        ggplot(aes(x = iter, y=scrabbleLoss,z = get(m2durian[j]), color=deconv_method)) +
        stat_summary_hex(fun=mean,bins=5,aes(color=deconv_method),size=2,name=deconv_method) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=2,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Greens",guide = guide_colourbar(title="Impute Err")) +
        scale_color_manual(values = myPal) +
        facet_grid(cols=vars(deconv_method))+
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Grouping") +
        theme(
        axis.title=element_text(size = 25,face="bold"),,
        panel.border=element_rect(colour="black",size=1,fill=NA),
        axis.text = element_text(size = 15,face="bold"),
        strip.text.x = element_text(size = 20, face = "bold", angle = 0),
        legend.title = element_text(face="bold",size=18), #change legend title font size
        legend.text = element_text(face="bold",size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
    p2 = filter(df,iter>0,!is.na(deconv_rmse)) %>%
        ggplot(aes(x = iter, y=scrabbleLoss,z = get(m2deconv[jj]), color=deconv_method)) +
        stat_summary_hex(fun=mean,bins=5,aes(color=deconv_method),size=2,name=deconv_method) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=2,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Blues",guide = guide_colourbar(title="Deconv Err")) +
        scale_color_manual(values = myPal) +
        facet_grid(cols=vars(deconv_method))+
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Grouping") +
        theme(
        axis.title=element_text(size = 25,face="bold"),
        panel.border=element_rect(colour="black",size=1,fill=NA),
        axis.text = element_text(size = 15,face="bold"),
        strip.text.x = element_text(size = 20, face = "bold", angle = 0),
        legend.title = element_text(face="bold",size=18), #change legend title font size
        legend.text = element_text(face="bold",size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
    pgrid1 = grid.arrange(p1, p2, nrow = 2)
    ggsave(plot = pgrid1,file.path(plotdir2,paste0("du_",m2durian[j],",dc_",m2deconv[jj],",facet.pdf")),width = 8, height = 6.7)


    p3 = filter(df,iter>0,!is.na(deconv_rmse)) %>%
        ggplot(aes(x = iter, y=scrabbleLoss,z = get(m2durian[j]), color=mean_dropout)) +
        stat_summary_hex(fun=mean,bins=5,aes(color=mean_dropout),size=2,name=sparsityparam) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=2,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Greens",guide = guide_colourbar(title="Impute Err")) +
        scale_color_manual(values = myPal) +
        facet_grid(~mean_dropout) +
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Grouping") +
        theme(
        axis.title=element_text(size = 25,face="bold"),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 15,face="bold"),
        panel.border=element_rect(colour="black",size=1,fill=NA),
        strip.text.x = element_text(size = 20, face = "bold", angle = 0),
        legend.title = element_text(face="bold",size=18), #change legend title font size
        legend.text = element_text(face="bold",size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

    p4 = filter(df,iter>0,!is.na(deconv_rmse)) %>%
        ggplot(aes(x = iter, y=scrabbleLoss,z = get(m2deconv[jj]), color=mean_dropout)) +
        stat_summary_hex(fun=mean,bins=5,aes(color=mean_dropout),size=2,name=sparsityparam) +
        geom_smooth(aes(x=iter,y=scrabbleLoss),method=loess,linetype="dashed",se=FALSE,size=2,color="grey",alpha=0.5) +
        scale_fill_distiller(palette="Blues",guide = guide_colourbar(title="Deconv Err")) +
        scale_color_manual(values = myPal) +
        facet_grid(~mean_dropout) +
        theme_bw() +
        xlab("DURIAN Iterations") +
        ylab("Imputation Objective") +
        labs(color = "Grouping") +
        theme(
        axis.title=element_text(size = 25,face="bold"),
        axis.title.y=element_blank(),
        panel.border=element_rect(colour="black",size=1,fill=NA),
        axis.text = element_text(size = 15,face="bold"),
        strip.text.x = element_text(size = 20, face = "bold", angle = 0),
        legend.title = element_text(face="bold",size=18), #change legend title font size
        legend.text = element_text(face="bold",size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
    pgrid2 = grid.arrange(p3, p4, nrow = 2)
    ggsave(plot = pgrid2,file.path(plotdir2,paste0("du_",m2durian[j],",dc_",m2deconv[jj],",facet.pdf")),width = 8, height = 6.7)

    # p = grid.arrange(p1,p3,p2,p4, nrow = 2,
    # left = textGrob("Imputation Objective", hjust = 0.5,gp=gpar(fontsize=20,fontface="bold"),rot=90),
    # bottom = textGrob("DURIAN Iterations", vjust = 0.5,gp=gpar(fontsize=20,fontface="bold")))
    p = grid.arrange(p1,p3,p2,p4, nrow = 2)

    ggsave(plot = p,file.path(plotdir,paste0("du_",m2durian[j],",dc_",m2deconv[jj],",grid.pdf")),width = 18, height = 8.2)
  }
}