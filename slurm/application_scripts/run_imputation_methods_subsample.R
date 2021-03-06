# run_imputation_methods_clusterMetrics.R

# .libPaths(c("/tmp/mkarikom/mylibs","/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0"))
# .libPaths(c(Sys.getenv("LOCAL_R_LIBS_USER"),Sys.getenv("R_LIBS_USER")))
.libPaths(c(Sys.getenv("R_LIBS_USER")))

library(Seurat,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(MuSiC,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")
library(Biobase)
library(xbioc)
library(CellChat)
library(ComplexHeatmap)
library(reshape2, include.only = c("melt"))
library(umap,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")
library(mtSCRABBLE,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")
library(fpc,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")
library(DURIAN,lib.loc="/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0")

### arguments
datapath = Sys.getenv("DATAPATH")
sourcepath = Sys.getenv("SOURCEPATH")
prefix = strsplit(Sys.getenv("PREFIXTUPLE"),";")[[1]]
imethod = Sys.getenv("IMPUTE_METHOD")
deconv_method = Sys.getenv("DECONVMETHOD")
error_out_threshold = as.numeric(Sys.getenv("ERR_OUT_THRESH"))
error_inner_threshold = as.numeric(Sys.getenv("ERR_IN_THRESH"))
nEM = as.integer(Sys.getenv("nEM"))
durianEps=as.numeric(Sys.getenv("durianEps"))

scrgenethresh=as.numeric(Sys.getenv("SCRGENETHRESH"))
deconvgenethresh=as.numeric(Sys.getenv("DECONVGENETHRESH"))

DunIterOuter = as.integer(Sys.getenv("DunIterOuter"))
DunIterInner = as.integer(Sys.getenv("DunIterInner"))
DunSDCIters = as.integer(Sys.getenv("DunSDCIters"))
DurianAlpha=as.numeric(Sys.getenv("DurianAlpha"))
DurianBeta=as.numeric(Sys.getenv("DurianBeta"))
DurianGamma=as.numeric(Sys.getenv("DurianGamma"))
useirlba=as.logical(Sys.getenv("USEIRLBA"))

print(paste0("value of RUNOUTERSTATS=",Sys.getenv("RUNOUTERSTATS")))
if(nchar(Sys.getenv("RUNOUTERSTATS"))>1){
    runouterstats=as.logical(Sys.getenv("RUNOUTERSTATS"))
}else{
    runouterstats=FALSE
}

print(paste0("value of RUNSTABILITY=",Sys.getenv("RUNSTABILITY")))
if(nchar(Sys.getenv("RUNSTABILITY"))>1){
    runstability=as.logical(Sys.getenv("RUNSTABILITY"))
}else{
    runstability=FALSE
}

ScrnIterOuter = as.integer(Sys.getenv("ScrnIterOuter"))
ScrnIterInner = as.integer(Sys.getenv("ScrnIterInner"))
ScrnSDCIters = as.integer(Sys.getenv("ScrnSDCIters"))
ScrabbleAlpha=as.numeric(Sys.getenv("ScrabbleAlpha"))
ScrabbleBeta=as.numeric(Sys.getenv("ScrabbleBeta"))
ScrabbleGamma=as.numeric(Sys.getenv("ScrabbleGamma"))
simrep=as.integer(Sys.getenv("SIMREP"))
        
### external code
source(Sys.getenv("ETCLIB"))
source("slurm/scrabble_helper_functions/library_cluster_metrics.R")
source("slurm/scrabble_helper_functions/library_extra_plots.R")

### additional methods
g2s3_fn = Sys.getenv("G2S3LIB")
cmfimpute_fn = Sys.getenv("CMFLIB")

fn_T = file.path(sourcepath,paste0(prefix[2],"_T.csv"))
fn_pDataC = file.path(sourcepath,paste0(prefix[1],"_pDataC.csv"))
fn_C = file.path(sourcepath,paste0(prefix[1],"_C.csv"))
fn_trueC = file.path(sourcepath,paste0(prefix[1],"_trueC.csv"))
fn_trueP = file.path(sourcepath,paste0(prefix[2],"_trueP.csv"))
print(paste0("exists T=",file.exists(fn_T)))
print(paste0("exists pdatac=",file.exists(fn_pDataC)))
print(paste0("exists c=",file.exists(fn_C)))
print(paste0("exists truec=",file.exists(fn_trueC)))
print(paste0("exists truep=",file.exists(fn_trueP)))

if(file.exists(fn_trueC)){
    print("case 1: with ground truth, with prefix")
    C = read.csv(fn_C,row.names=1)
    T = read.csv(fn_T,row.names=1)
    pDataC = read.csv(fn_pDataC,row.names=1)
    trueC = read.csv(fn_trueC,row.names=1)
    trueP = read.csv(fn_trueP,row.names=1)
}else if(file.exists(file.path(sourcepath,"trueC.csv"))){
    print("case 2: with ground truth, no prefix")
    C = read.csv(file.path(sourcepath,"C.csv"),row.names=1)
    T = read.csv(file.path(sourcepath,"T.csv"),row.names=1)
    pDataC = read.csv(file.path(sourcepath,"pDataC.csv"),row.names=1)
    trueC = read.csv(file.path(sourcepath,"trueC.csv"),row.names=1)
    trueP = read.csv(file.path(sourcepath,"trueP.csv"),row.names=1)

    # remove prefix
    fn_T = file.path(sourcepath,"T.csv")
    fn_pDataC = file.path(sourcepath,"pDataC.csv")
    fn_C = file.path(sourcepath,"C.csv")
    fn_trueC = file.path(sourcepath,"trueC.csv")
    fn_trueP = file.path(sourcepath,"trueP.csv")
}else if(file.exists(fn_C)){
    print("case 3: no ground truth, with prefix")
    C = read.csv(fn_C,row.names=1)
    T = read.csv(fn_T,row.names=1)
    pDataC = read.csv(fn_pDataC,row.names=1)
    trueC = NULL
    trueP = NULL

    print(paste0("C (",nrow(C),",",ncol(C),") mean:",mean(as.matrix(C))))
    print(paste0("pDataC (",nrow(pDataC),",",ncol(pDataC),") class:",class(pDataC)))
    print(paste0("T (",nrow(T),",",ncol(T),") mean:",mean(as.matrix(T))))

}else if(file.exists(file.path(sourcepath,"C.csv"))){
    print("case 4: no ground truth, no prefix")
    C = read.csv(file.path(sourcepath,"C.csv"),row.names=1)
    T = read.csv(file.path(sourcepath,"T.csv"),row.names=1)
    pDataC = read.csv(file.path(sourcepath,"pDataC.csv"),row.names=1)
    trueC = NULL
    trueP = NULL

    # remove prefix
    fn_T = file.path(sourcepath,"T.csv")
    fn_pDataC = file.path(sourcepath,"pDataC.csv")
    fn_C = file.path(sourcepath,"C.csv")
}

#### subsample the data
target = as.integer(Sys.getenv("SUBTARGETSIZE")) # 1500
mincells = as.integer(Sys.getenv("SUBMINCELLS")) # 5
generate = as.numeric(Sys.getenv("SUBGENERATE"))
set.seed(simrep)
sampn = min(ncol(C),target)
C = as.data.frame(t(as.data.frame(t(C)) %>% sample_n(sampn)))
C = subsetsc(scremoutlier(C),generate=generate,return_obj=TRUE,nsd=3)               
pDataC = pDataC[colnames(C),]

cttab = table(pDataC$cellType)
toinclude = which(pDataC$cellType %in% names(which(cttab >= mincells)))

pDataC = pDataC[toinclude,]
C = C[,rownames(pDataC)]
T = T[rownames(C),]

print(paste0("subsampled C (",nrow(C),",",ncol(C),") mean:",mean(as.matrix(C))))
print(paste0("subsampled T (",nrow(T),",",ncol(T),") mean:",mean(as.matrix(T))))

print("environment:")
print(Sys.getenv())

if(imethod=="dropout"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod,",simrep_",simrep))
    dir.create(savepath,recursive=TRUE)
    print("running null model")
    set.seed(simrep)

    impresult=run_null(
        path=savepath,
        scdata=C,
        imputebenchmark = trueC)

    # save cluster metrics
    logdf <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(0)))
    write.csv(logdf,file.path(savepath,paste0(imethod,"_logdf.csv")))
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}else if(imethod=="DrImpute"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod,",simrep_",simrep))
    dir.create(savepath,recursive=TRUE)
    print("running drimpute")
    set.seed(simrep)
    
    logdf0 <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(0)))

    Start=Sys.time()
    impresult=run_drimpute(
        path=savepath,
        scdata=C,
        imputebenchmark = trueC)
    End=Sys.time()
    Start_POSIX = as.POSIXct(as.numeric(Start), origin="1970-01-01")
    End_POSIX = as.POSIXct(as.numeric(End), origin="1970-01-01")
    totaltime = difftime(End_POSIX,Start_POSIX,units="mins")

    # save cluster metrics
    logdf <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(totaltime)))
    write.csv(rbind(logdf0,logdf),file.path(savepath,paste0(imethod,"_logdf.csv")))
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}else if(imethod=="SCRABBLE"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod,",simrep_",simrep,",sgene_",scrgenethresh,",sA_",ScrabbleAlpha,",sB_",ScrabbleBeta,",sG_",ScrabbleGamma,",errIn_",error_inner_threshold,",errOut_",error_out_threshold,",nIterOut_",ScrnIterOuter,",nIterIn_",ScrnIterInner,",nSDC_",ScrnSDCIters))
    dir.create(savepath,recursive=TRUE)
    print("running scrabble with outer stats")
    set.seed(simrep)

    logdf0 <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(0)))

    Start=Sys.time()
    impresult_list=run_scrabble(
        path=savepath,
        scrabble_parameters = c(ScrabbleAlpha,ScrabbleBeta,ScrabbleGamma),
        error_out_threshold = error_out_threshold,
        error_inner_threshold = error_inner_threshold,
        scdata=C,
        bulkdata=T,
        nIter_outer = ScrnIterOuter,
        nIter_inner = ScrnIterInner,
        nSDCIters = ScrnSDCIters,
        outerStats = runouterstats,
        metadata=pDataC,
        imputebenchmark = trueC,
        runstability = runstability,
        useIrlba=useirlba)
    End=Sys.time()
    Start_POSIX = as.POSIXct(as.numeric(Start), origin="1970-01-01")
    End_POSIX = as.POSIXct(as.numeric(End), origin="1970-01-01")
    totaltime = difftime(End_POSIX,Start_POSIX,units="mins")

    impresult = impresult_list[["C"]]
    # save cluster metrics
    logdf <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(totaltime)))
    write.csv(rbind(logdf0,logdf),file.path(savepath,paste0(imethod,"_logdf.csv")))
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}else if(imethod=="mtSCRABBLE"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod,",simrep_",simrep,",sgene_",scrgenethresh,",sA_",ScrabbleAlpha,",sB_",ScrabbleBeta,",sG_",ScrabbleGamma,",errIn_",error_inner_threshold,",errOut_",error_out_threshold,",nIterOut_",ScrnIterOuter,",nIterIn_",ScrnIterInner,",nSDC_",ScrnSDCIters))
    dir.create(savepath,recursive=TRUE)
    print("running mtscrabble with outer stats")
    set.seed(simrep)

    logdf0 <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(0)))

    Start=Sys.time()
    impresult_list=run_scrabble_m(
        path=savepath,
        scrabble_parameters = c(ScrabbleAlpha,ScrabbleBeta,ScrabbleGamma),
        error_out_threshold = error_out_threshold,
        error_inner_threshold = error_inner_threshold,
        scdata=C,
        metadata=pDataC,
        bulkdata=T,
        nIter_outer = ScrnIterOuter,
        nIter_inner = ScrnIterInner,
        nSDCIters = ScrnSDCIters,
        thetahat = t(trueP),
        outerStats = runouterstats,
        imputebenchmark = trueC,
        runstability = runstability,
        useIrlba=useirlba)
    End=Sys.time()
    Start_POSIX = as.POSIXct(as.numeric(Start), origin="1970-01-01")
    End_POSIX = as.POSIXct(as.numeric(End), origin="1970-01-01")
    totaltime = difftime(End_POSIX,Start_POSIX,units="mins")

    impresult = impresult_list[["C"]]

    # save cluster metrics
    logdf <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(totaltime)))
    write.csv(rbind(logdf0,logdf),file.path(savepath,paste0(imethod,"_logdf.csv")))
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}else if(imethod=="DURIAN"){
    print("running durian with outer stats")

    if(deconv_method=="dsLDA"){
        print("setup juliacall")
        library(JuliaCall)
        julia_setup(JULIA_HOME = Sys.getenv("JULIA_HOME"),verbose=TRUE,rebuild = TRUE,install=FALSE)
        julia_library("DistributedStwdLDA")
        julia_library("Distributed")
        nwrk_target = as.integer(Sys.getenv("nCoresAvail"))
        julia_call("procsN",nwrk_target)
        # julia_command("@everywhere using Revise")
        julia_command("@everywhere using DistributedStwdLDA")
    }
    savepath = file.path(datapath,paste0("imputemodel_",imethod,".",deconv_method,",simrep_",simrep,",dgene_",deconvgenethresh,",sgene_",scrgenethresh,",sA_",DurianAlpha,",sB_",DurianBeta,",sG_",DurianGamma,",errIn_",error_inner_threshold,",errOut_",error_out_threshold,",nIterOut_",DunIterOuter,",nIterIn_",DunIterInner,",nSDC_",DunSDCIters,",duEps_",durianEps,",nEM_",nEM))
    dir.create(savepath,recursive=TRUE)
    set.seed(simrep)

    impresult_list=run_durian(
        path = savepath,
        scrabble_parameters = c(DurianAlpha,DurianBeta,DurianGamma),
        error_out_threshold = error_out_threshold,
        error_inner_threshold = error_inner_threshold,
        nEM = nEM,
        scdata = C,
        metadata = pDataC,
        bulkdata = T,
        deconv_method = deconv_method,
        imputebenchmark = trueC,
        deconvbenchmark = trueP,
        MCNITER = as.integer(Sys.getenv("MCNITER")),
        MINCELLSTOPICCORP = as.integer(Sys.getenv("MINCELLSTOPICCORP")),
        MCNPARTITIONS = as.integer(Sys.getenv("MCNPARTITIONS")),
        MCNCHAINS = as.integer(Sys.getenv("MCNCHAINS")),
        MCTHINNING = as.integer(Sys.getenv("MCTHINNING")),
        MCBURNRATE = as.numeric(Sys.getenv("MCBURNRATE")),
        LDARUNQC = as.logical(Sys.getenv("LDARUNQC")),
        emDiag = as.logical(Sys.getenv("EMDIAG")),
        nIter_outer = DunIterOuter,
        nIter_inner = DunIterInner,
        nSDCIters = DunSDCIters,
        saveDeconvolutionLog = as.logical(Sys.getenv("SUMMARIZEDECONV")),
        DECONVGENETHRESH=deconvgenethresh,
        SCRGENETHRESH=scrgenethresh,
        LDASCALEFACBLK = as.numeric(Sys.getenv("LDASCALEFACBLK")),
        LDASCALESC=Sys.getenv("LDASCALESC"),
        LDASCALEBLK=Sys.getenv("LDASCALEBLK"),
        outerStats = runouterstats,
        durianEps=durianEps,
        saveImputedStep=TRUE,
        runstability = runstability,
        useIrlba=useirlba)
    impresult = impresult_list[["C"]]
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}else if(imethod=="CMFImpute"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod,",simrep_",simrep))
    dir.create(savepath,recursive=TRUE)
    print("running CMF-Impute")
    set.seed(simrep)

    fn_C = paste0(savepath,"/C.csv")
    write.csv(C,fn_C)

    logdf0 <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(0)))

    Start=Sys.time()
    impresult=run_cmfimpute(
        path=savepath,
        scdata=fn_C,
        imputebenchmark=trueC,
        lib=cmfimpute_fn)

    print(paste0("impresult (",nrow(impresult),",",ncol(impresult),") class:",class(impresult)))
    print(impresult[1:5,1:5])

    End=Sys.time()
    Start_POSIX = as.POSIXct(as.numeric(Start), origin="1970-01-01")
    End_POSIX = as.POSIXct(as.numeric(End), origin="1970-01-01")
    totaltime = difftime(End_POSIX,Start_POSIX,units="mins")

    # save cluster metrics
    logdf <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(totaltime)))
    write.csv(rbind(logdf0,logdf),file.path(savepath,paste0(imethod,"_logdf.csv")))
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}else if(imethod=="G2S3"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod,",simrep_",simrep))
    dir.create(savepath,recursive=TRUE)
    print("running G2S3")
    set.seed(simrep)

    fn_C = paste0(savepath,"/C.csv")
    write.csv(C,fn_C)

    logdf0 <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(0)))

    Start=Sys.time()
    impresult=run_g2s3(
        path=savepath,
        scdata=fn_C,
        imputebenchmark=trueC,
        lib=g2s3_fn)

    End=Sys.time()
    Start_POSIX = as.POSIXct(as.numeric(Start), origin="1970-01-01")
    End_POSIX = as.POSIXct(as.numeric(End), origin="1970-01-01")
    totaltime = difftime(End_POSIX,Start_POSIX,units="mins")

    # save cluster metrics
    logdf <- data.frame(
        iter = as.integer(c(NA)),
        ldaMeanRhat = as.numeric(c(NA)),
        scrabbleLoss = as.numeric(c(NA)),
        converged=as.integer(c(1)),
        status=c("converged"),
        wallclock=as.numeric(c(totaltime)))
    write.csv(rbind(logdf0,logdf),file.path(savepath,paste0(imethod,"_logdf.csv")))
    run_cluster_plots(imputedC=impresult,pdataC=pDataC,trueC=trueC,savepath=file.path(savepath,"cluster_plots"))
    write.csv(pDataC,file.path(savepath,"pDataC.csv"))
}

# output the environment info
sessionInfo()