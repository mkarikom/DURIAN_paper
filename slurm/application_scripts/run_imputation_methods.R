library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(MuSiC)
library(Biobase)
library(xbioc)

### arguments
durianlib = Sys.getenv("DURIANLIB")
datapath = Sys.getenv("DATAPATH")
sourcepath = Sys.getenv("SOURCEPATH")

prefix = strsplit(Sys.getenv("PREFIXTUPLE"),",")[[1]]

etclib = Sys.getenv("ETCLIB")

nsimulations = as.integer(Sys.getenv("NSIMINTERNAL")) # the number of times to regenerate a set of npseudobulk pseudobulk samples from the given sc source data
imethod = Sys.getenv("IMPUTE_METHOD")

### external code
source(dconvlib)
source(durianlib)
source(etclib)

# no slurm array for dispatching simulations, use internal numbering

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
    C = read.csv(fn_C,row.names=1)
    T = read.csv(fn_T,row.names=1)
    pDataC = read.csv(fn_pDataC,row.names=1)
    trueC = read.csv(fn_trueC,row.names=1)
    trueP = read.csv(fn_trueP,row.names=1)
}else if(file.exists(file.path(sourcepath,"trueC.csv"))){
    C = read.csv(file.path(sourcepath,"C.csv"),row.names=1)
    T = read.csv(file.path(sourcepath,"T.csv"),row.names=1)
    pDataC = read.csv(file.path(sourcepath,"pDataC.csv"),row.names=1)
    trueC = read.csv(file.path(sourcepath,"trueC.csv"),row.names=1)
    trueP = read.csv(file.path(sourcepath,"trueP.csv"),row.names=1)
}else if(file.exists(fn_C)){
    C = read.csv(fn_C,row.names=1)
    T = read.csv(fn_T,row.names=1)
    pDataC = read.csv(fn_pDataC,row.names=1)
    trueC = NULL
    trueP = NULL
}else if(file.exists(file.path(sourcepath,"C.csv"))){
    C = read.csv(file.path(sourcepath,"C.csv"),row.names=1)
    T = read.csv(file.path(sourcepath,"T.csv"),row.names=1)
    pDataC = read.csv(file.path(sourcepath,"pDataC.csv"),row.names=1)
    trueC = NULL
    trueP = NULL
}

print("environment:")
print(Sys.getenv())

if(imethod=="dropout"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running null model")
    impresult=run_null(
        path=savepath,
        scdata=C,
        imputebenchmark = trueC)
}else if(imethod=="DrImpute"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running drimpute")
    impresult=run_drimpute(
        path=savepath,
        scdata=C,
        imputebenchmark = trueC)
}else if(imethod=="SCRABBLE"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running scrabble")
    impresult=run_scrabble(
        path=savepath,
        scrabble_parameters = c(as.numeric(Sys.getenv("ScrabbleAlpha")), as.numeric(Sys.getenv("ScrabbleBeta")), as.numeric(Sys.getenv("ScrabbleGamma"))),
        scdata=C,
        bulkdata=T,
        imputebenchmark = trueC,
        nIter_outer = as.integer(Sys.getenv("ScrnIterOuter")),
        nIter_inner = as.integer(Sys.getenv("ScrnIterInner")),
        nSDCIters = as.integer(Sys.getenv("ScrnSDCIters")))
}else if(imethod=="mtSCRABBLE"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running mtscrabble")
    impresult=run_scrabble_m(
        path=savepath,
        scrabble_parameters = c(as.numeric(Sys.getenv("ScrabbleAlpha")), as.numeric(Sys.getenv("ScrabbleBeta")), as.numeric(Sys.getenv("ScrabbleGamma"))),
        scdata=C,
        metadata=pDataC,
        bulkdata=T,
        imputebenchmark = trueC,
        nIter_outer = as.integer(Sys.getenv("ScrnIterOuter")),
        nIter_inner = as.integer(Sys.getenv("ScrnIterInner")),
        nSDCIters = as.integer(Sys.getenv("ScrnSDCIters")),
        thetahat = t(trueP))
}else if(imethod=="DURIAN"){
    if(Sys.getenv("DECONVMETHOD")=="dsLDA"){
        print("setup juliacall")
        library(JuliaCall)
        julia_setup(JULIA_HOME = Sys.getenv("JULIA_HOME"),verbose=TRUE,rebuild = TRUE,install=FALSE)
        julia_library("DistributedTopicModels")
        julia_library("Distributed")
        nwrk_target = as.integer(Sys.getenv("nCoresAvail"))
        julia_call("procsN",nwrk_target)
        # julia_command("@everywhere using Revise")
        julia_command("@everywhere using DistributedTopicModels")
    }
    savepath = file.path(datapath,paste0("imputemodel_",imethod,".",Sys.getenv("DECONVMETHOD")))
    dir.create(savepath,recursive=TRUE)
    print("running durian")
    impresult=run_durian(
        path = savepath,
        scrabble_parameters = c(as.numeric(Sys.getenv("DurianAlpha")), as.numeric(Sys.getenv("DurianBeta")), as.numeric(Sys.getenv("DurianGamma"))),
        nEM = as.integer(Sys.getenv("nEM")),
        scdata = C,
        metadata = pDataC,
        bulkdata = T,
        deconv_method = Sys.getenv("DECONVMETHOD"),
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
        nIter_outer = as.integer(Sys.getenv("DunIterOuter")),
        nIter_inner = as.integer(Sys.getenv("DunIterInner")),
        nSDCIters = as.integer(Sys.getenv("DunSDCIters")),
        summarizeDeconv = as.logical(Sys.getenv("SUMMARIZEDECONV")),
        DECONVGENETHRESH=as.numeric(Sys.getenv("DECONVGENETHRESH")),
        SCRGENETHRESH=as.numeric(Sys.getenv("SCRGENETHRESH")),
        LDASCALEFACBLK = as.numeric(Sys.getenv("LDASCALEFACBLK")),
        LDASCALESC=Sys.getenv("LDASCALESC"),
        LDASCALEBLK=Sys.getenv("LDASCALEBLK"),
        durianEps=as.numeric(Sys.getenv("durianEps")))
}

# output the environment info
sessionInfo()