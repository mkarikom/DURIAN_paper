# run_imputation_methods_signaling.R

library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(MuSiC)
library(Biobase)
library(xbioc)
library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(reshape2, include.only = c("melt"))

### arguments
datapath = Sys.getenv("DATAPATH")
sourcepath = Sys.getenv("SOURCEPATH")
prefix = strsplit(Sys.getenv("PREFIXTUPLE"),",")[[1]]
imethod = Sys.getenv("IMPUTE_METHOD")

### external code
source(Sys.getenv("DURIANLIB"))
source(Sys.getenv("ETCLIB"))
source(Sys.getenv("SIGNALINGLIB"))

# no slurm array for dispatching simulations, use internal numbering

fn_T = file.path(sourcepath,paste0(prefix[2],"_T.csv"))
fn_pDataC = file.path(sourcepath,paste0(prefix[1],"_pDataC.csv"))
fn_C = file.path(sourcepath,paste0(prefix[1],"_C.csv"))
print(paste0("exists T=",file.exists(fn_T)))
print(paste0("exists pdatac=",file.exists(fn_pDataC)))
print(paste0("exists c=",file.exists(fn_C)))

C = read.csv(fn_C,row.names=1)
T = read.csv(fn_T,row.names=1)
pDataC = read.csv(fn_pDataC,row.names=1)
CellChatDB = get(Sys.getenv("CELLCHATDB")) # load the object, given its name

print("environment:")
print(Sys.getenv())

if(imethod=="dropout"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running null model")
    impresult=run_null(
        path=savepath,
        scdata=C)
    run_signaling(savepath,as.matrix(impresult),pDataC,C,CellChatDB)
}else if(imethod=="DrImpute"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running drimpute")
    impresult=run_drimpute(
        path=savepath,
        scdata=C)
    run_signaling(savepath,as.matrix(impresult),pDataC,C,CellChatDB)
}else if(imethod=="SCRABBLE"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running scrabble")
    impresult_list=run_scrabble(
        path=savepath,
        scrabble_parameters = c(as.numeric(Sys.getenv("ScrabbleAlpha")), as.numeric(Sys.getenv("ScrabbleBeta")), as.numeric(Sys.getenv("ScrabbleGamma"))),
        scdata=C,
        bulkdata=T,
        nIter_outer = as.integer(Sys.getenv("ScrnIterOuter")),
        nIter_inner = as.integer(Sys.getenv("ScrnIterInner")),
        nSDCIters = as.integer(Sys.getenv("ScrnSDCIters")))
    impresult = impresult_list[["C"]]
    run_signaling(savepath,as.matrix(impresult),pDataC,C,CellChatDB)
}else if(imethod=="mtSCRABBLE"){
    savepath = file.path(datapath,paste0("imputemodel_",imethod))
    dir.create(savepath,recursive=TRUE)
    print("running mtscrabble")
    impresult_list=run_scrabble_m(
        path=savepath,
        scrabble_parameters = c(as.numeric(Sys.getenv("ScrabbleAlpha")), as.numeric(Sys.getenv("ScrabbleBeta")), as.numeric(Sys.getenv("ScrabbleGamma"))),
        scdata=C,
        metadata=pDataC,
        bulkdata=T,
        nIter_outer = as.integer(Sys.getenv("ScrnIterOuter")),
        nIter_inner = as.integer(Sys.getenv("ScrnIterInner")),
        nSDCIters = as.integer(Sys.getenv("ScrnSDCIters")))
    impresult = impresult_list[["C"]]
    run_signaling(savepath,as.matrix(impresult),pDataC,C,CellChatDB)
}else if(imethod=="DURIAN"){
    if(Sys.getenv("DECONVMETHOD")=="dsLDA"){
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
    savepath = file.path(datapath,paste0("imputemodel_",imethod,".",Sys.getenv("DECONVMETHOD")))
    dir.create(savepath,recursive=TRUE)
    print("running durian")
    impresult_list=run_durian(
        path = savepath,
        scrabble_parameters = c(as.numeric(Sys.getenv("DurianAlpha")), as.numeric(Sys.getenv("DurianBeta")), as.numeric(Sys.getenv("DurianGamma"))),
        nEM = as.integer(Sys.getenv("nEM")),
        scdata = C,
        metadata = pDataC,
        bulkdata = T,
        deconv_method = Sys.getenv("DECONVMETHOD"),
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
    impresult = impresult_list[["C"]]
    run_signaling(savepath,as.matrix(impresult),pDataC,C,CellChatDB)
}

# output the environment info
sessionInfo()