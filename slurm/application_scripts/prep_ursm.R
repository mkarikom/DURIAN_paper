# Code adapted from https://github.com/favilaco/deconv_benchmark:
# 
# Copyright 2019 Francisco Avila Cobos, 2021 Matt Karikomi
library(doParallel)

print("preparing ursm files")

### arguments
datapath = Sys.getenv("DATAPATH")
etclib = Sys.getenv("ETCLIB")
output_directory = file.path(Sys.getenv("DATAPATH"),"imputemodel_URSM")

to_remove = Sys.getenv("CTREMOVE") # 'none', 'help', or something from unique(full_phenoData$cellType)
nsimulations = as.integer(Sys.getenv("NSIMINTERNAL")) # the number of times to regenerate a set of npseudobulk pseudobulk samples from the given sc source data
npseudobulk = as.integer(Sys.getenv("NPBULK")) # the number of pseudobulk samples to generate for each simulation
ctmincells = as.integer(Sys.getenv("CTMINCELLS")) # eg 50, minimum number of cells after qc for a cell type to be retained
poolsize = as.integer(Sys.getenv("POOLSZ")) # number of cells to randomly select from sc data to construct pseudobulk
scalesingle = Sys.getenv("SCALESC") # how to scale the single cell data, eg "column", "TMM"
scalepbulk = Sys.getenv("SCALEPBULK") # how to scale the pseudobulk, eg "LogNormalize"
julia_cores = as.integer(Sys.getenv("NJULIACORES"))
impute_methodlist = strsplit(Sys.getenv("IMPUTE_METHODLIST"),",")[[1]]

### external code
source(etclib)

fn_T = file.path(datapath,"T.csv")
fn_pDataC = file.path(datapath,"pDataC.csv")
fn_C = file.path(datapath,"C.csv")
fn_trueC = file.path(datapath,"trueC.csv")

C = read.csv(fn_C,row.names=1)
T = read.csv(fn_T,row.names=1)
pDataC = read.csv(fn_pDataC,row.names=1)
trueC = read.csv(fn_trueC,row.names=1)

dir.create(output_directory,recursive=TRUE)

save_ursm_tmp_files(
    C_fn=fn_C,
    T_fn=fn_T,
    pDataC_fn=fn_pDataC,
    sctmpfn = file.path(output_directory,"ursmsc.csv"),
    bktmpfn = file.path(output_directory,"ursmbulk.csv"),
    cttmpfn = file.path(output_directory,"ursmcelltype.csv"))
