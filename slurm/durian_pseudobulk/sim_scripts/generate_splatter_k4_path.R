library(SingleCellExperiment)
library(splatter)
library(Matrix)
library(dplyr)


### arguments
ngenes = as.integer(Sys.getenv("NGENE"))
ncells = as.integer(Sys.getenv("NCELL"))
dropout_mid = as.numeric(Sys.getenv("DROPOUTMID"))
seedval = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
nbatch = as.integer(Sys.getenv("NBATCH"))
nbatch_pb = as.integer(Sys.getenv("NBATCHPB"))
nbatch_sc = as.integer(Sys.getenv("NBATCHSC"))
dtype = Sys.getenv("DTYPE")
gProbs = as.numeric(strsplit(Sys.getenv("gProbs"),"-")[[1]])
batchFacLoc = as.numeric(Sys.getenv("batchFacLoc"))
batchFacScale = as.numeric(Sys.getenv("batchFacScale"))
deFacLoc = as.numeric(Sys.getenv("deFacLoc"))
deFacScale = as.numeric(Sys.getenv("deFacScale"))
deProbs = as.numeric(strsplit(Sys.getenv("deProbs"),"-")[[1]])

### output
outputdir = Sys.getenv("DATAPATH")

nbatch=nbatch_pb+nbatch_sc
floorsize = as.integer(floor(ncells/nbatch))
floorsizes = rep(floorsize,nbatch)
remsize = as.integer(ncells-sum(floorsizes))
floorsizes[length(floorsizes)] = as.integer(floorsizes[length(floorsizes)] + remsize)
sum(floorsizes)==ncells

sim.true <- splatSimulate(
    batchCells = floorsizes,
    nGenes = ngenes,
    group.prob = gProbs,
    de.prob = deProbs,
    de.facLoc = deFacLoc,
    de.facScale = deFacScale,
    batch.facLoc = batchFacLoc,
    batch.facScale = batchFacScale,
    path.from = c(0,1,2,3),
    batch.rmEffect = TRUE,
    method = "paths", 
    verbose = FALSE,
    dropout.type = "experiment",
    dropout.shape = -1,
    dropout.mid =dropout_mid,
    seed = seedval)

sim.batch <- splatSimulate(
    batchCells = floorsizes,
    nGenes = ngenes,
    group.prob = gProbs,
    de.prob = deProbs,
    de.facLoc = deFacLoc,
    de.facScale = deFacScale,
    batch.facLoc = batchFacLoc,
    batch.facScale = batchFacScale,
    path.from = c(0,1,2,3),
    batch.rmEffect = FALSE,
    method = "paths", 
    verbose = FALSE,
    dropout.type = "experiment",
    dropout.shape = -1,
    dropout.mid =dropout_mid,
    seed = seedval)

pb_batch_inds = 1:nbatch_pb
sc_batch_inds = (nbatch_pb+1):(nbatch_pb+nbatch_sc)
pb_batches = paste0("Batch",pb_batch_inds)
sc_batches = paste0("Batch",sc_batch_inds)

indsc = which(sim.true@colData$Batch %in% sc_batches)

pDataC = data.frame(cellID=sim.batch@colData$Cell,cellType=sim.batch@colData$Group,sampleID=sim.batch@colData$Batch) 

dataTrueC = sim.true@assays@data$TrueCounts

dataC = as.matrix(counts(sim.batch))

pblist = list()
plist = list()
celltypes = as.character(sort(unique(sim.batch@colData$Group)))

for(i in 1:length(pb_batches)){
    indpb = which(sim.batch@colData$Batch %in% pb_batches[i])
    cellids = pDataC$cellID[indpb]
    pblist[[i]] = rowMeans(dataTrueC[,cellids])
    plist[[i]] = (table(pDataC$cellType[indpb])/length(pDataC$cellType[indpb]))[celltypes]
}

T = as.data.frame(dplyr::bind_cols(pblist))
rownames(T) = rownames(dataTrueC)
colnames(T) = paste0("mix",1:length(pb_batch_inds))
trueP = t(as.data.frame(dplyr::bind_rows(plist)))
colnames(trueP) = paste0("mix",1:length(pb_batches))
C = dataC[,indsc]
pDataC = pDataC[indsc,]
trueC = dataTrueC[,indsc]

write.csv(pDataC,file=file.path(outputdir,"pDataC.csv"))
write.csv(trueC,file=file.path(outputdir,"trueC.csv"))
write.csv(C,file=file.path(outputdir,"C.csv"))
write.csv(T,file=file.path(outputdir,"T.csv"))
write.csv(trueP,file=file.path(outputdir,"trueP.csv"))