######## SCRABBLE OG ##############
run_scrabble <- function(
                      path, 
                      scrabble_parameters = c(1, 1e-06, 1e-04),
                      nIter_outer = 20,
                      nIter_inner = 20,
                      nSDCIters = 10000,
                      C_fn="ldaC.csv.gz",
                      T_fn="T.csv.gz",
                      scdata=NULL,
                      bulkdata=NULL,
                      scrcellids=NULL,
                      scrgeneids=NULL,
                      imputebenchmark=NULL,
                      allScr=TRUE,
                      error_out_threshold=1e-7,
                      error_inner_threshold=1e-5,
                      outerStats = FALSE,
                      metadata=NULL,
                      useIrlba=TRUE,
                      runstability=FALSE){
  dir.create(path,recursive = TRUE)
  print(paste0("library_scrabble_clusterMetrics: run_scrabble outerStats =",outerStats))

  #############################################################################
  # read in source data
  #############################################################################
  if(is.null(scdata)){
    bulkdata = read.csv(T_fn,row.names=1)
    scdata = read.csv(C_fn,row.names=1)
  }
  
  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  if(is.null(scrcellids)){
    scrcellids = colnames(scdata)
  }
  if(is.null(scrgeneids)){
    scrgeneids = rownames(scdata)
  }
  # genezero = c(names(which(rowSums(scdata)==0)),names(which(rowSums(bulkdata)==0)))  
  # scrgeneids = setdiff(scrgeneids,genezero)
  scrgeneids = intersect(scrgeneids,rownames(bulkdata))
  if(outerStats){
    result_list = scrabble_admm(
      scdata=scdata,
      bulkdata=bulkdata,
      scrgeneids=scrgeneids,
      scrcellids=scrcellids,
      scrabble_parameters=scrabble_parameters,
      nIter_outer=nIter_outer,
      nIter_inner=nIter_inner,
      nSDCIters=nSDCIters,
      error_out_threshold=error_out_threshold,
      error_inner_threshold=error_inner_threshold,
      outerStatsPath = file.path(path,"outerStats"),
      metadata = metadata,
      useIrlba=useIrlba)
  }else{
    result_list = scrabble_admm(
      scdata=scdata,
      bulkdata=bulkdata,
      scrgeneids=scrgeneids,
      scrcellids=scrcellids,
      scrabble_parameters=scrabble_parameters,
      nIter_outer=nIter_outer,
      nIter_inner=nIter_inner,
      nSDCIters=nSDCIters,
      error_out_threshold=error_out_threshold,
      error_inner_threshold=error_inner_threshold,
      metadata = metadata,
      useIrlba=useIrlba)
  }


  result=result_list[["imputed"]]
  scrabble_loss=result_list[["error"]]

  scdata[scrgeneids,scrcellids] = result

  rownames(result) = rownames(scdata[scrgeneids,scrcellids])
  colnames(result) = colnames(scdata[scrgeneids,scrcellids])


  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    if(allScr){
      true_sub = imputebenchmark
      result_sub = scdata
    }else{
      true_sub = imputebenchmark[result_genes,result_cells]
      result_sub = scdata[result_genes,result_cells]
    }

    impute_metrics = getmetrics(result_sub,true_sub)

    dropout_rate = getdroprate(result_sub,true_sub)
    sparsity_true = getsparsity(true_sub)
    sparsity_obs = getsparsity(result_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("SCRABBLE"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
      sparsity_true=c(sparsity_true),
      sparsity_obs=c(sparsity_obs),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }

  write.csv(scdata,file.path(path,"imputed_C.csv"))
  return(list(C=scdata))
}

######## SCRABBLE Multitask ##############
run_scrabble_m <- function(
                      path, 
                      scrabble_parameters = c(1, 1e-06, 1e-04),
                      nIter_outer = 20,
                      nIter_inner = 20,
                      nSDCIters = 10000,
                      C_fn="ldaC.csv",
                      T_fn="T.csv",
                      pDataC_fn="pDataC.csv",
                      scdata=NULL,
                      bulkdata=NULL,
                      metadata=NULL,
                      thetahat=NULL,
                      imputebenchmark=NULL,
                      scrcellids=NULL,
                      scrgeneids=NULL,
                      allScr=TRUE,
                      error_out_threshold=1e-7,
                      error_inner_threshold=1e-5,
                      outerStats = FALSE,
                      useIrlba=TRUE,
                      runstability=FALSE){
  dir.create(path,recursive = TRUE)
  print(paste0("library_scrabble_clusterMetrics: run_scrabble_m outerStats =",outerStats))

  #############################################################################
  # read in source data
  #############################################################################
  if(is.null(scdata)){
    metadata = read.csv(pDataC_fn,row.names = 1)
    bulkdata = read.csv(T_fn,row.names=1)
    scdata = read.csv(C_fn,row.names=1)
  }
  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  if(is.null(scrcellids)){
    scrcellids = colnames(scdata)
  }
  if(is.null(scrgeneids)){
    scrgeneids = rownames(scdata)
  }
  scrgeneids = intersect(scrgeneids,rownames(bulkdata))
  write.csv(thetahat,file.path(path,"thetahat.csv"))
  if(outerStats){
    result_list = mtscrabble_admm(
      scdata=scdata,
      bulkdata=bulkdata,
      metadata=metadata,
      scrgeneids=scrgeneids,
      scrcellids=scrcellids,
      scrabble_parameters=scrabble_parameters,
      nIter_outer=nIter_outer,
      nIter_inner=nIter_inner,
      nSDCIters=nSDCIters,
      error_out_threshold=error_out_threshold,
      error_inner_threshold=error_inner_threshold,
      thetahat=thetahat,
      outerStatsPath = file.path(path,"outerStats"),
      useIrlba=useIrlba)
  }else{
    result_list = mtscrabble_admm(
      scdata=scdata,
      bulkdata=bulkdata,
      metadata=metadata,
      scrgeneids=scrgeneids,
      scrcellids=scrcellids,
      scrabble_parameters=scrabble_parameters,
      nIter_outer=nIter_outer,
      nIter_inner=nIter_inner,
      nSDCIters=nSDCIters,
      error_out_threshold=error_out_threshold,
      error_inner_threshold=error_inner_threshold,
      thetahat=thetahat,
      useIrlba=useIrlba)
  }

  
  result=result_list[["imputed"]]
  scrabble_loss=result_list[["error"]]
  # update the result of imputation
  scdata[scrgeneids,scrcellids] = result
  # add annotation to result for backup
  colnames(result) = colnames(scdata[scrgeneids,scrcellids])
  rownames(result) = rownames(scdata[scrgeneids,scrcellids])

  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    if(allScr){
      true_sub = imputebenchmark
      result_sub = scdata
    }else{
      true_sub = imputebenchmark[result_genes,result_cells]
      result_sub = scdata[result_genes,result_cells]
    }

    impute_metrics = getmetrics(result_sub,true_sub)

    dropout_rate = getdroprate(result_sub,true_sub)
    sparsity_true = getsparsity(true_sub)
    sparsity_obs = getsparsity(result_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("mtSCRABBLE"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
      sparsity_true=c(sparsity_true),
      sparsity_obs=c(sparsity_obs),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }
  
  write.csv(scdata,file.path(path,"imputed_C.csv"))
  return(list(C=scdata))
}

# 1. find the correlation matrices for true and obs data:  
#     eg the correlation matrix for obs is a symmetric matrix giving correlation between each obs column and all the other obs columns
# 2. return the correlation between the vectorized obs and vectorized true
getmetrics <- function(obs.orig,true.orig){
  # remove cols or rows that are all zero
  obs.nzcol = which(colSums(obs.orig)>0)
  obs.nzrow = which(rowSums(obs.orig)>0)
  true.nzcol = which(colSums(true.orig)>0)
  true.nzrow = which(rowSums(true.orig)>0)
  nzcol = intersect(obs.nzcol,true.nzcol)
  nzrow = intersect(obs.nzrow,true.nzrow)
  true = true.orig[nzrow,nzcol]
  obs = obs.orig[nzrow,nzcol]  
  
  print("find rmse")
  rmse = sqrt(mean(as.matrix(obs-true)^2))
  # errnorm = norm(log10(as.matrix(obs) + 1) - log10(as.matrix(true) + 1), type = "2")
  errnorm = norm(log(as.matrix(obs) + 1) - log(as.matrix(true)+1), type = "2")
  if(!foreach::getDoParRegistered()){
    if(is.na(as.integer(Sys.getenv("NCPUS")))){
      Sys.setenv(NCPUS=4)
    }
    ncores = as.integer(Sys.getenv("NCPUS"))
    doParallel::registerDoParallel(cores = ncores)
  }

  # foreach version
  # excellent post here: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
  print("find mean column-wise correlation")
  mean_col = foreach::foreach(
    i = 1:ncol(obs), 
    .combine = 'c'
  ) %dopar% {
      cor(obs[,i],true[,i])
  }

  print("find mean row-wise correlation")
  mean_row = foreach::foreach(
    i = 1:nrow(obs), 
    .combine = 'c'
  ) %dopar% {
      cor(as.vector(t(obs[i,])),as.vector(t(true[i,])))
  }

  mean_cor_row = mean(mean_row[!is.na(mean_row)])
  mean_cor_col = mean(mean_col[!is.na(mean_col)])

  print("find column-wise correlation")
  obs_col = cor(as.matrix(obs))
  true_col = cor(as.matrix(true))
  cor_col = calculate_similarity(obs_col,true_col)

  print("find row-wise correlation")
  obs_row = cor(t(as.matrix(obs)))
  true_row = cor(t(as.matrix(true)))
  cor_row = calculate_similarity(obs_row,true_row)

  print("returning metrics")
  return(list(rmse=rmse,errnorm=errnorm,row=cor_row,col=cor_col,mean_row=mean_cor_row,mean_col=mean_cor_col))
}

# return the ratio of zero entries in `obs` that are non-zero in `orig` 
getdroprate <- function(obs,orig){
  nzorig = which(as.vector(as.logical(as.matrix(orig))))
  drvec = as.vector(as.logical(as.matrix(obs)))[nzorig] - as.vector(as.logical(as.matrix(orig)))[nzorig]
  drind = which(drvec == -1)
  length(drind) / length(drvec)
}

# return the percentage of nonzero values in the data 
getsparsity <- function(x){
  round(Matrix::nnzero(x == 0, na.counted = NA)/
                             (dim(x)[1]*dim(x)[2])*100)
}

calculate_similarity <- function(data1,data2){

  d = cor(c(data1[lower.tri(data1)]),c(data2[lower.tri(data2)]))

  return(d)

}

# remove outlier cells by library size
scremoutlier <- function(data,nsd=3){
  sums = colSums(data)
  msum = mean(sums)
  sdsum = sd(sums)
  inds = which(as.logical((sums >= msum-nsd*sdsum)*(sums <= msum+nsd*sdsum)))
  data[,inds]
}

# get the cell and gene ids corresponding to gene expressed in at least `generate * ncells` cells and for those genes, the cells that have nonzero total expression
subsetsc <- function(x=NULL,generate=0.05,geneids=NULL,return_obj=FALSE,nsd=NULL){
   # remove outlier cells
   if(!is.null(nsd)){
     xtmp = scremoutlier(x,nsd=nsd)
   }else{
     xtmp = x
   }
   # make sure genes are expressed at least `generate`
   if(is.null(geneids)){
     geneids = names(which(rowSums(xtmp > 0) > dim(xtmp)[2]*generate))
   }
   # make sure all cells express something
   cellids = names(which(colSums(xtmp[geneids,]) > 0))
   if(return_obj){
     return(xtmp[geneids,cellids])
   }else{
     return(list(gene=geneids,cell=cellids))
   }
}

