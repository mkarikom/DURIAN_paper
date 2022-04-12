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

