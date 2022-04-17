######## the null model ##############
run_null <- function(
  path=NULL,
  scdata=NULL,
  imputebenchmark=NULL){
  
  # write the data
  write.csv(scdata,file.path(path, "imputed_C.csv"))

  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  # compute the imputation loss
  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    true_sub = imputebenchmark
    result_sub = scdata

    impute_metrics = getmetrics(result_sub,true_sub)
    dropout_rate = getdroprate(result_sub,true_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("dropout"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }
  return(scdata)
}

######## DrImpute ##############
run_drimpute <- function(
  scdata=NULL,
  path=NULL,
  imputebenchmark=NULL){
    
  # impute the data using DrImpute  
  exdata <- DrImpute::DrImpute(Matrix::Matrix(as.matrix(scdata)))
  
  # write the data
  write.csv(exdata,file.path(path, "imputed_C.csv"))

  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  # compute the imputation loss
  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    true_sub = imputebenchmark
    result_sub = exdata

    impute_metrics = getmetrics(result_sub,true_sub)

    dropout_rate = getdroprate(result_sub,true_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("DrImpute"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }

  return(as.data.frame(exdata))  
}

######## URSM ##############
save_ursm_tmp_files <- function(
  C_fn="/dfs6/pub/mkarikom/code/URSM/demo/demo_data/C.csv",
  T_fn="/dfs6/pub/mkarikom/code/URSM/demo/demo_data/T.csv",
  pDataC_fn="/dfs6/pub/mkarikom/code/URSM/demo/demo_data/pDataC.csv",
  sctmpfn = NULL,
  bktmpfn = NULL,
  cttmpfn = NULL,
  scdata=NULL,
  bulkdata=NULL,
  metadata=NULL){

  # ensure the following was run:
  # conda create --name ursm27 python=2.7
  # conda activate ursm27
  # conda install --channel "conda-forge" pypolyagamma

  ursmScript = "/dfs6/pub/mkarikom/code/DTMwork/slurm/URSM/scUnif_LinuxEnv.py" # main ursm script
  indoffset = 0

  # load the data
  if(is.null(scdata)){
    scdata = read.csv(C_fn,row.names=1)
    bulkdata = read.csv(T_fn,row.names=1)
    metadata = read.csv(pDataC_fn,row.names=1)
  }

  ursmsc = t(scdata)
  ursmbulk = t(bulkdata)
  ursmcelltype = as.integer(as.factor(metadata$cellType))-1
  write.table(ursmsc,sctmpfn,sep=",",col.names=F,row.names=F)
  write.table(ursmbulk,bktmpfn,sep=",",col.names=F,row.names=F)
  write.table(ursmcelltype,cttmpfn,sep=",",col.names=F,row.names=F)

}

analyze_ursm_tmp_files <- function(
  pbdir=NULL,
  outprefix = "gemout_"){
  
  imputebenchmark = file.path(pbdir,"trueC.csv")
  C_fn = file.path(pbdir,"C.csv")
  sctmpfn = file.path(pbdir,"imputemodel_URSM","ursmsc.csv")
  cttmpfn = file.path(pbdir,"imputemodel_URSM","ursmcelltype.csv")
  estAfn = file.path(pbdir,"imputemodel_URSM",paste0(outprefix,"est_A.csv"))
  expSfn = file.path(pbdir,"imputemodel_URSM",paste0(outprefix,"exp_S.csv"))

  ursmcelltype = as.vector(read.csv(cttmpfn,header=FALSE)[,1])
  estA = read.csv(estAfn,header=FALSE)
  expS = read.csv(expSfn,header=FALSE)
  ursmsc = read.csv(sctmpfn,header=FALSE)
  C = read.csv(C_fn,row.names=1)
  trueC = read.csv(imputebenchmark,row.names=1)

  if(!is.null(trueC)){
    # com_genes = sort(intersect(rownames(trueC),rownames(ursmsc)))
    # com_cells = sort(intersect(colnames(trueC),colnames(ursmsc)))
    # dropout_rate_orig = getdroprate(trueC[com_genes,com_cells],ursmsc[com_genes,com_cells])
    dropout_rate_orig = getdroprate(trueC,ursmsc)
  }

  inddrop = which(signif(expS,10) >= 1.0,arr.ind=TRUE)
  indgene = inddrop[,2]
  indcell = inddrop[,1]
  indsignature = ursmcelltype[indcell] + 1
  celldepth = rowSums(expS)
  for(i in 1:length(indsignature)){
    # print(paste0("updating cell:",indcell[i],",gene:",indgene[i]))
    ursmsc[indcell[i],indgene[i]] = estA[indgene[i],indsignature[i]] * celldepth[indcell[i]]
  }
  ursmsc = t(ursmsc)
  colnames(ursmsc) = colnames(C)
  rownames(ursmsc) = rownames(C)

  impute_metrics = getmetrics(trueC,ursmsc)
  dropout_rate = getdroprate(ursmsc,trueC)

  impute_rmse = impute_metrics[["rmse"]]
  cor_gene = impute_metrics[["row"]]
  cor_cell = impute_metrics[["col"]]
  errnorm = impute_metrics[["errnorm"]]
  mean_cor_gene = impute_metrics[["mean_row"]]
  mean_cor_cell = impute_metrics[["mean_col"]]

  write.csv(data.frame(
    modelname=c("URSM"),
    ENORM=c(errnorm),
    RMSE=c(impute_rmse),
    Gene=c(cor_gene),
    Cell=c(cor_cell),
    MeanGene=c(mean_cor_gene),
    MeanCell=c(mean_cor_cell),
    Dropout=c(dropout_rate),
    dropout_rate=c(dropout_rate_orig)),file.path(pbdir,"imputemodel_URSM","imputation_loss.csv"))
  write.csv(ursmsc,file=file.path(pbdir,"imputemodel_URSM","imputed_C.csv"))
  return(ursmsc)
}

######## ALRA ##############
run_alra <- function(
  scdata=NULL,
  path=NULL,
  imputebenchmark=NULL){
    
  # impute the data using ALRA

  A_norm <- normalize_data(t(scdata))
  k_choice <- choose_k(A_norm)
  exdata <- Matrix::Matrix(t(alra(A_norm,k=k_choice$k)[[3]]))

  # write the data
  write.csv(exdata,file.path(path, "imputed_C.csv"))

  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  # compute the imputation loss
  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    true_sub = imputebenchmark
    result_sub = exdata

    impute_metrics = getmetrics(result_sub,true_sub)

    dropout_rate = getdroprate(result_sub,true_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("DrImpute"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }

  return(as.data.frame(exdata))
}

######## CMF-Impute ##############
run_cmfimpute <- function(
  scdata=NULL,
  path=NULL,
  imputebenchmark=NULL,
  lib=NULL){

  # impute the data using CMF-Impute
  Sys.setenv(MYSAVEPATH=path)
  Sys.setenv(MYDATA=scdata)

  myfile = paste0("\'",lib,"\'")
  system(paste0('matlab -nodesktop -nodisplay -nosplash -r "run(',myfile,');exit;"'))

  exdata = t(read.csv(paste0(path,"/imputed.csv"),header=FALSE))
  cellids = read.csv(paste0(path,"/cellids.csv"))
  geneids = read.csv(paste0(path,"/geneids.csv"))
  colnames(exdata) = cellids$cellID
  rownames(exdata) = geneids$geneID

  # write the data
  write.csv(exdata,file.path(path, "imputed_C.csv"))

  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  # compute the imputation loss
  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    true_sub = imputebenchmark
    result_sub = exdata

    impute_metrics = getmetrics(result_sub,true_sub)

    dropout_rate = getdroprate(result_sub,true_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("DrImpute"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }

  return(as.data.frame(exdata))
}

######## G2S3 ##############
run_g2s3 <- function(
  scdata=NULL,
  path=NULL,
  imputebenchmark=NULL,
  lib=NULL){
    
  # impute the data using CMF-Impute
  Sys.setenv(MYSAVEPATH=path)
  Sys.setenv(MYDATA=scdata)

  myfile = paste0("\'",lib,"\'")
  system(paste0('matlab -nodesktop -nodisplay -nosplash -r "run(',myfile,');exit;"'))
  system(paste0('matlab -nodesktop -nodisplay -nosplash -r "run(',myfile,');exit;"'))
  system(paste0('matlab -nodesktop -nodisplay -nosplash -r "run(',myfile,');exit;"'))

  exdata = read.csv(paste0(path,"/imputed.csv"))
  cellids = read.csv(paste0(path,"/cellids.csv"))
  geneids = read.csv(paste0(path,"/geneids.csv"))
  colnames(exdata) = cellids$cellID
  rownames(exdata) = geneids$geneID

  # write the data
  write.csv(exdata,file.path(path, "imputed_C.csv"))

  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  # compute the imputation loss
  if(!is.null(imputebenchmark)){
    # compare scdata with the benchmark true data
    print("intersect original data with result")
    result_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    result_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    true_sub = imputebenchmark
    result_sub = exdata

    impute_metrics = getmetrics(result_sub,true_sub)

    dropout_rate = getdroprate(result_sub,true_sub)

    impute_rmse = impute_metrics[["rmse"]]
    cor_gene = impute_metrics[["row"]]
    cor_cell = impute_metrics[["col"]]
    mean_cor_gene = impute_metrics[["mean_row"]]
    mean_cor_cell = impute_metrics[["mean_col"]]
    errnorm = impute_metrics[["errnorm"]]
    write.csv(data.frame(
      modelname=c("DrImpute"),
      ENORM=c(errnorm),
      RMSE=c(impute_rmse),
      Gene=c(cor_gene),
      Cell=c(cor_cell),
      MeanGene=c(mean_cor_gene),
      MeanCell=c(mean_cor_cell),
      Dropout=c(dropout_rate),
    dropout_rate=c(dropout_rate_orig)),file.path(path,"imputation_loss.csv"))
  }

  return(as.data.frame(exdata))  
}
