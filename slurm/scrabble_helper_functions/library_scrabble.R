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
                      allScr=TRUE){
  dir.create(path,recursive = TRUE)
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
  result_list = scrabble_admm(
    scdata=scdata,
    bulkdata=bulkdata,
    scrgeneids=scrgeneids,
    scrcellids=scrcellids,
    scrabble_parameters=scrabble_parameters,
    nIter_outer=nIter_outer,
    nIter_inner=nIter_inner,
    nSDCIters=nSDCIters,
    error_out_threshold=1e-7,
    error_inner_threshold = 1e-5)


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
                      allScr=TRUE){
  dir.create(path,recursive = TRUE)

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
    error_out_threshold=1e-7,
    error_inner_threshold = 1e-5,
    thetahat=thetahat)
  
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

run_durian <- function(path,
                      scrabble_parameters = c(1, 1e-06, 1e-04),
                      nIter_outer = 20,
                      nIter_inner = 20,
                      nSDCIters = 10000,
                      emDiag=FALSE,
                      nEM=10,
                      C_fn=NULL,
                      T_fn=NULL,
                      pDataC_fn=NULL,
                      scdata=NULL,
                      bulkdata=NULL,
                      metadata=NULL,
                      deconv_method="dsLDA",
                      MCNITER=5000,
                      MCNPARTITIONS=5,
                      MCNCHAINS=2,
                      MCBLOCKSIZE=5,
                      MCPLOT=TRUE,
                      LDAPHILATENT=2,
                      LDABETAPSEUDO=0.0,
                      LDABETAEPS=0.0,
                      LDAALPHA=0.1,
                      LDANSCTHRESH=3.0,
                      SCRNSCTHRESH=3.0,
                      DECONVGENETHRESH=0.01,
                      SCRGENETHRESH=0.3,
                      MINCELLSTOPICCORP=5,
                      SCRSCALESC="ident",
                      SCRSCALEFACSC=1e4,
                      LDASCALESC="ident",
                      LDASCALEBLK="ident",
                      LDASCALEFACSC=1.0,
                      LDASCALEFACBLK=10000,
                      LDAINITFLAVOR="unif",
                      LDARUNQC=FALSE,
                      LDAVERBOSE=FALSE,
                      MCTHINNING=5,
                      MCRMCHAINS=TRUE,
                      MCBURNRATE=0.25,
                      deconvbenchmark=NULL, # this is a legacy option which will be removed
                      imputebenchmark=NULL,
                      protectedgenes=c("dummygene"),
                      summarizeDeconv=FALSE,
                      ZPOSTERIOR=FALSE, # whether to cache z and do autocorrelation etc
                      maxT=0,
                      initscrabble=FALSE,
                      allScr=TRUE,
                      limrhat=10.0,
                      durianEps=1e-6){
  #############################################################################
  # read in source data
  #############################################################################
  if(is.null(scdata) & is.null(metadata) & is.null(bulkdata)){
    metadata = read.csv(pDataC_fn)[,c("cellID","cellType","sampleID")]  
    rownames(metadata) = metadata$cellID
    # truncate the bulk sample list for debugging
    bulkdata = read.csv(T_fn,row.names=1)
    if(maxT>0){
      bulkdata = bulkdata[,1:min(dim(bulkdata)[2],maxT)]
    }else{
      bulkdata = bulkdata
    }
    scdata = read.csv(C_fn,row.names=1)
  }
  if(nchar(Sys.getenv("SUBSETCELLTYPES"))>1){
    metadata = metadata %>% dplyr::filter(cellType %in% strsplit(Sys.getenv("SUBSETCELLTYPES"),"-")[[1]])
    scdata = scdata[,metadata$cellID]
  }

  if(!is.null(imputebenchmark)){
    com_genes = sort(intersect(rownames(imputebenchmark),rownames(scdata)))
    com_cells = sort(intersect(colnames(imputebenchmark),colnames(scdata)))
    dropout_rate_orig = getdroprate(imputebenchmark[com_genes,com_cells],scdata[com_genes,com_cells])
  }

  #############################################################################
  # create the diagnostic table
  #############################################################################
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
  }
  
  logdf <- data.frame(
    iter = as.integer(c(0)),
    ldaMeanRhat = as.numeric(c(NA)),
    scrabbleLoss = as.numeric(c(NA)),
    converged=as.integer(c(0)),
    status=c("running"),
    wallclock=as.numeric(c(0)))

  if(!is.null(imputebenchmark)){
    # define the diagnostic output table
    logdf$errnorm = as.numeric(c(errnorm))
    logdf$durian_rmse = as.numeric(c(impute_rmse))
    logdf$durian_genecor = as.numeric(c(cor_gene))
    logdf$durian_cellcor = as.numeric(c(cor_cell))
    logdf$durian_mean_genecor = as.numeric(c(mean_cor_gene))
    logdf$durian_mean_cellcor = as.numeric(c(mean_cor_cell))
    logdf$dropout_rate = as.numeric(c(dropout_rate))
    logdf$sparsity_true = as.numeric(c(sparsity_true))
    logdf$sparsity_obs = as.numeric(c(sparsity_obs))
  }

  if(!is.null(deconvbenchmark)){
    # define the diagnostic output table
    logdf$deconv_rmse = as.numeric(c(NA))
    logdf$deconv_cor_celltype = as.numeric(c(NA))
    logdf$deconv_cor_bulksample = as.numeric(c(NA))
    logdf$deconv_cor_mean_celltype = as.numeric(c(NA))
    logdf$deconv_cor_mean_bulksample = as.numeric(c(NA))
  }

  #############################################################################
  # run sclda on source files
  #############################################################################
  itercount=0
  emstatus=1
  iter0_path = file.path(path,paste0("emIter_",itercount))
  dir.create(iter0_path,recursive = TRUE)
  
  #############################################################################
  #
  #
  # main durian iteration
  #
  #
  #############################################################################
  # preserve scdata_orig for possible downstream comparison, use scdata for durian iterations
  #############################################################################
  scdata_orig = scdata
  deconvgeneids = subsetsc(x=scdata,generate=DECONVGENETHRESH)[["gene"]]
  deconvgeneids = deconvgeneids[!is.na(match(deconvgeneids,row.names(bulkdata)))]
  scrgeneids = subsetsc(x=scdata,generate=SCRGENETHRESH)[["gene"]]
  deconvcellids = subsetsc(x=scdata,generate=DECONVGENETHRESH)[["cell"]]
  scrcellids = subsetsc(x=scdata,generate=SCRGENETHRESH)[["cell"]]
  scrgeneids = intersect(scrgeneids,rownames(bulkdata))
  bulkids = colnames(bulkdata)
  meanrhat = NA
  #############################################################################
  # allow initial scrabble imputation 
  #############################################################################
  if(initscrabble){
    result_list = scrabble_admm(
      scdata=scdata,
      bulkdata=bulkdata,
      scrgeneids=scrgeneids,
      scrcellids=scrcellids,
      scrabble_parameters=scrabble_parameters,
      nIter_outer=nIter_outer,
      nIter_inner=nIter_inner,
      nSDCIters=nSDCIters,
      error_out_threshold=1e-7,
      error_inner_threshold = 1e-5)
    result=result_list[["imputed"]]
    scdata[scrgeneids,scrcellids] = result
  }

  Start=Sys.time()
  while(emstatus==1){
    itercount=itercount+1
    iterPath = file.path(path,paste0("emIter_",itercount))
    dir.create(iterPath, showWarnings = FALSE,recursive = TRUE)

    if(deconv_method=="MuSiC"){
        rownames(metadata) = metadata$cellID
        if(length(unique(metadata$sampleID))){
          metadata$sampleID = sample(paste0(unique(metadata$sampleID),c(".a",".b")),nrow(metadata),replace=TRUE)
        }
        C.eset <- Biobase::ExpressionSet(assayData = as.matrix(scdata[deconvgeneids,deconvcellids]),phenoData = Biobase::AnnotatedDataFrame(metadata[deconvcellids,]))
        T.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkdata[deconvgeneids,]))
        thetahat = MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',markers = NULL, normalize = FALSE, samples = 'sampleID', verbose = F)$Est.prop.weighted
    }

    if(deconv_method=="dsLDA"){
      meanrhat=limrhat+1
      niterset = MCNITER
      while(meanrhat > limrhat){
        jl_output = julia_call("DistributedTopicModels.dsLDA_E_step",
                          t(as.matrix(scdata[deconvgeneids,deconvcellids])), # the name of the single cell data file
                          deconvgeneids,
                          deconvcellids,
                          metadata[which(deconvcellids %in% metadata$cellID),], # the name of the single cell metadata file
                          t(as.matrix(bulkdata[deconvgeneids,])),
                          deconvgeneids,
                          bulkids,
                          c(""),
                          nparts=MCNPARTITIONS,
                          runqc=LDARUNQC,
                          ldagenethresh=DECONVGENETHRESH,
                          scrgenethresh=SCRGENETHRESH,
                          ldanscthresh=LDANSCTHRESH,
                          scrnscthresh=SCRNSCTHRESH,
                          minCellsTopicCorp=MINCELLSTOPICCORP,
                          scalebulk=LDASCALEBLK,
                          bulkfactor=LDASCALEFACBLK,
                          scalesc=LDASCALESC,
                          betapseudo=LDABETAPSEUDO,
                          scfactor=LDASCALEFACSC,
                          betaeps=LDABETAEPS,
                          nchains=MCNCHAINS,
                          alpha=LDAALPHA,
                          philatent=LDAPHILATENT,
                          blocksize=MCBLOCKSIZE,
                          niter=niterset,
                          initflavor=LDAINITFLAVOR,
                          verbose=LDAVERBOSE,
                          burn=MCBURNRATE,
                          thinning=MCTHINNING,
                          rmchains=MCRMCHAINS,
                          mcplot=MCPLOT,
                          protectedgenes=protectedgenes,
                          zPosterior=ZPOSTERIOR)
      
        #############################################################################
        # update the cell and gene ids based on qc
        #############################################################################
        thetahat=jl_output[[1]]
        meanrhat=jl_output[[2]]
        rownames(thetahat) = colnames(bulkdata)
        colnames(thetahat) = jl_output[[7]]

        deconvgeneids = jl_output[[3]]
        scrgeneids = jl_output[[4]]
        deconvcellids = jl_output[[5]]
        scrcellids = jl_output[[6]]

        if(meanrhat > limrhat){
          print(paste0("mean r-hat > ",limrhat,", doubling iterations..."))
          niterset = as.integer(niterset * 2)
        }
      }
    }

    write.csv(thetahat,file.path(path,"thetahat.csv"))

    if(summarizeDeconv){
        print("finished deconvolution, plotting summary...")

        #############################################################################
        # summarize deconvolution
        #############################################################################

        outputdir = file.path(iterPath,"deconv_plots")
        dir.create(outputdir,recursive=TRUE)

        Sc_refactor <- metadata %>% 
        group_by(sampleID, cellType) %>% 
        summarise(count = n()) %>% 
        mutate(perc = count/sum(count))
        Sc_refactor$cellType = factor(Sc_refactor$cellType,levels=sort(unique(Sc_refactor$cellType)))
        # Sc_refactor$status = as.factor(Sc_refactor$cellType)
        # Sc_refactor$sampleID = as.factor(Sc_refactor$cellType)

        print(paste0("levels are ",levels(Sc_refactor$cellType)))
        myColors = palette(rainbow(length(levels(Sc_refactor$cellType))))

        names(myColors) = levels(droplevels(Sc_refactor$cellType))[1:length(myColors)]
        colScale = scale_fill_manual(name = "celltype",values = myColors)

        ggplot(Sc_refactor, aes(x = factor(sampleID), y = perc*100, fill = cellType)) +
        geom_bar(stat="identity", width = 0.7) +
        labs(x = "donor ID", y = "percent", fill = "cell type") +
        theme_minimal(base_size = 14) +
            colScale

        ggsave(file.path(outputdir,"Sc_refactor.pdf"))


        #########################################

        # Bulk deconvoluted ratios
        premelt = as.data.frame(thetahat)
        write.csv(premelt,file.path(outputdir,"P.csv"))
        premelt$name = rownames(thetahat)
        Bulk_refactor = reshape2::melt(premelt,id.vars="name")
        colnames(Bulk_refactor) = c("sampleID","cellType","perc")

        Bulk_refactor$cellType = factor(Bulk_refactor$cellType,levels=sort(unique(Bulk_refactor$cellType)))

        ggplot(Bulk_refactor, aes(x = factor(sampleID), y = perc*100, fill = cellType)) +
        geom_bar(stat="identity", width = 0.7) +
        labs(x = "donor ID", y = "percent", fill = "cell type") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 90)) +
            colScale

        ggsave(file.path(outputdir,"Bulk_refactor.pdf"))

    }

    print("finished deconvolution, creating A matrix...")
    
    #############################################################################
    # run mtscrabble (create deconvolution coefficients)
    #############################################################################

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
      error_out_threshold=1e-7,
      error_inner_threshold = 1e-5,
      thetahat=thetahat)
        
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
        print("find imputation rmse")
        impute_rmse = impute_metrics[["rmse"]]
        cor_gene = impute_metrics[["row"]]
        cor_cell = impute_metrics[["col"]]
        mean_cor_gene = impute_metrics[["mean_row"]]
        mean_cor_cell = impute_metrics[["mean_col"]]
        errnorm = impute_metrics[["errnorm"]]
        sparsity_true = getsparsity(true_sub)
        sparsity_obs = getsparsity(result_sub)
    }

    if(!is.null(deconvbenchmark)){
        # compare scdata with the benchmark true data
        observed_values = as.matrix(thetahat)
        expected_values = t(deconvbenchmark[colnames(observed_values),])

        deconv_metrics = getmetrics(observed_values,expected_values)
        deconv_rmse = deconv_metrics[["rmse"]]
        deconv_cor_celltype = deconv_metrics[["row"]] # the mean (over pseudobulk samples) correlation between the predicted celltype percentages and the true
        deconv_cor_bulksample = deconv_metrics[["col"]] # not used
        deconv_cor_mean_celltype = deconv_metrics[["mean_col"]]
        deconv_cor_mean_bulksample = deconv_metrics[["mean_row"]]
        
        if(summarizeDeconv){
          outputdir = file.path(iterPath,"deconv_tables")
          dir.create(outputdir,recursive=TRUE)
          write.csv(observed_values,file.path(outputdir,"observed.csv"))
          write.csv(expected_values,file.path(outputdir,"expected.csv"))
        }
    }

    # calculate iteration run time
    End=Sys.time()
    Start_POSIX <- as.POSIXct(as.numeric(Start), origin="1970-01-01")
    End_POSIX <- as.POSIXct(as.numeric(End), origin="1970-01-01")
    totaltime=difftime(End_POSIX,Start_POSIX)

    if(!emDiag && itercount == 1){
      print("cond 1")
      # first iteration E-step succeeds, continue
      iteroutput <- c(as.integer(itercount),meanrhat,scrabble_loss,as.integer(0),"running",totaltime)
      names(iteroutput) = c("iter","ldaMeanRhat","scrabbleLoss","converged","status","wallclock")
      if(!is.null(imputebenchmark)){
        iteroutput$durian_rmse = impute_rmse
        iteroutput$durian_genecor = cor_gene
        iteroutput$durian_cellcor = cor_cell
        iteroutput$durian_mean_genecor = mean_cor_gene
        iteroutput$durian_mean_cellcor = mean_cor_cell
        iteroutput$dropout_rate = dropout_rate
        iteroutput$errnorm = errnorm
        iteroutput$sparsity_true = sparsity_true
        iteroutput$sparsity_obs = sparsity_obs
      }
      if(!is.null(deconvbenchmark)){
        iteroutput$deconv_rmse = deconv_rmse
        iteroutput$deconv_cor_celltype = deconv_cor_celltype
        iteroutput$deconv_cor_bulksample = deconv_cor_bulksample
        iteroutput$deconv_cor_mean_celltype = deconv_cor_mean_celltype
        iteroutput$deconv_cor_mean_bulksample = deconv_cor_mean_bulksample
      }
      logdf[itercount+1,] = iteroutput[colnames(logdf)]
      write.csv(logdf,file.path(path,"durian_logdf.csv"))
    }else if(!emDiag && is.na(scrabble_loss)){
      print("cond 2")
      # E-step succeeds but mtscrabble increases, stop
      emstatus=0
      logdf$converged[itercount] = "next step terminated with NA loss"
      # do not write sc data (EM succeeds, use output of nth-1 iteration)
    # }else if(!emDiag && scrabble_loss > as.numeric(logdf$scrabbleLoss[itercount])){
    # }else if(!emDiag && abs(scrabble_loss - as.numeric(logdf$scrabbleLoss[itercount])) > durianEps*as.numeric(logdf$scrabbleLoss[itercount])){
    }else if(!emDiag && abs(scrabble_loss - as.numeric(logdf$scrabbleLoss[itercount])) <= durianEps){
      print("cond 3")
      # E-step succeeds but mtscrabble increases, stop
      emstatus=0
      # logdf$converged[itercount-1] = "change(loss) < eps*loss(t-1)"
      # do not write sc data (EM succeeds, use output of nth-1 iteration)
      iteroutput <- c(as.integer(itercount),meanrhat,scrabble_loss,as.integer(0),"change(scrabble loss) <= eps",totaltime)
      names(iteroutput) = c("iter","ldaMeanRhat","scrabbleLoss","converged","status","wallclock")
      if(!is.null(imputebenchmark)){
        iteroutput$durian_rmse = impute_rmse
        iteroutput$durian_genecor = cor_gene
        iteroutput$durian_cellcor = cor_cell
        iteroutput$durian_mean_genecor = mean_cor_gene
        iteroutput$durian_mean_cellcor = mean_cor_cell
        iteroutput$dropout_rate = dropout_rate
        iteroutput$errnorm = errnorm
        iteroutput$sparsity_true = sparsity_true
        iteroutput$sparsity_obs = sparsity_obs
      }
      if(!is.null(deconvbenchmark)){
        iteroutput$deconv_rmse = deconv_rmse
        iteroutput$deconv_cor_celltype = deconv_cor_celltype
        iteroutput$deconv_cor_bulksample = deconv_cor_bulksample
        iteroutput$deconv_cor_mean_celltype = deconv_cor_mean_celltype
        iteroutput$deconv_cor_mean_bulksample = deconv_cor_mean_bulksample
      }
      logdf[itercount+1,] = iteroutput[colnames(logdf)]
      write.csv(logdf,file.path(path,"durian_logdf.csv"))
    }else if(itercount >= nEM){
      print("cond 4")
      # E-step succeeds, M-step succeeds, but iteration limit reached, stop
      emstatus=0
      iteroutput <- c(as.integer(itercount),meanrhat,scrabble_loss,as.integer(0),"reached iteration limit",totaltime)
      names(iteroutput) = c("iter","ldaMeanRhat","scrabbleLoss","converged","status","wallclock")
      if(!is.null(imputebenchmark)){
        iteroutput$durian_rmse = impute_rmse
        iteroutput$durian_genecor = cor_gene
        iteroutput$durian_cellcor = cor_cell
        iteroutput$durian_mean_genecor = mean_cor_gene
        iteroutput$durian_mean_cellcor = mean_cor_cell
        iteroutput$dropout_rate = dropout_rate
        iteroutput$errnorm = errnorm
        iteroutput$sparsity_true = sparsity_true
        iteroutput$sparsity_obs = sparsity_obs
      }
      if(!is.null(deconvbenchmark)){
        iteroutput$deconv_rmse = deconv_rmse
        iteroutput$deconv_cor_celltype = deconv_cor_celltype
        iteroutput$deconv_cor_bulksample = deconv_cor_bulksample
        iteroutput$deconv_cor_mean_celltype = deconv_cor_mean_celltype
        iteroutput$deconv_cor_mean_bulksample = deconv_cor_mean_bulksample
      }
      logdf[itercount+1,] = iteroutput[colnames(logdf)]
      write.csv(logdf,file.path(path,"durian_logdf.csv"))
    }else{
      print("cond 5")
      # E-step succeeds, M-step succeeds, continue
      iteroutput <- c(as.integer(itercount),meanrhat,scrabble_loss,as.integer(0),"running",totaltime)
      names(iteroutput) = c("iter","ldaMeanRhat","scrabbleLoss","converged","status","wallclock")
      if(!is.null(imputebenchmark)){
        iteroutput$durian_rmse = impute_rmse
        iteroutput$durian_genecor = cor_gene
        iteroutput$durian_cellcor = cor_cell
        iteroutput$durian_mean_genecor = mean_cor_gene
        iteroutput$durian_mean_cellcor = mean_cor_cell
        iteroutput$dropout_rate = dropout_rate
        iteroutput$errnorm = errnorm
        iteroutput$sparsity_true = sparsity_true
        iteroutput$sparsity_obs = sparsity_obs
      }
      if(!is.null(deconvbenchmark)){
        iteroutput$deconv_rmse = deconv_rmse
        iteroutput$deconv_cor_celltype = deconv_cor_celltype
        iteroutput$deconv_cor_bulksample = deconv_cor_bulksample
        iteroutput$deconv_cor_mean_celltype = deconv_cor_mean_celltype
        iteroutput$deconv_cor_mean_bulksample = deconv_cor_mean_bulksample
      }
      logdf[itercount+1,] = iteroutput[colnames(logdf)]
      write.csv(logdf,file.path(path,"durian_logdf.csv"))
    }
    # on each iteration, update the EM iteration output log 
  }
  #############################################################################
  # write imputed sc matrix
  #############################################################################
  
  if(!is.null(imputebenchmark)){
    write.csv(data.frame(
      modelname=c(paste0("DURIAN.",deconv_method)),
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
  return(list(P=thetahat,C=scdata))
}

scrabble_admm <- function(
  scdata=NULL,
  bulkdata=NULL,
  scrgeneids=NULL,
  scrcellids=NULL,
  scrabble_parameters=NULL,
  nIter_outer=NULL,
  nIter_inner=NULL,
  nSDCIters=NULL,
  error_out_threshold=1e-7,
  error_inner_threshold = 1e-5){

  data1 = list()
  # single cell data with dropped out values
  data1[[1]] = as.matrix(scdata[scrgeneids,scrcellids])+0.0 # Genes x Cells
  data1[[2]] = as.matrix(bulkdata[scrgeneids,])+0.0 # Genes x Bulk
  # run scrabble (result is the updated data$data_dropout: Genes x Cells)
  

  # data1[[3]] = data1[[3]]+matrix(rnorm(length(data1[[3]])),ncol=dim(data1[[3]])[2])*noisemult  
  # run scrabble (result is the updated data$data_dropout: Genes x Cells)

  result_list = SCRABBLE::scrabble(data1,
                          parameter = scrabble_parameters,
                          nIter_outer = nIter_outer,
                          nIter_inner = nIter_inner,
                          nSDCIters = nSDCIters,
                          error_out_threshold = error_out_threshold,
                          error_inner_threshold = error_inner_threshold)
  return(result_list)

}

mtscrabble_admm <- function(
  scdata=NULL,
  bulkdata=NULL,
  metadata=NULL,
  scrgeneids=NULL,
  scrcellids=NULL,
  scrabble_parameters=NULL,
  nIter_outer=NULL,
  nIter_inner=NULL,
  nSDCIters=NULL,
  error_out_threshold=1e-7,
  error_inner_threshold = 1e-5,
  thetahat=NULL){

  ncells = dim(scdata[scrgeneids,scrcellids])[2]
  nbulk = dim(bulkdata[scrgeneids,])[2]
  celltypes = unique(metadata$cellType)
  A = matrix(0,nbulk,ncells)
  for(i in 1:nbulk){
    for(j in 1:ncells){
      celltype=which(metadata$cellType[j] %in% celltypes)
      A[i,j]=thetahat[i,celltype]/ncells
    }
  }
  
  data1 = list()

  data1[[1]] = as.matrix(scdata[scrgeneids,scrcellids])+0.0 # Genes x Cells
  data1[[2]] = as.matrix(bulkdata[scrgeneids,])+0.0 # Genes x Bulk
  data1[[3]] = A # Bulk x Cells
  # run scrabble (result is the updated data$data_dropout: Genes x Cells)
  result_list = SCRABBLE::scrabble(data1,
                          parameter = scrabble_parameters,
                          nIter_outer = nIter_outer,
                          nIter_inner = nIter_inner,
                          nSDCIters = nSDCIters,
                          error_out_threshold = 1e-7,
                          error_inner_threshold = 1e-5,
                          multi=TRUE)
  return(result_list)
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
    doParallel::registerDoParallel(corese = ncores)
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