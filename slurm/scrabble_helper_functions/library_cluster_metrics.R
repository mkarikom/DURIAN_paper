run_fpc <- function(scdata,
                    pDataC,
                    n_samp_cell=800){
    n_samp_cell = min(n_samp_cell,ncol(scdata))
    umap_res_cell = umap(as.data.frame(t(scdata)) %>% sample_n(n_samp_cell))
    cstats_cell = fpc::cluster.stats(dist(umap_res_cell$layout),as.integer(as.factor(pDataC$cellType)))
}

get_model_params <- function(paramstring,sep=",",sep2="_"){
    params = strsplit(paramstring,split=sep)[[1]]
    paramlist = list()
    for(i in 1:length(params)){
        param = strsplit(params[i],split=sep2)[[1]]
        paramlist[[param[1]]] = parseargs(param[2])
    }
    paramlist
}

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

# return the percentage of nonzero values in the data 
getsparsity <- function(x){
  round(Matrix::nnzero(x == 0, na.counted = NA)/
                             (dim(x)[1]*dim(x)[2])*100)
}
