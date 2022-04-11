get_cellchat_genesets <- function(
                            CellChatDB,
                            limngene=10){

    interaction_input <- CellChatDB$interaction
    complex_input <- CellChatDB$complex
    cofactor_input <- CellChatDB$cofactor
    geneInfo <- CellChatDB$geneInfo

    # limngene is the minimum number of genes required in the geneset across all interacting species
    pathways = sort(unique(interaction_input$pathway_name))
    genesets = list()

    for(i in 1:length(pathways)){

        inames = interaction_input$interaction_name[interaction_input$pathway_name==pathways[i]]
        gnames = unique(unlist(sapply(inames,strsplit,"_")))
        
        coI = unique(interaction_input$co_I_receptor[interaction_input$pathway_name==pathways[i]])
        coIgenes = as.vector(as.matrix(cofactor_input[coI,]))[sapply(as.vector(as.matrix(cofactor_input[coI,])),nchar) > 0]        
        
        coA = unique(interaction_input$co_A_receptor[interaction_input$pathway_name==pathways[i]])
        coAgenes = as.vector(as.matrix(cofactor_input[coA,]))[sapply(as.vector(as.matrix(cofactor_input[coA,])),nchar) > 0]

        antagon = unique(interaction_input$antagonist[interaction_input$pathway_name==pathways[i]])
        antagongenes = as.vector(as.matrix(cofactor_input[antagon,]))[sapply(as.vector(as.matrix(cofactor_input[antagon,])),nchar) > 0]

        agon = unique(interaction_input$agonist[interaction_input$pathway_name==pathways[i]])
        agongenes = as.vector(as.matrix(cofactor_input[agon,]))[sapply(as.vector(as.matrix(cofactor_input[agon,])),nchar) > 0]

        tmpgenes = c(gnames,coIgenes,coAgenes,antagongenes,agongenes)
        genesets[[i]] = tmpgenes[!is.na(tmpgenes)]

    }

    genesets[which(lapply(genesets,length)>limngene)]
}

get_web_genesets <- function(
                            dirname,
                            limngene=50){

    # limngene is the minimum number of genes required in the geneset across all interacting species
    fnames_last = list.files(dirname)
    fnames_full = list.files(dirname,full.names=TRUE)
    genesets = list()

    inds_gene = grep("gene",fnames_last)
    inds_name = grep("name",fnames_last)

    for(i in 1:length(fnames_full[inds_gene])){
        genes = read.csv(fnames_full[inds_gene[i]],header=FALSE)[,1]
        name = read.csv(fnames_full[inds_name[i]],header=FALSE)[,1]
        genesets[[name]] = genes
    }
    genesets[which(lapply(genesets,length)>limngene)]
}

get_pathway_corr <- function(geneset,scdata,nperm = 100){
    
    set.seed(42)
    allcorr = c()

    inds = match(geneset,rownames(scdata))

    inds = sample(1:nrow(scdata),length(geneset))
    setcorr = cor(t(scdata[geneset,]))
    setcorr = abs(setcorr[upper.tri(setcorr)])
    meansetcorr0 = mean(setcorr[!is.na(setcorr)])
    for(i in 1:nperm){
        inds = sample(1:nrow(scdata),length(geneset))
        setcorr = cor(t(scdata[inds,]))
        setcorr = abs(setcorr[upper.tri(setcorr)])
        meansetcorr = mean(setcorr[!is.na(setcorr)])
        allcorr = c(allcorr,meansetcorr)
    }
    c(meansetcorr0,mean(allcorr))

}

get_pref <- function(string){
    sapply(strsplit(string,":")[[1]],strsplit,",")
}

plot_clustermetrics <- function(df_orig,outputdir,yvars = c("dPGCS")){

    nmodels = length(unique(df_orig$method))

    palette1 = scales::hue_pal()(nmodels)
    names(palette1) = sort(unique(df_orig$method))

    plist = list()
    for(j in 1:length(yvars)){
        df = filter(df_orig,metric==yvars[j])
        p = ggplot(df,
                aes(
                    x=as.factor(method), 
                    y=log(x+1),
                    fill=as.factor(method))) + 
                geom_boxplot(
                aes(
                    x=as.factor(method),
                    y=log(x+1),
                    fill=as.factor(method),
                    alpha=as.factor(dbname))) +
                scale_alpha_discrete(
                range = c(0.1, 0.8),
                guide = guide_legend(override.aes = list(fill = "black")),name="dbname") +
                facet_grid(~dataset,scale="free") +
                stat_summary(fun=mean, geom="point", shape=18, size=1,
                                aes(group=as.factor(dbname),x=as.factor(method), y=log(x+1)),
                                position = position_dodge(.75),stroke=1.5)+
                                guides(fill = FALSE)  
        # set_palette(p, palette1)
        p = p +
            theme_bw() +
                rotate_x_text(45) +
            theme(axis.text=element_text(size=25,face="bold"),
                # legend.position = "none",
                axis.title.x=element_blank(),
                panel.border=element_rect(colour="black",size=1,fill=NA),
                axis.title.y = element_text(size = 25, face = "bold"),
                strip.text.x = element_text(size = 25, face = "bold", angle = 0),
                legend.key.size = unit(3, 'cm'), #change legend key size
                legend.key.height = unit(1, 'cm'), #change legend key height
                legend.key.width = unit(1, 'cm'), #change legend key width
                legend.title = element_text(size=20,face="bold"), #change legend title font size
                legend.text = element_text(size=20),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()
                ) +
        ylab(yvars[j])
        ggsave(plot=p,file=file.path(outputdir,paste0("metrics_",yvars[j],".pdf")),width=100,height=30,units="cm")
        plist[[yvars[j]]] = p
    }
    return(plist)
}

run_clusterMetrics_final_nested_real <- function(
                        savepath_master,
                        outputmaster,
                        plotwidth=10,
                        plotheight=10,
                        excludemodels=c(),
                        maxwallclock=20,
                        subdir="plots",
                        exclude.sa=c(),
                        exclude.sb=c(),
                        exclude.sg=c()){

    metrics.orig = list()
    paramsets = list.files(outputmaster,full.names=TRUE,include.dirs=TRUE)
    pref_inds = grep("pref_",paramsets)
    for(paramind in pref_inds){
        dropsets = list.files(paramsets[paramind],full.names=TRUE,include.dirs=TRUE)
        dropsets_short = list.files(paramsets[paramind],include.dirs=TRUE)
        
        dataparams = get_model_params(tail(strsplit(paramsets[paramind],"/")[[1]],n=1))

        modeldirs = list.files(file.path(paramsets[paramind],"output_fit"),full.names=TRUE)
        modeldirs = modeldirs[grep("imputemodel_",modeldirs)]
        imethods = sapply(modeldirs,strsplit,"/") %>% sapply(last) %>% unname

        for(simind in 1:length(imethods)){
            myparams = get_model_params(imethods[simind])
            imputedir = modeldirs[simind]
            modelstring = tail(strsplit(imputedir,"/")[[1]],n=1)
            imethod = myparams[["imputemodel"]]
            ldfname = file.path(imputedir,paste0(imethod,"_logdf.csv"))
            scldfname = file.path(imputedir,"outerStats",paste0(imethod,"_logdf.csv")) # if we are running scrabble and things did not finish
            impname = file.path(imputedir,"imputed_C.csv")
            
            emiter_inds = grep("emIter_",list.files(imputedir,full.names=TRUE))

            if(file.exists(ldfname) & file.exists(impname)){
                print(paste0("processing converged metrics for simind=",simind,":",length(imethods)))
                ldf = read.csv(ldfname,row.names=1)

                data_myparams = as.data.frame(myparams) %>% dplyr::slice(rep(1:n(), each = nrow(ldf)))

                metrics = cbind(ldf,data_myparams,as.data.frame(dataparams))
                metrics$modelstring = modelstring
                metrics.orig[[length(metrics.orig)+1]] = tail(metrics,n=1)
            }else if(file.exists(scldfname)){
                # if(length(grep("DURIAN",imethod))>0){
                #     if(length(emiter_inds)>0){
                #         for(emi in emiter_inds){
                #             emiter_scrabble = file.path(list.files(imputedir,full.names=TRUE)[emi],"outerStats","SCRABBLE_logdf.csv")
                #             if(file.exists(emiter_scrabble)){
                #                 cmetrics = read.csv(emiter_scrabble,row.names=1)
                #                 cmetrics$modelname=imethod
                #                 cmetrics$simname=imethods[simind]
                #                 cmetrics$dropout_rate = droprate_backup

                #                 data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(cmetrics)))
                #                 cmetrics = cbind(cmetrics,data_cmetrics)
                #                 metrics.orig[[length(metrics.orig)+1]] = cmetrics
                #             }
                #         }
                #     }
                # }
                # print(paste0("processing not-yet-converged scrabble metrics for simind=",simind,":",length(imethods)))
                # cmetrics = read.csv(scldfname,row.names=1)
                # cmetrics$modelname=imethod
                # cmetrics$simname=imethods[simind]
                # cmetrics$dropout_rate = droprate_backup

                # data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(cmetrics)))
                # cmetrics = cbind(cmetrics,data_cmetrics)
                # metrics.orig[[length(metrics.orig)+1]] = cmetrics
                
            }else{
                print(paste0("imputation or log output is missing from: ",imputedir))
            }
        }

        
    }
    # metrics_df.orig = do.call(rbind,metrics.orig)

    metrics_df.orig = do.call(plyr::rbind.fill,metrics.orig)        

    metrics_df0 = metrics_df.orig %>% filter(!imputemodel %in% excludemodels & wallclock <= maxwallclock)
    metrics_df0 = metrics_df.orig %>% filter(!sA %in% exclude.sa & !sB %in% exclude.sb & !sG %in% exclude.sg)
    

    metrics_df0 = metrics_df0 %>% group_by(pref) %>% mutate(mean_dropout_cell_dunn = mean(cell_dunn[imputemodel=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(pref) %>% mutate(mean_dropout_gene_dunn = mean(gene_dunn[imputemodel=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(pref) %>% mutate(mean_dropout_cell_connectivity = mean(cell_connectivity[imputemodel=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(pref) %>% mutate(mean_dropout_gene_connectivity = mean(gene_connectivity[imputemodel=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(pref) %>% mutate(mean_dropout_cell_silhouette = mean(cell_silhouette[imputemodel=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(pref) %>% mutate(mean_dropout_gene_silhouette = mean(gene_silhouette[imputemodel=="dropout"]))
    
    dir.create(file.path(savepath_master,subdir),recursive=TRUE)
    fn = file.path(savepath_master,subdir,"mean_plot_point.pdf")

    metrics_df0
    for(cn in 1:length(colnames(metrics_df0))){
        if(length(grep("stab",colnames(metrics_df0)[cn]))>0){
            print("found")
            colnames(metrics)[cn] = stringr::str_sub(colnames(metrics)[cn], start=6L, end=-1L)
        }
    }
    metrics_prefix = c("dunn","connectivity","silhouette","apn","ad","adm","fom")
    metrics_prefix = c("dunn","connectivity","silhouette")
    for(iii in 1:length(metrics_prefix)){

        meanvar_cell = paste0("mean_dropout_cell_",metrics_prefix[iii])   
        meanvar_gene = paste0("mean_dropout_gene_",metrics_prefix[iii])   
        yvar_cell = paste0("gene_",metrics_prefix[iii])
        yvar_gene = paste0("cell_",metrics_prefix[iii])
        dir.create(file.path(savepath_master,subdir,"scatter"),recursive=TRUE)
        dir.create(file.path(savepath_master,subdir,"violin"),recursive=TRUE)
        dir.create(file.path(savepath_master,subdir,"box"),recursive=TRUE)
        print(paste0("plotting scatter prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"scatter",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=wallclock,y=get(yvar_cell),color=imputemodel)) + 
            geom_point() +
            geom_hline(
                aes(yintercept = get(meanvar_cell), 
                colour = imputemodel), 
                filter(metrics_df0,imputemodel=="dropout")) +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            facet_grid(rows=vars(pref)) +
            ggtitle(paste0("cell ",metrics_prefix[iii]," vs wallclock"))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting scatter prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"scatter",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=wallclock,y=get(yvar_gene),color=imputemodel)) + 
            geom_point() +
            geom_hline(
                aes(yintercept = get(meanvar_gene), 
                colour = imputemodel), 
                filter(metrics_df0,imputemodel=="dropout")) +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(pref)) +
            ggtitle(paste0("gene ",metrics_prefix[iii]," vs wallclock"))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting box prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"box",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=imputemodel,y=get(yvar_cell),color=imputemodel)) + 
            geom_boxplot() +
            ylab("Dunn Index") +
            xlab("Model") +
            rotate_x_text(45) +
            facet_grid(rows=vars(pref)) +
            theme(axis.text.x = element_text(size = 5, angle = 90, face = "plain")) +
            ggtitle(paste0("cell ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting box prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"box",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=imputemodel,y=get(yvar_gene),color=imputemodel)) + 
            geom_boxplot() +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(pref)) +
            theme(axis.text.x = element_text(size = 5, angle = 90, face = "plain")) +
            ggtitle(paste0("gene ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting violin prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"violin",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=imputemodel,y=get(yvar_cell),color=imputemodel)) + 
            geom_violin() +
            ylab("Dunn Index") +
            xlab("Model") +
            rotate_x_text(45) +
            facet_grid(rows=vars(pref)) +
            ggtitle(paste0("cell ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting violin prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"violin",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=imputemodel,y=get(yvar_gene),color=imputemodel)) + 
            geom_violin() +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(pref)) +
            ggtitle(paste0("gene ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)
    }
    return(metrics_df0)

}


run_clusterMetrics_final_nested_sim <- function(
                        savepath_master,
                        outputmaster,
                        plotwidth=10,
                        plotheight=10,
                        excludemodels=c(),
                        maxwallclock=20,
                        subdir="plots",
                        sparsityparam="dmid",
                        exclude.sa=c(),
                        exclude.sb=c(),
                        exclude.sg=c()){

    metrics.orig = list()
    paramsets = list.files(outputmaster,full.names=TRUE,include.dirs=TRUE)
    for(paramind in 1:length(paramsets)){
        dropsets = list.files(paramsets[paramind],full.names=TRUE,include.dirs=TRUE)
        dropsets_short = list.files(paramsets[paramind],include.dirs=TRUE)
        
        fit_inds = grep("output_fit",dropsets)

    
        for(fit_ind in fit_inds){
            dataparams = get_model_params(dropsets_short[fit_ind])
            simdirs = list.files(dropsets[fit_ind],full.names=TRUE)
            for(simdir in simdirs){
                modeldirs = list.files(simdir,full.names=TRUE)
                modeldirs = modeldirs[grep("imputemodel_",modeldirs)]
                imethods = sapply(modeldirs,strsplit,"/") %>% sapply(last) %>% unname
                ind_dropmodel = grep("imputemodel_dropout",modeldirs)
                if(file.exists(file.path(modeldirs[ind_dropmodel],"imputation_loss.csv"))){
                    droprate_backup = read.csv(file.path(modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$sparsity_obs
                    for(simind in 1:length(imethods)){
                        myparams = get_model_params(imethods[simind])
                        imputedir = modeldirs[simind]
                        modelstring = tail(strsplit(imputedir,"/")[[1]],n=1)
                        imethod = myparams[["imputemodel"]]
                        methodlossname = file.path(modeldirs[simind],"imputation_loss.csv")
                        ldfname = file.path(imputedir,paste0(imethod,"_logdf.csv"))
                        scldfname = file.path(imputedir,"outerStats",paste0(imethod,"_logdf.csv")) # if we are running scrabble and things did not finish
                        impname = file.path(imputedir,"imputed_C.csv")
                        
                        emiter_inds = grep("emIter_",list.files(imputedir,full.names=TRUE))

                        if(file.exists(ldfname) & file.exists(impname) & file.exists(methodlossname)){
                            print(paste0("processing converged metrics for simind=",simind,":",length(imethods)))
                            ldf = read.csv(ldfname,row.names=1)

                            methodloss = read.csv(methodlossname,row.names=1)
                            data_methodloss = methodloss %>% dplyr::slice(rep(1:n(), each = nrow(ldf)))

                            data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(ldf)))
                            
                            data_myparams = as.data.frame(myparams) %>% dplyr::slice(rep(1:n(), each = nrow(ldf)))

                            metrics = cbind(ldf,data_cmetrics,data_methodloss,data_myparams)
                            metrics$Dropout = NULL
                            metrics$dropout_rate = NULL
                            metrics$modelstring = modelstring
                            metrics.orig[[length(metrics.orig)+1]] = tail(metrics,n=1)
                        }else if(file.exists(scldfname)){
                            # if(length(grep("DURIAN",imethod))>0){
                            #     if(length(emiter_inds)>0){
                            #         for(emi in emiter_inds){
                            #             emiter_scrabble = file.path(list.files(imputedir,full.names=TRUE)[emi],"outerStats","SCRABBLE_logdf.csv")
                            #             if(file.exists(emiter_scrabble)){
                            #                 cmetrics = read.csv(emiter_scrabble,row.names=1)
                            #                 cmetrics$modelname=imethod
                            #                 cmetrics$simname=imethods[simind]
                            #                 cmetrics$dropout_rate = droprate_backup

                            #                 data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(cmetrics)))
                            #                 cmetrics = cbind(cmetrics,data_cmetrics)
                            #                 metrics.orig[[length(metrics.orig)+1]] = cmetrics
                            #             }
                            #         }
                            #     }
                            # }
                            # print(paste0("processing not-yet-converged scrabble metrics for simind=",simind,":",length(imethods)))
                            # cmetrics = read.csv(scldfname,row.names=1)
                            # cmetrics$modelname=imethod
                            # cmetrics$simname=imethods[simind]
                            # cmetrics$dropout_rate = droprate_backup

                            # data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(cmetrics)))
                            # cmetrics = cbind(cmetrics,data_cmetrics)
                            # metrics.orig[[length(metrics.orig)+1]] = cmetrics
                            
                        }else{
                            print(paste0("imputation or log output is missing from: ",imputedir))
                        }
                    }
                }

            }
        }
    }
    # metrics_df.orig = do.call(rbind,metrics.orig)

    metrics_df.orig = do.call(plyr::rbind.fill,metrics.orig)        

    metrics_df0 = metrics_df.orig %>% filter(!modelname %in% excludemodels & wallclock <= maxwallclock)
    metrics_df0 = metrics_df.orig %>% filter(!sA %in% exclude.sa & !sB %in% exclude.sb & !sG %in% exclude.sg)

    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_sparsity = mean(sparsity_obs[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_dunn = mean(cell_dunn[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_dunn = mean(gene_dunn[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_connectivity = mean(cell_connectivity[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_connectivity = mean(gene_connectivity[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_silhouette = mean(cell_silhouette[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_silhouette = mean(gene_silhouette[modelname=="dropout"]))

    dir.create(file.path(savepath_master,subdir),recursive=TRUE)
    fn = file.path(savepath_master,subdir,"mean_plot_point.pdf")

    metrics_df0
    for(cn in 1:length(colnames(metrics_df0))){
        if(length(grep("stab",colnames(metrics_df0)[cn]))>0){
            print(paste0("found ",colnames(metrics_df0)[cn]))
            colnames(metrics_df0)[cn] <- stringr::str_sub(colnames(metrics_df0)[cn], start=6L, end=-1L)
            print(paste0("fixed ",colnames(metrics_df0)[cn]))
        }
    }
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_gene_dunn = max(gene_dunn))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_cell_dunn = max(gene_dunn))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(min_gene_connectivity = max(gene_connectivity))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(min_cell_connectivity = max(gene_connectivity))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_gene_silhouette = max(gene_silhouette))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_cell_silhouette = max(cell_silhouette))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_gene_apn = max(gene_apn))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_cell_apn = max(cell_apn))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_gene_ad = max(gene_ad))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_cell_ad = max(cell_ad))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_gene_adm = max(gene_adm))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_cell_adm = max(cell_adm))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_gene_fom = max(gene_fom))
    metrics_df0 %>% group_by(get(sparsityparam), modelname) %>% summarise(max_cell_fom = max(cell_fom))

    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_apn = mean(cell_apn[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_apn = mean(gene_apn[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_ad = mean(cell_ad[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_ad = mean(gene_ad[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_adm = mean(cell_adm[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_adm = mean(gene_adm[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_fom = mean(cell_fom[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_fom = mean(gene_fom[modelname=="dropout"]))
    
    metrics_prefix = c("dunn","connectivity","silhouette","apn","ad","adm","fom")
    # metrics_prefix = c("dunn","connectivity","silhouette")
    for(iii in 1:length(metrics_prefix)){
        meanvar_cell = paste0("mean_dropout_cell_",metrics_prefix[iii])   
        meanvar_gene = paste0("mean_dropout_gene_",metrics_prefix[iii])   
        yvar_cell = paste0("gene_",metrics_prefix[iii])
        yvar_gene = paste0("cell_",metrics_prefix[iii])
        dir.create(file.path(savepath_master,subdir,"scatter"),recursive=TRUE)
        dir.create(file.path(savepath_master,subdir,"violin"),recursive=TRUE)
        dir.create(file.path(savepath_master,subdir,"box"),recursive=TRUE)
        dir.create(file.path(savepath_master,subdir,"bar"),recursive=TRUE)
        print(paste0("plotting scatter prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"scatter",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=wallclock,y=get(yvar_cell),color=modelname)) + 
            geom_point() +
            geom_hline(
                aes(yintercept = get(meanvar_cell), 
                colour = modelname), 
                filter(metrics_df0,modelname=="dropout")) +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            ggtitle(paste0("cell ",metrics_prefix[iii]," vs wallclock"))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting scatter prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"scatter",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=wallclock,y=get(yvar_gene),color=modelname)) + 
            geom_point() +
            geom_hline(
                aes(yintercept = get(meanvar_gene), 
                colour = modelname), 
                filter(metrics_df0,modelname=="dropout")) +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            ggtitle(paste0("gene ",metrics_prefix[iii]," vs wallclock"))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting box prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"box",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelstring,y=get(yvar_cell),color=modelname)) + 
            geom_boxplot() +
            ylab("Dunn Index") +
            xlab("Model") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            theme(axis.text.x = element_text(size = 5, angle = 90, face = "plain")) +
            ggtitle(paste0("cell ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting box prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"box",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelstring,y=get(yvar_gene),color=modelname)) + 
            geom_boxplot() +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            theme(axis.text.x = element_text(size = 5, angle = 90, face = "plain")) +
            ggtitle(paste0("gene ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting violin prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"violin",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelname,y=get(yvar_cell),color=modelname)) + 
            geom_violin() +
            ylab("Dunn Index") +
            xlab("Model") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            ggtitle(paste0("cell ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting violin prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"violin",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelname,y=get(yvar_gene),color=modelname)) + 
            geom_violin() +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            ggtitle(paste0("gene ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting bar prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"bar",paste0("max_",metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelname,y=get(yvar_cell),color=modelname)) + 
            geom_bar(stat="identity") +
            ylab("Dunn Index") +
            xlab("Model") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            ggtitle(paste0("max_cell ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting bar prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"bar",paste0("max_",metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelname,y=get(yvar_gene),color=modelname)) + 
            geom_bar(stat="identity") +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            rotate_x_text(45) +
            facet_grid(rows=vars(mean_dropout_sparsity)) +
            ggtitle(paste0("max_gene ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)
    }

    return(metrics_df0)
}

run_clusterMetrics_final_nested_clvalid <- function(
                        savepath_master,
                        outputmaster,
                        plotwidth=10,
                        plotheight=10,
                        excludemodels=c(),
                        maxwallclock=20,
                        n_samp_cell=800,
                        n_samp_gene=800,
                        subdir="plots",
                        sparsityparam="dmid",
                        k=3,
                        clMethod="hierarchical"){
    # paired_priority: a string match to the first level in a pair eg "NL" for non-lesion, where second level is "DM"
    # min_celltype_rate: eg "0.005" keep celltypes that make up at least 0.5% of all cells
    # typetable_list = 
    # lapply(cellchats_master,function(x){
    #     as.data.frame(t(table(x@idents)/sum(table(x@idents)) >= min_celltype_rate))
    # })
    # typetable = do.call(rbind,typetable_list)
    # typekeep = matrixStats::colProds(as.matrix(typetable))
    # names(typetable)[which(typekeep>0)]


    metrics.orig = list()
    paramsets = list.files(outputmaster,full.names=TRUE,include.dirs=TRUE)
    for(paramind in 1:length(paramsets)){
        dropsets = list.files(paramsets[paramind],full.names=TRUE,include.dirs=TRUE)
        dropsets_short = list.files(paramsets[paramind],include.dirs=TRUE)
        
        fit_inds = grep("output_fit",dropsets)

    
        for(fit_ind in fit_inds){
            dataparams = get_model_params(dropsets_short[fit_ind])
            simdirs = list.files(dropsets[fit_ind],full.names=TRUE)
            for(simdir in simdirs){
                modeldirs = list.files(simdir,full.names=TRUE)
                modeldirs = modeldirs[grep("imputemodel_",modeldirs)]
                imethods = sapply(modeldirs,strsplit,"/") %>% sapply(last) %>% unname
                ind_dropmodel = grep("imputemodel_dropout",modeldirs)
                if(file.exists(file.path(modeldirs[ind_dropmodel],"imputation_loss.csv"))){
                    droprate_backup = read.csv(file.path(modeldirs[ind_dropmodel],"imputation_loss.csv"),row.names=1)$sparsity_obs
                    for(simind in 1:length(imethods)){
                        myparams = get_model_params(imethods[simind])
                        imputedir = modeldirs[simind]
                        imethod = myparams[["imputemodel"]]
                        ldfname = file.path(imputedir,paste0(imethod,"_logdf.csv"))
                        scldfname = file.path(imputedir,"outerStats",paste0(imethod,"_logdf.csv")) # if we are running scrabble and things did not finish
                        impname = file.path(imputedir,"imputed_C.csv")
                        implossname = file.path(imputedir,"imputation_loss.csv")
                        
                        emiter_inds = grep("emIter_",list.files(imputedir,full.names=TRUE))

                        if(file.exists(ldfname) & file.exists(impname) & file.exists(implossname)){
                            print(paste0("processing converged metrics for simind=",simind,":",length(imethods)))
                            impdata = read.csv(impname,row.names=1)
                            implossdata = read.csv(implossname,row.names=1)
                            ldfdata = read.csv(ldfname,row.names=1)
                            n_samp_cell = min(n_samp_cell,ncol(impdata))
                            n_samp_gene = min(n_samp_gene,nrow(impdata))
                            if(length(k) > 1 & length(clMethod) > 1){
                                print(paste0("multiple k/clmethod detected, running optimal clustering"))
                                print(paste0("running optimal clustering on internal"))
                                intern_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), k, clMethods=clMethod, validation="internal",maxitems=10000)                            
                                clmetrics_cell = log(clValid::optimalScores(intern_cell)$Score)
                                names(clmetrics_cell) = rownames(clValid::optimalScores(intern_cell))
                                intern_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), k, clMethods=clMethod, validation="internal",maxitems=10000)                            
                                clmetrics_gene = log(clValid::optimalScores(intern_gene)$Score)
                                names(clmetrics_gene) = rownames(clValid::optimalScores(intern_gene))

                                print(paste0("running optimal clustering on stability"))
                                stab_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), 3, clMethods="kmeans", validation="stability",maxitems=10000)                            
                                clmetrics_cell_stab = log(clValid::optimalScores(stab_cell)$Score)
                                names(clmetrics_cell_stab) = rownames(clValid::optimalScores(stab_cell))
                                stab_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), 3, clMethods="kmeans", validation="stability",maxitems=10000)                            
                                clmetrics_gene_stab = log(clValid::optimalScores(stab_gene)$Score)
                                names(clmetrics_gene_stab) = rownames(clValid::optimalScores(stab_gene))

                            }else{
                                print(paste0("single k/clmethod detected, running specified clustering"))
                                print(paste0("running specified clustering on internal"))
                                intern_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), k, clMethods=clMethod, validation="internal",maxitems=10000)                            
                                clmetrics_cell = intern_cell@measures[,as.character(k),clMethod]
                                intern_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), k, clMethods=clMethod, validation="internal",maxitems=10000)                            
                                clmetrics_gene = intern_gene@measures[,as.character(k),clMethod]

                                print(paste0("running specified clustering on stability"))
                                stab_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), k, clMethods=clMethod, validation="stability",maxitems=10000)                            
                                clmetrics_cell_stab = log(clValid::optimalScores(stab_cell)$Score)
                                names(clmetrics_cell_stab) = rownames(clValid::optimalScores(stab_cell))
                                stab_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), k, clMethods=clMethod, validation="stability",maxitems=10000)                            
                                clmetrics_gene_stab = log(clValid::optimalScores(stab_gene)$Score)
                                names(clmetrics_gene_stab) = rownames(clValid::optimalScores(stab_gene))
                           }

                            names(clmetrics_cell) = paste("cell_",tolower(names(clmetrics_cell)),sep="")
                            names(clmetrics_gene) = paste("gene_",tolower(names(clmetrics_gene)),sep="")

                            names(clmetrics_cell_stab) = paste("cell_",tolower(names(clmetrics_cell_stab)),sep="")
                            names(clmetrics_gene_stab) = paste("gene_",tolower(names(clmetrics_gene_stab)),sep="")

                            cmetrics = cbind(tail(implossdata, n = 1),tail(ldfdata, n = 1),t(clmetrics_cell),t(clmetrics_gene),t(clmetrics_cell_stab),t(clmetrics_gene_stab),as.data.frame(dataparams))
                            cmetrics$simname=imethods[simind]
                            cmetrics$sparsity_true = droprate_backup
                            cmetrics$dropout_rate = NULL

                            # data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(cmetrics)))
                            # cmetrics = cbind(cmetrics,data_cmetrics)
                            metrics.orig[[length(metrics.orig)+1]] = cmetrics
                        # }else if(file.exists(scldfname)){
                        #     print(paste0("processing not-yet-converged scrabble metrics for simind=",simind,":",length(imethods)))
                        #     cmetrics = read.csv(scldfname,row.names=1)
                        #     cmetrics$modelname=imethod
                        #     cmetrics$simname=imethods[simind]
                        #     cmetrics$dropout_rate = droprate_backup

                        #     data_cmetrics = as.data.frame(dataparams) %>% dplyr::slice(rep(1:n(), each = nrow(cmetrics)))
                        #     cmetrics = cbind(cmetrics,data_cmetrics)
                        #     metrics.orig[[length(metrics.orig)+1]] = cmetrics
                            
                        }else{
                            print(paste0("imputation or log output is missing from: ",imputedir))
                        }
                    }
                }
            }
        }
    }
    metrics_df.orig = do.call(plyr::rbind.fill,metrics.orig)
    write.csv(metrics_df.orig,file=file.path(savepath_master,"metrics_df.orig.csv"))

    metrics_df0 = metrics_df.orig %>% filter(!modelname %in% excludemodels & wallclock <= maxwallclock)

    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_dunn = mean(cell_dunn[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_dunn = mean(gene_dunn[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_connectivity = mean(cell_connectivity[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_connectivity = mean(gene_connectivity[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_cell_silhouette = mean(cell_silhouette[modelname=="dropout"]))
    metrics_df0 = metrics_df0 %>% group_by(get(sparsityparam)) %>% mutate(mean_dropout_gene_silhouette = mean(gene_silhouette[modelname=="dropout"]))

    dir.create(file.path(savepath_master,subdir),recursive=TRUE)
    
    metrics_prefix = c("dunn","connectivity","silhouette","apn","ad","adm","fom")
    for(iii in 1:length(metrics_prefix)){
        meanvar_cell = paste0("mean_dropout_cell_",metrics_prefix[iii])   
        meanvar_gene = paste0("mean_dropout_gene_",metrics_prefix[iii])   
        yvar_cell = paste0("gene_",metrics_prefix[iii])
        yvar_gene = paste0("cell_",metrics_prefix[iii])
        dir.create(file.path(savepath_master,subdir,"scatter"),recursive=TRUE)
        dir.create(file.path(savepath_master,subdir,"box"),recursive=TRUE)
        print(paste0("plotting scatter prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"scatter",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=wallclock,y=get(yvar_cell),color=modelname)) + 
            geom_point() +
            geom_hline(
                aes(yintercept = get(meanvar_cell), 
                colour = modelname), 
                filter(metrics_df0,modelname=="dropout")) +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            facet_grid(rows=vars(sparsity_true)) +
            ggtitle(paste0("cell ",metrics_prefix[iii]," vs wallclock"))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting scatter prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"scatter",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=wallclock,y=get(yvar_gene),color=modelname)) + 
            geom_point() +
            geom_hline(
                aes(yintercept = get(meanvar_gene), 
                colour = modelname), 
                filter(metrics_df0,modelname=="dropout")) +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            facet_grid(rows=vars(sparsity_true)) +
            ggtitle(paste0("gene ",metrics_prefix[iii]," vs wallclock"))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting box prefix:",metrics_prefix[iii],", cell"))
        fn = file.path(savepath_master,subdir,"box",paste0(metrics_prefix[iii],"_cell.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelname,y=get(yvar_cell),color=modelname)) + 
            geom_boxplot() +
            ylab("Dunn Index") +
            xlab("Model") +
            facet_grid(rows=vars(sparsity_true)) +
            ggtitle(paste0("cell ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

        print(paste0("plotting box prefix:",metrics_prefix[iii],", gene"))
        fn = file.path(savepath_master,subdir,"box",paste0(metrics_prefix[iii],"_gene.pdf"))
        p = ggplot(
            data=metrics_df0,
            aes(x=modelname,y=get(yvar_gene),color=modelname)) + 
            geom_boxplot() +
            ylab("Dunn Index") +
            xlab("Wallclock") +
            facet_grid(rows=vars(sparsity_true)) +
            ggtitle(paste0("gene ",metrics_prefix[iii]))
            # ggsave(plot=p,file=file.path(savepath_master,paste0("ds_",dsname,",nIn_",nIn,",nOut_",nOut,",nSDC_",nSDC,"_conv.pdf")),width=plotwidth,height=plotheight)        
            ggsave(plot=p,file=fn,width=plotwidth,height=plotheight)

    }
    return(metrics_df0)
}

run_clvalid <- function(
                        impdata,
                        n_samp_cell=800,
                        n_samp_gene=800,
                        k=8,
                        clMethod="kmeans",
                        runstability=FALSE){
    n_samp_cell = min(n_samp_cell,ncol(impdata))
    n_samp_gene = min(n_samp_gene,nrow(impdata))
    
    if(length(k) > 1 & length(clMethod) > 1){
        print(paste0("multiple k/clmethod detected"))
        print(paste0("running optimal clustering on internal"))
        intern_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), k, clMethods=clMethod, validation="internal",maxitems=10000)
        clmetrics_cell = clValid::optimalScores(intern_cell)$Score
        names(clmetrics_cell) = tolower(paste0("cell_",rownames(clValid::optimalScores(intern_cell))))

        intern_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), k, clMethods=clMethod, validation="internal",maxitems=10000)
        clmetrics_gene = clValid::optimalScores(intern_gene)$Score
        names(clmetrics_gene) = tolower(paste0("gene_",rownames(clValid::optimalScores(intern_gene))))
        if(runstability){
            print(paste0("running optimal clustering on stability"))
            stab_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), 3, clMethods="kmeans", validation="stability",maxitems=10000)
            clmetrics_cell_stab = clValid::optimalScores(stab_cell)$Score
            names(clmetrics_cell_stab) = tolower(paste0("stab_cell_",rownames(clValid::optimalScores(stab_cell))))

            stab_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), 3, clMethods="kmeans", validation="stability",maxitems=10000)
            clmetrics_gene_stab = clValid::optimalScores(stab_gene)$Score
            names(clmetrics_gene_stab) = tolower(paste0("stab_gene_",rownames(clValid::optimalScores(stab_gene))))
        }

    }else{
        print(paste0("single k/clmethod detected"))
        print(paste0("running specified clustering on internal"))
        intern_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), k, clMethods=clMethod, validation="internal",maxitems=10000)
        clmetrics_cell = intern_cell@measures[,as.character(k),clMethod]
        names(clmetrics_cell) = tolower(paste0("cell_",names(clmetrics_cell)))
        
        intern_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), k, clMethods=clMethod, validation="internal",maxitems=10000)
        clmetrics_gene = intern_gene@measures[,as.character(k),clMethod]
        names(clmetrics_gene) = tolower(paste0("gene_",names(clmetrics_gene)))

        if(runstability){
            print(paste0("running specified clustering on stability"))
            stab_cell <- clValid::clValid(as.data.frame(t(impdata)) %>% sample_n(n_samp_cell), k, clMethods=clMethod, validation="stability",maxitems=10000)
            clmetrics_cell_stab = stab_cell@measures[,as.character(k),clMethod]
            names(clmetrics_cell_stab) = tolower(paste0("stab_cell_",names(clmetrics_cell_stab)))

            stab_gene <- clValid::clValid(as.data.frame(impdata) %>% sample_n(n_samp_gene), k, clMethods=clMethod, validation="stability",maxitems=10000)
            clmetrics_gene_stab = stab_gene@measures[,as.character(k),clMethod]
            names(clmetrics_gene_stab) = tolower(paste0("stab_gene_",names(clmetrics_gene_stab)))
        }
    }
    if(runstability){
        paste0("returning internal and stability metrics")
        return(c(clmetrics_cell,clmetrics_gene,clmetrics_cell_stab,clmetrics_gene_stab))
    }else{
        paste0("returning internal metrics")
        # return(cbind(t(clmetrics_cell),t(clmetrics_gene)))
        return(c(clmetrics_cell,clmetrics_gene))
    }
}

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
