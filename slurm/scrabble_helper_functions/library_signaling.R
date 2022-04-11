run_signaling <- function(
                        savepath,
                        imputedC,
                        pDataC,
                        origC,
                        CellChatDB){
    genesIntersect = intersect(rownames(origC),rownames(imputedC))
    cells_imputed = colnames(imputedC)
    #############################################################################
    # select disease subsets and normalize
    #############################################################################
    data.orig = normalizeData(as.matrix(origC[genesIntersect,cells_imputed]))
    data.imputed = normalizeData(as.matrix(imputedC[genesIntersect,cells_imputed]))
    colnames(data.imputed) = paste0("imputed_",colnames(data.imputed))

    meta.orig = pDataC[,c("cellType","sampleID")]
    rownames(meta.orig) = pDataC$cellID
    meta.orig = meta.orig[cells_imputed,]

    meta.imputed = pDataC[,c("cellType","sampleID")]
    rownames(meta.imputed) = paste0("imputed_",pDataC$cellID)
    meta.imputed = meta.imputed[paste0("imputed_",cells_imputed),]
    #############################################################################
    # remove the data wont use anymore (cellchat needs lots of ram)
    #############################################################################
    rm(list=c("origC","imputedC"))
    gc()

    print("creating cellchat object")
    cellchat.orig <- createCellChat(object = data.orig, meta = meta.orig, group.by = "cellType")
    cellchat.imputed <- createCellChat(object = data.imputed, meta = meta.imputed, group.by = "cellType")
    print("cellchat object successfully created")
    #############################################################################
    # remove the data wont use anymore (cellchat needs lots of ram)
    #############################################################################
    rm(list=c("data.orig","data.imputed"))
    gc()

    #############################################################################
    # set ligrec database
    #############################################################################

    cellchat.orig@DB <- CellChatDB
    cellchat.imputed@DB <- CellChatDB

    #############################################################################
    # preprocess the data
    #############################################################################
    cellchat.orig <- subsetData(cellchat.orig)
    cellchat.orig <- identifyOverExpressedGenes(cellchat.orig)
    cellchat.orig <- identifyOverExpressedInteractions(cellchat.orig)

    cellchat.imputed <- subsetData(cellchat.imputed)
    cellchat.imputed <- identifyOverExpressedGenes(cellchat.imputed)
    cellchat.imputed <- identifyOverExpressedInteractions(cellchat.imputed)

    # cellchat.orig <- projectData(cellchat.orig, PPI.mouse)
    # cellchat.imputed <- projectData(cellchat.imputed, PPI.mouse)
    gc()
    #############################################################################
    # compute probs without preprocessing
    #############################################################################
    cellchat.orig <- computeCommunProb(cellchat.orig,raw.use = TRUE)
    cellchat.imputed <- computeCommunProb(cellchat.imputed,raw.use = TRUE)


    #############################################################################
    # compute pathway probs
    #############################################################################
    cellchat.orig <- computeCommunProbPathway(cellchat.orig)
    cellchat.imputed <- computeCommunProbPathway(cellchat.imputed)

    cellchat.orig <- aggregateNet(cellchat.orig)
    cellchat.imputed <- aggregateNet(cellchat.imputed)

    #############################################################################
    #
    # Plotting
    #
    #############################################################################

    dir.create(file.path(savepath,"pData"))
    dir.create(file.path(savepath,"rdsData"))
    dir.create(file.path(savepath,"plots","heatmap"),recursive = TRUE)
    dir.create(file.path(savepath,"plots","circle"),recursive = TRUE)
    dir.create(file.path(savepath,"plots","violin"),recursive = TRUE)
    dir.create(file.path(savepath,"plots","bubble"))
    dir.create(file.path(savepath,"plots","aggregate_path_single_source"))
    dir.create(file.path(savepath,"plots","single_path_source_groups"))
    dir.create(file.path(savepath,"plots","centrality_heatmap"))

    cellchats = list(
    wt.orig=cellchat.orig,
    wt.imputed=cellchat.imputed
    )

    saveRDS(cellchats, file=file.path(savepath,"rdsData","cellchats.RDS"))

    #############################################################################
    # remove the data wont use anymore (cellchat needs lots of ram)
    #############################################################################
    rm(list=c("cellchat.orig","cellchat.imputed"))
    gc()
    
    # #############################################################################
    # # quantitative network-centrality based analysis: imputed vs orig
    # #############################################################################

    # change in centrality-based scores
    cellchats$wt.orig <- netAnalysis_computeCentrality(cellchats$wt.orig, slot.name = "netP")
    cellchats$wt.imputed <- netAnalysis_computeCentrality(cellchats$wt.imputed, slot.name = "netP")

    # centrality-based population roles, imputed vs orig 
    pathnames = cellchats$wt.orig@netP$pathways
    for(j in 1:length(pathnames)){
    # find out if sig interactions are available
    orig_sig = !any(is.na(match(pathnames[j],names(slot(cellchats$wt.orig,"netP")$centr))))
    imputed_sig = !any(is.na(match(pathnames[j],names(slot(cellchats$wt.imputed,"netP")$centr))))
    # plot any significant interactions
    if(orig_sig){
        print(paste0("detected significant interactions of orig: ",pathnames[j]))
        netA.orig = netAnalysis_signalingRole_network(cellchats$wt.orig, signaling = c(pathnames[j]));
        orig_ok = !any(is.nan(netA.orig$obj)) && !is.null(netA.orig$obj) # sometimes netAnalysis_computeCentrality generates NaN
        if(orig_ok){
        fn = file.path(savepath,"plots","centrality_heatmap",paste0("roles_",pathnames[j],"_orig.pdf"))
        print(paste0("plotting orig: ",pathnames[j]," to: ",fn))
        pdf(file=fn, width=9,height=3.5)
            pp = ComplexHeatmap::Heatmap(netA.orig$obj);
            draw(pp)
        dev.off()
        }else{
        print(paste0("NaN detected in orig: ",pathnames[j]))
        }
    }
    if(imputed_sig){
        print(paste0("detected significant interactions of imputed: ",pathnames[j]))
        netA.imputed = netAnalysis_signalingRole_network(cellchats$wt.imputed, signaling = c(pathnames[j]));
        imputed_ok = !any(is.nan(netA.imputed$obj)) && !is.null(netA.imputed$obj) # sometimes netAnalysis_computeCentrality generates NaN
        if(imputed_ok){
        fn = file.path(savepath,"plots","centrality_heatmap",paste0("roles_",pathnames[j],"_imputed.pdf"))
        print(paste0("plotting imputed: ",pathnames[j]," to: ",fn))
        pdf(file=fn, width=9,height=3.5)
            pp = ComplexHeatmap::Heatmap(netA.imputed$obj);
            draw(pp)
        dev.off()    
        }else{
        print(paste0("NaN detected in imputed: ",pathnames[j]))
        }
    }
    
    if(orig_sig && imputed_sig && orig_ok && imputed_ok){
        pdf(file=file.path(savepath,"plots","centrality_heatmap",paste0("roles_",pathnames[j],"_compare.pdf")), width=9,height=3.5)
        pp = ComplexHeatmap::Heatmap(netA.orig$obj) + ComplexHeatmap::Heatmap(netA.imputed$obj) + ComplexHeatmap::Heatmap(netA.imputed$obj-netA.orig$obj);
        draw(pp)
        dev.off()
    }
    }


    #############################################################################
    # generate circle plots and heatmaps
    #############################################################################

    saveRDS(cellchats, file=file.path(savepath,"rdsData","cellchats.RDS"))

    for(i in 1:length(cellchats)){
    status=strsplit(names(cellchats)[i],'\\.')[[1]][1]
    impute=strsplit(names(cellchats)[i],'\\.')[[1]][2]

    print(paste0("status: ",status,", impute: pathway, single source"))  
    groupSize <- as.numeric(table(cellchats[[i]]@idents))

    print("aggregate pathway, single source")
    mat = cellchats[[i]]@net$weight
    for (j in 1:nrow(mat)) {
        print(paste0("single pop ",rownames(mat)[j]))
        # print(rownames(mat)[i])
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[j, ] <- mat[j, ]
        # netvis circle
        fname = file.path(savepath,"plots","aggregate_path_single_source",
            paste("adjacency","circle",rownames(mat)[j],impute,status,"pdf",collapse="",sep="."))
        fname = gsub(" ", "", fname, fixed = TRUE)
        print(paste0("filename=",fname))
        pdf(file=fname,
            width=6,height=3.5)
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j])
        dev.off()
    }

    pathnames = cellchats[[i]]@netP$pathways
    for(j in 1:length(pathnames)){
        write.csv(cellchats[[i]]@netP$prob[,,j],
            file=file.path(savepath,"pData",
                            paste(status,impute,pathnames[j],"csv",collapse="",sep=".")))

        print(paste0("pathway ",pathnames[j]))

        print(paste(status,impute,pathnames[j]))    
        res <- extractEnrichedLR(cellchats[[i]], signaling = pathnames[j], 
                                geneLR.return = TRUE, enriched.only = TRUE)
        
        ngenes.plot = length(res$geneLR)
        ntypes.plot = length(unique(cellchats[[i]]@idents))
        hviolin = 20*ngenes.plot + 60
        wviolin = 10*ntypes.plot + 40
        hheatm = 10*ntypes.plot + 40
        wheatm = 10*ntypes.plot + 60
        hcirc = .5*ntypes.plot
        wcirc = .5*ntypes.plot
        
        print("adjacency circle")
        pdf(file=file.path(savepath,"plots","circle",
                        paste("adjacency","circle",status,impute,pathnames[j],"pdf",collapse="",sep=".")),width=wcirc,height=hcirc)
        netVisual_aggregate(cellchats[[i]], signaling = pathnames[j], layout="circle");
        dev.off()
        
        print("adjacency heatmap")
        melted = melt(cellchats[[i]]@netP$prob[,,j])
        colnames(melted) = c("sending","receiving","magnitude")
        gg.heatmap = ggplot(melted, aes(sending, receiving, fill= magnitude)) + 
        geom_tile() +
        theme(axis.text.x = 
                element_text(angle = -60,hjust = 0))
        ggsave(
        file.path(savepath,"plots","heatmap",
                    paste("adjacency","heatmap",status,impute,pathnames[j],"pdf",collapse="",sep=".")),
        plot = gg.heatmap,
        device = NULL,
        path = NULL,
        scale = 1,
        width = wheatm,
        height = hheatm,
        units = "mm",
        dpi = 300,
        limitsize = TRUE)

        print("violin plots single gene")    
        gg.violin = plotGeneExpression(cellchats[[i]], signaling = pathnames[j])
        ggsave(
        file.path(savepath,"plots","violin",
                    paste("expression","violin",status,impute,pathnames[j],"pdf",collapse="",sep=".")),
        plot = gg.violin,
        device = NULL,
        path = NULL,
        # scale = 1,
        width = wviolin,
        height = hviolin,
        units = "mm",
        dpi = 300,
        limitsize = TRUE)

        print("bubble plots single pathway")
        netVisual_bubble(cellchats[[i]], signaling = pathnames[j], remove.isolate = FALSE)
        ggsave(
        file.path(savepath,"plots","bubble",
                    paste("bubble",status,impute,pathnames[j],"pdf",collapse="",sep=".")))    

        print("single pathway, source groups")
        target_groups = list(c("FIB-A","FIB-P"))
        for(tg in 1:length(target_groups)){
        tg_ind = match(target_groups[[tg]],rownames(cellchats[[1]]@netP[[2]]))
        tg_ok = tg_ind[!is.na(tg_ind)]
        if(length(tg_ok)>0){
            dir.create(file.path(savepath,"plots","single_path_source_groups",paste(target_groups[[tg]],collapse="|")))
            fname = file.path(savepath,"plots","single_path_source_groups",paste(target_groups[[tg]],collapse="|"),
                        paste("single_path_source_groups",status,impute,pathnames[j],"pdf",collapse="",sep="."))
            pdf(file=fname,
                width=6,height=3.5)
            netVisual_aggregate(cellchats[[i]], 
                                targets.use = target_groups[[tg]], 
                                signaling = pathnames[j],
                                layout = "circle")
            dev.off()
        }
        }
    }
    }

    saveRDS(cellchats, file=file.path(savepath,"rdsData","cellchats_post.RDS"))

}

