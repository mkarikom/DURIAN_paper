run_signaling_multiset <- function(
                        savepath_master,
                        outputmaster,
                        sourcepath,
                        CellChatDB,
                        prefixes,
                        paired_priority=NULL,
                        common,
                        min_celltype_rate=0.005){

    # paired_priority: a string match to the first level in a pair eg "NL" for non-lesion, where second level is "DM"
    # min_celltype_rate: eg "0.005" keep celltypes that make up at least 0.5% of all cells

    typetable_list = list()
    for(prefixind in 1:length(prefixes)){
        prefix = prefixes[[prefixind]]
        pDataC = read.csv(file.path(sourcepath,paste0(prefix[1],"_pDataC.csv")),row.names=1)
        typetable_list[[prefixind]] = as.data.frame(t(table(pDataC$cellType)/sum(table(pDataC$cellType)) >= min_celltype_rate))
    }
    typetable0 = do.call(rbind,typetable_list)
    typematrix0 = matrixStats::colProds(as.matrix(typetable0))
    typekeep0 = names(typetable0)[which(typematrix0>0)]


    methodlist0 = c()
    for(prefixind in 1:length(prefixes)){
        prefix = prefixes[[prefixind]]

        outputdirs = list.files(outputmaster,full.names=TRUE,include.dirs=TRUE)
        fileind = grep(prefix[1],outputdirs)

        modeldirs = list.files(file.path(outputdirs[fileind],"output_fit"))
        modeldirs[grep("imputemodel",modeldirs)]
        methodlist0 = c(methodlist0,modeldirs[grep("imputemodel",modeldirs)])
    }
    # methodlist = unique(methodlist0) %>% strsplit(",") %>% sapply( "[", 1 )
    methodlist = unique(methodlist0)
    # methodlist = methodlist[-grep("dropout",methodlist)]
    methodlist = methodlist[grep("DURIAN.MuSiC",methodlist)]

    for(imethod in methodlist){
        print(paste0("running imethod:",imethod))
        cellchats_master.orig = list()
        for(prefixind in 1:length(prefixes)){
            prefix = prefixes[[prefixind]]

            outputdirs = list.files(outputmaster,full.names=TRUE,include.dirs=TRUE)
            fileind = grep(prefix[1],outputdirs)

            imputedir = file.path(outputdirs[fileind],"output_fit",imethod)

            imputedC = read.csv(file.path(imputedir,"imputed_C.csv"),row.names=1)
            origC = read.csv(file.path(sourcepath,paste0(prefix[1],"_C.csv")),row.names=1)
            pDataC = read.csv(file.path(sourcepath,paste0(prefix[1],"_pDataC.csv")),row.names=1)

            pDataC = filter(pDataC,cellType %in% typekeep0)
            common_cells = intersect(colnames(imputedC),rownames(pDataC))

            pDataC = pDataC[common_cells,]
            origC = origC[,common_cells]
            imputedC = imputedC[,common_cells]

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

            cellchat.orig <- netAnalysis_computeCentrality(cellchat.orig, slot.name = "netP")
            cellchat.imputed <- netAnalysis_computeCentrality(cellchat.imputed, slot.name = "netP")

            cellchats_master.orig[[paste0(prefix[1],".orig")]] = cellchat.orig
            cellchats_master.orig[[paste0(prefix[1],".imputed")]] = cellchat.imputed
        }
        if(length(prefixes)==1){
            print("single prefix detected")
            # make sure "orig" is first and "imputed" is second in all (eg differential) comparisons
            sortingdf = data.frame(
                name=names(cellchats_master.orig),
                imputed=stringr::str_detect(names(cellchats_master.orig),"orig"))
            cellchats_master = cellchats_master.orig[arrange(sortingdf,desc(imputed))[,"name"]]

            cellchats_pair = list()
            cellchats_pair[[paste0("orig:imputed,",prefixes[[1]][1])]] = cellchats_master
            cellchats_merged = list()
            cellchats_merged[[1]] = mergeCellChat(cellchats_pair[[1]], add.names = names(cellchats_pair[[1]]))
            names(cellchats_merged) = names(cellchats_pair)
        }else if(length(prefixes)==2){
            print("two prefixes detected")
            # when applicable make sure "orig" is first and "imputed" is second in all (eg differential) comparisons
            # when applicable make sure data name containing `paired_priority` is first in all (eg differential) comparisons
            sortingdf = data.frame(
                name=names(cellchats_master.orig),
                pair=stringr::str_detect(names(cellchats_master.orig),paired_priority),
                imputed=stringr::str_detect(names(cellchats_master.orig),"orig"))
            cellchats_master = cellchats_master.orig[arrange(sortingdf,desc(pair),desc(imputed))[,"name"]]

            name_priority = unlist(lapply(strsplit(names(cellchats_master)[grep("imputed",names(cellchats_master))],common),function(x){x[1]}))

            name_priority = unlist(lapply(strsplit(names(cellchats_master)[grep("orig",names(cellchats_master))],".orig"),function(x){x[1]}))

            cellchats_pair = list()
            cellchats_pair[[paste0("orig:imputed,",prefixes[[1]][1])]] = cellchats_master[grep(prefixes[[1]][1],names(cellchats_master))]
            cellchats_pair[[paste0("orig:imputed,",prefixes[[2]][1])]] = cellchats_master[grep(prefixes[[2]][1],names(cellchats_master))]
            cellchats_pair[[paste0("orig,",name_priority[1],":",name_priority[2])]] = cellchats_master[grep("orig",names(cellchats_master))]
            cellchats_pair[[paste0("imputed,",name_priority[1],":",name_priority[2])]] = cellchats_master[grep("imputed",names(cellchats_master))]
            cellchats_merged = list()
            cellchats_merged[[1]] = mergeCellChat(cellchats_pair[[1]], add.names = names(cellchats_pair[[1]]))
            cellchats_merged[[2]] = mergeCellChat(cellchats_pair[[2]], add.names = names(cellchats_pair[[2]]))
            cellchats_merged[[3]] = mergeCellChat(cellchats_pair[[3]], add.names = names(cellchats_pair[[3]]))
            cellchats_merged[[4]] = mergeCellChat(cellchats_pair[[4]], add.names = names(cellchats_pair[[4]]))
            names(cellchats_merged) = names(cellchats_pair)
        }else{
            print("error: multiple prefixes detected")
            return()
        }

        #############################################################################
        #
        # Plotting
        #
        #############################################################################


        typetable_list = 
        lapply(cellchats_master,function(x){
            as.data.frame(t(table(x@idents)/sum(table(x@idents)) >= min_celltype_rate))
        })
        typetable = do.call(rbind,typetable_list)
        typematrix = matrixStats::colProds(as.matrix(typetable))
        typekeep = names(typetable)[which(typematrix>0)]
        
        for(i in 1:length(cellchats_merged)){
            savepath = file.path(savepath_master,paste0("analysis_cellchat,given_imputemodel,",imethod,",compare_",names(cellchats_merged)[i]))
            # dir.create(file.path(savepath,"pData"))
            # dir.create(file.path(savepath,"rdsData"))
            # dir.create(file.path(savepath,"plots","heatmap"),recursive = TRUE)
            # dir.create(file.path(savepath,"plots","circle"),recursive = TRUE)
            # dir.create(file.path(savepath,"plots","violin"),recursive = TRUE)
            # dir.create(file.path(savepath,"plots","bubble"))
            # dir.create(file.path(savepath,"plots","single_path_source_groups"))
            # dir.create(file.path(savepath,"plots","centrality_heatmap"))
            cellchats = cellchats_pair[[i]]
            merged = cellchats_merged[[i]]

            # single-celltype target
            groupSize <- as.numeric(table(cellchats[[1]]@idents))
            
            ########################################################################
            ########################################################################
            # enumerate inter-gene for pathways discovered in second data set
            ########################################################################
            ########################################################################
            gained_pathways = setdiff(cellchats[[2]]@netP$pathway,cellchats[[1]]@netP$pathway)
            if(length(gained_pathways)>0){
                for(j in 1:length(gained_pathways)){
                    dir.create(file.path(savepath,"plots","single_pathway.new",gained_pathways[j]),recursive=TRUE)
                    fname_p2 = file.path(savepath,"plots","single_pathway.new",gained_pathways[j],paste0(names(cellchats)[2],".pdf"))
                    pdf(file=fname_p2,width=6,height=3.5)
                        netVisual_chord_gene(cellchats[[2]], signaling = gained_pathways[j],title.name=paste0(gained_pathways[j],"\n",names(cellchats)[2]))
                    dev.off()
                }
            }
            ########################################################################
            ########################################################################
            # enumerate inter-gene for pathways common to first and second data set
            ########################################################################
            ########################################################################
            shared_pathways = intersect(cellchats[[2]]@netP$pathway,cellchats[[1]]@netP$pathway)
            if(length(shared_pathways)>0){
                for(j in 1:length(shared_pathways)){
                    dir.create(file.path(savepath,"plots","single_pathway.shared",shared_pathways[j]),recursive=TRUE)
                    fname_p1 = file.path(savepath,"plots","single_pathway.shared",shared_pathways[j],paste0(names(cellchats)[1],".pdf"))
                    fname_p2 = file.path(savepath,"plots","single_pathway.shared",shared_pathways[j],paste0(names(cellchats)[2],".pdf"))
                    pdf(file=fname_p1,width=6,height=3.5)
                        netVisual_chord_gene(cellchats[[1]], signaling = shared_pathways[j],title.name=paste0(shared_pathways[j],"\n",names(cellchats)[1]))
                    dev.off()
                    pdf(file=fname_p2,width=6,height=3.5)
                        netVisual_chord_gene(cellchats[[2]], signaling = shared_pathways[j],title.name=paste0(shared_pathways[j],"\n",names(cellchats)[2]))
                    dev.off()
                }
            }

            ########################################################################
            ########################################################################
            # enumerate inter-celltype by weight
            ########################################################################
            ########################################################################
            mat_p1 = cellchats[[1]]@net$weight # priority data set
            mat_p2 = cellchats[[2]]@net$weight # secondary data set

            ########################################################################
            # single target
            ########################################################################
            # rows are sources, columns are targets
            for (j in 1:ncol(mat_p1)) {
                pathdiff = all(as.logical(mat_p1[,j])==as.logical(mat_p2[,j]))
                if(!pathdiff){
                    dir.create(file.path(savepath,"plots","single_target.weight"),recursive=TRUE)
                    mat_p1_target <- matrix(0, nrow = nrow(mat_p1), ncol = ncol(mat_p1), dimnames = dimnames(mat_p1))
                    mat_p2_target <- matrix(0, nrow = nrow(mat_p1), ncol = ncol(mat_p1), dimnames = dimnames(mat_p1))                
                    mat_p1_target[,j] <- mat_p1[,j]
                    mat_p2_target[,j] <- mat_p2[,j]
                    dir.create(file.path(savepath,"plots","single_target.weight",colnames(mat_p1)[j]),recursive=TRUE)
                    dir.create(file.path(savepath,"plots","single_target.weight",colnames(mat_p2)[j]),recursive=TRUE)
                    fname_p1 = file.path(savepath,"plots","single_target.weight",colnames(mat_p1)[j],paste0(names(cellchats)[1],".pdf"))
                    fname_p2 = file.path(savepath,"plots","single_target.weight",colnames(mat_p2)[j],paste0(names(cellchats)[2],".pdf"))
                    pdf(file=fname_p1,width=6,height=3.5)
                        netVisual_circle(mat_p1_target, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_p1_target), title.name = colnames(mat_p1_target)[j])
                    dev.off()
                    pdf(file=fname_p2,width=6,height=3.5)
                        netVisual_circle(mat_p2_target, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_p2_target), title.name = colnames(mat_p2_target)[j])
                    dev.off()
                }
            }
            ########################################################################
            # single source
            ########################################################################
            # rows are sources, columns are targets
            for (j in 1:nrow(mat_p1)) {
                pathdiff = all(as.logical(mat_p1[j,])==as.logical(mat_p2[j,]))
                if(!pathdiff){
                    dir.create(file.path(savepath,"plots","single_source.weight"),recursive=TRUE)
                    mat_p1_target <- matrix(0, nrow = nrow(mat_p1), ncol = ncol(mat_p1), dimnames = dimnames(mat_p1))
                    mat_p2_target <- matrix(0, nrow = nrow(mat_p1), ncol = ncol(mat_p1), dimnames = dimnames(mat_p1))                
                    mat_p1_target[j,] <- mat_p1[j,]
                    mat_p2_target[j,] <- mat_p2[j,]
                    dir.create(file.path(savepath,"plots","single_source.weight",rownames(mat_p1)[j]),recursive=TRUE)
                    dir.create(file.path(savepath,"plots","single_source.weight",rownames(mat_p2)[j]),recursive=TRUE)
                    fname_p1 = file.path(savepath,"plots","single_source.weight",rownames(mat_p1)[j],paste0(names(cellchats)[1],".pdf"))
                    fname_p2 = file.path(savepath,"plots","single_source.weight",rownames(mat_p2)[j],paste0(names(cellchats)[2],".pdf"))
                    pdf(file=fname_p1,width=6,height=3.5)
                        netVisual_circle(mat_p1_target, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_p1_target), title.name = rownames(mat_p1_target)[j])
                    dev.off()
                    pdf(file=fname_p2,width=6,height=3.5)
                        netVisual_circle(mat_p2_target, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_p2_target), title.name = rownames(mat_p2_target)[j])
                    dev.off()
                }
            }


            ########################################################################
            ########################################################################
            # differential expression
            ########################################################################
            ########################################################################
            if(sum(as.logical(cellchats[[2]]@net$count) - as.logical(cellchats[[1]]@net$count)) > 0){
                dir.create(file.path(savepath,"plots","differential_circle"),recursive=TRUE)
                weight.max <- getMaxWeight(cellchats, attribute = c("idents","count"))
                fname = file.path(savepath,"plots","differential_circle","count_orig.pdf")
                pdf(file=fname)
                netVisual_circle(
                    cellchats[[1]]@net$count, 
                    weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                    title.name = paste0("Number of interactions \n", names(cellchats)[1]))
                dev.off()
                fname = file.path(savepath,"plots","differential_circle","count_imputed.pdf")
                pdf(file=fname)
                netVisual_circle(
                    cellchats[[2]]@net$count, 
                    weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                    title.name = paste0("Number of interactions \n", names(cellchats)[2]))
                dev.off()

                weight.max <- getMaxWeight(cellchats, attribute = c("idents","count"))
                fname = file.path(savepath,"plots","differential_circle","weight_orig.pdf")
                pdf(file=fname)
                netVisual_circle(
                    cellchats[[1]]@net$weight, 
                    weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                    title.name = paste0("Strength of interactions \n", names(cellchats)[1]))
                dev.off()
                fname = file.path(savepath,"plots","differential_circle","weight_imputed.pdf")
                pdf(file=fname)
                netVisual_circle(
                    cellchats[[2]]@net$weight, 
                    weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                    title.name = paste0("Strength of interactions \n", names(cellchats)[2]))
                dev.off()
            }

            dir.create(file.path(savepath,"plots","differential_heatmap"),recursive=TRUE)
            fname = file.path(savepath,"plots","differential_heatmap","diff_counts.pdf")
            pdf(file=fname)
                p = netVisual_heatmap(cellchats_merged[[i]],font.size = 15,font.size.title = 20)
                draw(p)
            dev.off()
            fname = file.path(savepath,"plots","differential_heatmap","diff_weights.pdf")
            pdf(file=fname,height=5,width=5)
                p = netVisual_heatmap(cellchats_merged[[i]],measure="weight",font.size = 15,font.size.title = 20)
                draw(p)
            dev.off()
        }
    }
}

netAnalysis_signalingRole_network_retlist <- function(object, signaling, slot.name = "netP", measure = c("outdeg","indeg","flowbet","info"), measure.name = c("Sender","Receiver","Mediator","Influencer"),
                                    color.use = NULL, color.heatmap = "BuGn",
                                    width = 6.5, height = 1.4, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE) {
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  for(i in 1:length(centr)) {
    centr0 <- centr[[i]]
    mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
    if (!is.null(measure)) {
      mat <- mat[measure,]
      if (!is.null(measure.name)) {
        rownames(mat) <- measure.name
      }
    }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)

    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

    df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "column",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Importance",
                  bottom_annotation = col_annotation,
                  cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = paste0(names(centr[i]), " signaling pathway network"),column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 45,
                  heatmap_legend_param = list(title = "Importance", title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, at = c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1)),
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
    )
    draw(ht1)
  }
  list(obj=mat,hm=ht1)
}