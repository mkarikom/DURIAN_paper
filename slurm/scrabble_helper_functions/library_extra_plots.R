logratio_plot <- function(
    imputedC,pDataC,trueC,
    mthresh=c(-2,2),ylims=c(-7,7),plottitle){
    A = 0.5*rowMeans(log2(trueC+1) + log2(imputedC+1))
    M = rowMeans(log2(trueC+1) - log2(imputedC+1))
    mrate = abs(M/ylims[1])
    plotdf = data.frame(A=A,M=M)
    inbounds=c()
    for(i in 1:length(M)){
        if(M[i] <= mthresh[1]){
            inbounds = c(inbounds,2)
        }else if(M[i] >= mthresh[2]){
            inbounds = c(inbounds,3)
        }else{
            inbounds = c(inbounds,1)
        }
    }
    plotdf$inbounds = as.factor(inbounds)
    p_scr=ggplot(plotdf,aes(x=A,y=M,color=inbounds,alpha=mrate))+
            geom_point()+
            geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5)+
            geom_hline(yintercept=2, linetype="dashed",color = "red", size=1)+
            geom_hline(yintercept=-2, linetype="dashed",color = "blue", size=1)+
            ylim(ylims)+
            scale_colour_manual(values = c("grey30", "blue", "red"))+
            theme_bw()+
            theme(
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(size=22,face="bold"),
                axis.title=element_text(size=22,face="bold"))+
            guides(color = "none", alpha = "none")+
            ggtitle(plottitle)+
            xlab("Average counts")+
            ylab("Log ratio")
    p_scr
}

umap_plot <- function(dat,meta,plottitle,col="cellType",ptsize=5){
  umap <- umap(t(dat))

#   df <- data.frame(x = umap$layout[,1],
#                     y = umap$layout[,2],
#                     cellType = as.factor(meta[colnames(dat),"cellType"]),
#                     sampleID = as.factor(meta[colnames(dat),"sampleID"]))

    df = cbind(umap$layout,meta)
    colnames(df)=c("u1","u2",colnames(meta))

    ggplot(df, aes(u1, u2, colour = get(col))) +
    geom_point(size=ptsize) +
    ggtitle(plottitle)+
    xlab("U1")+
    ylab("U2")+
    theme_bw()+
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=22,face="bold"),
        axis.title=element_text(size=22,face="bold"))
}

umap_alpha_plot <- function(dat,meta,plottitle,col="cellType",inverted=FALSE,ptsize=5){
  umap <- umap(t(dat))


    df = cbind(umap$layout,meta)
    colnames(df)=c("u1","u2",colnames(meta))


    if(inverted){
        p = ggplot(df, aes(u1, u2, colour = get(col),alpha=1/get(col))) +
        geom_point(size=ptsize) +
        ggtitle(plottitle)+
        xlab("U1")+
        ylab("U2")+
        theme_bw()+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=22,face="bold"),
            axis.title=element_text(size=22,face="bold"))
    }else{
        p = ggplot(df, aes(u1, u2, colour = get(col),alpha=get(col))) +
        geom_point(size=ptsize) +
        ggtitle(plottitle)+
        xlab("U1")+
        ylab("U2")+
        theme_bw()+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=22,face="bold"),
            axis.title=element_text(size=22,face="bold"))
    }
}

umap_viridis_plot <- function(dat,meta,plottitle,col="cellType",inverted=FALSE,ptsize=5){
  umap <- umap(t(dat))


    df = cbind(umap$layout,meta)
    colnames(df)=c("u1","u2",colnames(meta))


    if(inverted){
        p = ggplot(df, aes(u1, u2, colour = get(col))) +
        geom_point(size=ptsize) +
        ggtitle(plottitle)+
        xlab("U1")+
        ylab("U2")+
        theme_bw()+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=22,face="bold"),
            axis.title=element_text(size=22,face="bold")) +
        scale_fill_viridis()
    }else{
        p = ggplot(df, aes(u1, u2, colour = get(col))) +
        geom_point(size=ptsize) +
        ggtitle(plottitle)+
        xlab("U1")+
        ylab("U2")+
        theme_bw()+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=22,face="bold"),
            axis.title=element_text(size=22,face="bold"))+
        scale_color_viridis()
    }
}


tsne_plot <- function(dat,meta,plottitle,col="cellType",ptsize=5){
    
    # browser()
    rownames(meta) = meta$cellID

    perplexity = max(1,floor((ncol(dat) - 1) /3))
    tsne <- Rtsne::Rtsne(t(dat),perplexity=perplexity,check_duplicates=FALSE)

    # df <- data.frame(x = tsne$Y[,1],
    #                 y = tsne$Y[,2],
    #                 cellType = as.factor(meta[colnames(dat),"cellType"]),
    #                 sampleID = as.factor(meta[colnames(dat),"sampleID"]))

    df = cbind(tsne$Y,meta)
    colnames(df)=c("t1","t2",colnames(meta))

    ggplot(df, aes(t1, t2, colour = get(col))) +
    geom_point(size=ptsize) +
    ggtitle(plottitle)+
    xlab("T1")+
    ylab("T2")+
    theme_bw()+
    theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=22,face="bold"),
    axis.title=element_text(size=22,face="bold"))
}

run_cluster_plots <- function(imputedC,pdataC,savepath,trueC=NULL){
    dir.create(savepath,recursive=TRUE)
    if(!is.null(trueC)){
        mergeC = trueC
        mergeC[rownames(imputedC),colnames(imputedC)] = imputedC
        p_lr=logratio_plot(imputedC=mergeC,pDataC=pDataC,trueC=trueC,plottitle="MA Plot")
        p_umap=umap_plot(dat=mergeC,meta=pDataC,plottitle="UMAP Plot",ptsize=2)
        p_tsne=tsne_plot(dat=mergeC,meta=pDataC,plottitle="tSNE Plot",col="sampleID",ptsize=2)

        p_lr_control=logratio_plot(imputedC=trueC,pDataC=pDataC,trueC=trueC,plottitle="MA Plot (True)")
        p_umap_control=umap_plot(dat=trueC,meta=pDataC,plottitle="UMAP Plot (True)",ptsize=2)
        p_tsne_control=tsne_plot(dat=trueC,meta=pDataC,plottitle="tSNE Plot (True)",col="sampleID",ptsize=2)

        ggsave(plot=p_lr,filename=file.path(savepath,"ma.pdf"),width=5,height=5)
        ggsave(plot=p_umap,filename=file.path(savepath,"umap.pdf"),width=5,height=5)
        ggsave(plot=p_tsne,filename=file.path(savepath,"tsne.pdf"),width=5,height=5)

        ggsave(plot=p_lr_control,filename=file.path(savepath,"ma_control.pdf"),width=5,height=5)
        ggsave(plot=p_umap_control,filename=file.path(savepath,"umap_control.pdf"),width=5,height=5)
        ggsave(plot=p_tsne_control,filename=file.path(savepath,"tsne_control.pdf"),width=5,height=5)
    }else{
        p_umap=umap_plot(dat=imputedC,meta=pDataC,plottitle="UMAP Plot",ptsize=2)
        p_tsne=tsne_plot(dat=imputedC,meta=pDataC,plottitle="tSNE Plot",col="sampleID",ptsize=2)

        ggsave(plot=p_umap,filename=file.path(savepath,"umap.pdf"),width=5,height=5)
        ggsave(plot=p_tsne,filename=file.path(savepath,"tsne.pdf"),width=5,height=5)
    }
}

get_source_data <- function(sourcepath,prefix){
    fn_T = file.path(sourcepath,paste0(prefix[2],"_T.csv"))
    fn_pDataC = file.path(sourcepath,paste0(prefix[1],"_pDataC.csv"))
    fn_C = file.path(sourcepath,paste0(prefix[1],"_C.csv"))
    fn_trueC = file.path(sourcepath,paste0(prefix[1],"_trueC.csv"))
    fn_trueP = file.path(sourcepath,paste0(prefix[2],"_trueP.csv"))
    print(paste0("exists T=",file.exists(fn_T)))
    print(paste0("exists pdatac=",file.exists(fn_pDataC)))
    print(paste0("exists c=",file.exists(fn_C)))
    print(paste0("exists truec=",file.exists(fn_trueC)))
    print(paste0("exists truep=",file.exists(fn_trueP)))

    if(file.exists(fn_trueC)){
        C = read.csv(fn_C,row.names=1)
        T = read.csv(fn_T,row.names=1)
        pDataC = read.csv(fn_pDataC,row.names=1)
        trueC = read.csv(fn_trueC,row.names=1)
        trueP = read.csv(fn_trueP,row.names=1)
    }else if(file.exists(file.path(sourcepath,"trueC.csv"))){
        C = read.csv(file.path(sourcepath,"C.csv"),row.names=1)
        T = read.csv(file.path(sourcepath,"T.csv"),row.names=1)
        pDataC = read.csv(file.path(sourcepath,"pDataC.csv"),row.names=1)
        trueC = read.csv(file.path(sourcepath,"trueC.csv"),row.names=1)
        trueP = read.csv(file.path(sourcepath,"trueP.csv"),row.names=1)
    }else if(file.exists(fn_C)){
        C = read.csv(fn_C,row.names=1)
        T = read.csv(fn_T,row.names=1)
        pDataC = read.csv(fn_pDataC,row.names=1)
        trueC = NULL
        trueP = NULL
    }else if(file.exists(file.path(sourcepath,"C.csv"))){
        C = read.csv(file.path(sourcepath,"C.csv"),row.names=1)
        T = read.csv(file.path(sourcepath,"T.csv"),row.names=1)
        pDataC = read.csv(file.path(sourcepath,"pDataC.csv"),row.names=1)
        trueC = NULL
        trueP = NULL
    }
    return(list(C=C,T=T,pDataC=pDataC,trueC=trueC,trueP=trueP))
}