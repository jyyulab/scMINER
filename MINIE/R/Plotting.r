##Function based MINIE
##Author:chenxi.qian@stjude.org
##Stjude.YuLab

#' @export
MICAplot<-function(input_eset=eset,label= metaName, visualize= "UMAP",title.size=5,title.name="",pct=0.5){

  if(!label%in%colnames(pData(input_eset))){stop("Label name not contained in phenotype data!","\n")}
  if(visualize=="UMAP") (stop("Must specify 2D embedding methods!","\n"))

  p <- ggplot(data=pData(input_eset),aes_string(x=X,y=Y,color = label))+
       geom_point(size=pct)+
       labs(title=title.name)+
       theme(plot.title = element_text(size = title.size, face = "bold"),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 15))+
       guides(colour = guide_legend(override.aes = list(size = 10)))+
       xlab(label = paste0(toString(visualize),"_1"))+
       ylab(label = paste0(toString(visualize),"_2"))

  return(p)

}


#' @export
gene_highlighting<-function(input_eset=eset,target,
                            ylabel="Expression",title.size=5){
  # change it to expr is ok
  input<-as.matrix(exprs(input_eset))
  indx<-which(rownames(input)%in%target)
  projection<-pData(input_eset)[colnames(input),c("tSNE_1","tSNE_2")]
  #gene expression visualized as columns
  if (length(indx)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    }else
      {target_values <- input[indx,]}#indx = 1

   	proj_target <- cbind(projection,target_values)
    proj_target_melt <- reshape2::melt(proj_target, id.vars=c("tSNE_1", "tSNE_2"))

      p<- ggplot(proj_target_melt, aes(tSNE_1, tSNE_2)) +
        geom_point(aes(colour=value),size=0.5) +
        facet_wrap(~variable,scales = "free")+
        labs(title="Markers highlighting",scales = "free") +
        scale_colour_gradientn(colors=rev(brewer.pal(11,"RdYlBu")))   +
        theme(plot.title = element_text(size = title.size, face = "bold"),
              axis.title = element_text(size = 10),
              legend.title = element_text(size = 10))+
        labs(x="Tsne_1",y="Tsne_2",color=ylabel)
  return(p)
}

###Function3: Gene violinplot from eset###
#' @export
gene_vlnplot <- function(eset,group_tag="celltype",target,
                         ylabel="Expression",ncol=3,
                         drawquantiles=FALSE,title.size=5){

  # extract input information
  input <- exprs(eset)
  target <- intersect(target,rownames(input))
  label <- as.factor(pData(eset)[,group_tag]);names(label) <- rownames(pData(eset))
  # Gene expression visualized as columns
  if (length(target)!=1) {
    target_values <- t(as.matrix(input[target,]))
    label<-label[rownames(target_values)]
    df <- data.frame(target_values,cluster=label)
  }else {
    target_values<-input[target,]
    label<-label[names(target_values)]
    df <- data.frame(target_values,cluster=label)
    colnames(df)[1]<-target
  }

  df_melt <- reshape2::melt(df, id.vars="cluster")

  p <- ggplot(df_melt, aes(x=cluster, y=value,fill=cluster))

  if (drawquantiles) p <- p + geom_violin(trim=TRUE,scale="width",weight=0.1,draw_quantiles = c(0.25,0.5,0.75),na.rm = TRUE)
  else p <- p + geom_violin(trim=TRUE,scale="width",na.rm = TRUE,size=0.1,width=0.5)

  p <- p + facet_wrap(~variable,scales = "free",ncol = ncol) +
    labs(x="Cluster",y=ylabel)+
    theme(axis.text.x = element_text(size=10),
          plot.title = element_text(size = title.size, face = "bold"),
          strip.background = element_rect(fill="#FFFFFF"))

  if (ylabel=="Activity") {p <- p + geom_boxplot(width=0.2,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)}

  return(p)
}

###Function6: Heatmap visualization from eset###needs editting!!!!!!!!!
#' @export
gene_heatmap <- function(eset,target,
                         group_tag="label",name="log2_expression",
                         save_plot=TRUE,width=4,height=8,
                         plot_name="GeneHeatmap.png"){

  target<-intersect(target,featureNames(eset))
  exp<-exprs(eset)[target,]
  lab<-pData(eset)[,group_tag];names(lab) <- sampleNames(eset)

  #re-order expressionmatrix and label
  ranks<-names(sort(lab,decreasing = FALSE))
  exp.ordered<-as.matrix(exp[,ranks])
  lab.ordered<-lab[ranks]
  df<-data.frame(scMINER=lab.ordered)

  #Define color annotations
  n<-length(unique(lab.ordered))
  ncols <- hue_pal()(n)
  names(ncols) <- unique(lab.ordered)
  myanndf = HeatmapAnnotation(df = df,col=list(scMINER = ncols))
  mycolors = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

  hmp <- Heatmap(exp.ordered, col = mycolors, name = name,
          show_row_names = TRUE,
          show_column_names = FALSE,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          top_annotation = myanndf)

  if(save_plot){
    png(filename = plot_name, width=width,height=height,units ="in",res = 300)
    draw(hmp)
    dev.off()
  }

  return(hmp)
}

###Function7: AssignCelltypes###
#' @export
AssignCellTypes.Hmp<-function(ref = NULL,eset = eset.demo,
                              save_plot = FALSE,
                              width=8.5, height=6.5,
                              plot_name="AnnotationHeatmap.png"){
  #start from eset
  #z-normalize each sample
  exp<-apply(exprs(eset),2,std)
  #filter reference marker sets
  ref<-filter(ref,markers%in%rownames(exp))
  celltypes<-unique(ref$celltype)

  ac<-matrix(NA,nrow=ncol(exp),ncol=length(celltypes),dimnames = list(colnames(exp),celltypes))

  for(i in 1:length(celltypes)){
    cat(i,"\n")
    ref.sel<-filter(ref,celltype==celltypes[i])
    n <- length(ref.sel$markers)

    if(n>1){
      mat<-t(exp[ref.sel$markers,])%*%as.vector(ref.sel$weight)
      ac[,i]<-mat[,1]/n
    }else if (n==1){
      ac[,i]<-exp[ref.sel$markers,]
    }
  }

  df<-data.frame(label=eset$label,ac);
  df<-df[,colSums(is.na(df))<nrow(df)];#remove NA columns
  df<-aggregate(.~label,df,mean)

  input<-t(df[,-1])

  colnames(input)<-1:length(unique(eset$label))
  myanndf<-data.frame(row.names=1:length(unique(eset$label)),scMINER=as.factor(1:length(unique(eset$label))))

  hmp<-pheatmap::pheatmap(input,kmeans_k = NA,
                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(30),
                          cluster_rows = FALSE,cluster_cols = FALSE,
                          scale = "row",
                          show_rownames = TRUE,
                          show_colnames = FALSE,
                          annotation_col= myanndf)
  hmp
  if(save_plot){ggsave(hmp,filename = plot_name ,device="png",width = width,height = height,dpi = 300)}
  return(hmp)
}
