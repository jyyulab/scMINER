##Function based MINIE
##Author:chenxi.qian@stjude.org
##Stjude.YuLab

#' @title plot MICA
#' @description This function help to generate a ggplot object for phenotypic visualization
#' @param input_eset ExpressionSet that include visualization coordinates in phenotype data
#' @param label Coloring criteria of data points on
#' @param visualize character, name of visualization method, this will be used on x or y axis title
#' @param X character, column name of x axis
#' @param Y character, column name of y axis
#' @param title.size numerical, size of plot title, default as 10
#' @param title.name character, title of plot, default as NULL
#' @param pct numerical, size of point, default as 0.5
#'
#' @export
MICAplot<-function(input_eset=eset,label= metaName,visualize=NULL,
                   X=NULL,Y=NULL,
                   title.size=10,title.name="",pct=0.5){

  if(!label%in%colnames(pData(input_eset))){stop("Label name not contained in phenotype data!","\n")}

  p <- ggplot(data=pData(input_eset),aes_string(x=X,y=Y,color = label))+
       geom_point(size=pct)+
       labs(title=title.name)+
       theme_classic()+
       theme(plot.title = element_text(size = title.size, face = "bold"),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 15))+
       guides(colour = guide_legend(override.aes = list(size = 10)))+
       xlab(label = paste0(toString(visualize),"_1"))+
       ylab(label = paste0(toString(visualize),"_2"))

  return(p)

}


#' @title Visualize gene expression level on scRNA-seq data
#' @description This plot will visualiz feature info in scatter plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param ylabel a characterm, title of y axis
#' @param title.size numerical, default as 5
#' @param x cordinates for x axis
#' @param y cordinates for y axis
#'
#' @export
feature_highlighting<-function(input_eset=eset,target,
                               feature="geneSymbol",
                               x="tSNE_1",y="tSNE_2",
                               ylabel="Expression",
                               title.size=5){

  # change it to expr is ok
  input<-as.matrix(exprs(input_eset))
  indx<-which(fData(input_eset)[,feature]%in%target)
  gn<-fData(input_eset)[,feature][indx]
  projection<-pData(input_eset)[colnames(input),c(x,y)]

  #gene expression visualized as columns
  if (length(indx)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    colnames(target_values)<-gn
    }else
      {target_values <- input[indx,]}#indx = 1

   	proj_target <- cbind(projection,target_values)
    proj_target_melt <- reshape2::melt(proj_target, id.vars=c(x, y))

    p<- ggplot(proj_target_melt, aes_string(x, y)) +
        theme_classic()+
        geom_point(aes(colour=value),size=0.5) +
        facet_wrap(~variable,scales = "free")+
        labs(title="",scales = "free") +
        scale_colour_gradientn(colors=rev(brewer.pal(11,"RdYlBu")))   +
        theme(plot.title = element_text(size = title.size, face = "bold"),
              axis.title = element_text(size = 10),
              legend.title = element_text(size = 10))+
        labs(x="Tsne_1",y="Tsne_2",color=ylabel)
  return(p)
}


#' @title Visualize gene expression level on scRNA-seq data
#' @description This plot will visualize feature info in violin plot by outputing a ggplot object
#' @param eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param ylabel a character, title of y axis
#' @param boxplot logical, whether to plot boxplot on violinplot
#' @param group_tag character, which group info
#' @param title.size numerical, default as 5
#'
#' @export
feature_vlnplot <- function(eset,group_tag="celltype",
                         target,feature="geneSymbol",
                         ylabel="Expression",ncol=3,
                         boxplot=FALSE,title.size=5){

  # extract input information
  input <- exprs(eset)
  indx<-which(fData(eset)[,feature]%in%target)
  gn<-fData(eset)[,feature][indx]

  label <- as.factor(pData(eset)[,group_tag])
  # Gene expression visualized as columns
  if (length(target)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    colnames(target_values)<-gn
    df <- data.frame(target_values,cluster=label)
  }else {
    target_values<-input[indx,]
    df <- data.frame(target_values,cluster=label)
    colnames(df)[1]<-gn
  }

  df_melt <- reshape2::melt(df, id.vars="cluster")

  p <- ggplot(df_melt, aes(x=cluster, y=value,fill=cluster))+
       theme_classic()+
       geom_violin(trim=TRUE,scale="width",na.rm = TRUE,size=0.1,width=0.5)+
       stat_summary(fun.y=median, geom="point", size=1.2, color="black")

  if(boxplot) p <- p + geom_boxplot(fill="white",width=0.1,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)

  p <- p + facet_wrap(~variable,scales = "free",ncol = ncol) +
    labs(x="Cluster",y=ylabel)+
    theme(axis.text.x = element_text(size=10),
          plot.title = element_text(size = title.size, face = "bold"),
          strip.background = element_rect(fill="#FFFFFF"))

  if (ylabel=="Activity") { p <- p + geom_boxplot(width=0.2,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)}

  return(p)
}


#' @title Visualize gene expression level on scRNA-seq data
#' @description This plot will visualiz feature info in scatter plot by outputing a ggplot object
#' @param eset Input expression set
#' @param group_tag a character, label to visualize on the top of heatmap
#' @param feature character, which feature to visualize
#' @param name character, name of value visualized in color scale
#' @param cluster_rows logical, if or not cluster rows
#' @param plot_name character, name of heamap
#' @param save_plot logical, whether to save plots or not
#' @param width numerical
#' @param height numerical
#'
#' @export
feature_heatmap <- function(eset,target,feature="geneSymbol",
                         group_tag="label",name="log2Exp",
                         save_plot=TRUE,width=4,height=8,
                         cluster_rows=FALSE,
                         plot_name="GeneHeatmap.png",
                         ...){

  input <- exprs(eset)
  indx<-which(fData(eset)[,feature]%in%target)
  gn<-fData(eset)[,feature][indx]

  exp<-exprs(eset)[indx,]
  rownames(exp)<-gn
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
          cluster_rows = cluster_rows,
          cluster_columns = FALSE,
          top_annotation = myanndf,
          ...)

  if(save_plot){
    png(filename = plot_name, width=width,height=height,units ="in",res = 300)
    draw(hmp)
    dev.off()
  }

  return(hmp)
}


#' Inner function for simple bubbleplots
#' @param df re-structured data.frame for bubble plots
#' @param xlab string
#' @param ylab string
#' @param clab string
#' @param slab string
#' @param low.col string,default as "#004C99"
#' @param high.col string, default as "CC0000
#' @param plot.title string
#'
#' @export
draw.bubblePlot2<-function(df=NULL,xlab,ylab,clab,slab,
                           low.col="#004C99",high.col="#CC0000",plot.title=NULL){

  p <- ggplot(df, aes_string(x= xlab, y= ylab)) +

    geom_point(aes_string(color=clab, size= slab))+

    scale_color_gradient2(low=low.col,high=high.col)+

    theme_minimal()+ # minimal theme

    scale_x_discrete(limits=levels(df[,xlab]))+

    scale_y_discrete(limits=levels(df[,ylab]))+

    theme(axis.text.x = element_text(size = 8),

          axis.text.y = element_text(size = 8))

    labs(x = xlab, y = ylab, title = plot.title)

  return(p)
}





#' draw.scRNAseq.QC
#' @title draw.scRNAseq.QC
#' @description generated a scRNA-seq quality control report in html with Rmarkdown
#' @param SparseEset an SparseEset generated by CreateSparseEset
#' @param project.name a character, project name to print on report
#' @param plot.dir a character, output directory for QC reports
#'
#' @return an R markdown QC report
#' @export
draw.scRNAseq.QC<-function(SparseEset,project.name,plot.dir="./QC/",output.cutoff=TRUE,group="group"){
  if(!dir.exists(plot.dir)) {dir.create(plot.dir)}

  #Calcualte Cutoffs
  pd<-pData(SparseEset)
  cfs<-list(nCell_cutoff = max(floor(0.005 * dim(SparseEset)[2]), 1),
            umi_cf_lo = max(floor(exp(median(log(pd$nUMI.total)) - 3 * mad(log(pd$nUMI.total)))),100),
            umi_cf_hi = ceiling(exp(median(log(pd$nUMI.total)) + 3 * mad(log(pd$nUMI.total)))),
            nGene_cf = max(floor(exp(median(log(pd$nGene)) - 3 * mad(log(pd$nGene)))),50),
            ERCC_cf = round(median(pd$percent.spikeIn) + 3 * mad(pd$percent.spikeIn),4),
            mito_cf = round(median(pd$percent.mito) + 3 * mad(pd$percent.mito),4))

  render(input=system.file("rmd", "SparseEset_QC_report.Rmd", package = "scMINER"),
           output_dir = plot.dir,
           output_file = paste0(project.name,"_scRNAseq_QC.html"),
           clean=TRUE,
           quiet =FALSE,
           params=list(
             Obj=SparseEset,
             projectName=project.name,
             cutoffs=cfs,
             plot.dir=plot.dir,
             output.cutoff=output.cutoff,
             group=group))

  return(cfs)
}

