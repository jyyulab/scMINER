##Function based MINIE
##Author:chenxi.qian@stjude.org
##Stjude.YuLab

#' @export
MICAplot<-function(input_eset=eset,label= metaName,visualize=NULL,
                   X=NULL,Y=NULL,
                   title.size=5,title.name="",pct=0.5){

  if(!label%in%colnames(pData(input_eset))){stop("Label name not contained in phenotype data!","\n")}

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
                         cluster_rows=FALSE,
                         plot_name="GeneHeatmap.png",
                         ...){

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

#inner function to draw pretty bubbleplot
draw.bubblePlot<-function(bb_size=NULL,bb_color=NULL,pdf_file=NULL,
                          gs_cex=1,driver_cex=1,
                          plot_target_size=FALSE,Z_val=NULL,main=""
                          )
{

  f_mat1<-bb_size[,-1]
  f_mat2<-bb_color[,-1]
  nr <- ncol(f_mat1)
  nc <- nrow(f_mat1)

  ## modified from NetBID2, credit to Sherry Dong
  gsWidth  <- max(strwidth(colnames(f_mat1),'inches',cex=gs_cex))
  gsHeight <- max(strheight(colnames(f_mat1),'inches',cex=gs_cex)*nrow(f_mat1))
  driverWidth  <- max(strwidth(show_label[colnames(f_mat1)],'inches',cex=gs_cex))
  driverHeight <- max(strheight(show_label[colnames(f_mat1)],'inches',cex=gs_cex)*ncol(f_mat1))

  ww <- (1+nc)*0.5+ gsWidth + 1.5 + 2
  hh <- (4+nr)*0.5+ driverWidth + 2+1

  ## output to pdf
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=ww,height=hh)

  layout(1);par(mai=c(driverWidth+2,gsWidth+1.5,1,2))

  plot(1,bty='n',col='white',
       xlim=c(0,nc+1),ylim=c(-2,nr+1),
       xaxt='n',yaxt='n',xlab="",ylab="",main=main)

  segments(x0=0,x1=nc,y0=0:nr,y1=0:nr,col='dark grey',xpd=TRUE)
  segments(x0=0:nc,x1=0:nc,y0=0,y1=nr,col='dark grey',xpd=TRUE)
  segments(x0=0:nc,x1=0:nc,y0=0,y1=-3,col='grey',xpd=TRUE)

  #
  text(-0.5,1:nr-0.5,colnames(f_mat1),xpd=TRUE,srt=0,adj=1,cex=gs_cex) ## sig pathways
  if(is.null(mark_gene)==TRUE){
    text(1:nc-0.5,-2.5,show_label[rownames(f_mat1)],xpd=TRUE,srt=90,adj=1,cex=driver_cex) ## sig regulators
  }else{
    bc <- rep('black',length.out=length(rownames(f_mat1)))
    bc[which(rownames(f_mat1) %in% mark_gene)] <- 'red'
    print(table(bc))
    text(1:nc-0.5,-2.5,show_label[rownames(f_mat1)],xpd=TRUE,srt=90,adj=1,col=bc,cex=driver_cex) ## sig regulators
  }

  ## draw circle
  max_size <- max(f_mat1,na.rm=TRUE)
  f_mat1 <- f_mat1/max_size

  cc_r <- matrix(num2col(f_mat2,n_len=30,threshold=0.01),ncol=ncol(f_mat2),byrow = FALSE)

  for(i in 1:nrow(f_mat1)){
    for(j in 1:ncol(f_mat1)){
      draw.circle(i-0.5,j-0.5,radius=f_mat1[i,j]/20,col=cc_r[i,j])
    }#library(plotrix)
  }
  ## draw circle legend
  legend_size <- unique(round(seq(1,max_size,length.out=5)))
  for(i in 1:length(legend_size)){
    draw.circle(length(legend_size)-i+1.5,nr+0.5,radius=0.5*legend_size[i]/max_size)
    text(length(legend_size)-i+1.5,nr+1,legend_size[i])
  }
  text(0.5,nr+0.5,'Size ')

  ## draw p-value legend
  {
  p_label <- c(1,0.1,0.05,0.01,0.001,1e-4,1e-10)
  p_thre <- qnorm(1-c(1,0.1,0.05,0.01,0.001,0.0001,1e-10))

  p_col  <- num2col(p_thre,n_len=30,sig_thre=qnorm(1-0.1))
  p_col_m  <- num2col(-p_thre,n_len=30,sig_thre=qnorm(1-0.1))
  ybottom <- seq(nr-6,nr-1,length.out=1+length(p_thre))[1:length(p_thre)]
  ytop <- seq(nr-6,nr-1,length.out=1+length(p_thre))[2:(length(p_thre)+1)]

  for(i in 1:length(p_thre)){
    rect(nc+0.5,ybottom[i],nc+1.5,ytop[i],col=p_col[i],xpd=TRUE)
    text(nc+1.7,(ybottom[i]+ytop[i])/2,adj=0,p_label[i],xpd=TRUE)
  }

  ybottom <- seq(nr-11,nr-6,length.out=1+length(p_thre))[1:length(p_thre)]
  ytop <- seq(nr-11,nr-6,length.out=1+length(p_thre))[2:(length(p_thre)+1)]
  for(i in 2:length(p_thre)){
    rect(nc+0.5,ybottom[i],nc+1.5,ytop[i],col=rev(p_col_m)[i-1],xpd=TRUE)
    text(nc+1.7,(ybottom[i]+ytop[i])/2,adj=0,rev(p_label)[i-1],xpd=TRUE)
  }

  text(nc+1.5,nr-0.5,'P-Value',xpd=TRUE)
  }
  # target size
  if(plot_target_size){
    max_target_size <- max(unlist(lapply(target_gene,function(x)max(unlist(lapply(x,length))))))
    ori_size <- unlist(lapply(target_gene,function(x)length(x[[1]])))
    pro_size <- unlist(lapply(target_gene,function(x)length(x[[2]])))

    rect(0:(nc-1)+0.15,-2,0:(nc-1)+0.45,1.5*ori_size/max_target_size-2,col='blue')
    rect(0:(nc-1)+0.55,-2,0:(nc-1)+0.85,1.5*pro_size/max_target_size-2,col='green')
    segments(y0=-2,y1=-0.5,x0=0,x1=0)
    segments(y0=seq(-2,-0.5,length.out=3),y1=seq(-2,-0.5,length.out=3),x0=-0.25,x1=0)

    text(-0.3,seq(-2,-0.5,length.out=3),round(seq(0,max_target_size,length.out=3)),adj=1,xpd=TRUE)
    legend(-5,-0.5,c('target_size','target_size\n(protein_coding)'),fill=c('blue','green'),border=NA,bty='n',cex=0.8,xpd=TRUE)}

  # add sig color
  if(!is.null(Z_val)){
  sig_col <- num2col(Z_val[driver_list],n_len=30)
  rect(0:(nc-1)+0.35,-0.4,0:(nc-1)+0.65,-0.1,col=sig_col)
  }
  # add driver type !!!
  if(is.null(driver_type)==FALSE){
    cc_tmp <- get.class.color(unique(driver_type[driver_list]))
    points(x=0:(nc-1)+0.5,y=rep(-2.2,length.out=nc),col=cc_tmp[driver_type[driver_list]],pch=16,cex=2.5)
    legend(nc+0.5,0.5,names(cc_tmp),fill=cc_tmp,border=NA,bty='n',cex=0.8,xpd=TRUE)
  }
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)

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


