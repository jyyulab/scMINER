#'@import Biobase ggplot2 kableExtra knitr limma
#'@importFrom reshape2 melt
#'@importFrom rhdf5 h5createFile h5write H5close
#'@importFrom ComplexHeatmap HeatmapAnnotation Heatmap draw
#'@importFrom igraph graph_from_data_frame set_edge_attr E
#'@importFrom rmarkdown render pandoc_available html_document
#'@importFrom utils read.delim write.table
#'@importFrom Matrix rowSums colSums readMM t
#'@importFrom dplyr filter select starts_with arrange desc left_join
#'@importFrom plyr ddply
#'@importFrom methods as new
#'@importFrom stats IQR aggregate as.dendrogram as.dist cutree density dist fisher.test gaussian glm mad t.test hclust kmeans ks.test lm median model.matrix na.omit order.dendrogram p.adjust pchisq pnorm prcomp pt qnorm quantile sd splinefun complete.cases
#'@importFrom utils read.delim write.table read.table
#'@importFrom grDevices colorRampPalette dev.off png
#'@importFrom RColorBrewer brewer.pal
#'@importFrom scales hue_pal
#'
#'
#############
# library(scales)
# library(ggplot2)
# library(reshape2)
# library(ComplexHeatmap)
# library(RColorBrewer)## for colors coding
# library(grDevices)
# library(anndata)
#
# library(Matrix)
# library(stats) ##  t.test
# library(methods)
# library(dplyr) #for easy filtering and apply function
# library(plyr) #for ddply
#
# library(Biobase) ## basic functions for bioconductor
#
# library(rmarkdown)
# library(kableExtra) # for Rmarkdown
# library(knitr)      # for Rmarkdown
#############
check_param <- function(para_name,envir){
  if(base::exists(para_name,envir=envir)==FALSE){message(sprintf('%s missing !',para_name));return(0)}
  if(is.null(base::get(para_name,envir=envir))==TRUE){message(sprintf('%s is NULL !',para_name));return(0)}
  return(1)
}

check_options <- function(para_name,option_list,envir){
  if(!base::get(para_name,envir=envir) %in% option_list){
    message(sprintf('Only accept %s set at: %s !',para_name,base::paste(option_list,collapse=';')));return(0)
  }
  return(1)
}

#############
#' @title SparseExpressionSet
#' @exportClass SparseExpressionSet
#' @importFrom Biobase ExpressionSet
setClass( "SparseExpressionSet",
          contains = "ExpressionSet",
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c(classVersion("ExpressionSet"), SparseExpressionSet = "1.0.0" ))))

###Activity inner function###
std<-function(x){
  x<-x[!is.na(x)]
  (x-mean(x))/sd(x)
}

es <- function(z,es.method="mean"){
  if(es.method=="maxmean"){
    n<-length(z)
    m1<-ifelse(sum(z>0)>0,sum(z[z>0])/n,0)
    m2<-ifelse(sum(z<0)>0,sum(z[z<0])/n,0)
    if(m1>-m2) es<-m1
    else es<-m2
  }
  else if(es.method=='absmean'){
    es<-mean(abs(z))
  }
  else if(es.method == 'mean'){
    es<-mean(z)
  }
  else{
    stop('Unsupported method!\n')
  }
  return(es)
}

#' DAG_ttest
DAG_ttest<-function(d,group){

  #d.tmp<-unlist(d[1,-1])
  d.tmp<-unlist(d[,-1])
  d.tmp<-data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)

  d.sel<-d.tmp[complete.cases(d.tmp),]
  d.sel$group<-as.factor(d.sel$group)

  n.cases<-nlevels(d.sel$group)

  cat("=")
  if(n.cases==2){
    res <- t.test(dplyr::filter(d.tmp,group=="Aim")$acs,dplyr::filter(d.tmp,group=="Ctrl")$acs)
    pval.t <-res$p.value
    t.stat <- res$statistic
    #extract results
    rs.t <- c(id=d$acs.id,
              t.stat,
              pval = pval.t,
              res$parameter,#df
              CI.low = res$conf.int[1],
              CI.high = res$conf.int[2],
              MeanAct.Aim = unname(res$estimate[1]),
              MeanAct.Ctrl = unname(res$estimate[2]))

  }else{
    rs.t <- c(id = d$acs.id,
              t=NA,
              pval = NA,
              df = NA,
              CI.low = NA,
              CI.high = NA,
              MeanAct.Aim = mean(filter(d.tmp,group=="Aim")$acs),
              MeanAct.Ctrl = mean(filter(d.tmp,group=="Ctrl")$acs))
  }

  rs.t<-c(rs.t, log2FC=rs.t["MeanAct.Aim"]-rs.t["MeanAct.Ctrl"])

  return(rs.t)
}

#' SJARACNe_filter
#' @description This is the inner function to help generate SJARACNe input for scRNA-seq data,
#' all non-informative (zero genes) will be filtered in by this function
#' @param eset.sel ExpressionSet to generate SJaracne input
#' @param tf.ref A vector of reference transcription factors
#' @param sig.ref A vector of reference signaling genes
#' @param wd.src path to store SjAracne input
#' @param grp.tag name of group for identification
#' @details Non-expressed genes in subgroups are filtered.
#' tf.ref should be coordinate with featureNames(eset.sel).
#' @return A folder with picked master regulator and filtered gene expression matrix
SJARACNe_filter<-function(eset.sel,tf.ref,sig.ref,wd.src,grp.tag){
  cat(grp.tag,'\n')

  #fData(eset.sel)$IQR<-apply(exprs(eset.sel),1,IQR)
  # exclude genes with all zero
  eset.sel<-eset.sel[apply(exprs(eset.sel),1,function(xx){sum(xx)!=0}),]

  fData(eset.sel)$geneNames<-fData(eset.sel)$geneSymbol

  ni<-nrow(eset.sel);ni
  ns<-ncol(eset.sel);ns
  ng<-nlevels(factor(fData(eset.sel)$geneNames));ng
  tag<-paste(grp.tag,nrow(eset.sel),ng,ncol(eset.sel),sep='_');tag

  dir.cur<-file.path(wd.src,tag);dir.cur
  dir.create(dir.cur,recursive = T)

  #write exp data to exp format
  expdata<-data.frame(cbind(isoformId=featureNames(eset.sel),geneSymbol=fData(eset.sel)$geneSymbol,as.matrix(exprs(eset.sel))),stringsAsFactors = FALSE)
  f.exp<-file.path(dir.cur,paste(grp.tag,"_",ni,"_",ng,"_",ns,".exp",sep=''));f.exp
  write.table(expdata,file=f.exp,sep="\t",row.names=FALSE,quote=FALSE)

  if (!is.null(tf.ref)){
    dir.create(file.path(dir.cur,'tf'),recursive = T)
    tf.eset.sel<-subset(eset.sel,fData(eset.sel)$geneSymbol%in%tf.ref)
    dim(tf.eset.sel)
    f.tf<-file.path(dir.cur,'tf',paste(grp.tag,"_",nrow(tf.eset.sel),"_",nlevels(factor(fData(tf.eset.sel)$geneSymbol)),"_",ns,"_tf.txt",sep=''));f.tf
    cat(featureNames(tf.eset.sel),file=f.tf,sep='\n')
  }

  if (!is.null(sig.ref)){
    dir.create(file.path(dir.cur,'sig'),recursive = T)
    sig.eset.sel<-subset(eset.sel,fData(eset.sel)$geneSymbol%in%sig.ref)
    dim(sig.eset.sel)
    f.sig<-file.path(dir.cur,'sig',paste(grp.tag,"_",nrow(sig.eset.sel),"_",nlevels(factor(fData(sig.eset.sel)$geneSymbol)),"_",ns,"_sig.txt",sep=''));f.sig
    cat(featureNames(sig.eset.sel),file=f.sig,sep='\n')
  }
  return(list(dir.cur,grp.tag))
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
#' @return a ggplot object
draw.bubblePlot2<-function(df=NULL,xlab,ylab,clab,slab,
                           low.col="#004C99",high.col="#CC0000",plot.title=NULL){

  p <- ggplot(df, aes_string(x= xlab, y= ylab)) +
    theme_classic()+
    geom_point(aes_string(fill=clab, size= slab),color="black",pch=21)+
    scale_fill_gradient2(low=low.col,high=high.col)+
    scale_x_discrete(limits=levels(df[,xlab]))+
    scale_y_discrete(limits=levels(df[,ylab]))+
    theme(panel.grid.major= element_line(colour = "grey",size=0.3),
          panel.grid.minor = element_line(colour = "grey",size=0.3),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  labs(x = xlab, y = ylab, title = plot.title)

  return(p)
}

