##Inner function of network visualization
##Author: xinran.dong@stjude.org
##Stjude.YuLab

#' getGSC
#'
#' @param sig signaling genes network table
#' @param tf transcription factors network table
#'
#' @return a list of master regulators and their targets
#' @export
getGSC<-function(sig=NULL,tf=NULL){
  sig.gsc<-NULL
  tf.gsc<-NULL

  if(!is.null(sig)){
    src <-unique(sig$source)
    n<-length(src)

    sig.gsc <-vector("list",n)
    for(i in 1:n){
      tag<-src[i]
      tmp<-filter(sig,sig$source==tag)
      names(sig.gsc)[i] <-paste(tag,"SIG",sep = ".")
      sig.gsc[[i]] <-as.character(tmp$target)}
  }

  if (!is.null(tf)){
    src <-unique(tf$source)
    n<-length(src)
    tf.gsc <-vector("list",n)
    for(i in 1:n){
      tag<-src[i]
      tmp<-filter(tf,tf$source==tag)
      names(tf.gsc)[i] <-paste(tag,"TF",sep = ".")
      tf.gsc[[i]] <-as.character(tmp$target)}}

  gsc<-c(sig.gsc,tf.gsc)
  return(gsc)
}




###Function: Inner function of Network rewiring visualization####
##credit to: Sherry Dong
#' Target network structure plot for two drivers.
#'
#' \code{draw.targetNet.TWO} will draw the network structure for the selected two drivers and their target genes.
#'
#' This is a function to draw target network structure for the selected two drivers.
#' The color bar represents the positive (red) or negative (blue) regulation with line width showing the strength.
#'
#' @param source1_label the label for the first(left) driver displayed on the plot.
#' @param source2_label the label for the second(right) driver displayed on the plot.
#' @param source1_z numeric, the Z statistic for the first driver, used to color the driver point.
#' If NULL, the driver will be colored in grey. Default is NULL.
#' @param source2_z numeric, the Z statistic for the second driver, used to color the driver point.
#' If NULL, the driver will be colored in grey. Default is NULL.
#' @param edge_score1 a vector of numeric values, indicating the correlation between the first driver and the target genes.
#' The value ranges from -1 to 1, with positive value indicating postivie regulation and negative value indicating negative correlation.
#' The names for the vector is the gene labels displayed on the plot.
#' @param edge_score2 a vector of numeric values, indicating the correlation between the second driver and the target genes.
#' Similar with \code{edge_score1}
#' @param label_cex numeric, \code{cex} for the target genes displayed on the plot. Default is 0.7.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param total_possible_target numeric or a vector of characters. If input numeric, will be the total number of possible targets.
#' If input a vector of characters, will be the background list of all possible target genes.
#' This parameter will be passed to function \code{test.targetNet.overlap} to test whether the target genes of the two drivers are significantly intersected.
#' If NULL, will do not perform this test. Default is NULL.
#' @param show_test logical, indicating whether the testing results will be printed and returned. Default is FALSE.
#'
#' @return if \code{show_test}==FALSE, will return logical value indicating whether the plot has been successfully generated, otherwise will return the statistics of testing.
#'
#' @examples
#' source1_label <- 'test1'
#' source1_z <- 1.96
#' edge_score1 <- (sample(1:160,size=80,replace=TRUE)-80)/80
#' names(edge_score1) <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' source2_z <- -2.36
#' edge_score2 <- (sample(1:240,size=120,replace=TRUE)-120)/120
#' names(edge_score2) <- sample(paste0('G',1:1000),size=120)
#' draw.targetNet.TWO(source1_label=source1_label,
#'                source2_label=source2_label,
#'                source1_z=source1_z,source2_z=source2_z,
#'                edge_score1=edge_score1,edge_score2=edge_score2,
#'                total_possible_target=paste0('G',1:1000),
#'                show_test=TRUE,label_cex=0.6)
#' \dontrun{
#' source1_label <- 'test1'
#' source1_z <- 1.96
#' edge_score1 <- (sample(1:160,size=100,replace=TRUE)-80)/80
#' names(edge_score1) <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' source2_z <- -2.36
#' edge_score2 <- (sample(1:240,size=100,replace=TRUE)-120)/120
#' names(edge_score2) <- sample(paste0('G',1:1000),size=120)
#' analysis.par <- list()
#' analysis.par$out.dir.PLOT <- getwd()
#' draw.targetNet.TWO(source1_label=source1_label,
#'                source2_label=source2_label,
#'                source1_z=source1_z,source2_z=source2_z,
#'                edge_score1=edge_score1,edge_score2=edge_score2,
#'                total_possible_target=paste0('G',1:1000),show_test=TRUE,
#'                pdf_file=sprintf('%s/targetNetTWO.pdf',
#'                analysis.par$out.dir.PLOT))
#' }
#' @export
draw.targetNet.TWO <- function(source1_label="",source2_label="",
                               source1_z=NULL,source2_z=NULL,
                               edge_score1=NULL,edge_score2=NULL,
                               label_cex=0.7,pdf_file=NULL,
                               total_possible_target=NULL,show_test=FALSE){
  tmp1 <- sapply(unique(names(edge_score1)),function(x){
    x1 <- edge_score1[which(names(edge_score1)==x)];x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score1));edge_score1 <- tmp1
  tmp1 <- sapply(unique(names(edge_score2)),function(x){
    x1 <- edge_score2[which(names(edge_score2)==x)];x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score2));edge_score2 <- tmp1
  edge_score1<- sort(edge_score1,decreasing = FALSE)
  edge_score2<- sort(edge_score2,decreasing = TRUE)
  g12 <- intersect(names(edge_score1),names(edge_score2))
  g1  <- setdiff(names(edge_score1),names(edge_score2))
  g2  <- setdiff(names(edge_score2),names(edge_score1))
  ec1 <- z2col(edge_score1*100,sig_thre=0,n_len=length(edge_score1));names(ec1) <- names(edge_score1)
  ec2 <- z2col(edge_score2*100,sig_thre=0,n_len=length(edge_score2));names(ec2) <- names(edge_score2)
  ew1 <- 2*label_cex*(abs(edge_score1)-min(abs(edge_score1)))/(max(abs(edge_score1))-min(abs(edge_score1)))+label_cex/2; names(ew1) <- names(edge_score1)
  ew2 <- 2*label_cex*(abs(edge_score2)-min(abs(edge_score2)))/(max(abs(edge_score2))-min(abs(edge_score2)))+label_cex/2; names(ew2) <- names(edge_score2)
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  geneWidth <- max(strwidth(g1,'inches',cex=label_cex))
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=10+4*geneWidth,height=8+2*geneWidth)
  plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab="",ylab="",xaxt='n',yaxt='n')
  par(mai=c(1,1,1,1))
  if(length(g1)>0){
    tt <- seq(-0.225,0.225,length.out=length(g1));init.angle <- -180;p1<-t2xy(tt,radius=0.8);
    for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g1[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    #p1<-t2xy(tt,radius=0.78);
    #p2<-t2xy(tt,radius=0.76);
    p1<-t2xy(tt,radius=0.8-label_cex/15);
    p2<-t2xy(tt,radius=0.8-label_cex/20);
    arrows(x0=-0.2,y0=0,x1=p1$x,y1=p1$y,col=ec1[g1],lwd=ew1[g1],angle=10,length=0.1*label_cex);
    points(p1$x,p1$y,pch=16,col='dark grey',cex=label_cex)
  }
  if(length(g2)>0){
    tt <- seq(-0.225,0.225,length.out=length(g2));init.angle <- 0;
    p1<-t2xy(tt,radius=0.8);
    for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g2[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    #p1<-t2xy(tt,radius=0.78);
    #p2<-t2xy(tt,radius=0.76);
    p1<-t2xy(tt,radius=0.8-label_cex/15);
    p2<-t2xy(tt,radius=0.8-label_cex/20);
    arrows(x0=0.2,y0=0,x1=p1$x,y1=p1$y,col=ec2[g2],lwd=ew2[g2],angle=10,length=0.1*label_cex)
    points(p1$x,p1$y,pch=16,col='dark grey',cex=label_cex)
  }
  if(length(g12)>0){
    tt <- seq(min(0.1*length(g12),0.7),-min(0.1*length(g12),0.7),length.out=length(g12));
    #segments(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12])
    #segments(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12])
    arrows(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12],angle=10,length=0.1*label_cex)
    arrows(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12],angle=10,length=0.1*label_cex)
    boxtext(0,tt,labels=g12,col.bg=get_transparent('light grey',0.3),cex=label_cex)
    #text(0,tt,g12,adj=0.5,cex=label_cex)
  }
  geneWidth1 <- strwidth(source1_label,'inches',cex=0.8)
  geneWidth2 <- strwidth(source2_label,'inches',cex=0.8)

  if(is.null(source2_z)==TRUE)
    points(0.2,0,col='light grey',cex=geneWidth1*18,pch=16)
  else
    points(0.2,0,col=z2col(source2_z),cex=geneWidth1*18,pch=16)

  text(0.2,0,source2_label,adj=0.5,cex=0.8)
  if(is.null(source1_z)==TRUE)
    points(-0.2,0,col='light grey',cex=geneWidth1*18,pch=16)
  else
    points(-0.2,0,col=z2col(source1_z),cex=geneWidth1*18,pch=16)

  text(-0.2,0,source1_label,adj=0.5,cex=0.8)
  # fisher test for target
  if(is.null(total_possible_target)==FALSE & show_test==TRUE){
    res <- test.targetNet.overlap(source1_label,source2_label,names(edge_score1),names(edge_score2),total_possible_target)
    if(is.null(pdf_file)==FALSE) dev.off()
    return(res)
  }
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}




#' Test for the target genes' intersection between two drivers.
#'
#' \code{test.targetNet.overlap} will test whether the target genes of two drivers are significantly intersected.
#'
#' This is a function to perform Fisher's Exact Test for the intersection of the target genes from two drivers.
#'
#' @param source1_label character, the label for the first driver.
#' @param source2_label character, the label for the second driver.
#' @param target1 a vector of characters, the list of target genes for the first driver.
#' @param target2 a vector of characters, the list of target genes for the second driver.
#' @param total_possible_target numeric or a vector of characters. If input numeric, will be the total number of possible targets.
#' If input a vector of characters, will be the background list of all possible target genes.
#'
#' @return Return the statistics of testing, including the \code{P.Value}, \code{Odds_Ratio} and \code{Intersected_Nuumber}.
#'
#' @examples
#' source1_label <- 'test1'
#' target1 <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' target2 <- sample(paste0('G',1:1000),size=120)
#' test.targetNet.overlap(source1_label=source1_label,source2_label=source2_label,
#'                target1=target1,target2=target2,
#'                total_possible_target=paste0('G',1:1000))
#' \dontrun{
#' }
#' @export
test.targetNet.overlap <- function(source1_label=NULL,source2_label=NULL,
                                   target1=NULL,target2=NULL,
                                   total_possible_target=NULL){
  t1  <- unique(target1)
  t2  <- unique(target2)
  print(sprintf('%s has %d unique targets !',source1_label,length(t1)))
  print(sprintf('%s has %d unique targets !',source2_label,length(t2)))
  n11 <- length(intersect(t1,t2))
  n12 <- length(setdiff(t1,t2))
  n21 <- length(setdiff(t2,t1))
  if(class(total_possible_target) %in% c('integer','numeric')){
    n22 <- total_possible_target-length(union(t1,t2))
  }else{
    n22 <- length(setdiff(total_possible_target,c(t1,t2)))
  }
  mm  <- cbind(c(n11,n21),c(n12,n22))
  ft  <- fisher.test(mm)$p.value
  or  <- n11/n12/(n21/n22)
  rownames(mm) <- c(sprintf('In %s target',source1_label),sprintf('Not in %s target',source1_label))
  colnames(mm) <- c(sprintf('In %s target',source2_label),sprintf('Not in %s target',source2_label))
  print(mm)
  res <- c('P.Value'=ft,'Odds_Ratio'=or,'Intersected_Number'=mm[1,1])
  return(res)
}
##
####

# plot_pattern: no, random_one, all
get.spByGene <- function(igraph_obj=NULL,driver_list=NULL,target_list=NULL,mode='all',transfer_tab=NULL,plot_pattern='random_one',...){
  ori_driver_list <- driver_list
  ori_target_list <- target_list
  if(is.null(driver_list) & length(target_list)==1){
    driver_list <- unlist(lapply(neighborhood(igraph_obj,nodes=target_list,mode='in'),names))
  }
  if(is.null(target_list) & length(driver_list)==1){
    target_list <- unlist(lapply(neighborhood(igraph_obj,nodes=driver_list,mode='out'),names))
  }
  driver_list <- unique(intersect(driver_list,names(V(igraph_obj))))
  target_list <- unique(intersect(target_list,names(V(igraph_obj))))
  if(length(driver_list)==0 & length(target_list)==0){
    message('No gene included, please check and re-try!');return(FALSE)
  }
  #r1 <- all_shortest_paths(igraph_obj,from=use_list,to=use_list,mode=mode)$res
  #r2 <- unique(unlist(lapply(r1,names)))
  print(driver_list)
  print(target_list)

  if(length(driver_list)==0) driver_list <- target_list
  if(length(target_list)==0) target_list <- driver_list

  all_r1 <- list()
  for(d in driver_list){
    for(t in target_list){
      if(d!=t) all_r1[[sprintf('%s %s',d,t)]] <- all_shortest_paths(igraph_obj,from=d,to=t,mode=mode)$res
    }
  }
  #print(all_r1)

  if(plot_pattern!='all') r2 <- unique(unlist(lapply(all_r1,function(x)names(x[[1]]))))
  if(plot_pattern=='all') r2 <- unique(unlist(lapply(all_r1,function(x)names(unlist(x)))))

  all_sp <- induced_subgraph(igraph_obj,r2)

  ## get layout matrix,x,y
  net <- all_sp
  if(plot_pattern!='no'){
    plot.spByGene(net,driver_list=ori_driver_list,target_list=ori_target_list,transfer_tab=transfer_tab,...)
  }
  return(list(all_sp=all_r1,use_sp=net))
}

### plot igraph
plot.spByGene <- function(net,driver_list=NULL,target_list=NULL,transfer_tab=NULL,vertex.size=25,vertex.label.cex=1,
                          vertex.frame.color="white",vertex.label.color="black",edge.color='black',
                          layout='auto',...){
  n1 <- names(V(net))
  cc <- brewer.pal(11,'Set3')
  c1 <- ifelse(n1 %in% driver_list,cc[3],ifelse(n1 %in% target_list,cc[1],cc[9]))
  names(c1) <- n1
  d1 <- degree(net,mode='out')
  n2 <- n1[order(d1,decreasing = TRUE)]
  y <- ifelse(d1[n2]>0,2,1)
  x <- c(seq(0,1,length.out=length(y[which(y==2)])),seq(0,1,length.out=length(y[which(y==1)])))
  l <- cbind(x,y)
  use_label <- names(V(net))
  if(is.null(transfer_tab)==FALSE){
    transfer_tab <- as.matrix(transfer_tab)
    tt1 <- rbind(cbind(paste0(transfer_tab[,1],'_TF'),paste0(transfer_tab[,2],'_TF')),
                 cbind(paste0(transfer_tab[,1],'_SIG'),paste0(transfer_tab[,2],'_SIG')),
                 transfer_tab[,1:2])
    tt1 <- unique(tt1)
    rownames(tt1) <- tt1[,1]
    w1 <- which(use_label %in% rownames(tt1))
    use_label[w1] <- tt1[use_label[w1],2]
  }
  if(class(layout) != 'function'){
    plot.igraph(net,vertex.label=use_label,vertex.color=c1,layout=l,vertex.size=vertex.size,vertex.frame.color=vertex.frame.color,
                vertex.label.color=vertex.label.color,vertex.label.cex=vertex.label.cex,
                edge.color=edge.color,...)
  }else{
    plot.igraph(net,vertex.label=use_label,vertex.color=c1,layout=layout,vertex.size=vertex.size,vertex.frame.color=vertex.frame.color,
                vertex.label.color=vertex.label.color,vertex.label.cex=vertex.label.cex,
                edge.color=edge.color,...)
  }
  return(TRUE)
}



#' get target list by input network data
#' @export
get_net2target_list <- function(net_dat=NULL) {
  all_source <- unique(net_dat$source)
  all_target <- lapply(all_source, function(x) {
    n1 <- net_dat[which(net_dat$source == x), c('target', 'MI', 'spearman')]
    rownames(n1) <- n1$target
    return(n1)
  })
  names(all_target) <- all_source
  return(all_target)
}
