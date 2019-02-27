##Function based MINIE
##Author:chenxi.qian@stjude.org
##Stjude.YuLab

###Function 1:read MICA result & MICA input, output as a eset###
#' @export
readMICAoutput<-function(input_file, output_file,load_clust_label=TRUE){

	cat("Reading Input...","\n")
	d <- read.table(input_file,header = TRUE,stringsAsFactors = FALSE,quote = "")
	input<-t(d[,-1]);colnames(input)<-d$ID;rm(d)

	MICA_Result <- read.table( output_file, # MICA output text file
                          	   header = TRUE,
                          	   stringsAsFactors = FALSE)

	fd<-data.frame(row.names=rownames(input),
					geneNames=rownames(input),stringsAsFactors = FALSE)

	pd<-data.frame(row.names=colnames(input),
					cellNames=colnames(input),
					tSNE_1=MICA_Result[,2],
					tSNE_2=MICA_Result[,3],stringsAsFactors=FALSE)

	if(load_clust_label){pd$label <- MICA_Result$label;cat("Clustering info is under 'label' slot.","\n")}

	eset<-new("ExpressionSet",phenoData= new("AnnotatedDataFrame",pd),
          featureData=new("AnnotatedDataFrame",fd), annotation="",exprs=as.matrix(input))
  cat("ExpressionSet Generated!","\n")
	return(eset)
}



###Network related functions###
###Funciton2: Generate SJAaracne input using scRNAseq data###
#' @export
SJARACNeInput_scRNAseq<-function(eset.sel,tf.ref,wd.src,grp.tag){

  cat(grp.tag,'\n')

  #fData(eset.sel)$IQR<-apply(exprs(eset.sel),1,IQR)

  # exclude genes with all zero
  eset.sel<-eset.sel[apply(exprs(eset.sel),1,function(xx){sum(xx)!=0}),]

  fData(eset.sel)$geneSymbol<-fData(eset.sel)$geneNames

  ni<-nrow(eset.sel);ni
  ns<-ncol(eset.sel);ns
  ng<-nlevels(factor(fData(eset.sel)$geneNames));ng
  tag<-paste(grp.tag,nrow(eset.sel),ng,ncol(eset.sel),sep='_');tag

  dir.cur<-file.path(wd.src,tag);dir.cur
  dir.create(dir.cur,recursive = T)
  dir.create(file.path(dir.cur,'tf'),recursive = T)

  #write exp data to exp format
  expdata<-data.frame(cbind(isoformId=featureNames(eset.sel),geneSymbol=fData(eset.sel)$geneSymbol,exprs(eset.sel)))
  f.exp<-file.path(dir.cur,paste(grp.tag,"_",ni,"_",ng,"_",ns,".exp",sep=''));f.exp
  write.table(expdata,file=f.exp,sep="\t",row.names=FALSE,quote=FALSE)

  tf.eset.sel<-subset(eset.sel,fData(eset.sel)$geneSymbol%in%tf.ref)
  dim(tf.eset.sel)
  f.tf<-file.path(dir.cur,'tf',paste(grp.tag,"_",nrow(tf.eset.sel),"_",nlevels(factor(fData(tf.eset.sel)$geneSymbol)),"_",ns,"_tf.txt",sep=''));f.tf
  cat(featureNames(tf.eset.sel),file=f.tf,sep='\n')

}

###Function8: Wrap up function for generate SJARACNe input###
#' @export
generateSJARACNeInput<-function(eset,tf.ref,wd.src,group_tag){

  if (!dir.exists(wd.src)) dir.create(wd.src,recursive = T)

  if (group_tag%in%colnames(pData(eset))){
    groups <- unique(pData(eset)[,"group_tag"])
    for (i in 1:length(groups)){
      grp.tag<-groups[i]; eset[,which(pData(eset)[,"group_tag"]==grp.tag)] -> eset.sel
      SJARACNeInput_scRNAseq(eset.sel=eset.sel,tf.ref=tf.ref,wd.src=wd.src,grp.tag=grp.tag)
    }#end for
  }else{
    stop("Lack of group info, please check your group_tag.","\n")
  }#end if
}#end function





###scRNA-seq data preprocess(before running clustering)
#' @export
pre.MICA <- function(data.input=NULL, #data matrix that have unique colnames and geneSymbol as rownames
            gene_filter=TRUE, 
            cell_filter=TRUE,
            cell_percentage=0.005, 
            ERCC_filter=TRUE,
            Mito_filter=TRUE, # one way filtering
            UMI_filter="both", #three way filtering 
            plotting=TRUE,
            plot.dir=".",
            norm=10e6, 
            sampleID="PBMC_12k",
            log2Transform=TRUE)
{
  cat("Pre-QC expression matrix dimention: ", dim(data.input),"\n")
  # 
  cells_per_gene <- rowSums(sign(as.matrix(data.input)))
  nGene <- sum(cells_per_gene > 0)
  cat("# of non-zero gene:", nGene, "\n")

  # Count the cells with >=1 identified gene(s)
  genes_per_cell <- colSums(sign(as.matrix(data.input)))
  nCell <- sum(genes_per_cell > 0)
  cat("# of non-zero cell:", nCell, "\n")

  #remove genes by percentage

  nCell_cutoff <- floor(cell_percentage * dim(data.input)[2])
  cut.gene <- names(which(cells_per_gene < nCell_cutoff))


  p.gene.qc <- ggplot(data=data.frame(ncells=cells_per_gene), aes(log10(ncells+1))) + 
      geom_histogram(bins = 100) + 
      labs(x = "LOG10 (Cell counts + 1)", 
           y = "Gene counts", 
           title = "Histogram: Number of cells expressed by each genes") + 
      theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 12, face = "bold"), 
          plot.title = element_text(size = 10, face = "bold"))
 

  data <- data.input[!rownames(data.input) %in% cut.gene, ]
  
  cat("Gene expressed in less than ", nCell_cutoff,
      "(",length(cut.gene)*100/dim(data.input)[1],"%) cells were filtered","\n",
      "Filtered expression matrix dimension:",dim(data),"\n")
  
  # extract mito-gene and spike-in genes
  mito.genes <- grep(pattern = "^mt-|^MT-", x = rownames(data), value = TRUE)
  spikeIn.genes <- grep(pattern = "^ERCC-", x = rownames(data), value = TRUE)
 
  pd <- data.frame(nUMI.total = colSums(data),
             nGene = genes_per_cell,
             percent.mito= Matrix::colSums(data[mito.genes, ]) / Matrix::colSums(data),
             percent.spikeIn = Matrix::colSums(data[spikeIn.genes, ]) / Matrix::colSums(data),
             sampleID=sampleID,
             stringsAsFactors = FALSE)

  #Compute filtering criteria
  umi_cf_lo <- floor(exp(median(log(pd$nUMI.total)) - 3 * mad(log(pd$nUMI.total))))
  umi_cf_hi <- floor(exp(median(log(pd$nUMI.total)) + 3 * mad(log(pd$nUMI.total))))
  ERCC_cf <- median(pd$percent.ERCC) + 3 * mad(pd$percent.ERCC)
  mito_cf <- median(pd$percent.mito) + 3 * mad(pd$percent.mito)


  if(UMI_filter=="low") {umi_cf_hi<-Inf
    }else if(UMI_filter=="high") {umi_cf_lo <- 0}
  
  cut.cell <- which((pd$nUMI.total < umi_cf_lo) | (pd$nUMI.total > umi_cf_hi))
  
  if(ERCC_filter) cut.cell <- union(cut.cell,
                  which(pd$percent.ERCC > ERCC_cf))

  if(Mito_filter) cut.cell <- union(cut.cell,
                  which(pd$percent.mito > mito_cf))


  # Visulize the distribution of Identified Genes, UMIs and %Mitochondrial genes among cells
  p1 <- ggplot(pd,aes(x=sampleID,y=nUMI.total,fill=sampleID)) + 
      geom_violin(scale="width", na.rm=TRUE)+
      geom_boxplot(width=0.3,outlier.size = 0.25, outlier.stroke = 0.5,fill="white")+
      labs(x=" ",y="",title="Total UMI count for each cell") + 
      theme(legend.position="none",plot.title=element_text(size=8))
  
  p2 <- ggplot(pd,aes(x=sampleID,y=nGene,fill=sampleID)) + 
      geom_violin(scale="width", na.rm=TRUE)+
      geom_boxplot(width=0.3,outlier.size = 0.25, outlier.stroke = 0.5,fill="white")+
      labs(x=" ",y="",title="Total # of gene expressed in each cell") + 
      theme(legend.position="none",plot.title=element_text(size=8))
  
  plot.list <- list(p1,p2)

  if(length(spikeIn.genes)!=0){
  p3 <- ggplot(data=pd, aes(x=nUMI.total,y=percent.ERCC))+ geom_point(na.rm=TRUE)+
      labs(x=" ",y="",title="Percentage of spike-in gene expression VS Total UMI counts in each cell") + 
      theme(legend.position="none",plot.title=element_text(size=8))
  
  plot.list$p3<-p3}
 
  if(length(mito.genes)!=0){
  p4 <- ggplot(data=pd, aes(x=nUMI.total,y=percent.mito))+ geom_point(na.rm=TRUE)+
      labs(x=" ",y="",title="Percentage of MT-gene expression VS Total UMI counts in each cell") + 
      theme(legend.position="none",plot.title=element_text(size=8))
  
  plot.list$p4<-p4}

  p.cell.qc <- cowplot::plot_grid(plotlist=plot.list,nrow=1)
  
 ##############################
 # do cell filtering 
  data<-data[,-c(cut.cell)]
  pd<-pd[data@Dimnames[[2]], ] 
 ###############################

  cat("A total of ",length(cut.cell), "(",length(cut.cell)*100/dim(data.input)[1],"%) cells were filtered","\n",
      "Filtered expression matrix dimension:",dim(data),"\n")
  cat("==============","\n")


  if(plotting){
  
  ggsave(plot = p.gene.qc, filename = "1_1_Gene_QCmetrics_before_filtering.png", 
        width = 5, height = 5, units = "in", dpi = 300)
  dev.off()
  
  png("1_2_Cell_QCmetrics_before_filtering.png",width = 14, height = 6, units = "in", res = 300)
  print(p.cell.qc)
  dev.off()

  message("Plots are generated and saved at",normalizePath(plot.dir),"\n")
  }

  # Normalization and log2 transformation #########
  # CPM 100k
  cellSum <- colSums(data)
  if(!is.null(norm)) {
    data <- sweep(data, 2, norm/unname(pd$nUMI.total), '*');
    cat("Data was normalized!","\n")}

  # log transformation
  if(log2Transform) {
    data <- log( data + 1, base=2)
    cat("Data was log2Transformed!","\n")}

  return(data)
}





