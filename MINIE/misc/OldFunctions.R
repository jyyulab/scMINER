

#' @title pre.MICA
#' @description This function helps to conduct scRNA-seq data preprocessing and quality
#' control for MICA input.
#' @param raw_data input matrix
#' @param projectName project name used for Rmarkdown report title
#' @param sampleID a string or a vector to indicate sample info
#' @param output_rmd logical, whether or not output Rmd QC report
#' @param plot.dir output path for QC report
#' @param gene_filter logical
#' @param cell_filter logical
#' @param cell_percentage numerical, default as 0.005
#' @param ERCC_filter logical
#' @param Mito_filter logical
#' @param nGene_filter logical
#' @param nUMI_filter logical
#' @param plotting logical
#' @param norm numerical, default to 10e6
#' @param logTransform logical,default as TRUE
#' @param base numerical, log(n+1,base=base)
#'
#' @export
pre.MICA <- function(raw_data=NULL,  #data matrix that have unique colnames and geneSymbol as rownames
                     projectName="SAMPLE",
                     sampleID="SAMPLE",
                     output_rmd=TRUE,
                     plot.dir=".",
                     gene_filter=TRUE,
                     cell_filter=TRUE,
                     cell_percentage=0.005,
                     ERCC_filter=TRUE,
                     Mito_filter=TRUE, # one way filtering
                     nGene_filter=TRUE,
                     nUMI_filter="both", #three way filtering
                     plotting=TRUE,
                     norm=10e6,
                     logTransform=TRUE,
                     base=NULL
)
{
  if(!class(raw_data)[1]%in%c("dgCMatrix","dgTMatrix","matrix")){
    stop("Input format should %in% c( 'matrix','dgTMatrix','dgCMatrix')","\n")
  }

  cat("Running QC...","\n",
      "Pre-QC expression matrix dimention: ", dim(raw_data),"\n")

  if(!dir.exists(plot.dir)) {dir.create(plot.dir)}
  if(output_rmd) {render(input=system.file("rmd", "Preprocessing.Rmd", package = "MINIE"),
                         output_dir = plot.dir,
                         output_file = paste0(projectName,"_scRNAseq_preprocessing.html"),
                         clean=TRUE,
                         quiet =TRUE,
                         params=list(
                           d=raw_data,
                           gene_filter=gene_filter,
                           projectName=projectName,
                           cell_filter=cell_filter,
                           cell_percentage=cell_percentage,
                           ERCC_filter=ERCC_filter,
                           Mito_filter=Mito_filter,
                           nGene_filter=nGene_filter,
                           nUMI_filter=nUMI_filter,
                           plot.dir=plot.dir,
                           norm=norm,
                           sampleID=sampleID,
                           logTransform=logTransform,
                           base=base))

    cat("Running QC...","\n",
        "After QC expression matrix dimention: ", dim(data),"\n")
     }

  else{
    d<-raw_data;rm(raw_data)

    cells_per_gene <- Matrix::rowSums(d!=0)
    nGene <- sum(cells_per_gene > 0)
    cat("# of non-zero gene:", nGene, "\n")

    # Count the cells with >=1 identified gene(s)
    genes_per_cell <- Matrix::colSums(d!=0)
    nCell <- sum(genes_per_cell > 0)
    cat("# of non-zero cell:", nCell, "\n")

    #Compute filtering criteria for genes
    nCell_cutoff <- max(floor(cell_percentage * dim(d)[2]),1) #make sure that all 0 genes are filtered

    gene <- unname((which(cells_per_gene >= nCell_cutoff)))

    d.tmp<-d[gene,]#use index instead of gene names to avoid mis-select problem

    cat("Defining meta data...","\n")
    # extract mito-gene and spike-in genes
    mito.genes <- grep(pattern = "^mt-|^MT-", x = rownames(d.tmp), value = TRUE)
    spikeIn.genes <- grep(pattern = "^ERCC-|^Ercc", x = rownames(d.tmp), value = TRUE)
    if(length(mito.genes)==0) Mito_filter=FALSE
    if(length(spikeIn.genes)==0) ERCC_filter=FALSE

    pd <- data.frame(nUMI.total = Matrix::colSums(d),
                     nGene = Matrix::colSums(sign(d)),  #need adjustment!
                     percent.mito = round(Matrix::colSums(d.tmp[mito.genes, ]) / Matrix::colSums(d.tmp),8),
                     percent.spikeIn = round(Matrix::colSums(d.tmp[spikeIn.genes, ]) / Matrix::colSums(d.tmp),8),
                     sampleID=sampleID,
                     stringsAsFactors = FALSE)

    #Compute filtering criteria for cells
    umi_cf_lo <- max(floor(exp(median(log(pd$nUMI.total)) - 3 * mad(log(pd$nUMI.total)))),100)
    umi_cf_hi <- ceiling(exp(median(log(pd$nUMI.total)) + 3 * mad(log(pd$nUMI.total))))
    nGene_cf <- max(floor(exp(median(log(pd$nGene)) - 3 * mad(log(pd$nGene)))),50)
    ERCC_cf <- round(median(pd$percent.spikeIn) + 3 * mad(pd$percent.spikeIn),4)
    mito_cf <- round(median(pd$percent.mito) + 3 * mad(pd$percent.mito),4)


    if(plotting){

      p.gene.qc <- ggplot(data = data.frame(ncells = cells_per_gene), aes(log10(ncells+1))) +
        geom_histogram(bins = 100) +
        labs(x = "Log10 (Cell counts + 1)",
             y = "Gene counts",
             title = "Histogram: Number of cells expressed by each genes") +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"),
              plot.title = element_text(size = 10, face = "bold"))+
        geom_vline(xintercept = log10 (nCell_cutoff + 1), data = as.data.frame(nCell_cutoff),
                   color="blue",size = 1,alpha = 0.3)

      ggsave(plot = p.gene.qc, filename = file.path(plot.dir,"1_1_Gene_QCmetrics_before_filtering.png"),
             width = 5, height = 4, units = "in", dpi = 300)
      dev.off()


      # Visulize the distribution of Identified Genes, UMIs and %Mitochondrial genes among cells
      p1 <- ggplot(pd,aes(x=sampleID,y=nUMI.total,fill=sampleID)) +
        geom_violin(scale="width", na.rm=TRUE)+
        geom_boxplot(width=0.3,outlier.size = 0.25, outlier.stroke = 0.5,fill="white")+
        labs(x=" ",y="",title="Total UMI count for each cell") +
        theme(legend.position="none",plot.title=element_text(size=8))+
        geom_hline(yintercept = umi_cf_lo, data = as.data.frame(umi_cf_lo),
                   color="blue",size=2,alpha = 0.3) +
        geom_hline(yintercept = umi_cf_hi, data = as.data.frame(umi_cf_hi),
                   color="blue",size=2,alpha = 0.3)


      p2 <- ggplot(pd,aes(x=sampleID,y=nGene,fill=sampleID)) +
        geom_violin(scale="width", na.rm=TRUE)+
        geom_boxplot(width=0.3,outlier.size = 0.25, outlier.stroke = 0.5,fill="white")+
        labs(x=" ",y="",title="Total # of gene expressed in each cell") +
        theme(legend.position="none",plot.title=element_text(size=8))+
        geom_hline(yintercept = nGene_cf, data = as.data.frame(nGene_cf),
                   color="blue",size=2,alpha = 0.3)

      plot.list <- list(p1,p2)

      if(length(spikeIn.genes)!=0){

        p3 <- ggplot(data=pd, aes(x=nUMI.total,y=percent.spikeIn))+
          geom_point(na.rm=TRUE,size=0.8,alpha=0.5)+
          labs(x="",y="",title="Percentage of spike-in gene expression VS Total UMI counts in each cell") +
          theme(legend.position="none",plot.title=element_text(size=8))+
          geom_hline(yintercept = ERCC_cf,data=as.data.frame(ERCC_cf),
                     size=2,color="blue",alpha = 0.3)
        plot.list$p3 <- p3}

      if(length(mito.genes)!=0){
        p4 <- ggplot(data=pd, aes(x=nUMI.total,y=percent.mito))+
          geom_point(na.rm=TRUE,size=0.8,alpha=0.5)+
          labs(x="",y="",title="Percentage of MT-gene expression VS Total UMI counts in each cell") +
          theme(legend.position="none",plot.title=element_text(size=8))+
          geom_hline(yintercept = mito_cf, data=as.data.frame(mito_cf),
                     size=2,color="blue",alpha = 0.3)
        plot.list$p4 <- p4}

      p.cell.qc <- cowplot::plot_grid(plotlist=plot.list,nrow=1)

      png(file.path(plot.dir,"1_2_Cell_QCmetrics_before_filtering.png"),width = 14, height = 6, units = "in", res = 300)
      print(p.cell.qc)
      dev.off()

      message("Plots are generated and saved at",normalizePath(plot.dir),"\n")
    }


    if(nUMI_filter=="low") {umi_cf_hi = Inf
    }else if(nUMI_filter=="high") {umi_cf_lo = 0}
    cell <- which((pd$nUMI.total > umi_cf_lo) & (pd$nUMI.total < umi_cf_hi))
    if(nGene_filter) cell <- intersect(cell, which(pd$nGene > nGene_cf))
    if(ERCC_filter) cell <- intersect(cell, which(pd$percent.spikeIn < ERCC_cf))
    if(Mito_filter) cell <- intersect(cell, which(pd$percent.mito < mito_cf))

    #gene filtering
    if(gene_filter){
      data <- d.tmp
      cat("Gene expressed in less than ", nCell_cutoff,
          "cells (",(dim(d)[1]-length(gene))*100/dim(d)[1],"% genes) were filtered","\n",
          "Filtered expression matrix dimension:",dim(data),"\n")
    }else{data <- d}

    #cell filtering
    if(cell_filter){
      data <- data[,cell]
      pd <- pd[colnames(data), ]

      cat("A total of ",dim(d)[2]-length(cell),
          "(",(dim(d)[2]-length(cell))*100/dim(d)[2],"%) cells were filtered","\n",
          "Filtered expression matrix dimension:",dim(data),"\n")}

    cat("==============","\n")

    # Normalization and log2 transformation
    # CPM 100k
    cellSum <- Matrix::colSums(data)
    if(!is.null(norm)) {
      data <- sweep(data, 2, norm/unname(pd$nUMI.total), '*');
      cat("Data was normalized!","\n")}

    # log transformation
    if(logTransform) {
      if(is.null(base)) data <- log( data + 1)
      else {data <- log(data+1,base=base)}
      cat(paste0("Data was log",base,"Transformed!"),"\n")}
  }
  return(data)
}#end pre.MICA


preMICA.QA<-function(raw.data=NULL,
                       feature.data=NULL,
                       meta.data=NULL,
                       project.name="SAMPLE",
                       sample.group="SAMPLE",
                       output_rmd=TRUE,
                       plot.dir="./QC/",
                       cell_percentage=0.005,
                       plotting=TRUE)
{
  if(!class(raw.data)[1]%in%c("dgCMatrix","dgTMatrix","matrix")){
    stop("Input format should %in% c( 'matrix','dgTMatrix','dgCMatrix')","\n")
  }

  if(is.null(feature.data)){
    warning("Will take rownames as geneSymbol","\n")
  }else if(any(rownames(raw.data)!=rownames(feature.data))){
      stop("Row names of feature data doesn't match with expression data row names!","\n")
  }

  if(!is.null(meta.data)){
    if(any(colnames(raw.data)!=rownames(meta.data))){
      stop("Row names of meta data doesnt match with expression data Column names!","\n")
    }
  }

  if(!is.null(sample.group)){
    if(length(sample.group)!=dim(raw.data)[2]&length(sample.group)!=1){
      stop("length of sample.group information should be 1 or equals to number of sample!","\n")
    }
  }

  if(any(duplicated(row.names(raw.data)))){
    cat("Found duplicated rownames, convert row.names of raw.data to an arbitrary vector!","\n")
    cat("Old feature name stored in feature.data$old.name..","\n")
    if(is.null(feature.data)){
      feature.data<-data.frame(old.name=rownames(raw.data))
    }else{
      feature.data$old.name=rownames(raw.data)
    }
    row.names(raw.data)<- NULL
  }

  cat("Passed sanity check..","\n",
      "Running Quality Assessment...","\n",
      "Expression matrix dimention: ", dim(raw.data),"\n")

  if(!dir.exists(plot.dir)) {dir.create(plot.dir)}

  ###DO QA
  if(output_rmd) {
    render(input=system.file("rmd", "Preprocessing_QA.Rmd", package = "scMINER"),
           output_dir = plot.dir,
           output_file = paste0(project.name,"_scRNAseq_QA.html"),
           clean=TRUE,
           quiet =TRUE,
           params=list(
             d=raw.data,
             feature.data=feature.data,
             meta.data=meta.data,
             projectName=project.name,
             cell_percentage=cell_percentage,
             plot.dir=plot.dir,
             sampleID=sample.group))

    cat("QA Done!","\n")
  }
  #without rmarkdown
  else{
    d<-raw.data;rm(raw.data)
    cells_per_gene <- Matrix::rowSums(d!=0)

    if(is.null(feature.data)){
      feature.data<-data.frame(nCells=cells_per_gene,
                               geneSymbol=rownames(d),
                               stringsAsFactors = FALSE)}

    else feature.data$nCells<-cells_per_gene

    nGene <- sum(cells_per_gene > 0)
    cat("# of non-zero gene:", nGene, "\n")

    # Count the cells with >=1 identified gene(s)
    genes_per_cell <- Matrix::colSums(d!=0)
    nCell <- sum(genes_per_cell > 0)
    cat("# of non-zero cell:", nCell, "\n")

    #Compute filtering criteria for genes
    nCell_cutoff <- max(floor(cell_percentage * dim(d)[2]),1) #make sure that all 0 genes are filtered
    gene <- unname((which(cells_per_gene >= nCell_cutoff)))

    d.tmp<-d[gene,]#use index instead of gene names to avoid mis-select problem

    cat("Defining meta data...","\n")
    # extract mito-gene and spike-in genes

    if ("geneSymbol" %in% colnames(feature.data)){
      mito.genes <- grep(pattern = "^mt-|^MT-", x = feature.data$geneSymbol)
      spikeIn.genes <- grep(pattern = "^ERCC-|^Ercc", x = feature.data$geneSymbol)
    }else{
      cat("Mitochondrial genes and spike-in genes were not found due to lack of geneSymbol information...","\n")
      mito.genes<-integer(0)
      spikeIn.genes<-integer(0)
    }

    if(length(mito.genes)==0) Mito_filter=FALSE
    if(length(spikeIn.genes)==0) ERCC_filter=FALSE

    pd <- data.frame(row.names = colnames(d),
                     nUMI.total = Matrix::colSums(d),
                     nGene = Matrix::colSums(sign(d)),
                     percent.mito = round(Matrix::colSums(d[mito.genes, ]) / Matrix::colSums(d),8),
                     percent.spikeIn = round(Matrix::colSums(d[spikeIn.genes, ]) / Matrix::colSums(d),8),
                     group=sample.group,
                     stringsAsFactors = FALSE)


    #Compute filtering criteria for cells
    umi_cf_lo <- max(floor(exp(median(log(pd$nUMI.total)) - 3 * mad(log(pd$nUMI.total)))),100)
    umi_cf_hi <- ceiling(exp(median(log(pd$nUMI.total)) + 3 * mad(log(pd$nUMI.total))))
    nGene_cf <- max(floor(exp(median(log(pd$nGene)) - 3 * mad(log(pd$nGene)))),50)
    ERCC_cf <- round(median(pd$percent.spikeIn) + 3 * mad(pd$percent.spikeIn),4)
    mito_cf <- round(median(pd$percent.mito) + 3 * mad(pd$percent.mito),4)

    if(plotting){

      p.gene.qc <- ggplot(data = data.frame(ncells = cells_per_gene), aes(log10(ncells+1))) +
        theme_classic()+
        geom_histogram(bins = 100) +
        labs(x = "Log10 (Cell counts + 1)",
             y = "Gene counts",
             title = "Histogram: Number of cells expressed by each genes") +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"),
              plot.title = element_text(size = 10, face = "bold"))+
        geom_vline(xintercept = log10 (nCell_cutoff + 1), data = as.data.frame(nCell_cutoff),
                   color="blue",size = 1,alpha = 0.3)

      ggsave(plot = p.gene.qc, filename = file.path(plot.dir,"1_1_Gene_QCmetrics_before_filtering.png"),
             width = 10, height = 6, units = "in", dpi = 300)

      # Visulize the distribution of Identified Genes, UMIs and %Mitochondrial genes among cells
      p1 <- ggplot(pd,aes(x=group,y=nUMI.total,fill=group)) +
        theme_classic()+
        geom_violin(scale="width", na.rm=TRUE)+
        geom_boxplot(width=0.3,outlier.size = 0.25, outlier.stroke = 0.5,fill="white")+
        labs(x=" ",y="nUMI",title="") +
        theme(legend.position="none",plot.title=element_text(size=8))+
        geom_hline(yintercept = umi_cf_lo, data = as.data.frame(umi_cf_lo),
                   color="blue",size=2,alpha = 0.3) +
        geom_text(aes(0.5,umi_cf_lo,label = umi_cf_lo, vjust = -1),size=3)+
        geom_hline(yintercept = umi_cf_hi, data = as.data.frame(umi_cf_hi),
                   color="blue",size=2,alpha = 0.3) +
        geom_text(aes(0.5,umi_cf_hi,label = umi_cf_hi, vjust = -1),size=3)

      p2 <- ggplot(pd,aes(x=group,y=nGene,fill=group)) +
        geom_violin(scale="width", na.rm=TRUE)+
        theme_classic()+
        geom_boxplot(width=0.3,outlier.size = 0.25, outlier.stroke = 0.5,fill="white")+
        labs(x=" ",y="nGene",title="") +
        theme(legend.position="none",plot.title=element_text(size=8))+
        geom_hline(yintercept = nGene_cf, data = as.data.frame(nGene_cf),
                   color="blue",size=2,alpha = 0.3)+
        geom_text(aes(0.5,nGene_cf,label = nGene_cf, vjust = -1),size=3)

      p.cell.1 <- cowplot::plot_grid(p1,p2,nrow=1,ncol=2)
      png(file.path(plot.dir,"2_1_Distribution_of_nGene_and_nUMI_preQC.png"),width = 10, height = 8, units = "in", res = 300)
      print(p.cell.1)
      dev.off()

      if(length(spikeIn.genes)!=0){
        p3 <- ggplot(data=pd, aes(x=nUMI.total,y=percent.spikeIn))+
          theme_classic()+
          geom_point(na.rm=TRUE,size=0.8,alpha=0.5)+
          labs(x="nUMI",y="ERCC fraction",
               title="") +
          theme(legend.position="none",plot.title=element_text(size=8))+
          geom_hline(yintercept = ERCC_cf,data=as.data.frame(ERCC_cf),
                     size=2,color="blue",alpha = 0.3)+
          geom_text(aes(0.5,ERCC_cf,label = ERCC_cf, vjust = -1,hjust=-6),size=3)
      }else {
        text<-"\n  No Spike-in genes detected. \n"
        p3 <- ggplot() +
              theme_void() +
              annotate("text", x = 4, y = 25, size=5, label = text)+
              labs(x="",y="")}

      if(length(mito.genes)!=0){
        p4 <- ggplot(data=pd, aes(x=nUMI.total,y=percent.mito))+
          theme_classic()+
          geom_point(na.rm=TRUE,size=0.8,alpha=0.5)+
          labs(x="nUMI",y="MT-gene fraction",
               title="") +
          theme(legend.position="none",plot.title=element_text(size=8))+
          geom_hline(yintercept = mito_cf, data=as.data.frame(mito_cf),
                     size=2,color="blue",alpha = 0.3)+
          geom_text(aes(0.5,mito_cf,label = mito_cf, vjust = -1,hjust=-6),size=3)
      }else{
        text<-"\n  No Mitochondrial genes detected. \n"
        p4 <- ggplot() +
              theme_void() +
              annotate("text", x = 4, y = 25, size=5, label = text)+
              labs(x="",y="")}

      p.cell.2 <- cowplot::plot_grid(p3,p4,nrow=1,ncol=2)
      png(file.path(plot.dir,"2_2_Fraction_of_MT_and_SpikeIn_VS_nUMI_preQC.png"),width = 10, height = 6, units = "in", res = 300)
      print(p.cell.2)
      dev.off()

      message("Plots are generated and saved at", normalizePath(plot.dir),"\n")
    }

    if(!is.null(meta.data)) pd<-cbind(pd,meta.data)
    Obj<-list(raw.data=d,
              meta.data=pd,
              feature.data=feature.data,
              plot.dir=normalizePath(plot.dir),
              cal.cutoffs=list(nCell_cutoff = nCell_cutoff,
                               umi_cf_lo=umi_cf_lo,
                               umi_cf_hi=umi_cf_hi,
                               nGene_cf = nGene_cf,
                               ERCC_cf = ERCC_cf,
                               mito_cf = mito_cf))
    cat("Computed cut-off:","\n")
    print(Obj$cal.cutoffs)
  }
  return(Obj)
}
