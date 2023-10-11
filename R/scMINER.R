## ---------------------------
##
## Script name: scMINER v1.0
##
## Purpose of script: designed for preprocessing, QC, clustering,
##    and hidden driver analysis of single-cell RNA-seq data
##
## Author: The scMINER software is developed and maintained by the Yu Laboratory @ St. Jude
##
## Date Created: 2023-03-24
##
## ---------------------------
##
## Notes: scMINER enables mutual information-based cell clustering,
##    cell-type-specific gene regulatory network (GRN) reverse engineering and
##    protein activity inference, to identify hidden transcriptional factors (TFs)
##    and signaling factors (SIGs) driving cellular lineage differentiation and
##    tissue specific specification.
##
##    scMINER software consists of three components:
##
##      (1) MICA: Mutual Information based Clustering analysis (https://github.com/jyyulab/MICA).
##      (2) SJARACNe:a scalable solution of ARACNe that improves the computational performance to reconstruct the regulatory network (https://github.com/jyyulab/SJARACNe).
##      (3) MINIE: Mutual Information-based Network Inference Engine (MINIE.R).
##
##    scMINER provides supporting functions to prepare datasets for running MICA and SJARACNe and import results into R, including following functions:
##    scMINER.dir.create()
##    CreateSparseEset()
##    readscRNAseqData()
##    draw.scRNAseq.QC()
##    draw.group.barplot()
##    draw.marker.bbp()
##    feature_vlnplot()
##    feature_heatmap()
##    feature_highlighting()
##    preMICA.filtering()
##    generateMICAinput()
##    readMICAoutput()
##    MICAplot()
##    generateSJARACNeInput()
## ---------------------------

#' Manipulation of Working Directories for scMINER pipeline
#'
#' \code{scMINER.dir.create} is used to help users create an organized working directory for the network construction step in scMINER analysis.
#' However, it is not essential for the analysis.
#' It creates a hierarchical working directory and returns a list contains this directory information.
#'
#' This function needs users to define the main working directory and the project's name.
#' It creates a main working directory with a subdirectory of the project.
#' It also automatically creates five subfolders (DATA, SJAR, MICA, QC, PLOT) within the project folder.
#' DATA/, storing data files;
#' SJAR/, storing files needed for running SJAracne;
#' MICA/, storing files needed for running MICA;
#' QC/,   storing Quality Control related plots;
#' PLOT/, storing plot files;
#' This function also returns a list object (example, \code{scminer.par} in the demo) with directory information wrapped inside.
#' @param project_main_dir character, name or absolute path of the main working directory.
#' @param project_name character, name of the project folder.
#'
#' @return \code{scMINER.dir.create} returns a list object, containing main.dir (path of the main working directory),
#' project.name (project name), out.dir (path of the project folder).
#'
#' @examples
#' project_main_dir <- '../test'
#' project_name <- 'PBMC14KDS'
#' scminer.par <- scMINER::scMINER.dir.create(project_main_dir = project_main_dir,
#'                                         project_name = project_name)
#' \dontrun{
#' # Creating a main working directory under the current working directory by folder name
#' scminer.par <- scMINER.dir.create("MyMainDir","MyProject")
#' # Or creating a main working directory under the current working directory by relative path
#' scminer.par <- scMINER.dir.create("./MyMainDir","MyProject")
#' # Or creating a main working directory to a specific path by absolute path
#' scminer.par <- scMINER.dir.create("~/Desktop/MyMainDir","MyProject")
#' }
#' @export
scMINER.dir.create <- function(project_main_dir=NULL,project_name=NULL){
  #
  if(base::exists('scminer.par')==TRUE){
    stop('scminer.par is occupied in the current session,please manually run: rm(scminer.par) and re-try, otherwise will not change !');
  }
  #
  all_input_para <- c('project_main_dir','project_name')
  check_res <- sapply(all_input_para,function(x)check_param(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  scminer.par <- list()
  scminer.par$main.dir <- project_main_dir
  scminer.par$project.name <- project_name
  scminer.par$out.dir <- sprintf('%s/%s',scminer.par$main.dir,scminer.par$project.name)
  # create output directory
  if (!dir.exists(project_main_dir)) {
    dir.create(project_main_dir, recursive = TRUE)
  }
  if (!dir.exists(scminer.par$out.dir)) {
    dir.create(scminer.par$out.dir, recursive = TRUE)
  }
  scminer.par$out.dir.MICA <- paste0(scminer.par$out.dir, '/MICA/')
  if (!dir.exists(scminer.par$out.dir.MICA)) {
    dir.create(scminer.par$out.dir.MICA, recursive = TRUE) ## directory for MICA
  }
  scminer.par$out.dir.PLOT <- paste0(scminer.par$out.dir, '/PLOT/')
  if (!dir.exists(scminer.par$out.dir.PLOT)) {
    dir.create(scminer.par$out.dir.PLOT, recursive = TRUE) ## directory for PLOT
  }
  scminer.par$out.dir.QC <- paste0(scminer.par$out.dir, '/QC/')
  if (!dir.exists(scminer.par$out.dir.QC)) {
    dir.create(scminer.par$out.dir.QC, recursive = TRUE) ## directory for QC
  }
  scminer.par$out.dir.DATA <- paste0(scminer.par$out.dir, '/DATA/')
  scminer.par$out.dir.DATA_par <- paste0(scminer.par$out.dir, '/DATA/scminer.par.RData')
  if (!dir.exists(scminer.par$out.dir.DATA)) {
    dir.create(scminer.par$out.dir.DATA, recursive = TRUE) ## directory for DATA
  }
  scminer.par$out.dir.SJAR <- paste0(scminer.par$out.dir, '/SJAR/')
  if (!dir.exists(scminer.par$out.dir.SJAR)) {
    dir.create(scminer.par$out.dir.SJAR, recursive = TRUE) ## directory for SJARACNe
  }
  message(sprintf('Project space created, please check %s',scminer.par$out.dir))
  return(scminer.par)
}

#' CreateSparseEset
#' @description Create a S4 class which utilize 'ExpressionSet' template yet compatible with sparseMatrix type of assaydata
#' @param data Sparse expression data, could be from either of these class:c('matrix','dgTMatrix','dgCMatrix').Required
#' @param meta.data phenotype data which rownames should be the same as data colnames; Optional; Default as NULL
#' @param feature.data feature data which rownames should be the same as data rownames; Optional; Default as NULL
#' @param add.meta logical; Whether or not calculate extra pheonotype info including total number of UMI,
#' number of non-zero gene for each cell, mitochondrial percentage and spike-in gene expression percentage and store them in Biobase::pData
#'
#' @return A customized S4 class using ExpressionSet class as prototype
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' pbmc.14k.DS.eset.log2.update <- CreateSparseEset(data = exprs(pbmc.14k.DS.eset.log2),
#'                          meta.data = Biobase::pData(pbmc.14k.DS.eset.log2),
#'                          feature.data = Biobase::fData(pbmc.14k.DS.eset.log2),
#'                          add.meta = F)
#' \dontrun{
#' meta.data<-SeuratObject@meta.data
#' feature.data<-data.frame(rownames(SeuratObject@assays$RNA@data))
#' colnames(feature.data)<-"geneSymbol"
#' rownames(feature.data)<-feature.data$geneSymbol
#' eset<-CreateSparseEset(data=SeuratObject@assays$RNA@data,meta.data = meta.data,
#'                      feature.data = feature.data,add.meta = F)
#' @export
CreateSparseEset<-function(data=NULL,meta.data=NULL,feature.data=NULL,add.meta=TRUE){

  if(!class(data)[1]%in%c("dgCMatrix","dgTMatrix","matrix","dgeMatrix")){
    stop("Input format should %in% c( 'matrix','dgTMatrix','dgCMatrix','dgeMatrix,'Matrix')","\n")
    if(class(data)%in%"matrix"){
      data<-as(data,"dgCMatrix")
    }
  }

  if(!is.null(feature.data)){
    if(any(rownames(data)!=rownames(feature.data))){
      stop("Row names of feature data doesn't match with expression data row names!","\n")
    }
  } else {
    cat("Will take rownames as geneSymbol","\n")
    feature.data<-data.frame(geneSymbol=rownames(data), stringsAsFactors=F)

    if(any(duplicated(row.names(data)))){
      cat("Found duplicated rownames, convert row.names of data to an arbitrary vector!","\n")
      row.names(data)<- NULL
      row.names(feature.data)<-NULL
    }else{
      rownames(feature.data)<-rownames(data)
    }
  }


  if(!is.null(meta.data)){
    if(any(colnames(data)!=rownames(meta.data))){
      stop("Row names of meta data doesnt match with expression data Column names!","\n")
    }
  } else {
    meta.data=data.frame(row.names=colnames(data),
                         cellName=colnames(data),stringsAsFactors=F)
  }
  cat("Passed sanity check..","\n")

  if(add.meta){

    cat("Adding phenotype data based on data ..","\n")

    d<-data;rm(data)
    cells_per_gene <- Matrix::rowSums(d!=0)
    feature.data$nCells<-cells_per_gene;

    nGene <- sum(cells_per_gene > 0)
    cat("# of non-zero gene:", nGene, "\n")

    # Count the cells with >=1 identified gene(s)
    genes_per_cell <- Matrix::colSums(d!=0)
    nCell <- sum(genes_per_cell > 0)
    cat("# of non-zero cell:", nCell, "\n")

    cat("Assessing mitochondrial and spike-in genes ..","\n")

    if ("geneSymbol" %in% colnames(feature.data)){
      mito.genes <- grep(pattern = "^mt-|^MT-", x = feature.data$geneSymbol)
      spikeIn.genes <- grep(pattern = "^ERCC-|^Ercc", x = feature.data$geneSymbol)
    }else{
      cat("Mitochondrial genes and spike-in genes were not found due to lack of geneSymbol information...","\n")
      mito.genes<-integer(0)
      spikeIn.genes<-integer(0)
    }

    pd <- data.frame(row.names = colnames(d),
                     nUMI.total = Matrix::colSums(d),
                     nGene = Matrix::colSums(sign(d)),
                     percent.mito = round(Matrix::colSums(d[mito.genes, ]) / Matrix::colSums(d),8),
                     percent.spikeIn = round(Matrix::colSums(d[spikeIn.genes, ]) / Matrix::colSums(d),8),
                     stringsAsFactors = FALSE)

    meta.data<-cbind(meta.data,pd)

    cat("Creating SparseMatrix Object..","\n")
    Obj <- new( "SparseExpressionSet",
                assayData = assayDataNew( "environment", exprs=d),
                phenoData= new("AnnotatedDataFrame",data=meta.data),
                featureData= new("AnnotatedDataFrame",data=feature.data))
  }else{

    cat("Creating SparseMatrix Object..","\n")
    Obj <- new( "SparseExpressionSet",
                assayData = assayDataNew("environment", exprs=data),
                phenoData= new("AnnotatedDataFrame",data=meta.data),
                featureData= new("AnnotatedDataFrame",data=feature.data))
  }

  return(Obj)
}

#' readscRNAseqData
#' @description read scRNA-seq data, a wrapper of conventional data reading (read.delim) and 10x genomics data standarad output reading
#'
#' @param file data path to 10x genomics output folder, which normally contains 3 files (matrix.mtx, gene or feature.tsv and barcode.csv),
#'  or data path to data txt/csv/tsv file
#' @param is.10x logical, whether or not inputs are from CellRanger standard output
#' @param CreateSparseEset logical, whether or not create sparse matrix incorporated expression set
#' @param add.meta logical, whether or not calculate metadata info from expression matrixm, this is not suggested before merging/downsampling your data
#' @param ... parameters pass to read.delim if is.10x = FALSE
#' @return A list or sparse matrix expression set
#' @examples
#' demo_dir <- system.file('PBMC14KDS_DemoDataSet/DATA/10X/',package = "scMINER")
#' pbmc.14k.DS.eset <- scMINER::readscRNAseqData(file = demo_dir,
#'                                               is.10x = TRUE,
#'                                               CreateSparseEset = TRUE,
#'                                               add.meta = TRUE)
#'
#' @export
readscRNAseqData <- function(file, is.10x=TRUE, CreateSparseEset=TRUE, add.meta=F, sep=','){

  if(is.10x){

    data.path <- file
    cat("Reading 10X genomcis data in",data.path, "\n")

    f<-list.files(path=data.path,full.names = TRUE)
    f<-f[grep(".gz$",f)]

    if (length(grep(".gz",f))!=0){
      cat("Decompressing .gz files.","\n")
      for (i in f){
        i<-normalizePath(i)
        system("gunzip ",i)
      }
    }

    if (file.exists(file.path(data.path,"matrix.mtx"))){
      data.raw <- Matrix::readMM(file.path(data.path,"matrix.mtx"))
    }else{
      stop("Matrix file not found!","\n")
    }

    if(file.exists(file.path(data.path,"barcodes.tsv"))){
      barcodes <- read.delim(file.path(data.path,"barcodes.tsv"), header=FALSE,stringsAsFactors = FALSE,sep = "\t")
      colnames(barcodes)[1]<-"CellNames"
      rownames(barcodes)<-barcodes[,1]
    }else{
      stop("Barcodes file not found!","\n")
    }


    if(file.exists(file.path(data.path,"genes.tsv"))){
      genes <- read.delim(file.path(data.path,"genes.tsv"),header=FALSE,stringsAsFactors = FALSE, sep = "\t")
      colnames(genes)<-c("ensembl","geneSymbol")
      rownames(genes)<-genes[,1]
    }else if(file.exists(file.path(data.path, "features.tsv"))){
      genes <- read.delim(file.path(data.path, "features.tsv"), header=FALSE,stringsAsFactors=FALSE,sep = "\t")
      colnames(genes)<-c("ensembl","geneSymbol","biotype")[1:ncol(genes)]
      rownames(genes)<-genes[,1]
    }else{
      cat("Genes/features file not found!","\n")
    }

    dimnames(data.raw)[[1]] <- rownames(genes)
    dimnames(data.raw)[[2]] <- rownames(barcodes)

    if (CreateSparseEset){
      d <- CreateSparseEset(data=data.raw,meta.data=barcodes,
                            feature.data=genes,add.meta = add.meta)
      cat("Sparse expression set generated!","\n")
    }else{
      d = list(raw.data=data.raw,
               meta.data=barcodes,
               feature.data=genes)}
  }
  else{
    d<- read.delim(file=file, header=TRUE, row.names=1, sep=sep)
  }
  return(d)
}

#' draw.scRNAseq.QC
#' @title draw.scRNAseq.QC
#' @description generated a scRNA-seq quality control report in html with Rmarkdown
#' @param SparseEset an SparseEset generated by CreateSparseEset
#' @param project.name a character, project name to print on report
#' @param plot.dir a character, output directory for QC reports
#' @param output.cutoff logical, whether or not return a list of suggested thresholds for filtering
#' @param only.cutoff logical, whether or not only return cuttoff without plotting QC
#' @param group a character, a variable name indicate groupping information (stored in Biobase::pData) to
#' help generate violin plots
#'
#' @return an R markdown QC report and a list of suggested threshold (if specify)
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' cutoffs <- scMINER::draw.scRNAseq.QC(SparseEset = pbmc.14k.DS.eset.log2,
#'                project.name = 'test',
#'                plot.dir = '.',
#'                group = "group",
#'                output.cutoff = TRUE,only.cutoff=TRUE)
#' @export
draw.scRNAseq.QC<-function(SparseEset,project.name,
                           plot.dir="./QC/",
                           output.cutoff=TRUE,group="group",only.cutoff=FALSE){
  if(!dir.exists(plot.dir)) {dir.create(plot.dir)}

  #Calcualte Cutoffs
  pd<-Biobase::pData(SparseEset)
  cfs<-list(nCell_cutoff = max(floor(0.005 * dim(SparseEset)[2]), 1),
            umi_cf_lo = max(floor(exp(median(log(pd$nUMI.total)) - 3 * mad(log(pd$nUMI.total)))),100),
            umi_cf_hi = ceiling(exp(median(log(pd$nUMI.total)) + 3 * mad(log(pd$nUMI.total)))),
            nGene_cf = max(floor(exp(median(log(pd$nGene)) - 3 * mad(log(pd$nGene)))),50),
            ERCC_cf = round(median(pd$percent.spikeIn) + 3 * mad(pd$percent.spikeIn),4),
            mito_cf = round(median(pd$percent.mito) + 3 * mad(pd$percent.mito),4))
  if(only.cutoff==TRUE) return(cfs)
  rmarkdown::render(input=system.file("rmd", "SparseEset_QC_report.Rmd", package = "scMINER"),
                    output_dir = plot.dir,
                    output_file = paste0(project.name,"_scRNAseq_QC.html"),
                    clean=TRUE,
                    quiet =FALSE,
                    params=list(
                      Obj=SparseEset,
                      projectName=project.name,
                      cfs=cfs,
                      output.cutoff=output.cutoff,
                      group=group))

  return(cfs)
}


#' Draw barplot for composition study
#' @param input_eset ExpressionSet that include group information in phenotype data
#' @param group_by Group criteria for bars, should be a variable stored in Biobase::pData(input_eset)
#' @param color_by Coloring criteria of bar fractions, should be a variable stored in Biobase::pData(input_eset)
#' @param colors color values to feed in scale_fill_manual, default as NULL; If NULL, then default color for ggplot will be used
#'
#' @return a ggplot object
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' draw.group.barplot(pbmc.14k.DS.eset.log2,
#'                    group_by = 'celltype',
#'                    color_by = 'celltype')
#' @export
#'
draw.group.barplot<-function(input_eset,
                             group_by,
                             color_by,
                             colors = NULL){

  input<-Biobase::pData(input_eset)
  if(group_by%in%colnames(input)){
    if(color_by%in%colnames(input)){
      input[,group_by]<-as.factor(input[,group_by])
      input[,color_by]<-as.factor(input[,color_by])
    }else{
      stop ('color_by name not found in Biobase::pData!','\n')
    }
  }else{
    stop ('group_by name not found in Biobase::pData!','\n')
  }


  p <- ggplot(data=input, aes_string(group_by))+
    geom_bar(aes_string(fill=color_by), position="fill")+
    xlab(label="Condition")+
    ylab(label="Fraction")+
    theme_classic()+
    guides(fill=guide_legend(size=15,title=""))+
    theme(plot.title = element_text(size=25),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14,face="bold"))

  if (!is.null(colors)){
    p <- p + scale_fill_manual(values = colors)
  }

  return(p)
}


#' Visualize marker score of different cell types on bubbleplot
#'
#' @title Generate visualization for marker scores via bubble plot
#' @description  Marker visualizatoin from known markers/signatures, requires knowledge-based marker list as input
#' @param ref reference dataframe, includes positive or negative markers for different cell types;
#' Specify first column as different cell types, second columns as markers, third columns as weight (postive or negative marker)
#' @param input_eset expressionSet/SparseExpressionSet object with clustering membership stored in Biobase::pData
#' @param group_name a character, the variable containing clustering label in Biobase::pData(eset); or any other group information stored in Biobase::pData(eset)
#' @param save_plot logical, whether or not save your plot; if TRUE, plot will be saved as plot_name
#' @param width default as 8, inch as unit
#' @param height default as 5, inch as unit
#' @param plot_name plot name, please include plot type
#' @param feature feature type from second column of your reference , should be in colnames(Biobase::fData(eset))
#' @return A ggplot object
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' markers_file <- system.file('PBMC14KDS_DemoDataSet/DATA/',
#'  'Immune_signatures.xlsx',
#'   package = "scMINER")
#' markers <- openxlsx::read.xlsx(markers_file)
#' draw.marker.bbp(ref = markers, input_eset = pbmc.14k.DS.eset.log2,
#'                width = 6, height = 4, feature = "geneSymbol",
#'                group_name = "ClusterRes", save_plot = FALSE)
#' \dontrun{
#' }
#'
#' @export
draw.marker.bbp<-function(ref = NULL,input_eset,
                          feature='geneSymbol',group_name="ClusterRes",
                          save_plot = FALSE,
                          width=8, height=5,
                          plot_name="AnnotationBubbleplot.png"){

  #exp<-apply(exprs(eset),2,std)
  #filter reference marker sets

  if (!feature%in%colnames(Biobase::fData(input_eset))) stop('Please check your feature!')
  colnames(ref)<-c("celltype","markers","weight")
  ref<-dplyr::filter(ref,markers%in%Biobase::fData(input_eset)[,feature])
  indx<-which(Biobase::fData(input_eset)[,feature]%in%ref$markers)
  if(length(indx)==0) stop("No genes from the reference list could be found in data!","\n")

  exp<-as.matrix(exprs(input_eset))[indx,]
  rownames(exp)<-Biobase::fData(input_eset)[,feature][indx]

  celltypes<-unique(ref$celltype)

  ac<-matrix(NA,nrow=ncol(exp),ncol=length(celltypes),dimnames = list(colnames(exp),celltypes))
  for(i in 1:length(celltypes)){
    cat(i,"\n")
    ref.sel<-dplyr::filter(ref,celltype==celltypes[i])
    n <- length(unique(ref.sel$markers))

    if(n>1){
      mat<-t(exp[ref.sel$markers,])%*%as.vector(ref.sel$weight)
      ac[,i]<-mat[,1]/n
    }else if (n==1){
      ac[,i]<-exp[ref.sel$markers,]
    }
  }

  ac_norm<-apply(ac,2,scale) #column normalization

  n_mtx<-(ac>0.5)
  df_n<-data.frame(label=Biobase::pData(input_eset)[,group_name],n_mtx)
  tmp1 <- aggregate(df_n,by=list(df_n$label),mean);tmp1$label<-tmp1$Group.1;df_n <- tmp1[,-1]
  #df_n<-aggregate(.~label,df_n,mean)
  library(reshape2)
  df_n_melt<-melt(df_n,id.vars = "label")
  df<-data.frame(label=Biobase::pData(input_eset)[,group_name],ac_norm);
  df<-df[,colSums(is.na(df))<nrow(df)];#remove NA columns
  #df<-aggregate(.~label,df,mean)
  tmp1 <- aggregate(df,by=list(df$label),mean);tmp1$label<-tmp1$Group.1;df <- tmp1[,-1]
  #
  input<-t(apply(df[,-1],1,scale))#row normalization
  input<-as.data.frame(cbind(df[,1],input))
  rownames(input)<-rownames(df)
  colnames(input)<-colnames(df)
  df_melt<-melt(input,id.vars = "label")
  ##
  df_n_melt<-melt(df_n[,colnames(input)],id.vars = "label") # 20230724

  if(all(df_melt[,c(1,2)]==df_n_melt[,c(1,2)])){

    d<-cbind(df_melt,df_n_melt[,3])
    colnames(d)<-c("Cluster", "CellType", "MarkerScore","ExpressionPercentage")
    d$Cluster<-as.factor(d$Cluster)
    d$MarkerScore <- as.numeric(d$MarkerScore)
    p<- draw.bubblePlot2(df=d, xlab="Cluster",ylab="CellType",
                        clab="MarkerScore",slab="ExpressionPercentage",
                        plot.title="Cell type annotation for each cluster")
  }

  if(save_plot){ggsave(plot = p, filename = plot_name , units="in",
                       width = width,height = height,dpi = 300)}
  print(d)
  return(p)
}

#' @title feature_vlnplot
#' @description This plot will visualize feature info in violin plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param stat a character, whether to plot median or mean as a black dot on violinplot
#' @param group_by character, which group info to visualize as x axis
#' @param color_by character, which group info to define color, if NULL, then violin plots will be colored by 'group_by'
#' @param colors character vector, default as NULL, will use ggplot default color palette
#' @param ylabel a character, title of y axis
#' @param boxplot logical, whether to plot boxplot on violinplot
#' @param title.size numerical, default as 5
#' @param ncol cordinates for y axis
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' genes_of_interest <-c("CD3D", "CD27", "IL7R",
#'                       "SELL", "CCR7", "IL32",
#'                       "GZMA", "GZMK", "DUSP2",
#'                       "CD8A", "GZMH", "GZMB",
#'                       "CD79A", "CD79B", "CD86", "CD14")
#' scMINER::feature_vlnplot(input_eset = pbmc.14k.DS.eset.log2,
#'                          target = genes_of_interest,
#'                          feature = "geneSymbol",
#'                          group_by = "ClusterRes",
#'                          ylabel = "log2Exp", ncol = 4)
#' @export
feature_vlnplot <- function(input_eset,
                            target=NULL,feature="geneSymbol",
                            group_by="celltype",ylabel="Expression",
                            color_by=NULL,colors=NULL,
                            ncol=3,stat="median",
                            boxplot=FALSE,title.size=5){

  if(!group_by%in% colnames(Biobase::pData(input_eset))) stop('Please check your group_by information!','\n')
  if(!feature%in% colnames(Biobase::fData(input_eset))) stop('Please check your feature information!','\n')

  # extract input information
  input <- exprs(input_eset)
  indx<-which(Biobase::fData(input_eset)[,feature]%in%target)
  gn<-Biobase::fData(input_eset)[,feature][indx]

  if(length(indx)==0) stop('No target feature found in data!','\n')

  label <- as.factor(Biobase::pData(input_eset)[,group_by])
  if (is.null(color_by)) color_by=group_by
  condition<- as.factor(Biobase::pData(input_eset)[,color_by])

  # Gene expression visualized as columns
  if (length(target)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    colnames(target_values)<-gn
    df <- data.frame(target_values,cluster=label,condition=condition)
  }else {
    target_values<-input[indx,]
    df <- data.frame(target_values,cluster=label,condition=condition)
    colnames(df)[1]<-gn
  }

  df_melt <- reshape2::melt(df, id.vars=c("cluster","condition"))

  p <- ggplot(df_melt, aes(x=cluster, y=value, fill=condition))+
    theme_classic()+
    geom_violin(trim=TRUE,scale="width",na.rm = TRUE,size=0.4)

  if(!is.null(stat)){
    if (stat=="median") p <- p + stat_summary(fun.y=median, geom="point", size=1.2, color="black",position=position_dodge(width=1))
    else if (stat=="mean") p <- p + stat_summary(fun.y=mean, geom="point", size=1.2, color="black",position=position_dodge(width=1))
    else cat("Stat not supported, please check your spelling.","\n")}

  if(boxplot) p <- p + geom_boxplot(fill="white",width=0.1,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)

  p <- p + facet_wrap(~variable,scales = "free",ncol = ncol) +
    labs(x=group_by,y=ylabel)+
    theme(axis.text.x = element_text(size=10),
          plot.title = element_text(size = title.size, face = "bold"),
          strip.background = element_rect(fill="#FFFFFF"))

  if (!is.null(colors)) p <- p+ scale_fill_manual(values=colors)

  if (ylabel=="Activity") { p <- p + geom_boxplot(width=0.2,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)}

  return(p)
}

#' @title Visualize gene expression level on scRNA-seq data via heatmap
#' @description This plot will visualiz feature info in scatter plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature a character, which feature to visualize
#' @param target a character or a character vector indicating feature names
#' @param group_name a character, label to visualize on the top of heatmap
#' @param name character, name of value visualized in color scale
#' @param cluster_rows logical, if or not cluster rows
#' @param colors color palette
#' @param plot_name character, name of heamap
#' @param save_plot logical, whether to save plots or not
#' @param width numerical
#' @param height numerical
#' @param ... parameter to be passed to ComplexHeatmap::Heatmap
#' @return a ggplot object
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' genes_of_interest <-c("CD3D", "CD27", "IL7R",
#'                       "SELL", "CCR7", "IL32",
#'                       "GZMA", "GZMK", "DUSP2",
#'                       "CD8A", "GZMH", "GZMB",
#'                       "CD79A", "CD79B", "CD86", "CD14")
#' scMINER::feature_heatmap(input_eset = pbmc.14k.DS.eset.log2,
#'                      target = genes_of_interest,
#'                      group_name = "ClusterRes",
#'                      save_plot = FALSE,
#'                      width = 6, height = 6,
#'                      name = "log2Exp")
#' @export
feature_heatmap <- function(input_eset,target,feature="geneSymbol",
                            group_name="label",name="log2Exp",
                            save_plot=TRUE,width=4,height=8,
                            cluster_rows = FALSE,
                            colors = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)),
                            plot_name="GeneHeatmap.png",
                            ...){

  input <- exprs(input_eset)
  gn<-intersect(target,Biobase::fData(input_eset)[,feature])
  indx<-match(gn,Biobase::fData(input_eset)[,feature])

  exp<-exprs(input_eset)[indx,]
  rownames(exp)<-gn
  lab<-Biobase::pData(input_eset)[,group_name];names(lab) <- sampleNames(input_eset)

  #re-order expressionmatrix and label
  ranks<-names(sort(lab,decreasing = FALSE))
  exp.ordered<-as.matrix(exp[,ranks])
  lab.ordered<-lab[ranks]
  df<-data.frame(scMINER=lab.ordered)

  #Define color annotations
  n<-length(unique(lab.ordered))
  ncols <- scales::hue_pal()(n)
  names(ncols) <- unique(lab.ordered)
  myanndf = ComplexHeatmap::HeatmapAnnotation(df = df,col=list(scMINER = ncols))
  mycolors = colors

  hmp <- ComplexHeatmap::Heatmap(exp.ordered, col = mycolors, name = name,
                 show_row_names = TRUE,
                 show_column_names = FALSE,
                 cluster_rows = cluster_rows,
                 cluster_columns = FALSE,
                 top_annotation = myanndf,
                 ...)

  if(save_plot){
    png(filename = plot_name, width=width,height=height,units ="in",res = 300)
    ComplexHeatmap::draw(hmp)
    dev.off()
  }

  return(hmp)
}

#' @title feature_highlighting
#' @description This plot will visualize feature info on scatter plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param ylabel a characterm, title of y axis
#' @param title.size numerical, default as 5
#' @param x cordinates for x axis
#' @param y cordinates for y axis
#' @param wrap_by character, variable to wrap plot with
#' @param ncol cordinates for y axis
#' @param alpha numerical, default as 0.8
#' @param colors color palette for feature highlighting
#' @param pct.size numrical, point size
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' genes_of_interest <-c("CD3D", "CD27", "IL7R",
#'                       "SELL", "CCR7", "IL32",
#'                       "GZMA", "GZMK", "DUSP2",
#'                       "CD8A", "GZMH", "GZMB",
#'                       "CD79A", "CD79B", "CD86", "CD14")
#' scMINER::feature_highlighting(input_eset = pbmc.14k.DS.eset.log2,
#'                       target = genes_of_interest,
#'                       feature = "geneSymbol",
#'                       ylabel = "log2Exp",
#'                       x = "X", y = "Y",
#'                       pct.size = 0.5)
#' @export
feature_highlighting<-function(input_eset,target=NULL,
                               feature="geneSymbol",
                               x="X",y="Y",
                               wrap_by=NULL,
                               ylabel="Expression",pct.size=0.8,
                               title.size=15,ncol=4, alpha=0.8,
                               colors=colorRampPalette(c("#E3E3E3", "#BCA2FC","#4900FE"),interpolate="linear")(8)){

  # change it to expr is ok
  input<-as.matrix(exprs(input_eset))
  indx<-which(Biobase::fData(input_eset)[,feature]%in%target)
  if(length(indx)==0) stop("Target feature not found.")

  gn<-Biobase::fData(input_eset)[,feature][indx]
  id.vars<-c(x,y,wrap_by)
  projection<-Biobase::pData(input_eset)[colnames(input),id.vars]

  #gene expression visualized as columns
  if (length(indx)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    colnames(target_values)<-gn

    proj_target <- cbind(projection,target_values)
    proj_target_melt <- reshape2::melt(proj_target, id.vars=id.vars)

    p<- ggplot(proj_target_melt, aes_string(x, y)) +
      theme_classic()+
      facet_wrap(c("variable",wrap_by),scales = "free",ncol = ncol)
    labs(title="")

  }else{
    target_values <- input[indx,]
    proj_target <- cbind(projection,target=target_values)
    proj_target_melt <- reshape2::melt(proj_target, id.vars=id.vars)

    p<- ggplot(proj_target_melt, aes_string(x, y)) +
      theme_classic()+
      labs(title=target,scales = "free")

    if(!is.null(wrap_by)) p <- p + facet_wrap(c(wrap_by),scales = "free",ncol = ncol)
  }#indx = 1

  p<- p + geom_point(aes(colour=value),size=pct.size,alpha=alpha) +
    scale_colour_gradientn(colors=colors)   +
    theme(plot.title = element_text(size = title.size, face = "bold"),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10))+
    labs(color=ylabel)

  return(p)
}


#' preMICA.filtering
#' @title preMICA.filtering
#' @description scRNA-seq filtering function
#' @param SparseEset the sparseEset object outputted from draw.scRNAseq.QC
#' @param cutoffs a list outputted by draw.scRNAseq.QC, if NULL, manual input will be required
#' @param gene_filter logical; or a numerical number, indicating lower threshold for gene filtering based on how many non-zero cells each gene expressed in
#' @param ERCC_filter logical; a numerical number, indicating upper threshold put on ERCC percentage for cell filtering
#' @param Mito_filter logical;a numerical number, indicating upper threshold put on Mitochondrial gene expression fraction for cell filtering
#' @param nGene_filter logical; a numerical number, indicating lower threshold put on number of gene expression in each cell for cell filtering
#' @param nUMI_filter logical;a vector of two numerical number, indicating lower threshold and upper threshold put on number of total UMI for cell filtering
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' cutoffs <- draw.scRNAseq.QC(SparseEset = pbmc.14k.DS.eset.log2,
#'                                      project.name = project_name,
#'                                      plot.dir = 'test/',
#'                                      group = "group",
#'                                      output.cutoff = TRUE,
#'                                      only.cutoff = TRUE)
#' pbmc.14k.DS.eset.filter <- preMICA.filtering(SparseEset = pbmc.14k.DS.eset.log2,
#'       cutoffs = cutoffs)
#' @return A Sparse expression set
#' @export
preMICA.filtering <- function(SparseEset,
                              cutoffs,
                              gene_filter=TRUE,
                              nGene_filter=TRUE,
                              nUMI_filter=TRUE,
                              ERCC_filter=TRUE,
                              Mito_filter=TRUE)
{
  cat("Pre-filtering dimension: ",dim(SparseEset),"\n")

  if(is.null(cutoffs)){
    if(all(is.logical(gene_filter),is.logical(ERCC_filter),
           is.logical(Mito_filter),is.logical(nGene_filter),is.logical(nUMI_filter))){
      stop("No filtering will be conducted due to lack of numerical threshold input.","\n",
           "You can take default cutoff calculated by function draw.scRNAseq.QC and feed them into cutoffs,
             or input them manually in each function parameters.")
    }else{
      if(any(isTRUE(gene_filter),isTRUE(ERCC_filter),isTRUE(Mito_filter),isTRUE(nGene_filter),isTRUE(nUMI_filter)))
        stop("When cutoffs=NULL, please indicate numerical threshold instead of 'TRUE'
                 if you want to do filtering on that particular criteria !","\n")
      cfs<-list()
      if (gene_filter) cfs$nCell_cutoff=gene_filter
      if (nGene_filter) cfs$nGene_cf=nGene_filter
      if (!isFALSE(nUMI_filter)) cfs$umi_cf_lo=nUMI_filter[1];cfs$umi_cf_hi=nUMI_filter[2]
      if (ERCC_filter) cfs$ERCC_cf=ERCC_filter
      if (Mito_filter) cfs$mito_cf=Mito_filter
    }
  }else{cfs<-cutoffs}

  if(isFALSE(nUMI_filter)){cfs$umi_cf_lo=0; cfs$umi_cf_hi=Inf}
  cell <- which((SparseEset$nUMI.total > cfs$umi_cf_lo) & (SparseEset$nUMI.total < cfs$umi_cf_hi))
  if(nGene_filter) cell <- intersect(cell, which(SparseEset$nGene > cfs$nGene_cf))
  if(ERCC_filter&(cfs$ERCC_cf!=0)) cell <- intersect(cell, which(SparseEset$percent.spikeIn < cfs$ERCC_cf))
  if(Mito_filter&(cfs$mito_cf!=0)) cell <- intersect(cell, which(SparseEset$percent.mito < cfs$mito_cf))

  cat("Below threshold were used:","\n")
  print(cfs)

  #gene filtering
  if(gene_filter){
    gene <- unname((which(Biobase::fData(SparseEset)$nCells >= cfs$nCell_cutoff)))
    eset.sel<-SparseEset[gene,]

    cat("Gene expressed in less than ", cfs$nCell_cutoff,
        "cells (",(dim(SparseEset)[1]-length(gene))*100/dim(SparseEset)[1],"% genes) were filtered","\n",
        "-Filtered expression matrix dimension:",dim(eset.sel),"\n")
  }else { eset.sel <- SparseEset }

  eset.sel2 <- eset.sel[,cell]
  cat("A total of ",dim(SparseEset)[2]-length(cell),
      "(",(dim(SparseEset)[2]-length(cell))*100/dim(SparseEset)[2],"%) cells were filtered","\n",
      "-Filtered expression matrix dimension:",dim(eset.sel2),"\n")

  cat("Data filtering done!","\n")
  cat("====================","\n")

  return(eset.sel2)
}#end preMICA.filtering

#' Generate MICA input accepted txt or h5ad file
#' @description A utility function that helps generate MICA input from a data matrix with rownames and colnames
#' @usage generateMICAinput(d, filename="project_name_MICAinput.h5")
#' @param d matrix with colnames as cell/sample info, rownames as gene/feature info
#' @param filename filename of your MICA input file, supported format: txt or h5
#' @param sampleN integer, number of cells that should be sampled for computational efficiency purpose, default is 50000. If sampleN=NULL, no sampling will be performed.
#' @param seed integer, the random seed for sampling, default is 1.
#' @param scminer.par list for the parameter settings in scMINER pipeline, optional.
#' @return A txt file or a h5 file that could be read in MICA
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' MICA.cmd <- generateMICAinput(eset = pbmc.14k.DS.eset.log2 ,
#'                       filepath = 'test.h5ad')
#' @export
generateMICAinput <- function(eset, filepath, sampleN=50000, seed=1, scminer.par=NULL){
  recommend_cmd <- NULL
  exp_mat <- as.matrix(exprs(eset))
  obs_dat <- Biobase::pData(eset)
  var_dat <- Biobase::fData(eset)
  cell_size <- ncol(exp_mat)
  if(is.null(sampleN)==FALSE){
    if(cell_size>sampleN){
      set.seed(seed=seed)
      s1 <- sample(1:cell_size,size=sampleN,replace = FALSE)
      exp_mat <- exp_mat[,s1,drop=F]
      obs_dat <- obs_dat[s1,,drop=F]
      message(sprintf('original %s cells are downsampled to %s cells for MICA input',cell_size,length(s1)))
      }
  }
  if (length(grep(".txt$", filepath))!=0){
    mica.input <- as.data.frame(t(exp_mat))
    cat("Writing MICA input to a .txt file...", "\n")
    write.table(mica.input, file = filepath, sep = "\t",
                row.names = TRUE, col.names = NA, quote = FALSE)
  }else if(length(grep(".h5ad$", filepath))!=0){
    cat("Writing MICA input to a .h5ad file...", "\n")
    ad <- anndata::AnnData(
                          X = t(exp_mat),
                          obs = obs_dat,
                          var = var_dat)
    anndata::write_h5ad(ad, filepath)
  }else{
    stop("Your filepath should be ended with .txt or .h5ad", "\n")
  }
  cell_size <- ncol(exp_mat)
  if(cell_size>=5000){
    cat('For dataset with more than 5k cells, MICA GE mode is recommended.\n')
    if(is.null(scminer.par)==FALSE){
      recommend_cmd <- sprintf('mica ge -i %s -o %s -pn %s -ar 3.0 -ss 0.2 -nw 1 -nnm 80',scminer.par$out.dir.MICA_input,scminer.par$out.dir.MICA,scminer.par$project.name)
      cat(sprintf("Suggested command line is:\n\n%s\n
                  Where options represent:
                  -ar determines the maximum size of the communities (default: 3.0);
                  -ss is the step size to sweep resolutions in a range (default: 0.2);
                  -nw specifies the number of workers to run in parallel (default: 1);
                  -nnm is the number of neighbors to build mutual information-based nearest neighbor graph (default 80).
                  Use mica ge -h to find more options of MICA GE mode.",recommend_cmd))
    }
  }else{
    cat('For dataset with less than 5k cells, MICA MDS mode is recommended.\n')
    if(is.null(scminer.par)==FALSE){
      recommend_cmd <- sprintf('mica mds -i %s -o %s -pn %s -nc 4 5 6 7 8 9 10 -dk 19',scminer.par$out.dir.MICA_input,scminer.par$out.dir.MICA,scminer.par$project.name)
      cat(sprintf("Suggested command line is:\n\n%s\n
                -pn specifies a project name for naming the output files;
                -nc is an array of integers delimited by a single space, where each integer specifies a k to perform a k-mean clustering;
                -dk can be an integer or an array of integers delimited by a single space (default is 19), it specifies the number of dimensions used in k-mean clusterings.
                Use mica mds -h to see more options of MICA MDS mode.",recommend_cmd))
    }
  }
  cat("Done.","\n")
  return(recommend_cmd)
}


#' @title readMICAoutput
#' @description Read MICA input and output to create an expressionSet for downstream analysis
#'
#' @param eset a SparseMatrix Eset
#' @param input_file input expression txt file of MICA pipeline
#' @param output_file output ClusterMem.txt file from MICA pipeline
#' @param load_ClusterRes logical, if TRUE, clustering results will be store at Biobase::pData(eset)$label
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' MICA_output <- system.file('PBMC14KDS_DemoDataSet/MICA/clustering_umap_euclidean_19.txt',package = "scMINER")
#' pbmc.14k.DS.eset.log2 <- scMINER::readMICAoutput(eset = pbmc.14k.DS.eset.log2,
#'                           load_ClusterRes = TRUE,
#'                           output_file = MICA_output)
#' @return A sparse expressionSet object
#' @export
readMICAoutput<-function(eset=NULL, input_file, output_file,load_ClusterRes =TRUE){

  res <- read.table( output_file, # MICA output text file
                     header = TRUE,
                     stringsAsFactors = FALSE)

  if (!is.null(eset)){
    if(!all(colnames(eset)==res[,1])) stop("Check your eset, Output file doesnt match with input list object.","\n")
    else{
      eset$X=res[,2]
      eset$Y=res[,3]
    }
  }else{
    cat("Reading input...","\n")
    d <- read.delim(input_file,header = FALSE,
                    stringsAsFactors = FALSE,
                    quote = "",
                    check.names = FALSE)

    gn<-unname(unlist(d[1,-1]))
    input<-t(d[-1,-1]);colnames(input)<-d$V1[-1]
    class(input)<-"numeric"
    rownames(input)<-gn;

    pd<-data.frame(row.names=colnames(input),cellNames=colnames(input),
                   X=res[,2],Y=res[,3],stringsAsFactors=FALSE)

    eset<-CreateSparseEset(data=input,meta.data=pd,add.meta = F)
    cat("SparseMatrix Expression Set Generated!","\n")
  }

  if(load_ClusterRes){eset$ClusterRes <- as.factor(res$label);
  cat("Clustering info is under 'ClusterRes' slot.","\n")}
  return(eset)
}


#' @title plot MICA clustering results or other meta variables
#' @description This function helps to generate a ggplot object for phenotypic visualization
#' @param input_eset ExpressionSet that include visualization coordinates in phenotype data
#' @param color_by Coloring criteria of data points, should be a variable stored in Biobase::pData(input_eset)
#' @param colors character, color values, if NULL then use ggplot default color
#' @param X character, column name of x axis
#' @param Y character, column name of y axis
#' @param show_label logical, whether or not to show label on tSNE plot
#' @param label.size numerical, size of label ploted on figure
#' @param title.size numerical, size of plot title, default as 10
#' @param title.name character, title of plot, default as NULL
#' @param pct numerical, size of point, default as 0.5
#' @param aplha numerical, indicate point transparency
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' scMINER::MICAplot(input_eset = pbmc.14k.DS.eset.log2, X = "X", Y = "Y",
#'                   color_by = "ClusterRes", pct = 0.5)
#' \dontrun{
#' }
#' @export
MICAplot<-function(input_eset,
                   color_by= "ClusterRes", colors= NULL,
                   X=NULL,Y=NULL,show_label=FALSE,label.size=10,
                   title.size=20,title.name="",pct=0.5,alpha=1){

  input<-Biobase::pData(input_eset)
  if(!color_by%in%colnames(input)){stop("Label name not found in phenotype data!","\n")}
  if(!X%in%colnames(input)|!Y%in%colnames(input)){stop("Please check your x or y axis assignment!","\n")}
  p <- ggplot(data=input,aes_string(x=X, y=Y,color = color_by))+
    geom_point(size=pct, alpha=alpha)+
    labs(title=title.name, x = X, y= Y)+
    theme_classic() +
    theme(plot.title = element_text(size=title.size),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14,face="bold"))

  if (show_label){
    loc <- stats::aggregate(input[,c(X,Y)],by=list(Cluster=input[,color_by]),mean)
    p <- p + geom_text(data = loc, aes_string(x=X,y=Y,label="Cluster"),
                       check_overlap = TRUE, na.rm = FALSE,size=label.size,
                       show.legend = F, inherit.aes = FALSE)
  }# add number on tsne plot if TRUE

  # add number of cells in each catergory on legend
  if (class(input$label)!="numerical"){
    input$label<-input[,color_by]
    lgnd<-paste0(names(table(input$label))," (",table(input$label),")")

    if (!is.null(colors)){
      p <- p + scale_color_manual(values=colors, labels=lgnd)
    }else{
      p <- p + scale_color_discrete(labels=lgnd)
    }
    p <- p + guides(col=guide_legend(override.aes = list(size=10),
                                     title=paste0(color_by,"\n(",dim(input)[1] ,")")),
                    ncol=ceiling(length(unique(input$label))/10))
  }else{
    if (!is.null(colors)) p <- p + scale_color_gradientn(values=colors)
    p<- p + guides(colour = guide_legend(override.aes = list(size = 10)))
  }
  return(p)
}

#' generateSJARACNeInput
#'
#' @title Generate SJARACNE input with designed folder structure
#' @description This function helps to generate appropriate input files for SJARACNe pipeline.
#' It can take transcription factor/signaling gene reference from internal(stored in package) or external (manual define)
#'
#' @param input_eset An expressionSet
#' @param ref c("hg", "mm"), could be a manually defined geneSymbol vector
#' @param funcType c("TF","SIG", NULL), if NULL then both TF and SIG will be considered
#' @param input_driver, a list of drivers for calculation. If NULL, curated driver list in the R package will be used.
#' @param sampleN integer, number of cells that should be sampled per group for computational efficiency purpose, default is 1000. If sampleN=NULL, no sampling will be performed.
#' @param seed integer, the random seed for sampling, default is 1.
#' @param wd.src output path
#' @param group_name name of group for sample identification
#' @param symbolColumnName name for the column that save gene symbols.
#' @return SJARACNe input files for each subgroups
#' @keywords SJARACNe
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' SJAR.cmd.tf <- generateSJARACNeInput(
#'                input_eset = pbmc.14k.DS.eset.log2,
#'                funcType = "TF",
#'                ref = "hg",  # human
#'                wd.src = 'test/',  # output directory
#'                group_name = "celltype")
#' \dontrun{
#' }
#' @export
generateSJARACNeInput<-function(input_eset,
                                ref=NULL,funcType=NULL,
                                input_driver=NULL,
                                sampleN=1000,
                                seed=1,
                                wd.src,group_name,symbolColumnName='geneSymbol'){
  SJAR.cmd <- NULL
  if (!dir.exists(wd.src)) dir.create(wd.src,recursive = TRUE)
  # input_driver, 20230919
  if(is.null(input_driver)==TRUE | length(input_driver)<1){
    if (ref%in%c("hg","mm")){
      ref_file<-system.file("RData",paste0("tf_sigs_",ref,".RData"),package = "scMINER")
      load(ref_file)
      cat("Using references from: ", ref_file,"\n")
      sig.ref <- NULL;tf.ref <- NULL
      tf.ref<- dplyr::filter(tf_sigs, isTF==TRUE)$geneSymbol
      sig.ref<- dplyr::filter(tf_sigs, isSIG==TRUE)$geneSymbol
      if (!is.null(funcType)){
        if (funcType=="TF") sig.ref <- NULL
        else if (funcType=="SIG") tf.ref <- NULL
      }
    }else {
      if (funcType=="TF") tf.ref <- ref
      else if (funcType=="SIG") sig.ref <- ref
      else warning("Activity calculations will not be supported!","\n")
    }
  }else{
    sig.ref <- NULL;tf.ref <- NULL
    ref <- input_driver
    if (funcType=="TF") tf.ref <- ref
    else if (funcType=="SIG") sig.ref <- ref
  }
  if(group_name%in%colnames(Biobase::pData(input_eset))){
    groups <- unique(Biobase::pData(input_eset)[,group_name])
    all.SJAR.cmd.1 <- c()
    all.SJAR.cmd.2 <- c()
    for (i in 1:length(groups)){
      grp.tag<-groups[i];
      eset.sel <- input_eset[,which(Biobase::pData(input_eset)[,group_name]==grp.tag)]
      ########### 20230913, downsample eset.sel
      if(is.null(sampleN)==FALSE){
        cell_size <- ncol(Biobase::exprs(eset.sel))
        if(cell_size>sampleN){
          set.seed(seed=seed)
          s1 <- sample(1:cell_size,size=sampleN,replace = FALSE)
          eset.sel <- eset.sel[,s1]
          message(sprintf('original %s cells are downsampled to %s cells for SJARACNe input',cell_size,length(s1)))
        }
      }
      res <- SJARACNe_filter(eset.sel=eset.sel,tf.ref=tf.ref,
                             sig.ref=sig.ref,wd.src=wd.src,grp.tag=grp.tag,
                             symbolColumnName=symbolColumnName)
      print(res)
      ###########
      #return(c(dir.cur,grp.tag,f.exp,f.tfsig))
      ## suggested command for SJARACNe
      SJAR.cmd.1 <- sprintf('sjaracne lsf -e %s -g %s -o %s_final -n 100 -pc 1e-5 -tmp %s/tmp/ -j ../SJARACNe/SJARACNe/config/config_cwlexec.json',
                            res[3],res[4],res[1],res[1])
      SJAR.cmd.2 <- sprintf('sjaracne local -e %s -g %s -o %s_final -n 100 -pc 1e-5 -tmp %s/tmp/',
                            res[3],res[4],res[1],res[1])
      all.SJAR.cmd.1 <- c(all.SJAR.cmd.1,SJAR.cmd.1)
      all.SJAR.cmd.2 <- c(all.SJAR.cmd.2,SJAR.cmd.2)
    }#end for
    SJAR.note <- sprintf('-e: input expression file;
                            -g: input candidate driver file;
                            -o: output directory;
                            -n: the number of bootstraps;
                           -pc: the consensus p-value threshold;
                            -j: json config file.')
    cat(sprintf("Suggested command line for lsf is:\n\n%s\n
                  Suggested command line for local is:\n\n%s\n
                  Where options represent:%s",
                paste(all.SJAR.cmd.1,collapse='\n'),
                paste(all.SJAR.cmd.2,collapse='\n'),SJAR.note))
  }else{
    stop("Lack group info, please check your group_name.","\n")
  }#end if
  save(input_eset, file=file.path(wd.src,"Input.eset")) # save input file as expressionSet
  return(list(lsf=all.SJAR.cmd.1,local=all.SJAR.cmd.2))
}#end function

### inner functions
#' SJARACNe_filter
#' @description This is the inner function to help generate SJARACNe input for scRNA-seq data,
#' all non-informative (zero genes) will be filtered in by this function
#' @param eset.sel ExpressionSet to generate SJaracne input
#' @param tf.ref A vector of reference transcription factors
#' @param sig.ref A vector of reference signaling genes
#' @param wd.src path to store SjAracne input
#' @param grp.tag name of group for identification
#' @param symbolColumnName name for the column that save gene symbols.
#' @details Non-expressed genes in subgroups are filtered.
#' tf.ref should be coordinate with featureNames(eset.sel).
#' @return A folder with picked master regulator and filtered gene expression matrix
SJARACNe_filter<-function(eset.sel,tf.ref,sig.ref,wd.src,grp.tag,symbolColumnName='geneSymbol'){
  cat(grp.tag,'\n')

  #Biobase::fData(eset.sel)$IQR<-apply(exprs(eset.sel),1,IQR)
  # exclude genes with all zero
  eset.sel<-eset.sel[apply(exprs(eset.sel),1,function(xx){sum(xx)!=0}),]
  ## 20231011
  if(!'geneNames' %in% colnames(Biobase::fData(eset.sel))){
    if(symbolColumnName %in% colnames(Biobase::fData(eset.sel))){
      Biobase::fData(eset.sel)$geneNames<-Biobase::fData(eset.sel)[,symbolColumnName]
    }
  }

  ni<-nrow(eset.sel);ni
  ns<-ncol(eset.sel);ns
  ng<-nlevels(factor(Biobase::fData(eset.sel)$geneNames));ng
  tag<-paste(grp.tag,nrow(eset.sel),ng,ncol(eset.sel),sep='_');tag

  dir.cur<-file.path(wd.src,tag);dir.cur
  dir.create(dir.cur,recursive = TRUE)

  #write exp data to exp format
  expdata<-data.frame(cbind(isoformId=featureNames(eset.sel),geneSymbol=Biobase::fData(eset.sel)$geneSymbol,as.matrix(exprs(eset.sel))),stringsAsFactors = FALSE)
  f.exp<-file.path(dir.cur,paste(grp.tag,"_",ni,"_",ng,"_",ns,".exp",sep=''));f.exp
  write.table(expdata,file=f.exp,sep="\t",row.names=FALSE,quote=FALSE)

  if (!is.null(tf.ref)){
    dir.create(file.path(dir.cur,'tf'),recursive = TRUE)
    tf.eset.sel<-subset(eset.sel,Biobase::fData(eset.sel)$geneSymbol%in%tf.ref)
    dim(tf.eset.sel)
    f.tf<-file.path(dir.cur,'tf',paste(grp.tag,"_",nrow(tf.eset.sel),"_",nlevels(factor(Biobase::fData(tf.eset.sel)$geneSymbol)),"_",ns,"_tf.txt",sep=''));f.tf
    cat(featureNames(tf.eset.sel),file=f.tf,sep='\n')
    f.tfsig<-f.tf
    dir.cur.tfsig <- sprintf('%s/tf',dir.cur)
  }

  if (!is.null(sig.ref)){
    dir.create(file.path(dir.cur,'sig'),recursive = TRUE)
    sig.eset.sel<-subset(eset.sel,Biobase::fData(eset.sel)$geneSymbol%in%sig.ref)
    dim(sig.eset.sel)
    f.sig<-file.path(dir.cur,'sig',paste(grp.tag,"_",nrow(sig.eset.sel),"_",nlevels(factor(Biobase::fData(sig.eset.sel)$geneSymbol)),"_",ns,"_sig.txt",sep=''));f.sig
    cat(featureNames(sig.eset.sel),file=f.sig,sep='\n')
    f.tfsig<-f.sig
    dir.cur.tfsig <- sprintf('%s/sig',dir.cur)
  }
  return(c(dir.cur.tfsig,as.character(grp.tag),f.exp,f.tfsig))
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
                           low.col="#004C99",high.col="#CC0000",plot.title=NULL,
                           xlab_angle=0,xlab_hjust=0.5){

  p <- ggplot(df, aes_string(x= xlab, y= ylab)) +
    theme_classic()+
    geom_point(aes_string(fill=clab, size= slab),color="black",pch=21)+
    scale_fill_gradient2(low=low.col,high=high.col)+
    scale_x_discrete(limits=levels(df[,xlab]))+
    scale_y_discrete(limits=levels(df[,ylab]))+
    theme(panel.grid.major= element_line(colour = "grey",size=0.3),
          panel.grid.minor = element_line(colour = "grey",size=0.3),
          axis.text.x = element_text(size = 12,angle = xlab_angle, hjust = xlab_hjust),
          axis.text.y = element_text(size = 12))
  labs(x = xlab, y = ylab, title = plot.title)
  return(p)
}
#############
#' @title SparseExpressionSet
#' @exportClass SparseExpressionSet
#' @importFrom Biobase ExpressionSet
setClass( "SparseExpressionSet",
          contains = "ExpressionSet",
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c(classVersion("ExpressionSet"), SparseExpressionSet = "1.0.0" ))))

##
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
