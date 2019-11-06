#'@import Biobase ggplot2 kableExtra knitr limma
#'@importFrom reshape2 melt
#'@importFrom ComplexHeatmap HeatmapAnnotation Heatmap
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
# library(kableExtra) #for Rmarkdown
# library(knitr) #for Rmarkdown

#####
#' @title SparseExpressionSet
#' @exportClass SparseExpressionSet
#' @importFrom Biobase ExpressionSet
setClass( "SparseExpressionSet",
            contains = "ExpressionSet",
            prototype = prototype( new( "VersionedBiobase",
              versions = c(classVersion("ExpressionSet"), SparseExpressionSet = "1.0.0" ))))


#' CreateSparseEset
#' @description Create a S4 class which utilize 'ExpressionSet' template yet compatible with sparseMatrix type of assaydata
#' @param data Sparse expression data, could be from either of these class:c('matrix','dgTMatrix','dgCMatrix').Required
#' @param meta.data phenotype data which rownames should be the same as data colnames; Optional; Default as NULL
#' @param feature.data feature data which rownames should be the same as data rownames; Optional; Default as NULL
#' @param add.meta logical; Whether or not calculate extra pheonotype info including total number of UMI,
#' number of non-zero gene for each cell, mitochondrial percentage and spike-in gene expression percentage and store them in pData
#'
#' @return A customized S4 class using ExpressionSet class as prototype
#' @export
CreateSparseEset<-function(data=NULL,meta.data=NULL,feature.data=NULL,add.meta=T){

  if(!class(data)[1]%in%c("dgCMatrix","dgTMatrix","matrix","dgeMatrix")){
    stop("Input format should %in% c( 'matrix','dgTMatrix','dgCMatrix','dgeMatrix,'Matrix')","\n")
    if(class(data)%in%"matrix"){
      data<-as(data,"sparseMatrix")
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
#'
#' @return A list or sparse matrix expression set
#' @export
readscRNAseqData <- function(file,is.10x=TRUE,CreateSparseEset=TRUE, add.meta=F,...){

  if(is.10x){

    data.path <- file
    cat("Reading 10X genomcis data in",data.path, "\n")

    f<-list.files(path=data.path,full.names = T)
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
        colnames(genes)<-c("ensembl","geneSymbol","biotype")
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
    d<- read.delim(file=file,...)
  }
  return(d)
}



#' generateMICAinput
#'
#' @description A utility function that helps generate MICA input from a data matrix with rownames and colnames
#' @param d matrix with colnames as cell/sample info, rownames as gene/feature info
#' @param filename filename of your MICA input file, supported format: txt
#'
#' @export
generateMICAinput <- function(d,filename){

  mica.input <- as.data.frame(as.matrix(Matrix::t(d)))
  mica.input$ID <- colnames(d)
  mica.input <- mica.input[,c("ID",rownames(d))]

  if(length(grep("txt",filename))==0){
    filename<-paste0(filename,".txt")
  }

  cat("Writing MICA input to table...","\n")
  write.table(mica.input,file = filename,sep = "\t",
              row.names = FALSE,col.names = TRUE,quote = FALSE)
  cat("Done.","\n")

}


#' Generate MICA command script for job submission on local or LSF
#'
#' @description Generate command for running MICA locally or on LSF.
#' MICA is a high-performance clustering analysis method that was implemented in python and cwl.
#'
#' @param save_sh_at path to save your MICA command file;
#' @param input_file path to actual input file for MICA, usually use \code{generateMICAinput} to generete file with desired format.
#' @param project_name character, name of your project, will be used for naming output data
#' @param num_cluster a vector or a numerical number, the number of clusters
#' @param output_path path to MICA output file
#' @param host character, whether you want to run MICA pipeline on "lsf" or "local"
#' @param queue character. If host="lsf", which queue to submit your run_MICA.sh, default as NULL
#' @param config_file path to your customized config file where you could define queue, requeseted memory and processor; if use default file, then keep as NULL; defualt as NULL.
#' @param bootstrap number of iterations of k-means process, default as 10.
#' @param dim_reduction_method character, default as "mds". Other supported methods include "pca" and "lpl".
#' @param visualization character, visualization method used for visualize final clustering results. default as "tsne", other supported methods includes "umap"
#' @param perplexity numeical. parameter used in tsne visualization to indicate plot density
#' @param min_dist numerical. parameter used in umap visualization to indicate the distance between neighbours.
#' @param slice_size numerical. MICA sliced original data matrix into smaller pieces to improve computational efficiency,this parameter was used to indicate how many cells/samples should be grouped together as one downsized sub-matrix.
#' Default as 1000.
#'
#' @return A .sh script with runnable command for MICA execution
#' @examples
#' \dontrun{
#' generateMICAcmd<-function(  save_sh_at="./", #path to save shell script
#'                             input_file="MICA_input.txt", #your MICA input file
#'                             project_name="test",
#'                             num_cluster=c(3,4,5), #a vector of numerical number
#'                             output_path="./",
#'                             host="local", #or local
#'                             dim_reduction_method="MDS",
#'                             visualization="tsne")
#'}
#' @export
generateMICAcmd<-function(save_sh_at,
                            input_file,
                            project_name="test",
                            num_cluster=c(3,4,5),
                            output_path,
                            host="lsf",
                            queue="standard",
                            config_file=NULL,
                            bootstrap=10,
                            dim_reduction_method="MDS",
                            visualization="tsne",
                            perplexity=30,
                            min_dist=0.01,
                            slice_size=1000){

  #sanity check
  if(!dir.exists(save_sh_at)) dir.create(save_sh_at)
  if(!dir.exists(output_path)) dir.create(output_path)
  if(!file.exists(input_file)) stop("Input file not found!")
  file.sh<-file.path(save_sh_at,paste0("01_run_MICA_",project_name,'.sh'))
  if(file.exists(file.sh)) stop("File already existed!")

  #add host specific attributes
  if (tolower(host)=="lsf"){
    cat("For lsf usage, if you need to specify your own configuration file,
    please take https://github.com/jyyulab/MICA/blob/master/MICA/config/config_cwlexec.json as your reference.","\n")

    project<-ifelse(is.null(project_name),'',paste('#BSUB -P ',project_name,' \n'))
    queue.bash<-ifelse(is.null(queue),'',paste('#BSUB -q ', queue,' \n'))
    config.bash<-ifelse(is.null(config_file),
                        paste0("-j ", system.file("config_template", "config_cwlexec.json", package = "scMINER"), " "),
                        paste0("-j ", normalizePath(config_file), " "))

    job<-ifelse(is.null(project_name),
          '#BSUB -J MICA',
          paste0('#BSUB -J ','MICA_',project_name,'\n'))

    sh.scminer<-paste0(
      '#!/bin/env bash\n',
      project,
      job,
      '#BSUB -oo ',project_name,'.sh.out \n',
      '#BSUB -eo ',project_name,'.sh.err \n',
      '#BSUB -R \"rusage[mem=2000]\" \n',
      queue.bash,
      "mica lsf ",
      config.bash)

  } else if (tolower(host)=="local") {
    sh.scminer<-paste0(
      '#!/bin/env bash\n' ,
      "mica local ",
      ifelse(is.null(bootstrap),"",paste0("-b ",bootstrap," ")))
  }

  #add generic attributes
  sh.scminer<-paste0(sh.scminer,
    paste0("-i ", normalizePath(input_file), " "),
    paste0("-p ", project_name, " "),
    paste0("-k ",  paste0(num_cluster, collapse="", " ")),
    paste0("-o ",  normalizePath(output_path), " "),
    ifelse(is.null(bootstrap),"",paste0("-b ",bootstrap," ")),
    ifelse(is.null(dim_reduction_method),"",paste0("-dr ",dim_reduction_method," ")),
    ifelse(is.null(visualization),"",paste0("-v ",visualization, " ")),
    ifelse(is.null(perplexity),"",paste0("-pp ",perplexity," ")),
    ifelse(is.null(min_dist),"",paste0("-d ",min_dist," ")),
    ifelse(is.null(slice_size),"",paste0("-sn ", slice_size," "))
  )

  sink(file.sh)
  cat(sh.scminer)
  sink()

  cat(basename(file.sh),'is generated!\n')

}



#' @title readMICAoutput
#'
#' @description Read MICA input and output to create an expressionSet for downstream analysis
#'
#' @param eset a SparseMatrix Eset
#' @param input_file input expression txt file of MICA pipeline
#' @param output_file output ClusterMem.txt file from MICA pipeline
#' @param load_ClusterRes logical, if TRUE, clustering results will be store at pData(eset)$label
#'
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
    cat("Reading Input...","\n")
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



#' generateSJARACNeInput
#'
#' @title Generate SJARACNE input with designed folder structure
#' @description This function helps to generate appropriate input files for SJARACNe pipeline.
#' It can take transcription factor/signaling gene reference from internal(stored in package) or external (manual define)
#'
#'
#' @param input_eset An expressionSet
#' @param ref c("hg", "mm"), could be a manually defined geneSymbol vector
#' @param funcType c("TF","SIG", NULL), if NULL then both TF and SIG will be considered
#' @param wd.src output path
#' @param group_tag name of group for sample identification
#' @return SJARACNe input files for each subgroups
#' @keywords SJARACNe
#' @examples
#' \dontrun{
#' generateSJARACNeInput(input_eset = eset ,ref = "hg",funcType="TF",
#' wd.src = "./",group_tag = "celltype")
#' }
#' @export
generateSJARACNeInput<-function(input_eset,ref=NULL,funcType=NULL,wd.src,group_tag){

  if (!dir.exists(wd.src)) dir.create(wd.src,recursive = T)

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

  if(group_tag%in%colnames(pData(input_eset))){
    groups <- unique(pData(input_eset)[,group_tag])
    for (i in 1:length(groups)){
      grp.tag<-groups[i]; input_eset[,which(pData(input_eset)[,group_tag]==grp.tag)] -> eset.sel
      SJARACNe_filter(eset.sel=eset.sel,tf.ref=tf.ref,sig.ref=sig.ref,wd.src=wd.src,grp.tag=grp.tag)
    }#end for
  }else{
    stop("Lack of group info, please check your group_tag.","\n")
  }#end if

  save(input_eset, file=file.path(wd.src,"Input.eset")) # save input file as expressionSet
}#end function




#' draw.marker.bbp
#'
#' @title Generate visualization for marker scores via bubble plot
#' @description  Marker visualizatoin from known markers/signatures, requires knowledge-based marker list as input
#' @param ref reference dataframe, includes positive or negative markers for different cell types
#' @param input_eset expressionSet/SparseExpressionSet object with clustering membership stored in pData
#' @param group_tag a character, the variable containing clustering label in pData(eset)
#' @param save_plot logical, whether or not save your plot
#' @param width default as 8, inch as unit
#' @param height default as 5, inch as unit
#' @param plot_name plot name, please include plot type
#' @param feature feature type from your reference, should be in colnames(fData(eset))
#' @return A ggplot object
#'
#' @export
draw.marker.bbp<-function(ref = NULL,input_eset,
                              feature='geneSymbol',group_tag="ClusterRes",
                              save_plot = FALSE,
                              width=8, height=5,
                              plot_name="AnnotationBubbleplot.png"){

  #exp<-apply(exprs(eset),2,std)
  #filter reference marker sets

  if (!feature%in%colnames(fData(input_eset))) stop('Please check your feature!')
  colnames(ref)<-c("celltype","markers","weight")
  ref<-dplyr::filter(ref,markers%in%fData(input_eset)[,feature])
  indx<-which(fData(input_eset)[,feature]%in%ref$markers)
  if(length(indx)==0) stop("No genes from the reference list could be found in data!","\n")

  exp<-as.matrix(exprs(input_eset))[indx,]
  rownames(exp)<-fData(input_eset)[,feature][indx]

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
  df_n<-data.frame(label=pData(input_eset)[,group_tag],n_mtx)
  df_n<-aggregate(.~label,df_n,mean)
  library(reshape2)
  df_n_melt<-melt(df_n,id.vars = "label")

  df<-data.frame(label=pData(input_eset)[,group_tag],ac_norm);
  df<-df[,colSums(is.na(df))<nrow(df)];#remove NA columns
  df<-aggregate(.~label,df,mean)
  input<-t(apply(df[,-1],1,scale))#row normalization
  input<-as.data.frame(cbind(df[,1],input))
  rownames(input)<-rownames(df)
  colnames(input)<-colnames(df)
  df_melt<-melt(input,id.vars = "label")


  if(all(df_melt[,c(1,2)]==df_n_melt[,c(1,2)])){

    d<-cbind(df_melt,df_n_melt[,3])
    colnames(d)<-c("Cluster", "CellType", "MarkerScore","ExpressionPercentage")
    d$Cluster<-as.factor(d$Cluster)

    p<-draw.bubblePlot2(df=d, xlab="Cluster",ylab="CellType",
                        clab="MarkerScore",slab="ExpressionPercentage",
                        plot.title="Cell type annotation for each cluster")
  }

  if(save_plot){ggsave(plot = p, filename = plot_name , units="in",
                       width = width,height = height,dpi = 300)}
  return(p)
}



##################
#' GetActivityFromSJARACNe
#'
#' @description Allocate network information from SJARACNe and calculate activity score for each hub genes.
#'
#' @param SJARACNe_output_path Path to SJARACNe output folder(s)
#' @param SJARACNe_input_eset Expressionset that you generate input from
#' @param group_tag a string, group name stored in pData that defines expression matrix separation
#' @param activity.method c("weighted,unweighted), default to "unweighted"
#' @param activity.norm logical, default to TRUE.
#' @param functype character c("tf","sig"); If NULL, both activity from TF and SIG network will be calculated;
#' default as NULL
#' @param save_network_file logical, default to FALSE
#' @param save_path Path to save network file
#'
#' @return An expressionset with activity values
#'
#' @examples
#' \dontrun{
#' acs.12k <- GetActivityFromSJARACNe(
#'              SJARACNe_output_path ="./",
#'              SJARACNe_input_eset = eset.12k,
#'              activity.method="unweighted",
#'              activity.norm=TRUE,
#'              group_tag = "celltype",
#'              save_network_file=TRUE, #default as false, but recommended to be TRUE
#'              save_path="./networks")
#'}
#'
#' @keywords GetActivity
#' @author Chenxi Qian, \email{chenxi.qian@stjude.org}
#'
#'
#' @export
GetActivityFromSJARACNe<-function(SJARACNe_output_path=NA,
							   SJARACNe_input_eset=NA,
							   functype="tf",
							   group_tag=NA,
							   activity.method="unweighted",
							   activity.norm=TRUE,
							   save_network_file=FALSE,
							   save_path=NULL){

  eset<-SJARACNe_input_eset;

  if(!group_tag%in%colnames(pData(eset))){
	  stop('Check your group_tag please.','\n')
	}

  if(!activity.method%in%c("weighted", "unweighted")){
    stop('Check your group_tag please.','\n')
  }

  #retrieve networks
	output.files<-list.files(path=SJARACNe_output_path,
						pattern="consensus_network_ncol_.txt",recursive = TRUE,full.names = TRUE)

	if(length(output.files)==0) stop ("Please check your SJARACNe output path!",'\n')

	if(!is.null(functype)) {
    if (!functype%in%c("tf","sig")) stop("Only accept functype %in% c('tf','sig')","\n")
    else output.files<-output.files[grep(paste0("/",functype,"/"),output.files)]}

  net.names<-gsub(SJARACNe_output_path,"",output.files)
	net.names<-gsub("\\_.*","",net.names);
	net.names<-gsub("[/]","",net.names)
	celltypes<-unique(net.names)

  #initialize actiivty list
	acs_master<-data.frame(ID=NA,stringsAsFactors=FALSE)
	deg_master<-data.frame(ID=NA,stringsAsFactors=FALSE)


	for( i in 1:length(celltypes)){

	    net<-celltypes[i]
      cat("Retrieve Network from ",i,net,"\n")

      TF.table<-NULL
      SIG.table<-NULL

      f<-output.files[grep(paste0("/",net,"_"),output.files)]

      if (length(grep("/tf/",f)!=0))
        {TF.table<-get.network.scMINER(network_file = f[grep("/tf/",f)])}
      if(length(grep("/sig/",f)!=0))
        {SIG.table<-get.network.scMINER(network_file= f[grep("/sig/",f)])}

      if(save_network_file){
        if(!is.null(TF.table)) save(TF.table,file=file.path(save_path,paste0(net,".TF.network")))
        if(!is.null(SIG.table)) save(SIG.table,file=file.path(save_path,paste0(net,".TF.network")))
        cat("Network saved for ", net,"\n")
      }

  	  cat("Calculate Activity for ",net,"!",'\n')
  	  eset.sel<-eset[,pData(eset)[,group_tag]==net]

      acs1<-get_activity(Net = TF.table$network_dat,tag = "TF",normalize=activity.norm,
    					   eset = eset.sel, activity.method = activity.method)

      acs2<-get_activity(Net = SIG.table$network_dat,tag = "SIG",normalize=activity.norm,
                         eset = eset.sel, activity.method = activity.method)

      acs<-t(cbind(acs1,acs2));rm(acs1,acs2)

 	    #update full gene list
 	    acs.ID <- sapply(strsplit(rownames(acs),"_"),"[",1)

	    acs.deg <- data.frame(ID=acs.ID,
	  						Degree=as.numeric(sapply(strsplit(rownames(acs),"_"),"[",2)),
	  						stringsAsFactors=FALSE)

  	  acs.tmp <- acs; rownames(acs.tmp)<-acs.ID

  	  acs_master<-merge(acs_master,acs.tmp,by.x="ID",by.y="row.names",all=TRUE)
  	  deg_master<-merge(deg_master,acs.deg,by="ID",all=TRUE)

  	  colnames(deg_master)[i+1]<-paste0("degree_",net)

  	  cat("Activity Done!!","\n")
  	  rm(acs)

  	  cat('==============================================',"\n")
  	  gc()

	}#end for

	# generate acs expression set
	deg_master<-filter(deg_master,!is.na(ID))
	fd <- data.frame(ID=deg_master$ID,
	                 fn=sapply(strsplit(deg_master$ID,"\\."),"[",1),
	                 FuncType=sapply(strsplit(deg_master$ID,"\\."),"[",2),
	                 deg_master[,-1],stringsAsFactors = FALSE)

  fd <-merge(fd,fData(eset),by.x="fn",by.y="row.names")
  rownames(fd)<-fd$ID

  pd <- pData(eset)

	acs.mtx <- as.matrix(acs_master[,-1])
	rownames(acs.mtx)<- acs_master$ID

	acs.mtx<-acs.mtx[,rownames(pd)]
  acs.mtx<-acs.mtx[-which(is.na(rownames(acs.mtx))),]

  acs.eset<-new("ExpressionSet",phenoData= new("AnnotatedDataFrame",pd),
          featureData=new("AnnotatedDataFrame",fd), annotation="",exprs=as.matrix(acs.mtx))

  return(acs.eset)
}#end activity function



#############################
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


####inner function to calculate activity
get_activity<-function(Net,eset,tag,use.symbol=FALSE,
                       es.method="mean", activity.method="weighted",
                       normalize=TRUE,test=FALSE,sep.symbol="."){
  library(dplyr)
  if(is.null(Net)) {return(NULL)}

  if(!(use.symbol)) src<-unique(Net$source)
  else{ src<-unique(Net$source.symbol)
        if(!"geneSymbol"%in%colnames(fData(eset))) stop("please check your geneSymbol in fData!")
        else { cat("Using geneSymbol to retrieve network")}
  }

  gsc<-vector("list", length(src))
  names(gsc)<-paste(src, tag, sep = sep.symbol)

  exp<-exprs(eset)
  ac<-matrix(NA, nrow=ncol(eset), ncol=length(gsc), dimnames=list(colnames(eset), names(gsc)))

  #z-normalize each sample
  if(normalize){
    exp<-apply(exp,2,std)
    cat("normalized!\n")
  }
  else{
    cat("Non_normalized!\n")
  }


  for(i in 1:length(gsc)){
    #NetBID based geneset
    if(use.symbol){
      tmp<-filter(Net, Net$source.symbol==src[i]);tmp<-tmp[!duplicated(tmp$target.symbol),]
      gsc[[i]]<-unique(as.character(tmp$target.symbol))
    }
    else{
      tmp<-filter(Net, Net$source==src[i]);tmp<-tmp[!duplicated(tmp$target),]
      gsc[[i]]<-unique(as.character(tmp$target))#network from file
    }

    #update the overlap between NetBID based geneset and original expression data
    if(length(intersect(featureNames(eset),gsc[[i]]))==0) stop("Please check your featureNames!")
      else{
        eset.sel<-eset[featureNames(eset)%in%gsc[[i]],]
        gsc[[i]]<-featureNames(eset.sel)
      }

    #n=degree
    n<-length(gsc[[i]])
    name<-paste(colnames(ac)[i], n, sep='_')

    if(n>1){

      if(activity.method == 'unweighted') ac[,i]<-apply(exp[gsc[[i]],], 2, es, es.method)

      else if (activity.method == 'weighted'){

        fd.sel<-data.frame(fn=featureNames(eset.sel),
                           geneSymbol=fData(eset.sel)$geneSymbol,
                           stringsAsFactors=FALSE)

        tmp<-merge(fd.sel,tmp,by.x="fn",by.y="target",sort = FALSE)

        #if(!all(rownames(exp[gsc[[i]],])==fd.sel$fn)) stop("MI is not coordinate with Feature names! \n")

        tmp$p.sign<-sign(as.numeric(tmp$spearman))
        tmp$p.sign[tmp$p.sign == 0]<-1
        tmp$MI.sign<-as.numeric(tmp$MI)*(tmp$p.sign)

        mat<-t(exp[gsc[[i]],])%*%(tmp$MI.sign)
        MI.sum<-sum(tmp$MI)
        ac[,i]<-mat[,1]/MI.sum

        if (test)
        {ac[,i]<-mat[,1]}
      }
      colnames(ac)[i]<-name
    }
  }
  return(ac)
}


#' DAG_ttest
#' @description Inner function to do t.test(pairwise/2case) from activity matrix
#' @param d A vector of gene expression
#' @param group A vector of group information
#'
#' @export
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


#' @title Find differential activity genes from activity matrix
#'
#' @description  \code{get.DA} is a wraper of (\code{DAG_test}, and \code{getDE.limma}),
#'  which helps to conduct two_sided t.test on all genes in specific group VS Others
#'  to find differential activity genes, a table with essential statistics will be outputted.
#'
#' @param input_eset ExpressionSet that stores group information in pData
#' @param group_tag a character string, column name in pData(input_eset) that indicates group info
#' @param group_case NULL(If do get.DA for all group vs others) or a character string (one specific group vs others) of
#'  column name in pData(input_eset) that indicates group info
#' @param group_ctrl NULL(If one vs Others); a character indicate case group if do pairwise analysis
#' @param method a character from c("t.test", "limma"), which method will be used to identify differential activity gene
#' @return output would be a data.frame containing: t.statistics, p.value, log2FC, z.score, and mean Activity value
#'
#' @seealso DAG_ttest; getDE.limma
#' @examples
#' \dontrun{get.DA(eset, group_tag="group")}
#'  # Try find DAG for only group 1
#'  \dontrun{get.DA(eset, group_tag="group",group_case="group1")}
#' @export

get.DA<-function(input_eset=NULL,group_tag="celltype",group_case=NULL, group_ctrl=NULL, method="t.test"){

  d<-data.frame(id = featureNames(input_eset), exprs(input_eset), stringsAsFactors=FALSE)
  rs <- fData(input_eset);rs$id <- d$id

  if(!group_tag%in%colnames(pData(input_eset))) {
    stop('Please check your group_tag.',"\n")}


  if(!is.null(group_case)){
    if(!group_case%in%pData(input_eset)[,group_tag]){
      stop('Please check your group_case.',"\n")
    }

    if(!is.null(group_ctrl)){
      if(!group_case%in%pData(input_eset)[,group_tag]){
        stop('Please check your group_ctrl',"\n")
        }
      input_eset<-input_eset[,which(pData(input_eset)[,group_tag]%in%c(group_case,group_ctrl))]
    }else{
      group_ctrl<-"Others"
    }

    cat("Find differential activity genes for ", group_case ," vs ",group_ctrl, "only!","\n")
    input_eset$da_group <- ifelse(pData(input_eset)[,group_tag]==group_case,"Aim","Ctrl") #label info

    if(method=="t.test"){
      da <- plyr::ddply(d,'id','DAG_ttest',group=input_eset$da_group)
      rs <- merge(rs,da,by="id")}
    else{
      da <- getDE.limma(eset=input_eset,
                                    G1_name=group_case,G0_name = "Others",
                                    G1=colnames(input_eset[,which(input_eset$da_group=="Aim")]),
                                    G0=colnames(input_eset[,which(input_eset$da_group=="Ctrl")]),
                                    verbose=FALSE)
    }
  }else{

    cat('\n','Find Differential Activity TF for all groups!','\n')

    #do all group cases in all
    if (method=="t.test"){

      da.list <- lapply(unique(pData(input_eset)[,group_tag]),function(xx){

        d.label <- ifelse(pData(input_eset)[,group_tag] == xx, "Aim", "Ctrl") #label info

        da <- plyr::ddply(d,'id','DAG_ttest',group=d.label)

        da$pval<-sapply(da$pval,function(xx){ifelse(xx!=0, xx,.Machine$double.xmin)})

        da$Z<-abs(qnorm(da$pval)/2)*sign(da$t)

        da <- da[,setdiff(colnames(da),"MeanAct.Ctrl")]

        colnames(da)[-1] <- paste0(colnames(da)[-1],'_',xx)

        return(da)})

      rs.tmp <- as.data.frame(da.list,stringsAsFactors=FALSE)
      rs.full <- merge(rs,rs.tmp,by='id');rm(rs.tmp)

      rs <- dplyr::select(rs.full,
                        geneSymbol,
                        id:FuncType,
                        starts_with("degree"),
                        starts_with("t"),
                        starts_with("pval"),
                        starts_with("Z"),
                        starts_with("MeanAct"),
                        starts_with("log2FC"))
    }else {
    #use limma
      da.list <- lapply(unique(pData(input_eset)[,group_tag]),function(xx){
        da <- getDE.limma(eset=input_eset,
                                      G1_name=xx,G0_name = "Others",
                                      G1=colnames(input_eset[,which(pData(input_eset)[,group_tag]==xx)]),
                                      G0=colnames(input_eset[,which(pData(input_eset)[,group_tag]!=xx)]),
                                      verbose=FALSE)
        indx<-match(rs$id, da$ID)
        da<-da[indx,]
        colnames(da)[-1]<-paste0(colnames(da)[-1],"_",xx,"VSothers") # rm ID column
      return(da)})

      rs.tmp <- as.data.frame(da.list,stringsAsFactors=FALSE)
      rs.full <- merge(rs,rs.tmp,by.x='id',by.y="ID");rm(rs.tmp);rm(da.list)

      rs <- dplyr::select(rs.full,
                          geneSymbol,
                          id:FuncType,
                          starts_with("degree"),
                          starts_with("Z.statistics"),
                          starts_with("t_"),
                          starts_with("logFC"),
                          starts_with("P.Val"),
                          starts_with("adj.P.Val"),
                          starts_with("Ave."))

    }
  }#end else
  return(rs)
}



#' @title get.Topdrivers
#' @description Help quick pick top master regulators from previous
#' differential activity analysis results
#' @param DAG_result Output table from function FindDAG
#' @param n threshold to pick top master regulators(top n)
#' @param degree_filter filter out drivers with target number less than certain value
#' @param celltype character, output top hits are from which celltype
#' @return A list of top master regulators among different groups
#' @export
get.Topdrivers <- function(DAG_result= DAG_result, n=5, degree_filter=c(50,500), celltype=NULL){

  rownames(DAG_result)<-DAG_result$id
  cols <- colnames(DAG_result)
  if(is.null(celltype)) celltypes<-gsub("^degree_","",cols[grep("^degree_",cols)])
  else { celltypes<-celltype }

  cat("Output top regulators for ",celltypes ,"\n")
  res<-NULL

  for (i in celltypes){
    cat("Top MR for:" ,i,"\n")
    if(!is.null(degree_filter)) {

      DAG_result<-DAG_result[which(DAG_result[,paste0("degree_",i)]> degree_filter[1] & DAG_result[,paste0("degree_",i)] <degree_filter[2]),]
    }

    tcol<-grep(paste0("^t_",i),colnames(DAG_result))

    topDriver<-DAG_result$id[sort(DAG_result[,tcol],decreasing=TRUE,index.return=TRUE,na.last=TRUE)$ix][1:n]

    D2print<-DAG_result[topDriver,c(1,grep(i,colnames(DAG_result)))]


    print(D2print)

    res<- c(res,topDriver)
    cat("==============","\n")
  }
  cat("Done!")
  return(res)
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
#' @export
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
#'
#' @return A Sparse expression set
#'
#' @export
preMICA.filtering <- function(SparseEset,
                              cutoffs,
                              gene_filter=T,
                              nGene_filter=T,
                              nUMI_filter=T,
                              ERCC_filter=T,
                              Mito_filter=T)
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
      gene <- unname((which(fData(SparseEset)$nCells >= cfs$nCell_cutoff)))
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



##################################################################################################
#' @title plot MICA clustering results or other meta variables
#' @description This function helps to generate a ggplot object for phenotypic visualization
#' @param input_eset ExpressionSet that include visualization coordinates in phenotype data
#' @param label Coloring criteria of data points, should be a variable stored in pData(input_eset)
#' @param visualize character, name of visualization method, this will be used as x or y axis label
#' @param X character, column name of x axis
#' @param Y character, column name of y axis
#' @param show_label logical, whether or not to show label on tSNE plot
#' @param title.size numerical, size of plot title, default as 10
#' @param title.name character, title of plot, default as NULL
#' @param pct numerical, size of point, default as 0.5
#'
#' @export
MICAplot<-function(input_eset,
                   label= NULL,visualize=NULL,
                   X=NULL,Y=NULL,show_label=TRUE,
                   title.size=10,title.name="",pct=0.5){

  if(!label%in%colnames(pData(input_eset))){stop("Label name not found in phenotype data!","\n")}

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


#' @title feature_highlighting
#' @description This plot will visualize feature info on scatter plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param ylabel a characterm, title of y axis
#' @param title.size numerical, default as 5
#' @param x cordinates for x axis
#' @param y cordinates for y axis
#' @param visualize character, name of visualization method, this will be used as x or y axis label
#' @param wrap_by character, variable to wrap plot with
#' @param ncol cordinates for y axis
#' @param alpha numerical, default as 0.8
#' @param colors color palette for feature highlighting
#' @param pct.size numrical, point size
#'
#' @export
feature_highlighting<-function(input_eset,target=NULL,
                               feature="geneSymbol",
                               x="X",y="Y",visualize="tSNE",wrap_by=NULL,
                               ylabel="Expression",pct.size=0.8,
                               title.size=15,ncol=4, alpha=0.8,
                               colors=colorRampPalette(c("#E3E3E3", "#BCA2FC","#4900FE"),interpolate="linear")(8)){

  # change it to expr is ok
  input<-as.matrix(exprs(input_eset))
  indx<-which(fData(input_eset)[,feature]%in%target)
  if(length(indx)==0) stop("Target feature not found.")

  gn<-fData(input_eset)[,feature][indx]
  id.vars<-c(x,y,wrap_by)
  projection<-pData(input_eset)[colnames(input),id.vars]

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
       xlab(label = paste0(toString(visualize),"_1"))+
       ylab(label = paste0(toString(visualize),"_2"))+
       labs(color=ylabel)

   return(p)
}


#' @title feature_vlnplot
#' @description This plot will visualize feature info in violin plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param stat a character, whether to plot median or mean by a black dot on violinplot
#' @param group_tag character, which group info to visualize as y axis
#' @param color_by character, which group info to define color
#' @param colors character vector, default as NULL, will use ggplot default color palette
#' @param ylabel a character, title of y axis
#' @param boxplot logical, whether to plot boxplot on violinplot
#' @param title.size numerical, default as 5
#' @param ncol cordinates for y axis
#'
#' @export
feature_vlnplot <- function(input_eset, group_tag="celltype",
                         target=NULL,feature="geneSymbol",color_by="cluster",
                         ylabel="Expression",ncol=3,stat="median",
                         colors=NULL,
                         boxplot=FALSE,title.size=5){

  # extract input information
  input <- exprs(input_eset)
  indx<-which(fData(input_eset)[,feature]%in%target)
  gn<-fData(input_eset)[,feature][indx]

  label <- as.factor(pData(input_eset)[,group_tag])

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

  p <- ggplot(df_melt, aes(x=cluster, y=value))+
       theme_classic()+
       geom_violin(aes_string(fill=color_by),trim=TRUE,scale="width",na.rm = TRUE,size=0.1,width=0.5)

  if(!is.null(stat)){
    if (stat=="median") p <- p + stat_summary(fun.y=median, geom="point", size=1.2, color="black")
    else if (stat=="mean") p <-p + stat_summary(fun.y=mean, geom="point", size=1.2, color="black")
    else cat("Stat not supported, please check your spelling.","\n")}

  if(boxplot) p <- p + geom_boxplot(fill="white",width=0.1,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)

  p <- p + facet_wrap(~variable,scales = "free",ncol = ncol) +
    labs(x="Cluster",y=ylabel)+
    theme(axis.text.x = element_text(size=10),
          plot.title = element_text(size = title.size, face = "bold"),
          strip.background = element_rect(fill="#FFFFFF"))

  if (ylabel=="Activity") { p <- p + geom_boxplot(width=0.2,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)}

  return(p)
}


#' @title Visualize gene expression level on scRNA-seq data via heatmap
#' @description This plot will visualiz feature info in scatter plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature a character, which feature to visualize
#' @param target a character or a character vector indicating feature names
#' @param group_tag a character, label to visualize on the top of heatmap
#' @param name character, name of value visualized in color scale
#' @param cluster_rows logical, if or not cluster rows
#' @param colors color palette
#' @param plot_name character, name of heamap
#' @param save_plot logical, whether to save plots or not
#' @param width numerical
#' @param height numerical
#' @param ... parameter to be passed to ComplexHeatmap::Heatmap
#' @return a ggplot object
#'
#' @export
feature_heatmap <- function(input_eset,target,feature="geneSymbol",
                         group_tag="label",name="log2Exp",
                         save_plot=TRUE,width=4,height=8,
                         cluster_rows=FALSE,
                         colors=rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)),
                         plot_name="GeneHeatmap.png",
                         ...){

  input <- exprs(input_eset)
  gn<-intersect(target,fData(input_eset)[,feature])
  indx<-match(gn,fData(input_eset)[,feature])

  exp<-exprs(input_eset)[indx,]
  rownames(exp)<-gn
  lab<-pData(input_eset)[,group_tag];names(lab) <- sampleNames(input_eset)

  #re-order expressionmatrix and label
  ranks<-names(sort(lab,decreasing = FALSE))
  exp.ordered<-as.matrix(exp[,ranks])
  lab.ordered<-lab[ranks]
  df<-data.frame(scMINER=lab.ordered)

  #Define color annotations
  n<-length(unique(lab.ordered))
  ncols <- scales::hue_pal()(n)
  names(ncols) <- unique(lab.ordered)
  myanndf = HeatmapAnnotation(df = df,col=list(scMINER = ncols))
  mycolors = colors

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
#' @return a ggplot object
#' @export
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



#' draw.scRNAseq.QC
#' @title draw.scRNAseq.QC
#' @description generated a scRNA-seq quality control report in html with Rmarkdown
#' @param SparseEset an SparseEset generated by CreateSparseEset
#' @param project.name a character, project name to print on report
#' @param plot.dir a character, output directory for QC reports
#' @param output.cutoff logical, whether or not return a list of suggested thresholds for filtering
#' @param group a character, a variable name indicate groupping information (stored in pData) to
#' help generate violin plots
#'
#' @return an R markdown QC report and a list of suggested threshold (if specify)
#' @export
draw.scRNAseq.QC<-function(SparseEset,project.name,
                           plot.dir="./QC/",
                           output.cutoff=TRUE,group="group"){
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
             cfs=cfs,
             output.cutoff=output.cutoff,
             group=group))

  return(cfs)
}

