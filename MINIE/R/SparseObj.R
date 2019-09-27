
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

  if(!class(data)[1]%in%c("dgCMatrix","dgTMatrix","matrix")){
    stop("Input format should %in% c( 'matrix','dgTMatrix','dgCMatrix')","\n")
    if(class(data)%in%"matrix")
      data<-as(data,"sparseMatrix")
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

    if(length(mito.genes)==0) Mito_filter=FALSE
    if(length(spikeIn.genes)==0) ERCC_filter=FALSE

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
