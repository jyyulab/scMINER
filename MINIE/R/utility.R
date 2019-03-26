#' @title readMICAoutput
#' @description Read MICA input and output to create an expressionSet for downstream analysis
#' @param input_file input txt file for MICA pipeline
#' @param output_file output ggplot.txt file from MICA pipeline
#' @param load_clust_label logical, if TRUE, clustering results will be store at pData(eset)$label
#' @return An expressionSet
#' @export
readMICAoutput<-function(input_file, output_file,load_clust_label=TRUE){

  cat("Reading Input...","\n")
  d <- read.table(input_file,header = FALSE,
                  stringsAsFactors = FALSE,
                  quote = "",
                  check.names = FALSE)

  MICA_Result <- read.table( output_file, # MICA output text file
                             header = TRUE,
                             stringsAsFactors = FALSE)

  gn<-unname(unlist(d[1,-1]))
  input<-t(d[-1,-1]);colnames(input)<-d$V1[-1]
  class(input)<-"numeric"

  if(length(which(duplicated(gn)))!=0){
    cat("Colnames has duplicates, assigned new rownames. Gene symbols stored in fData.","\n")
  }else{
    rownames(input)<-gn;
    cat("Gene Symbol as rownames.","\n")}

  fd<-data.frame(row.names=rownames(input),
                 geneSymbol=gn,stringsAsFactors = FALSE)

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


#' generateMICAinput
#' @description function to generate MICA input from pre.MICA output
#' @param d matrix with colnames as cell info, rownames as gene info
#' @param filename filename of your MICA input file, supported format: txt
#'
#' @export
generateMICAinput <- function(d,filename){

  mica.input <- as.data.frame(t(d))
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






