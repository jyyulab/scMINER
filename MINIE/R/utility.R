##read data, a wrapper of conventional data reading (read.delim2) and 10x genomics data reading
#' @export
readscRNAseqData <- function(file,is.10x=TRUE,...){

  if(is.10x){
    data.path <- file
    data.raw <- Matrix::readMM(file.path(data.path,"matrix.mtx"))


    barcodes <- read.table(file.path(data.path,"barcodes.tsv"),header=FALSE,stringsAsFactors = FALSE)
    
    if(file.exists(file.path(data.path,"genes.tsv"))){
      genes <- read.table(file.path(data.path,"genes.tsv"),header=FALSE,stringsAsFactors = FALSE)
      colnames(genes)<-c("ensembl","geneSymbol")
    }else if(file.exists(file.path(data.path, "feature.tsv"))){
      genes <- read.table(file.path(data.path, "feature.tsv"),header=FALSE,stringsAsFactors=FALSE)
      colnames(genes)<-c("ensembl","geneSymbol","biotype")
    }

    dimnames(data.raw)[[1]] <- genes$V1
    dimnames(data.raw)[[2]] <- barcodes$V1

    d=list(raw.data=data.raw,
           meta.data=barcodes,
           feature.data=genes)
  }
  else{
    d<- read.delim2(file=file,...)
  }
  return(d)
}




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

  res <- read.table( output_file, # MICA output text file
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
                 X=res[,2],
                 Y=res[,3],stringsAsFactors=FALSE)

  if(load_clust_label){pd$label <- as.factor(res$label);
  cat("Clustering info is under 'label' slot.","\n")}

  eset<-new("ExpressionSet",phenoData= new("AnnotatedDataFrame",pd),
            featureData=new("AnnotatedDataFrame",fd), annotation="",exprs=as.matrix(input))
  cat("ExpressionSet Generated!","\n")
  return(eset)
}



#' generateMICAinput
#' @description This function helps generate MICA input from a data matrix with rownames and colnames
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
#' @description Generate command for running MICA locally or on LSF.
#' MICA is a high-performance clustering analysis method that was implemented in python and cwl.
#' @usage
#' @param save_sh_at character, path to save your MICA command file
#' @param input_file character, path to actual input file for MICA, usually use \code{generateMICAinput} to generete file with desired format.
#' @param project_name character, name of your project, will be used for naming of output data
#' @param num_cluster a vector or a numerical number, the number of clusters
#' @param output_path character, path to MICA output file
#' @param host character, whether you want to run MICA pipeline on "lsf" or "local"
#' @param queue character. If host="lsf", which queue to submit your job, default as NULL
#' @param memory a vector of numerical number, default as NULL
#' @param threads number of threads for pooling in clustering step, default as 10.
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
#' generate_MICA_cmd<-function("./cmd", #path to save shell script
#'                             "MICA_input.txt", #your MICA input file
#'                             "test",
#'                             c(3,4,5), #a vector of numerical number
#'                             "./",
#'                             host="lsf", #or local
#'                             queue=NULL, #your queue to submit the job
#'                             memory=NULL, #specify if you use LSF
#'                             dim_reduction_method="MDS",
#'                             visualization="tsne")
#'}
#' @export
generate_MICA_cmd<-function(save_sh_at,
                            input_file,
                            project_name,
                            num_cluster,
                            output_path,
                            host="lsf",
                            queue="standard",
                            memory=NULL,
                            threads=10,
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

  sink(file.sh)

  if (tolower(host)=="lsf"){
    project<-ifelse(is.null(project_name),'',paste('#BSUB -P ',project_name,' \n'))
    queue<-ifelse(is.null(queue),'',paste('#BSUB -q ', queue,' \n'))
    
    job<-ifelse(is.null(project_name),
          '#BSUB -J MICA',
          paste('#BSUB -J ','MICA_',project_name,',\n'))

    sh.scminer<-paste0(
      '#!/bin/env bash\n',
      project,
      job,
      '#BSUB -oo ',project_name,'.sh.out \n',
      '#BSUB -eo ',project_name,'.sh.err \n',
      '#BSUB -R \"rusage[mem=2000]\" \n',
      queue,
      "mica lsf ",
      ifelse(is.null(memory),"",paste0("-r ",paste0(memory, collapse = ""), " ")))
      ifelse(is.null(queue),"",paste0("-q ", queue, " "))
  } else if (tolower(host)=="local") {
    sh.scminer<-paste0(
      '#!/bin/env bash\n' ,
      "mica local ")
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
    ifelse(is.null(slice_size),"",paste0("-sn ", slice_size," ")),
    ifelse(is.null(threads),"",paste0("-t ", threads," "))
  )
  cat(sh.scminer)
  sink()

  cat(basename(file.sh),'is generated!\n')

}



###Function7: AssignCelltypes###
#' AssignCellTypes.bbp
#' @description  Cell type annotation from known markers/signatures
#'
#' @param ref reference dataframe, includes positive or negative markers for different cell types
#' @param eset expressionSet with clustering membership stored in pData
#' @param save_plot logical
#' @param width default as 8
#' @param height default as 5
#' @param plot_name plot name
#'
#' @return A ggplot object
#'
#' @export
AssignCellTypes.bbp<-function(ref = NULL,eset = eset.demo,
                              save_plot = FALSE,
                              width=8, height=5,
                              plot_name="AnnotationBubbleplot.png"){

  #exp<-apply(exprs(eset),2,std)
  #filter reference marker sets
  exp<-exprs(eset)
  ref<-filter(ref,markers%in%rownames(exp))
  celltypes<-unique(ref$celltype)

  ac<-matrix(NA,nrow=ncol(exp),ncol=length(celltypes),dimnames = list(colnames(exp),celltypes))
  for(i in 1:length(celltypes)){
    cat(i,"\n")
    ref.sel<-filter(ref,celltype==celltypes[i])
    n <- length(ref.sel$markers)

    if(n>1){
      mat<-t(exp[ref.sel$markers,])%*%as.vector(ref.sel$weight)
      ac[,i]<-mat[,1]/n
    }else if (n==1){
      ac[,i]<-exp[ref.sel$markers,]
    }
  }

  ac_norm<-apply(ac,2,scale) #column normalization

  n_mtx<-(ac>0.5)
  df_n<-data.frame(label=eset$label,n_mtx)
  df_n<-aggregate(.~label,df_n,mean)
  library(reshape2)
  df_n_melt<-melt(df_n,id.vars = "label")

  df<-data.frame(label=eset$label,ac_norm);
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

  if(save_plot){ggsave(plot = p, filename = plot_name ,
                       device="png",width = width,height = height,dpi = 300)}
  return(p)
}


