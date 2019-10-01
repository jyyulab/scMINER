#' readscRNAseqData
#' @description read scRNA-seq data, a wrapper of conventional data reading (read.delim2) and 10x genomics data standarad output reading
#'
#' @param file data path to 10x genomics output folder, which normally contains 3 files (matrix.mtx, gene or feature.tsv and barcode.csv),
#'  or data path to data txt/csv/tsv file
#' @param is.10x logical, whether or not inputs are from CellRanger standard output
#' @param CreateSparseEset logical, whether or not create sparse matrix incorporated expression set
#' @param add.meta logical, whether or not calculate metadata info from expression matrixm, this is not suggested before merging/downsampling your data
#' @param ... paramters pass to read.delim if is.10x = FALSE
#'
#' @return A list or sparse matrix expression set
#' @export
readscRNAseqData <- function(file,is.10x=TRUE,CreateSparseEset=TRUE, add.meta=F,...){

  if(is.10x){
    data.path <- file
    data.raw <- Matrix::readMM(file.path(data.path,"matrix.mtx"))

    barcodes <- read.table(file.path(data.path,"barcodes.tsv"),header=FALSE,stringsAsFactors = FALSE)
    colnames(barcodes)[1]<-"CellNames"
    rownames(barcodes)<-barcodes[,1]

    if(file.exists(file.path(data.path,"genes.tsv"))){
      genes <- read.table(file.path(data.path,"genes.tsv"),header=FALSE,stringsAsFactors = FALSE)
      colnames(genes)<-c("ensembl","geneSymbol")
      rownames(genes)<-genes[,1]
    }else if(file.exists(file.path(data.path, "feature.tsv"))){
      genes <- read.table(file.path(data.path, "feature.tsv"),header=FALSE,stringsAsFactors=FALSE)
      colnames(genes)<-c("ensembl","geneSymbol","biotype")
      rownames(genes)<-genes[,1]
    }


    dimnames(data.raw)[[1]] <- rownames(genes)
    dimnames(data.raw)[[2]] <- rownames(barcodes)

    if (CreateSparseEset){
      d <- CreateSparseEset(data=data.raw,meta.data=barcodes,
        feature.data=genes,add.meta = add.meta)
    }else{
      d = list(raw.data=data.raw,
           meta.data=barcodes,
           feature.data=genes)}
  }
  else{
    d<- read.delim2(file=file,...)
  }
  return(d)
}



#' generateMICAinput
#'
#' @description This utility function helps generate MICA input from a data matrix with rownames and colnames
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
#' @param save_sh_at character, path to save your MICA command file;
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
#' generateMICAcmd<-function("./cmd", #path to save shell script
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
generateMICAcmd<-function(save_sh_at,
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
    queue.bash<-ifelse(is.null(queue),'',paste('#BSUB -q ', queue,' \n'))

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
      ifelse(is.null(memory),"",paste0("-r ",paste0(memory, collapse = ""), " ")),
      ifelse(is.null(queue),"",paste0("-q ", queue, " ")))
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



#' @title readMICAoutput
#'
#' @description Read MICA input and output to create an expressionSet for downstream analysis
#'
#'
#' @param eset a SparseMatrix Eset
#' @param input_file input expression txt file of MICA pipeline
#' @param output_file output ClusterMem.txt file from MICA pipeline
#' @param load_clust_label logical, if TRUE, clustering results will be store at pData(eset)$label
#' @param NewSparseEset logical, if TRUE, return a eset obj
#'
#'
#' @return An expressionSet
#' @export
readMICAoutput<-function(eset=NULL, input_file, output_file,load_ClusterRes=TRUE){

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
    d <- read.table(input_file,header = FALSE,
                    stringsAsFactors = FALSE,
                    quote = "",
                    check.names = FALSE)

    gn<-unname(unlist(d[1,-1]))
    input<-t(d[-1,-1]);colnames(input)<-d$V1[-1]
    class(input)<-"numeric"
    rownames(input)<-gn;

    pd<-data.frame(row.names=colnames(input),cellNames=colnames(input),
                   X=res[,2],Y=res[,3],stringsAsFactors=FALSE)

    eset<-CreateSparseEset(data=input,meta.data=pd)
    cat("SparseMatrix Expression Set Generated!","\n")
  }

  if(load_clust_label){eset$ClusterRes <- as.factor(res$label);
  cat("Clustering info is under 'ClusterRes' slot.","\n")}

}



#' generateSJARACNeInput
#'
#'
#' @description This function helps to generate appropriate input files for SJARACNe pipeline.
#' It can take transcription factor/signaling gene reference from internal(stored in package) or external (manual define)
#'
#'
#' @param eset scRNA-seq ExpressionSet
#' @param ref c("hg", "mm"), could be a manually defined geneSymbol vector
#' @param funcType c("TF","SIG", NULL), if NULL then both TF and SIG will be considered
#' @param wd.src output path
#' @param group_tag name of group for sample identification
#' @return SJARACNe input files for each subgroups
#' @keywords SJARACNe
#' @examples
#' \dontrun{
#' generateSJARACNeInput(eset = eset.demo,ref = "hg",wd.src = "./",group_tag = "celltype")}
#' @export
generateSJARACNeInput<-function(eset,ref=NULL,funcType=NULL,wd.src,group_tag){

  if (!dir.exists(wd.src)) dir.create(wd.src,recursive = T)

  if (ref%in%c("hg","mm")){
    ref_file<-system.file("RData",paste0("tf_sigs_",ref,".RData"),package = "MINIE")
    load(ref_file)
    cat("Using references from: ", ref_file,"\n")
    sig.ref <- NULL;tf.ref <- NULL
    tf.ref<- filter(tf_sigs, isTF==TRUE)$geneSymbol
    sig.ref<- filter(tf_sigs, isSIG==TRUE)$geneSymbol
    if (!is.null(funcType)){
      if (funcType=="TF") sig.ref <- NULL
      else if (funcType=="SIG") tf.ref <- NULL
    }
  }else {
    if (funcType=="TF") tf.ref <- ref
    else if (funcType=="SIG") sig.ref <- ref
    else warning("Activity calculations will not be supported!","\n")
  }

  if(group_tag%in%colnames(pData(eset))){
    groups <- unique(pData(eset)[,group_tag])
    for (i in 1:length(groups)){
      grp.tag<-groups[i]; eset[,which(pData(eset)[,group_tag]==grp.tag)] -> eset.sel
      SJARACNe_filter(eset.sel=eset.sel,tf.ref=tf.ref,sig.ref=sig.ref,wd.src=wd.src,grp.tag=grp.tag)
    }#end for
  }else{
    stop("Lack of group info, please check your group_tag.","\n")
  }#end if

  save(eset, file=file.path(wd.src,"Input.eset")) # save input file as expressionSet
}#end function




#' @title marker_bbplot
#' @description  Marker visualizatoin from known markers/signatures, requires knowledge-based marker list as input
#' @param ref reference dataframe, includes positive or negative markers for different cell types
#' @param eset expressionSet/SparseExpressionSet object with clustering membership stored in pData
#' @param save_plot logical, whether or not save your plot
#' @param width default as 8, inch as unit
#' @param height default as 5, inch as unit
#' @param plot_name plot name, please include plot type
#'
#' @return A ggplot object
#'
#' @export
marker_bbplot<-function(ref = NULL,eset = eset.demo,feature='geneSymbol',
                              save_plot = FALSE,
                              width=8, height=5,
                              plot_name="AnnotationBubbleplot.png"){

  #exp<-apply(exprs(eset),2,std)
  #filter reference marker sets

  if (!feature%in%colnames(fData(eset))) stop('Please check your feature!')

  ref<-filter(ref,markers%in%fData(eset)[,feature])
  indx<-which(fData(eset)[,feature]%in%ref$markers)

  exp<-exprs(eset[indx,])
  rownames(exp)<-fData(eset)[,feature][indx]

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

  if(save_plot){ggsave(plot = p, filename = plot_name , unit="in",
                       width = width,height = height,dpi = 300)}
  return(p)
}
