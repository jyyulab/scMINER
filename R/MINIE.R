## ---------------------------
##
## Script name: MINIE (Mutual Information-based Network Inference Engine)
##
## Purpose of script: to reconstruct cell-type-specific GRNs for driver activity inference and target network rewiring analysis
##
## Author: The scMINER software is developed and maintained by the Yu Laboratory @ St. Jude
##
## Date Created: 2023-03-24
##
## ---------------------------
##
## Notes:
##    (1) MINIE takes inputs of a gene expression profile and cell cluster labels.
##    (2) Then MINIE invokes SJARACNe to reconstruct cell-type-specific transcriptional factor and signaling networks.
##    (3) Finally, with the predicted targets of a driver for a cell cluster, MINIE calculates the driver activity by performing a column-wise normalization to ensure each cell is on a similar expression level, followed by averaging the expression of the driver's target genes.
##
## Include following functions:
## get.network.scMINER()
## ConvertNet2List()
## GetActivityFromSJARACNe()
## get_activity()
## getDE.limma()
## get.DA()
## get.Topdrivers()
## combinePvalVector()
## ---------------------------

#' Read SJARACNe Network Result and Return it as List Object(adapted from NetBID2)
#'
#' \code{get.network.scMINER} reads SJARACNe network construction result and returns a list object
#' with network data frame, driver-to-target list and igraph object wrapped inside.
#'
#' In the demo, "consensus_network_ncol_.txt" file will be read and convert into a list object.
#' This list contains three elements, \code{network_data}, \code{target_list} and \code{igraph_obj}.
#' \code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed.
#' \code{target_list} is a driver-to-target list object. Please check details in \code{get_net2target_list}.
#' \code{igraph_obj} is an igraph object used to save this directed and weighted network.
#' Each edge of the network has two attributes, \code{weight} and \code{sign}.
#' \code{weight} is the "MI (mutual information)" value and \code{sign} is the sign of the spearman
#' correlation coefficient (1, positive regulation; -1, negative regulation).
#'
#' @param network_file character, the path for storing network file. For the output of SJAracne, the name of the network file will be "consensus_network_ncol_.txt" under the output directory.
#'
#' @return Return a list containing three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#'
#' @examples
#' TF_file <- system.file('PBMC14KDS_DemoDataSet/SJAR/Bcell_11546_11539_77/tf_final/consensus_network_ncol_.txt',
#'             package = "scMINER")
#' tf.network  <- get.network.scMINER(network_file=TF_file)
#' \dontrun{
#' }
#' @export
get.network.scMINER<- function(network_file=NULL){
  all_input_para <- c('network_file')
  check_res <- sapply(all_input_para,function(x)check_param(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  net_dat      <- read.delim(file=network_file,stringsAsFactors = FALSE)
  target_list  <- ConvertNet2List(net_dat)
  igraph_obj   <- igraph::graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE) ## add edge weight ???
  if('MI' %in% colnames(net_dat)) igraph_obj   <- igraph::set_edge_attr(igraph_obj,'weight',index=igraph::E(igraph_obj),value=net_dat[,'MI'])
  if('spearman' %in% colnames(net_dat)) igraph_obj   <- igraph::set_edge_attr(igraph_obj,'sign',index=igraph::E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}

#' Convert Pairwise Network Data Frame to Driver-to-Target List
#' \code{ConvertNet2List} is a helper function in the \code{get.SJAracne.network}.
#' But if users have their own pairwise gene network files, they can convert it to driver-to-target list object.
#'
#' @param net_dat data.frame, must contain two columns with column names "source" (driver) and "target" (target genes).
#' "MI" (mutual information) and "spearman" (spearman correlation coefficient) columns are optional, but strongly suggested to use.
#' If "MI" and "spearman" columns are missing, errors may occur in some following steps (e.g. es.method='weightedmean' in \code{cal.Activity}).
#'
#' @return Return a list. The names of the list elements are drivers.
#' Each element is a data frame, contains three columns. "target", target gene names;
#' "MI", mutual information; "spearman", spearman correlation coefficient.
#'
#' @examples
#' TF_file <- system.file('PBMC14KDS_DemoDataSet/SJAR/Bcell_11546_11539_77/tf_final/consensus_network_ncol_.txt',
#'             package = "scMINER")
#' tf.network  <- get.network.scMINER(network_file=TF_file)
#' network_list <- ConvertNet2List(tf.network$network_dat)
#' \dontrun{
#' }
#' @export
ConvertNet2List <- function(net_dat=NULL) {
  all_input_para <- c('net_dat')
  check_res <- sapply(all_input_para,function(x)check_param(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  all_source <- base::unique(net_dat$source)
  all_target <- lapply(all_source, function(x) {
    n1 <- net_dat[which(net_dat$source == x), base::intersect(c('target', 'MI', 'spearman'),colnames(net_dat))]
    if(class(n1)=='character') n1 <- data.frame('target'=n1,'MI'=1,'spearman'=1,stringsAsFactors=F)
    n1 <- unique(n1)
    if(length(unique(n1$target))!=length(n1$target)){ ## multiple
      t1 <- table(n1$target)
      w1 <- names(which(t1==1)); w2 <- names(which(t1>1))
      w21 <- n1[which(n1$target %in% w1),]
      w22 <- do.call(rbind,lapply(w2,function(x){
        x1 <- n1[which(n1$target==x),]
        x1 <- x1[which.max(x1$MI),]
      }))
      n1 <- rbind(w21,w22)
    }
    rownames(n1) <- n1$target
    return(n1)
  })
  names(all_target) <- all_source
  return(all_target)
}

#' GetActivityFromSJARACNe
#'
#' @description Allocate network information from SJARACNe and calculate activity score for each hub genes.
#'
#' @param SJARACNe_output_path Path to SJARACNe output folder(s)
#' @param SJARACNe_input_eset Expressionset that you generate input from
#' @param group_name a string, group name stored in Biobase::pData that defines expression matrix separation
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
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' out.dir.SJAR <- system.file('PBMC14KDS_DemoDataSet/SJAR/',
#'                         package = "scMINER")
#' acs.14k.tf <- GetActivityFromSJARACNe(
#'               SJARACNe_output_path = out.dir.SJAR,
#'               SJARACNe_input_eset = pbmc.14k.DS.eset.log2,
#'               activity.method="unweighted",
#'               activity.norm=TRUE,
#'               group_name =  "celltype",
#'               save_network_file=FALSE,
#'               functype="tf",
#'               save_path='.')
#' \dontrun{
#'}
#'
#' @keywords GetActivity
#' @author Chenxi Qian, \email{chenxi.qian@stjude.org}
#' @export
GetActivityFromSJARACNe<-function(SJARACNe_output_path=NA,
                                  SJARACNe_input_eset=NA,
                                  functype="tf",
                                  group_name=NA,
                                  activity.method="unweighted",
                                  activity.norm=TRUE,
                                  save_network_file=FALSE,
                                  save_path=NULL){

  eset<-SJARACNe_input_eset;
  if(!group_name%in%colnames(Biobase::pData(eset))){
    stop('Check your group_name please.','\n')
  }

  if(!activity.method%in%c("weighted", "unweighted")){
    stop('Check your group_name please.','\n')
  }
  #retrieve networks
  output.files<-list.files(path=SJARACNe_output_path,
                           pattern="consensus_network_ncol_.txt",
                           recursive = TRUE,full.names = TRUE)

  if(length(output.files)==0) stop ("Please check your SJARACNe output path!",'\n')
  print(output.files)
  # TODO: "/tf/" is not in sjaracne's output directory anymore. Should not format output.files based on that
  if(!is.null(functype))
    if (!functype%in%c("tf","sig")) stop("Only accept functype %in% c('tf','sig')","\n")
  #else output.files<-output.files[grep(paste0("/",functype,"/"),output.files)]}
  output.files_use <- output.files[grep(functype,output.files)]
  if(length(output.files_use)>0) output.files <- output.files_use
  if(length(output.files)==0) stop ("Please check your SJARACNe output path!",'\n')
  ##
  net.names<-gsub(SJARACNe_output_path,"",output.files)
  #net.names<-gsub("\\_.*","",net.names);
  net.names<-gsub("^/","",net.names);
  net.names <- gsub('^(.*)_\\d+_\\d+_\\d+/.*','\\1',net.names)
  #net.names<-gsub("[/]","",net.names)
  celltypes<-unique(net.names)
  print(celltypes)

  #initialize activity list
  acs_master<-data.frame(ID=NA,stringsAsFactors=FALSE)
  deg_master<-data.frame(ID=NA,stringsAsFactors=FALSE)

  for( i in 1:length(celltypes)){

    net<-celltypes[i]
    cat("Retrieve Network from ",i,net,"\n")

    TF.table<-NULL
    SIG.table<-NULL
    f<-output.files[grep(paste0(net,"_"),output.files)]

    #if (length(grep("/tf/",f)!=0))
    if(functype=="tf"){
      TF.table<-get.network.scMINER(network_file = f[grep("consensus*",f)])
      SIG.table<-NULL
    }
    #if(length(grep("/sig/",f)!=0))
    if(functype=="sig"){
      SIG.table<-get.network.scMINER(network_file= f[grep("consensus*",f)])
      TF.table<-NULL
    }

    if(save_network_file){
      if(!is.null(TF.table)) save(TF.table,file=file.path(save_path,paste0(net,".TF.network.RData")))
      if(!is.null(SIG.table)) save(SIG.table,file=file.path(save_path,paste0(net,".SIG.network.RData")))
      cat("Network saved for ", net,"\n")
    }

    cat("Calculate Activity for ",net,"!",'\n')
    eset.sel<-eset[,which(Biobase::pData(eset)[,group_name]==celltypes[i])]

    acs1<-get_activity(Net = TF.table$network_dat,tag = "TF",normalize=activity.norm,
                       eset = eset.sel, activity.method = activity.method, use.symbol=TRUE)

    acs2<-get_activity(Net = SIG.table$network_dat,tag = "SIG",normalize=activity.norm,
                       eset = eset.sel, activity.method = activity.method, use.symbol=TRUE)

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

  # Was not getting any values when "row.names", changed to "geneSymbol" since
  # "fn" in fd is equivalent to "geneSymbol" in Biobase::fData(eset) -- unsure if this change
  # is 100% correct.
  fd <-merge(fd,Biobase::fData(eset),by.x="fn",by.y="geneSymbol")
  rownames(fd)<-fd$ID

  pd <- Biobase::pData(eset)

  acs.mtx <- as.matrix(acs_master[,-1])
  rownames(acs.mtx)<- acs_master$ID

  acs.mtx<-acs.mtx[,rownames(pd)]
  acs.mtx<-acs.mtx[-which(is.na(rownames(acs.mtx))),]

  acs.eset<-new("ExpressionSet",phenoData= new("AnnotatedDataFrame",pd),
                featureData=new("AnnotatedDataFrame",fd), annotation="",exprs=as.matrix(acs.mtx))

  return(acs.eset)
}#end activity function


#' Calculate activity from network file or gene list
#'
#' @param Net Network data frame
#' @param eset ExpressionSet/SparseExpressionSet with expression data
#' @param tag If network is TF network or SIG network
#' @param genelist A list of signature gene list
#' @param use.symbol logical, in network file, use geneSymbol or use geneID
#' @param feature character, use which feature as ID in Biobase::fData(eset)
#' @param es.method character, which method to use to calculate actiivty value ("mean","maxmean")
#' @param activity.method character, which method to use to estimate activity ("weighted","unweighted")
#' @param normalize logical, if normalize or not
#' @param sep.symbol which symbol to sparate name and tag
#'
#' @details
#' If network object was loaded by get.network.scMINER function, then network dataframe is could be retrieved under network_dat slot.
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#'                         package = "scMINER")
#' load(demo_file)
#' TF_file <- system.file('PBMC14KDS_DemoDataSet/SJAR/Bcell_11546_11539_77/tf_final/consensus_network_ncol_.txt',
#'             package = "scMINER")
#' tf.network  <- get.network.scMINER(network_file=TF_file)
#' acs1<-get_activity(Net = tf.network$network_dat,tag = "TF",
#'          normalize=T,
#'          eset = pbmc.14k.DS.eset.log2,
#'          activity.method = 'unweighted',
#'          use.symbol=TRUE)
#' @return
#' @export
get_activity<-function(Net=NULL,
                       eset,
                       tag=NULL,
                       genelist=NULL,
                       use.symbol=FALSE,
                       feature='geneSymbol',
                       es.method="mean",
                       activity.method="weighted",
                       normalize=TRUE,
                       sep.symbol="."){

  if(!feature%in%colnames(Biobase::fData(eset))) stop("please check your geneSymbol in Biobase::fData!","\n")
  else{cat("Using",feature,"to match targets!","\n")}

  if(!is.null(Net)) {
    if(!(use.symbol)) src<-unique(Net$source)
    else{src<-unique(Net$source.symbol)}
  }else if(!is.null(genelist)){
    if(any(duplicated(names(genelist)))) stop("Please check the names of your gene lists!","\n")
    src<-names(genelist)
  }else {return(NULL)}

  gsc<-vector("list", length(src))
  if(!is.null(tag)) {names(gsc)<-paste(src, tag, sep = sep.symbol)}
  else {names(gsc)<-src}

  exp<-exprs(eset)
  ac<-matrix(NA, nrow=ncol(eset), ncol=length(gsc), dimnames=list(colnames(eset), names(gsc)))
  #z-normalize each sample
  if(normalize){exp<-apply(exp,2,std);cat("normalized!\n")}
  else{cat("Not normalized!", "\n")}

  for(i in 1:length(gsc)){
    #NetBID based geneset
    if(!is.null(Net)){
      if(use.symbol){
        tmp<-filter(Net, Net$source.symbol==src[i]);tmp<-tmp[!duplicated(tmp$target.symbol),]
        gsc[[i]]<-unique(as.character(tmp$target.symbol))
      }
      else{
        tmp<-filter(Net, Net$source==src[i]);tmp<-tmp[!duplicated(tmp$target),]
        gsc[[i]]<-unique(as.character(tmp$target))#network from file
      }
    }else if (!is.null(genelist)){
      gsc[[i]]<-unique(genelist[[i]])
    }

    #update the overlap between NetBID based geneset and original expression data
    if(length(intersect(Biobase::fData(eset)[,feature],gsc[[i]]))==0){
      cat("Genelist",names(gsc)[i], "has no overlap with eset feature names.","\n")
      next
    }
    else{
      eset.sel<-eset[Biobase::fData(eset)[,feature]%in%gsc[[i]],]
      gsc[[i]]<-featureNames(eset.sel)
    }

    #n=degree
    n<-length(gsc[[i]])
    name<-paste(colnames(ac)[i], n, sep='_')

    if(n>1){

      if(activity.method == 'unweighted') ac[,i]<-apply(exp[gsc[[i]],], 2, es, es.method)

      else if (activity.method == 'weighted'){

        fd.sel<-data.frame(fn=featureNames(eset.sel),
                           geneSymbol=Biobase::fData(eset.sel)$geneSymbol,
                           stringsAsFactors=FALSE)

        tmp<-merge(fd.sel,tmp,by.x="fn",by.y="target",sort = FALSE)

        #if(!all(rownames(exp[gsc[[i]],])==fd.sel$fn)) stop("MI is not coordinate with Feature names! \n")

        tmp$p.sign<-sign(as.numeric(tmp$spearman))
        tmp$p.sign[tmp$p.sign == 0]<-1
        tmp$MI.sign<-as.numeric(tmp$MI)*(tmp$p.sign)

        mat<-t(exp[gsc[[i]],])%*%(tmp$MI.sign)
        MI.sum<-sum(tmp$MI)
        ac[,i]<-mat[,1]/MI.sum
      }
      colnames(ac)[i]<-name
    }
  }
  return(ac)
}



#' Differential Expression Analysis and Differential Activity Analysis Between 2 Sample Groups Using Limma
#'
#' \code{getDE.limma} is a function performs differential gene expression analysis and differential driver activity analysis
#' between control group (parameter G0) and experimental group (parameter G1), using limma related functions.
#'
#' @param eset ExpressionSet class object, contains gene expression data or driver activity data.
#' @param G1 a vector of characters, the sample names of experimental group.
#' @param G0 a vecotr of characters, the sample names of control group.
#' @param G1_name character, the name of experimental group (e.g. "Male"). Default is "G1".
#' @param G0_name character, the name of control group (e.g. "Female"). Default is "G0".
#' @param verbose logical, if TRUE, sample names of both groups will be printed. Default is TRUE.
#' @param random_effect a vector of characters, vector or factor specifying a blocking variable.
#' Default is NULL, no random effect will be considered.
#'
#' @return
#' Return a data frame. Rows are genes/drivers, columns are "ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Z-statistics", "Ave.G1" and "Ave.G0".
#' Names of the columns may vary from different group names. Sorted by P-values.
#' @examples
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' each_subtype <- 'G4'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma <- getDE.limma(eset=analysis.par$cal.eset,
#'                                 G1=G1,G0=G0,
#'                                 G1_name=each_subtype,
#'                                 G0_name='other')
#' DA_driver_limma <- getDE.limma(eset=analysis.par$merge.ac.eset,
#'                                 G1=G1,G0=G0,
#'                                 G1_name=each_subtype,
#'                                 G0_name='other')
#' }
getDE.limma <- function(eset=NULL, G1=NULL, G0=NULL,G1_name=NULL,G0_name=NULL,verbose=TRUE,random_effect=NULL) {
  #
  all_input_para <- c('eset','G1','G0','verbose')
  check_res <- sapply(all_input_para,function(x)check_param(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_options('verbose',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  exp_mat <- Biobase::exprs(eset)
  G1 <- base::intersect(G1,colnames(exp_mat))
  G0 <- base::intersect(G0,colnames(exp_mat))
  if(verbose==TRUE){
    print(sprintf('G1:%s', base::paste(G1, collapse = ';')))
    print(sprintf('G0:%s', base::paste(G0, collapse = ';')))
  }
  if(base::length(G1)==0 | base::length(G0)==0){
    message('Too few samples, please check the sample name of G1, G0 and samples in eset !');return(FALSE);
  }
  #
  all_samples <- colnames(Biobase::exprs(eset))
  use_samples <- c(G0, G1)

  new_eset<-eset[,use_samples]
  new_mat  <- Biobase::exprs(new_eset)
  ##
  design.mat <-as.data.frame(matrix(NA, nrow = base::length(use_samples), ncol = 1))
  rownames(design.mat) <-use_samples
  colnames(design.mat) <- 'group'
  design.mat[base::intersect(G0, use_samples), 'group'] <- 'G0'
  design.mat[base::intersect(G1, use_samples), 'group'] <- 'G1'
  #  design <- model.matrix( ~ group + 0, design.mat)
  group <- factor(design.mat$group)
  design <- model.matrix(~0+group);
  colnames(design) <- levels(group); rownames(design) <- colnames(new_mat)

  if(is.null(random_effect)==TRUE){
    fit <- limma::lmFit(new_mat,design)
  }else{
    random_effect <- random_effect[colnames(Biobase::exprs(new_eset))]
    corfit <- limma::duplicateCorrelation(new_eset,design,block=random_effect)
    fit <- limma::lmFit(new_mat,design,block=random_effect,correlation=corfit$consensus)
  }
  contrasts <- limma::makeContrasts(G1-G0,levels=design)
  fit2 <- limma::contrasts.fit(fit,contrasts=contrasts)
  fit2 <- limma::eBayes(fit2,trend=TRUE)
  #summary(decideTests(fit2, method="global"))
  ##
  tT <- limma::topTable(fit2, adjust.method = "fdr", number = Inf,coef=1)
  if(nrow(tT)==1){
    rownames(tT) <- rownames(new_mat)
  }
  tT <- base::cbind(ID=rownames(tT),tT,stringsAsFactors=FALSE)
  tT <- tT[rownames(new_mat),,drop=FALSE]
  exp_G1 <- base::rowMeans(new_mat[,G1,drop=FALSE]);
  exp_G0 <- base::rowMeans(new_mat[,G0,drop=FALSE]);
  w1 <- which(tT$P.Value<=0);
  if(base::length(w1)>0) tT$P.Value[w1] <- .Machine$double.xmin;
  z_val <- sapply(tT$P.Value*sign(tT$logFC),function(x)combinePvalVector(x,twosided = TRUE)[1])
  if(is.null(random_effect)==TRUE){
    tT <- base::cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  }else{
    tT <- base::cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1,
                      'Ave.G0_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G0'],
                      'Ave.G1_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G1'])
  }
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  tT <- tT[order(tT$P.Value, decreasing = FALSE), ]
  return(tT)
}

#' @title Find differential activity genes from activity matrix
#'
#' @description  \code{get.DA} is a wraper of (\code{DAG_test}, and \code{getDE.limma}),
#'  which helps to conduct two_sided t.test on all genes in specific group VS Others
#'  to find differential activity genes, a table with essential statistics will be outputted.
#'
#' @param input_eset ExpressionSet that stores group information in Biobase::pData
#' @param group_name a character string, column name in Biobase::pData(input_eset) that indicates group info
#' @param group_case NULL(If do get.DA for all group vs others) or a character string (one specific group vs others) of
#'  column name in Biobase::pData(input_eset) that indicates group info
#' @param group_ctrl NULL(If one vs Others); a character indicate case group if do pairwise analysis
#' @param method a character from c("t.test", "limma"), which method will be used to identify differential activity gene
#' @return output would be a data.frame containing: t.statistics, p.value, log2FC, z.score, and mean Activity value
#'
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/SJAR/DATA/celltype_Activity.RData',
#'                                              package = "scMINER")
#' load(demo_file)
#' DAG_result_tf <- get.DA(input_eset = AC_eset$AC.TF,
#'                        group_name = "celltype")
#' @export
get.DA<-function(input_eset=NULL,group_name="celltype",group_case=NULL, group_ctrl=NULL, method="t.test"){

  d<-data.frame(id = featureNames(input_eset), exprs(input_eset), stringsAsFactors=FALSE)
  rs <- fData(input_eset);rs$id <- d$id

  if(!group_name%in%colnames(pData(input_eset))) {
    stop('Please check your group_name.',"\n")}


  if(!is.null(group_case)){
    if(!group_case%in%pData(input_eset)[,group_name]){
      stop('Please check your group_case.',"\n")
    }

    if(!is.null(group_ctrl)){
      if(!group_case%in%pData(input_eset)[,group_name]){
        stop('Please check your group_ctrl',"\n")
      }
      input_eset<-input_eset[,which(pData(input_eset)[,group_name]%in%c(group_case,group_ctrl))]
    }else{
      group_ctrl<-"Others"
    }

    cat("Find differential activity genes for ", group_case ," vs ",group_ctrl, "only!","\n")
    input_eset$da_group <- ifelse(pData(input_eset)[,group_name]==group_case,"Aim","Ctrl") #label info

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

      da.list <- lapply(unique(pData(input_eset)[,group_name]),function(xx){

        d.label <- ifelse(pData(input_eset)[,group_name] == xx, "Aim", "Ctrl") #label info

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
      da.list <- lapply(unique(pData(input_eset)[,group_name]),function(xx){
        da <- getDE.limma(eset=input_eset,
                          G1_name=xx,G0_name = "Others",
                          G1=colnames(input_eset[,which(pData(input_eset)[,group_name]==xx)]),
                          G0=colnames(input_eset[,which(pData(input_eset)[,group_name]!=xx)]),
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
#' @examples
#' demo_file <- system.file('PBMC14KDS_DemoDataSet/SJAR/DATA/celltype_Activity.RData',
#'                                              package = "scMINER")
#' load(demo_file)
#' DAG_result_tf <- get.DA(input_eset = AC_eset$AC.TF,
#'                        group_name = "celltype")
#' celltype <- levels(Biobase::pData(AC_eset$AC.TF)[,"celltype"])
#' TF_list <- get.Topdrivers(DAG_result = DAG_result_tf,
#'                           celltype = celltype,
#'                           n = 5, degree_filter = c(50, 600))
#' @export
get.Topdrivers <- function(DAG_result= DAG_result, n=5, degree_filter=c(50,500), celltype=NULL){

  rownames(DAG_result)<-DAG_result$id
  cols <- colnames(DAG_result)
  if(is.null(celltype)){
    celltypes<-gsub("^degree_","",cols[grep("^degree_",cols)])
  } else {
    celltypes<-celltype
  }
  ori_DAG_result<-DAG_result
  cat("Output top regulators for ",celltypes ,"\n")
  res<-NULL
  for (i in celltypes){
    cat("Top MR for:" ,i,"\n")
    if(!is.null(degree_filter)) {
      DAG_result<-ori_DAG_result[which(ori_DAG_result[,paste0("degree_",i)]> degree_filter[1] & ori_DAG_result[,paste0("degree_",i)] <degree_filter[2]),] ### change to ori_DAG_result
    }
    tcol<-grep(paste0("^t_",i),colnames(DAG_result))
    topDriver<-DAG_result$id[sort(DAG_result[,tcol],decreasing=TRUE,index.return=TRUE,na.last=TRUE)$ix][1:n]

    D2print<-DAG_result[topDriver,c(1,grep(i,colnames(DAG_result)))]
    D2print <- D2print[which(is.na(D2print[,1])==F),,drop=F]
    print(D2print)

    res<- c(res,topDriver[complete.cases(topDriver)])
    cat("==============","\n")
  }
  cat("Done!")
  return(res)
}

#' Combine P Values Using Fisher's Method or Stouffer's Method
#'
#' \code{combinePvalVector} is a function to combine multiple comparison's P values using Fisher's method or Stouffer's method.
#'
#' @param pvals a vector of numerics, the P values from multiple comparison need to be combined.
#' @param method character, users can choose between "Stouffer" and "Fisher". Default is "Stouffer".
#' @param signed logical, if TRUE, will give a sign to the P value to indicate the direction of testing.
#' Default is TRUE.
#' @param twosided logical, if TRUE, P value is calculated in a one-tailed test.
#' If FALSE, P value is calculated in a two-tailed test, and it falls within the range 0 to 0.5.
#' Default is TRUE.
#' @return Return a vector contains the "Z-statistics" and "P.Value".
#' @examples
#' combinePvalVector(c(0.1,1e-3,1e-5))
#' combinePvalVector(c(0.1,1e-3,-1e-5))
#' @export
combinePvalVector <-function(pvals,
                             method = 'Stouffer',
                             signed = TRUE,
                             twosided = TRUE) {

  #remove NA pvalues
  pvals <- pvals[!is.na(pvals) & !is.null(pvals)]
  pvals[which(abs(pvals)<=0)] <- .Machine$double.xmin
  if (sum(is.na(pvals)) >= 1) {
    stat <- NA
    pval <- NA
  } else{
    if (twosided & (sum(pvals > 1 | pvals < -1) >= 1))
      stop('pvalues must between 0 and 1!\n')
    if (!twosided & (sum(pvals > 0.5 | pvals < -0.5) >= 1))
      stop('One-sided pvalues must between 0 and 0.5!\n')

    if (!signed) {
      pvals <- abs(pvals)
    }

    signs <- sign(pvals)
    signs[signs == 0] <- 1

    if (grepl('Fisher', method, ignore.case = TRUE)) {
      if (twosided & signed) {
        neg.pvals <- pos.pvals <- abs(pvals) / 2
        pos.pvals[signs < 0] <- 1 - pos.pvals[signs < 0]
        neg.pvals[signs > 0] <- 1 - neg.pvals[signs > 0]
      } else{
        neg.pvals <- pos.pvals <- abs(pvals)
      }
      pvals <-
        c(1, -1) * c(
          pchisq(
            -2 * sum(log(as.numeric(pos.pvals))),
            df = 2 * base::length(pvals),
            lower.tail = FALSE
          ) / 2,
          pchisq(
            -2 * sum(log(as.numeric(neg.pvals))),
            df = 2 * base::length(pvals),
            lower.tail = FALSE
          ) / 2
        )
      pval <- base::min(abs(pvals))[1]
      #if two pvals are equal, pick up the first one
      stat <-
        sign(pvals[abs(pvals) == pval])[1] * qnorm(pval, lower.tail = F)[1]
      pval <- 2 * pval
    }
    else if (grepl('Stou', method, ignore.case = TRUE)) {
      if (twosided) {
        zs <- signs * qnorm(abs(pvals) / 2, lower.tail = FALSE)
        stat <- sum(zs) / sqrt(base::length(zs))
        pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
      }
      else{
        zs <- signs * qnorm(abs(pvals), lower.tail = FALSE)
        stat <- sum(zs) / sqrt(base::length(zs))
        pval <- pnorm(abs(stat), lower.tail = FALSE)
      }
    }
    else{
      stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
    }
  }
  return(c(`Z-statistics` = stat, `P.Value` = pval))
}


###Activity inner function###
std <- function(x){
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
#' @export
DAG_ttest<-function(d,group){

  #d.tmp<-unlist(d[1,-1])
  d.tmp<-unlist(d[,-1])
  d.tmp<-data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)

  d.sel<-d.tmp[complete.cases(d.tmp),]
  d.sel$group<-as.factor(d.sel$group)

  n.cases<-nlevels(d.sel$group)

  #cat("=")
  if(n.cases==2){
    res <- t.test(dplyr::filter(d.tmp,group=="Aim")$acs,
                  dplyr::filter(d.tmp,group=="Ctrl")$acs)
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
