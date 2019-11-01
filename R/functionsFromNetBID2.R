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


#' Convert Pairwise Network Data Frame to Driver-to-Target List
#'
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
#' tf.network.file <- sprintf('%s/demo1/network/SJAR/project_2019-02-14/%s/%s',
#'                    system.file(package = "NetBID2"),
#'                    'output_tf_sjaracne_project_2019-02-14_out_.final',
#'                    'consensus_network_ncol_.txt')
#' net_dat      <- read.delim(file=tf.network.file,stringsAsFactors = FALSE)
#' target_list  <- ConvertNet2List(net_dat)
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
#' \dontrun{
#' tf.network  <- get.network.scMINER(network_file="consensus_network_ncol_.txt")
#' }
#' @export
get.network.scMINER<- function(network_file=NULL){
  all_input_para <- c('network_file')
  check_res <- sapply(all_input_para,function(x)check_param(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  net_dat      <- read.delim(file=network_file,stringsAsFactors = FALSE)
  target_list  <- ConvertNet2List(net_dat)
  igraph_obj   <- igraph::graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE) ## add edge weight ???
  if('MI' %in% colnames(net_dat)) igraph_obj   <- igraph::set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  if('spearman' %in% colnames(net_dat)) igraph_obj   <- igraph::set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
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
#'
#'
#' @examples
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
#' @export
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
  phe <- as.data.frame(Biobase::pData(eset)[use_samples, ,drop=FALSE]);
  rownames(phe) <- use_samples
  new_eset <- generate.eset(exp_mat=Biobase::exprs(eset)[, use_samples,drop=F],phenotype_info=phe,
                            feature_info=Biobase::fData(eset))
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
    rownames(tΩT) <- rownames(new_mat)
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
