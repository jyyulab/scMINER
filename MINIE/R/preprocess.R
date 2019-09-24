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
  expdata<-data.frame(cbind(isoformId=featureNames(eset.sel),geneSymbol=fData(eset.sel)$geneSymbol,exprs(eset.sel)))
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
    f.sig<-file.path(dir.cur,'sig',paste(grp.tag,"_",nrow(sig.eset.sel),"_",nlevels(factor(fData(tf.eset.sel)$geneSymbol)),"_",ns,"_sig.txt",sep=''));f.sig
    cat(featureNames(sig.eset.sel),file=f.sig,sep='\n')
  }

}

#' generateSJARACNeInput
#' @description This function helps to generate appropriate input files for SJARACNe pipeline.
#' It can take transcription factor/signaling gene reference from internal(stored in package) or external (manual define)
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



#' preMICA.filtering
#' @title preMICA.filtering
#' @description Qualtiy control function to do the actual filtering
#' @param Obj the list object outputted by preMICA.QA
#' @param project.name character
#' @param gene_filter logical
#' @param cell_filter logical
#' @param ERCC_filter logical
#' @param Mito_filter logical
#' @param nGene_filter logical
#' @param nUMI_filter "high","low" or "both"
#'
#' @return a list contains raw data, filtered data, feature data, meta data, and thresholds.
#' @export
#'
preMICA.filtering <- function(SparseEset,
                        project.name=NULL,
                        gene_filter=TRUE,
                        ERCC_filter=TRUE,
                        Mito_filter=TRUE,
                        nGene_filter=TRUE,
                        nUMI_filter="low"){

    cfs<-Obj$cal.cutoffs
    d<-Obj$raw.data

    if(nUMI_filter=="low") {cfs$umi_cf_hi = Inf
    }else if(nUMI_filter=="high") {cfs$umi_cf_lo = 0}

    cell <- which((SparseEset$nUMI.total > cfs$umi_cf_lo) & (SparseEset$nUMI.total < cfs$umi_cf_hi))
    if(nGene_filter) cell <- intersect(cell, which(SparseEset$nGene > cfs$nGene_cf))
    if(ERCC_filter&cfs$ERCC_cf!=0) cell <- intersect(cell, which(Obj$meta.data$percent.spikeIn < cfs$ERCC_cf))
    if(Mito_filter&cfs$ERCC_cf!=0) cell <- intersect(cell, which(Obj$meta.data$percent.mito < cfs$mito_cf))

    #gene filtering
    if(gene_filter){
      gene <- unname((which(Obj$feature.data$nCells >= Obj$cal.cutoffs$nCell_cutoff)))
      data <- d[gene,]
      fd<-Obj$feature.data[gene,]
      cat("Gene expressed in less than ", Obj$cal.cutoffs$nCell_cutoff,
          "cells (",(dim(d)[1]-length(gene))*100/dim(d)[1],"% genes) were filtered","\n",
          "Filtered expression matrix dimension:",dim(data),"\n")
      }else { data <- d; fd <- Obj$feature.data}

    #cell filtering
    if(cell_filter){
      data <- data[,cell]
      pd <- Obj$meta.data[cell,]
      cat("A total of ",dim(d)[2]-length(cell),
          "(",(dim(d)[2]-length(cell))*100/dim(d)[2],"%) cells were filtered","\n",
          "Filtered expression matrix dimension:",dim(data),"\n")
      }else {pd<-Obj@meta.data}

    cat("Data filtering done!","\n")
    cat("====================","\n")

  return(Obj.new)
}#end preMICA.filtering







