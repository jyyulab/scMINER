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
#' @description scRNA-seq filtering function
#' @param Obj the list object outputted by draw.scRNAseq.QC
#' @param cutoffs a list outputted by draw.scRNAseq.QC, if NULL, manual input will be required
#' @param gene_filter logical; or a numerical number, indicating lower threshold for gene filtering based on how many non-zero cells each gene expressed in
#' @param ERCC_filter logical; a numerical number, indicating upper threshold put on ERCC percentage for cell filtering
#' @param Mito_filter logical;a numerical number, indicating upper threshold put on Mitochondrial gene expression fraction for cell filtering
#' @param nGene_filter logical; a numerical number, indicating lower threshold put on number of gene expression in each cell for cell filtering
#' @param nUMI_filter logical;a vector of two numerical number, indicating lower threshold and upper threshold put on number of total UMI for cell filtering
#'
#' @return A Sparse expression set
#'
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

    if(is.null(cutoffs)){
      if(all(is.logical(gene_filter,ERCC_filter,Mito_filter,nGene_filter,nUMI_filter))){
        stop("No filtering will be conducted due to lack of numerical threshold input.","\n",
             "You can take default cutoff calculated by function draw.scRNAseq.QC and feed them into cutoffs,
             or input them manually in each function parameters.")
        }else{
          if(any(isTRUE(gene_filter),isTRUE(ERCC_filter),isTRUE(Mito_filter),is.TRUE(nGene_filter),is.TRUE(nUMI_filter)))
            stop("When cutoffs=NULL, please indicate numerical threshold instead of 'TRUE' if you want to do filtering on that particular criteria !","\n")
          cfs<-list()
          if (gene_filter) cfs$nCell_cutoff=gene_filter
          if (nGene_filter) cdf$nGene_cf=nGene_filter
          if (nUMI_filter) cdf$umi_cf_lo=nUMI_filter[1];cdf$umi_cf_hi=nUMI_filter[2]
          if (is.null(ERCC_filter)) ERCC_filter=FALSE else cdf$ERCC_cf=ERCC_filter
        }
    }else{cfs<-cutoffs}
    if(!nUMI_filter) cdf$umi_cf_lo=0; cdf$umi_cf_hi=Inf
    cell <- which((SparseEset$nUMI.total > cfs$umi_cf_lo) & (SparseEset$nUMI.total < cfs$umi_cf_hi))
    if(nGene_filter) cell <- intersect(cell, which(SparseEset$nGene > cfs$nGene_cf))
    if(ERCC_filter&cfs$ERCC_cf!=0) cell <- intersect(cell, which(SparseEset$percent.spikeIn < cfs$ERCC_cf))
    if(Mito_filter&cfs$Mito_filter!=0) cell <- intersect(cell, which(SparseEset$percent.mito < cfs$mito_cf))

    #gene filtering
    if(gene_filter){
      gene <- unname((which(fData(SparseEset)$nCells >= cfs$nCell_cutoff)))
      eset.sel<-SparseEset[gene,]

      cat("Gene expressed in less than ", cfs$nCell_cutoff,
          "cells (",(dim(SparseEset)[1]-length(gene))*100/dim(SparseEset)[1],"% genes) were filtered","\n",
          "Filtered expression matrix dimension:",dim(data),"\n")
      }else { eset.sel <- SparseEset }

    #cell filtering
    if(cell_filter){
      eset.sel2 <- eset.sel[,cell]
      cat("A total of ",dim(SparseEset)[2]-length(cell),
          "(",(dim(SparseEset)[2]-length(cell))*100/dim(SparseEset)[2],"%) cells were filtered","\n",
          "Filtered expression matrix dimension:",dim(eset.sel2),"\n")
      }else {eset.sel2<-eset.sel}

    cat("Data filtering done!","\n")
    cat("====================","\n")

  return(eset.sel2)
}#end preMICA.filtering

