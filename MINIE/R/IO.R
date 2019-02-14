###Function 1:read MICA result & MICA input, output as a eset###
#' @export
readMICAoutput<-function(input_file, output_file,load_clust_label=TRUE){

	cat("Reading Input...","\n")
	d <- read.table(input_file,header = TRUE,stringsAsFactors = FALSE,quote = "")
	input<-t(d[,-1]);colnames(input)<-d$ID;rm(d)

	MICA_Result <- read.table( output_file, # MICA output text file
                          	   header = TRUE,
                          	   stringsAsFactors = FALSE)

	fd<-data.frame(row.names=rownames(input),
					geneNames=rownames(input),stringsAsFactors = FALSE)

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



###Network related functions###
###Funciton2: Generate SJAaracne input using scRNAseq data###
#' @export
SJARACNeInput_scRNAseq<-function(eset.sel,tf.ref,wd.src,grp.tag){

  cat(grp.tag,'\n')

  #fData(eset.sel)$IQR<-apply(exprs(eset.sel),1,IQR)

  # exclude genes with all zero
  eset.sel<-eset.sel[apply(exprs(eset.sel),1,function(xx){sum(xx)!=0}),]

  fData(eset.sel)$geneSymbol<-fData(eset.sel)$geneNames

  ni<-nrow(eset.sel);ni
  ns<-ncol(eset.sel);ns
  ng<-nlevels(factor(fData(eset.sel)$geneNames));ng
  tag<-paste(grp.tag,nrow(eset.sel),ng,ncol(eset.sel),sep='_');tag

  dir.cur<-file.path(wd.src,tag);dir.cur
  dir.create(dir.cur,recursive = T)
  dir.create(file.path(dir.cur,'tf'),recursive = T)

  #write exp data to exp format
  expdata<-data.frame(cbind(isoformId=featureNames(eset.sel),geneSymbol=fData(eset.sel)$geneSymbol,exprs(eset.sel)))
  f.exp<-file.path(dir.cur,paste(grp.tag,"_",ni,"_",ng,"_",ns,".exp",sep=''));f.exp
  write.table(expdata,file=f.exp,sep="\t",row.names=FALSE,quote=FALSE)

  tf.eset.sel<-subset(eset.sel,fData(eset.sel)$geneSymbol%in%tf.ref)
  dim(tf.eset.sel)
  f.tf<-file.path(dir.cur,'tf',paste(grp.tag,"_",nrow(tf.eset.sel),"_",nlevels(factor(fData(tf.eset.sel)$geneSymbol)),"_",ns,"_tf.txt",sep=''));f.tf
  cat(featureNames(tf.eset.sel),file=f.tf,sep='\n')

}

###Function8: Wrap up function for generate SJARACNe input###
#' @export
generateSJARACNeInput<-function(eset,tf.ref,wd.src,group_tag){

  if (!dir.exists(wd.src)) dir.create(wd.src,recursive = T)

  if (group_tag%in%colnames(pData(eset))){
    groups <- unique(pData(eset)[,"group_tag"])
    for (i in 1:length(groups)){
      grp.tag<-groups[i]; eset[,which(pData(eset)[,"group_tag"]==grp.tag)] -> eset.sel
      SJARACNeInput_scRNAseq(eset.sel=eset.sel,tf.ref=tf.ref,wd.src=wd.src,grp.tag=grp.tag)
    }#end for
  }else{
    stop("Lack of group info, please check your group_tag.","\n")
  }#end if
}#end function

