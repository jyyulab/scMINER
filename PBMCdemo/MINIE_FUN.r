##Function based MINIE
##Author:chenxi.qian@stjude.org
##Stjude.YuLab

require(Biobase)
require(reshape2)
require(ggplot2)
require(RColorBrewer)
require(scales)
require(ComplexHeatmap)
require(pheatmap)
require(car)
require(dplyr)
require(openxlsx)

###Function 1:read MICA result & MICA input, output as a eset###
readMICAoutput<-function(input_file="PBMC_Demo_MICA_input_mini.txt", load_clust_label=TRUE,
						output_file="scMINER_PBMCdemo/scMINER_MICA/scMINER_PBMCdemo_MDS_4/scMINER_MICA_out/PBMCdemo.ggplot.txt"){
	
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

###Function2: Gene highlighting from eset###
gene_highlighting<-function(input_eset=eset,target,
                            ylabel="Expression",title.size=5){
  # change it to expr is ok
  input<-as.matrix(exprs(input_eset))
  indx<-which(rownames(input)%in%target)
  projection<-pData(input_eset)[colnames(input),c("tSNE_1","tSNE_2")]
  #gene expression visualized as columns
  if (length(indx)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    }else 
      {target_values <- input[indx,]}#indx = 1 
    
   	proj_target <- cbind(projection,target_values)
    proj_target_melt <- melt(proj_target, id.vars=c("tSNE_1", "tSNE_2"))
      
      p<- ggplot(proj_target_melt, aes(tSNE_1, tSNE_2)) +
        geom_point(aes(colour=value),size=0.5) + 
        facet_wrap(~variable,scales = "free")+
        labs(title="Markers highlighting",scales = "free") +
        scale_colour_gradientn(colors=rev(brewer.pal(11,"RdYlBu")))   +
        theme(plot.title = element_text(size = title.size, face = "bold"),
              axis.title = element_text(size = 10),
              legend.title = element_text(size = 10))+
        labs(x="Tsne_1",y="Tsne_2",color=ylabel)
  return(p)
}

###Function3: Gene violinplot from eset###
gene_vlnplot <- function(eset,group_tag="celltype",target,
                         ylabel="Expression",ncol=3,
                         drawquantiles=FALSE,title.size=5){
  
  # extract input information
  input <- exprs(eset)
  target <- intersect(target,rownames(input))
  label <- as.factor(pData(eset)[,group_tag]);names(label) <- rownames(pData(eset))
  # Gene expression visualized as columns
  if (length(target)!=1) {
    target_values <- t(as.matrix(input[target,]))
    label<-label[rownames(target_values)]
    df <- data.frame(target_values,cluster=label)
  }else {
    target_values<-input[target,]
    label<-label[names(target_values)]
    df <- data.frame(target_values,cluster=label)
    colnames(df)[1]<-target
  }
  
  df_melt <- melt(df, id.vars="cluster")
  
  p <- ggplot(df_melt, aes(x=cluster, y=value,fill=cluster))
  
  if (drawquantiles) p <- p + geom_violin(trim=TRUE,scale="width",weight=0.1,draw_quantiles = c(0.25,0.5,0.75),na.rm = TRUE)
  else p <- p + geom_violin(trim=TRUE,scale="width",na.rm = TRUE,size=0.1,width=0.5)
  
  p <- p + facet_wrap(~variable,scales = "free",ncol = ncol) +
    labs(x="Cluster",y=ylabel)+
    theme(axis.text.x = element_text(size=10),
          plot.title = element_text(size = title.size, face = "bold"),
          strip.background = element_rect(fill="#FFFFFF"))
  
  if (ylabel=="Activity") {p <- p + geom_boxplot(width=0.2,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)}
  
  return(p)
}

###Function6: Heatmap visualization from eset###needs editting!!!!!!!!!
gene_heatmap <- function(eset,target,
                         group_tag="label",name="log2_expression",
                         save_plot=TRUE,width=4,height=8,
                         plot_name="GeneHeatmap.png"){

  target<-intersect(target,featureNames(eset))
  exp<-exprs(eset)[target,]
  lab<-pData(eset)[,group_tag];names(lab) <- sampleNames(eset)
  
  #re-order expressionmatrix and label
  ranks<-names(sort(lab,decreasing = FALSE))
  exp.ordered<-as.matrix(exp[,ranks])
  lab.ordered<-lab[ranks]
  df<-data.frame(scMINER=lab.ordered)
  
  #Define color annotations
  n<-length(unique(lab.ordered))
  ncols <- hue_pal()(n)
  names(ncols) <- unique(lab.ordered)
  myanndf = HeatmapAnnotation(df = df,col=list(scMINER = ncols))
  mycolors = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
  
  hmp <- Heatmap(exp.ordered, col = mycolors, name = name, 
          show_row_names = TRUE,
          show_column_names = FALSE,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          top_annotation = myanndf) 
   
  if(save_plot){
    png(filename = plot_name, width=width,height=height,units ="in",res = 300)
    draw(hmp)
    dev.off()
  }
  
  return(hmp)
}

###Function7: AssignCelltypes###
AssignCellTypes.Hmp<-function(ref = NULL,eset = eset.demo,
                              save_plot = FALSE,
                              width=8.5, height=6.5, 
                              plot_name="AnnotationHeatmap.png"){
  #start from eset
  #z-normalize each sample
  exp<-apply(exprs(eset),2,std)
  #filter reference marker sets
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
  
  df<-data.frame(label=eset$label,ac);
  df<-df[,colSums(is.na(df))<nrow(df)];#remove NA columns
  df<-aggregate(.~label,df,mean)
  
  input<-t(df[,-1]) 
  
  colnames(input)<-1:length(unique(eset$label))
  myanndf<-data.frame(row.names=1:length(unique(eset$label)),scMINER=as.factor(1:length(unique(eset$label))))
  
  hmp<-pheatmap::pheatmap(input,kmeans_k = NA,
                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(30),
                          cluster_rows = FALSE,cluster_cols = FALSE,
                          scale = "row",
                          show_rownames = TRUE,
                          show_colnames = FALSE,
                          annotation_col= myanndf)
  hmp
  if(save_plot){ggsave(hmp,filename = plot_name ,device="png",width = width,height = height,dpi = 300)}
  return(hmp)
}

###Network related functions###
###Funciton9: Generate SJAaracne input using scRNAseq data###
SJAracneInput_scRNAseq<-function(eset.sel,tf.ref,wd.src,grp.tag){
  
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

###Function8: Wrap up function for generate sjaracne input###
generateSJAracneInput<-function(eset,tf.ref,wd.src,group_tag){
  
  if (!dir.exists(wd.src)) dir.create(wd.src,recursive = T)
  
  if (group_tag%in%colnames(pData(eset))){
    groups <- unique(pData(eset)[,"group_tag"])
    for (i in 1:length(groups)){
      grp.tag<-groups[i]; eset[,which(pData(eset)[,"group_tag"]==grp.tag)] -> eset.sel
      SJAracneInput_scRNAseq(eset.sel=eset.sel,tf.ref=tf.ref,wd.src=wd.src,grp.tag=grp.tag)
    }#end for
  }else{
    stop("Lack of group info, please check your group_tag.","\n")
  }#end if
}#end function

###Function9: Activity Function(from sjaracne output)###
## calculate activity and networks from network files
## path to where you want to calculate activity
GetActivityFromSJARANCE<-function(SJaracne_output_path=NA,
							   SJaracne_input_eset=NA,
							   group_tag=NA,
							   activity.method="unweighted",
							   activity.norm=TRUE,
							   save_network_file=FALSE,
							   save_path=NA)
{
  eset<-SJaracne_input_eset; 
	if(!group_tag%in%colnames(pData(eset))){
	  stop('Check your group_tag please.','\n')
	}
  
  #retrieve networks 
	output.files<-list.files(path=SJaracne_output_path,
						pattern="consensus_network_3col_.txt",recursive = TRUE,full.names = TRUE)
	
  print(output.files)
	#initialize actiivty list
	acs_master<-data.frame(geneSymbol=NA,stringsAsFactors=FALSE)
	deg_master<-data.frame(geneSymbol=NA,stringsAsFactors=FALSE)
	
	eset<-SJaracne_input_eset
	for( i in 1:length(output.files)){
      net.name<-gsub(SJaracne_output_path,"",output.files[i])
  	  net.name<-gsub("\\_.*","",net.name);
  	  net.name<-gsub("[/]","",net.name)
  
  	  cat("Retrieve Network from ",i,net.name,"\n")
      TF.table<-read.table(output.files[i],header = TRUE,
  						stringsAsFactors = FALSE,check.names = FALSE)
  
  	  if(save_network_file) 
  		{ gsc <- getGSC(tf = TF.table)
  	 	  save(gsc,file=file.path(save_path,paste0("gsc.",netname)))}
    
  	  cat("Calculate Activity for ",net.name,"!",'\n')  
  	  eset.sel<-eset[,pData(eset)[,group_tag]==net.name]
  	  
  	  fData(eset.sel)$geneSymbol<- fData(eset.sel)$geneNames
      acs.tmp<-get_activity(Net = TF.table,tag = "TF",normalize=activity.norm,
    					   eset = eset.sel, activity.method = activity.method)
      acs<-t(cbind(acs.tmp));rm(acs.tmp)

 	  #update full gene list
 	  acs.ID <- sapply(strsplit(rownames(acs),"_"),"[",1)
	  
	  acs.deg <- data.frame(geneSymbol=acs.ID,
	  						Degree=as.numeric(sapply(strsplit(rownames(acs),"_"),"[",2)),
	  						stringsAsFactors=FALSE)

  	  acs.tmp <- acs; rownames(acs.tmp)<-acs.ID

  	  acs_master<-merge(acs_master,acs.tmp,by.x="geneSymbol",by.y="row.names",all=TRUE)
  	  deg_master<-merge(deg_master,acs.deg,by="geneSymbol",all=TRUE)

  	  colnames(deg_master)[i+1]<-paste0("degree_",net.name)
  	  
  	  cat("Activity Done!!","\n")
  	  rm(acs)
  	
  	  cat('==============================================',"\n")
  	  gc()
	
	}#end for
	
	# generate acs expression set
	deg_master<-filter(deg_master,!is.na(geneSymbol))
	fd <- data.frame(row.names=deg_master$geneSymbol,
	                 geneSymbol=sapply(strsplit(deg_master$geneSymbol,"\\."),"[",1),
	                 FuncType=sapply(strsplit(deg_master$geneSymbol,"\\."),"[",2),
	                 deg_master[,-1],stringsAsFactors = FALSE)
	pd <- pData(eset)

	acs.mtx <- as.matrix(acs_master[,-1])
	rownames(acs.mtx)<- acs_master$geneSymbol

	acs.mtx<-acs.mtx[,rownames(pd)]
  acs.mtx<-acs.mtx[-which(is.na(rownames(acs.mtx))),]
	
  acs.eset<-new("ExpressionSet",phenoData= new("AnnotatedDataFrame",pd),
          featureData=new("AnnotatedDataFrame",fd), annotation="",exprs=as.matrix(acs.mtx))

  return(acs.eset)
}#end activity function 

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

get_activity<-function(Net,eset,tag,exp.match=NULL, match.method=NULL, 
                       es.method="mean", activity.method="weighted", 
                       normalize=TRUE,test=FALSE,sep.symbol="."){
  library(dplyr)  
  if(is.null(exp.match)) src<-unique(Net$source)
  else src<-unique(Net$source.symbol)
  
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
    if(is.null(exp.match)){
      tmp<-filter(Net, Net$source==src[i]);tmp<-tmp[!duplicated(tmp$target),]
      gsc[[i]]<-unique(as.character(tmp$target))#network from file
      
    }
    else{
      tmp<-filter(Net, Net$source.symbol==src[i]);tmp<-tmp[!duplicated(tmp$target.symbol),]
      gsc[[i]]<-unique(as.character(tmp$target.symbol))
      gsc[[i]]<-as.character(exp.match[toupper(exp.match$symbol.new) %in% toupper(gsc[[i]]), "exp.rowname"])
    }
    
    #update the overlap between NetBID based geneset and real expression data
    if(is.null(eset@featureData@data$geneSymbol))stop("Please check your geneSymbol!")
    
    else{
      eset.sel<-eset[eset@featureData@data$geneSymbol%in%gsc[[i]],]
      gsc[[i]]<-featureNames(eset.sel)
    }
    
    #n=degree
    n<-length(gsc[[i]])
    name<-paste(colnames(ac)[i], n, sep='_')
    
    if(n>1){
      
      if(activity.method == 'unweighted') ac[,i]<-apply(exp[gsc[[i]],], 2, es, es.method)
      
      else if (activity.method == 'weighted' && is.null(exp.match)){
        
        fd.sel<-data.frame(fn=featureNames(eset.sel),geneSymbol=fData(eset.sel)$geneSymbol,stringsAsFactors=FALSE)
        tmp<-merge(fd.sel,tmp,by.x="geneSymbol",by.y="target",sort = FALSE)
        
        #if(!all(rownames(exp[gsc[[i]],])==fd.sel$fn)) stop("MI is not coordinate with Feature names! \n")
        
        tmp$p.sign<-sign(as.numeric(tmp$spearman))
        tmp$p.sign[tmp$p.sign == 0]<-1
        tmp$MI.sign<-as.numeric(tmp$MI)*(tmp$p.sign)
        
        mat<-t(exp[gsc[[i]],])%*%(tmp$MI.sign)
        #mat<-t(exp[gsc[[i]],])%*%tmp$MI
        MI.sum<-sum(tmp$MI)
        ac[,i]<-mat[,1]/MI.sum
        
        if (test)
        {ac[,i]<-mat[,1]}
      }
      
      
      else if(activity.method == 'weighted' && !is.null(exp.match) && match.method == 'topMI'){
        tmp<-arrange(tmp, target.symbol, desc(MI))
        tmp<-tmp[!duplicated(tmp$target.symbol),]
        tmp$p.sign<-sign(as.numeric(tmp$spearman))
        tmp$p.sign[tmp$p.sign == 0]<-1
        tmp$MI.sign<-as.numeric(tmp$MI)*tmp$p.sign
        match.mat<-exp.match[exp.match$exp.rowname %in% gsc[[i]],]
        match.mat<-left_join(match.mat, tmp, by =c('symbol.new'='target.symbol'))
        match.mat<-na.omit(match.mat)
        gsc[[i]] <- as.character(match.mat$exp.rowname)
        mat<-t(exp[gsc[[i]],])%*%match.mat$MI.sign
        ac[,i]<-mat[,1]/length(gsc[[i]])
      }
      else if(activity.method == 'weighted' && !is.null(exp.match) && match.method == 'avg.signed.MI'){
        
        tmp$p.sign<-sign(as.numeric(tmp$spearman))
        tmp$MI.sign<-as.numeric(tmp$MI)*tmp$p.sign
        tmp$p.sign[tmp$p.sign == 0]<-1
        tmp<-aggregate(MI.sign ~ target.symbol, data = tmp, mean)
        #generate MI.sign vector match to the exp rows
        match.mat<-exp.match[exp.match$exp.rowname %in% gsc[[i]],]
        match.mat<-left_join(match.mat, tmp, by =c('symbol.new'='target.symbol'))
        match.mat<-na.omit(match.mat)
        gsc[[i]] <- as.character(match.mat$exp.rowname)
        mat<-t(exp[gsc[[i]],])%*%match.mat$MI.sign
        ac[,i]<-mat[,1]/length(gsc[[i]])
      }
      colnames(ac)[i]<-name
    }
  }
  return(ac)
}

###Differential Activity analysis#######

###Anova####
HVG_Anova<-function(d,group){
  
  d.tmp<-unlist(d[,-1])
  d.tmp<-data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)
  d.sel<-d.tmp[complete.cases(d.tmp),]
  d.sel$group<-as.factor(d.sel$group)
  
  n.cases<-nlevels(d.sel$group)
  
  #print(n.cases)
  if(n.cases>1){
    rs.lm<-lm(data = d.sel,formula = acs~group)
    rs.aov<-suppressWarnings(Anova(rs.lm,Type="II",white.adjust=TRUE))
    rs.aov<-c(id=d$acs.id,F.value = (rs.aov$F)[1],pval = (rs.aov$`Pr(>F)`)[1])
  }else{
    rs.aov<-c(id=d$acs.id,F.value = NA, pval= NA)  
  }
  
  rs.acs<-tapply(d.tmp$acs, list(d.tmp$group), mean)
  names(rs.acs)<-paste0("MeanAct_",names(rs.acs))
  rs<-c(rs.aov,rs.acs)
  return(rs)
}

###Wrapper of Anova###
FindHVG<-function(eset = acs.demo,group_tag = "celltype",
                   print_screen = TRUE){
  
  if(!group_tag%in%colnames(pData(eset))) {
    cat('Please check your group_tag!',"\n")}
  
  d <- data.frame(id = featureNames(eset),exprs(eset),stringsAsFactors = FALSE)
  
  rs <- fData(eset); rs$id <- d$id
  
  
  zz <- file("all.Rout", open = "wt")
  sink(zz,type="message")
  da <- plyr::ddply(d,'id','HVG_Anova',group=pData(eset)[,group_tag])
  sink()
  unlink("all.Rout")
  
  rs <-merge(rs,da,by = "id")
  
  if(print_screen){
    
    cat("Top 10 highly variable TF...")
    
    indx<-sort.int(da$F.value, decreasing = TRUE,na.last = NA,index.return=TRUE)$ix[1:10]
    
    cat(da$id[indx],"\n")
  }
  
  return(rs)
}

###t.test(pairwise/2case)###
DAG_ttest <- function(d,group){
  
  #d.tmp<-unlist(d[1,-1])
  d.tmp <- unlist(d[,-1])
  d.tmp <- data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)
  d.sel <- d.tmp[complete.cases(d.tmp),]
  d.sel$group <- as.factor(d.sel$group)
  
  n.cases <- nlevels(d.sel$group)
  
  #print(n.cases)
  if(n.cases == 2){
    res <- t.test(dplyr::filter(d.tmp,group=="Aim")$acs,dplyr::filter(d.tmp,group=="Ctrl")$acs)
    pval.t <- res$p.value
    t.stat <- res$statistic
    #extract results
    rs.t <- c(id= d$acs.id,
              t.stat,
              pval = pval.t, 
              res$parameter,#df
              CI.low = res$conf.int[1],
              CI.high = res$conf.int[2],
              MeanAct.Aim = unname(res$estimate[1]),
              MeanAct.Ctrl = unname(res$estimate[2]))
    	if (pval.t == 0) {rs.t <- c(rs.t, z = sign(t.stat)*40)}
    	 else {rs.t <- c(rs.t, z  = sign(pval.t) * abs(qnorm(pval.t/2)))}
  
  }else{
    rs.t <- c(id = d$acs.id,
              t.stat = NA,
              pval = NA,
              z  = NA,
              df = NA,
              CI.low = NA,
              CI.high = NA,
              MeanAct.Aim = mean(filter(d.tmp,group=="Aim")$acs),
              MeanAct.Ctrl = mean(filter(d.tmp,group=="Ctrl")$acs))
  }
  return(rs.t)
}

###Wrapper of T-test###
FindDAG <- function(eset=acs.demo,group_tag="celltype",group_case=NULL){
  
  d <- data.frame(id = featureNames(eset), exprs(eset), stringsAsFactors=FALSE)
  rs <- fData(eset);rs$id <- d$id
  
  #3 cases 1) all 1vsOther 
  if(!group_tag%in%colnames(pData(eset))) {
    stop('Please check your group_tag.',"\n")}
  
  if(!is.null(group_case)){
    if(!group_case%in%pData(eset)[,"group_tag"]){
      stop('Please check your group_case.',"\n")
    }
    
    d.label <- ifelse(pData(eset)[,group_tag]==group_case,"Aim","Ctrl") #label info
    
    da <- plyr::ddply(d,'id','DAG_ttest',group=d.label)
    
    rs <- merge(rs,da,by="id")
    
  }else{
    
    cat('Find Differential Activity TF for all groups!')
      #do all group cases in all
    da.list <- lapply(unique(pData(eset)[,group_tag]),function(xx){
      
      d.label <- ifelse(pData(eset)[,group_tag] == xx, "Aim", "Ctrl") #label info
      
      da <- plyr::ddply(d,'id','DAG_ttest',group=d.label)
      
      da <- da[,setdiff(colnames(da),"MeanAct.Ctrl")]
    
      colnames(da)[-1] <- paste0(colnames(da)[-1],'_',xx)      
      
      return(da)})
    
    rs.tmp <- as.data.frame(da.list,stringsAsFactors=FALSE)  
    rs.full <- merge(rs,rs.tmp,by='id');rm(rs.tmp)
    
    rs <- dplyr::select(rs.full,
                            id:FuncType,
                            starts_with("degree"),
                            starts_with("t"),
                            starts_with("pval"),
                            starts_with("z"),
                            starts_with("MeanAct"))
    }
  return(rs)
}



