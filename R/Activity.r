##Author:chenxi.qian@stjude.org
##Stjude.YuLab
#'
#' GetActivityFromSJARACNe
#'
#' @description Allocate network information from SJARACNe and calculate activity score for each hub genes.
#'
#' @param SJARACNe_output_path Path to SJARACNe output folder(s)
#' @param SJARACNe_input_eset Expressionset that you generate input from
#' @param group_tag a string, group name stored in pData that defines expression matrix separation
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
#' \dontrun{
#' acs.12k <- GetActivityFromSJARACNe(
#'              SJARACNe_output_path ="./",
#'              SJARACNe_input_eset = eset.12k,
#'              activity.method="unweighted",
#'              activity.norm=TRUE,
#'              group_tag = "celltype",
#'              save_network_file=TRUE, #default as false, but recommended to be TRUE
#'              save_path="./networks")
#'}
#'
#' @keywords GetActivity
#' @author Chenxi Qian, \email{chenxi.qian@stjude.org}
#'
#'
#' @export
GetActivityFromSJARACNe<-function(SJARACNe_output_path=NA,
							   SJARACNe_input_eset=NA,
							   group_tag=NA,
							   activity.method="unweighted",
							   activity.norm=TRUE,
							   save_network_file=FALSE,
							   functype="tf",
							   save_path=NA){
  eset<-SJARACNe_input_eset;

  if(!group_tag%in%colnames(pData(eset))){
	  stop('Check your group_tag please.','\n')
	}

  #retrieve networks
	output.files<-list.files(path=SJARACNe_output_path,
						pattern="consensus_network_3col_.txt",recursive = TRUE,full.names = TRUE)

	if(length(output.files)==0) stop ("Please check your SJARACNe output path!",'\n')

	if(!is.null(functype)) {
    if (!functype%in%c("tf","sig")) stop("Only accept functype %in% c('tf','sig')","\n")
    else output.files<-output.files[grep(paste0("/",functype,"/"),output.files)]}

  net.names<-gsub(SJARACNe_output_path,"",output.files)
	net.names<-gsub("\\_.*","",net.names);
	net.names<-gsub("[/]","",net.names)
	celltypes<-unique(net.names)

  #initialize actiivty list
	acs_master<-data.frame(ID=NA,stringsAsFactors=FALSE)
	deg_master<-data.frame(ID=NA,stringsAsFactors=FALSE)


	for( i in 1:length(celltypes)){

	    net<-celltypes[i]
      cat("Retrieve Network from ",i,net,"\n")

      TF.table<-NULL
      SIG.table<-NULL

      f<-output.files[grep(paste0("/",net,"_"),output.files)]

      if (length(grep("/tf/",f)!=0))
        {TF.net<-NetBID2::get.SJAracne.network(network_file = f[grep("/tf/",f)])}
      if(length(grep("/sig/",f)!=0))
        {SIG.table<-NetBID2::get.SJAracne.network(network_file= f[grep("/sig/",f)])}

      if(save_network_file){
        if(!is.null(TF.table)) save(TF.net,file=file.path(save_path,paste0(net,".TF.network")))
        if(!is.null(SIG.table)) save(SIG.table,file=file.path(save_path,paste0(net,".TF.network")))
        cat("Network saved for ", net,"\n")
      }

  	  cat("Calculate Activity for ",net,"!",'\n')
  	  eset.sel<-eset[,pData(eset)[,group_tag]==net]

      acs1<-get_activity(Net = TF.net$network_dat,tag = "TF",normalize=activity.norm,
    					   eset = eset.sel, activity.method = activity.method)

      acs2<-get_activity(Net = SIG.table$network_dat,tag = "SIG",normalize=activity.norm,
                         eset = eset.sel, activity.method = activity.method)

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

  fd <-merge(fd,fData(eset),by.x="fn",by.y="row.names")
  rownames(fd)<-fd$ID

  pd <- pData(eset)

	acs.mtx <- as.matrix(acs_master[,-1])
	rownames(acs.mtx)<- acs_master$ID

	acs.mtx<-acs.mtx[,rownames(pd)]
  acs.mtx<-acs.mtx[-which(is.na(rownames(acs.mtx))),]

  acs.eset<-new("ExpressionSet",phenoData= new("AnnotatedDataFrame",pd),
          featureData=new("AnnotatedDataFrame",fd), annotation="",exprs=as.matrix(acs.mtx))

  return(acs.eset)
}#end activity function



#############################
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


####inner function to calculate activity
get_activity<-function(Net,eset,tag,exp.match=NULL, match.method=NULL,
                       es.method="mean", activity.method="weighted",
                       normalize=TRUE,test=FALSE,sep.symbol="."){
  library(dplyr)
  if(is.null(Net)) {return(NULL)}

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
    if(length(intersect(featureNames(eset),gsc[[i]]))==0) stop("Please check your featureNames!")
      else{
        eset.sel<-eset[featureNames(eset)%in%gsc[[i]],]
        gsc[[i]]<-featureNames(eset.sel)
      }

    #n=degree
    n<-length(gsc[[i]])
    name<-paste(colnames(ac)[i], n, sep='_')

    if(n>1){

      if(activity.method == 'unweighted') ac[,i]<-apply(exp[gsc[[i]],], 2, es, es.method)

      else if (activity.method == 'weighted' && is.null(exp.match)){

        fd.sel<-data.frame(fn=featureNames(eset.sel),
                           geneSymbol=fData(eset.sel)$geneSymbol,
                           stringsAsFactors=FALSE)
        tmp<-merge(fd.sel,tmp,by.x="fn",by.y="target",sort = FALSE)

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


#' DAG_ttest
#' @description inner function to do t.test(pairwise/2case) from activity matrix
#' @param d A vector of gene expression
#' @param group A vector of group information
#'
#' @export
DAG_ttest<-function(d,group){

  #d.tmp<-unlist(d[1,-1])
  d.tmp<-unlist(d[,-1])
  d.tmp<-data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)

  d.sel<-d.tmp[complete.cases(d.tmp),]
  d.sel$group<-as.factor(d.sel$group)

  n.cases<-nlevels(d.sel$group)

  cat("=")
  if(n.cases==2){
    res <- t.test(dplyr::filter(d.tmp,group=="Aim")$acs,dplyr::filter(d.tmp,group=="Ctrl")$acs)
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




#' Find differential activity genes from activity matrix
#'
#' @description  This function is a wraper of (/code {DAG_test}),
#'  which helps to conduct two_sided t.test on all genes in specific group VS Others
#'  to find differential activity genes, a table with essential statistics will be outputted.
#'
#' @param eset ExpressionSet that stores group information in pData
#' @param group_tag a character string, column name in pData(eset) that indicates group info
#' @param group_case NULL(for 1vsOthers) or a character string(pairwise),
#'  column name in pData(eset) that indicates group info
#'
#' @return output would be a table containing: t.statistics, p.value, z.score, and mean Activity value
#'
#' @seealso DAG_ttest
#' @examples
#' \dontrun{FindDAG(eset, group_tag="group")}
#'  # Try find DAG for only group 1
#'  \dontrun{FindDAG(eset, group_tag="group",group_case="group1")}
#' @export

FindDAG<-function(eset=NULL,group_tag="celltype",group_case=NULL){

  d<-data.frame(id = featureNames(eset), exprs(eset), stringsAsFactors=FALSE)
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

    cat('\n','Find Differential Activity TF for all groups!','\n')
    #do all group cases in all
    da.list <- lapply(unique(pData(eset)[,group_tag]),function(xx){

      d.label <- ifelse(pData(eset)[,group_tag] == xx, "Aim", "Ctrl") #label info

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
                        starts_with("log2FC"),
    )
  }
  return(rs)
}


##
#' @title TopMasterRegulator
#' @description Help quick pick top master regulators from previous
#' differential activity analysis results
#' @param DAG_result Output table from function FindDAG
#' @param n threshold to pick top master regulators(top n)
#' @param degree_filter filter out drivers with target number less than certain value
#' @param celltype character, output top hits are from which celltype
#' @return A list of top master regulators among different groups
#' @export
TopMasterRegulator <- function(DAG_result=res,n=5,degree_filter=c(50,500),celltype=NULL){

  rownames(DAG_result)<-DAG_result$id
  cols <- colnames(DAG_result)
  if(is.null(celltype)) celltypes<-gsub("^degree_","",cols[grep("^degree_",cols)])
  else { celltypes<-celltype }

  cat("Output top regulators for ",celltypes ,"\n")
  res<-NULL

  for (i in celltypes){

    if(!is.null(degree_filter)) {

    DAG_result<-DAG_result[which(DAG_result[,paste0("degree_",i)]> degree_filter[1] & DAG_result[,paste0("degree_",i)] <degree_filter[2]),]}

    topDriver<-DAG_result$id[sort(DAG_result[,paste0("t_",i)],decreasing=TRUE,index.return=TRUE,na.last=TRUE)$ix][1:n]

    D2print<-DAG_result[topDriver,c(1,grep(i,colnames(DAG_result)))]

    cat("Top MR for:" ,i,"\n")
    print(D2print)

    res<- c(res,topDriver)
  }
  cat("Done!")
  return(res)
}


