##Author:chenxi.qian@stjude.org
##Stjude.YuLab


#' GetActivityFromSJARACNe
#'
#' @description Allocate network information from SJARACNe and calculate activity score for each hub genes.
#'
#' @param SJARACNe_output_path Path to SJARACNe output folder(s)
#' @param SJARACNe_input_eset Expressionset that you generate input from
#' @param group_tag a string, group name stored in pData that defines expression matrix separation
#' @param activity.method c("weighted,unweighted), default to "unweighted"
#' @param activity.norm logical, default to TRUE.
#' @param save_network_file logical, default to FALSE
#' @param save_path Path to save network file
#'
#' @return An expressionset with activity values
#'
#' @examples
#' acs.12k <- GetActivityFromSJARACNe(
#'              SJARACNe_output_path ="./",
#'              SJARACNe_input_eset = eset.12k,
#'              activity.method="unweighted",
#'              activity.norm=TRUE,
#'              group_tag = "celltype",
#'              save_network_file=TRUE, #default as false, but recommended to be TRUE
#'              save_path="./networks")
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
							   save_path=NA)
{
  eset<-SJARACNe_input_eset;
	if(!group_tag%in%colnames(pData(eset))){
	  stop('Check your group_tag please.','\n')
	}

  #retrieve networks
	output.files<-list.files(path=SJARACNe_output_path,
						pattern="consensus_network_3col_.txt",recursive = TRUE,full.names = TRUE)

  print(output.files)
	#initialize actiivty list
	acs_master<-data.frame(geneSymbol=NA,stringsAsFactors=FALSE)
	deg_master<-data.frame(geneSymbol=NA,stringsAsFactors=FALSE)

	eset<-SJARACNe_input_eset
	for( i in 1:length(output.files)){
      net.name<-gsub(SJARACNe_output_path,"",output.files[i])
  	  net.name<-gsub("\\_.*","",net.name);
  	  net.name<-gsub("[/]","",net.name)

  	  cat("Retrieve Network from ",i,net.name,"\n")
      TF.table<-read.table(output.files[i],header = TRUE,
  						stringsAsFactors = FALSE,check.names = FALSE)


      if(save_network_file)
  		{ gsc <- getGSC(tf = TF.table, sig=SIG.table)
  	 	  save(gsc,file=file.path(save_path,paste0("gsc.",netname)))
        cat("Network saved for ", net.name,"\n")}

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
