##Function based MINIE
##Author:chenxi.qian@stjude.org
##Stjude.YuLab

###Differential Activity analysis#######
###Anova####
#' @export
HVG_Anova<-function(d,group){

  d.tmp<-unlist(d[,-1])
  d.tmp<-data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)
  d.sel<-d.tmp[complete.cases(d.tmp),]
  d.sel$group<-as.factor(d.sel$group)

  n.cases<-nlevels(d.sel$group)

  #print(n.cases)
  if(n.cases>1){
    rs.lm<-lm(data = d.sel,formula = acs~group)
    rs.aov<-Anova(rs.lm,Type="II",white.adjust=TRUE)
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
#' @export
FindHVG<-function(eset=acs.demo,group_tag="celltype",
                   print_screen=TRUE){

  if(!group_tag%in%colnames(pData(eset))) {
    cat('Please check your group_tag!',"\n")}

  d<-data.frame(id=featureNames(eset),exprs(eset),stringsAsFactors=FALSE)

  rs <- fData(eset);rs$id<-d$id


  zz <- file("all.Rout", open = "wt")
  sink(zz,type="message")
  da <- plyr::ddply(d,'id','HVG_Anova',group=pData(eset)[,group_tag])
  sink()
  unlink("all.Rout")

  rs <-merge(rs,da,by="id")

  if(print_screen){
    cat("Top 10 highly variable TF...")
    indx<-sort.int(da$F.value, decreasing = TRUE,na.last = NA,index.return=TRUE)$ix[1:10]
    cat(da$id[indx],"\n")
  }
  return(rs)
}


###t.test(pairwise/2case)###
#' @export
DAG_ttest<-function(d,group){

  #d.tmp<-unlist(d[1,-1])
  d.tmp<-unlist(d[,-1])
  d.tmp<-data.frame(acs=d.tmp,group=group,stringsAsFactors = FALSE)

  d.sel<-d.tmp[complete.cases(d.tmp),]
  d.sel$group<-as.factor(d.sel$group)

  n.cases<-nlevels(d.sel$group)

  #print(n.cases)
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
              MeanAct.Ctrl = unname(res$estimate[2])
              )

  }else{
    rs.t <- c(id = d$acs.id,
              t.stat=NA,
              pval = NA,
              df = NA,
              CI.low = NA,
              CI.high = NA,
              MeanAct.Aim = mean(filter(d.tmp,group=="Aim")$acs),
              MeanAct.Ctrl = mean(filter(d.tmp,group=="Ctrl")$acs))
  }
  return(rs.t)
}

###Wrapper of T-test###
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
                            starts_with("MeanAct")
                            )
    }
  return(rs)
}

##Pick top master regulator from t.test result
#' @export
TopMasterRegulator <- function(DAG_result=res,n=5,degree_filter= NULL){

  cols <- colnames(res)
  celltypes<-gsub("^degree_","",cols[grep("^degree_",cols)])
  cat("Output top regulators for",length(celltypes) ,"clusters.","\n")

  for (i in celltypes){

    if(!is.null(degree_filter)) {res<-res[which(res[,paste0("degree_",i)]> degree_filter),]}

    topMR<-res$id[sort(res[,paste0("t_",i)],decreasing=TRUE,na.last=TRUE,index.return=TRUE)$ix][1:n]

    MR2print<-filter(res[,c(1,grep(i,colnames(res)))],id%in%topMR)

    cat("Top MR for:" ,i,"\n")
    print(MR2print)
  }
  cat("Done!")
}







