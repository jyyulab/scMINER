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



#' DAG_ttest
#' @description t.test(pairwise/2case)
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

#
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
                            starts_with("z"),
                            starts_with("MeanAct")
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







