##Inner function of network visualization
##Author: xinran.dong@stjude.org
##Stjude.YuLab

#' getGSC
#'
#' @param sig signaling genes network table
#' @param tf transcription factors network table
#'
#' @return a list of master regulators and their targets
#' @export
getGSC<-function(sig=NULL,tf=NULL){
  sig.gsc<-NULL
  tf.gsc<-NULL

  if(!is.null(sig)){
    src <-unique(sig$source)
    n<-length(src)

    sig.gsc <-vector("list",n)
    for(i in 1:n){
      tag<-src[i]
      tmp<-filter(sig,sig$source==tag)
      names(sig.gsc)[i] <-paste(tag,"SIG",sep = ".")
      sig.gsc[[i]] <-as.character(tmp$target)}
  }

  if (!is.null(tf)){
    src <-unique(tf$source)
    n<-length(src)
    tf.gsc <-vector("list",n)
    for(i in 1:n){
      tag<-src[i]
      tmp<-filter(tf,tf$source==tag)
      names(tf.gsc)[i] <-paste(tag,"TF",sep = ".")
      tf.gsc[[i]] <-as.character(tmp$target)}}

  gsc<-c(sig.gsc,tf.gsc)
  return(gsc)
}


