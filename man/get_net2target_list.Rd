% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_analysis.R
\name{get_net2target_list}
\alias{get_net2target_list}
\title{Convert a txt network file to a list}
\usage{
get_net2target_list(net_dat = NULL)
}
\arguments{
\item{net_dat}{a data frame generated from SJARANCe network file. These columns are needed: "\code{source}", "\code{target}", "\code{MI}" and "\code{spearman}".}
}
\value{
A list of source genes, for each source gene, the value is a data frame with "\code{target}", "\code{MI}" and "\code{spearman}" as the columns
}
\description{
This function is used to convert a data frame generated from the SJARACNe network file (by \code{read.delim()} or others) to a list that the names are the source genes and the values are
a data frame of target information of each source gene. This is an Iinternal function for \code{cal_Activity()}, \code{getActivity_individual()}, \code{getActivity_inBatch()}.
}
\examples{
net_data <-  read.table(system.file("extdata/demo_pbmc14k/SJARACNe/B/TF/bt100_pc001/consensus_network_ncol_.txt", package = "scMINER"),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
                        colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric"))
target_list <- get_net2target_list(net_data)

}
